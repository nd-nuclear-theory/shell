/****************************************************************

  relative_xform.cpp

  Patrick J. Fasano, University of Notre Dame.

****************************************************************/


#include "relative_xform.h"

#include "fmt/format.h"
#include "mcutils/eigen.h"
#include "quadrature/mesh.h"
#include "quadrature/oscillator_wf.h"
#include "spline/spline.h"

namespace relative {

  ////////////////////////////////////////////////////////////////
  // oscillator radial matrix
  ////////////////////////////////////////////////////////////////


  void PopulateRelativeXformMatrix(
      std::size_t source_dim, std::size_t target_dim, int L,
      double b_ratio, int num_steps,
      basis::OperatorBlock<double>& xform_matrix
    )
  // Construct transformation matrix for dilation between oscillator bases.
  //
  // Arguments:
  //   source_dim, target_dim (input): dimensions of radial subspaces
  //   L (input): orbital angular momentum quantum number
  //   b_ratio (input): ratio of target to source oscillator length parameter
  //
  // Returns:
  //   transformation matrix
  {
    // short circuit on non-dilation transform (b_ratio = 1)
    if (b_ratio == 1.0)
    {
      xform_matrix = basis::OperatorBlock<double>::Identity(source_dim,target_dim);
      return;
    }

    // symmetrically distribute ratio between source and target
    assert(b_ratio>0);
    const double source_b = 1./std::sqrt(b_ratio);
    const double target_b = std::sqrt(b_ratio);

    // convert number of steps to number of points
    const int num_size = num_steps+1;

    const Eigen::ArrayXd z_grid = shell::UniformMesh(num_size);
    const Eigen::ArrayXd r_grid = shell::TransformedCoordinate<+1>(1.0, z_grid);
    const Eigen::ArrayXd jacobian_grid =
        shell::TransformedCoordinateJacobian<+1>(1.0, z_grid);

    const auto source_oscillator_functions = shell::OscillatorWaveFunctions<+1>(
        source_dim-1, L, source_b, r_grid, 0
      );
    const auto target_oscillator_functions = shell::OscillatorWaveFunctions<+1>(
        target_dim-1, L, target_b, r_grid, 0
      );

    // initialize matrix
    xform_matrix = basis::OperatorBlock<double>::Zero(source_dim,target_dim);

    // populate nonzero entries
    #pragma omp parallel for collapse(2) default(none)                  \
      shared(                                                           \
          xform_matrix, z_grid, jacobian_grid, source_dim, target_dim,  \
          source_oscillator_functions, target_oscillator_functions      \
        )
    for (std::size_t source_n=0; source_n<source_dim; ++source_n)
      for (std::size_t target_n=0; target_n<target_dim; ++target_n)
      {
        // evaluate radial integral
        Eigen::ArrayXd integrand =
          source_oscillator_functions[source_n]
          * target_oscillator_functions[target_n]
          * jacobian_grid;
        shell::FixEndpointSingularitiesToZero(integrand);
        double radial_integral = spline::CubicIntegrate(z_grid, integrand);
        if (std::isnan(radial_integral))
          throw std::runtime_error("radial_integral is nan");
        xform_matrix(source_n,target_n) = radial_integral;
      }
  }

  void TransformRelativeOperatorLSJT(
      const basis::RelativeOperatorParametersLSJT& source_operator_parameters,
      const basis::RelativeSpaceLSJT& source_relative_space,
      const std::array<basis::RelativeSectorsLSJT,3>& source_relative_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& source_relative_component_blocks,
      const basis::RelativeOperatorParametersLSJT& target_operator_parameters,
      const basis::RelativeSpaceLSJT& target_relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& target_relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& target_relative_component_blocks,
      double b_ratio,
      int num_steps,
      bool verbose
    )
  {
    // write diagnostics
    if (verbose)
    {
      std::cout << "Transforming relative operator..." << std::endl;
      std::cout
        << "  Source operator properties:"
        << " J0 " << source_operator_parameters.J0
        << " g0 " << source_operator_parameters.g0
        << " T0_min " << source_operator_parameters.T0_min
        << " T0_max " << source_operator_parameters.T0_max
        << " symmetry " << int(source_operator_parameters.symmetry_phase_mode)
        << std::endl
        << "  Source truncation:"
        << " Nmax " << source_operator_parameters.Nmax
        << " Jmax " << source_operator_parameters.Jmax
        << std::endl;
      std::cout
        << "  Target operator properties:"
        << " J0 " << target_operator_parameters.J0
        << " g0 " << target_operator_parameters.g0
        << " T0_min " << target_operator_parameters.T0_min
        << " T0_max " << target_operator_parameters.T0_max
        << " symmetry " << int(target_operator_parameters.symmetry_phase_mode)
        << std::endl
        << "  Target truncation:"
        << " Nmax " << target_operator_parameters.Nmax
        << " Jmax " << target_operator_parameters.Jmax
        << std::endl;
    }

    // check compatibility of source and target parameters
    assert(target_operator_parameters.J0==source_operator_parameters.J0);
    assert(target_operator_parameters.g0==source_operator_parameters.g0);
    assert(target_operator_parameters.T0_max<=source_operator_parameters.T0_max);
    assert(target_operator_parameters.T0_min>=source_operator_parameters.T0_min);
    assert(target_operator_parameters.Jmax<=source_operator_parameters.Jmax);

    // write diagnostics
    if (verbose)
    {
      std::cout << "  Source matrix elements:";
      for (int T0=source_operator_parameters.T0_min; T0<=source_operator_parameters.T0_max; ++T0)
        std::cout << " " << basis::UpperTriangularEntries(source_relative_component_sectors[T0]);
      std::cout << std::endl;
      std::cout << "  Allocated:";
      for (int T0=source_operator_parameters.T0_min; T0<=source_operator_parameters.T0_max; ++T0)
        std::cout << " " << basis::AllocatedEntries(source_relative_component_blocks[T0]);
      std::cout << std::endl;
    }

    // construct target indexing and zero initialize matrices
    basis::ConstructZeroOperatorRelativeLSJT(
        target_operator_parameters,
        target_relative_space,
        target_relative_component_sectors,
        target_relative_component_blocks
      );

    // transformation matrix storage
    std::unordered_map<std::size_t,std::tuple<std::size_t,basis::OperatorBlock<double>>> xform_matrix_map;

    // precompute transformation matrices
    if (verbose)
    {
      std::cout << "  Constructing overlap matrices..." << std::flush;
    }
    for (std::size_t target_subspace_index=0; target_subspace_index<target_relative_space.size(); ++target_subspace_index)
    {
      const auto& target_subspace = target_relative_space.GetSubspace(target_subspace_index);
      const std::size_t source_subspace_index = source_relative_space.LookUpSubspaceIndex(target_subspace.labels());

      // default initialize xform matrix
      basis::OperatorBlock<double> xform_matrix;

      // populate xform matrix if source subspace exists
      if (source_subspace_index != basis::kNone)
      {
        const auto& source_subspace = source_relative_space.GetSubspace(source_subspace_index);
        PopulateRelativeXformMatrix(
            source_subspace.size(), target_subspace.size(), target_subspace.L(),
            b_ratio, num_steps, xform_matrix
          );
      }
      xform_matrix_map[target_subspace_index] = {source_subspace_index, xform_matrix};
      std::cout << "." << std::flush;
    }
    if (verbose)
    {
      std::cout << "done." << std::endl;
    }


    if (verbose)
    {
      std::cout << "  Applying transformation..." << std::flush;
    }
    for (int T0 = target_operator_parameters.T0_min; T0 <= target_operator_parameters.T0_max; ++T0)
    {
      // select T0 component of both source and target
      const basis::RelativeSectorsLSJT& source_sectors = source_relative_component_sectors[T0];
      const basis::OperatorBlocks<double>& source_blocks = source_relative_component_blocks[T0];
      const basis::RelativeSectorsLSJT& target_sectors = target_relative_component_sectors[T0];
      basis::OperatorBlocks<double>& target_blocks = target_relative_component_blocks[T0];

      for (std::size_t target_sector_index=0; target_sector_index<target_sectors.size(); ++target_sector_index)
      {
        // set up alias for indexing and get target subspace indices
        const auto& target_sector = target_sectors.GetSector(target_sector_index);
        const auto& target_bra_subspace = target_sector.bra_subspace();
        const auto& target_ket_subspace = target_sector.ket_subspace();
        const auto target_bra_subspace_index = target_sector.bra_subspace_index();
        const auto target_ket_subspace_index = target_sector.ket_subspace_index();

        // get source subspace indices and xform matrices
        const auto& [source_bra_subspace_index, bra_xform_matrix]
          = xform_matrix_map.at(target_bra_subspace_index);
        const auto& [source_ket_subspace_index, ket_xform_matrix]
          = xform_matrix_map.at(target_ket_subspace_index);

        // get source sector
        const auto source_sector_index = source_sectors.LookUpSectorIndex(
            source_bra_subspace_index, source_ket_subspace_index
          );

        // print progress indicator
        if (verbose)
          std::cout << "." << std::flush;
        // skip if source sector does not exist
        if (source_sector_index==basis::kNone)
          continue;

        // set up alias for indexing and get target subspace indices
        const auto& source_sector = source_sectors.GetSector(source_sector_index);
        const auto& source_bra_subspace = source_sector.bra_subspace();
        const auto& source_ket_subspace = source_sector.ket_subspace();

        // apply transformation and store
        basis::OperatorBlock<double> source_block = source_blocks[source_sector_index];
        if (source_sector.IsDiagonal())
          mcutils::CompleteLowerTriangle(source_block);
        // std::cout << source_block << std::endl;
        target_blocks[target_sector_index]
          = bra_xform_matrix.transpose() * source_block * ket_xform_matrix;
        // std::cout << target_blocks[target_sector_index] << std::endl;
        // print progress indicator
      }
    }
    if (verbose)
      std::cout << "done." << std::endl;

    // write diagnostics
    if (verbose)
    {
      std::cout << "  Target matrix elements:";
      for (int T0=target_operator_parameters.T0_min; T0<=target_operator_parameters.T0_max; ++T0)
        std::cout << " " << basis::UpperTriangularEntries(target_relative_component_sectors[T0]);
      std::cout << std::endl;
      std::cout << "  Allocated:";
      for (int T0=target_operator_parameters.T0_min; T0<=target_operator_parameters.T0_max; ++T0)
        std::cout << " " << basis::AllocatedEntries(target_relative_component_blocks[T0]);
      std::cout << std::endl;
    }

  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
}  // namespace relative
