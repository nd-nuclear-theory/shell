/****************************************************************
  tbme_mapping.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "tbme/tbme_mapping.h"

#include <sstream>
#include <utility>
#include <tuple>

#include "cppformat/format.h" // for debugging


namespace shell {

  TwoBodyMapping::TwoBodyMapping(
      const basis::OrbitalSpacePN& source_orbital_space,
      const basis::TwoBodySpaceJJJPN& source_space,
      const basis::OrbitalSpacePN& target_orbital_space,
      const basis::TwoBodySpaceJJJPN& target_space
    )
  {
    // initialize informational flags for mapping properties
    domain_subspaces_covered = true;
    domain_states_covered = true;
    range_subspaces_covered = true;
    range_states_covered = true;

    // set up orbital mapping
    assert(source_orbital_space.size()==2);
    assert(target_orbital_space.size()==2);
    orbital_mapping.resize(2);

    // iterate over orbital subspaces
    for (int orbital_subspace_index=0; orbital_subspace_index<2; ++orbital_subspace_index)
      {
        // set up subspace aliases
        const basis::OrbitalSubspacePN& source_orbital_subspace
          = source_orbital_space.GetSubspace(orbital_subspace_index);
        const basis::OrbitalSubspacePN& target_orbital_subspace
          = target_orbital_space.GetSubspace(orbital_subspace_index);

        // set up vector to hold index mappings for this subspace
        orbital_mapping[orbital_subspace_index].resize(source_orbital_subspace.size());

        // find mappings for source orbitals
        for (
            int source_orbital_index=0;
            source_orbital_index<source_orbital_subspace.size();
            ++source_orbital_index
          )
          {
            // look up corresponding target state index -- old version for OrDie lookup
            // int target_orbital_index = TwoBodyMapping::kNone;
            // const typename basis::OrbitalSubspacePN::StateLabelsType state_labels
            //   = source_orbital_subspace.GetStateLabels(source_orbital_index);
            // bool target_orbital_found = target_orbital_subspace.ContainsState(state_labels);
            // if (target_orbital_found)
            //     target_orbital_index = target_orbital_subspace.LookUpStateIndex(state_labels);

            // look up corresponding target state index
            //
            // The target index value will be basis::kNone if no
            // target matches the source labels.
            int target_orbital_index = target_orbital_subspace.LookUpStateIndex(
              source_orbital_subspace.GetStateLabels(source_orbital_index)
              );
            orbital_mapping[orbital_subspace_index][source_orbital_index] = target_orbital_index;
          }
      }

    // set up two_body mapping
    subspace_mapping.resize(source_space.size());
    state_mapping.resize(source_space.size());
    int num_target_subspaces_found = 0;  // used to track range coverage

    // iterate over two-body subspaces
    for (int source_subspace_index=0; source_subspace_index<source_space.size(); ++source_subspace_index)
      {
        // set up source subspace alias
        const basis::TwoBodySubspaceJJJPN& source_subspace = source_space.GetSubspace(source_subspace_index);

        // look up target subspace -- old version for "OrDie" lookup
        // bool target_subspace_found = target_space.ContainsSubspace(source_subspace.labels());
        // if (!target_subspace_found)
        //   // handle missing target subspace
        //   {
        //     orbital_mapping.emplace_back();
        //     continue;
        //   }
        // int target_subspace_index = target_space.LookUpSubspaceIndex(source_subspace.labels());
        // const basis::TwoBodySubspaceJJJPN& target_subspace = source_space.GetSubspace(target_subspace_index);

        // look up corresponding target subspace index
        //
        // The target index value will be basis::kNone if no
        // target matches the source labels.
        int target_subspace_index = target_space.LookUpSubspaceIndex(source_subspace.labels());
        subspace_mapping[source_subspace_index] = target_subspace_index;
        // std::cout << fmt::format("Subspace {} -> {}",source_subspace_index,target_subspace_index) << std::endl;

        // do diagnostic record keeping
        if (target_subspace_index==basis::kNone)
          {
            domain_subspaces_covered = false;
            domain_states_covered = false;
          }
        else
          ++num_target_subspaces_found;

        // short circuit if target subspace missing
        if (target_subspace_index==basis::kNone)
            continue;

        // set up alias for target subspace
        const basis::TwoBodySubspaceJJJPN& target_subspace = target_space.GetSubspace(target_subspace_index);

        // set up vector to hold index mappings for this subspace
        state_mapping[source_subspace_index].resize(source_subspace.size());
        int num_target_states_found = 0;  // used to track range coverage

        // find mappings for source two-body states
        for (
            int source_state_index=0;
            source_state_index<source_subspace.size();
            ++source_state_index
          )
          {
            // remap state labels to use target orbital indexing
            const basis::TwoBodyStateJJJPN source_state(source_subspace,source_state_index);
            const basis::OrbitalStatePN source_orbital1 = source_state.GetOrbital1();
            const basis::OrbitalStatePN source_orbital2 = source_state.GetOrbital2();
            int target_orbital1_index =
              orbital_mapping[static_cast<int>(source_orbital1.orbital_species())][source_orbital1.index()];
            int target_orbital2_index =
              orbital_mapping[static_cast<int>(source_orbital2.orbital_species())][source_orbital2.index()];

            // canonicalize orbital indices for new space
            int relative_phase = 1;
            if (
              (source_orbital1.orbital_species() == source_orbital2.orbital_species())
              &&
              (target_orbital1_index > target_orbital2_index)
            )
            {
              std::swap(target_orbital1_index, target_orbital2_index);
              // see Suhonen (8.30)
              // TODO is this phase correct for non-scalar operators?
              int grel = static_cast<int>(source_orbital1.j() + source_orbital2.j()) + source_state.J();
              relative_phase = std::pow(-1, grel + 1);
            }

            // construct labels in target subspace
            const typename basis::TwoBodyStateJJJPN::StateLabelsType target_state_labels(
              target_orbital1_index, target_orbital2_index
              );

            // look up corresponding target state index
            //
            // "Missing" target orbitals will be indexed by basis::kNone.
            int target_state_index = target_subspace.LookUpStateIndex(target_state_labels);
            if (target_state_index==basis::kNone)
              relative_phase = 0;
            state_mapping[source_subspace_index][source_state_index]
              = std::tuple<int,int>(target_state_index, relative_phase);

            // do diagnostic record keeping
            if (target_state_index==basis::kNone)
              domain_states_covered = false;
            else
              ++num_target_states_found;

            // std::cout << fmt::format(
            //     "  subspace {} source {} target {}",
            //     source_subspace_index,source_state_index,target_state_index
            //   )
            //           << std::endl;

          }

        // check if all target states covered for this subspace
        range_states_covered &= (num_target_states_found==target_subspace.size());
      }

    // check target subspace coverage
    range_subspaces_covered &= (num_target_subspaces_found==target_space.size());
    range_states_covered &= range_subspaces_covered;

  }

  std::string TwoBodyMapping::DebugStr() const
  {
    std::ostringstream os;

    // dump orbital mapping
    os << "Orbitals" << std::endl;
    for (int orbital_subspace_index=0; orbital_subspace_index<2; ++orbital_subspace_index)
        for (
            int source_orbital_index=0;
            source_orbital_index<orbital_mapping[orbital_subspace_index].size();
            ++source_orbital_index
          )
          {
            int target_orbital_index = orbital_mapping[orbital_subspace_index][source_orbital_index];
            os << fmt::format(
                "  subspace {:3} source {:3} target {:3}",
                orbital_subspace_index,source_orbital_index,target_orbital_index
              )
               << std::endl;
          }

    // dump two-body state mapping
    os << "Two-body states" << std::endl;
    for (int source_subspace_index=0; source_subspace_index<state_mapping.size(); ++source_subspace_index)
      {
        // indicate subspace mapping
        int target_subspace_index = subspace_mapping[source_subspace_index];
        os << fmt::format(
            "  subspace {:3} => {:3}",
            source_subspace_index,target_subspace_index
          )
           << std::endl;

        // skip unmapped source subspace
        if (target_subspace_index==basis::kNone)
          continue;

        // dump mappings within subspace
        for (
            int source_state_index=0;
            source_state_index<state_mapping[source_subspace_index].size();
            ++source_state_index
          )
          {
            int target_state_index, relative_phase;
            std::tie(target_state_index, relative_phase)
              = state_mapping[source_subspace_index][source_state_index];
            os << fmt::format(
                "  subspace {:3} source {:3} target {:3} relative phase {:3}",
                source_subspace_index,source_state_index,target_state_index,relative_phase
              )
               << std::endl;
          }
      }

    // dump status flags
    os << fmt::format(
        "Diagnostics:\n"
        "  domain_subspaces_covered {} domain_states_covered {}\n"
        "  range_subspaces_covered {} range_states_covered {}",
        domain_subspaces_covered, domain_states_covered, range_subspaces_covered, range_states_covered
      )
              <<std::endl;

    return os.str();

  }

Eigen::MatrixXd
RemappedMatrixJJJPN(
      const basis::TwoBodySectorsJJJPN::SectorType& source_sector,
      const basis::TwoBodySectorsJJJPN::SectorType& target_sector,
      const shell::TwoBodyMapping& two_body_mapping,
      const Eigen::MatrixXd& source_matrix
    )
{
  // Note: Could restructure loop to remove redundant lookups if needed.

  // zero initialize target
  //
  // Target matrix might not be covered by source matrix elements, so zero padding is important.
  Eigen::MatrixXd target_matrix
    =Eigen::MatrixXd::Zero(target_sector.bra_subspace().size(),target_sector.ket_subspace().size());

  // copy matrix elements
  for (int source_bra_index=0; source_bra_index<source_sector.bra_subspace().size(); ++source_bra_index)
    for (int source_ket_index=0; source_ket_index<source_sector.ket_subspace().size(); ++source_ket_index)
          {
            // look up target matrix entry indices
            int remapped_bra_index, bra_relative_phase;
            std::tie(remapped_bra_index, bra_relative_phase)
              = two_body_mapping.state_mapping[source_sector.bra_subspace_index()][source_bra_index];
            if (remapped_bra_index == basis::kNone)
              continue;
            int remapped_ket_index, ket_relative_phase;
            std::tie(remapped_ket_index, ket_relative_phase)
              = two_body_mapping.state_mapping[source_sector.ket_subspace_index()][source_ket_index];
            if (remapped_ket_index == basis::kNone)
              continue;

            // copy entry
            target_matrix(remapped_bra_index,remapped_ket_index)
              = (bra_relative_phase*ket_relative_phase) * source_matrix(source_bra_index,source_ket_index);
          }

  return target_matrix;
}


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
