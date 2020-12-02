/****************************************************************
  tbme_mapping.h

  Defines mapping of JJJPN-scheme two-body matrix state indexings and
  remapping copy function for two-body matrices.

  Mark A. Caprio
  University of Notre Dame

  + 10/21/16 (mac): Created, as two_body_mapping.
  + 11/4/16 (mac):
    - Add storage of subspace_mapping.
    - Add RemappedMatrixJJJPN.
    - Rename to tbme_mapping.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.
  + 11/25/20 (pjf): Add OpenMP parallelization to RemappedMatrixJJJPN.

****************************************************************/

#ifndef TBME_MAPPING_H_
#define TBME_MAPPING_H_

#include <cstddef>
#include <vector>
#include <limits>

#include "eigen3/Eigen/Dense"

#include "basis/nlj_orbital.h"
#include "basis/jjjpn_scheme.h"

namespace shell {

  struct TwoBodyMapping
  // Container for mapping two-body state indices from one set of
  // JJJPN-scheme subspaces to another:
  //
  //   (source subspace index, source state index) -> (target state index, relative phase)
  //
  // We do not keep track of the target subspace index, since this is
  // generally already known in situations where we are using this
  // mapping.  If no target subspace exists with labels matching the
  // source subspace, the state indexing vector for that source
  // subspace is left empty.
  //
  // This may also be used as a mask on the source two-body states, to
  // select only source two-body states which fall into some "target"
  // truncation, by simply considering the boolean condition that the
  // source state map to *some* state in the target indexing.
  {

    // constructors

    TwoBodyMapping() = default;
    // default constructor

    TwoBodyMapping(
        const basis::OrbitalSpacePN& source_orbital_space,
        const basis::TwoBodySpaceJJJPN& source_space,
        const basis::OrbitalSpacePN& target_orbital_space,
        const basis::TwoBodySpaceJJJPN& target_space
      );
    // Construct mapping based on given source and target orbitals and
    // two-body spaces.

    // mapping data

    static constexpr std::size_t kNone = basis::kNone;
    // Flag value for missing target.

    std::vector<std::vector<std::size_t>> orbital_mapping;
    // Nested vectors providing mapping of orbital indices from
    // (source_subspace_index, state_index) to target_state_index.
    //
    // Precondition: It is assumed that both p and n subspaces will be
    // present.

    std::vector<std::size_t> subspace_mapping;
    // Vectors providing mapping of two-body subspace indices from
    // source_subspace_index to target_subspace_index.
    //
    // Note: In practice this information may not be useful, since
    // commonly we will be working in a context in which we have
    // already been given a target sector (bra and ket subspaces) to
    // construct, and will have looked up the source sector (bra and
    // ket subspaces) which needs to be mapped to it, so we do not
    // need to look up the indices of the subspaces.  Nonetheless, the
    // subspace mapping can be informative diagnostic information.

    std::vector<std::vector<std::tuple<std::size_t,std::size_t>>> state_mapping;
    // Nested vectors providing mapping of two-body state indices from
    // (source_subspace_index, state_index) to
    // (target_state_index, relative_phase).

    // informational flags for mapping properties

    bool domain_subspaces_covered, domain_states_covered,
      range_subspaces_covered, range_states_covered;
    // Whether or not all source and target two-body states are used.

    // debugging

    std::string DebugStr() const;

  };

Eigen::MatrixXd
RemappedMatrixJJJPN(
      const basis::TwoBodySectorsJJJPN::SectorType& source_sector,
      const basis::TwoBodySectorsJJJPN::SectorType& target_sector,
      const shell::TwoBodyMapping& two_body_mapping,
      const Eigen::MatrixXd& source_matrix
  );
// Copy matrix elements from source block to target block, with
// remapping of state indexing.


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
