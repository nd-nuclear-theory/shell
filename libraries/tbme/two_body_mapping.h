/****************************************************************
  two_body_mapping.h                       

  Defines mapping and masking of JJJPN-scheme two-body matrix element
  indexings.
                                  
  Mark A. Caprio
  University of Notre Dame

  10/21/16 (mac): Created.
  11/4/16 (mac): Add storage of subspace_mapping.
     
****************************************************************/

#ifndef TWO_BODY_MAPPING_H_
#define TWO_BODY_MAPPING_H_

#include <vector>

#include "basis/nlj_orbital.h"
#include "basis/jjjpn_scheme.h"

namespace shell {

  struct TwoBodyMapping
  // Container for mapping two-body state indices from one set of
  // JJJPN-scheme subspaces to another:
  //
  //   (source subspace index, source state index) -> target state index
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

    static const int kNone = -1;
    // Flag value for missing target.

    std::vector<std::vector<int>> orbital_mapping;
    // Nested vectors providing mapping of orbital indices from
    // (source_subspace_index, state_index) to target_state_index.
    // 
    // Precondition: It is assumed that both p and n subspaces will be
    // present.

    std::vector<int> subspace_mapping;
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

    std::vector<std::vector<int>> state_mapping;
    // Nested vectors providing mapping of two-body state indices from
    // (source_subspace_index, state_index) to target_state_index.
    
    // informational flags for mapping properties

    bool domain_subspaces_covered, domain_states_covered,
      range_subspaces_covered, range_states_covered;
    // Whether or not all source and target two-body states are used.
    
    // debugging

    std::string DebugStr() const;

  };

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
