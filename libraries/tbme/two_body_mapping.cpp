/****************************************************************
  two_body_mapping.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "tbme/two_body_mapping.h"

#include <sstream>

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
            //
            // "Missing" target orbitals will be indexed by basis::kNone.
            const basis::TwoBodyStateJJJPN source_state(source_subspace,source_state_index);
            const typename basis::TwoBodyStateJJJPN::StateLabelsType target_state_labels(
                orbital_mapping[int(source_state.GetOrbital1().orbital_species())][source_state.index1()],
                orbital_mapping[int(source_state.GetOrbital2().orbital_species())][source_state.index2()]
              );

            // look up corresponding target state index
            int target_state_index = target_subspace.LookUpStateIndex(
                source_subspace.GetStateLabels(source_state_index)
              );
            state_mapping[source_subspace_index][source_state_index] = target_state_index;

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
            int target_state_index = state_mapping[source_subspace_index][source_state_index];
            os << fmt::format(
                "  subspace {:3} source {:3} target {:3}",
                source_subspace_index,source_state_index,target_state_index
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


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
