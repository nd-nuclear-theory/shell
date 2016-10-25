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
    domain_covered = true;
    range_covered = true;

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
    state_mapping.resize(source_space.size());

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
        // std::cout << fmt::format("Subspace {} -> {}",source_subspace_index,target_subspace_index) << std::endl;
        if (target_subspace_index==basis::kNone)
            continue;
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
              domain_covered = false;
            else
              ++num_target_states_found;

            // std::cout << fmt::format(
            //     "  subspace {} source {} target {}",
            //     source_subspace_index,source_state_index,target_state_index
            //   )
            //           << std::endl;

          }

        // check if all target states covered for this subspace
        range_covered &= (num_target_states_found==target_subspace.size());
      }
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
        // skip unmapped source subspace
        if (state_mapping[source_subspace_index].size()==0)
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
        "Flags: domain_covered {} range_covered {}",
        domain_covered,range_covered
      )
              <<std::endl;

    return os.str();

  }


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
