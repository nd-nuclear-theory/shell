/****************************************************************
  h2_io.h

  Patrick J. Fasano
  University of Notre Dame

  + 08/31/12 (mac): Adapted from mfdn_io.h (as mfdn_h2).
****************************************************************/

#ifndef PARTITIONUTILS_GROUP_READ_H_
#define PARTITIONUTILS_GROUP_READ_H_

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <array>
#include <execution>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>

#include "basis/nlj_orbital.h"
#include "mcutils/fortran_io.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"
#include "parallel_hashmap/phmap.h"

// namespace
// {

#ifndef MAX_A
#  define MAX_A 20
#endif

#define SYMMETRIC_PARTITIONS


using replacements_map_t = std::map<uint16_t, uint16_t>;
using groupid_t = std::array<uint16_t, MAX_A>;
using value_t = uint32_t;
namespace std
{
template<> struct hash<groupid_t>
{
  std::size_t operator()(groupid_t const& g) const
  {
    return phmap::Hash<decltype(std::tuple_cat(g))>()(std::tuple_cat(g));
  }
};
}  // namespace std

using hash_map_t = phmap::parallel_flat_hash_map<
    groupid_t,
    value_t,
    phmap::priv::hash_default_hash<groupid_t>,
    phmap::priv::hash_default_eq<groupid_t>,
    phmap::priv::Allocator<std::pair<const groupid_t, value_t>>,  // alias for std::allocator
    12,
    std::mutex
  >;

struct MBGroupsMetadata
{
  int32_t version_number, num_classes, num_particles, reserved3;
  int32_t parity, twoM, max_Nmax, num_diag;
  int32_t num_groupids, num_states, reserved10, num_blocks;
  int32_t reserved12, reserved13, reserved14, reserved15;
};

struct MFDnPartitioning
{
  std::vector<int> proton_partitions, neutron_partitions;
};

struct MFDnSMWFInfo
{
  int version;
  int Z, N, twoM;
  std::string uniqueID;
  basis::OrbitalSpacePN orbital_space;
  int num_proton_states, num_neutron_states;
  MFDnPartitioning partitioning;
  int parity;
  float weight_max;
  std::size_t dimension;
  int num_diag;
};

std::size_t total_dimension(const hash_map_t& groupid_map)
{
  return std::reduce(
      std::execution::par,
      std::cbegin(groupid_map),
      std::cend(groupid_map),
      std::size_t{0},
      [](const auto& lhs, const auto& rhs) -> std::size_t { return std::size_t(lhs) + std::size_t(rhs.second); }
    );
}

auto group_map_stats(const hash_map_t& groupid_map)
{
  const auto& [largest_group, largest_dim] = *std::max_element(
      std::execution::par,
      groupid_map.begin(),
      groupid_map.end(),
      [](const auto& lhs, const auto& rhs) { return (lhs.second < rhs.second); }
    );

  return make_tuple(largest_group, largest_dim);
}

MFDnSMWFInfo ReadMFDnSMWFInfo(std::string filename)
{
  MFDnSMWFInfo info;
  std::string line;
  int line_count = 0;
  mcutils::FileExistCheck(
      filename, /*exit_on_nonexist=*/true, /*warn_on_overwrite=*/false
    );
  auto stream = std::ifstream(filename, std::ios_base::in);

  // line 1: version
  {
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> info.version;
    mcutils::ParsingCheck(line_stream, line_count, line);
    assert(info.version == 15200);
  }

  // line 2: Z, N, 2Mj
  {
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> info.Z >> info.N >> info.twoM;
    mcutils::ParsingCheck(line_stream, line_count, line);
    assert((info.Z + info.N) % 2 == info.twoM % 2);
  }

  // line 3: unique identifier
  {
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> info.uniqueID;
    mcutils::ParsingCheck(line_stream, line_count, line);
  }

  // line 4+: number of proton, neutron orbitals, orbital listing body
  {
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    int num_orbitals_p, num_orbitals_n;
    line_stream >> num_orbitals_p >> num_orbitals_n;
    mcutils::ParsingCheck(line_stream, line_count, line);

    std::size_t num_orbitals = num_orbitals_p + num_orbitals_n;
    std::string orbital_info_str;
    for (std::size_t orbital_line_count = 0; orbital_line_count < num_orbitals;
         ++orbital_line_count)
    {
      mcutils::GetLine(stream, line, line_count);
      orbital_info_str.append(line);
      orbital_info_str.append("\n");  // need to restore newline to input line
    }
    std::istringstream orbital_info_stream(orbital_info_str);
    auto orbital_list = basis::ParseOrbitalPNStream(
        orbital_info_stream,
        /*standalone=*/false,
        basis::MFDnOrbitalFormat::kVersion15200
      );
    info.orbital_space = basis::OrbitalSpacePN(orbital_list);
  }

  // derive number of proton and neutron states
  info.num_proton_states = 0;
  info.num_neutron_states = 0;
  const auto& proton_subspace = info.orbital_space.GetSubspace(0);
  for (int index = 0; index < proton_subspace.size(); ++index)
  {
    info.num_proton_states += TwiceValue(proton_subspace.GetState(index).j()) + 1;
  }
  const auto& neutron_subspace = info.orbital_space.GetSubspace(0);
  for (int index = 0; index < neutron_subspace.size(); ++index)
  {
    info.num_neutron_states +=
        TwiceValue(neutron_subspace.GetState(index).j()) + 1;
  }


  // partitioning
  {
    int num_proton_partitions, num_neutron_partitions;
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> num_proton_partitions >> num_neutron_partitions;
    mcutils::ParsingCheck(line_stream, line_count, line);

    auto& partitioning = info.partitioning;
    partitioning.proton_partitions.reserve(num_proton_partitions);
    for (int i = 0; i < num_proton_partitions; ++i)
    {
      int partition;
      stream >> partition;
      partitioning.proton_partitions.push_back(partition);
    }
    partitioning.neutron_partitions.reserve(num_neutron_partitions);
    for (int i = 0; i < num_neutron_partitions; ++i)
    {
      int partition;
      stream >> partition;
      partitioning.neutron_partitions.push_back(info.num_proton_states + partition);
    }
  }

  // basis information
  {
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> info.parity >> info.weight_max >> info.dimension
        >> info.num_diag;
    mcutils::ParsingCheck(line_stream, line_count, line);
  }

  return info;
}

std::vector<uint16_t> GenerateSPBinPartitionTable(const MFDnSMWFInfo& info)
{
  std::vector<uint16_t> sp_bin_partition_table{};
  sp_bin_partition_table.reserve(info.num_proton_states + info.num_neutron_states);

  {
    auto partition_it = info.partitioning.proton_partitions.cbegin();
    auto partition_end = info.partitioning.proton_partitions.cend();
    for (int i = 1; i <= info.num_proton_states; ++i)
    {
      if ((std::next(partition_it) != partition_end) && (i >= *std::next(partition_it)))
        partition_it = std::next(partition_it);
      sp_bin_partition_table.push_back(*partition_it);
    }
  }

  {
    auto partition_it = info.partitioning.neutron_partitions.cbegin();
    auto partition_end = info.partitioning.neutron_partitions.cend();
    for (int i = info.num_proton_states + 1;
         i <= info.num_proton_states + info.num_neutron_states;
         ++i)
    {
      if ((std::next(partition_it) != partition_end) && (i >= *std::next(partition_it)))
        partition_it = std::next(partition_it);
      sp_bin_partition_table.push_back(*partition_it);
    }
  }

  return sp_bin_partition_table;
}

MFDnPartitioning GenerateAutomaticPartitioning(const basis::OrbitalSpacePN& orbital_space)
{
  MFDnPartitioning partitioning{};
  const auto& proton_subspace = orbital_space.GetSubspace(0);
  int partition = 1;
  for (int index = 0; index < proton_subspace.size(); ++index)
  {
    partitioning.proton_partitions.push_back(partition);
    partition += TwiceValue(proton_subspace.GetState(index).j()) + 1;
  }
  const auto& neutron_subspace = orbital_space.GetSubspace(0);
  for (int index = 0; index < neutron_subspace.size(); ++index)
  {
    partitioning.neutron_partitions.push_back(partition);
    partition += TwiceValue(proton_subspace.GetState(index).j()) + 1;
  }

  return partitioning;
}

void ReadMBGroups(
    std::string filename_pattern,
    const MFDnSMWFInfo& smwf_info,
    hash_map_t& groupid_map,
    bool verbose = false
  )
{
  // convenience mode variable
  static constexpr std::ios_base::openmode mode_argument =
      std::ios_base::in | std::ios_base::binary;

  // get sp_bin
  const auto spbin_partition_table = GenerateSPBinPartitionTable(smwf_info);

  // read metadata
  {
    // open stream
    std::string filename = fmt::format(filename_pattern, 1);
    mcutils::FileExistCheck(
        filename, /*exit_on_nonexist=*/true, /*warn_on_overwrite=*/false
      );
    auto stream = std::ifstream(filename, mode_argument);
    mcutils::StreamCheck(
        bool(stream), filename, "Failure opening groups file for input"
      );

    const auto metadata_vector =
        mcutils::ReadFortranRecord<MBGroupsMetadata>(stream);
    const auto metadata = metadata_vector[0];
    assert(metadata.num_particles <= MAX_A);
    assert(metadata.num_particles == smwf_info.Z + smwf_info.N);
    assert(metadata.twoM == smwf_info.twoM);
    groupid_map.reserve(metadata.num_groupids * metadata.num_diag);
  }


  // read group files and build hash map
  double avg_time;
  for (std::size_t i = 1; i <= smwf_info.num_diag; ++i)
  {
    mcutils::SteadyTimer timer{};
    const std::string filename = fmt::format(filename_pattern, i);
    if (verbose)
      fmt::print("reading file {} ({}/{})\n", filename, i, smwf_info.num_diag);
    auto stream = std::ifstream(filename, mode_argument);
    mcutils::StreamCheck(
        bool(stream), filename, "Failure opening groups file for input"
      );

    if (verbose)
      timer.Start();
    const auto metadata_vector =
        mcutils::ReadFortranRecord<MBGroupsMetadata>(stream);
    const auto& metadata = metadata_vector[0];
    if (verbose)
      fmt::print("  num_groups: {}\n", metadata.num_groupids);
    mcutils::SkipFortranRecord(stream);  // nblksNm
    const auto Mstateptr = mcutils::ReadFortranRecord<int32_t>(stream);
    const auto groupidlist = mcutils::ReadFortranRecord<int16_t>(stream);

    assert(metadata.num_groupids * metadata.num_particles == groupidlist.size());
    assert(Mstateptr.size() == metadata.num_groupids + 1);

#pragma omp parallel for
    for (std::size_t j = 0; j < metadata.num_groupids; ++j)
    {
      groupid_t groupid{};
      // build groupid, translating to partition numbers
      //#pragma omp simd
      for (int p = 0; p < metadata.num_particles; ++p)
        groupid[p] =
            spbin_partition_table[groupidlist[metadata.num_particles * j + p] - 1];
      assert(Mstateptr[j + 1] - Mstateptr[j] > 0);
      groupid_map.lazy_emplace_l(
          groupid,
          [](hash_map_t::value_type&) {},
          [&](
              const hash_map_t::constructor& ctor
            ) { ctor(groupid, Mstateptr[j + 1] - Mstateptr[j]); }
        );
    }
    if (verbose)
    {
      timer.Stop();
      fmt::print(
          "  size:{}  load_factor: {}\n",
          groupid_map.size(),
          groupid_map.load_factor()
        );
      fmt::print(
          "elapsed time: {:>8.2f}  average time: {:>8.2f}\n\n",
          timer.ElapsedTime(),
          avg_time = (avg_time * (i - 1) + timer.ElapsedTime()) / i
        );
      fflush(stdout);
    }
  }

  std::size_t dimension = total_dimension(groupid_map);

  if (verbose)
    fmt::print("dimension: {:d}\n\n", dimension);

  assert(smwf_info.dimension == dimension);
}

template<
    typename T,
    std::enable_if_t<std::is_same_v<typename T::key_type, typename T::mapped_type>>* = nullptr
  >
typename T::mapped_type recurse_map(const T& map, const typename T::key_type& key)
{
  typename T::mapped_type value = map.at(key);
  if (map.count(value))
    return recurse_map(map, value);
  return value;
}

std::map<uint16_t, uint16_t> GetPartitionReplacements(
    const MFDnSMWFInfo& smwf_info, bool verbose = false
  )
{
  auto auto_partitioning = GenerateAutomaticPartitioning(smwf_info.orbital_space);
  if (verbose)
  {
    fmt::print(
        "automatic partitions_p: {:>4d}\n",
        fmt::join(auto_partitioning.proton_partitions, " ")
      );
    fmt::print(
        "automatic partitions_n: {:>4d}\n",
        fmt::join(auto_partitioning.neutron_partitions, " ")
      );
  }

  std::map<uint16_t, uint16_t> partition_replacements{};
#ifdef SYMMETRIC_PARTITIONS
  const auto& auto_partitions = auto_partitioning.proton_partitions;
  const auto& partitions = smwf_info.partitioning.proton_partitions;
#else
  for (auto&& [auto_partitions, partitions] :
       {std::make_pair(
            std::ref(auto_partitioning.proton_partitions),
            std::ref(smwf_info.partitioning.proton_partitions)
          ),
        std::make_pair(
            std::ref(auto_partitioning.neutron_partitions),
            std::ref(smwf_info.partitioning.neutron_partitions)
          )})
#endif
    for (int i = 0; i < auto_partitions.size(); ++i)
    {
      auto start = std::upper_bound(
          partitions.begin(), partitions.end(), auto_partitions[i]
        );
      auto end = partitions.end();
      if ((i + 1) < auto_partitions.size())
        end = std::lower_bound(
            partitions.begin(), partitions.end(), auto_partitions[i + 1]
          );
      if (std::distance(start, end) > 0)
      {
        for (auto it = --end; it != start; --it)
        {
          partition_replacements[*it] = *std::prev(it);
        }
        partition_replacements[*start] = auto_partitions[i];
      }
    }

  if (verbose)
    fmt::print(
        "partition replacements: {}\n", fmt::join(partition_replacements, " ")
      );

  return partition_replacements;
}


template<typename T>
typename T::mapped_type atdefault(
    const T& map,
    const typename T::key_type& key,
    const typename T::mapped_type& default_val = {}
  )
{
  if (map.count(key))
    return map.at(key);
  return default_val;
}

auto ReplacePartitions(
    const hash_map_t& groupid_map,
    const replacements_map_t& replacements,
    bool verbose = false
  )
{
  hash_map_t new_groupid_map{};
  if (replacements.size() > 0)
  {
    replacements_map_t transformed_replacements{};
    std::transform(
        replacements.begin(),
        replacements.end(),
        std::inserter(transformed_replacements, transformed_replacements.end()),
        [&](auto v) { return std::pair{v.first, recurse_map(replacements, v.first)}; }
      );
    if (verbose)
      fmt::print("replacements: {}\n", fmt::join(transformed_replacements, " "));

    auto replace =
        [&transformed_replacements](auto i) { return atdefault(transformed_replacements, i, i); };
    auto replaced_groupid = [&replace](groupid_t g) { std::transform(g.begin(), g.end(), g.begin(), replace); return g; };

    new_groupid_map.reserve(groupid_map.size());
    groupid_map.for_each(
        std::execution::par,
        [&](const auto& el) {
            groupid_t new_g = replaced_groupid(el.first);
            new_groupid_map.try_emplace_l(
                new_g,
                [&](auto& v){v.second += el.second;},
                el.second
              );
          }
      );
  }
  else
  {
    new_groupid_map = groupid_map;
  }

  return new_groupid_map;
}

static double rpm_avg_time{};
static int rpm_count{};
std::optional<std::pair<replacements_map_t, uint32_t>> ReplacePartitionsMain(
    const MFDnSMWFInfo& smwf_info,
    const hash_map_t& groupid_map,
    const replacements_map_t& replacements,
    const replacements_map_t& active_replacements = {}
  )
{
  mcutils::SteadyTimer timer{};
  timer.Start();
  const auto& new_groupid_map =
      ReplacePartitions(groupid_map, active_replacements, false);
  const auto new_dimension  = total_dimension(new_groupid_map);
  assert(new_dimension == smwf_info.dimension);
  const auto& [largest_group, largest_dim] = group_map_stats(new_groupid_map);
  timer.Stop();
  ++rpm_count;

  fmt::print(
      "replacements: {}\n"
      "  num groupids: {:d}\n"
      "  largest group: {}\n"
      "  largest group dimension: {:d}\n"
      "elapsed time: {:>8.2f}  average time: {:>8.2f}\n",
      active_replacements,
      groupid_map.size(),
      largest_group,
      largest_dim,
      timer.ElapsedTime(),
      rpm_avg_time = (rpm_avg_time * (rpm_count - 1) + timer.ElapsedTime()) / rpm_count
    );
  fflush(stdout);

  if (largest_dim > std::numeric_limits<int16_t>::max())
  {
    fmt::print(
        "largest group is larger than std::numeric_limits<int16_t>::max()\n\n"
      );
    fflush(stdout);
    return std::nullopt;
  }
  else
  {
    fmt::print("\n");
  }

  replacements_map_t best_replacements = active_replacements;
  uint32_t min_number_groups = groupid_map.size();
  for (auto&& it = replacements.begin(); it != replacements.end(); ++it)
  {
    replacements_map_t new_replacements = active_replacements;
    new_replacements.insert(*it);
#ifdef SYMMETRIC_PARTITIONS
    new_replacements.insert({
        it->first+smwf_info.num_proton_states,
        it->second+smwf_info.num_proton_states
      });
#endif
    auto stats = ReplacePartitionsMain(
        smwf_info, new_groupid_map, replacements_map_t{next(it), replacements.end()}, new_replacements
      );
    if (stats && (stats->second < min_number_groups))
    {
      std::tie(best_replacements, min_number_groups) = *stats;
    }
  }

  return {{best_replacements, min_number_groups}};
}

MFDnPartitioning GenerateSimplifiedPartitioning(
    MFDnPartitioning partitioning, const replacements_map_t& replacements
  )
{
  for (std::vector<int>& partitions :
       {std::ref(partitioning.proton_partitions),
        std::ref(partitioning.neutron_partitions)})
  {
    std::for_each(
        partitions.begin(), partitions.end(),
        [&](auto& p){ if (replacements.count(p)) p = recurse_map(replacements, p); }
      );
    auto last = std::unique(partitions.begin(), partitions.end());
    partitions.erase(last, partitions.end());
  }

  return partitioning;
}

int main(int argc, char* argv[])
{
  hash_map_t groupid_map{};
  MBGroupsMetadata metadata{};
  const auto smwf_info = ReadMFDnSMWFInfo("mfdn_smwf.info");
  fmt::print(
      "partitions_p: {:>4d}\n",
      fmt::join(smwf_info.partitioning.proton_partitions, " ")
    );
  fmt::print(
      "partitions_n: {:>4d}\n",
      fmt::join(smwf_info.partitioning.neutron_partitions, " ")
    );
  auto replacements = GetPartitionReplacements(smwf_info);
  fmt::print("dimension: {:d}\n", smwf_info.dimension);
  fflush(stdout);
  ReadMBGroups("mfdn_MBgroups{:03d}", smwf_info, groupid_map, false);
  fmt::print("number of groups: {:d}\n", groupid_map.size());
  fflush(stdout);

  const auto& stats = ReplacePartitionsMain(smwf_info, groupid_map, replacements);

  if (stats)
  {
    const auto& [best_replacements, min_groups] = *stats;
    fmt::print(
        "\n\n"
        "best replacements: {}\n"
        "minimum number of groups: {:d}\n"
        "initial number of groups: {:d}\n"
        "reduction factor: {:4.1f}x",
        best_replacements,
        min_groups,
        groupid_map.size(),
        double(groupid_map.size())/double(min_groups)
      );

    const auto& simplified_partitioning =
        GenerateSimplifiedPartitioning(smwf_info.partitioning, best_replacements);
    fmt::print("\nmfdn_partitioning.info:\n{:->64s}\n", "");
    fmt::print(
        "  {:>4d} {:>4d}\n",
        simplified_partitioning.proton_partitions.size(),
        simplified_partitioning.neutron_partitions.size()
      );
    fmt::print(
        "  {:>4d}\n", fmt::join(simplified_partitioning.proton_partitions, " ")
      );
    std::vector<int> shifted_neutron_partitions{};
    std::transform(
        simplified_partitioning.neutron_partitions.begin(),
        simplified_partitioning.neutron_partitions.end(),
        std::back_inserter(shifted_neutron_partitions),
        [&](
            auto& p
          ) { return p - simplified_partitioning.neutron_partitions.at(0) + 1; }
      );
    fmt::print("  {:>4d}\n", fmt::join(shifted_neutron_partitions, " "));
  }
  else
  {
    fmt::print("\n\nNo valid partitionings found!\n");
  }
  return 0;
}

// }  // namespace
#endif  // PARTITIONUTILS_GROUP_READ_H_
