#include "qudrip.hpp"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string_view>

using namespace qudrip;

namespace {

template <typename Build>
auto measure(std::string_view label, Build &&build) {
  const auto start = std::chrono::steady_clock::now();
  auto sector = build();
  const auto stop = std::chrono::steady_clock::now();
  const std::chrono::duration<double> elapsed = stop - start;
  std::cout << label << ": dimension=" << sector.range()
            << ", seconds=" << elapsed.count()
            << ", direct_rank=" << sector.has_direct_rank() << '\n';
  return sector;
}

template <typename Build>
void measure_matrix(std::string_view label, Build &&build) {
  const auto start = std::chrono::steady_clock::now();
  auto matrix = build();
  const auto stop = std::chrono::steady_clock::now();
  const std::chrono::duration<double> elapsed = stop - start;
  std::cout << label << ": nonzeros=" << matrix.nonZeros()
            << ", seconds=" << elapsed.count() << '\n';
}

}  // namespace

int main(int argc, char **argv) {
  const int spin_sites = argc > 1 ? std::atoi(argv[1]) : 20;
  const int hubbard_sites = argc > 2 ? std::atoi(argv[2]) : 10;
  const int spin_particles =
      argc > 3 ? std::atoi(argv[3]) : spin_sites / 2;
  const int hubbard_particles =
      argc > 4 ? std::atoi(argv[4]) : hubbard_sites / 2;
  const bool direct_rank = argc > 5 && std::atoi(argv[5]) != 0;
  if (spin_sites <= 0 || hubbard_sites <= 0 || spin_particles < 0 ||
      spin_particles > spin_sites || hubbard_particles < 0 ||
      hubbard_particles > hubbard_sites) {
    std::cerr << "usage: constraint_bench [spin-sites] [Hubbard-sites] "
                 "[spin-particles] [particles-per-spin] [direct-rank]\n";
    return 2;
  }

  auto legacy_spin = getQbits(spin_sites);
  auto legacy_spin_number = getNset(legacy_spin);
  auto legacy_spin_sector = measure("spin legacy", [&] {
    return Constrain(legacy_spin, legacy_spin_number = spin_particles);
  });

  auto clustered_spin = getQbits(spin_sites);
  auto clustered_spin_number = getNset(clustered_spin);
  auto spin_layout = getSpinClusters(clustered_spin);
  spin_layout.useDirectRank(direct_rank);
  auto clustered_spin_sector = measure("spin cluster", [&] {
    return Constrain(spin_layout, clustered_spin_number = spin_particles);
  });

  auto legacy_spin_hop = getHopBoseGate(legacy_spin);
  auto clustered_spin_hop = getHopBoseGate(clustered_spin);
  measure_matrix("spin legacy hop", [&] {
    return legacy_spin_hop(1, 0) >> legacy_spin_sector;
  });
  measure_matrix("spin cluster hop", [&] {
    return clustered_spin_hop(1, 0) >> clustered_spin_sector;
  });
  auto legacy_spin_sum = getEmptyOperator();
  auto clustered_spin_sum = getEmptyOperator();
  for (int site = 0; site + 1 < spin_sites; ++site) {
    legacy_spin_sum =
        legacy_spin_sum + legacy_spin_hop(site + 1, site) +
        legacy_spin_hop(site, site + 1);
    clustered_spin_sum =
        clustered_spin_sum + clustered_spin_hop(site + 1, site) +
        clustered_spin_hop(site, site + 1);
  }
  measure_matrix("spin legacy sum", [&] {
    return legacy_spin_sum >> legacy_spin_sector;
  });
  measure_matrix("spin cluster sum", [&] {
    return clustered_spin_sum >> clustered_spin_sector;
  });

  auto legacy_hubbard = getQbits(2 * hubbard_sites);
  auto legacy_up = getNup(legacy_hubbard, hubbard_sites);
  auto legacy_down = getNdo(legacy_hubbard, hubbard_sites);
  auto legacy_hubbard_sector = measure("Hubbard legacy", [&] {
    return Constrain(legacy_hubbard, legacy_up = hubbard_particles,
                     legacy_down = hubbard_particles);
  });

  auto clustered_hubbard = getQbits(2 * hubbard_sites);
  auto clustered_up = getNup(clustered_hubbard, hubbard_sites);
  auto clustered_down = getNdo(clustered_hubbard, hubbard_sites);
  auto site_layout =
      getSiteClusters(clustered_hubbard, hubbard_sites);
  site_layout.useDirectRank(direct_rank);
  auto clustered_hubbard_sector = measure("Hubbard cluster", [&] {
    return Constrain(site_layout, clustered_up = hubbard_particles,
                     clustered_down = hubbard_particles);
  });

  auto legacy_hubbard_hop = getHopFermiGate(legacy_hubbard);
  auto clustered_hubbard_hop = getHopFermiGate(clustered_hubbard);
  measure_matrix("Hubbard legacy hop", [&] {
    return legacy_hubbard_hop(1, 0) >> legacy_hubbard_sector;
  });
  measure_matrix("Hubbard cluster hop", [&] {
    return clustered_hubbard_hop(1, 0) >> clustered_hubbard_sector;
  });
  auto legacy_hubbard_sum = getEmptyOperator();
  auto clustered_hubbard_sum = getEmptyOperator();
  for (int site = 0; site + 1 < hubbard_sites; ++site) {
    for (const int offset : {0, hubbard_sites}) {
      legacy_hubbard_sum =
          legacy_hubbard_sum +
          legacy_hubbard_hop(offset + site + 1, offset + site) +
          legacy_hubbard_hop(offset + site, offset + site + 1);
      clustered_hubbard_sum =
          clustered_hubbard_sum +
          clustered_hubbard_hop(offset + site + 1, offset + site) +
          clustered_hubbard_hop(offset + site, offset + site + 1);
    }
  }
  measure_matrix("Hubbard legacy sum", [&] {
    return legacy_hubbard_sum >> legacy_hubbard_sector;
  });
  measure_matrix("Hubbard cluster sum", [&] {
    return clustered_hubbard_sum >> clustered_hubbard_sector;
  });
}
