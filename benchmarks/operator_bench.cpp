#include "qudrip.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace qudrip;

namespace {

using clock_type = std::chrono::steady_clock;

struct Timing {
  double minimum;
  double median;
};

template <typename Function>
Timing time_trials(int trials, Function &&function) {
  std::vector<double> samples;
  samples.reserve(trials);
  for (int trial = 0; trial < trials; ++trial) {
    const auto begin = clock_type::now();
    function();
    const auto end = clock_type::now();
    samples.push_back(
        std::chrono::duration<double>(end - begin).count());
  }
  std::sort(samples.begin(), samples.end());
  return {samples.front(), samples[samples.size() / 2]};
}

value_type matrix_checksum(const SpMatrix &matrix) {
  value_type checksum = 0.0;
  for (Eigen::Index outer = 0; outer < matrix.outerSize(); ++outer) {
    for (SpMatrix::InnerIterator entry(matrix, outer); entry; ++entry) {
      const value_type weight(
          1.0 + static_cast<double>(entry.row()),
          0.5 + static_cast<double>(entry.col()));
      checksum += weight * entry.value();
    }
  }
  return checksum;
}

value_type vector_checksum(const Matrix &vector) {
  value_type checksum = 0.0;
  for (Eigen::Index row = 0; row < vector.rows(); ++row) {
    const value_type weight(1.0 + static_cast<double>(row),
                            0.25 * static_cast<double>(row + 1));
    checksum += weight * vector(row, 0);
  }
  return checksum;
}

template <typename RAW, typename SECTOR>
void run_spin_case(const std::string &label, RAW &raw, SECTOR &sector,
                   int sites, int trials) {
  auto hop = getHopBoseGate(raw);
  auto hamiltonian = getEmptyOperator();
  for (int site = 0; site + 1 < sites; ++site) {
    hamiltonian =
        hamiltonian + hop(site + 1, site) + hop(site, site + 1);
  }

  auto source = getState(sector, 1);
  auto destination = getState(sector, 1);
  for (idxv_type i = 0; i < sector.range(); ++i) {
    source.data()(i, 0) =
        value_type(std::sin(0.17 * double(i + 1)),
                   std::cos(0.11 * double(i + 1)));
  }
  source.data().col(0).normalize();

  SpMatrix matrix = hamiltonian >> sector;
  hamiltonian.map(destination, 0, source, 0);

  const Timing assembly = time_trials(trials, [&] {
    matrix = hamiltonian >> sector;
  });
  const Timing lazy_apply = time_trials(trials, [&] {
    hamiltonian.map(destination, 0, source, 0);
  });

  const Matrix expected = matrix * source.data().col(0);
  const double map_error =
      (destination.data().col(0) - expected).norm();
  const value_type assembly_checksum = matrix_checksum(matrix);
  const value_type map_checksum =
      vector_checksum(destination.data().col(0));

  std::cout << std::setprecision(12)
            << "operator_bench case=" << label
            << " operation=assembly"
            << " dimension=" << sector.range()
            << " terms=" << 2 * (sites - 1)
            << " nonzeros=" << matrix.nonZeros()
            << " trials=" << trials
            << " seconds_min=" << assembly.minimum
            << " seconds_median=" << assembly.median
            << " checksum_real=" << assembly_checksum.real()
            << " checksum_imag=" << assembly_checksum.imag() << '\n';

  std::cout << std::setprecision(12)
            << "operator_bench case=" << label
            << " operation=map"
            << " dimension=" << sector.range()
            << " terms=" << 2 * (sites - 1)
            << " trials=" << trials
            << " seconds_min=" << lazy_apply.minimum
            << " seconds_median=" << lazy_apply.median
            << " checksum_real=" << map_checksum.real()
            << " checksum_imag=" << map_checksum.imag()
            << " matrix_error=" << map_error << '\n';
}

void run_dense_bit_case(int trials) {
  auto index = getQbits(10);
  auto gate0 = getBoseGate(index);
  auto gate1 = getBoseGate(index);
  Matrix u0(2, 2);
  Matrix u1(2, 2);
  u0 << value_type(0.5, 0.25), value_type(-0.75, 0.5),
      value_type(1.25, -0.25), value_type(0.375, 0.125);
  u1 << value_type(-0.25, 0.75), value_type(0.5, -0.125),
      value_type(0.875, 0.25), value_type(-1.0, 0.375);
  gate0 << u0;
  gate1 << u1;
  auto op = gate0(2) * gate1(7) * gate0(4) * gate1(1);

  SpMatrix matrix = op >> index;
  const Timing assembly = time_trials(trials, [&] {
    matrix = op >> index;
  });
  const auto checksum = matrix_checksum(matrix);
  std::cout << std::setprecision(12)
            << "operator_bench case=dense_bit"
            << " operation=assembly"
            << " dimension=" << index.range()
            << " factors=4"
            << " nonzeros=" << matrix.nonZeros()
            << " trials=" << trials
            << " seconds_min=" << assembly.minimum
            << " seconds_median=" << assembly.median
            << " checksum_real=" << checksum.real()
            << " checksum_imag=" << checksum.imag() << '\n';
}

void run_dense_mode_case(int trials) {
  auto mode = getSingleMode(8);
  auto gate0 = getModeOp(mode);
  auto gate1 = getModeOp(mode);
  Matrix u0(8, 8);
  Matrix u1(8, 8);
  for (Eigen::Index row = 0; row < 8; ++row) {
    for (Eigen::Index column = 0; column < 8; ++column) {
      u0(row, column) =
          value_type(0.03 * double(row + 1) - 0.02 * double(column),
                     0.01 * double(row - column));
      u1(row, column) =
          value_type(-0.015 * double(row) +
                         0.025 * double(column + 1),
                     0.005 * double(row + column));
    }
  }
  gate0 << u0;
  gate1 << u1;
  auto op = gate0(0) * gate1(0) * gate0(0);

  SpMatrix matrix = op >> mode;
  const Timing assembly = time_trials(trials, [&] {
    matrix = op >> mode;
  });
  const auto checksum = matrix_checksum(matrix);
  std::cout << std::setprecision(12)
            << "operator_bench case=dense_mode"
            << " operation=assembly"
            << " dimension=" << mode.range()
            << " factors=3"
            << " nonzeros=" << matrix.nonZeros()
            << " trials=" << trials
            << " seconds_min=" << assembly.minimum
            << " seconds_median=" << assembly.median
            << " checksum_real=" << checksum.real()
            << " checksum_imag=" << checksum.imag() << '\n';
}

}  // namespace

int main(int argc, char **argv) {
  const int sites = argc > 1 ? std::atoi(argv[1]) : 16;
  const int particles = argc > 2 ? std::atoi(argv[2]) : sites / 2;
  const int trials = argc > 3 ? std::atoi(argv[3]) : 3;
  if (sites < 2 || sites >= BIT_LIMIT || particles < 0 ||
      particles > sites || trials <= 0) {
    std::cerr
        << "usage: operator_bench [spin-sites>=2] [particles] [trials]\n";
    return 2;
  }

  auto hash_raw = getQbits(sites);
  auto hash_number = getNset(hash_raw);
  auto hash_layout = getSpinClusters(hash_raw);
  auto hash_sector =
      Constrain(hash_layout, hash_number = particles);

  auto direct_raw = getQbits(sites);
  auto direct_number = getNset(direct_raw);
  auto direct_layout = getSpinClusters(direct_raw);
  direct_layout.useDirectRank();
  auto direct_sector =
      Constrain(direct_layout, direct_number = particles);

  if (hash_sector.range() == 0 ||
      hash_sector.range() != direct_sector.range())
    throw std::runtime_error("invalid benchmark sectors");

  run_spin_case("spin_hash", hash_raw, hash_sector, sites, trials);
  run_spin_case("spin_direct", direct_raw, direct_sector, sites, trials);
  run_dense_bit_case(trials);
  run_dense_mode_case(trials);
}
