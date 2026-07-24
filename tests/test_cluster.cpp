#include "qudrip.hpp"

#include <climits>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

using namespace qudrip;

namespace {

int cluster_checks = 0;

#define CLUSTER_CHECK(cond)                                               \
  do {                                                                    \
    ++cluster_checks;                                                     \
    if (!(cond)) {                                                        \
      std::cerr << "CLUSTER FAILED (line " << __LINE__ << "): " #cond    \
                << "\n";                                                  \
      std::exit(1);                                                       \
    }                                                                     \
  } while (0)

template <typename SUB, typename RAW>
std::vector<idxv_type> raw_mapping(SUB &sub, RAW &raw) {
  std::vector<idxv_type> result;
  result.reserve(sub.range());
  for (idxv_type i = 0; i < sub.range(); ++i) {
    sub[i];
    result.push_back(idxv_type(raw));
    CLUSTER_CHECK(idxv_type(sub) == i);
  }
  return result;
}

Matrix permutation_from_raw(const std::vector<idxv_type> &legacy,
                            const std::vector<idxv_type> &cluster) {
  Matrix permutation = Matrix::Zero(legacy.size(), cluster.size());
  for (size_t compact = 0; compact < cluster.size(); ++compact) {
    const auto it =
        std::find(legacy.begin(), legacy.end(), cluster[compact]);
    CLUSTER_CHECK(it != legacy.end());
    permutation(it - legacy.begin(), compact) = 1.0;
  }
  return permutation;
}

}  // namespace

int run_cluster_tests() {
  // One-bit clusters preserve the existing raw integer order.
  {
    auto raw = getQbits(4);
    auto number = getNset(raw);
    auto mutable_layout = getSpinClusters(raw);
    mutable_layout.useDirectRank();
    const auto layout = mutable_layout;
    auto sector = Constrain(layout, number = 2);
    const std::vector<idxv_type> expected{3, 5, 6, 9, 10, 12};

    CLUSTER_CHECK(sector.has_direct_rank());
    CLUSTER_CHECK(raw_mapping(sector, raw) == expected);
    auto &stripped = strip(sector);
    CLUSTER_CHECK(&stripped == &raw);

    for (idxv_type full = 0; full < raw.range(); ++full) {
      raw[full];
      const auto found = std::find(expected.begin(), expected.end(), full);
      const idxv_type expected_rank =
          found == expected.end() ? sector.range() : found - expected.begin();
      CLUSTER_CHECK(idxv_type(sector) == expected_rank);
    }
  }

  // The charge dimension is part of the backend type. Exercise more
  // components than the old small benchmark cases and verify every raw rank.
  {
    auto raw = getQbits(6);
    Nset<> bit0(raw, 0, 0);
    Nset<> bit1(raw, 1, 1);
    Nset<> bit2(raw, 2, 2);
    Nset<> bit3(raw, 3, 3);
    auto layout = getSpinClusters(raw);
    layout.useDirectRank();
    auto sector = Constrain(layout, bit0 = 1, bit1 = 0, bit2 = 1,
                            bit3 = 0);
    const std::vector<idxv_type> expected{5, 21, 37, 53};

    CLUSTER_CHECK(sector.has_direct_rank());
    CLUSTER_CHECK(raw_mapping(sector, raw) == expected);
    for (idxv_type full = 0; full < raw.range(); ++full) {
      raw[full];
      const auto found = std::find(expected.begin(), expected.end(), full);
      const idxv_type expected_rank =
          found == expected.end() ? sector.range() : found - expected.begin();
      CLUSTER_CHECK(idxv_type(sector) == expected_rank);
    }
    raw[raw.range()];
    CLUSTER_CHECK(idxv_type(sector) == sector.range());
  }

  // A Hubbard site is a four-state primitive. Cluster order differs from the
  // legacy spin-major raw order, while membership stays identical.
  {
    constexpr int sites = 3;

    auto legacy_raw = getQbits(2 * sites);
    auto legacy_up = getNup(legacy_raw, sites);
    auto legacy_down = getNdo(legacy_raw, sites);
    auto legacy =
        Constrain(legacy_raw, legacy_up = 1, legacy_down = 1);

    auto cluster_raw = getQbits(2 * sites);
    auto cluster_up = getNup(cluster_raw, sites);
    auto cluster_down = getNdo(cluster_raw, sites);
    auto layout = getSiteClusters(cluster_raw, sites);
    layout.useDirectRank();
    auto clustered =
        Constrain(layout, cluster_up = 1, cluster_down = 1);

    const std::vector<idxv_type> expected_legacy{
        9, 10, 12, 17, 18, 20, 33, 34, 36};
    const std::vector<idxv_type> expected_cluster{
        9, 10, 17, 18, 12, 20, 33, 34, 36};

    const auto legacy_map = raw_mapping(legacy, legacy_raw);
    const auto cluster_map = raw_mapping(clustered, cluster_raw);
    CLUSTER_CHECK(legacy_map == expected_legacy);
    CLUSTER_CHECK(cluster_map == expected_cluster);
    auto sorted_cluster = cluster_map;
    std::sort(sorted_cluster.begin(), sorted_cluster.end());
    CLUSTER_CHECK(sorted_cluster == legacy_map);
    CLUSTER_CHECK(clustered.has_direct_rank());
    for (idxv_type full = 0; full < cluster_raw.range(); ++full) {
      cluster_raw[full];
      const auto found =
          std::find(cluster_map.begin(), cluster_map.end(), full);
      const idxv_type expected_rank =
          found == cluster_map.end() ? clustered.range()
                                     : found - cluster_map.begin();
      CLUSTER_CHECK(idxv_type(clustered) == expected_rank);
    }

    auto reverse_raw = getQbits(2 * sites);
    auto reverse_up = getNup(reverse_raw, sites);
    auto reverse_down = getNdo(reverse_raw, sites);
    auto reverse_layout = getSiteClusters(reverse_raw, sites);
    reverse_layout.useDirectRank();
    auto reversed =
        Constrain(reverse_layout, reverse_down = 1, reverse_up = 1);
    CLUSTER_CHECK(raw_mapping(reversed, reverse_raw) == expected_cluster);

    auto legacy_hop = getHopFermiGate(legacy_raw);
    auto cluster_hop = getHopFermiGate(cluster_raw);
    auto legacy_number_gate = getBoseGate(legacy_raw);
    auto cluster_number_gate = getBoseGate(cluster_raw);
    Matrix occupation = Matrix::Zero(2, 2);
    occupation(1, 1) = 1.0;
    legacy_number_gate << occupation;
    cluster_number_gate << occupation;
    auto legacy_hamiltonian =
        legacy_hop(1, 0) + legacy_hop(0, 1) +
        legacy_hop(sites + 1, sites) + legacy_hop(sites, sites + 1) +
        value_type(0.37) * legacy_number_gate(0) +
        value_type(0.19) *
            (legacy_hop(1, 0) * legacy_hop(sites + 1, sites) +
             legacy_hop(0, 1) * legacy_hop(sites, sites + 1));
    auto cluster_hamiltonian =
        cluster_hop(1, 0) + cluster_hop(0, 1) +
        cluster_hop(sites + 1, sites) + cluster_hop(sites, sites + 1) +
        value_type(0.37) * cluster_number_gate(0) +
        value_type(0.19) *
            (cluster_hop(1, 0) * cluster_hop(sites + 1, sites) +
             cluster_hop(0, 1) * cluster_hop(sites, sites + 1));

    const Matrix old_matrix = Matrix(legacy_hamiltonian >> legacy);
    const Matrix cluster_matrix =
        Matrix(cluster_hamiltonian >> clustered);
    const Matrix permutation =
        permutation_from_raw(legacy_map, cluster_map);
    CLUSTER_CHECK(cluster_matrix.isApprox(
        permutation.adjoint() * old_matrix * permutation, 1e-12));

    SpMatrix full_cluster_matrix = cluster_hamiltonian >> cluster_raw;
    CLUSTER_CHECK(Matrix(restrict(cluster_raw, clustered,
                                  full_cluster_matrix))
                      .isApprox(cluster_matrix, 1e-12));

    auto source = getState(clustered, 1);
    auto destination = getState(clustered, 1);
    for (idxv_type i = 0; i < clustered.range(); ++i)
      source.data()(i, 0) =
          value_type(double(i + 1), -0.25 * double(i));
    destination(0) = cluster_hamiltonian * source(0);
    CLUSTER_CHECK(destination.data().col(0).isApprox(
        cluster_matrix * source.data().col(0), 1e-12));

    auto raw_docc = getQbits(2 * sites);
    auto nup = getNup(raw_docc, sites);
    auto ndo = getNdo(raw_docc, sites);
    auto nd = getNd(raw_docc, sites);
    auto docc_layout = getSiteClusters(raw_docc, sites);
    docc_layout.useDirectRank();
    auto docc =
        Constrain(docc_layout, nup = 1, ndo = 1, nd = 1);
    CLUSTER_CHECK(raw_mapping(docc, raw_docc) ==
                  std::vector<idxv_type>({9, 18, 36}));

    auto raw_holes = getQbits(2 * sites);
    auto holes_up = getNup(raw_holes, sites);
    auto holes_down = getNdo(raw_holes, sites);
    auto nh = getNh(raw_holes, sites);
    auto holes_layout = getSiteClusters(raw_holes, sites);
    holes_layout.useDirectRank();
    auto holes =
        Constrain(holes_layout, holes_up = 1, holes_down = 1, nh = 1);
    CLUSTER_CHECK(raw_mapping(holes, raw_holes) ==
                  std::vector<idxv_type>({10, 17, 12, 20, 33, 34}));

    // CQContainer is intentionally opaque to the first compiler. Its fallback
    // still scans in site-cluster order rather than sorting raw labels.
    auto opaque_raw = getQbits(2 * sites);
    auto opaque_up = getNup(opaque_raw, sites);
    auto opaque_down = getNdo(opaque_raw, sites);
    auto opaque_total = opaque_up + opaque_down;
    auto opaque_layout = getSiteClusters(opaque_raw, sites);
    auto opaque = Constrain(opaque_layout, opaque_total = 2);
    const std::vector<idxv_type> opaque_expected{
        9, 3, 10, 17, 24, 18, 5, 12, 6, 20, 33, 40, 34, 48, 36};
    CLUSTER_CHECK(!opaque.has_direct_rank());
    CLUSTER_CHECK(raw_mapping(opaque, opaque_raw) == opaque_expected);
  }

  // A DP-budget fallback still enumerates in cluster order and uses the
  // materialized inverse map.
  {
    auto raw = getQbits(6);
    auto up = getNup(raw, 3);
    auto down = getNdo(raw, 3);
    auto layout = getSiteClusters(raw, 3);
    const auto old_budget = clusterDPMemoryBudget();
    setClusterDPMemoryBudget(0);
    auto fallback = Constrain(layout, up = 1, down = 1);
    setClusterDPMemoryBudget(old_budget);

    const std::vector<idxv_type> expected{
        9, 10, 17, 18, 12, 20, 33, 34, 36};
    CLUSTER_CHECK(!fallback.has_direct_rank());
    CLUSTER_CHECK(raw_mapping(fallback, raw) == expected);
  }

  // A safely representable impossible target is rejected from cluster maxima
  // without allocating a DP table or scanning the raw space.
  {
    auto raw = getQbits(2);
    auto number = getNset(raw);
    auto layout = getSpinClusters(raw);
    layout.useDirectRank();
    auto empty = Constrain(layout, number = 3);
    CLUSTER_CHECK(empty.range() == 0);
    CLUSTER_CHECK(empty.has_direct_rank());
    for (idxv_type full = 0; full < raw.range(); ++full) {
      raw[full];
      CLUSTER_CHECK(idxv_type(empty) == empty.range());
    }
  }

  // Targets outside the exact, int-representable DP domain retain the
  // historical Evaluate() round-trip semantics through the fallback.
  {
    auto raw = getQbits(2);
    auto number = getNset(raw);
    auto layout = getSpinClusters(raw);
    const idxv_type wrapped_target = idxv_type(1) << 32;
    auto wrapped = Constrain(layout, number = wrapped_target);
    CLUSTER_CHECK(!wrapped.has_direct_rank());
    CLUSTER_CHECK(raw_mapping(wrapped, raw) ==
                  std::vector<idxv_type>({0}));
  }

  // Tensor-product bosons use their existing mixed-radix order.
  {
    auto mode0 = getSingleMode(2);
    auto mode1 = getSingleMode(3);
    auto mode2 = getSingleMode(4);
    auto modes = mode0 * mode1 * mode2;
    auto layout = getModeClusters(modes);
    layout.useDirectRank();
    auto number = getTotalOccupation(layout);
    auto sector = Constrain(layout, number = 3);
    const std::vector<idxv_type> expected{3, 6, 9, 14, 17, 20};

    CLUSTER_CHECK(sector.has_direct_rank());
    CLUSTER_CHECK(raw_mapping(sector, modes) == expected);
    for (idxv_type full = 0; full < modes.range(); ++full) {
      modes[full];
      const auto found = std::find(expected.begin(), expected.end(), full);
      const idxv_type expected_rank =
          found == expected.end() ? sector.range() : found - expected.begin();
      CLUSTER_CHECK(idxv_type(sector) == expected_rank);
    }

    auto copied = sector;
    copied[2];
    ++copied;
    CLUSTER_CHECK(idxv_type(copied) == 3);
    auto moved = std::move(copied);
    ++moved;
    CLUSTER_CHECK(idxv_type(moved) == 4);

    auto legacy_mode0 = getSingleMode(2);
    auto legacy_mode1 = getSingleMode(3);
    auto legacy_mode2 = getSingleMode(4);
    auto legacy_modes = legacy_mode0 * legacy_mode1 * legacy_mode2;
    auto legacy_layout = getModeClusters(legacy_modes);
    auto legacy_number = getTotalOccupation(legacy_layout);
    auto legacy =
        Constrain(legacy_modes, legacy_number = 3);
    CLUSTER_CHECK(raw_mapping(legacy, legacy_modes) == expected);

    auto create0 = getModeOp(mode0);
    auto destroy1 = getModeOp(mode1);
    create0 << LadderMatrix(mode0.range(), Mode::Upper);
    destroy1 << LadderMatrix(mode1.range(), Mode::Lower);
    auto cluster_hop = create0(0) * destroy1(0);

    auto old_create0 = getModeOp(legacy_mode0);
    auto old_destroy1 = getModeOp(legacy_mode1);
    old_create0 << LadderMatrix(legacy_mode0.range(), Mode::Upper);
    old_destroy1 << LadderMatrix(legacy_mode1.range(), Mode::Lower);
    auto legacy_hop = old_create0(0) * old_destroy1(0);

    CLUSTER_CHECK(Matrix(cluster_hop >> sector)
                      .isApprox(Matrix(legacy_hop >> legacy), 1e-12));
  }

  // The same mode adapter also handles the one-cluster tensor-product base
  // case.
  {
    auto mode = getSingleMode(5);
    auto layout = getModeClusters(mode);
    auto number = getTotalOccupation(layout);
    auto sector = Constrain(layout, number = 3);
    CLUSTER_CHECK(raw_mapping(sector, mode) ==
                  std::vector<idxv_type>({3}));
  }

  // A bosonic quantity is bound to both the raw product and its mode
  // dimensions. Reinterpreting the same raw range with different radices must
  // select the exhaustive predicate rather than a mismatched DP.
  {
    auto mode0 = getSingleMode(2);
    auto mode1 = getSingleMode(6);
    auto modes = mode0 * mode1;
    auto original_layout = getModeClusters(modes);
    auto number = getTotalOccupation(original_layout);
    ClusterLayout<decltype(modes)> alternate_layout(
        modes, ClusterLayoutKind::BosonModes, {3, 4});
    auto sector = Constrain(alternate_layout, number = 3);
    CLUSTER_CHECK(!sector.has_direct_rank());
    CLUSTER_CHECK(raw_mapping(sector, modes) ==
                  std::vector<idxv_type>({3, 8}));
  }

  // The built-in bosonic adapter rejects tensor products containing a
  // non-mode leaf.
  {
    auto spin = getQbits(1);
    auto mode = getSingleMode(3);
    auto mixed = spin * mode;
    bool mixed_rejected = false;
    try {
      auto invalid_layout = getModeClusters(mixed);
      (void)invalid_layout;
    } catch (const std::invalid_argument &) {
      mixed_rejected = true;
    }
    CLUSTER_CHECK(mixed_rejected);
  }

  // Public layout construction enforces the encoding-specific local
  // dimensions, so encode/decode and DP rank cannot disagree.
  {
    auto spin_raw = getQbits(2);
    bool spin_rejected = false;
    try {
      ClusterLayout<QbitsIndex<>> malformed(
          spin_raw, ClusterLayoutKind::SpinBits, {4});
      (void)malformed;
    } catch (const std::invalid_argument &) {
      spin_rejected = true;
    }
    CLUSTER_CHECK(spin_rejected);

    auto site_raw = getQbits(4);
    bool sites_rejected = false;
    try {
      ClusterLayout<QbitsIndex<>> malformed(
          site_raw, ClusterLayoutKind::HubbardSites, {4, 4}, 1);
      (void)malformed;
    } catch (const std::invalid_argument &) {
      sites_rejected = true;
    }
    CLUSTER_CHECK(sites_rejected);
  }

  // Characterize the existing unsigned-difference behavior in Evaluate().
  {
    auto raw = getQbits(2);
    auto number = getNset(raw);
    auto below = LooseConstrain<1, true>(raw, number = 2);
    CLUSTER_CHECK(raw_mapping(below, raw) ==
                  std::vector<idxv_type>({1, 2, 3}));

    auto above = LooseConstrain<1, true>(raw, number = 0);
    CLUSTER_CHECK(raw_mapping(above, raw) ==
                  std::vector<idxv_type>({0, 1, 2}));
  }

  // Tensor-product range multiplication rejects overflow before a layout or
  // DP sees a wrapped range.
  {
    auto large0 = getSingleMode(INT_MAX);
    auto large1 = getSingleMode(INT_MAX);
    auto large2 = getSingleMode(INT_MAX);
    auto product01 = large0 * large1;
    bool overflow_rejected = false;
    try {
      auto overflow = product01 * large2;
      (void)overflow;
    } catch (const std::overflow_error &) {
      overflow_rejected = true;
    }
    CLUSTER_CHECK(overflow_rejected);
  }

  return cluster_checks;
}
