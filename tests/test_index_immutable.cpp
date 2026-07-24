#include "qudrip.hpp"

#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace qudrip;

namespace {

int checks = 0;

#define INDEX_CHECK(condition)                                                 \
  do {                                                                         \
    ++checks;                                                                  \
    if (!(condition)) {                                                        \
      std::cerr << "FAILED (line " << __LINE__ << "): " #condition << '\n';    \
      std::exit(1);                                                            \
    }                                                                          \
  } while (false)

template <typename Basis> void check_identity_basis(const Basis &basis) {
  for (idxv_type compact = 0; compact < basis.range(); ++compact) {
    INDEX_CHECK(basis.raw_at(compact) == compact);
    INDEX_CHECK(basis.rank(compact) == compact);
  }
  INDEX_CHECK(basis.raw_at(basis.range()) == basis.range());
  INDEX_CHECK(basis.raw_at(std::numeric_limits<idxv_type>::max()) ==
              basis.range());
  INDEX_CHECK(basis.rank(basis.range()) == basis.range());
  INDEX_CHECK(basis.rank(std::numeric_limits<idxv_type>::max()) ==
              basis.range());
}

template <typename Basis> void check_constrained_basis(const Basis &basis) {
  for (idxv_type compact = 0; compact < basis.range(); ++compact)
    INDEX_CHECK(basis.rank(basis.raw_at(compact)) == compact);
  INDEX_CHECK(basis.raw_at(basis.range()) == basis.range());
  INDEX_CHECK(basis.raw_at(std::numeric_limits<idxv_type>::max()) ==
              basis.range());
}

} // namespace

int main() {
  {
    bool negative_qbits_rejected = false;
    bool oversized_qbits_rejected = false;
    bool empty_mode_rejected = false;
    try {
      (void)getQbits(-1);
    } catch (const std::invalid_argument &) {
      negative_qbits_rejected = true;
    }
    try {
      (void)getQbits(std::numeric_limits<idx_size_t>::digits);
    } catch (const std::overflow_error &) {
      oversized_qbits_rejected = true;
    }
    try {
      (void)getSingleMode(0);
    } catch (const std::invalid_argument &) {
      empty_mode_rejected = true;
    }
    INDEX_CHECK(negative_qbits_rejected);
    INDEX_CHECK(oversized_qbits_rejected);
    INDEX_CHECK(empty_mode_rejected);
  }

  {
    auto bits = getQbits(4);
    bits[9];
    check_identity_basis(bits);
    INDEX_CHECK(idxv_type(bits) == 9);

    const auto bindings = detail::collect_raw_leaf_bindings(bits);
    INDEX_CHECK(bindings.size() == 1);
    INDEX_CHECK(bindings[0].identity == std::addressof(bits));
    INDEX_CHECK(bindings[0].kind == detail::RawLeafKind::Qbits);
    INDEX_CHECK(bindings[0].range == 16);
    INDEX_CHECK(bindings[0].stride == 1);
    INDEX_CHECK(bindings[0].extract(11) == 11);
    INDEX_CHECK(bindings[0].replace(11, 11, 3) == 3);
    INDEX_CHECK(idxv_type(bits) == 9);
  }

  {
    auto mode = getSingleMode(5);
    mode[3];
    check_identity_basis(mode);
    INDEX_CHECK(idxv_type(mode) == 3);

    const auto bindings = detail::collect_raw_leaf_bindings(mode);
    INDEX_CHECK(bindings.size() == 1);
    INDEX_CHECK(bindings[0].identity == std::addressof(mode));
    INDEX_CHECK(bindings[0].kind == detail::RawLeafKind::SingleMode);
    INDEX_CHECK(bindings[0].range == 5);
    INDEX_CHECK(bindings[0].stride == 1);
  }

  {
    auto mode2 = getSingleMode(2);
    auto mode3 = getSingleMode(3);
    auto mode4 = getSingleMode(4);
    auto product = mode2 * mode3 * mode4;
    product[17];

    check_identity_basis(product);
    INDEX_CHECK(idxv_type(product) == 17);
    INDEX_CHECK(idxv_type(mode2) == 1);
    INDEX_CHECK(idxv_type(mode3) == 1);
    INDEX_CHECK(idxv_type(mode4) == 1);

    const auto bindings = detail::collect_raw_leaf_bindings(product);
    INDEX_CHECK(bindings.size() == 3);
    const std::vector<const void *> identities{
        std::addressof(mode2), std::addressof(mode3), std::addressof(mode4)};
    const std::vector<idxv_type> ranges{2, 3, 4};
    const std::vector<idxv_type> strides{12, 4, 1};
    for (size_t leaf = 0; leaf < bindings.size(); ++leaf) {
      INDEX_CHECK(bindings[leaf].identity == identities[leaf]);
      INDEX_CHECK(bindings[leaf].range == ranges[leaf]);
      INDEX_CHECK(bindings[leaf].stride == strides[leaf]);
    }

    for (idxv_type raw = 0; raw < product.range(); ++raw) {
      for (size_t leaf = 0; leaf < bindings.size(); ++leaf) {
        const auto old_local = bindings[leaf].extract(raw);
        INDEX_CHECK(old_local == (raw / strides[leaf]) % ranges[leaf]);
        for (idxv_type replacement = 0; replacement < ranges[leaf];
             ++replacement) {
          const auto changed =
              bindings[leaf].replace(raw, old_local, replacement);
          INDEX_CHECK(bindings[leaf].extract(changed) == replacement);
          for (size_t other = 0; other < bindings.size(); ++other) {
            if (other == leaf)
              continue;
            INDEX_CHECK(bindings[other].extract(changed) ==
                        bindings[other].extract(raw));
          }
        }
      }
    }

    INDEX_CHECK(idxv_type(product) == 17);
    INDEX_CHECK(idxv_type(mode2) == 1);
    INDEX_CHECK(idxv_type(mode3) == 1);
    INDEX_CHECK(idxv_type(mode4) == 1);

    auto tail = mode3 * mode4;
    auto right_associated = mode2 * tail;
    const auto right_bindings =
        detail::collect_raw_leaf_bindings(right_associated);
    INDEX_CHECK(right_bindings.size() == bindings.size());
    for (size_t leaf = 0; leaf < bindings.size(); ++leaf) {
      INDEX_CHECK(right_bindings[leaf].identity == bindings[leaf].identity);
      INDEX_CHECK(right_bindings[leaf].range == bindings[leaf].range);
      INDEX_CHECK(right_bindings[leaf].stride == bindings[leaf].stride);
    }
  }

  {
    auto bits = getQbits(2);
    auto mode = getSingleMode(3);
    auto bits_then_mode = bits * mode;
    auto mode_then_bits = mode * bits;

    const auto left = detail::collect_raw_leaf_bindings(bits_then_mode);
    INDEX_CHECK(left.size() == 2);
    INDEX_CHECK(left[0].identity == std::addressof(bits));
    INDEX_CHECK(left[0].stride == 3);
    INDEX_CHECK(left[1].identity == std::addressof(mode));
    INDEX_CHECK(left[1].stride == 1);

    const auto right = detail::collect_raw_leaf_bindings(mode_then_bits);
    INDEX_CHECK(right.size() == 2);
    INDEX_CHECK(right[0].identity == std::addressof(mode));
    INDEX_CHECK(right[0].stride == 4);
    INDEX_CHECK(right[1].identity == std::addressof(bits));
    INDEX_CHECK(right[1].stride == 1);
  }

  {
    auto mode = getSingleMode(3);
    auto ambiguous = mode * mode;
    const auto bindings = detail::collect_raw_leaf_bindings(ambiguous);
    INDEX_CHECK(bindings.size() == 2);
    INDEX_CHECK(detail::find_unique_raw_leaf_binding(
                    bindings, std::addressof(mode)) == nullptr);
    int absent = 0;
    INDEX_CHECK(detail::find_unique_raw_leaf_binding(
                    bindings, std::addressof(absent)) == nullptr);
  }

  {
    auto raw = getQbits(3);
    using raw_reference = decltype(strip(raw));
    SubIndex<raw_reference, true> sector(strip(raw),
                                         std::vector<idxv_type>{1, 2, 4});

    check_constrained_basis(sector);
    INDEX_CHECK(sector.raw_at(0) == 1);
    INDEX_CHECK(sector.raw_at(1) == 2);
    INDEX_CHECK(sector.raw_at(2) == 4);
    INDEX_CHECK(sector.rank(3) == sector.range());

    sector[1];
    INDEX_CHECK(idxv_type(sector) == 1);
    (void)sector.raw_at(2);
    (void)sector.rank(4);
    INDEX_CHECK(idxv_type(sector) == 1);
    INDEX_CHECK(idxv_type(raw) == 2);

    raw[6];
    INDEX_CHECK(idxv_type(raw) == 6);

    const auto bindings = detail::collect_raw_leaf_bindings(sector);
    INDEX_CHECK(bindings.size() == 1);
    INDEX_CHECK(bindings[0].identity == std::addressof(raw));
    INDEX_CHECK(bindings[0].range == raw.range());
    INDEX_CHECK(bindings[0].stride == 1);
    INDEX_CHECK(idxv_type(raw) == 6);

    SubIndex<raw_reference, true> empty(strip(raw), std::vector<idxv_type>{});
    INDEX_CHECK(empty.range() == 0);
    INDEX_CHECK(empty.raw_at(0) == empty.range());
    INDEX_CHECK(empty.rank(0) == empty.range());
    const auto empty_bindings = detail::collect_raw_leaf_bindings(empty);
    INDEX_CHECK(empty_bindings.size() == 1);
    INDEX_CHECK(empty_bindings[0].identity == std::addressof(raw));
  }

  {
    auto hash_raw = getQbits(4);
    auto hash_number = getNset(hash_raw);
    auto hash_layout = getSpinClusters(hash_raw);
    auto hash_sector = Constrain(hash_layout, hash_number = 2);
    check_constrained_basis(hash_sector);
    hash_raw[7];
    INDEX_CHECK(hash_sector.rank(7) == hash_sector.range());
    INDEX_CHECK(idxv_type(hash_raw) == 7);

    auto direct_raw = getQbits(4);
    auto direct_number = getNset(direct_raw);
    auto direct_layout = getSpinClusters(direct_raw);
    direct_layout.useDirectRank();
    auto direct_sector = Constrain(direct_layout, direct_number = 2);
    check_constrained_basis(direct_sector);
    direct_raw[7];
    INDEX_CHECK(direct_sector.rank(7) == direct_sector.range());
    INDEX_CHECK(idxv_type(direct_raw) == 7);

    const auto direct_bindings =
        detail::collect_raw_leaf_bindings(direct_sector);
    INDEX_CHECK(direct_bindings.size() == 1);
    INDEX_CHECK(direct_bindings[0].identity == std::addressof(direct_raw));
    INDEX_CHECK(idxv_type(direct_raw) == 7);
  }

  {
    detail::RawLeafBinding binding{nullptr, 3, 4};
    INDEX_CHECK(binding.extract(11) == 2);
    INDEX_CHECK(binding.replace(11, 2, 0) == 3);
    INDEX_CHECK(binding.replace(3, 0, 2) == 11);

    bool rejected = false;
    try {
      (void)binding.replace(3, 1, 0);
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    INDEX_CHECK(rejected);

    rejected = false;
    try {
      (void)binding.replace(3, 0, 3);
    } catch (const std::out_of_range &) {
      rejected = true;
    }
    INDEX_CHECK(rejected);

    detail::RawLeafBinding invalid{nullptr, 2, 0};
    rejected = false;
    try {
      (void)invalid.extract(0);
    } catch (const std::logic_error &) {
      rejected = true;
    }
    INDEX_CHECK(rejected);
  }

  std::cout << checks << " immutable-index checks passed\n";
}
