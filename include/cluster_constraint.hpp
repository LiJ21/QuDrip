#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace qudrip {

enum class ClusterLayoutKind { SpinBits, HubbardSites, BosonModes };

template <typename RawIndex>
  requires (!std::is_const_v<RawIndex>)
class ClusterLayout {
 public:
  struct is_cluster_layout {};
  using raw_index_type = RawIndex;

 private:
  RawIndex *raw_;
  ClusterLayoutKind kind_;
  std::vector<idxv_type> local_dims_;
  std::vector<idxv_type> local_strides_;
  size_t site_count_;
  bool direct_rank_ = false;

 public:
  ClusterLayout(RawIndex &raw, ClusterLayoutKind kind,
                std::vector<idxv_type> local_dims, size_t site_count = 0)
      : raw_(&raw),
        kind_(kind),
        local_dims_(std::move(local_dims)),
        local_strides_(local_dims_.size()),
        site_count_(site_count) {
    idxv_type product = 1;
    for (size_t cluster = local_dims_.size(); cluster-- > 0;) {
      const auto dim = local_dims_[cluster];
      if (dim == 0)
        throw std::invalid_argument("cluster local dimension must be positive");
      local_strides_[cluster] = product;
      product = detail::checked_index_product(product, dim);
    }
    switch (kind_) {
      case ClusterLayoutKind::SpinBits:
        if (std::any_of(local_dims_.begin(), local_dims_.end(),
                        [](idxv_type dim) { return dim != 2; }))
          throw std::invalid_argument(
              "spin clusters must all have dimension two");
        break;
      case ClusterLayoutKind::HubbardSites:
        if (site_count_ != local_dims_.size() ||
            site_count_ > BIT_LIMIT / 2 ||
            std::any_of(local_dims_.begin(), local_dims_.end(),
                        [](idxv_type dim) { return dim != 4; }))
          throw std::invalid_argument(
              "Hubbard clusters require one four-state cluster per site");
        break;
      case ClusterLayoutKind::BosonModes:
        break;
    }
    if (product != raw_->range())
      throw std::invalid_argument(
          "cluster dimensions do not match the raw index range");
  }

  RawIndex &raw_index() const noexcept { return *raw_; }
  ClusterLayoutKind kind() const noexcept { return kind_; }
  size_t cluster_count() const noexcept { return local_dims_.size(); }
  idxv_type local_state_count(size_t cluster) const {
    return local_dims_.at(cluster);
  }
  const std::vector<idxv_type> &local_dimensions() const noexcept {
    return local_dims_;
  }
  size_t site_count() const noexcept { return site_count_; }
  idxv_type range() const noexcept { return raw_->range(); }
  bool uses_direct_rank() const noexcept { return direct_rank_; }

  ClusterLayout &useDirectRank(bool enabled = true) noexcept {
    direct_rank_ = enabled;
    return *this;
  }

  idxv_type encode(std::span<const idxv_type> local_states) const {
    if (local_states.size() != local_dims_.size())
      throw std::invalid_argument("wrong number of cluster-local states");

    idxv_type raw = 0;
    switch (kind_) {
      case ClusterLayoutKind::SpinBits:
        for (size_t cluster = 0; cluster < local_states.size(); ++cluster) {
          if (local_states[cluster] >= 2)
            throw std::out_of_range("invalid spin-cluster state");
          raw = (raw << 1) | local_states[cluster];
        }
        break;

      case ClusterLayoutKind::HubbardSites:
        for (size_t cluster = 0; cluster < local_states.size(); ++cluster) {
          const auto local = local_states[cluster];
          if (local >= 4)
            throw std::out_of_range("invalid Hubbard-cluster state");
          const size_t site = site_count_ - cluster - 1;
          if (local & 1ULL) raw |= idxv_type(1) << site;
          if (local & 2ULL) raw |= idxv_type(1) << (site_count_ + site);
        }
        break;

      case ClusterLayoutKind::BosonModes:
        for (size_t cluster = 0; cluster < local_states.size(); ++cluster) {
          if (local_states[cluster] >= local_dims_[cluster])
            throw std::out_of_range("invalid bosonic cluster state");
          raw = detail::checked_index_product(raw, local_dims_[cluster]);
          if (raw > std::numeric_limits<idxv_type>::max() -
                        local_states[cluster])
            throw std::overflow_error("cluster encoding overflow");
          raw += local_states[cluster];
        }
        break;
    }
    return raw;
  }

  bool decode(idxv_type raw, std::span<idxv_type> local_states) const {
    if (raw >= range() || local_states.size() != local_dims_.size())
      return false;

    switch (kind_) {
      case ClusterLayoutKind::SpinBits:
        for (size_t cluster = 0; cluster < local_states.size(); ++cluster) {
          const size_t bit = local_states.size() - cluster - 1;
          local_states[cluster] = (raw >> bit) & 1ULL;
        }
        break;

      case ClusterLayoutKind::HubbardSites:
        for (size_t cluster = 0; cluster < local_states.size(); ++cluster) {
          const size_t site = site_count_ - cluster - 1;
          const idxv_type up = (raw >> site) & 1ULL;
          const idxv_type down = (raw >> (site_count_ + site)) & 1ULL;
          local_states[cluster] = up | (down << 1);
        }
        break;

      case ClusterLayoutKind::BosonModes:
        for (size_t cluster = local_states.size(); cluster-- > 0;) {
          local_states[cluster] = raw % local_dims_[cluster];
          raw /= local_dims_[cluster];
        }
        break;
    }
    return true;
  }

  idxv_type local_state(idxv_type raw, size_t cluster) const {
    if (raw >= range() || cluster >= local_dims_.size())
      throw std::out_of_range("invalid raw or cluster index");

    switch (kind_) {
      case ClusterLayoutKind::SpinBits: {
        const size_t bit = local_dims_.size() - cluster - 1;
        return (raw >> bit) & 1ULL;
      }
      case ClusterLayoutKind::HubbardSites: {
        const size_t site = site_count_ - cluster - 1;
        const idxv_type up = (raw >> site) & 1ULL;
        const idxv_type down = (raw >> (site_count_ + site)) & 1ULL;
        return up | (down << 1);
      }
      case ClusterLayoutKind::BosonModes:
        return (raw / local_strides_[cluster]) % local_dims_[cluster];
    }
    return 0;
  }

  idxv_type local_raw_fragment(size_t cluster, idxv_type local_state) const {
    if (cluster >= local_dims_.size() ||
        local_state >= local_dims_[cluster])
      throw std::out_of_range("invalid cluster-local state");

    switch (kind_) {
      case ClusterLayoutKind::SpinBits: {
        const size_t bit = local_dims_.size() - cluster - 1;
        return local_state << bit;
      }
      case ClusterLayoutKind::HubbardSites: {
        const size_t site = site_count_ - cluster - 1;
        idxv_type raw = 0;
        if (local_state & 1ULL) raw |= idxv_type(1) << site;
        if (local_state & 2ULL)
          raw |= idxv_type(1) << (site_count_ + site);
        return raw;
      }
      case ClusterLayoutKind::BosonModes:
        return local_state;
    }
    return 0;
  }

  idxv_type append_local_state(idxv_type raw_prefix, size_t cluster,
                               idxv_type local_state) const {
    if (cluster >= local_dims_.size() ||
        local_state >= local_dims_[cluster])
      throw std::out_of_range("invalid cluster-local state");

    switch (kind_) {
      case ClusterLayoutKind::SpinBits:
        return (raw_prefix << 1) | local_state;
      case ClusterLayoutKind::HubbardSites:
        return raw_prefix | local_raw_fragment(cluster, local_state);
      case ClusterLayoutKind::BosonModes: {
        raw_prefix =
            detail::checked_index_product(raw_prefix, local_dims_[cluster]);
        if (raw_prefix >
            std::numeric_limits<idxv_type>::max() - local_state)
          throw std::overflow_error("cluster encoding overflow");
        return raw_prefix + local_state;
      }
    }
    return raw_prefix;
  }
};

namespace detail {

inline size_t active_qbits(idxv_type range) {
  if (range == 0 || (range & (range - 1)) != 0)
    throw std::invalid_argument("QbitsIndex range must be a power of two");
  size_t bits = 0;
  while (range > 1) {
    range >>= 1;
    ++bits;
  }
  return bits;
}

}  // namespace detail

template <typename T>
auto getSpinClusters(QbitsIndex<T> &index) {
  const size_t bits = detail::active_qbits(index.range());
  return ClusterLayout<QbitsIndex<T>>(
      index, ClusterLayoutKind::SpinBits,
      std::vector<idxv_type>(bits, idxv_type(2)));
}

template <typename T>
auto getSiteClusters(QbitsIndex<T> &index, size_t site_count) {
  const size_t bits = detail::active_qbits(index.range());
  if (site_count > BIT_LIMIT / 2 || bits != 2 * site_count)
    throw std::invalid_argument(
        "Hubbard site layout requires exactly two bits per site");
  return ClusterLayout<QbitsIndex<T>>(
      index, ClusterLayoutKind::HubbardSites,
      std::vector<idxv_type>(site_count, idxv_type(4)), site_count);
}

template <typename IDX>
concept IndicesLike = requires {
  typename std::remove_reference_t<IDX>::is_indices;
};

template <typename T>
auto getModeClusters(SingleModeIndex<T> &index) {
  return ClusterLayout<SingleModeIndex<T>>(
      index, ClusterLayoutKind::BosonModes, {index.range()});
}

template <IndicesLike IDX>
  requires (!std::is_const_v<IDX>)
auto getModeClusters(IDX &index) {
  if (!index.all_single_mode_indices())
    throw std::invalid_argument(
        "bosonic mode layouts require only SingleModeIndex leaves");
  auto dimensions = index.local_ranges();
  return ClusterLayout<IDX>(index, ClusterLayoutKind::BosonModes,
                            std::move(dimensions));
}

template <typename RawIndex>
class TotalOccupation
    : public ConservedQuantityBase<TotalOccupation<RawIndex>, idx_size_t> {
  using base =
      ConservedQuantityBase<TotalOccupation<RawIndex>, idx_size_t>;
  friend class ConservedQuantityBase<TotalOccupation<RawIndex>, idx_size_t>;

  RawIndex *index_;
  std::vector<idxv_type> dimensions_;

  idx_size_t eval() {
    idxv_type raw = idxv_type(*index_);
    idx_size_t total = 0;
    for (size_t cluster = dimensions_.size(); cluster-- > 0;) {
      total += raw % dimensions_[cluster];
      raw /= dimensions_[cluster];
    }
    return total;
  }

 public:
  using quant_type = idx_size_t;
  using raw_index_type = RawIndex;
  using base::operator=;

  TotalOccupation(RawIndex &index, std::vector<idxv_type> dimensions)
      : index_(&index), dimensions_(std::move(dimensions)) {}

  bool bound_to(const RawIndex &candidate) const noexcept {
    return index_ == &candidate;
  }
  bool matches_dimensions(
      const std::vector<idxv_type> &candidate) const noexcept {
    return dimensions_ == candidate;
  }
};

template <typename RawIndex>
auto getTotalOccupation(ClusterLayout<RawIndex> &layout) {
  if (layout.kind() != ClusterLayoutKind::BosonModes)
    throw std::invalid_argument(
        "total occupation requires a bosonic mode layout");
  return TotalOccupation<RawIndex>(layout.raw_index(),
                                   layout.local_dimensions());
}

template <typename RawIndex>
auto getTotalOccupation(ClusterLayout<RawIndex> &&) = delete;

namespace detail {

template <typename T>
struct is_nset : std::false_type {};
template <typename T>
struct is_nset<Nset<T>> : std::true_type {};

template <typename T>
struct is_nd : std::false_type {};
template <typename T>
struct is_nd<Nd<T>> : std::true_type {};

template <typename T>
struct is_total_occupation : std::false_type {};
template <typename T>
struct is_total_occupation<TotalOccupation<T>> : std::true_type {};

inline size_t &cluster_dp_memory_budget_storage() {
  static size_t budget = size_t(256) * 1024 * 1024;
  return budget;
}

inline size_t checked_size_product(size_t lhs, size_t rhs) {
  if (lhs != 0 && rhs > std::numeric_limits<size_t>::max() / lhs)
    throw std::overflow_error("cluster DP table size overflow");
  return lhs * rhs;
}

inline bool try_size_product(size_t lhs, size_t rhs,
                             size_t &result) noexcept {
  if (lhs != 0 && rhs > std::numeric_limits<size_t>::max() / lhs)
    return false;
  result = lhs * rhs;
  return true;
}

inline bool try_size_add(size_t lhs, size_t rhs, size_t &result) noexcept {
  if (rhs > std::numeric_limits<size_t>::max() - lhs) return false;
  result = lhs + rhs;
  return true;
}

inline bool try_account_bytes(size_t count, size_t item_size,
                              size_t &total) noexcept {
  size_t bytes = 0;
  size_t next_total = 0;
  if (!try_size_product(count, item_size, bytes) ||
      !try_size_add(total, bytes, next_total))
    return false;
  total = next_total;
  return true;
}

inline idxv_type checked_count_add(idxv_type lhs, idxv_type rhs) {
  if (rhs > std::numeric_limits<idxv_type>::max() - lhs)
    throw std::overflow_error("constrained Hilbert-space dimension overflow");
  return lhs + rhs;
}

template <size_t M>
struct ClusterChargeProgram {
  static_assert(M > 0);
  using charge_type = std::array<idxv_type, M>;

  std::array<idxv_type, M> targets{};
  // Rows for cluster c are [cluster_offsets[c], cluster_offsets[c + 1]).
  std::vector<size_t> cluster_offsets;
  std::vector<charge_type> delta;

  std::array<size_t, M> radices{};
  std::array<size_t, M> weights{};
  std::vector<size_t> weighted_delta;
  std::vector<std::uint8_t> target_feasible;
  size_t charge_state_count = 0;
  size_t target_state = 0;

  size_t row(size_t cluster, idxv_type local) const noexcept {
    return cluster_offsets[cluster] + static_cast<size_t>(local);
  }
};

template <typename Layout, typename CQ>
bool supports_cluster_charge(const Layout &layout, CQ &cq) {
  using charge_type = std::remove_cvref_t<CQ>;
  using raw_type = typename std::remove_reference_t<Layout>::raw_index_type;

  if constexpr (is_nset<charge_type>::value) {
    if constexpr (!std::is_same_v<raw_type, QbitsIndex<idx_size_t>>) {
      return false;
    } else {
      if (!cq.bound_to(layout.raw_index()) ||
          (layout.kind() != ClusterLayoutKind::SpinBits &&
           layout.kind() != ClusterLayoutKind::HubbardSites))
        return false;
      return true;
    }
  } else if constexpr (is_nd<charge_type>::value) {
    if constexpr (!std::is_same_v<raw_type, QbitsIndex<idx_size_t>>) {
      return false;
    } else {
      if (!cq.bound_to(layout.raw_index()) ||
          layout.kind() != ClusterLayoutKind::HubbardSites ||
          cq.site_count() != static_cast<int>(layout.site_count()))
        return false;
      return true;
    }
  } else if constexpr (is_total_occupation<charge_type>::value) {
    if constexpr (!std::is_same_v<
                      raw_type,
                      typename charge_type::raw_index_type>) {
      return false;
    } else {
      if (!cq.bound_to(layout.raw_index()) ||
          layout.kind() != ClusterLayoutKind::BosonModes ||
          !cq.matches_dimensions(layout.local_dimensions()))
        return false;
      return true;
    }
  } else {
    return false;
  }
}

template <typename Layout, typename CQ>
idxv_type cluster_charge_contribution(const Layout &layout, CQ &cq,
                                      size_t cluster, idxv_type local) {
  using charge_type = std::remove_cvref_t<CQ>;

  if constexpr (is_nset<charge_type>::value) {
    const bits_type active_mask(layout.raw_index().range() - 1);
    const bits_type mask = cq.bit_mask() & active_mask;
    const bits_type fragment(layout.local_raw_fragment(cluster, local));
    return (fragment & mask).count();
  } else if constexpr (is_nd<charge_type>::value) {
    return cq.counts_holes() ? idxv_type(local == 0)
                             : idxv_type(local == 3);
  } else if constexpr (is_total_occupation<charge_type>::value) {
    return local;
  } else {
    return 0;
  }
}

template <typename Layout, typename... CQs>
auto try_compile_cluster_constraints(const Layout &layout, CQs &...cqs) {
  constexpr size_t component_count = sizeof...(CQs);
  ClusterChargeProgram<component_count> program;

  if (!(supports_cluster_charge(layout, cqs) && ...))
    return std::pair<bool, ClusterChargeProgram<component_count>>(
        false, std::move(program));

  size_t component = 0;
  ((program.targets[component++] =
        static_cast<idxv_type>(cqs.val())),
   ...);

  size_t layer_count = 0;
  if (!try_size_add(layout.cluster_count(), size_t(1), layer_count))
    return std::pair<bool, ClusterChargeProgram<component_count>>(
        false, std::move(program));

  size_t row_count = 0;
  for (size_t cluster = 0; cluster < layout.cluster_count(); ++cluster) {
    const auto local_count = layout.local_state_count(cluster);
    if (local_count > std::numeric_limits<size_t>::max() ||
        !try_size_add(row_count, static_cast<size_t>(local_count),
                      row_count))
      return std::pair<bool, ClusterChargeProgram<component_count>>(
          false, std::move(program));
  }

  size_t program_bytes = 0;
  if (layer_count > program.cluster_offsets.max_size() ||
      row_count > program.delta.max_size() ||
      row_count > program.weighted_delta.max_size() ||
      row_count > program.target_feasible.max_size() ||
      !try_account_bytes(layer_count, sizeof(size_t), program_bytes) ||
      !try_account_bytes(row_count,
                         sizeof(typename decltype(program)::charge_type),
                         program_bytes) ||
      !try_account_bytes(row_count, sizeof(size_t), program_bytes) ||
      !try_account_bytes(row_count, sizeof(std::uint8_t), program_bytes) ||
      program_bytes > cluster_dp_memory_budget_storage())
    return std::pair<bool, ClusterChargeProgram<component_count>>(
        false, std::move(program));

  program.cluster_offsets.resize(layer_count);
  size_t row_offset = 0;
  for (size_t cluster = 0; cluster < layout.cluster_count(); ++cluster) {
    program.cluster_offsets[cluster] = row_offset;
    row_offset +=
        static_cast<size_t>(layout.local_state_count(cluster));
  }
  program.cluster_offsets[layout.cluster_count()] = row_offset;
  program.delta.resize(row_count);

  component = 0;
  const auto populate = [&](auto &charge) {
    for (size_t cluster = 0; cluster < layout.cluster_count(); ++cluster)
      for (idxv_type local = 0;
           local < layout.local_state_count(cluster); ++local)
        program.delta[program.row(cluster, local)][component] =
            cluster_charge_contribution(layout, charge, cluster, local);
    ++component;
  };
  (populate(cqs), ...);

  return std::pair<bool, ClusterChargeProgram<component_count>>(
      true, std::move(program));
}

template <size_t M>
idxv_type reachable_charge_max(const ClusterChargeProgram<M> &program,
                               size_t component) {
  idxv_type reachable = 0;
  for (size_t cluster = 0; cluster + 1 < program.cluster_offsets.size();
       ++cluster) {
    idxv_type local_max = 0;
    for (size_t row = program.cluster_offsets[cluster];
         row < program.cluster_offsets[cluster + 1]; ++row)
      local_max = std::max(local_max, program.delta[row][component]);
    if (local_max >
        std::numeric_limits<idxv_type>::max() - reachable) {
      return std::numeric_limits<idxv_type>::max();
    }
    reachable += local_max;
  }
  return reachable;
}

template <size_t M>
bool targets_fit_exact_evaluate_domain(
    const ClusterChargeProgram<M> &program) {
  const auto max_exact =
      static_cast<idxv_type>(std::numeric_limits<int>::max());
  for (size_t component = 0; component < M; ++component)
    if (program.targets[component] > max_exact ||
        reachable_charge_max(program, component) > max_exact)
      return false;
  return true;
}

template <size_t M>
bool target_exceeds_reachable_max(
    const ClusterChargeProgram<M> &program) {
  for (size_t component = 0; component < M; ++component)
    if (program.targets[component] >
        reachable_charge_max(program, component))
      return true;
  return false;
}

template <typename Layout, size_t M>
class ClusterSectorDP {
  using program_type = ClusterChargeProgram<M>;
  using charge_type = typename program_type::charge_type;

  Layout layout_;
  program_type program_;
  size_t charge_state_count_;
  std::vector<idxv_type> ways_;
  idxv_type sector_size_;

  void decode_charge(size_t index, charge_type &charge) const noexcept {
    for (size_t component = M; component-- > 0;) {
      charge[component] = index % program_.radices[component];
      index /= program_.radices[component];
    }
  }

  bool transition(const charge_type &charge, size_t charge_state,
                  size_t row, charge_type &result,
                  size_t &result_state) const noexcept {
    if (!program_.target_feasible[row]) return false;
    const auto &delta = program_.delta[row];
    // The flattened integer comparison alone is insufficient because it
    // could borrow across mixed-radix charge components.
    for (size_t component = 0; component < M; ++component) {
      if (charge[component] < delta[component]) return false;
      result[component] = charge[component] - delta[component];
    }
    assert(charge_state >= program_.weighted_delta[row]);
    result_state = charge_state - program_.weighted_delta[row];
    return true;
  }

  idxv_type suffix_count(size_t cluster, size_t charge_state) const noexcept {
    return ways_[cluster * charge_state_count_ + charge_state];
  }

  void enumerate_from(size_t cluster, idxv_type raw_prefix,
                      const charge_type &remaining,
                      size_t remaining_state,
                      std::vector<idxv_type> &mapping) const {
    if (cluster == layout_.cluster_count()) {
      mapping.push_back(raw_prefix);
      return;
    }

    charge_type next;
    for (idxv_type local = 0; local < layout_.local_state_count(cluster);
         ++local) {
      const size_t row = program_.row(cluster, local);
      size_t next_state = 0;
      if (!transition(remaining, remaining_state, row, next, next_state))
        continue;
      if (suffix_count(cluster + 1, next_state) == 0) continue;
      enumerate_from(
          cluster + 1,
          layout_.append_local_state(raw_prefix, cluster, local),
          next, next_state, mapping);
    }
  }

  ClusterSectorDP(Layout layout, program_type program)
      : layout_(std::move(layout)),
        program_(std::move(program)),
        charge_state_count_(program_.charge_state_count),
        ways_(checked_size_product(program_.cluster_offsets.size(),
                                   charge_state_count_),
              0),
        sector_size_(0) {
    ways_[layout_.cluster_count() * charge_state_count_] = 1;

    charge_type charge;
    charge_type previous;
    for (size_t cluster = layout_.cluster_count(); cluster-- > 0;) {
      for (size_t state = 0; state < charge_state_count_; ++state) {
        decode_charge(state, charge);
        idxv_type count = 0;
        for (idxv_type local = 0;
             local < layout_.local_state_count(cluster); ++local) {
          const size_t row = program_.row(cluster, local);
          size_t previous_state = 0;
          if (!transition(charge, state, row, previous, previous_state))
            continue;
          count = checked_count_add(
              count, ways_[(cluster + 1) * charge_state_count_ +
                           previous_state]);
        }
        ways_[cluster * charge_state_count_ + state] = count;
      }
    }
    sector_size_ = suffix_count(0, program_.target_state);
  }

 public:
  ClusterSectorDP(const ClusterSectorDP &) = default;
  ClusterSectorDP(ClusterSectorDP &&) noexcept = default;
  ClusterSectorDP &operator=(const ClusterSectorDP &) = default;
  ClusterSectorDP &operator=(ClusterSectorDP &&) noexcept = default;

  static std::optional<ClusterSectorDP> create(Layout layout,
                                                program_type program) {
    size_t charge_states = 1;
    for (size_t component = M; component-- > 0;) {
      const auto target = program.targets[component];
      if (target == std::numeric_limits<idxv_type>::max()) return {};
      const idxv_type radix = target + 1;
      if (radix > std::numeric_limits<size_t>::max()) return {};
      program.radices[component] = static_cast<size_t>(radix);
      program.weights[component] = charge_states;
      size_t next_charge_states = 0;
      if (!try_size_product(charge_states, static_cast<size_t>(radix),
                            next_charge_states))
        return {};
      charge_states = next_charge_states;
    }
    program.charge_state_count = charge_states;
    program.target_state = charge_states - 1;

    size_t entries = 0;
    size_t total_bytes = 0;
    const auto account = [&](size_t count, size_t item_size) {
      size_t bytes = 0;
      size_t next_total = 0;
      if (!try_size_product(count, item_size, bytes) ||
          !try_size_add(total_bytes, bytes, next_total))
        return false;
      total_bytes = next_total;
      return true;
    };
    if (!try_size_product(program.cluster_offsets.size(), charge_states,
                          entries))
      return {};
    if (entries > std::vector<idxv_type>().max_size() ||
        !account(program.cluster_offsets.size(), sizeof(size_t)) ||
        !account(program.delta.size(), sizeof(charge_type)) ||
        !account(program.delta.size(), sizeof(size_t)) ||
        !account(program.delta.size(), sizeof(std::uint8_t)) ||
        !account(entries, sizeof(idxv_type)) ||
        total_bytes > cluster_dp_memory_budget_storage())
      return {};

    program.weighted_delta.resize(program.delta.size());
    program.target_feasible.assign(program.delta.size(), std::uint8_t(1));
    for (size_t row = 0; row < program.delta.size(); ++row) {
      size_t decrement = 0;
      for (size_t component = 0; component < M; ++component) {
        if (program.delta[row][component] > program.targets[component]) {
          program.target_feasible[row] = 0;
          break;
        }
        size_t term = 0;
        size_t next_decrement = 0;
        if (!try_size_product(
                static_cast<size_t>(program.delta[row][component]),
                program.weights[component], term) ||
            !try_size_add(decrement, term, next_decrement))
          return {};
        decrement = next_decrement;
      }
      if (program.target_feasible[row])
        program.weighted_delta[row] = decrement;
    }

    return ClusterSectorDP(std::move(layout), std::move(program));
  }

  idxv_type size() const noexcept { return sector_size_; }

  std::vector<idxv_type> materialize() const {
    std::vector<idxv_type> mapping;
    if (sector_size_ > mapping.max_size())
      throw std::length_error("constrained basis exceeds vector capacity");
    mapping.reserve(static_cast<size_t>(sector_size_));
    enumerate_from(0, 0, program_.targets, program_.target_state, mapping);
    if (mapping.size() != sector_size_)
      throw std::logic_error("cluster DP enumeration count mismatch");
    return mapping;
  }

  idxv_type rank(idxv_type raw) const {
    if (raw >= layout_.range()) return sector_size_;

    charge_type remaining = program_.targets;
    charge_type next;
    size_t remaining_state = program_.target_state;
    idxv_type result = 0;

    for (size_t cluster = 0; cluster < layout_.cluster_count(); ++cluster) {
      const auto actual = layout_.local_state(raw, cluster);
      for (idxv_type local = 0; local < actual; ++local) {
        const size_t row = program_.row(cluster, local);
        size_t next_state = 0;
        if (!transition(remaining, remaining_state, row, next, next_state))
          continue;
        const idxv_type block = suffix_count(cluster + 1, next_state);
        // These are disjoint lexicographic prefix blocks, so their partial
        // sum is bounded by the already validated sector size.
        assert(block <= sector_size_ &&
               result <= sector_size_ - block);
        result += block;
      }
      const size_t actual_row = program_.row(cluster, actual);
      size_t next_state = 0;
      if (!transition(remaining, remaining_state, actual_row, next,
                      next_state))
        return sector_size_;
      if (suffix_count(cluster + 1, next_state) == 0) return sector_size_;
      remaining = next;
      remaining_state = next_state;
    }

    if (remaining_state != 0) return sector_size_;
    return result < sector_size_ ? result : sector_size_;
  }
};

template <typename Layout, size_t M>
class ClusterLexMapping {
  using dp_type = ClusterSectorDP<Layout, M>;
  struct EmptyDirectRank {};
  struct DirectTag {};
  struct EmptyDirectTag {};
  using rank_state =
      std::variant<InverseHash, dp_type, EmptyDirectRank>;

  std::vector<idxv_type> forward_;
  rank_state inverse_;

  ClusterLexMapping(std::vector<idxv_type> forward, DirectTag, dp_type dp)
      : forward_(std::move(forward)),
        inverse_(std::in_place_type<dp_type>, std::move(dp)) {}

  ClusterLexMapping(std::vector<idxv_type> forward, EmptyDirectTag)
      : forward_(std::move(forward)),
        inverse_(std::in_place_type<EmptyDirectRank>) {}

 public:
  explicit ClusterLexMapping(std::vector<idxv_type> forward)
      : forward_(std::move(forward)),
        inverse_(std::in_place_type<InverseHash>, forward_) {}

  static ClusterLexMapping direct(std::vector<idxv_type> forward,
                                  dp_type dp) {
    return ClusterLexMapping(std::move(forward), DirectTag{},
                             std::move(dp));
  }

  static ClusterLexMapping empty_direct() {
    return ClusterLexMapping({}, EmptyDirectTag{});
  }

  idxv_type size() const noexcept { return forward_.size(); }

  idxv_type raw_at(idxv_type compact) const noexcept {
    return forward_[compact];
  }

  idxv_type rank(idxv_type raw) const {
    switch (inverse_.index()) {
      case 0:
        return std::get<InverseHash>(inverse_).rank(raw, size());
      case 1:
        return std::get<dp_type>(inverse_).rank(raw);
      case 2:
        return size();
      default:
        throw std::logic_error("invalid cluster rank backend");
    }
  }

  bool has_direct_rank() const noexcept {
    return !std::holds_alternative<InverseHash>(inverse_);
  }
};

template <size_t tight, typename Layout, typename CQ, typename... CQs>
void exhaustive_cluster_visit(Layout &layout, size_t cluster,
                              idxv_type raw_prefix,
                              std::vector<idxv_type> &mapping, CQ &cq,
                              CQs &...cqs) {
  if (cluster == layout.cluster_count()) {
    layout.raw_index()[raw_prefix];
    if (Evaluate<tight>(cq, cqs...)) mapping.push_back(raw_prefix);
    return;
  }

  for (idxv_type local = 0; local < layout.local_state_count(cluster);
       ++local)
    exhaustive_cluster_visit<tight>(
        layout, cluster + 1,
        layout.append_local_state(raw_prefix, cluster, local), mapping, cq,
        cqs...);
}

template <size_t tight, typename Layout, typename CQ, typename... CQs>
std::vector<idxv_type> exhaustive_cluster_mapping(Layout &layout, CQ &cq,
                                                  CQs &...cqs) {
  std::vector<idxv_type> mapping;

  if (layout.kind() == ClusterLayoutKind::HubbardSites) {
    exhaustive_cluster_visit<tight>(layout, 0, 0, mapping, cq, cqs...);
  } else {
    // Spin-bit and bosonic mixed-radix cluster order is raw index order.
    for (idxv_type raw = 0; raw < layout.range(); ++raw) {
      layout.raw_index()[raw];
      if (Evaluate<tight>(cq, cqs...)) mapping.push_back(raw);
    }
  }
  return mapping;
}

}  // namespace detail

inline size_t clusterDPMemoryBudget() {
  return detail::cluster_dp_memory_budget_storage();
}

inline void setClusterDPMemoryBudget(size_t bytes) {
  detail::cluster_dp_memory_budget_storage() = bytes;
}

template <size_t tight, bool P, typename Layout, typename CQ, typename... CQs>
  requires ClusterLayoutLike<Layout> &&
           std::is_lvalue_reference_v<Layout &&>
auto LooseConstrain(Layout &&layout, CQ &&cq, CQs &&...cqs) {
  using layout_type = std::remove_cvref_t<Layout>;
  using raw_reference = decltype((layout.raw_index()));
  constexpr size_t component_count = 1 + sizeof...(CQs);
  using dp_type = detail::ClusterSectorDP<layout_type, component_count>;
  using backend_type =
      detail::ClusterLexMapping<layout_type, component_count>;

  std::vector<idxv_type> mapping;
  std::optional<dp_type> direct_dp;
  bool empty_direct = false;
  bool used_dp = false;

  if constexpr (tight == 0) {
    auto [supported, program] =
        detail::try_compile_cluster_constraints(layout, cq, cqs...);
    if (supported && detail::targets_fit_exact_evaluate_domain(program)) {
      if (detail::target_exceeds_reachable_max(program)) {
        used_dp = true;
        empty_direct = layout.uses_direct_rank();
      } else {
        auto dp = dp_type::create(layout, std::move(program));
        if (dp) {
          used_dp = true;
          mapping = dp->materialize();
          if (layout.uses_direct_rank())
            direct_dp.emplace(std::move(*dp));
        }
      }
    }
  }

  if (!used_dp)
    mapping =
        detail::exhaustive_cluster_mapping<tight>(layout, cq, cqs...);

  auto backend = [&]() -> backend_type {
    if (empty_direct) return backend_type::empty_direct();
    if (direct_dp)
      return backend_type::direct(std::move(mapping),
                                  std::move(*direct_dp));
    return backend_type(std::move(mapping));
  }();
  auto sub = SubIndex<raw_reference, P, backend_type>(
      layout.raw_index(), std::move(backend));
  sub[0];
  return sub;
}

template <typename Layout, typename... CQs>
  requires ClusterLayoutLike<Layout> &&
           std::is_lvalue_reference_v<Layout &&>
auto Constrain(Layout &&layout, CQs &&...cqs) {
  return LooseConstrain<0, true>(std::forward<Layout>(layout),
                                 std::forward<CQs>(cqs)...);
}

}  // namespace qudrip
