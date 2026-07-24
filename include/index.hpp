#include <algorithm>
#include <array>
#include <bitset>
#include <complex>
#include <concepts>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#pragma once
namespace qudrip {
template <typename T>
concept IndexLike = requires {
  typename std::remove_reference_t<T>::is_index;
};

template <typename T>
concept ClusterLayoutLike = requires {
  typename std::remove_reference_t<T>::is_cluster_layout;
};

template <typename T>
concept PlainIndexLike = IndexLike<T> && (!ClusterLayoutLike<T>);

template <class T>
class QbitsIndex;

template <class T>
class SingleModeIndex;

template <class T1, class T2>
class Indices;

namespace detail {
class MaterializedMapping;
}

template <typename T, bool P = false,
          typename Backend = detail::MaterializedMapping>
class SubIndex;

namespace detail {

inline idx_size_t checked_index_product(idx_size_t lhs, idx_size_t rhs) {
  if (lhs != 0 && rhs > std::numeric_limits<idx_size_t>::max() / lhs)
    throw std::overflow_error("tensor-product index range overflow");
  return lhs * rhs;
}

inline idxv_type checked_raw_add(idxv_type lhs, idxv_type rhs) {
  if (rhs > std::numeric_limits<idxv_type>::max() - lhs)
    throw std::overflow_error("raw leaf replacement overflow");
  return lhs + rhs;
}

inline idxv_type checked_raw_subtract(idxv_type lhs, idxv_type rhs) {
  if (rhs > lhs)
    throw std::overflow_error("raw leaf replacement underflow");
  return lhs - rhs;
}

enum class RawLeafKind { Unknown, Qbits, SingleMode };

struct RawLeafBinding {
  const void *identity = nullptr;
  idxv_type range = 0;
  idxv_type stride = 0;
  RawLeafKind kind = RawLeafKind::Unknown;

  idxv_type extract(idxv_type raw) const {
    if (range == 0 || stride == 0)
      throw std::logic_error("invalid raw leaf binding");
    return (raw / stride) % range;
  }

  idxv_type replace(idxv_type raw, idxv_type old_local,
                    idxv_type new_local) const {
    if (old_local >= range || new_local >= range)
      throw std::out_of_range("raw leaf state outside local range");
    if (extract(raw) != old_local)
      throw std::invalid_argument(
          "raw leaf replacement does not match the encoded state");

    if (new_local > old_local) {
      const auto increment = checked_index_product(
          new_local - old_local, stride);
      return checked_raw_add(raw, increment);
    }
    if (old_local > new_local) {
      const auto decrement = checked_index_product(
          old_local - new_local, stride);
      return checked_raw_subtract(raw, decrement);
    }
    return raw;
  }
};

using RawLeafBindings = std::vector<RawLeafBinding>;

template <class T>
void append_raw_leaf_bindings(const QbitsIndex<T> &index,
                              idxv_type stride,
                              RawLeafBindings &bindings);

template <class T>
void append_raw_leaf_bindings(const SingleModeIndex<T> &index,
                              idxv_type stride,
                              RawLeafBindings &bindings);

template <class T1, class T2>
void append_raw_leaf_bindings(const Indices<T1, T2> &index,
                              idxv_type stride,
                              RawLeafBindings &bindings);

template <typename T, bool P, typename Backend>
void append_raw_leaf_bindings(const SubIndex<T, P, Backend> &index,
                              idxv_type stride,
                              RawLeafBindings &bindings);

template <typename IDX>
RawLeafBindings collect_raw_leaf_bindings(const IDX &index) {
  RawLeafBindings bindings;
  append_raw_leaf_bindings(index, 1, bindings);
  return bindings;
}

inline const RawLeafBinding *find_unique_raw_leaf_binding(
    const RawLeafBindings &bindings, const void *identity) noexcept {
  const RawLeafBinding *match = nullptr;
  for (const auto &binding : bindings) {
    if (binding.identity != identity) continue;
    if (match != nullptr) return nullptr;
    match = std::addressof(binding);
  }
  return match;
}
}  // namespace detail

template <size_t = 0, bool = false, typename IDX, typename CQ, typename... CQs>
  requires PlainIndexLike<IDX>
auto LooseConstrain(IDX &&index, CQ &&cq, CQs &&...cqs);

template <typename IDX, typename... CQs>
  requires PlainIndexLike<IDX>
auto Constrain(IDX &&index, CQs &&...cqs);

template <size_t = 0, bool = false, typename Layout, typename CQ,
          typename... CQs>
  requires ClusterLayoutLike<Layout> &&
           std::is_lvalue_reference_v<Layout &&>
auto LooseConstrain(Layout &&layout, CQ &&cq, CQs &&...cqs);

template <typename Layout, typename... CQs>
  requires ClusterLayoutLike<Layout> &&
           std::is_lvalue_reference_v<Layout &&>
auto Constrain(Layout &&layout, CQs &&...cqs);

template <typename IDX, typename std::remove_reference_t<IDX>::ear::type = 0>
decltype(auto) strip(IDX &&index);

template <typename IDX,
          typename std::remove_reference_t<IDX>::is_earless::type = 0>
decltype(auto) strip(IDX &&index);

template <typename IDX, typename std::remove_reference_t<IDX>::ear::type,
          typename = typename std::enable_if<IDX::pretend_ == true>::type>
decltype(auto) strip(IDX &&index);

using bits_type = std::bitset<BIT_LIMIT>;

// CRTP
template <class Derived, class T = idx_size_t>
class IndexBase {
  T range_;

 protected:
  T idx_value_;

 public:
  struct is_index {
    using type = int;
  };
  struct is_earless {
    using type = int;
  };

  IndexBase(T range) : range_(range), idx_value_(0) {}

  operator T() const { return idx_value_; }

  void operator[](T idx_val) { idx_value_ = idx_val; }

  template <class IDX>
  Derived &operator=(IDX &&idx) noexcept {
    static_cast<Derived *>(this)->t_idx(std::forward<IDX>(idx));
    return (*static_cast<Derived *>(this));
  }

  Derived &operator++() noexcept {
    ++idx_value_;
    return (*static_cast<Derived *>(this));
  }

  T range() const { return range_; }

  idxv_type raw_at(idxv_type compact) const noexcept {
    const auto missing = static_cast<idxv_type>(range_);
    return compact < missing ? compact : missing;
  }

  idxv_type rank(idxv_type raw) const noexcept {
    const auto missing = static_cast<idxv_type>(range_);
    return raw < missing ? raw : missing;
  }
};
//------------------------------------------------------------------

template <typename int_type>
int_type bit_convert(const bits_type &idx) {
  return idx.to_ulong();
}
template <>
inline unsigned long long bit_convert(const bits_type &idx) {
  return idx.to_ullong();
}

template <class T = idx_size_t>
class QbitsIndex : public IndexBase<QbitsIndex<T>, T> {
  using base = IndexBase<QbitsIndex<T>, T>;
  friend class IndexBase<QbitsIndex<T>, T>;

  static T checked_range(T bit_count) {
    if constexpr (std::is_signed_v<T>) {
      if (bit_count < 0)
        throw std::invalid_argument("Qbit count must be nonnegative");
    }
    if (static_cast<std::make_unsigned_t<T>>(bit_count) >=
        static_cast<std::make_unsigned_t<T>>(
            std::numeric_limits<T>::digits))
      throw std::overflow_error("Qbit index range is not representable");
    return T(1) << bit_count;
  }

  void t_idx(const bits_type &idx) { base::idx_value_ = bit_convert<T>(idx); }
  void t_idx(std::string idx) { t_idx(bits_type(idx)); }

  bits_type content() const { return bits_type(base::idx_value_); }

 public:
  struct is_qbits {};
  using IndexBase<QbitsIndex<T>, T>::operator=;

  QbitsIndex(T N) : IndexBase<QbitsIndex<T>, T>(checked_range(N)) {}
};
//------------------------------------------------------------------

template <class T = idx_size_t>
class SingleModeIndex : public IndexBase<SingleModeIndex<T>, T> {
  using base = IndexBase<SingleModeIndex<T>, T>;
  friend class IndexBase<SingleModeIndex<T>, T>;
  int n_mode_;
  static T checked_range(int n_mode) {
    if (n_mode <= 0)
      throw std::invalid_argument("mode range must be positive");
    if (static_cast<unsigned long long>(n_mode) >
        static_cast<unsigned long long>(std::numeric_limits<T>::max()))
      throw std::overflow_error("mode range is not representable");
    return static_cast<T>(n_mode);
  }
  void t_idx(T N) { base::idx_value_ = N; }
  T content() const { return base::idx_value_; }

 public:
  struct is_single_mode {};
  using IndexBase<SingleModeIndex<T>, T>::operator=;
  SingleModeIndex(int n_mode)
      : IndexBase<SingleModeIndex<T>, T>(checked_range(n_mode)),
        n_mode_(n_mode) {}
};
//------------------------------------------------------------------

namespace detail {

template <typename IDX>
void append_index_ranges(const IDX &index, std::vector<idx_size_t> &ranges) {
  if constexpr (requires { index.append_local_ranges(ranges); })
    index.append_local_ranges(ranges);
  else
    ranges.push_back(index.range());
}

template <typename IDX>
constexpr bool all_single_mode_indices(const IDX &index) {
  if constexpr (requires { typename IDX::is_single_mode; })
    return true;
  else if constexpr (requires { index.all_single_mode_indices(); })
    return index.all_single_mode_indices();
  else
    return false;
}

}  // namespace detail

//------------------------------------------------------------------

template <class T1, class T2>
class Indices {
  template <class L, class R>
  friend void detail::append_raw_leaf_bindings(
      const Indices<L, R> &index, idxv_type stride,
      detail::RawLeafBindings &bindings);

  using lhs_type = T1;
  using rhs_type = T2;

  lrHandle<lhs_type> LH;
  lrHandle<rhs_type> RH;
  typename std::remove_reference<lhs_type>::type &lhs_;
  typename std::remove_reference<rhs_type>::type &rhs_;
  idx_size_t range_;

 public:
  struct is_indices {
    using type = int;
  };
  struct is_index {
    using type = int;
  };
  struct is_earless {
    using type = int;
  };

  Indices(lhs_type lhs, rhs_type rhs)
      : LH(lhs),
        RH(rhs),
        lhs_(LH.val_),
        rhs_(RH.val_),
        range_(detail::checked_index_product(lhs_.range(), rhs_.range())) {}

  Indices(const Indices &idx)
      : LH(idx.LH),
        RH(idx.RH),
        lhs_(LH.val_),
        rhs_(RH.val_),
        range_(detail::checked_index_product(lhs_.range(), rhs_.range())) {}

  Indices(Indices &&idx)
      : LH(std::move(idx.LH)),
        RH(std::move(idx.RH)),
        lhs_(LH.val_),
        rhs_(RH.val_),
        range_(detail::checked_index_product(lhs_.range(), rhs_.range())) {}

  operator idx_size_t() const {
    return idx_size_t(lhs_) * (rhs_.range()) + idx_size_t(rhs_);
  }

  idx_size_t range() const { return range_; }

  idxv_type raw_at(idxv_type compact) const noexcept {
    return compact < range_ ? compact : range_;
  }

  idxv_type rank(idxv_type raw) const noexcept {
    return raw < range_ ? raw : range_;
  }

  void append_local_ranges(std::vector<idx_size_t> &ranges) const {
    detail::append_index_ranges(lhs_, ranges);
    detail::append_index_ranges(rhs_, ranges);
  }

  std::vector<idx_size_t> local_ranges() const {
    std::vector<idx_size_t> ranges;
    append_local_ranges(ranges);
    return ranges;
  }

  constexpr bool all_single_mode_indices() const {
    return detail::all_single_mode_indices(lhs_) &&
           detail::all_single_mode_indices(rhs_);
  }

  //=========================================================================================
  void operator[](idx_size_t i) {
    lhs_[i / rhs_.range()];
    rhs_[i % rhs_.range()];
  }

  //=========================================================================================
  Indices &operator++() {
    if (rhs_ < rhs_.range() - 1)
      ++rhs_;
    else {
      rhs_[0];
      ++lhs_;
    }
    return *this;
  }
};

//------------------------------------------------------------------
template <class T1, class T2,
          typename = typename std::remove_reference<T1>::type::is_index,
          typename = typename std::remove_reference<T2>::type::is_index>
Indices<T1, T2> operator*(T1 &&lhs, T2 &&rhs) {
  return Indices<T1, T2>(lhs, rhs);
}

//------------------------------------------------------------------
namespace detail {

class InverseHash {
  std::unordered_map<idxv_type, idxv_type> mapping_;

 public:
  InverseHash() = default;

  explicit InverseHash(const std::vector<idxv_type> &forward) {
    mapping_.reserve(forward.size());
    for (idxv_type compact = 0; compact < forward.size(); ++compact) {
      const auto inserted = mapping_.emplace(forward[compact], compact);
      if (!inserted.second)
        throw std::logic_error("duplicate raw index in constrained basis");
    }
  }

  idxv_type rank(idxv_type raw, idxv_type missing) const noexcept {
    const auto found = mapping_.find(raw);
    return found == mapping_.end() ? missing : found->second;
  }
};

class MaterializedMapping {
  std::vector<idxv_type> forward_;
  InverseHash inverse_;

 public:
  MaterializedMapping() = default;

  explicit MaterializedMapping(std::vector<idxv_type> forward)
      : forward_(std::move(forward)), inverse_(forward_) {}

  idxv_type size() const noexcept { return forward_.size(); }

  idxv_type raw_at(idxv_type compact) const noexcept {
    return forward_[compact];
  }

  idxv_type rank(idxv_type raw) const noexcept {
    return inverse_.rank(raw, size());
  }

  constexpr bool has_direct_rank() const noexcept { return false; }
};

}  // namespace detail

template <typename T, bool P, typename Backend>
class SubIndex {
  template <typename Wrapped, bool Pretend, typename MappingBackend>
  friend void detail::append_raw_leaf_bindings(
      const SubIndex<Wrapped, Pretend, MappingBackend> &index,
      idxv_type stride, detail::RawLeafBindings &bindings);

  template <typename IDX, typename std::remove_reference_t<IDX>::ear::type>
  friend decltype(auto) strip(IDX &&index);

  template <typename IDX,
            typename std::remove_reference_t<IDX>::is_earless::type>
  friend decltype(auto) strip(IDX &&index);

  template <typename IDX, typename std::remove_reference_t<IDX>::ear::type,
            typename>
  friend decltype(auto) strip(IDX &&index);

  using index_type = T;
  lrHandle<index_type> IH;
  typename std::remove_reference<index_type>::type &earless_idx_;
  Backend backend_;
  idxv_type idx_value_ = 0;

 public:
  const static bool pretend_ = P;
  struct is_index {
    using type = int;
  };
  struct ear {
    using type = int;
  };

  explicit SubIndex(index_type idx)
    requires std::default_initializable<Backend>
      : IH(idx), earless_idx_(IH.val_), backend_() {}

  SubIndex(index_type idx, Backend backend)
      : IH(idx),
        earless_idx_(IH.val_),
        backend_(std::move(backend)) {}

  SubIndex(index_type idx, std::vector<idxv_type> raw_mapping)
    requires std::same_as<Backend, detail::MaterializedMapping>
      : SubIndex(idx, Backend(std::move(raw_mapping))) {}

  SubIndex(const SubIndex &sub)
      : IH(sub.IH),
        earless_idx_(IH.val_),
        backend_(sub.backend_),
        idx_value_(sub.idx_value_) {}

  SubIndex(SubIndex &&sub)
      : IH(std::move(sub.IH)),
        earless_idx_(IH.val_),
        backend_(std::move(sub.backend_)),
        idx_value_(sub.idx_value_) {}

  operator idxv_type() const {
    const idxv_type raw = idxv_type(earless_idx_);
    if (idx_value_ < backend_.size() &&
        backend_.raw_at(idx_value_) == raw)
      return idx_value_;
    return backend_.rank(raw);
  }

  void operator[](idxv_type i) {
    idx_value_ = i;
    if (idx_value_ < backend_.size())
      earless_idx_[backend_.raw_at(idx_value_)];
    else
      earless_idx_[earless_idx_.range()];
  }

  SubIndex &operator++() {
    (*this)[idx_value_ + 1];
    return *this;
  }

  idxv_type range() const noexcept { return backend_.size(); }

  idxv_type raw_at(idxv_type compact) const
      noexcept(noexcept(backend_.raw_at(compact))) {
    return compact < backend_.size() ? backend_.raw_at(compact)
                                     : backend_.size();
  }

  idxv_type rank(idxv_type raw) const
      noexcept(noexcept(backend_.rank(raw))) {
    return backend_.rank(raw);
  }

  bool has_direct_rank() const noexcept {
    return backend_.has_direct_rank();
  }
};

//------------------------------------------------------------------
namespace detail {

template <typename Primitive>
void append_primitive_raw_leaf_binding(
    const Primitive &index, idxv_type stride,
    RawLeafKind kind, RawLeafBindings &bindings) {
  const auto local_range = static_cast<idxv_type>(index.range());
  if (local_range == 0)
    throw std::invalid_argument("raw leaf range must be positive");
  if (stride == 0)
    throw std::invalid_argument("raw leaf stride must be positive");

  // The largest contribution of this leaf must fit in the raw label.
  (void)checked_index_product(local_range - 1, stride);
  bindings.push_back(
      {static_cast<const void *>(std::addressof(index)), local_range,
       stride, kind});
}

template <class T>
void append_raw_leaf_bindings(const QbitsIndex<T> &index,
                              idxv_type stride,
                              RawLeafBindings &bindings) {
  append_primitive_raw_leaf_binding(
      index, stride, RawLeafKind::Qbits, bindings);
}

template <class T>
void append_raw_leaf_bindings(const SingleModeIndex<T> &index,
                              idxv_type stride,
                              RawLeafBindings &bindings) {
  append_primitive_raw_leaf_binding(
      index, stride, RawLeafKind::SingleMode, bindings);
}

template <class T1, class T2>
void append_raw_leaf_bindings(const Indices<T1, T2> &index,
                              idxv_type stride,
                              RawLeafBindings &bindings) {
  const auto rhs_range = static_cast<idxv_type>(index.rhs_.range());
  if (rhs_range == 0)
    throw std::invalid_argument("raw product leaf range must be positive");
  const auto lhs_stride = checked_index_product(stride, rhs_range);
  append_raw_leaf_bindings(index.lhs_, lhs_stride, bindings);
  append_raw_leaf_bindings(index.rhs_, stride, bindings);
}

template <typename T, bool P, typename Backend>
void append_raw_leaf_bindings(const SubIndex<T, P, Backend> &index,
                              idxv_type stride,
                              RawLeafBindings &bindings) {
  append_raw_leaf_bindings(index.earless_idx_, stride, bindings);
}

}  // namespace detail

//------------------------------------------------------------------
template <size_t start, typename T, size_t N, typename = typename T::is_index,
          typename std::enable_if<(start < N - 1), int>::type = 0>
constexpr auto product(std::array<T, N> &idxs) {
  return idxs[start] * product<start + 1>(idxs);
}

template <size_t start, typename T, size_t N, typename = typename T::is_index,
          typename std::enable_if<(start == N - 1), int>::type = 0>
constexpr T product(std::array<T, N> &idxs) {
  return idxs[start];
}

template <typename T, size_t N, typename = typename T::is_index>
constexpr auto product(std::array<T, N> &idxs) {
  return product<0>(idxs);
}

//------------------------------------------------------------------

template <typename IDX, typename std::remove_reference_t<IDX>::ear::type,
          typename>
decltype(auto) strip(IDX &&index) {
  return std::forward<IDX>(index);
}

template <typename IDX, typename std::remove_reference_t<IDX>::ear::type>
decltype(auto) strip(IDX &&index) {
  return strip(index.earless_idx_);
}

template <typename IDX, typename std::remove_reference_t<IDX>::is_earless::type>
decltype(auto) strip(IDX &&index) {
  return index;
}

//------------------------------------------------------------------
inline auto getQbits(int N) {
  if (N < 0)
    throw std::invalid_argument("Qbit count must be nonnegative");
  return QbitsIndex<>(static_cast<idx_size_t>(N));
}

inline auto getSingleMode(int Nph) { return SingleModeIndex<>(Nph); }

//------------------------------------------------------------------
inline bits_type getQbitsLabel(int N) { return bits_type(N); }

//
}  // namespace qudrip
