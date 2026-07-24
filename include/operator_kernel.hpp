#pragma once

#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace qudrip::detail {

static_assert(std::is_integral_v<idxv_type> && std::is_unsigned_v<idxv_type>,
              "operator kernels require an unsigned integral raw index");

inline constexpr std::size_t raw_index_bits =
    std::numeric_limits<idxv_type>::digits;

// Stateless paths use raw zero as an ordinary label. A dead branch is
// represented by not emitting a RawTransition.
struct RawTransition {
  idxv_type raw{};
  value_type amplitude{};
};

enum class FactorKind : std::uint8_t { BitMatrix, HopBit, ModeMatrix };

// Layout-independent immutable description of one built-in elementary
// operator. The matrix span borrows from the shared elementary operator, whose
// lifetime is owned by Operator::ops_ for the duration of a bound plan.
struct FactorDescriptor {
  FactorKind kind{FactorKind::ModeMatrix};
  const void *primitive_identity{nullptr};
  std::size_t input_range{};
  std::size_t output_range{};
  std::span<const value_type> matrix;
  std::size_t site{};
  std::size_t site2{};
  int string_value{1};

  value_type matrix_element(std::size_t input,
                            std::size_t output) const noexcept {
    return matrix[input * output_range + output];
  }
};

constexpr idxv_type low_bits_mask(std::size_t bit_count) noexcept {
  if (bit_count == 0)
    return idxv_type{0};
  if (bit_count >= raw_index_bits)
    return std::numeric_limits<idxv_type>::max();
  return (idxv_type{1} << bit_count) - idxv_type{1};
}

inline idxv_type checked_site_mask(std::size_t site, idxv_type active_mask) {
  if (site >= raw_index_bits)
    throw std::out_of_range("Qbit site exceeds raw-index width");
  const idxv_type mask = idxv_type{1} << site;
  if ((mask & active_mask) == 0)
    throw std::out_of_range("Qbit site is outside the active domain");
  return mask;
}

constexpr idxv_type higher_site_mask(std::size_t site,
                                     idxv_type active_mask) noexcept {
  if (site >= raw_index_bits - 1)
    return idxv_type{0};
  return active_mask & (std::numeric_limits<idxv_type>::max() << (site + 1));
}

constexpr idxv_type strictly_between_site_mask(std::size_t site_a,
                                               std::size_t site_b,
                                               idxv_type active_mask) noexcept {
  const std::size_t lower = site_a < site_b ? site_a : site_b;
  const std::size_t upper = site_a < site_b ? site_b : site_a;
  if (lower >= raw_index_bits - 1 || upper <= lower + 1)
    return idxv_type{0};
  return higher_site_mask(lower, active_mask) & low_bits_mask(upper);
}

// For an active-domain raw label x, this represents
//
//   valid(x) = (x & must_set) == must_set
//           && (x & must_clear) == 0
//   output(x) = x ^ flip
//   amplitude(x) = scale * (-1)^popcount(x & parity)
//
// All masks are expressed in the local, low-bit Qbits domain. Layout binding
// is deliberately kept outside this math-only type.
struct BitMaskKernel {
  idxv_type active_mask{};
  idxv_type must_set{};
  idxv_type must_clear{};
  idxv_type flip{};
  idxv_type parity{};
  value_type scale{1.0};
  bool identically_zero{false};
};

constexpr bool
masks_within_active_domain(const BitMaskKernel &kernel) noexcept {
  const idxv_type used =
      kernel.must_set | kernel.must_clear | kernel.flip | kernel.parity;
  return (used & ~kernel.active_mask) == 0;
}

inline unsigned parity_count(idxv_type raw, idxv_type mask) noexcept {
  return std::popcount(raw & mask) & 1U;
}

inline int parity_sign(idxv_type raw, idxv_type mask) noexcept {
  return parity_count(raw, mask) == 0 ? 1 : -1;
}

inline BitMaskKernel zero_bit_mask_kernel(idxv_type active_mask) noexcept {
  BitMaskKernel result;
  result.active_mask = active_mask;
  result.scale = value_type{};
  result.identically_zero = true;
  return result;
}

// Put a kernel into its unique useful form: conflicts and an exact-zero scale
// become the zero kernel, while parity on fixed inputs becomes a constant
// phase folded into scale.
inline BitMaskKernel canonicalize_bit_mask_kernel(BitMaskKernel kernel) {
  if (!masks_within_active_domain(kernel))
    throw std::invalid_argument(
        "bit-mask kernel uses bits outside its active domain");

  if (kernel.identically_zero || kernel.scale == value_type{} ||
      (kernel.must_set & kernel.must_clear) != 0)
    return zero_bit_mask_kernel(kernel.active_mask);

  if (parity_count(kernel.must_set, kernel.parity) != 0)
    kernel.scale = -kernel.scale;
  kernel.parity &= ~(kernel.must_set | kernel.must_clear);
  kernel.parity &= kernel.active_mask;
  return kernel;
}

// Compose A after B. This direction matches the mathematical product A * B:
// B sees the input first and A sees B's output.
inline BitMaskKernel compose_after(const BitMaskKernel &after,
                                   const BitMaskKernel &before) {
  if (after.active_mask != before.active_mask)
    throw std::invalid_argument(
        "cannot compose bit-mask kernels from different active domains");

  const BitMaskKernel a = canonicalize_bit_mask_kernel(after);
  const BitMaskKernel b = canonicalize_bit_mask_kernel(before);
  const idxv_type active = a.active_mask;

  if (a.identically_zero || b.identically_zero)
    return zero_bit_mask_kernel(active);

  const idxv_type unchanged_by_b = active & ~b.flip;

  BitMaskKernel result;
  result.active_mask = active;
  result.must_set =
      b.must_set | (a.must_set & unchanged_by_b) | (a.must_clear & b.flip);
  result.must_clear =
      b.must_clear | (a.must_clear & unchanged_by_b) | (a.must_set & b.flip);
  result.flip = b.flip ^ a.flip;
  result.parity = b.parity ^ a.parity;
  result.scale = a.scale * b.scale;
  if (parity_count(b.flip, a.parity) != 0)
    result.scale = -result.scale;

  return canonicalize_bit_mask_kernel(std::move(result));
}

inline bool bit_mask_input_in_domain(const BitMaskKernel &kernel,
                                     idxv_type raw) noexcept {
  return (raw & ~kernel.active_mask) == 0;
}

inline bool bit_mask_accepts(const BitMaskKernel &kernel,
                             idxv_type raw) noexcept {
  return !kernel.identically_zero && bit_mask_input_in_domain(kernel, raw) &&
         (raw & kernel.must_set) == kernel.must_set &&
         (raw & kernel.must_clear) == 0;
}

// This applies structural, exact-zero pruning only. The caller owns the
// terminal eps policy and compact-basis rank.
inline bool apply_bit_mask(const BitMaskKernel &kernel, idxv_type raw_in,
                           value_type amplitude_in,
                           RawTransition &output) noexcept {
  if (!bit_mask_accepts(kernel, raw_in))
    return false;

  value_type amplitude = amplitude_in * kernel.scale;
  if (parity_count(raw_in, kernel.parity) != 0)
    amplitude = -amplitude;
  if (amplitude == value_type{})
    return false;

  output.raw = raw_in ^ kernel.flip;
  output.amplitude = amplitude;
  return true;
}

// A deterministic table has at most one exact structural output per local
// input. output_range is also its dead-output sentinel.
struct DeterministicTable {
  idxv_type input_range{};
  idxv_type output_range{};
  std::vector<idxv_type> output;
  std::vector<value_type> amplitude;
};

enum class DeterministicTarget : std::uint8_t { QbitSite, Mode };

enum class ParityEvaluation : std::uint8_t { Input, Output };

// A phase is evaluated on either the state entering its stage or the state
// emitted by that stage. Multiple entries remain in their original order.
// odd_multiplier is normally -1 for fermions and 1 for hard-core bosons, but
// it can also retain a custom elementary-operator string value.
struct OrderedParityPhase {
  ParityEvaluation evaluation{ParityEvaluation::Input};
  idxv_type active_mask{};
  idxv_type parity_mask{};
  value_type odd_multiplier{1.0};
};

struct DeterministicStage {
  DeterministicTarget target{DeterministicTarget::Mode};
  const void *primitive_identity{nullptr};

  // QbitSite: the bit site within the referenced Qbits primitive.
  // Mode: unused and kept zero.
  std::size_t coordinate{};

  DeterministicTable table;
  std::vector<OrderedParityPhase> phases;
};

enum class KernelKind : std::uint8_t { BitMask, Deterministic, Dfs, Legacy };

} // namespace qudrip::detail
