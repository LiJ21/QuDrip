//------------------------------------------------------------------

#pragma once
#include "operator_kernel.hpp"

namespace qudrip {
using idx_tree_type = std::vector<idxv_type>;
using val_tree_type = std::vector<value_type>;
using idx_it_type = idx_tree_type::iterator;
using val_it_type = val_tree_type::iterator;

//=========================================================================

class Operator;
class OperatorSum;

namespace detail {
inline size_t checked_factor_matrix_size(size_t input, size_t output) {
  if (input != 0 &&
      output > std::numeric_limits<size_t>::max() / input)
    throw std::length_error("elementary operator matrix size overflow");
  return input * output;
}
}  // namespace detail

using U_out_type = std::vector<value_type>;
class elOpBase {
 protected:
  U_out_type U_;
  size_t n_input, n_output;

 public:
  elOpBase(size_t n_input, size_t n_output)
      : U_(detail::checked_factor_matrix_size(n_input, n_output)),
        n_input(n_input),
        n_output(n_output) {}
  virtual ~elOpBase() = default;

  // U_ is laid out as U_[in * n_output + out] = U(out, in), i.e. the
  // amplitude <out|U|in>, matching how branch() consumes it.
  virtual void operator<<(Matrix U) {
    assert(U.rows() == n_output);
    assert(U.cols() == n_input);
    for (size_t i = 0; i < n_input; ++i)
      for (size_t j = 0; j < n_output; ++j) U_[i * n_output + j] = U(j, i);
  }

  virtual void invert() {
    assert(n_input == n_output);
    Matrix M(n_input, n_output);
    for (size_t i = 0; i < n_input; ++i)
      for (size_t j = 0; j < n_output; ++j) M(j, i) = U_[i * n_output + j];
    operator<<(M.inverse().eval());
  }

  virtual void feed_idx(idx_size_t idx) = 0;
  virtual idx_size_t get_idx() = 0;
  virtual void branch(idx_it_type &, val_it_type &, idx_it_type &,
                      val_it_type &) = 0;
  virtual bool describe_factor(detail::FactorDescriptor &) const noexcept {
    return false;
  }
  size_t range() const noexcept { return n_output; }
};

using ops_type = std::vector<std::shared_ptr<elOpBase>>;
//=========================================================================
// auxiliary function
inline size_t tree_size(const ops_type &ops) {
  size_t res = 1;
  for (const auto &op : ops) {
    const auto width = op->range();
    if (width != 0 &&
        res > (std::numeric_limits<size_t>::max() - 1) / width)
      throw std::length_error("operator expansion tree size overflow");
    res = 1 + res * width;
  }
  return res;
}

namespace detail {

inline constexpr idxv_type eigen_sparse_storage_limit() noexcept {
  using storage_index_type = SpMatrix::StorageIndex;
  static_assert(std::is_integral_v<storage_index_type>);
  static_assert(std::numeric_limits<storage_index_type>::max() >= 0);
  return static_cast<idxv_type>(
      std::numeric_limits<storage_index_type>::max());
}

inline Eigen::Index checked_eigen_sparse_index(
    idxv_type value, const char *message) {
  constexpr auto eigen_index_limit = static_cast<idxv_type>(
      std::numeric_limits<Eigen::Index>::max());
  if (value > eigen_sparse_storage_limit() ||
      value > eigen_index_limit)
    throw std::length_error(message);
  return static_cast<Eigen::Index>(value);
}

struct LegacyTreeWorkspace {
  idx_tree_type indices;
  val_tree_type amplitudes;

  explicit LegacyTreeWorkspace(const ops_type &ops) {
    const auto size = tree_size(ops);
    indices.resize(size);
    amplitudes.resize(size);
  }
};

struct BoundFactor {
  FactorDescriptor descriptor;
  RawLeafBinding binding;
  size_t scratch_offset{};
  bool direct_leaf{};
};

struct BoundDeterministicStep {
  DeterministicStage stage;
  RawLeafBinding binding;
  FactorKind kind{FactorKind::ModeMatrix};
  size_t site2{};
  int string_value{1};
  bool direct_leaf{};
};

struct BoundDfsStep {
  std::variant<BoundFactor, BoundDeterministicStep> operation;
  size_t scratch_offset{};
  size_t output_capacity{};
};

struct BoundOperatorPlan {
  KernelKind kind{KernelKind::Legacy};
  BitMaskKernel bit_mask;
  std::vector<value_type> bit_scales;
  std::vector<BoundDeterministicStep> deterministic;
  std::vector<BoundFactor> factors;
  std::vector<BoundDfsStep> dfs;
  size_t scratch_size{};
};

struct DfsFrame {
  idxv_type raw{};
  value_type amplitude{};
  size_t next_child{};
  size_t child_count{};
  bool expanded{};
};

inline constexpr size_t recursive_dfs_depth_limit = 32;

struct OperatorWorkspace {
  std::vector<RawTransition> scratch;
  std::vector<DfsFrame> frames;

  explicit OperatorWorkspace(const BoundOperatorPlan &plan)
      : scratch(plan.scratch_size),
        frames(plan.kind == KernelKind::Dfs &&
                       plan.dfs.size() > recursive_dfs_depth_limit
                   ? plan.dfs.size()
                   : 0) {}
};

inline size_t checked_size_add(size_t lhs, size_t rhs,
                               const char *message) {
  if (rhs > std::numeric_limits<size_t>::max() - lhs)
    throw std::length_error(message);
  return lhs + rhs;
}

inline bool matrix_descriptor_valid(const FactorDescriptor &factor) noexcept {
  if (factor.input_range == 0 || factor.output_range == 0) return false;
  if (factor.input_range >
      std::numeric_limits<size_t>::max() / factor.output_range)
    return false;
  return factor.matrix.size() ==
         factor.input_range * factor.output_range;
}

inline bool finite_value(value_type value) noexcept {
  return std::isfinite(value.real()) && std::isfinite(value.imag());
}

inline bool factor_is_deterministic(
    const FactorDescriptor &factor) noexcept {
  if (factor.kind == FactorKind::HopBit) return true;
  if (!matrix_descriptor_valid(factor)) return false;
  for (const auto entry : factor.matrix)
    if (!finite_value(entry)) return false;
  for (size_t input = 0; input < factor.input_range; ++input) {
    size_t live = 0;
    for (size_t output = 0; output < factor.output_range; ++output) {
      if (factor.matrix_element(input, output) != value_type{} && ++live > 1)
        return false;
    }
  }
  return true;
}

inline bool try_bit_matrix_mask(const FactorDescriptor &factor,
                                idxv_type active_mask,
                                BitMaskKernel &kernel) {
  if (factor.kind != FactorKind::BitMatrix ||
      factor.input_range != 2 || factor.output_range != 2 ||
      !matrix_descriptor_valid(factor) ||
      (factor.string_value != 1 && factor.string_value != -1))
    return false;
  for (const auto entry : factor.matrix)
    if (!finite_value(entry)) return false;

  const auto site_mask = checked_site_mask(factor.site, active_mask);
  std::array<idxv_type, 2> output{2, 2};
  std::array<value_type, 2> amplitude{};
  size_t active_inputs = 0;
  for (size_t input = 0; input < 2; ++input) {
    size_t live = 0;
    for (size_t out = 0; out < 2; ++out) {
      const auto entry = factor.matrix_element(input, out);
      if (entry == value_type{}) continue;
      if (++live > 1) return false;
      output[input] = out;
      amplitude[input] = entry;
    }
    if (live != 0) ++active_inputs;
  }

  if (active_inputs == 0) {
    kernel = zero_bit_mask_kernel(active_mask);
    return true;
  }

  size_t first_input = output[0] != 2 ? 0 : 1;
  const idxv_type flip_value =
      static_cast<idxv_type>(first_input) ^ output[first_input];
  for (size_t input = 0; input < 2; ++input) {
    if (output[input] == 2) continue;
    if ((static_cast<idxv_type>(input) ^ output[input]) != flip_value)
      return false;
  }

  kernel = BitMaskKernel{};
  kernel.active_mask = active_mask;
  kernel.scale = amplitude[first_input];
  if (flip_value != 0) kernel.flip = site_mask;
  if (factor.string_value == -1)
    kernel.parity = higher_site_mask(factor.site, active_mask);

  if (active_inputs == 1) {
    if (first_input == 0)
      kernel.must_clear = site_mask;
    else
      kernel.must_set = site_mask;
  } else if (amplitude[1] == amplitude[0]) {
    // no input-dependent local sign
  } else if (amplitude[1] == -amplitude[0]) {
    kernel.parity ^= site_mask;
  } else {
    return false;
  }

  kernel = canonicalize_bit_mask_kernel(std::move(kernel));
  return true;
}

inline bool try_hop_bit_mask(const FactorDescriptor &factor,
                             idxv_type active_mask,
                             BitMaskKernel &kernel) {
  if (factor.kind != FactorKind::HopBit ||
      (factor.string_value != 1 && factor.string_value != -1))
    return false;

  const auto source = checked_site_mask(factor.site, active_mask);
  const auto destination = checked_site_mask(factor.site2, active_mask);
  kernel = BitMaskKernel{};
  kernel.active_mask = active_mask;
  kernel.must_set = source;
  if (factor.site == factor.site2) {
    kernel.scale = value_type(1.0);
    return true;
  }

  kernel.must_clear = destination;
  kernel.flip = source | destination;
  if (factor.string_value == -1) {
    kernel.parity = strictly_between_site_mask(
        factor.site, factor.site2, active_mask);
  }
  kernel.scale = value_type(1.0);
  return true;
}

inline bool try_factor_bit_mask(const FactorDescriptor &factor,
                                idxv_type active_mask,
                                BitMaskKernel &kernel) {
  if (factor.kind == FactorKind::BitMatrix)
    return try_bit_matrix_mask(factor, active_mask, kernel);
  if (factor.kind == FactorKind::HopBit)
    return try_hop_bit_mask(factor, active_mask, kernel);
  return false;
}

inline BoundDeterministicStep make_deterministic_step(
    const BoundFactor &factor) {
  BoundDeterministicStep result;
  result.binding = factor.binding;
  result.kind = factor.descriptor.kind;
  result.site2 = factor.descriptor.site2;
  result.string_value = factor.descriptor.string_value;
  result.direct_leaf = factor.direct_leaf;
  result.stage.primitive_identity = factor.descriptor.primitive_identity;

  if (factor.descriptor.kind == FactorKind::HopBit) {
    result.stage.target = DeterministicTarget::QbitSite;
    result.stage.coordinate = factor.descriptor.site;
    return result;
  }

  result.stage.target =
      factor.descriptor.kind == FactorKind::BitMatrix
          ? DeterministicTarget::QbitSite
          : DeterministicTarget::Mode;
  result.stage.coordinate = factor.descriptor.site;
  result.stage.table.input_range = factor.descriptor.input_range;
  result.stage.table.output_range = factor.descriptor.output_range;
  result.stage.table.output.assign(factor.descriptor.input_range,
                                   factor.descriptor.output_range);
  result.stage.table.amplitude.assign(factor.descriptor.input_range,
                                      value_type{});
  for (size_t input = 0; input < factor.descriptor.input_range; ++input) {
    for (size_t output = 0; output < factor.descriptor.output_range;
         ++output) {
      const auto entry = factor.descriptor.matrix_element(input, output);
      if (entry == value_type{}) continue;
      result.stage.table.output[input] = output;
      result.stage.table.amplitude[input] = entry;
      break;
    }
  }

  if (factor.descriptor.kind == FactorKind::BitMatrix &&
      factor.descriptor.string_value != 1) {
    const auto active_mask = factor.binding.range - 1;
    result.stage.phases.push_back(
        {ParityEvaluation::Input, active_mask,
         higher_site_mask(factor.descriptor.site, active_mask),
         value_type(factor.descriptor.string_value)});
  }
  return result;
}

template <typename Basis>
concept StatelessBasis = requires(
    const Basis &basis, RawLeafBindings &bindings, idxv_type raw) {
  { basis.raw_at(raw) } -> std::convertible_to<idxv_type>;
  { basis.rank(raw) } -> std::convertible_to<idxv_type>;
  append_raw_leaf_bindings(basis, idxv_type{1}, bindings);
};

template <typename Basis>
BoundOperatorPlan bind_operator_plan(const ops_type &ops,
                                     const Basis &basis) {
  BoundOperatorPlan plan;
  if constexpr (!StatelessBasis<Basis>) {
    return plan;
  } else {
    const auto bindings = collect_raw_leaf_bindings(basis);
    plan.factors.reserve(ops.size());
    for (auto op = ops.rbegin(); op != ops.rend(); ++op) {
      FactorDescriptor descriptor;
      if (!(*op)->describe_factor(descriptor)) return BoundOperatorPlan{};
      const auto *binding = find_unique_raw_leaf_binding(
          bindings, descriptor.primitive_identity);
      if (binding == nullptr) return BoundOperatorPlan{};

      switch (descriptor.kind) {
      case FactorKind::ModeMatrix:
        if (binding->kind != RawLeafKind::SingleMode ||
            !matrix_descriptor_valid(descriptor) ||
            descriptor.input_range != binding->range ||
            descriptor.output_range != binding->range)
          return BoundOperatorPlan{};
        break;
      case FactorKind::BitMatrix: {
        if (binding->kind != RawLeafKind::Qbits ||
            !matrix_descriptor_valid(descriptor) ||
            descriptor.input_range != 2 ||
            descriptor.output_range != 2)
          return BoundOperatorPlan{};
        const auto active_mask = binding->range - 1;
        try {
          (void)checked_site_mask(descriptor.site, active_mask);
        } catch (const std::out_of_range &) {
          return BoundOperatorPlan{};
        }
        break;
      }
      case FactorKind::HopBit: {
        if (binding->kind != RawLeafKind::Qbits ||
            descriptor.input_range != 2 ||
            descriptor.output_range != 1 ||
            !descriptor.matrix.empty())
          return BoundOperatorPlan{};
        const auto active_mask = binding->range - 1;
        try {
          (void)checked_site_mask(descriptor.site, active_mask);
          (void)checked_site_mask(descriptor.site2, active_mask);
        } catch (const std::out_of_range &) {
          return BoundOperatorPlan{};
        }
        break;
      }
      default:
        return BoundOperatorPlan{};
      }

      BoundFactor bound{descriptor, *binding, plan.scratch_size};
      bound.direct_leaf =
          bindings.size() == 1 && binding->stride == 1;
      plan.scratch_size = checked_size_add(
          plan.scratch_size, descriptor.output_range,
          "operator DFS scratch size overflow");
      plan.factors.push_back(std::move(bound));
    }

    bool all_deterministic = true;
    for (const auto &factor : plan.factors)
      all_deterministic &= factor_is_deterministic(factor.descriptor);

    bool mask_compatible = bindings.size() == 1 && !plan.factors.empty();
    BitMaskKernel combined;
    bool have_mask = false;
    if (mask_compatible) {
      const auto identity = plan.factors.front().descriptor.primitive_identity;
      const auto active_mask = bindings.front().range - 1;
      for (const auto &factor : plan.factors) {
        if (factor.descriptor.primitive_identity != identity ||
            factor.binding.stride != 1 ||
            (factor.descriptor.kind != FactorKind::BitMatrix &&
             factor.descriptor.kind != FactorKind::HopBit)) {
          mask_compatible = false;
          break;
        }
        BitMaskKernel local;
        if (!try_factor_bit_mask(factor.descriptor, active_mask, local)) {
          mask_compatible = false;
          break;
        }
        local = canonicalize_bit_mask_kernel(std::move(local));
        if (!local.identically_zero) {
          plan.bit_scales.push_back(local.scale);
          local.scale = value_type(1.0);
        }
        combined = have_mask ? compose_after(local, combined) : local;
        have_mask = true;
      }
    }

    if (mask_compatible && have_mask) {
      plan.kind = KernelKind::BitMask;
      plan.bit_mask = std::move(combined);
      plan.scratch_size = 0;
      return plan;
    }

    if (all_deterministic) {
      plan.kind = KernelKind::Deterministic;
      plan.deterministic.reserve(plan.factors.size());
      for (const auto &factor : plan.factors) {
        plan.deterministic.push_back(make_deterministic_step(factor));
      }
      plan.scratch_size = 0;
      return plan;
    }

    plan.kind = KernelKind::Dfs;
    plan.dfs.reserve(plan.factors.size());
    plan.scratch_size = 0;
    for (const auto &factor : plan.factors) {
      if (factor_is_deterministic(factor.descriptor)) {
        auto deterministic = make_deterministic_step(factor);
        plan.dfs.push_back(
            {std::move(deterministic), plan.scratch_size, 1});
        plan.scratch_size = checked_size_add(
            plan.scratch_size, 1,
            "operator DFS scratch size overflow");
      } else {
        auto dfs_factor = factor;
        dfs_factor.scratch_offset = plan.scratch_size;
        plan.dfs.push_back(
            {std::move(dfs_factor), plan.scratch_size,
             factor.descriptor.output_range});
        plan.scratch_size = checked_size_add(
            plan.scratch_size, factor.descriptor.output_range,
            "operator DFS scratch size overflow");
      }
    }
    return plan;
  }
}

inline value_type bit_string_phase(idxv_type qbits, size_t site,
                                   int string_value) noexcept {
  const auto higher =
      site >= raw_index_bits - 1 ? idxv_type{} : (qbits >> (site + 1));
  return (std::popcount(higher) & 1U) != 0
             ? value_type(string_value)
             : value_type(1.0);
}

inline size_t emit_bound_factor(const BoundFactor &factor, idxv_type raw_in,
                                value_type amplitude_in,
                                std::span<RawTransition> output) {
  if (amplitude_in == value_type{}) return 0;
  const auto &descriptor = factor.descriptor;
  const auto local =
      factor.direct_leaf ? raw_in : factor.binding.extract(raw_in);
  const auto encode_output = [&](idxv_type local_out) {
    return factor.direct_leaf
               ? local_out
               : factor.binding.replace(raw_in, local, local_out);
  };
  size_t count = 0;

  if (descriptor.kind == FactorKind::BitMatrix) {
    if (output.size() < descriptor.output_range) return 0;
    const auto input_bit = static_cast<size_t>((local >> descriptor.site) & 1);
    const auto phase =
        descriptor.string_value == 1
            ? value_type(1.0)
            : bit_string_phase(local, descriptor.site,
                               descriptor.string_value);
    const auto site_mask = idxv_type{1} << descriptor.site;
    for (size_t out = 0; out < descriptor.output_range; ++out) {
      const auto entry = descriptor.matrix_element(input_bit, out);
      const auto amplitude = phase * entry * amplitude_in;
      if (amplitude == value_type{}) continue;
      const auto local_out =
          (local & ~site_mask) | (static_cast<idxv_type>(out) << descriptor.site);
      output[count++] = {encode_output(local_out), amplitude};
    }
    return count;
  }

  if (descriptor.kind == FactorKind::ModeMatrix) {
    if (output.size() < descriptor.output_range) return 0;
    for (size_t out = 0; out < descriptor.output_range; ++out) {
      const auto entry = descriptor.matrix_element(local, out);
      const auto amplitude = entry * amplitude_in;
      if (amplitude == value_type{}) continue;
      output[count++] = {encode_output(out), amplitude};
    }
    return count;
  }

  if (descriptor.kind != FactorKind::HopBit || output.empty()) return 0;
  const auto source_mask = idxv_type{1} << descriptor.site;
  const auto destination_mask = idxv_type{1} << descriptor.site2;
  if (descriptor.site == descriptor.site2) {
    if ((local & source_mask) == 0) return 0;
    output[0] = {raw_in, amplitude_in};
    return 1;
  }
  if ((local & source_mask) == 0 || (local & destination_mask) != 0)
    return 0;

  auto local_out = local & ~source_mask;
  value_type phase(1.0);
  if (descriptor.string_value != 1) {
    phase = bit_string_phase(
        local, descriptor.site, descriptor.string_value);
    phase *= bit_string_phase(
        local_out, descriptor.site2, descriptor.string_value);
  }
  local_out |= destination_mask;
  const auto amplitude = phase * amplitude_in;
  if (amplitude == value_type{}) return 0;
  output[0] = {encode_output(local_out), amplitude};
  return 1;
}

inline bool apply_deterministic_step(const BoundDeterministicStep &step,
                                     RawTransition &state) {
  if (step.kind == FactorKind::HopBit) {
    FactorDescriptor descriptor;
    descriptor.kind = FactorKind::HopBit;
    descriptor.site = step.stage.coordinate;
    descriptor.site2 = step.site2;
    descriptor.string_value = step.string_value;
    BoundFactor bound{descriptor, step.binding, 0};
    bound.direct_leaf = step.direct_leaf;
    RawTransition output;
    if (emit_bound_factor(bound, state.raw, state.amplitude,
                          std::span<RawTransition>(&output, 1)) == 0)
      return false;
    state = output;
    return true;
  }

  const auto local =
      step.direct_leaf ? state.raw : step.binding.extract(state.raw);
  idxv_type table_input = local;
  idxv_type local_output = local;
  if (step.stage.target == DeterministicTarget::QbitSite) {
    table_input = (local >> step.stage.coordinate) & 1;
  }
  const auto &table = step.stage.table;
  if (table_input >= table.output.size()) return false;
  const auto table_output = table.output[table_input];
  if (table_output >= table.output_range) return false;

  if (step.stage.target == DeterministicTarget::QbitSite) {
    const auto site_mask = idxv_type{1} << step.stage.coordinate;
    local_output =
        (local & ~site_mask) | (table_output << step.stage.coordinate);
  } else {
    local_output = table_output;
  }

  value_type phase(1.0);
  for (const auto &ordered : step.stage.phases) {
    const auto phase_label =
        ordered.evaluation == ParityEvaluation::Input ? local : local_output;
    if (parity_count(phase_label, ordered.parity_mask) != 0)
      phase *= ordered.odd_multiplier;
  }

  state.amplitude *= phase * table.amplitude[table_input];
  if (state.amplitude == value_type{}) return false;
  state.raw = step.direct_leaf
                  ? local_output
                  : step.binding.replace(state.raw, local, local_output);
  return true;
}

template <typename Basis, typename Sink>
void emit_ranked_leaf(const Basis &basis, idxv_type source_compact,
                      idxv_type source_raw, const RawTransition &leaf,
                      Sink &&sink) {
  if (!(std::abs(leaf.amplitude) > eps)) return;
  const auto compact =
      leaf.raw == source_raw ? source_compact : basis.rank(leaf.raw);
  if (compact < basis.range()) sink(compact, leaf.amplitude);
}

inline size_t emit_dfs_step(const BoundDfsStep &step,
                            const RawTransition &input,
                            OperatorWorkspace &workspace) {
  auto outputs = std::span<RawTransition>(
      workspace.scratch.data() + step.scratch_offset,
      step.output_capacity);
  if (std::holds_alternative<BoundDeterministicStep>(
          step.operation)) {
    RawTransition state = input;
    if (!apply_deterministic_step(
            std::get<BoundDeterministicStep>(step.operation), state))
      return 0;
    outputs[0] = state;
    return 1;
  }
  return emit_bound_factor(
      std::get<BoundFactor>(step.operation), input.raw,
      input.amplitude, outputs);
}

template <typename Basis, typename Sink>
void evaluate_dfs_recursive(
    const BoundOperatorPlan &plan, const Basis &basis,
    idxv_type source_compact, idxv_type source_raw, size_t depth,
    const RawTransition &input, OperatorWorkspace &workspace, Sink &sink) {
  const auto &step = plan.dfs[depth];
  const auto child_count = emit_dfs_step(step, input, workspace);
  const auto child_begin = workspace.scratch.data() + step.scratch_offset;
  if (depth + 1 == plan.dfs.size()) {
    for (size_t child = 0; child < child_count; ++child)
      emit_ranked_leaf(
          basis, source_compact, source_raw, child_begin[child], sink);
    return;
  }
  for (size_t child = 0; child < child_count; ++child)
    evaluate_dfs_recursive(
        plan, basis, source_compact, source_raw, depth + 1,
        child_begin[child], workspace, sink);
}

template <typename Basis, typename Sink>
void evaluate_source(const BoundOperatorPlan &plan, const Basis &basis,
                     idxv_type source_compact, idxv_type source_raw,
                     value_type initial, OperatorWorkspace &workspace,
                     Sink &&sink) {
  if (initial == value_type{}) return;

  if (plan.kind == KernelKind::BitMask) {
    if (!bit_mask_accepts(plan.bit_mask, source_raw)) return;
    value_type amplitude = initial;
    for (const auto scale : plan.bit_scales) {
      amplitude *= scale;
      if (amplitude == value_type{}) return;
    }
    amplitude *= plan.bit_mask.scale;
    if (parity_count(source_raw, plan.bit_mask.parity) != 0)
      amplitude = -amplitude;
    if (amplitude == value_type{}) return;
    const RawTransition output{
        source_raw ^ plan.bit_mask.flip, amplitude};
    emit_ranked_leaf(basis, source_compact, source_raw, output,
                     std::forward<Sink>(sink));
    return;
  }

  if (plan.kind == KernelKind::Deterministic) {
    RawTransition state{source_raw, initial};
    for (const auto &step : plan.deterministic)
      if (!apply_deterministic_step(step, state)) return;
    emit_ranked_leaf(basis, source_compact, source_raw, state,
                     std::forward<Sink>(sink));
    return;
  }

  if (plan.kind != KernelKind::Dfs || plan.dfs.empty()) return;
  if (plan.dfs.size() <= recursive_dfs_depth_limit) {
    const RawTransition root{source_raw, initial};
    evaluate_dfs_recursive(
        plan, basis, source_compact, source_raw, 0, root, workspace, sink);
    return;
  }
  auto &root = workspace.frames[0];
  root = {source_raw, initial, 0, 0, false};
  size_t depth = 0;
  while (true) {
    auto &frame = workspace.frames[depth];
    const auto &step = plan.dfs[depth];
    if (!frame.expanded) {
      frame.child_count = emit_dfs_step(
          step, {frame.raw, frame.amplitude}, workspace);
      frame.next_child = 0;
      frame.expanded = true;
    }

    if (frame.next_child < frame.child_count) {
      const auto child =
          workspace.scratch[step.scratch_offset + frame.next_child++];
      if (depth + 1 == plan.dfs.size()) {
        emit_ranked_leaf(basis, source_compact, source_raw, child, sink);
      } else {
        ++depth;
        workspace.frames[depth] = {
            child.raw, child.amplitude, 0, 0, false};
      }
      continue;
    }

    frame.expanded = false;
    if (depth == 0) break;
    --depth;
  }
}

}  // namespace detail
//=========================================================================
template <typename OP, typename STATE>
class OPsiType {
  friend STATE;
  const OP &op;
  const STATE &psi;

 public:
  class isOPsi {};
  OPsiType(const OP &op, const STATE &psi) : op(op), psi(psi) {}
};
//=========================================================================
class Operator {
  friend class OperatorSum;

  ops_type ops_;

  template <typename trans_index_type, typename index_type>
  inline size_t legacy_full_branch(
      trans_index_type &trans_index, index_type &index, value_type init,
      detail::LegacyTreeWorkspace &workspace, idx_it_type &idx_b,
      val_it_type &val_b) const;

  template <typename index_type>
  void append_triplets(index_type &index, value_type output_scale,
                       triplets_type &triplets) const;

 public:
  explicit Operator(const std::shared_ptr<elOpBase> &op) : ops_{op} {}

  Operator operator*(const Operator &Op) const {
    Operator Opp(*this);
    Opp.ops_.insert(Opp.ops_.end(), Op.ops_.begin(), Op.ops_.end());
    return Opp;
  }

  template <typename index_type>
  void map_acc(State<index_type> &to_psi, int t1,
               const State<index_type> &from_psi, int t2,
               value_type = 1.0) const;

  template <typename index_type>
  void map(State<index_type> &to_psi, int t1,
           const State<index_type> &from_psi, int t2,
           value_type = 1.0) const;

  template <typename index_type>
  auto operator*(const State<index_type> &psi) const {
    return OPsiType<Operator, State<index_type>>(*this, psi);
  }

  template <typename index_type,
            typename = typename std::remove_reference_t<index_type>::is_index>
  SpMatrix operator>>(index_type &&index) const;

  template <
      typename index_type, typename trans_index_type,
      typename = typename std::remove_reference_t<index_type>::is_index,
      typename = typename std::remove_reference_t<trans_index_type>::is_index>
  SpMatrix operator()(index_type &&index, trans_index_type &&trans_index) const {
    return this->operator>>(index);
  }

  bool operator==(const Operator &op) const {
    if (ops_.size() != op.ops_.size()) return false;

    for (auto i : range(ops_.size())) {
      if (ops_[i] != op.ops_[i]) return false;
    }

    return true;
  }

#ifdef QUDRIP_TESTING
  template <typename index_type>
  detail::KernelKind selected_kernel_for_test(const index_type &index) const {
    return detail::bind_operator_plan(ops_, index).kind;
  }
#endif
};

//==========================================================================
template <typename trans_index_type, typename index_type>
inline size_t Operator::legacy_full_branch(
    trans_index_type &trans_index, index_type &index, value_type init,
    detail::LegacyTreeWorkspace &workspace, idx_it_type &idx_b,
    val_it_type &val_b) const {
  workspace.indices[0] = idx_size_t(trans_index) + 1;
  workspace.amplitudes[0] = init;

  auto idx_front_b = workspace.indices.begin();
  auto idx_front_e = workspace.indices.begin() + 1;
  auto val_front_b = workspace.amplitudes.begin();
  auto val_front_e = workspace.amplitudes.begin() + 1;

  for (auto op_it = ops_.rbegin(); op_it != ops_.rend(); ++op_it) {
    const auto &op = *op_it;
    const auto idx_next_b = idx_front_e;
    const auto val_next_b = val_front_e;
    auto idx_emit = idx_next_b;
    auto val_emit = val_next_b;
    auto idx_compact = idx_next_b;
    auto val_compact = val_next_b;

    auto idx_read = idx_front_b;
    auto val_read = val_front_b;
    for (; idx_read != idx_front_e; ++idx_read, ++val_read) {
      if (*idx_read == 0 || *val_read == value_type{}) continue;

      trans_index[*idx_read - 1];
      *idx_read = op->get_idx() + 1;

      const auto idx_emitted_b = idx_emit;
      const auto val_emitted_b = val_emit;
      auto idx_input = idx_read;
      auto val_input = val_read;
      op->branch(idx_input, val_input, idx_emit, val_emit);

      auto child_idx = idx_emitted_b;
      auto child_val = val_emitted_b;
      for (; child_idx != idx_emit; ++child_idx, ++child_val) {
        if (*child_idx == 0 || *child_val == value_type{}) continue;
        op->feed_idx(*child_idx - 1);
        *idx_compact = idxv_type(trans_index) + 1;
        *val_compact = *child_val;
        ++idx_compact;
        ++val_compact;
      }
    }

    idx_front_b = idx_next_b;
    idx_front_e = idx_compact;
    val_front_b = val_next_b;
    val_front_e = val_compact;
    if (idx_front_b == idx_front_e) break;
  }

  auto idx_write = idx_front_b;
  auto val_write = val_front_b;
  auto idx_read = idx_front_b;
  auto val_read = val_front_b;
  for (; idx_read != idx_front_e; ++idx_read, ++val_read) {
    if (*idx_read == 0 || !(std::abs(*val_read) > eps)) continue;
    trans_index[*idx_read - 1];
    const auto compact_index = idxv_type(index);
    if (compact_index >= index.range()) continue;
    *idx_write = compact_index + 1;
    *val_write = *val_read;
    ++idx_write;
    ++val_write;
  }

  idx_b = idx_front_b;
  val_b = val_front_b;
  return static_cast<size_t>(idx_write - idx_front_b);
}
//==========================================================================
template <typename index_type>
void Operator::map(State<index_type> &to_psi, int t1,
                   const State<index_type> &from_psi, int t2,
                   value_type coeff) const {
  to_psi(t1).mat().setZero();
  map_acc(to_psi, t1, from_psi, t2, coeff);
}
//==========================================================================
template <typename index_type>
void Operator::map_acc(State<index_type> &to_psi, int t1,
                       const State<index_type> &from_psi, int t2,
                       value_type coeff) const {
  auto &index = from_psi.get_index();
  const auto plan = detail::bind_operator_plan(ops_, index);
  if (plan.kind != detail::KernelKind::Legacy) {
    detail::OperatorWorkspace workspace(plan);
    auto output = to_psi(t1).mat();
    const auto ndim = index.range();
    for (idxv_type i = 0; i < ndim; ++i) {
      const auto raw = index.raw_at(i);
      detail::evaluate_source(
          plan, index, i, raw, coeff * from_psi.data()(i, t2), workspace,
          [&](idxv_type compact, value_type amplitude) {
            output(compact, 0) += amplitude;
          });
    }
    return;
  }

  if constexpr (std::is_const_v<index_type>) {
    throw std::invalid_argument(
        "legacy elementary operators require a mutable index");
  } else {
    auto &mutable_index = const_cast<index_type &>(index);
    auto &trans_index = strip(mutable_index);
    auto ndim = mutable_index.range();
    detail::LegacyTreeWorkspace workspace(ops_);
    idx_it_type idx_b;
    val_it_type val_b;
    for (idxv_type i = 0; i < ndim; ++i) {
      mutable_index[i];
      const auto live_count =
          legacy_full_branch(trans_index, mutable_index,
                             coeff * from_psi.data()(i, t2), workspace, idx_b,
                             val_b);
      for (size_t j = 0; j < live_count; ++j) {
        to_psi(t1).mat()(*(idx_b + j) - 1, 0) += *(val_b + j);
      }
    }
  }
}

//==========================================================================
template <typename index_type>
void Operator::append_triplets(index_type &index, value_type output_scale,
                               triplets_type &triplets) const {
  using storage_index_type = SpMatrix::StorageIndex;
  (void)detail::checked_eigen_sparse_index(
      index.range(), "operator dimension exceeds Eigen sparse index range");
  const auto plan = detail::bind_operator_plan(ops_, index);
  if (plan.kind != detail::KernelKind::Legacy) {
    detail::OperatorWorkspace workspace(plan);
    const auto ndim = index.range();
    for (idxv_type i = 0; i < ndim; ++i) {
      const auto raw = index.raw_at(i);
      detail::evaluate_source(
          plan, index, i, raw, value_type(1.0), workspace,
          [&](idxv_type compact, value_type amplitude) {
            triplets.emplace_back(
                static_cast<storage_index_type>(compact),
                static_cast<storage_index_type>(i),
                output_scale * amplitude);
          });
    }
    return;
  }

  if constexpr (std::is_const_v<std::remove_reference_t<index_type>>) {
    throw std::invalid_argument(
        "legacy elementary operators require a mutable index");
  } else {
    auto &trans_index = strip(index);
    auto ndim = index.range();
    detail::LegacyTreeWorkspace workspace(ops_);
    idx_it_type idx_b;
    val_it_type val_b;

#ifdef PRINT_TREE
    std::cout << "tree expansion..." << std::endl;
#endif
    for (idxv_type i = 0; i < ndim; ++i) {
      index[i];
      const auto live_count = legacy_full_branch(
          trans_index, index, 1.0, workspace, idx_b, val_b);

      for (size_t j = 0; j < live_count; ++j) {
        triplets.emplace_back(
            static_cast<storage_index_type>((*(idx_b + j)) - 1),
            static_cast<storage_index_type>(i),
            output_scale * *(val_b + j));
      }

#ifdef PRINT_TREE
      std::cout << std::endl;
#endif
    }
  }
}

template <typename index_type, typename>
SpMatrix Operator::operator>>(index_type &&index) const {
  const auto ndim = index.range();
  const auto matrix_dimension = detail::checked_eigen_sparse_index(
      ndim, "operator dimension exceeds Eigen sparse index range");
  SpMatrix mat(matrix_dimension, matrix_dimension);
  triplets_type triplets;

  if (ndim <= triplets.max_size())
    triplets.reserve(static_cast<size_t>(ndim));
  append_triplets(index, value_type(1.0), triplets);

  if (triplets.size() > detail::eigen_sparse_storage_limit())
    throw std::length_error(
        "operator triplet count exceeds Eigen sparse storage range");
  mat.setFromTriplets(triplets.begin(), triplets.end());

  return mat;
}

inline size_t count_bits(const idxv_type &bits) {
  return bits_type(bits).count();
}

}  // namespace qudrip

#include "operators/bit.hpp"
#include "operators/hopbit.hpp"
#include "operators/mode.hpp"
