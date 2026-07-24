#include "qudrip.hpp"

#include <unsupported/Eigen/KroneckerProduct>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

using namespace qudrip;

namespace {

int operator_checks = 0;

#define OPERATOR_CHECK(cond)                                                   \
  do {                                                                         \
    ++operator_checks;                                                         \
    if (!(cond)) {                                                             \
      std::cerr << "OPERATOR FAILED (line " << __LINE__ << "): " #cond         \
                << "\n";                                                       \
      std::exit(1);                                                            \
    }                                                                          \
  } while (0)

Matrix dense(const SpMatrix &matrix) { return Matrix(matrix); }

Matrix scaled_identity(size_t dimension, value_type scale) {
  return scale * Matrix::Identity(dimension, dimension);
}

template <typename RAW> auto make_spin_hamiltonian(RAW &raw, int sites) {
  auto hop = getHopBoseGate(raw);
  auto hamiltonian = getEmptyOperator();
  for (int site = 0; site + 1 < sites; ++site) {
    hamiltonian = hamiltonian + hop(site + 1, site) + hop(site, site + 1);
  }
  return hamiltonian;
}

// Deliberately implements only the original mutable elementary-operator API.
// In particular, it inherits elOpBase::describe_factor(), so it must remain on
// the compatibility kernel.
class LegacyOnlyBitOperator final : public elOpBase {
  QbitsIndex<> &index_;
  size_t site_;
  size_t branch_calls_ = 0;

public:
  LegacyOnlyBitOperator(QbitsIndex<> &index, size_t site)
      : elOpBase(2, 2), index_(index), site_(site) {}

  void feed_idx(idx_size_t idx) override { index_[idx]; }
  idx_size_t get_idx() override { return index_; }

  void branch(idx_it_type &idx_b, val_it_type &val_b, idx_it_type &idx_w,
              val_it_type &val_w) override {
    ++branch_calls_;
    if (*idx_b == 0) {
      for (size_t output = 0; output < 2; ++output) {
        *idx_w++ = 0;
        *val_w++ = value_type{};
      }
      return;
    }

    bits_type bits(*idx_b - 1);
    const auto input = static_cast<size_t>(bits[site_]);
    for (size_t output = 0; output < 2; ++output) {
      bits[site_] = output != 0;
      *idx_w++ = bit_convert<idxv_type>(bits) + 1;
      *val_w++ = U_[input * 2 + output] * *val_b;
    }
  }

  size_t branch_calls() const noexcept { return branch_calls_; }
};

class MalformedDescriptorOperator final : public elOpBase {
  QbitsIndex<> &index_;

public:
  explicit MalformedDescriptorOperator(QbitsIndex<> &index)
      : elOpBase(1, 1), index_(index) {}

  void feed_idx(idx_size_t idx) override { index_[idx]; }
  idx_size_t get_idx() override { return index_; }

  void branch(idx_it_type &idx_b, val_it_type &val_b, idx_it_type &idx_w,
              val_it_type &val_w) override {
    *idx_w++ = *idx_b;
    *val_w++ = *val_b;
  }

  bool describe_factor(detail::FactorDescriptor &factor) const
      noexcept override {
    factor.kind = static_cast<detail::FactorKind>(255);
    factor.primitive_identity = std::addressof(index_);
    factor.input_range = 1;
    factor.output_range = 0;
    factor.matrix = {};
    factor.site = 0;
    factor.site2 = 0;
    return true;
  }
};

} // namespace

int run_operator_tests() {
  {
    bool overflow_rejected = false;
    try {
      (void)detail::checked_factor_matrix_size(
          std::numeric_limits<size_t>::max(), size_t(2));
    } catch (const std::length_error &) {
      overflow_rejected = true;
    }
    OPERATOR_CHECK(overflow_rejected);
  }

  // Arbitrary local matrices retain their supplied orientation. Products are
  // evaluated right-to-left: A * B applies B first and therefore represents
  // the matrix product A B.
  {
    auto index = getQbits(1);
    auto gate_a = getBoseGate(index);
    auto gate_b = getBoseGate(index);
    Matrix matrix_a(2, 2);
    Matrix matrix_b(2, 2);
    matrix_a << value_type(1.0, 0.5), value_type(-2.0, 1.0),
        value_type(0.25, -0.75), value_type(3.0, 0.0);
    matrix_b << value_type(-0.5, 0.25), value_type(1.5, -0.5),
        value_type(2.0, 1.0), value_type(-1.0, 0.75);
    gate_a << matrix_a;
    gate_b << matrix_b;

    auto op_a = gate_a(0);
    auto op_b = gate_b(0);
    OPERATOR_CHECK(dense(op_a >> index).isApprox(matrix_a, 1e-13));
    OPERATOR_CHECK(
        dense((op_a * op_b) >> index).isApprox(matrix_a * matrix_b, 1e-13));
  }

  // An external subclass that implements only the original stateful virtual
  // interface stays on Legacy, while matrix construction and lazy application
  // retain the same numerical behavior.
  {
    auto index = getQbits(1);
    auto legacy = std::make_shared<LegacyOnlyBitOperator>(index, 0);
    Matrix local(2, 2);
    local << value_type(0.5, -0.25), value_type(1.25, 0.5),
        value_type(-0.75, 0.25), value_type(2.0, -0.5);
    *legacy << local;
    Operator op(std::static_pointer_cast<elOpBase>(legacy));

    OPERATOR_CHECK(op.selected_kernel_for_test(index) ==
                   detail::KernelKind::Legacy);
    const SpMatrix assembled = op >> index;
    OPERATOR_CHECK(dense(assembled).isApprox(local, 1e-13));
    OPERATOR_CHECK(legacy->branch_calls() > 0);

    auto builtin_gate = getBoseGate(index);
    builtin_gate << id;
    const SpMatrix mixed_sum = (op + builtin_gate(0)) >> index;
    OPERATOR_CHECK(
        dense(mixed_sum).isApprox(local + Matrix::Identity(2, 2), 1e-13));

    auto source = getState(index, 1);
    auto destination = getState(index, 1);
    source.data()(0, 0) = value_type(0.75, -0.5);
    source.data()(1, 0) = value_type(-1.25, 0.25);
    const Matrix before = source.data().col(0);
    const auto calls_before = legacy->branch_calls();
    destination(0) = op * source(0);
    OPERATOR_CHECK(destination.data().col(0).isApprox(local * before, 1e-13));
    OPERATOR_CHECK(source.data().col(0).isApprox(before, 0.0));
    OPERATOR_CHECK(legacy->branch_calls() > calls_before);

    const auto &const_index = index;
    bool rejected_const_legacy = false;
    try {
      (void)(op >> const_index);
    } catch (const std::invalid_argument &) {
      rejected_const_legacy = true;
    }
    OPERATOR_CHECK(rejected_const_legacy);
  }

  // The compact legacy frontier never expands an exactly dead parent, while
  // a nonzero sub-eps parent remains eligible for later amplification.
  {
    auto index = getQbits(1);
    auto zero = std::make_shared<LegacyOnlyBitOperator>(index, 0);
    auto amplify = std::make_shared<LegacyOnlyBitOperator>(index, 0);
    *zero << Matrix::Zero(2, 2);
    *amplify << scaled_identity(2, 4.0);
    const Operator zero_op(std::static_pointer_cast<elOpBase>(zero));
    const Operator amplify_op(std::static_pointer_cast<elOpBase>(amplify));

    const SpMatrix dead = (amplify_op * zero_op) >> index;
    OPERATOR_CHECK(dead.nonZeros() == 0);
    OPERATOR_CHECK(zero->branch_calls() == index.range());
    OPERATOR_CHECK(amplify->branch_calls() == 0);

    auto small = std::make_shared<LegacyOnlyBitOperator>(index, 0);
    *small << scaled_identity(2, eps / 2.0);
    const Operator small_op(std::static_pointer_cast<elOpBase>(small));
    const SpMatrix revived = (amplify_op * small_op) >> index;
    OPERATOR_CHECK(revived.nonZeros() == 2);
    OPERATOR_CHECK(dense(revived).isApprox(
        scaled_identity(2, 2.0 * eps), 1e-20));
  }

  // Malformed opt-in metadata cannot enter a stateless emitter with an
  // invalid scratch span; it is contained by the compatibility path.
  {
    auto index = getQbits(1);
    auto malformed = std::make_shared<MalformedDescriptorOperator>(index);
    const Operator op(std::static_pointer_cast<elOpBase>(malformed));
    OPERATOR_CHECK(op.selected_kernel_for_test(index) ==
                   detail::KernelKind::Legacy);
    OPERATOR_CHECK(
        dense(op >> index).isApprox(Matrix::Identity(2, 2), 0.0));
  }

  // Fermionic elementary gates retain the current Jordan-Wigner convention
  // and canonical anticommutation relations.
  {
    auto index = getQbits(3);
    auto annihilate = getFermiGate(index);
    auto create = getFermiGate(index);
    Matrix a(2, 2);
    Matrix adag(2, 2);
    a << 0.0, 1.0, 0.0, 0.0;
    adag << 0.0, 0.0, 1.0, 0.0;
    annihilate << a;
    create << adag;
    const Matrix identity = Matrix::Identity(index.range(), index.range());

    for (int lhs = 0; lhs < 3; ++lhs) {
      for (int rhs = 0; rhs < 3; ++rhs) {
        const Matrix c_lhs = dense(annihilate(lhs) >> index);
        const Matrix cd_rhs = dense(create(rhs) >> index);
        const Matrix anti =
            c_lhs * cd_rhs +
            dense(create(rhs) >> index) * dense(annihilate(lhs) >> index);
        if (lhs == rhs)
          OPERATOR_CHECK(anti.isApprox(identity, 1e-13));
        else
          OPERATOR_CHECK(anti.norm() < 1e-13);
      }
    }
  }

  // Hopping arguments are (destination, source), and site zero is the low
  // bit. Equal-site hopping is the corresponding number projector.
  {
    auto index = getQbits(2);
    auto hop = getHopBoseGate(index);
    Matrix expected_10 = Matrix::Zero(4, 4);
    Matrix expected_01 = Matrix::Zero(4, 4);
    Matrix expected_n0 = Matrix::Zero(4, 4);
    expected_10(2, 1) = 1.0;
    expected_01(1, 2) = 1.0;
    expected_n0(1, 1) = 1.0;
    expected_n0(3, 3) = 1.0;

    OPERATOR_CHECK(dense(hop(1, 0) >> index).isApprox(expected_10, 1e-13));
    OPERATOR_CHECK(dense(hop(0, 1) >> index).isApprox(expected_01, 1e-13));
    OPERATOR_CHECK(dense(hop(0, 0) >> index).isApprox(expected_n0, 1e-13));
  }

  // The public custom string multiplier is promoted before the two
  // Jordan-Wigner factors are multiplied on the legacy compatibility API.
  {
    auto index = getQbits(4);
    hopBitOp<idx_size_t> custom(index, std::numeric_limits<int>::max());
    (void)custom(2, 0);
    idx_tree_type input_indices{idxv_type(9) + 1};
    val_tree_type input_values{value_type(1.0)};
    idx_tree_type output_indices(1);
    val_tree_type output_values(1);
    auto idx_b = input_indices.begin();
    auto val_b = input_values.begin();
    auto idx_w = output_indices.begin();
    auto val_w = output_values.begin();
    custom.branch(idx_b, val_b, idx_w, val_w);
    const double string_value =
        static_cast<double>(std::numeric_limits<int>::max());
    OPERATOR_CHECK(output_indices[0] == idxv_type(12) + 1);
    OPERATOR_CHECK(
        std::abs(output_values[0] -
                 value_type(string_value * string_value)) /
            (string_value * string_value) <
        1e-13);
  }

  // Intermediate pruning is exact-zero only. A sub-eps intermediate can be
  // amplified above eps, while a final value exactly equal to eps is absent.
  {
    auto index = getQbits(1);
    auto small_gate = getBoseGate(index);
    auto amplify_gate = getBoseGate(index);
    auto zero_gate = getBoseGate(index);
    small_gate << scaled_identity(2, eps / 2.0);
    amplify_gate << scaled_identity(2, 4.0);
    zero_gate << Matrix::Zero(2, 2);

    const SpMatrix amplified = (amplify_gate(0) * small_gate(0)) >> index;
    OPERATOR_CHECK(amplified.nonZeros() == 2);
    OPERATOR_CHECK(
        dense(amplified).isApprox(scaled_identity(2, 2.0 * eps), 1e-13));

    const SpMatrix exact_zero = (amplify_gate(0) * zero_gate(0)) >> index;
    OPERATOR_CHECK(exact_zero.nonZeros() == 0);

    auto boundary_gate = getBoseGate(index);
    boundary_gate << scaled_identity(2, eps);
    const SpMatrix boundary = boundary_gate(0) >> index;
    OPERATOR_CHECK(boundary.nonZeros() == 0);

    auto accepted_gate = getBoseGate(index);
    accepted_gate << scaled_identity(2, 2.0 * eps);
    const SpMatrix accepted = accepted_gate(0) >> index;
    OPERATOR_CHECK(accepted.nonZeros() == 2);
  }

  // The final positive predicate rejects NaN and accepts infinity.
  {
    auto index = getQbits(1);
    auto nan_gate = getBoseGate(index);
    Matrix nan_matrix = Matrix::Zero(2, 2);
    nan_matrix(0, 0) =
        value_type(std::numeric_limits<double>::quiet_NaN(), 0.0);
    nan_gate << nan_matrix;
    const auto nan_op = nan_gate(0);
    OPERATOR_CHECK(nan_op.selected_kernel_for_test(index) ==
                   detail::KernelKind::Dfs);
    const SpMatrix nan_result = nan_op >> index;
    OPERATOR_CHECK(nan_result.nonZeros() == 0);

    auto inf_gate = getBoseGate(index);
    Matrix inf_matrix = Matrix::Zero(2, 2);
    inf_matrix(0, 0) = value_type(std::numeric_limits<double>::infinity(), 0.0);
    inf_gate << inf_matrix;
    const auto inf_op = inf_gate(0);
    OPERATOR_CHECK(inf_op.selected_kernel_for_test(index) ==
                   detail::KernelKind::Dfs);
    const SpMatrix inf_result = inf_op >> index;
    OPERATOR_CHECK(inf_result.nonZeros() == 1);
    OPERATOR_CHECK(std::isinf(inf_result.coeff(0, 0).real()));
  }

  // Matrix assembly thresholds an unscaled term, then applies its sum
  // coefficient. Lazy application includes the coefficient in its seed.
  {
    auto index = getQbits(1);
    auto small_gate = getBoseGate(index);
    small_gate << scaled_identity(2, eps / 2.0);
    auto small_op = small_gate(0);
    auto boosted = value_type(4.0) * small_op;

    const SpMatrix boosted_matrix = boosted >> index;
    OPERATOR_CHECK(boosted_matrix.nonZeros() == 0);

    auto boosted_source = getState(index, 1);
    auto boosted_destination = getState(index, 1);
    boosted_source.data()(0, 0) = 1.0;
    boosted_destination(0) = boosted * boosted_source(0);
    OPERATOR_CHECK(std::abs(boosted_destination.data()(0, 0) - 2.0 * eps) <
                   1e-22);

    auto large_gate = getBoseGate(index);
    large_gate << scaled_identity(2, 2.0 * eps);
    auto large_op = large_gate(0);
    auto damped = value_type(0.1) * large_op;

    const SpMatrix damped_matrix = damped >> index;
    OPERATOR_CHECK(damped_matrix.nonZeros() == 2);
    OPERATOR_CHECK(std::abs(damped_matrix.coeff(0, 0) - 0.2 * eps) < 1e-22);

    auto damped_source = getState(index, 1);
    auto damped_destination = getState(index, 1);
    damped_source.data()(0, 0) = 1.0;
    damped_destination(0) = damped * damped_source(0);
    OPERATOR_CHECK(damped_destination.data().col(0).norm() == 0.0);
  }

  // Compiled products retain chronological multiplication with the lazy
  // source seed. Pre-multiplying these two finite scales would overflow even
  // though the original staged evaluation has a finite result.
  {
    auto index = getQbits(1);
    auto first_gate = getBoseGate(index);
    auto second_gate = getBoseGate(index);
    first_gate << scaled_identity(2, 1e200);
    second_gate << scaled_identity(2, 1e200);
    const auto product = second_gate(0) * first_gate(0);
    OPERATOR_CHECK(product.selected_kernel_for_test(index) ==
                   detail::KernelKind::BitMask);

    auto source = getState(index, 1);
    auto destination = getState(index, 1);
    source.data()(0, 0) = value_type(1e-300);
    destination(0) = product * source(0);
    OPERATOR_CHECK(
        std::isfinite(destination.data()(0, 0).real()));
    OPERATOR_CHECK(
        std::abs(destination.data()(0, 0) - value_type(1e100)) /
            1e100 <
        1e-13);

    auto mode = getSingleMode(2);
    auto first_mode_gate = getModeOp(mode);
    auto second_mode_gate = getModeOp(mode);
    first_mode_gate << scaled_identity(2, 1e200);
    second_mode_gate << scaled_identity(2, 1e200);
    const auto mode_product = second_mode_gate(0) * first_mode_gate(0);
    OPERATOR_CHECK(mode_product.selected_kernel_for_test(mode) ==
                   detail::KernelKind::Deterministic);
    auto mode_source = getState(mode, 1);
    auto mode_destination = getState(mode, 1);
    mode_source.data()(0, 0) = value_type(1e-300);
    mode_destination(0) = mode_product * mode_source(0);
    OPERATOR_CHECK(
        std::isfinite(mode_destination.data()(0, 0).real()));
    OPERATOR_CHECK(
        std::abs(mode_destination.data()(0, 0) - value_type(1e100)) /
            1e100 <
        1e-13);

  }

  // Distinct terms that reach the same coordinate add. Current Eigen sparse
  // accumulation retains an explicit stored zero after exact cancellation.
  {
    auto index = getQbits(1);
    auto gate = getBoseGate(index);
    gate << id;
    auto first = gate(0);
    auto second = gate(0);

    const SpMatrix doubled = (first + second) >> index;
    OPERATOR_CHECK(doubled.nonZeros() == 2);
    OPERATOR_CHECK(
        dense(doubled).isApprox(2.0 * Matrix::Identity(2, 2), 1e-13));

    const auto cancelling_sum = first - second;
    const SpMatrix cancelled = cancelling_sum >> index;
    const SpMatrix cancelled_fused =
        cancelling_sum.fused_triplets_for_test(index);
    OPERATOR_CHECK(cancelled.nonZeros() == 2);
    OPERATOR_CHECK(cancelled_fused.nonZeros() == 2);
    OPERATOR_CHECK(dense(cancelled).norm() == 0.0);
    OPERATOR_CHECK(
        dense(cancelled).isApprox(dense(cancelled_fused), 0.0));

    const value_type complex_scale(0.25, -0.75);
    const SpMatrix complex_scaled = (complex_scale * first) >> index;
    OPERATOR_CHECK(
        dense(complex_scaled)
            .isApprox(complex_scale * Matrix::Identity(2, 2), 1e-13));

    const SpMatrix empty = getEmptyOperator() >> index;
    OPERATOR_CHECK(empty.rows() == 2 && empty.cols() == 2);
    OPERATOR_CHECK(empty.nonZeros() == 0);
  }

  // Dense mode matrices retain orientation and duplicate paths within a
  // product reduce to ordinary matrix multiplication.
  {
    auto mode = getSingleMode(3);
    auto gate_a = getModeOp(mode);
    auto gate_b = getModeOp(mode);
    Matrix matrix_a(3, 3);
    Matrix matrix_b(3, 3);
    matrix_a << value_type(1.0, 0.0), value_type(0.5, 0.25),
        value_type(-1.0, 0.5), value_type(0.25, -0.5), value_type(2.0, 0.0),
        value_type(0.75, 0.0), value_type(-0.5, 0.0), value_type(1.0, -0.25),
        value_type(0.0, 0.5);
    matrix_b << value_type(0.5, 0.0), value_type(1.0, -0.5),
        value_type(0.0, 0.25), value_type(-0.75, 0.0), value_type(0.25, 0.5),
        value_type(1.5, 0.0), value_type(1.0, 0.0), value_type(-0.5, 0.25),
        value_type(0.75, -0.5);
    gate_a << matrix_a;
    gate_b << matrix_b;

    auto op_a = gate_a(0);
    auto op_b = gate_b(0);
    OPERATOR_CHECK(dense(op_a >> mode).isApprox(matrix_a, 1e-13));
    OPERATOR_CHECK(
        dense((op_a * op_b) >> mode).isApprox(matrix_a * matrix_b, 1e-13));
  }

  // A four-state cyclic mode map is deterministic even though it is neither
  // a ladder nor an XOR permutation. Both matrix and lazy stateless paths
  // leave the mode cursor untouched.
  {
    auto mode = getSingleMode(4);
    auto gate = getModeOp(mode);
    Matrix cyclic = Matrix::Zero(4, 4);
    for (idxv_type input = 0; input < mode.range(); ++input) {
      cyclic((input + 1) % mode.range(), input) =
          value_type(double(input + 1), -0.125 * double(input));
    }
    gate << cyclic;
    auto op = gate(0);

    OPERATOR_CHECK(op.selected_kernel_for_test(mode) ==
                   detail::KernelKind::Deterministic);
    mode[2];
    const SpMatrix assembled = op >> mode;
    OPERATOR_CHECK(dense(assembled).isApprox(cyclic, 1e-13));
    OPERATOR_CHECK(idxv_type(mode) == 2);

    auto source = getState(mode, 1);
    auto destination = getState(mode, 1);
    for (idxv_type input = 0; input < mode.range(); ++input)
      source.data()(input, 0) =
          value_type(0.5 * double(input + 1), -0.2 * double(input));
    const Matrix before = source.data().col(0);
    mode[3];
    destination(0) = op * source(0);
    OPERATOR_CHECK(destination.data().col(0).isApprox(cyclic * before, 1e-13));
    OPERATOR_CHECK(source.data().col(0).isApprox(before, 0.0));
    OPERATOR_CHECK(idxv_type(mode) == 3);
  }

  // Unequal tensor-product modes preserve the existing mixed-radix order.
  {
    auto mode0 = getSingleMode(2);
    auto mode1 = getSingleMode(3);
    auto modes = mode0 * mode1;
    auto create0 = getModeOp(mode0);
    auto destroy1 = getModeOp(mode1);
    const Matrix upper = LadderMatrix(mode0.range(), Mode::Upper);
    const Matrix lower = LadderMatrix(mode1.range(), Mode::Lower);
    create0 << upper;
    destroy1 << lower;

    const Matrix actual = dense((create0(0) * destroy1(0)) >> modes);
    const Matrix expected = Eigen::kroneckerProduct(upper, lower).eval();
    OPERATOR_CHECK(actual.isApprox(expected, 1e-13));
  }

  // A bit factor embedded in a mixed-radix product cannot use a global XOR.
  // It binds to the Qbits digit and remains stateless through the deterministic
  // kernel; matrix and lazy evaluation preserve both primitive cursors.
  {
    auto bits = getQbits(2);
    auto mode = getSingleMode(3);
    auto mixed = bits * mode;
    auto gate = getBoseGate(bits);
    gate << pauli_x;
    auto op = gate(1);

    const auto selected = op.selected_kernel_for_test(mixed);
    OPERATOR_CHECK(selected != detail::KernelKind::BitMask);
    OPERATOR_CHECK(selected == detail::KernelKind::Deterministic);

    const Matrix bit_matrix = dense(op >> bits);
    const Matrix expected =
        Eigen::kroneckerProduct(bit_matrix,
                                Matrix::Identity(mode.range(), mode.range()))
            .eval();

    bits[2];
    mode[1];
    const SpMatrix assembled = op >> mixed;
    OPERATOR_CHECK(dense(assembled).isApprox(expected, 1e-13));
    OPERATOR_CHECK(idxv_type(bits) == 2);
    OPERATOR_CHECK(idxv_type(mode) == 1);

    auto source = getState(mixed, 1);
    auto destination = getState(mixed, 1);
    for (idxv_type input = 0; input < mixed.range(); ++input)
      source.data()(input, 0) =
          value_type(double(input + 1), -0.05 * double(input));
    const Matrix before = source.data().col(0);
    bits[3];
    mode[2];
    destination(0) = op * source(0);
    OPERATOR_CHECK(
        destination.data().col(0).isApprox(expected * before, 1e-12));
    OPERATOR_CHECK(source.data().col(0).isApprox(before, 0.0));
    OPERATOR_CHECK(idxv_type(bits) == 3);
    OPERATOR_CHECK(idxv_type(mode) == 2);
  }

  // Deterministic stages on both sides of one genuinely branching factor are
  // retained inside the DFS plan and applied in product order.
  {
    auto index = getQbits(1);
    auto left_gate = getBoseGate(index);
    auto dense_gate = getBoseGate(index);
    auto right_gate = getBoseGate(index);
    Matrix dense_factor(2, 2);
    dense_factor << value_type(1.0, 0.25), value_type(-0.5, 0.75),
        value_type(0.25, -1.0), value_type(1.5, 0.5);
    left_gate << pauli_x;
    dense_gate << dense_factor;
    right_gate << pauli_z;

    auto op = left_gate(0) * dense_gate(0) * right_gate(0);
    const Matrix expected = pauli_x * dense_factor * pauli_z;
    OPERATOR_CHECK(op.selected_kernel_for_test(index) ==
                   detail::KernelKind::Dfs);

    index[1];
    const SpMatrix assembled = op >> index;
    OPERATOR_CHECK(dense(assembled).isApprox(expected, 1e-13));
    OPERATOR_CHECK(idxv_type(index) == 1);

    auto source = getState(index, 1);
    auto destination = getState(index, 1);
    source.data()(0, 0) = value_type(0.5, -0.25);
    source.data()(1, 0) = value_type(-1.0, 0.75);
    const Matrix before = source.data().col(0);
    index[0];
    destination(0) = op * source(0);
    OPERATOR_CHECK(
        destination.data().col(0).isApprox(expected * before, 1e-13));
    OPERATOR_CHECK(source.data().col(0).isApprox(before, 0.0));
    OPERATOR_CHECK(idxv_type(index) == 0);
  }

  // Hash and direct-DP cluster sectors produce the same physical operator.
  // A single ladder that leaves a fixed-number sector is discarded.
  {
    constexpr int sites = 6;

    auto hash_raw = getQbits(sites);
    auto hash_number = getNset(hash_raw);
    auto hash_layout = getSpinClusters(hash_raw);
    auto hash_sector = Constrain(hash_layout, hash_number = sites / 2);

    auto direct_raw = getQbits(sites);
    auto direct_number = getNset(direct_raw);
    auto direct_layout = getSpinClusters(direct_raw);
    direct_layout.useDirectRank();
    auto direct_sector = Constrain(direct_layout, direct_number = sites / 2);

    OPERATOR_CHECK(!hash_sector.has_direct_rank());
    OPERATOR_CHECK(direct_sector.has_direct_rank());
    OPERATOR_CHECK(hash_sector.range() == direct_sector.range());

    auto hash_hop_gate = getHopBoseGate(hash_raw);
    auto direct_hop_gate = getHopBoseGate(direct_raw);
    auto hash_probe = hash_hop_gate(1, 0);
    auto direct_probe = direct_hop_gate(1, 0);
    OPERATOR_CHECK(hash_probe.selected_kernel_for_test(hash_sector) ==
                   detail::KernelKind::BitMask);
    OPERATOR_CHECK(direct_probe.selected_kernel_for_test(direct_sector) ==
                   detail::KernelKind::BitMask);

    hash_sector[4];
    direct_sector[4];
    const auto hash_raw_before_matrix = idxv_type(hash_raw);
    const auto direct_raw_before_matrix = idxv_type(direct_raw);
    const SpMatrix hash_probe_matrix = hash_probe >> hash_sector;
    const SpMatrix direct_probe_matrix = direct_probe >> direct_sector;
    OPERATOR_CHECK(
        dense(hash_probe_matrix).isApprox(dense(direct_probe_matrix), 1e-13));
    OPERATOR_CHECK(idxv_type(hash_sector) == 4);
    OPERATOR_CHECK(idxv_type(direct_sector) == 4);
    OPERATOR_CHECK(idxv_type(hash_raw) == hash_raw_before_matrix);
    OPERATOR_CHECK(idxv_type(direct_raw) == direct_raw_before_matrix);

    auto hash_hamiltonian = make_spin_hamiltonian(hash_raw, sites);
    auto direct_hamiltonian = make_spin_hamiltonian(direct_raw, sites);
    const SpMatrix hash_matrix = hash_hamiltonian >> hash_sector;
    const SpMatrix direct_matrix = direct_hamiltonian >> direct_sector;
    const SpMatrix hash_fused =
        hash_hamiltonian.fused_triplets_for_test(hash_sector);
    const SpMatrix direct_fused =
        direct_hamiltonian.fused_triplets_for_test(direct_sector);
    OPERATOR_CHECK(dense(hash_matrix).isApprox(dense(direct_matrix), 1e-13));
    OPERATOR_CHECK(dense(hash_matrix).isApprox(dense(hash_fused), 1e-13));
    OPERATOR_CHECK(
        dense(direct_matrix).isApprox(dense(direct_fused), 1e-13));
    OPERATOR_CHECK(idxv_type(hash_sector) == 4);
    OPERATOR_CHECK(idxv_type(direct_sector) == 4);
    OPERATOR_CHECK(idxv_type(hash_raw) == hash_raw_before_matrix);
    OPERATOR_CHECK(idxv_type(direct_raw) == direct_raw_before_matrix);

    Matrix create_matrix(2, 2);
    create_matrix << 0.0, 0.0, 1.0, 0.0;
    auto create = getBoseGate(hash_raw);
    create << create_matrix;
    OPERATOR_CHECK((create(0) >> hash_sector).nonZeros() == 0);

    auto hash_source = getState(hash_sector, 1);
    auto hash_destination = getState(hash_sector, 1);
    auto direct_source = getState(direct_sector, 1);
    auto direct_destination = getState(direct_sector, 1);
    for (idxv_type i = 0; i < hash_sector.range(); ++i) {
      const value_type value(double(i + 1), -0.125 * double(i));
      hash_source.data()(i, 0) = value;
      direct_source.data()(i, 0) = value;
    }
    const Matrix hash_before = hash_source.data().col(0);
    const Matrix direct_before = direct_source.data().col(0);
    hash_sector[7];
    direct_sector[7];
    const auto hash_raw_before_lazy = idxv_type(hash_raw);
    const auto direct_raw_before_lazy = idxv_type(direct_raw);
    hash_destination(0) = hash_hamiltonian * hash_source(0);
    direct_destination(0) = direct_hamiltonian * direct_source(0);

    OPERATOR_CHECK(hash_destination.data().col(0).isApprox(
        hash_matrix * hash_before, 1e-12));
    OPERATOR_CHECK(direct_destination.data().col(0).isApprox(
        direct_matrix * direct_before, 1e-12));
    OPERATOR_CHECK(hash_destination.data().col(0).isApprox(
        direct_destination.data().col(0), 1e-12));
    OPERATOR_CHECK(hash_source.data().col(0).isApprox(hash_before, 0.0));
    OPERATOR_CHECK(direct_source.data().col(0).isApprox(direct_before, 0.0));
    OPERATOR_CHECK(idxv_type(hash_sector) == 7);
    OPERATOR_CHECK(idxv_type(direct_sector) == 7);
    OPERATOR_CHECK(idxv_type(hash_raw) == hash_raw_before_lazy);
    OPERATOR_CHECK(idxv_type(direct_raw) == direct_raw_before_lazy);
  }

  // Exact structural recognition is based on the local matrix, not on a
  // named gate. Every single-entry matrix is a mask; unequal deterministic
  // amplitudes and many-to-one maps use tables; a multi-output column uses
  // DFS.
  {
    auto index = getQbits(1);
    std::vector<std::pair<Matrix, detail::KernelKind>> cases;
    cases.emplace_back(Matrix::Zero(2, 2), detail::KernelKind::BitMask);
    for (Eigen::Index row = 0; row < 2; ++row) {
      for (Eigen::Index column = 0; column < 2; ++column) {
        Matrix single = Matrix::Zero(2, 2);
        single(row, column) = value_type(1.25, -0.5);
        cases.emplace_back(std::move(single), detail::KernelKind::BitMask);
      }
    }

    Matrix unequal_diagonal = Matrix::Zero(2, 2);
    unequal_diagonal(0, 0) = value_type(2.0, 0.25);
    unequal_diagonal(1, 1) = value_type(-0.75, 0.5);
    cases.emplace_back(unequal_diagonal,
                       detail::KernelKind::Deterministic);

    Matrix unequal_antidiagonal = Matrix::Zero(2, 2);
    unequal_antidiagonal(1, 0) = value_type(0.5, -0.25);
    unequal_antidiagonal(0, 1) = value_type(3.0, 0.75);
    cases.emplace_back(unequal_antidiagonal,
                       detail::KernelKind::Deterministic);

    Matrix reset = Matrix::Zero(2, 2);
    reset(0, 0) = 1.0;
    reset(0, 1) = value_type(-2.0, 0.5);
    cases.emplace_back(reset, detail::KernelKind::Deterministic);

    Matrix branching = Matrix::Zero(2, 2);
    branching(0, 0) = 1.0;
    branching(1, 0) = value_type(0.5, 0.25);
    branching(1, 1) = -1.0;
    cases.emplace_back(branching, detail::KernelKind::Dfs);

    for (const auto &[local, expected_kind] : cases) {
      auto gate = getBoseGate(index);
      gate << local;
      auto op = gate(0);
      OPERATOR_CHECK(op.selected_kernel_for_test(index) == expected_kind);
      OPERATOR_CHECK(dense(op >> index).isApprox(local, 1e-13));
    }
  }

  // Product construction is linear in factor count; it no longer allocates
  // the old product-width expansion arena.
  {
    auto index = getQbits(1);
    auto dense_gate = getBoseGate(index);
    Matrix dense_factor(2, 2);
    dense_factor << 1.0, 2.0, 3.0, 4.0;
    dense_gate << dense_factor;
    const auto factor = dense_gate(0);
    auto product = factor;
    for (size_t depth = 1; depth < 80; ++depth)
      product = product * factor;
    OPERATOR_CHECK(product.selected_kernel_for_test(index) ==
                   detail::KernelKind::Dfs);

    auto identity_gate = getBoseGate(index);
    identity_gate << id;
    const auto identity_factor = identity_gate(0);
    auto deep_sparse_product = factor;
    for (size_t depth = 1; depth < 40; ++depth)
      deep_sparse_product = deep_sparse_product * identity_factor;
    OPERATOR_CHECK(deep_sparse_product.selected_kernel_for_test(index) ==
                   detail::KernelKind::Dfs);
    OPERATOR_CHECK(
        dense(deep_sparse_product >> index).isApprox(dense_factor, 1e-13));
  }

  // Mask composition has an intentionally directional phase correction.
  // Compare it with literal sequential application over the complete small
  // raw domain, including a flip-before-parity case.
  {
    constexpr idxv_type active = 0b111;
    std::vector<detail::BitMaskKernel> kernels{
        {active, 0, 0, 0b001, 0, value_type(1.0), false},
        {active, 0b001, 0b010, 0b011, 0b110,
         value_type(0.5, 0.25), false},
        {active, 0, 0b100, 0b100, 0b011,
         value_type(-1.0, 0.5), false},
        {active, 0b010, 0, 0, 0b101, value_type(2.0), false},
        detail::zero_bit_mask_kernel(active)};

    for (const auto &before : kernels) {
      for (const auto &after : kernels) {
        const auto composed = detail::compose_after(after, before);
        for (idxv_type raw = 0; raw <= active; ++raw) {
          const value_type seed(0.75, -0.125);
          detail::RawTransition intermediate;
          detail::RawTransition sequential;
          detail::RawTransition direct;
          const bool before_live =
              detail::apply_bit_mask(before, raw, seed, intermediate);
          const bool sequential_live =
              before_live &&
              detail::apply_bit_mask(after, intermediate.raw,
                                     intermediate.amplitude, sequential);
          const bool direct_live =
              detail::apply_bit_mask(composed, raw, seed, direct);
          OPERATOR_CHECK(sequential_live == direct_live);
          if (sequential_live) {
            OPERATOR_CHECK(sequential.raw == direct.raw);
            OPERATOR_CHECK(
                std::abs(sequential.amplitude - direct.amplitude) < 1e-13);
          }
        }
      }
    }
  }

  // Dispatch follows the exact structural hierarchy: flat Qbits masks,
  // arbitrary one-output deterministic tables, then genuinely branching DFS.
  {
    auto bits = getQbits(3);
    auto bose = getBoseGate(bits);
    bose << pauli_x;
    auto mask_op = bose(1);
    OPERATOR_CHECK(mask_op.selected_kernel_for_test(bits) ==
                   detail::KernelKind::BitMask);

    bits[5];
    const SpMatrix mask_matrix = mask_op >> bits;
    OPERATOR_CHECK(idxv_type(bits) == 5);
    const auto &const_bits = bits;
    OPERATOR_CHECK(
        dense(mask_op >> const_bits).isApprox(dense(mask_matrix), 0.0));
    OPERATOR_CHECK(
        dense((value_type(2.0) * mask_op) >> const_bits)
            .isApprox(2.0 * dense(mask_matrix), 1e-13));
    auto mask_source = getState(bits, 1);
    auto mask_destination = getState(bits, 1);
    for (idxv_type input = 0; input < bits.range(); ++input)
      mask_source.data()(input, 0) =
          value_type(double(input + 1), 0.125 * double(input));
    const Matrix mask_before = mask_source.data().col(0);
    bits[6];
    mask_destination(0) = mask_op * mask_source(0);
    OPERATOR_CHECK(mask_destination.data().col(0).isApprox(
        mask_matrix * mask_before, 1e-13));
    OPERATOR_CHECK(mask_source.data().col(0).isApprox(mask_before, 0.0));
    OPERATOR_CHECK(idxv_type(bits) == 6);

    Matrix reset(2, 2);
    reset << 1.0, 1.0, 0.0, 0.0;
    auto reset_gate = getBoseGate(bits);
    reset_gate << reset;
    auto reset_op = reset_gate(0);
    OPERATOR_CHECK(reset_op.selected_kernel_for_test(bits) ==
                   detail::KernelKind::Deterministic);

    Matrix dense_bit(2, 2);
    dense_bit << 1.0, 2.0, 3.0, 4.0;
    auto dense_gate = getBoseGate(bits);
    dense_gate << dense_bit;
    auto dense_op = dense_gate(0);
    OPERATOR_CHECK(dense_op.selected_kernel_for_test(bits) ==
                   detail::KernelKind::Dfs);

    auto mode = getSingleMode(4);
    auto mode_gate = getModeOp(mode);
    mode_gate << LadderMatrix(4, Mode::Upper);
    auto ladder = mode_gate(0);
    OPERATOR_CHECK(ladder.selected_kernel_for_test(mode) ==
                   detail::KernelKind::Deterministic);

    auto mixed = bits * mode;
    OPERATOR_CHECK(mask_op.selected_kernel_for_test(mixed) ==
                   detail::KernelKind::Deterministic);

    bitOp<idx_size_t> custom_string(bits, 2);
    Matrix create(2, 2);
    create << 0.0, 0.0, 1.0, 0.0;
    custom_string << create;
    auto custom_op = custom_string(0);
    OPERATOR_CHECK(custom_op.selected_kernel_for_test(bits) ==
                   detail::KernelKind::Deterministic);

    bool bad_site_rejected = false;
    try {
      (void)bose(3);
    } catch (const std::out_of_range &) {
      bad_site_rejected = true;
    }
    OPERATOR_CHECK(bad_site_rejected);

    auto hop = getHopBoseGate(bits);
    bad_site_rejected = false;
    try {
      (void)hop(3, 0);
    } catch (const std::out_of_range &) {
      bad_site_rejected = true;
    }
    OPERATOR_CHECK(bad_site_rejected);
  }

  return operator_checks;
}
