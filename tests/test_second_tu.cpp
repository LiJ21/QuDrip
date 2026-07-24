// Second translation unit including the full library header.
// Linking this together with test_core.cpp verifies that every definition
// in the headers is inline (no ODR / duplicate-symbol violations).
#include "qudrip.hpp"
#include "find_param.hpp"
#include "lattice.hpp"

using namespace qudrip;

Matrix other_tu_gate_matrix() {
  auto idx = getQbits(1);
  auto g = getBoseGate(idx);
  g << pauli_x;
  return Matrix(g(0) >> idx);
}
