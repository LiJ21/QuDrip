#include <algorithm>
#include <complex>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>

#pragma once
#define USE_EIGEN

namespace qudrip {
using value_type = std::complex<double>;
using idx_size_t = unsigned long long;
using range_type = std::vector<size_t>;
using idxv_type = unsigned long long;
static const double eps = 1e-10;
constexpr short BIT_LIMIT = 64;
const value_type II = value_type(0.0, 1.0);
//--------------------------------------------------------------------
// range functions
range_type range(const size_t N) {
  range_type r(N);
  std::iota(r.begin(), r.end(), 0);
  return r;
}

range_type range(const size_t M, const size_t N) {
  range_type r(N - M);
  std::iota(r.begin(), r.end(), M);
  return r;
}

range_type range(const size_t M, const size_t N, const size_t STEP) {
  range_type r((N - M + 1) / STEP);
  auto s = M - STEP;
  std::generate(r.begin(), r.end(), [&s, STEP]() { return (s += STEP); });
  return r;
}

//--------------------------------------------------------------------

// Handle for left/right containers
template <typename VAL>
struct lrHandle {
  VAL val_;
  lrHandle(VAL val) : val_(val) {}
};

template <typename VAL>
struct lrHandle<VAL &> {
  VAL &val_;
  lrHandle(VAL &val) : val_(val) {}
};

}  // namespace qudrip

#include "index.hpp"
#include "conservation.hpp"
#include "matrix.hpp"
#include "state.hpp"
#include "operator.hpp"
#include "operator_abel.hpp"
#include "sparse_algorithm.hpp"
#include "utils.hpp"
#include "mc_manager.hpp"
