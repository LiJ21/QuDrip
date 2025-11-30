#pragma once

#ifdef USE_EIGEN

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>

namespace qudrip {
using Matrix = Eigen::MatrixXcd;
using SpMatrix = Eigen::SparseMatrix<value_type>;
using t_type = Eigen::Triplet<value_type>;
using triplets_type = std::vector<t_type>;
template <typename T, int N>
using StaticVector = Eigen::Matrix<T, N, 1>;
template <int N>
using StaticIVector = StaticVector<int, N>;
template <int N>
using StaticDVector = StaticVector<double, N>;
auto MatrixExp(const Matrix& mat) { return mat.exp(); }

}  // namespace qudrip

#else
namespace qudrip {
class Matrix {};

class SpMatrix {};
}  // namespace qudrip
#endif
