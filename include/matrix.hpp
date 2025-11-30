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

auto MatrixExp(const Matrix& mat) { return mat.exp(); }

}  // namespace qudrip

#else
namespace qudrip {
class Matrix {};

class SpMatrix {};
}  // namespace qudrip
#endif
