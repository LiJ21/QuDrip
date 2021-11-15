#pragma once

#ifdef USE_EIGEN

//#define EIGEN_USE_LAPACKE
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>

namespace qudrip
{
using Matrix = Eigen::MatrixXcd;
using SpMatrix = Eigen::SparseMatrix<value_type>;
using t_type = Eigen::Triplet<value_type>;
using triplets_type = std::vector<t_type>;
//using Block = Eigen::Block<Matrix, Eigen::Dynamic, Eigen::Dynamic>;

auto MatrixExp(const Matrix & mat)
{
	return mat.exp();
}

} // namespace

#elif
namespace qudrip
{
class Matrix
{};

class SpMatrix
{};
} // namespace
#endif
