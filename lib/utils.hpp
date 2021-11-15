#pragma once
#include <random>
#include <chrono>
namespace qudrip
{
template<class IDX1, class IDX2>
SpMatrix restrict(IDX1 && idx1, IDX2 && idx2, SpMatrix & mat)
{
	auto ndim = idx2.range();
	SpMatrix res(ndim, ndim);
	res.reserve(ndim);
	for(idx_size_t k = 0; k < mat.outerSize(); ++k)
		for(SpMatrix::InnerIterator it(mat, k); it; ++it)
		{
			idx1[it.row()];
			idxv_type pp = idx2;
			idx1[it.col()];
			idxv_type ppp = idx2;
			if(pp != -1 && ppp != -1)
				res.insert(pp, ppp) = it.value();
		}

	return res;

}

template<class TIDX, class IDX>
void fill_state(State<TIDX> & psi, const std::vector<value_type> & fill_vec, IDX && index, idxv_type base_idx = 0, int tstp = 0)
{
	TIDX & tot_index = psi.get_index();
	auto ncopies = tot_index.range() / index.range();
	auto counter = 0;
	for(tot_index[0]; tot_index < tot_index.range() && counter < ncopies; ++tot_index)
	{
		if(index == base_idx)
		{
			value_type base_val = psi[tstp];
			for(index[0]; index < index.range(); ++index)
			{
				psi[tstp] = base_val * fill_vec[index];
			}
			index[base_idx];
			++counter;
		}
		
	}

}

template<class IDX1, class IDX2>
std::vector<double> state_to_p(int tstp, State<IDX1> & psi, IDX1 & tot_index, const IDX2 & rho_index)
{
	std::vector<double> p(rho_index.range());
	tot_index[0];
	for(auto i : range(tot_index.range()))
	{
		p[rho_index] += (std::conj(psi[tstp]) * psi[tstp]).real();
		++tot_index;
	}
	return p;
}

Matrix herm_exp(const Matrix & A, value_type alpha)
{
	assert(A.rows() == A.cols());
	Matrix expA(A.rows(), A.cols());
	Matrix temp(A.rows(), A.cols());
	temp.setZero();

	Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(A);
	if(eigensolver.info() != Eigen::Success) 
	{
		std::cout << "herm_exp: failed to solve the matrix" << std::endl;
		abort();
	}

	for(int i = 0; i < A.rows(); ++i)
	{
		temp(i, i) = std::exp(alpha * eigensolver.eigenvalues()(i));
	}

	expA.noalias() = eigensolver.eigenvectors() * temp * eigensolver.eigenvectors().adjoint();

	return expA;

}

template<typename Ir_R, typename Io_R, typename Ir_L, typename Io_L, typename State_r, typename State_o>
void reduce(int tstp, State_r & rho_r, const State_o & rho_o, Ir_R && idx_rR, Io_R && idx_oR, 
		Ir_L && idx_rL, Io_L && idx_oL)
{
	for(idx_rR[0]; idx_rR < idx_rR.range(); ++idx_rR)
	{
		for(idx_rL[0]; idx_rL < idx_rL.range(); ++idx_rL)
		{
			rho_r[tstp] = 0.0;
			for(idx_oR[0], idx_oL[0]; idx_oR < idx_oR.range(); ++idx_oR, ++idx_oL)
			{
				rho_r[tstp] += rho_o[tstp];
			}
		}
	}
}

template<typename RHO_IDX, typename L_IDX, typename R_IDX>
void state_to_rho(L_IDX & l_idx, R_IDX & r_idx, State<RHO_IDX> & rho, const State<R_IDX> & psi, int tstp = 0)
{
	for(l_idx[0]; l_idx < l_idx.range(); ++l_idx)
	{
		r_idx = l_idx;
		auto rval = psi[tstp];
		for(r_idx[0]; r_idx < r_idx.range(); ++r_idx)
		{
			auto lval = psi[tstp];
			rho[tstp] = std::conj(rval) * lval;
		}
	}
}

template<typename T>
void shuffle_state(State<T> & psi)
{
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dice(-1.0, 1.0);
	
	auto & index = psi.get_index();
	for(index[0]; index < index.range(); ++index)
		psi[0] = value_type(dice(mt), dice(mt));

	psi(0) /= std::sqrt(Hermitian::norm(psi(0)));
}

template<typename T>
class sec_clock
{
	std::chrono::time_point<std::chrono::high_resolution_clock> t_;
	public:
	sec_clock():
		t_(std::chrono::high_resolution_clock::now())
	{}

	auto note()
	{
		t_ = std::chrono::high_resolution_clock::now();
		return t_;
	}

	auto release()
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		return std::chrono::duration_cast<T>(t1 - t_).count() * 1e-3;
	}
};

template<typename T = std::chrono::milliseconds>
auto get_clock()
{
	return sec_clock<T>();
}

} // namespace
