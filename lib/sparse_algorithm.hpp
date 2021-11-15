#pragma once
#include<fstream>
namespace qudrip
{
namespace Arnoldi
{
template<typename T>
inline double norm(const T & psi)
{
	assert(psi.cols() == 1);
	return (psi.adjoint() * psi)(0,0).real();
}


template<typename T>
inline double norm(const State<T> & psi)
{
	return (psi.mat().adjoint() * psi.mat())(0,0).real();
}

template<bool herm = 0>
int expandBasis(const SpMatrix & H, Matrix & V, std::vector<value_type> & h, int start, int add_dim)
{
	double n = 10.;
	int idx = 0;
	for(idx = start + 1; idx <= start + add_dim; ++idx)
	{
		V.col(idx) = H * V.col(idx - 1);
		for(auto i = 0; i < idx; ++i)
		{
			value_type new_h = (V.col(i).adjoint() * V.col(idx))(0,0);
			V.col(idx) -= new_h * V.col(i);
			h.push_back(new_h);
		}
		n = std::sqrt(norm(V.col(idx)));
		V.col(idx) /= n;
		h.push_back(n);
		if(n < eps) break;
	}
	return idx - 1;
}

template<>
int expandBasis<1>(const SpMatrix & H, Matrix & V, std::vector<value_type> & h, int start, int add_dim)
{
	double n = 10.;
	int idx = 0;
	for(idx = start + 1; idx <= start + add_dim; ++idx)
	{

		int idx0 = idx % 3;
		int idxm1 = (idx - 1) % 3;
		int idxm2 = (idx - 2) % 3;
		V.col(idx0) = H * V.col(idxm1);

		for(auto i = 0; i < idx - 2; ++i) h.push_back(0.);
		if(idxm2 >= 0) 
		{
			value_type new_h = (V.col(idxm2).adjoint() * V.col(idx0))(0,0);
			V.col(idx0) -= new_h * V.col(idxm2);
			h.push_back(new_h);
		}

		value_type new_h = (V.col(idxm1).adjoint() * V.col(idx0))(0,0);
		V.col(idx0) -= new_h * V.col(idxm1);
		h.push_back(new_h);

		n = std::sqrt(norm(V.col(idx0)));
		V.col(idx0) /= n;

		h.push_back(n);
		if(n < eps) break;
	}

	return idx - 1;
}

template<bool herm = 0>
void initV0(Matrix & V, int ndim)
{
	V.resize(0, 0);
}

template<>
void initV0<1>(Matrix & V, int ndim)
{
	V.resize(ndim, 1);
}

template<bool herm = 0>
int get_d2(int kry_dim) {return kry_dim + 1;}

template<>
int get_d2<1>(int kry_dim) {return 3;}

template<int herm>
struct herm_type
{};

template<typename OP_MAT, typename PSI_MAT, typename EV_TYPE>
void getState(Matrix & V, const OP_MAT & H, const EV_TYPE & evs, PSI_MAT && psi, const Matrix & V0, int s_idx, herm_type<0> dummy)
{
	int act_dim = evs.cols();
	int ndim = psi.rows();
	psi = V.block(0,0,ndim,act_dim) * evs.col(s_idx);
}

template<typename OP_MAT, typename PSI_MAT, typename EV_TYPE>
void getState(Matrix & V, const OP_MAT & H, const EV_TYPE & evs, PSI_MAT && psi, const Matrix & V0, int s_idx, herm_type<1> dummy)
{
	double n = 10.;
	int idx = 0;
	auto act_dim = evs.cols();
	V.col(0) = V0;
	std::vector<value_type> h(0);
	psi = evs(0, s_idx) * V.col(0);


	for(int i = 1; i < act_dim; ++i)
	{
		expandBasis<1>(H, V, h, i - 1, 1);
		psi += evs(i, s_idx) * V.col(i % 3);
	}
	expandBasis<1>(H, V, h, act_dim - 1, 1);
}

template<typename T, bool herm = 0>
class ArnoldiSolver
{
	using index_type = T;
	private:
		int kry_dim_;
		int act_dim_;
		int ndim_;
		Matrix KR, diagkry, V0;
		State<T> psi_V;
		Matrix & V;
		double err_=0;
		value_type n_last = 0.;
		Eigen::SelfAdjointEigenSolver<Matrix> eigensolver;
		Matrix evs_;

	public:

		double error(){return err_;}
		int act_dim(){return act_dim_;}
		void convergeEV(int N, value_type iv = 0.)
		{
			assert(N < kry_dim_);
			evs_.resize(N, 1);
			for(int i = 0; i < N; ++i)
			{
				evs_(i, 0) = iv;
			}
		}

		ArnoldiSolver(index_type & index, int kry_dim):
			kry_dim_(kry_dim),
			ndim_(index.range()),
			psi_V(index, get_d2(kry_dim)),
			V(psi_V.data()),
			KR(kry_dim_, kry_dim_),
			diagkry(kry_dim_, kry_dim_),
			act_dim_(0),
			evs_(1,1)
	{		
		initV0<herm>(V0, ndim_);
	}

		template <typename PSI_MAT>
		void initKrylov(const PSI_MAT & psi)
		{
			KR.setZero();
			V.col(0) = psi / std::sqrt(norm(psi));
			if(V0.rows() > 0) V0 = V.col(0);
		}

		template <typename OP_MAT>
		void updateKrylovRepr(const OP_MAT & H, int add_dim = -1)noexcept
		{

			if(add_dim == -1 || kry_dim_ < act_dim_ + add_dim) add_dim = kry_dim_ - act_dim_;
			std::vector<value_type> h;
			h.reserve(add_dim);
			double n = 10.0;

			int start = act_dim_;

			act_dim_ = expandBasis<herm>(H, V, h, start, add_dim);

			size_t idx = 0;
			if(act_dim_ > start && start > 0) KR(start, start - 1) = n_last;
			for(auto i = start; i < act_dim_ - 1; ++i)
			{
				for(auto j = 0; j <= i + 1; ++j)
				{
					KR(j,i) = h[idx];
					++idx;
				}

			}
			for(auto j = 0; j <= act_dim_ - 1; ++j)
			{
				KR(j,act_dim_ - 1) = h[idx];
				++idx;
			}
			n_last = h[idx];

		}

		template<typename OP_MAT>
		void evolve(const OP_MAT & H, State<index_type> & psi, double h)noexcept
		{
			int tstp = psi.t();
			assert(tstp < psi.nt() - 1);
			getKrylovRepr(H, psi);
			
			diagkry = MatrixExp(-II * KR.block(0, 0, act_dim_, act_dim_) * h);
			psi(tstp + 1).mat() = (V.block(0,0,ndim_,act_dim_) * diagkry).block(0, 0, ndim_, 1);
		}

		template<typename OP_MAT>
		void getGroundState(const OP_MAT & H, State<index_type> & psi, int add_dim)noexcept
		{
			int tstp = psi.t();
			if(act_dim_ >= kry_dim_) return;
			if(act_dim_ == 0) initKrylov(psi.mat());
			updateKrylovRepr(H, add_dim);
			eigensolver.compute(KR.block(0,0,act_dim_,act_dim_));
			if(eigensolver.info() != Eigen::Success) 
			{
				std::cout << "failed to solve the Krylov matrix" << std::endl;
				abort();
			}
			getState(V, H, eigensolver.eigenvectors(), psi.mat(), V0, 0, herm_type<herm>());
			err_ = std::sqrt(norm((eigensolver.eigenvalues().block(0,0,std::min((int)evs_.rows(),act_dim_),1)
						       	- evs_.block(0,0,std::min((int)evs_.rows(),act_dim_),1))));
			evs_.block(0,0,std::min((int)evs_.rows(),act_dim_),1) = 
				eigensolver.eigenvalues().block(0,0,std::min((int)evs_.rows(),act_dim_),1);
		}

		void storeE(std::string fname)
		{
			std::ofstream ef(fname);
                        for(auto i : range(eigensolver.eigenvalues().rows()))
                        {
                                ef << eigensolver.eigenvalues()(i, 0) << std::endl;
                        }
                        ef.close();
		}
};



template<class T>
ArnoldiSolver<T> getSolver(T & index, int kry_dim)
{
	return ArnoldiSolver<T>(index, kry_dim);
}

template<class T>
ArnoldiSolver<T, 1> getLanczosSolver(T & index, int kry_dim)
{
	return ArnoldiSolver<T, 1>(index, kry_dim);
}

}//namespace arnoldi

namespace Hermitian = Arnoldi;
}//namespace qudrip

