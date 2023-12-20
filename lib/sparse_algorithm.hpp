#pragma once
#include<fstream>
namespace qudrip
{

namespace Hermitian
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
//====================================================================================================
template<typename OP_MAT, typename VEC_MAT, typename PSI_MAT>
double getKrylovReprAd(const OP_MAT & H, const PSI_MAT & psi, VEC_MAT & v, VEC_MAT & w, Matrix & KR) 
{
	assert(KR.cols() == KR.rows());
	auto kry_dim = KR.rows() - 1;
	v = psi / std::sqrt(norm(psi));
	std::vector<value_type> a,b;
	a.reserve(10);
	b.reserve(10);
	w.setZero();
	a.push_back((v.adjoint() * H * v)(0,0));
	w = H * v - a[0] * v;
	double n = norm(w);
	w /= std::sqrt(n);
	bool use_v = false;
	auto idx = 1;
	while(n > eps)
	{
		if(use_v)
		{
			a.push_back((v.adjoint() * H * v)(0,0));
			b.push_back(std::sqrt(n));

			w = H * v - a[a.size() - 1] * v - std::sqrt(n) * w;
			n = norm(w);
			w /= std::sqrt(n);

		}
		else
		{
			a.push_back((w.adjoint() * H * w)(0,0));
			b.push_back(std::sqrt(n));

			v = H * w - a[a.size() - 1] * w - std::sqrt(n) * v;
			n = norm(v);
			v /= std::sqrt(n);

		}
		use_v = !use_v;

		std::cout << "Lanczos idx, err = " << idx << ", " << std::sqrt(n) << std::endl;
		++idx;

		if(idx > kry_dim) break;
	}

	KR.setZero();
	KR(0,0) = a[0];
	for(auto i : range(1, a.size()))
	{
		KR(i,i) = a[i];
		KR(i-1,i) = b[i - 1];
		KR(i,i-1) = b[i - 1];
	}

	return std::sqrt(n);	
}

//====================================================================================================
template<typename OP_MAT, typename VEC_MAT, typename PSI_MAT>
Matrix getKrylovReprAd(const OP_MAT & H, const PSI_MAT & psi, VEC_MAT & v, VEC_MAT & w, size_t kry_dim, double & err) 
{
	v = psi / std::sqrt(norm(psi));
	std::vector<value_type> a,b;
	a.reserve(10);
	b.reserve(10);
	w.setZero();
	a.push_back((v.adjoint() * H * v)(0,0));
	w = H * v - a[0] * v;
	double n = norm(w);
	w /= std::sqrt(n);
	bool use_v = false;
	auto idx = 1;
	while(n > eps)
	{
		if(use_v)
		{
			a.push_back((v.adjoint() * H * v)(0,0));
			b.push_back(std::sqrt(n));

			w = H * v - a[a.size() - 1] * v - std::sqrt(n) * w;
			n = norm(w);
			w /= std::sqrt(n);

		}
		else
		{
			a.push_back((w.adjoint() * H * w)(0,0));
			b.push_back(std::sqrt(n));

			v = H * w - a[a.size() - 1] * w - std::sqrt(n) * v;
			n = norm(v);
			v /= std::sqrt(n);

		}

		std::cout << "Lanczos idx, err = " << idx << ", " << std::sqrt(n) << std::endl;
		use_v = !use_v;
		++idx;
		if(idx > kry_dim) break;
	}

	Matrix KR(a.size(), a.size());
	KR.setZero();
	KR(0,0) = a[0];
	for(auto i : range(1, a.size()))
	{
		KR(i,i) = a[i];
		KR(i-1,i) = b[i - 1];
		KR(i,i-1) = b[i - 1];
	}
	err = std::sqrt(n);

	return KR;	
}
//====================================================================================================
template<typename OP_MAT, typename PSI_MAT>
std::vector<Matrix> getKrylovRepr(const OP_MAT & H, const PSI_MAT & psi, size_t KryDim) 
{
	std::vector<Matrix> v;
	v.reserve(10);
	v.push_back(psi / std::sqrt(norm(psi)));

	std::vector<value_type> a,b;
	a.reserve(10);
	b.reserve(10);
	a.push_back((v[0].adjoint() * H * v[0])(0,0));

	v.push_back(H * v[0] - a[0] * v[0]);
	auto n = norm(v[1]);
	v[1] /= std::sqrt(n);

	auto idx = 1;
	while(n > eps)
	{
		a.push_back((v[idx].adjoint() * H * v[idx])(0,0));
		b.push_back(std::sqrt(n));

		v.push_back(H * v[idx] - a[idx] * v[idx] - b[idx - 1] * v[idx - 1]);
		n = norm(v[idx + 1]);
		v[idx + 1] /= std::sqrt(n);

		std::cout << "Lanczos idx, err = " << idx << ", " << std::sqrt(n) << std::endl;
		++idx;

		if(idx > KryDim) break;

	}
	auto kry_dim = a.size();
	Matrix KR(kry_dim, kry_dim);
	Matrix V(v[0].rows(), kry_dim);
	KR.setZero();
	KR(0,0) = a[0];
	V.col(0) = v[0];
	for(auto i = 1; i < kry_dim; ++i)
	{
		KR(i,i) = a[i];
		KR(i - 1, i) = b[i - 1];
		KR(i, i - 1) = KR(i - 1, i);

		V.col(i) = v[i];
	}


		
	return std::vector<Matrix>{KR, V};	
}
//====================================================================================================

template<typename OP_MAT, typename VEC_MAT, typename PSI_MAT, typename RES_MAT>
void getKrylovGS(const OP_MAT & H, const PSI_MAT & psi, const Matrix & GS, VEC_MAT & v, VEC_MAT & w, RES_MAT && psi_gs) 
{
	auto KryDim = GS.rows();
	v = psi / std::sqrt(norm(psi));
	w.setZero();
	auto alpha = 0.0;

	psi_gs = GS(0,0) * v;
	alpha  = (v.adjoint() * H * v)(0,0).real();
	w = H * v - alpha * v;
	double n = norm(w);
	w /= std::sqrt(n);
	for(auto idx = 1; idx < KryDim; ++idx)
	{
		if(!(idx % 2))
		{
			psi_gs += GS(idx, 0) * v;
			alpha  = (v.adjoint() * H * v)(0,0).real();

			w = H * v - alpha * v - std::sqrt(n) * w;
			n = norm(w);
			w /= std::sqrt(n);
		}
		else
		{
			psi_gs += GS(idx, 0) * w;
			alpha = (w.adjoint() * H * w)(0,0).real();

			v = H * w - alpha * w - std::sqrt(n) * v;
			n = norm(v);
			v /= std::sqrt(n);

		}
	}

}
//====================================================================================================

template<typename OP_MAT, typename VEC_MAT, typename PSI_MAT>
Matrix getKrylovGS(const OP_MAT & H, const PSI_MAT & psi, const Matrix & GS, VEC_MAT & v, VEC_MAT & w) 
{
	auto KryDim = GS.rows();
	v = psi / std::sqrt(norm(psi));
	w.setZero();
	auto alpha = 0.0;
	VEC_MAT psi_gs(psi.rows(), 1);
	psi_gs = GS(0,0) * v;
	alpha  = (v.adjoint() * H * v)(0,0).real();
	w = H * v - alpha * v;
	double n = norm(w);
	w /= std::sqrt(n);
	for(auto idx = 1; idx < KryDim; ++idx)
	{
		if(!(idx % 2))
		{
			psi_gs += GS(idx, 0) * v;
			alpha  = (v.adjoint() * H * v)(0,0).real();

			w = H * v - alpha * v - std::sqrt(n) * w;
			n = norm(w);
			w /= std::sqrt(n);
		}
		else
		{
			psi_gs += GS(idx, 0) * w;
			alpha = (w.adjoint() * H * w)(0,0).real();

			v = H * w - alpha * w - std::sqrt(n) * v;
			n = norm(v);
			v /= std::sqrt(n);

		}
	}

	return psi_gs;
}

template<typename T>
class LanczosSolver
{
	using index_type = T;
	private:
		int kry_dim_;
		int ndim_;
		Matrix helper1, helper2;
		double err;
		Eigen::SelfAdjointEigenSolver<Matrix> eigensolver;

	public:
		double error(){return err;}

		LanczosSolver(const index_type & index, int kry_dim):
			kry_dim_(kry_dim),
			ndim_(index.range()),
			helper1(ndim_, 1),
			helper2(ndim_, 1)
	{		
	}


		template<typename OP_MAT>
		void evolve(const OP_MAT & H, State<index_type> & psi, int tstp, double h)noexcept
		{
			assert(tstp < psi.nt() - 1);
			auto Mat = Hermitian::getKrylovReprAd(H, psi(tstp).mat(), helper1, helper2, kry_dim_, err);
			eigensolver.compute(Mat);
			if(eigensolver.info() != Eigen::Success) 
			{
				std::cout << "failed to solve the Krylov matrix" << std::endl;
				abort();
			}
			Matrix diagkry(Mat.rows(), Mat.cols());
			diagkry.setIdentity();

			for(auto j : range(Mat.rows())) 
				diagkry(j,j) = std::cos(eigensolver.eigenvalues()(j,0) * h) - II * std::sin(eigensolver.eigenvalues()(j,0) *h);

			Matrix kry_psi(psi.ndim(), 1);
			kry_psi = (eigensolver.eigenvectors() * diagkry * eigensolver.eigenvectors().adjoint()).block(0,0,psi.ndim(),1);

			Hermitian::getKrylovGS(H, psi(tstp).mat(), kry_psi, helper1, helper2, psi(tstp + 1).mat());
		}

		template<typename OP_MAT>
		void getGroundState(const OP_MAT & H, State<index_type> & psi, int tstp)noexcept
		{
			assert(tstp < psi.nt());
			auto Mat = Hermitian::getKrylovReprAd(H, psi(tstp), helper1, helper2, kry_dim_, err);
			eigensolver.compute(Mat);
			if(eigensolver.info() != Eigen::Success) 
			{
				std::cout << "failed to solve the Krylov matrix" << std::endl;
				abort();
			}
			psi(tstp) = Hermitian::getKrylovGS(H, psi(tstp), eigensolver.eigenvectors().col(0), helper1, helper2);
		
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
LanczosSolver<T> getSolver(const T & index, int kry_dim)
{
	return LanczosSolver<T>(index, kry_dim);
}

}//namespace hermitian

//====================================================================================================
namespace Arnoldi
{
template<typename OP_MAT, typename PSI_MAT>
std::vector<Matrix> getKrylovRepr(const OP_MAT & H, const PSI_MAT & psi, size_t KryDim) 
{
	std::vector<Matrix> v;
	v.reserve(10);
	v.push_back(psi / std::sqrt(Hermitian::norm(psi)));

	std::vector<value_type> h;
	h.reserve(10);
	auto idx = 1;
	auto n = 10.0;

	while(n > eps)
	{
		v.push_back(H * v[idx - 1]);
		for(auto i = 0; i < idx; ++i)
		{
			value_type new_h = (v[i].adjoint() * v[idx])(0,0);
			v[idx] -= new_h * v[i];
			h.push_back(new_h);
		}
		n = Hermitian::norm(v[idx]);
		v[idx] /= std::sqrt(n);
		h.push_back(std::sqrt(n));

		++idx;
		if(idx > KryDim) break;

	}
	auto kry_dim = v.size() - 1;
	Matrix KR(kry_dim, kry_dim);
	Matrix V(v[0].rows(), kry_dim);
	KR.setZero();
	idx = 0;
	for(auto i = 0; i < kry_dim - 1; ++i)
	{
		for(auto j = 0; j <= i + 1; ++j)
		{
			KR(j,i) = h[idx];
			++idx;
		}

		V.col(i) = v[i];
	}
	for(auto j = 0; j <= kry_dim - 1; ++j)
	{
		KR(j,kry_dim - 1) = h[idx];
		++idx;
	}

	V.col(kry_dim - 1) = v[kry_dim - 1];

	Matrix err(1,1);	
	err(0,0) = n;
	return std::vector<Matrix>{KR, V, err};	
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
		n = std::sqrt(Hermitian::norm(V.col(idx)));
		V.col(idx) /= n;
		h.push_back(n);
		if(n < eps) 
		{
			++idx;
			break;
		}
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

		n = std::sqrt(Hermitian::norm(V.col(idx0)));
		V.col(idx0) /= n;

		/*
		V.col(idx) = H * V.col(idx - 1);
		for(auto i = 0; i < idx - 2; ++i) h.push_back(0.);
		for(auto i = std::max(idx - 2,0); i < idx; ++i)
		{
			value_type new_h = (V.col(i).adjoint() * V.col(idx))(0,0);
			V.col(idx) -= new_h * V.col(i);
			h.push_back(new_h);
		}
		n = std::sqrt(Hermitian::norm(V.col(idx)));
		V.col(idx) /= n;

		*/
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

/*
template<bool herm = 0>
void initV(Matrix & V, int ndim, int kry_dim)
{
	V.resize(ndim, kry_dim + 1);
}

template<>
void initV<1>(Matrix & V, int ndim, int kry_dim)
{
	V.resize(ndim, 3);
}
*/

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
		//initV<herm>(psi_V, ndim_, kry_dim_);
		initV0<herm>(V0, ndim_);
	}

		template <typename PSI_MAT>
		void initKrylov(const PSI_MAT & psi)
		{
			KR.setZero();
			V.col(0) = psi / std::sqrt(Hermitian::norm(psi));
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
		void evolve(const OP_MAT & H, State<index_type> & psi, int tstp, double h, int kry_dim_rt = 0)noexcept
		{
			assert(tstp < psi.nt() - 1);
			if(kry_dim_rt == 0) kry_dim_rt = kry_dim_;
			//getKrylovRepr(H, psi(tstp).mat());

			act_dim_ = 0;
			initV0<herm>(V0, ndim_);
			initKrylov(psi(tstp).mat());
			updateKrylovRepr(H, kry_dim_rt - 1);
			
			diagkry = MatrixExp(-II * KR.block(0, 0, act_dim_, act_dim_) * h);
			psi(tstp + 1).mat() = (V.block(0,0,ndim_,act_dim_) * diagkry).block(0, 0, ndim_, 1);
		}

		template<typename OP_MAT>
		void getGroundState(const OP_MAT & H, State<index_type> & psi, int tstp, int add_dim)noexcept
		{
			if(act_dim_ >= kry_dim_) return;
			assert(tstp < psi.nt());
			//value_type E0 = psi.eval(tstp, H);
			//value_type E0 = eigensolver.eigenvalues()(0,0);
			if(act_dim_ == 0) initKrylov(psi(tstp).mat());
			updateKrylovRepr(H, add_dim);
			eigensolver.compute(KR.block(0,0,act_dim_,act_dim_));
			if(eigensolver.info() != Eigen::Success) 
			{
				std::cout << "failed to solve the Krylov matrix" << std::endl;
				abort();
			}
			getState(V, H, eigensolver.eigenvectors(), psi(tstp).mat(), V0, 0, herm_type<herm>());
			err_ = std::sqrt(Hermitian::norm((eigensolver.eigenvalues().block(0,0,std::min((int)evs_.rows(),act_dim_),1)
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
}//namespace exact

