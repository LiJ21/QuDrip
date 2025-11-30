#pragma once
#include<string>
#include<unordered_map>
#include<Eigen/Dense>

namespace qudrip
{

//=====================================================================================================================
using c1 = StaticIVector<1>;
using c2 = StaticIVector<2>;
using c3 = StaticIVector<3>;

template<typename T>
class params_type
{
	using param_val_type = T;
	std::unordered_map<std::string, param_val_type> params_;

	public:

	param_val_type & operator()(const std::string & name)
	{
		auto it = params_.find(name);
		if(it == params_.end())
		{
			auto nit = params_.insert(std::make_pair(name, param_val_type{}));
			return nit.first -> second;
		}
		return it -> second;
	}

	const param_val_type & operator()(const std::string & name)const
	{
		auto it = params_.find(name);
		if(it == params_.end()) throw(std::string("undefined params"));
		return it -> second;
	}

	std::string stringize()const
	{
		std::string ss = "";
		for(auto & el : params_)
		{
			ss += el.first + "=" + std::to_string(el.second) + ",";
		}

		ss[ss.size() - 1] = ';';
		return ss;
	}
};

template<typename T>
std::ostream& operator<< (std::ostream & stream, const params_type<T> & params)
{
	stream << params.stringize();
	return stream;
}

//=====================================================================================================================

template<int N>
int dispatch_getr(std::vector<StaticIVector<N>> & rpvec, const std::vector<StaticIVector<N>> & pvec)
{
	return 0;
}

template<int N>
int get_reci(std::vector<StaticIVector<N>> & rpvec, const std::vector<StaticIVector<N>> & pvec)
{
	assert(rpvec.size() == N);
	assert(pvec.size() == N);

	return dispatch_getr(rpvec, pvec);
}

template<>
int dispatch_getr<1>(std::vector<StaticIVector<1>> & rpvec, const std::vector<StaticIVector<1>> & pvec)
{
	rpvec[0](0) = pvec[0](0);
	return StaticInnerProduct(pvec[0], pvec[0]);
	//rpvec[0] *= 2. * M_PI / StaticInnerProduct(pvec[0], pvec[0]).real();
}

template<>
int dispatch_getr<2>(std::vector<StaticIVector<2>> & rpvec, const std::vector<StaticIVector<2>> & pvec)
{
	rpvec[0](0) = -pvec[1](1);
	rpvec[0](1) = pvec[1](0);
	//rpvec[0] *= 2. * M_PI / StaticInnerProduct(rpvec[0], pvec[0]).real();

	rpvec[1](0) = pvec[0](1);
	rpvec[1](1) = -pvec[0](0);
	//rpvec[1] *= 2. * M_PI / StaticInnerProduct(rpvec[1], pvec[1]).real();
	return StaticInnerProduct(rpvec[0], pvec[0]);
}

template<>
int dispatch_getr<3>(std::vector<StaticIVector<3>> & rpvec, const std::vector<StaticIVector<3>> & pvec)
{
	StaticCrossProduct(rpvec[0], pvec[1], pvec[2]);
	StaticCrossProduct(rpvec[1], pvec[2], pvec[0]);
	StaticCrossProduct(rpvec[2], pvec[0], pvec[1]);

	return StaticInnerProduct(pvec[0], StaticCrossProduct(pvec[1], pvec[2]));

	//rpvec[0] *= 2. * M_PI / norm;
	//rpvec[1] *= 2. * M_PI / norm;
	//rpvec[2] *= 2. * M_PI / norm;
}

//=====================================================================================================================

template<int NDIM>
class Lattice
{
	static const int ndim_ = NDIM;
	//using coord_type = std::array<int, ndim_>;
public:
	using coord_type = StaticIVector<ndim_>;
	using coords_type = std::vector<coord_type>;
	//using rcoord_type = StaticDVector<ndim_>;
	//using rcoords_type = std::vector<rcoord_type>;
private:

	size_t N_;
	int vol_;
	coords_type coord_;
	//std::unordered_map<coord_type, size_t> index_;
	//std::vector<double> norm_;
	coords_type pvec_, rpvec_;
	//rcoords_type rpvec_; //periodicity vectors and reciprocals

	//void create_index()
	//{
		//for(auto i = 0; i < coord_.size(); ++i)
		//{
		//	index_.insert(std::make_pair(coord_[i], i));
		//}
	//}
	size_t index(const coord_type & c)const
	{
		for(int i = 0; i < coord_.size(); ++i)
		{
			if(coord_[i] == c) return i;
		}
		return N_;
	}

public:

	Lattice(const coords_type & coord, const coords_type & pvec):
		coord_(coord), N_(coord.size()), 
		pvec_(pvec), rpvec_(ndim_)
	{
		assert(pvec_.size() == ndim_);
		vol_ = get_reci<ndim_>(rpvec_, pvec_);

		for(auto & p : coord_) for(int i = 0; i < p.rows(); ++i) 
			if(p(i) < 0) 
			{
				std::cout << "Warning: non-positive coordinates are invalid for periodicity reduction." << std::endl;
				return;
			}
		//create_index();
	}

	Lattice(const coords_type & coord): Lattice(coord, coords_type{})
	{}
	
	size_t N() const {return N_;}
	int vol() const {return vol_;}

	coord_type reduce(const coord_type & c)const
	{
		coord_type r = c;	
		if(pvec_.size() == 0) return r;
		//if(index(r) != N_) return r;

		for(int i = 0; i < ndim_; ++i)
		{
			auto prod = StaticInnerProduct(r, rpvec_[i]);
			auto ir = StaticInnerProduct(r, rpvec_[i]) / vol_;
			if(prod != 0 && prod < 0 != vol_ < 0) --ir;
			r -= ir * pvec_[i];
			//std::cout << i << "," << ir << std::endl;
			//std::cout << "c = (" << c.transpose() << "), prod = " + std::to_string(prod) + ", ir = " + std::to_string(ir) + ", r = (" << r.transpose() << ")" << std::endl;
		}
		return r;
	}
	
	const coord_type & operator()(size_t i) const {assert(i < N_); return coord_[i];}

	size_t operator()(const coord_type & c) const {return index(reduce(c));}
	const coord_type & pvec(size_t i) const {assert(i < ndim_); return pvec_[i];}
	//const rcoord_type & rpvec(size_t i) const {assert(i < ndim_); return pvec_[i];}
	const coord_type & rpvec(size_t i) const {assert(i < ndim_); return rpvec_[i];}
};

//=====================================================================================================================

template<typename T, int NDIM>
class LatticeFunction
{
	using val_type = T;
	static const int ndim_ = NDIM;
	using coord_type = StaticIVector<ndim_>;
	using arg_type = std::vector<int>;

	Lattice<ndim_> latt_;
	std::vector<val_type> data_;
	int argn_;

	int convert_idx(const arg_type & arg)
	{
		int idx = 0;
		for(int i = 0; i < arg.size(); ++i) 
		{
			idx *= latt_.N();
			idx += arg[i];
		}
		return idx;
	}

	void convert_arg(arg_type & arg, int idx)
	{
		for(int i = arg.size() - 1; i >=0; --i)
		{
			arg[i] = idx % latt_.N();
			idx /= latt_.N();
		}
	}

public:

	LatticeFunction(const Lattice<ndim_> & latt, int argn = 1):
	latt_(latt), argn_(argn), data_(std::pow(latt.N(), argn))
	{}

	const val_type & operator()(size_t i) const 
	{
		assert(argn_ == 1);
		return data_[i];
	}

	const val_type & operator()(const arg_type & arg) const 
	{
		return data_[convert_idx(arg)];
	}

	val_type & operator()(size_t i) 
	{
		return const_cast<val_type &>(
			const_cast<const LatticeFunction<T,NDIM> *>(this)->operator()(i)
			);
	}

	val_type & operator()(const arg_type & arg) 
	{
		return const_cast<val_type &>(
			const_cast<const LatticeFunction<T,NDIM> *>(this)->operator()(arg)
			);
	}

	val_type DFT(const coord_type & k) const //k is momentum (reciprocal lattice)
	{
		assert(argn_ == 1);
		val_type result;
		for(size_t i = 0; i < latt_.N(); ++i)
		{		
			//auto & x = this -> operator()(i);
			auto & x = latt_(i);
			double arg = 0.;
			for(int j = 0; j < ndim_; ++j) 
				arg += StaticInnerProduct(x, k(j) * latt_.rpvec(j)) / (double) latt_.vol();
			result += data_[i] * value_type(std::cos(arg), std::sin(arg));
		}
		result /= latt_.N();
		return result;
	}

	val_type DFT(const std::vector<coord_type> & ks) const //k is momentum (reciprocal lattice)
	{
		assert(argn_ == ks.size());
		val_type result;
		for(size_t i = 0; i < data_.size(); ++i)
		{		
			arg_type argi;
			convert_arg(argi, i);
			double arg = 0.;
			for(int ix = 0; ix < argn_; ++ix)
			{
				auto & x = latt_(argi[ix]);
				for(int j = 0; j < ndim_; ++j) 
					arg += StaticInnerProduct(x, ks[i](j) * latt_.rpvec(j)) 
						/ (double) latt_.vol();
			}
			result += data_[i] * value_type(std::cos(arg), std::sin(arg));
		}
		result /= data_.size();
		return result;
	}

};

//=====================================================================================================================

template<typename T, int NDIM>
class Model
{ 
	static const int ndim_ = NDIM;
	using param_val_type = T;

	std::vector<std::pair<int, int>> bonds_;
	std::vector<params_type<param_val_type>> params_;
	Lattice<ndim_> latt_;

	public:
	
	Model(const Lattice<ndim_> & latt):
		latt_(latt)
	{}

	void addBond(int s1, int s2, const params_type<param_val_type> & param)
	{
		assert(s1 < latt_.N() && s2 < latt_.N()); 
		bonds_.push_back(std::make_pair(s1, s2)); 
		params_.push_back(param);
	}

	//void addBond(const typename Lattice<NDIM>::coord_type & c1, const typename Lattice<NDIM>::coord_type & c2, const params_type & param)
	//{
	//	addBond(latt_(c1), latt_(c2), param);
	//}

	void print()const
	{
		std::cout << "Bonds: " << std::endl;
		for(int i = 0; i < bonds_.size(); ++i)
		{
			auto & site1 = latt_(bonds_[i].first);
			auto & site2 = latt_(bonds_[i].second);

			std::cout << "(" + std::to_string(site1(0));
			for(int j = 1; j < site1.rows(); ++j) std::cout << "," + std::to_string(site1(1));
			std::cout << ") and (" + std::to_string(site2(0));
			for(int j = 1; j < site1.rows(); ++j) std::cout << "," + std::to_string(site2(1));
			std::cout << "): " << params_[i] << std::endl;
		}

	}

	template<typename FUNC>
	SpMatrix getHamiltonianMatrix(FUNC && f, int ndim) 
	{
		SpMatrix ham_op(ndim, ndim);
		ham_op.setZero(); 
		for(int i = 0; i < bonds_.size(); ++i)
		{
			auto & ss = bonds_[i];
			auto & param = params_[i];
			int s1 = ss.first;
			int s2 = ss.second;
			ham_op = ham_op + f(s1, s2, param);
		}

		return ham_op;
	}

	template<typename FUNC>
	auto getHamiltonian(FUNC && f) 
	{
		auto ham_op = getEmptyOperator();
		for(int i = 0; i < bonds_.size(); ++i)
		{
			auto & ss = bonds_[i];
			auto & param = params_[i];
			int s1 = ss.first;
			int s2 = ss.second;
			ham_op = ham_op + f(s1, s2, param);
		}

		return ham_op;
	}
};

}//namespace
