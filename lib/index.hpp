#include<vector>
#include<iostream>
#include<array>
#include<algorithm>
#include<complex>
#include<bitset>
#include<unordered_map>
#include<tuple>


#pragma once
namespace qudrip
{
template<typename T, bool P = false>
class SubIndex;

template<size_t = 0, bool = false, typename IDX, typename CQ, typename... CQs>
	auto LooseConstrain(IDX && index, CQ && cq, CQs&&... cqs);

template<typename IDX, typename... CQs>
	auto Constrain(IDX && index, CQs&&... cqs);

template<typename IDX,
	typename std::remove_reference_t<IDX>::ear::type = 0>
decltype(auto) strip(IDX && index);

template<typename IDX,
	typename std::remove_reference_t<IDX>::is_earless::type = 0>
decltype(auto) strip(IDX && index);

template<typename IDX,
	typename std::remove_reference_t<IDX>::ear::type,
	typename = typename std::enable_if<IDX::pretend_ == true>::type>
decltype(auto) strip(IDX && index);


using bits_type = std::bitset<BIT_LIMIT>;

//CRTP
template<class Derived, class T = idx_size_t>
class IndexBase
{
	T range_;

	protected:
	T idx_value_;

	public:
	struct is_index {using type = int;};
	struct is_earless {using type = int;};

	IndexBase(T range) : range_(range)
	{}

	//bool exists() {return idx_value_ < range_;}

	operator T() const
	{
		return idx_value_;
	}

	void operator[](T idx_val)
	{
		idx_value_ = idx_val;

	}

	template <class IDX>
	Derived & operator=(IDX && idx) noexcept
	{
		static_cast<Derived *> (this) -> t_idx(std::forward <IDX> (idx));
		return (* static_cast<Derived *> (this));
	}

	Derived & operator++() noexcept
	{
		++idx_value_;
		return (* static_cast<Derived *> (this));
	}

	T range() const{return range_;}
	
};
//------------------------------------------------------------------

template<typename int_type>
int_type bit_convert(const bits_type & idx){return idx.to_ulong();}
template<>
unsigned long long bit_convert(const bits_type & idx) {return idx.to_ullong();}


template <class T = idx_size_t>
class QbitsIndex : public IndexBase<QbitsIndex<T>, T>
{
	using base = IndexBase<QbitsIndex<T>, T>;
	friend class IndexBase<QbitsIndex<T>, T>;

	void t_idx(const bits_type & idx) 
	{
		base::idx_value_ = bit_convert<T> (idx);

	}
	void t_idx(std::string idx) 
	{
		t_idx(bits_type(idx));
	}

	bits_type content()const
	{
		return bits_type(base::idx_value_);
	}

	public:
	using IndexBase<QbitsIndex<T>, T>::operator=;

	QbitsIndex(T N) : IndexBase<QbitsIndex<T>, T>(T(1) << N)
	{}

};
//------------------------------------------------------------------

template<class T = idx_size_t>
class SingleModeIndex : public IndexBase<SingleModeIndex<T>, T>
{
	using base = IndexBase<SingleModeIndex<T>, T>;
	friend class IndexBase<SingleModeIndex<T>, T>;
	int n_mode_;
	void t_idx(T N){base::idx_value_ = N;}
	T content()const{return base::idx_value_;}
	public:
	using IndexBase<SingleModeIndex<T>, T>::operator=;
	SingleModeIndex(int n_mode) : IndexBase<SingleModeIndex<T>, T>(n_mode),
		n_mode_(n_mode)
	{}

};
//------------------------------------------------------------------

template<class T1, class T2>
class Indices
{
	using lhs_type = T1;
	using rhs_type = T2;

	lrHandle<lhs_type> LH;
	lrHandle<rhs_type> RH;
	typename std::remove_reference<lhs_type>::type & lhs_;
	typename std::remove_reference<rhs_type>::type & rhs_;
	idx_size_t range_;

	public:
	struct is_indices {using type = int;};
	struct is_index {using type = int;};
	struct is_earless{using type = int;};

	Indices(lhs_type lhs, rhs_type rhs) : 
		LH(lhs), RH(rhs),
		lhs_(LH.val_), rhs_(RH.val_), range_(lhs_.range() * rhs_.range())
	{}

	Indices(const Indices & idx) :
		LH(idx.LH), RH(idx.RH),
		lhs_(LH.val_), rhs_(RH.val_), range_(lhs_.range() * rhs_.range())
	{}

	Indices(Indices && idx) :
		LH(std::move(idx.LH)), RH(std::move(idx.RH)),
		lhs_(LH.val_), rhs_(RH.val_), range_(lhs_.range() * rhs_.range())
	{}

	operator idx_size_t() const
	{
		return idx_size_t(lhs_) * (rhs_.range()) + idx_size_t(rhs_);
	}

	idx_size_t range() const {
		return range_;
	}
//=========================================================================================
	void operator[](idx_size_t i)
	{
		lhs_[i / rhs_.range()]; 
		rhs_[i % rhs_.range()];
	}

//=========================================================================================
	Indices & operator++()
	{
		if(rhs_ < rhs_.range() - 1) ++rhs_;
		else {rhs_[0]; ++lhs_;}
		return *this;
	}

};

//------------------------------------------------------------------
template<class T1, class T2, typename = typename std::remove_reference<T1>::type::is_index, typename = typename std::remove_reference<T2>::type::is_index>
Indices<T1, T2> operator*(T1 && lhs, T2 && rhs)
{
	return Indices<T1, T2>(lhs, rhs);
}

//------------------------------------------------------------------
template<typename T, bool P>
class SubIndex
{
	template<size_t, bool, typename IDX, typename CQ, typename... CQs>
		friend auto LooseConstrain(IDX && index, CQ && cq, CQs&&... cqs);
	template<typename IDX, typename... CQs>
		friend auto Constrain(IDX && index, CQs&&... cqs);

	template<typename IDX,
		typename std::remove_reference_t<IDX>::ear::type>
		friend decltype(auto) strip(IDX && index);

	template<typename IDX,
		typename std::remove_reference_t<IDX>::is_earless::type>
		friend decltype(auto) strip(IDX && index);

	template<typename IDX,
		typename std::remove_reference_t<IDX>::ear::type,
		typename>
		friend decltype(auto) strip(IDX && index);


	using index_type = T;
	lrHandle<index_type> IH;
	typename std::remove_reference<index_type>::type & earless_idx_;
	std::vector<idxv_type> mapping;
	std::unordered_map<idxv_type, idxv_type> inv_mapping;
	idxv_type idx_value_;

	public:
	const static bool pretend_ = P;
	struct is_index {using type = int;};
	struct ear{using type = int;};
	SubIndex(index_type idx) : IH(idx), earless_idx_(IH.val_)
	{}

	SubIndex(const SubIndex & sub) : 
		IH(sub.IH), earless_idx_(IH.val_),
		mapping(sub.mapping), inv_mapping(sub.inv_mapping)
	{}

	SubIndex(SubIndex && sub) : 
		IH(std::move(sub.IH)), earless_idx_(IH.val_),
		mapping(std::move(sub.mapping)), inv_mapping(std::move(sub.inv_mapping))
	{}

	operator idxv_type() const
	{
		auto it = inv_mapping.find(idxv_type(earless_idx_));
		if(it == inv_mapping.end())
			return range();
		else
			return it -> second;
		//return idx_value_;
	}

	//bool exists() { return (idx_value_ < range() && inv_mapping.find(idxv_type(earless_idx_)) != inv_mapping.end()); }

	void operator[](idxv_type i)
	{
		idx_value_ = i;
		if(idx_value_ < mapping.size())
			earless_idx_[mapping[idx_value_]];
		else
			//earless_idx_[mapping[0]];
			earless_idx_[earless_idx_.range()];
	}

	SubIndex & operator++()
	{
		++idx_value_;
		if(idx_value_ < mapping.size())
			earless_idx_[mapping[idx_value_]];
		else
			//earless_idx_[mapping[0]];
			earless_idx_[earless_idx_.range()];

		return *this;
	}

	idxv_type range()const {return mapping.size();}
};

//------------------------------------------------------------------
template <size_t start, typename T, size_t N, typename = typename T::is_index, 
	typename std::enable_if<(start < N - 1), int>::type = 0 >
constexpr auto product(std::array<T, N> & idxs)
{
	return idxs[start] * product<start + 1>(idxs);
}

template <size_t start, typename T, size_t N, typename = typename T::is_index, 
	typename std::enable_if<(start == N - 1), int>::type = 0 >
constexpr T product(std::array<T, N> & idxs)
{
	return idxs[start];
}

template <typename T, size_t N, typename = typename T::is_index>
constexpr auto product(std::array<T, N> & idxs)
{
	return product<0>(idxs);
}

//------------------------------------------------------------------

template<typename IDX,
	typename std::remove_reference_t<IDX>::ear::type,
	typename>
decltype(auto) strip(IDX && index)
{
	return std::forward<IDX>(index);
}

template<typename IDX,
	typename std::remove_reference_t<IDX>::ear::type>
decltype(auto) strip(IDX && index)
{
	return strip(index.earless_idx_);
}

template<typename IDX,
	typename std::remove_reference_t<IDX>::is_earless::type>
decltype(auto) strip(IDX && index)
{
	return index;
}

//------------------------------------------------------------------
auto getQbits(int N)
{
	return QbitsIndex<>(N);
}

auto getSingleMode(int Nph)
{
	return SingleModeIndex<>(Nph);
}

//------------------------------------------------------------------
bits_type getQbitsLabel(int N)
{
	return bits_type(N);
}

//
}//namespace
