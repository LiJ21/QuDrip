#pragma once
#include <type_traits>
#include <cmath>

namespace qudrip
{
//------------------------------------------------------------------------------

struct MUL {};
struct PLU {};
struct MIN {};
struct DIV {};

//------------------------------------------------------------------------------
template<typename VAL, typename OP>
struct dispatcher {};

template<typename VAL>
struct dispatcher<VAL, MUL>
{
VAL dispatch (const VAL & lhs, const VAL & rhs) {return lhs * rhs;}
};
template<typename VAL>
struct dispatcher<VAL, PLU>
{
VAL dispatch(const VAL & lhs, const VAL & rhs) {return lhs + rhs;}
};
template<typename VAL>
struct dispatcher<VAL, MIN>
{
VAL dispatch (const VAL & lhs, const VAL & rhs) {return lhs - rhs;}
};
template<typename VAL>
struct dispatcher<VAL, DIV>
{
VAL dispatch (const VAL & lhs, const VAL & rhs) {return lhs / rhs;}
};

//------------------------------------------------------------------------------
template <typename LHS, typename RHS, typename OP, typename VAL = idx_size_t>
class CQContainer
{
	public:
	using quant_type = VAL;
	private:
	using lhs_type = LHS;
	using rhs_type = RHS;

	lrHandle<LHS> LH;
	lrHandle<RHS> RH;
	struct dispatcher<quant_type, OP> dispatch;
	typename std::remove_reference<lhs_type>::type lhs_;
	typename std::remove_reference<rhs_type>::type rhs_;
	quant_type val_;

	public:
	struct is_cq{};
	CQContainer(lhs_type lhs, rhs_type rhs) : 
		LH(lhs), RH(rhs), lhs_(LH.val_), rhs_(RH.val_)
	{}

	operator quant_type()
	{
		return dispatch.dispatch (lhs_, rhs_);
	}

	CQContainer & operator=(quant_type val)
	{
		val_ = val;
		return *this;
	}

	quant_type val() {return val_;}
};
//CRTP
template <typename sub, typename T>
class ConservedQuantityBase
{
	public:
	using quant_type = T;
	private:
	quant_type val_;
	public:
	struct is_cq{};
	ConservedQuantityBase() {}

	sub & operator=(quant_type val)
	{
		val_ = val;
		return *static_cast<sub*>(this);
	}

	operator quant_type()
	{
		return static_cast<sub*>(this) -> eval();
	}

	quant_type val() {return val_;}
	
};
//------------------------------------------------------------------------------

template <typename T1 = idx_size_t>
class Nset : public ConservedQuantityBase<Nset<T1>, T1>
{
	public:
	using quant_type = T1;
	private:
	friend class ConservedQuantityBase<Nset<T1>, T1>;
	using base = ConservedQuantityBase<Nset<T1>, T1>;
	using index_type = QbitsIndex<idx_size_t>;

	const index_type & index;
	std::bitset<BIT_LIMIT> mask;

	quant_type eval()
	{
		return (std::bitset<BIT_LIMIT>(quant_type(index)) & mask).count();	
	}

	public:
	using base::operator=;

	Nset(const index_type & index, int start = -1, int end = -1) : 
		index(index), mask(0)
	{
		if(start == -1) start = 0;
		if(end == -1) end = bits_type(index.range() - 1).count();
		for(auto i : range(start, end + 1))
			mask[i] = 1;
	
	}
};

template <typename T1 = idx_size_t>
class Nd : public ConservedQuantityBase<Nd<T1>, T1>
{
	public:
	using quant_type = T1;
	private:
	int Nsite_;
	friend class ConservedQuantityBase<Nd<T1>, T1>;
	using base = ConservedQuantityBase<Nd<T1>, T1>;
	using index_type = QbitsIndex<idx_size_t>;

	const index_type & index;

	std::bitset<BIT_LIMIT> mask;
	bool hole_;

	quant_type eval()
	{
		auto biti = std::bitset<BIT_LIMIT>(quant_type(index));
		return (((biti >> Nsite_)^mask) & (biti^mask)).count();	
	}

	public:
	using base::operator=;

	Nd(const index_type & index, int Nsite, bool hole = 0) : 
		index(index), Nsite_(Nsite), hole_(hole), mask(0)
	{
		for(auto i : range(0, Nsite))
		{
			mask[i] = hole_;
		}
		/*
		for(auto i : range(0, Nsite))
		{
			mask1[i] = 1;
		}
		for(auto i : range(Nsite, 2 * Nsite))
		{
			mask2[i] = 1;
		}
		*/

	}
};
//------------------------------------------------------------------------------
template<typename LHS, typename RHS, 
	typename = typename std::remove_reference<LHS>::type::is_cq, 
	typename = typename std::remove_reference<RHS>::type::is_cq>
auto operator+(LHS && lhs, RHS && rhs)
{
	using NL = typename std::remove_reference<LHS>::type::quant_type;
	using NR = typename std::remove_reference<RHS>::type::quant_type;
	using QT = typename std::conditional<(sizeof(NL) > sizeof(NR)), NL, NR>::type;

	return CQContainer<LHS, RHS, PLU, QT> (lhs, rhs);
}
template<typename LHS, typename RHS, 
	typename = typename std::remove_reference<LHS>::type::is_cq, 
	typename = typename std::remove_reference<RHS>::type::is_cq>
auto operator-(LHS && lhs, RHS && rhs)
{
	using NL = typename std::remove_reference<LHS>::type::quant_type;
	using NR = typename std::remove_reference<RHS>::type::quant_type;
	using QT = typename std::conditional<(sizeof(NL) > sizeof(NR)), NL, NR>::type;

	return CQContainer<LHS, RHS, MIN, QT> (lhs, rhs);
}
template<typename LHS, typename RHS, 
	typename = typename std::remove_reference<LHS>::type::is_cq, 
	typename = typename std::remove_reference<RHS>::type::is_cq>
auto operator*(LHS && lhs, RHS && rhs)
{
	using NL = typename std::remove_reference<LHS>::type::quant_type;
	using NR = typename std::remove_reference<RHS>::type::quant_type;
	using QT = typename std::conditional<(sizeof(NL) > sizeof(NR)), NL, NR>::type;

	return CQContainer<LHS, RHS, MUL, QT> (lhs, rhs);
}
template<typename LHS, typename RHS, 
	typename = typename std::remove_reference<LHS>::type::is_cq, 
	typename = typename std::remove_reference<RHS>::type::is_cq>
auto operator/(LHS && lhs, RHS && rhs)
{
	using NL = typename std::remove_reference<LHS>::type::quant_type;
	using NR = typename std::remove_reference<RHS>::type::quant_type;
	using QT = typename std::conditional<(sizeof(NL) > sizeof(NR)), NL, NR>::type;

	return CQContainer<LHS, RHS, DIV, QT> (lhs, rhs);
}


//============================================================================
template<size_t tight, typename CQ, typename... CQs>
bool Evaluate(CQ && cq, CQs &&... cqs)
{
	return (std::abs(int(cq - cq.val())) <= tight) && Evaluate<tight>(std::forward<CQs>(cqs)...);
}

template<size_t tight, typename CQ, typename... CQs>
bool Evaluate(CQ && cq)
{
	return (std::abs(int(cq - cq.val())) <= tight);
}

template<size_t tight, bool P, typename IDX, typename CQ, typename... CQs>
auto LooseConstrain(IDX && index, CQ && cq, CQs&&... cqs)
{
	std::cout << "Constrain: original range " << index.range() << std::endl;
	auto & origin_index = strip(index);
	auto sub = SubIndex<decltype(strip(index)), P>(strip(index));
	for(index[0]; index < index.range(); ++index)
	{
		if(Evaluate<tight>(std::forward<CQ>(cq), std::forward<CQs>(cqs)...))
		{
			sub.mapping.push_back(origin_index);
			sub.inv_mapping.insert({origin_index, sub.mapping.size() - 1});
		}
	}
	sub[0];
	return sub;
}

template<typename IDX, typename... CQs>
auto Constrain(IDX && index, CQs&&... cqs)
{
	return LooseConstrain<0, true>(std::forward<IDX>(index), std::forward<CQs>(cqs)...);
}
//============================================================================
template<typename index_type>
auto getNset(index_type && index)
{
	return Nset<>(std::forward<index_type>(index), -1, -1);
}

template<typename index_type>
auto getNup(index_type && index, int Nsite)
{
	return Nset<>(std::forward<index_type>(index), 0, Nsite - 1);
}

template<typename index_type>
auto getNdo(index_type && index, int Nsite)
{
	return Nset<>(std::forward<index_type>(index), Nsite, 2 * Nsite - 1);
}

template<typename index_type>
auto getNd(index_type && index, int Nsite)
{
	return Nd<>(std::forward<index_type>(index), Nsite, 0);
}

template<typename index_type>
auto getNh(index_type && index, int Nsite)
{
	return Nd<>(std::forward<index_type>(index), Nsite, 1);
}

}//namespace

