//------------------------------------------------------------------

#pragma once
namespace qudrip
{
using idx_tree_type = std::vector<idxv_type>;
using val_tree_type = std::vector<value_type>;
using idx_it_type = idx_tree_type::iterator;
using val_it_type = val_tree_type::iterator;

//=========================================================================

class Operator;

using U_out_type = std::vector<value_type>;
class elOpBase
{
	protected:
		//std::vector<U_out_type> U_;
		U_out_type U_;
		size_t n_input, n_output;
	public:
		elOpBase(size_t n_input, size_t n_output):
			//U_(n_input, U_out_type(n_output)), 
			U_(n_input * n_output),
			n_input(n_input),
			n_output(n_output)
		{}

		virtual void operator<<(Matrix U) 
		{
			assert(U.rows() == n_output);
			assert(U.cols() == n_input);
			for(auto i = 0 ; i < n_input; ++i)
				//U_[i] = U.block(0, i, n_output, 1);
				for(auto j = 0; j < n_output; ++j)
					U_[i * n_output + j] = U(i, j);
		}

		virtual void invert()
		{
			assert(n_input == n_output);
			Matrix M(n_input,n_output);
			for(auto i = 0 ; i < n_input; ++i)
				for(auto j = 0; j < n_output; ++j)
					M(i,j) = U_[i * n_output + j];
			operator<< (M.inverse().eval());
		}

		virtual void feed_idx(idx_size_t idx) = 0;
		virtual idx_size_t get_idx() = 0;
		virtual void branch(idx_it_type & , val_it_type & , idx_it_type & , val_it_type & ) = 0;
		size_t range() {return n_output;}
};

using ops_type = std::vector<std::shared_ptr<elOpBase>>;
//=========================================================================
//auxiliary function
size_t tree_size(const ops_type & ops)
{
	size_t res = 1;
	for(auto op = ops.begin(); op != ops.end(); ++op)
	{
		res = 1 + res * ((*op) -> range());
	}
	return res;
}
//=========================================================================
template<typename OP, typename STATE>
class OPsiType
{
	friend STATE;
	OP & op;
	STATE & psi;
	public:
	class isOPsi
	{};
	OPsiType(OP & op, STATE & psi):
		op(op), psi(psi)
	{}

};
//=========================================================================
class Operator
{
	ops_type ops_; //raw pointer
	idx_tree_type idx_tree;
	val_tree_type val_tree;
	size_t n_leaves;

	public:

		Operator(const std::shared_ptr<elOpBase> & op)
		{
			ops_.push_back(op);
			idx_tree.resize(1 + op -> range());
			val_tree.resize(1 + op -> range());
			n_leaves = ops_[0] -> range();
		}

		Operator operator*(const Operator & Op) const
		{
			Operator Opp(*this);
			Opp.ops_.insert(Opp.ops_.end(), Op.ops_.begin(), 
					Op.ops_.end());
			auto size = tree_size(Opp.ops_);
			Opp.idx_tree.resize(size);	
			Opp.val_tree.resize(size);	
			Opp.n_leaves *= Op.n_leaves;

			return Opp;
		}

		template<typename trans_index_type, typename index_type>
		inline void full_branch(trans_index_type & trans_index, index_type & index, value_type init, idx_it_type &idx_b, val_it_type &val_b);

		template<typename index_type>
		void map_acc(State<index_type> & to_psi, int t1, State<index_type> & from_psi, int t2, value_type=1.0);

		template<typename index_type>
		void map(State<index_type> & to_psi, int t1, State<index_type> & from_psi, int t2, value_type=1.0);

		template<typename index_type>
		auto operator*(const State<index_type> & psi) const
		{
			return OPsiType<Operator, State<index_type>>(*this, psi);
		}

		template<typename index_type, typename = typename std::remove_reference_t<index_type>::is_index>
		SpMatrix operator>>(index_type && index);

		template<typename index_type, typename trans_index_type, 
			typename = typename std::remove_reference_t<index_type>::is_index, 
			typename = typename std::remove_reference_t<trans_index_type>::is_index>
		SpMatrix operator()(index_type && index, trans_index_type && trans_index)
		{
			return this -> operator>>(index);
		}

		bool operator==(const Operator & op)
		{
			if(ops_.size() != op.ops_.size()) return false;

			for(auto i : range(ops_.size()))
			{
				if (ops_[i] != op.ops_[i]) return false;
			}

			return true;
		}

};

//==========================================================================
template<typename trans_index_type, typename index_type>
inline void Operator::full_branch(trans_index_type & trans_index, index_type & index, value_type init, idx_it_type & idx_b, val_it_type & val_b)
{
#ifdef PRINT_TREE
	std::cout << std::bitset<6>(idxv_type(trans_index)) << std::endl;
	std::flush(std::cout);
#endif
	idx_tree[0] = idx_size_t(trans_index) + 1;
	val_tree[0] = init;

	//idx_it_type idx_b = idx_tree.begin();
	idx_b = idx_tree.begin();
	idx_it_type idx_e = idx_tree.begin() + 1;
	//val_it_type val_b = val_tree.begin();
	val_b = val_tree.begin();
	val_it_type val_e = val_tree.begin() + 1;

	for(int i_op = ops_.size() - 1; i_op >= 0; --i_op)
	{
		auto & op = ops_[i_op];
		auto idx_w = idx_e;
		auto val_w = val_e;

		for(; idx_b != idx_e; ++idx_b) //convert total index to "label index"
		{
#ifdef PRINT_TREE
			std::cout << std::bitset<6>(*idx_b - 1) << ": "; std::flush(std::cout);
#endif
			if(*idx_b != 0)
			{
				trans_index[*idx_b - 1];
				*idx_b = op -> get_idx() + 1;
			}

			auto idx_ww = idx_w;
			op -> branch(idx_b, val_b, idx_w, val_w);

			for(auto it = idx_ww; it != idx_w; ++it) //convert "label index" to the total index
			{
				if(*it != 0)
				{
					op -> feed_idx(*it - 1);
					(*it) = trans_index + 1;
				}
#ifdef PRINT_TREE
				std::cout << " " << std::bitset<6>((*it) - 1);std::flush(std::cout);
#endif

			}
#ifdef PRINT_TREE
			std::cout << std::endl << "->" << std::endl;
			std::flush(std::cout);
#endif
			++val_b;
		}

		idx_e = idx_w;
		val_e = val_w;
	}

	for(size_t j = 0; j < n_leaves; ++j)
	{
		//translate back to index
		if(*(idx_b + j) != 0)
		{
			trans_index[*(idx_b + j) - 1];
			//if(index.exists())
			if(index < index.range())
				*(idx_b + j) = index + 1;
			else
				*(idx_b + j) = 0;
		}
#ifdef PRINT_TREE
		std::cout << " " << std::bitset<6>((*(idx_b + j)) - 1); std::flush(std::cout);
#endif
	}

}
//==========================================================================
template<typename index_type>
void Operator::map(State<index_type> & to_psi, int t1, State<index_type> & from_psi, int t2, value_type coeff)
{
	to_psi(t2).mat().setZero();
	map(to_psi, t1, from_psi, t2, coeff);
}
//==========================================================================
template<typename index_type>
void Operator::map_acc(State<index_type> & to_psi, int t1, State<index_type> & from_psi, int t2, value_type coeff)
{
	auto & index = from_psi.get_index();
	auto & trans_index = strip(index);
	auto ndim = index.range();
	idx_it_type idx_b;
	val_it_type val_b;
	for(idxv_type i = 0; i < ndim; ++i)
	{
		index[i];
		full_branch(trans_index, index, coeff * from_psi(t2).mat()(i,0), idx_b, val_b);
		for(size_t j = 0; j < n_leaves; ++j)
		{
				if(*(idx_b + j) != 0 && std::abs(*(val_b + j)) > eps)
					to_psi(t1).mat()(*(idx_b + j) - 1, 0) += *(val_b +j);
		}
	}
}

//==========================================================================
template<typename index_type, 
	typename>
SpMatrix Operator::operator>>(index_type && index)
{
	auto & trans_index = strip(index);
	auto ndim = index.range();
	SpMatrix mat(ndim, ndim);
	triplets_type triplets;
	idx_it_type idx_b;
	val_it_type val_b;

	triplets.reserve(ndim);
	mat.reserve(ndim);

	//for(auto i : range(ndim))
#ifdef PRINT_TREE
	std::cout << "tree expansion..." << std::endl;
#endif
	for(idxv_type i = 0; i < ndim; ++i)
	{
		index[i];
		full_branch(trans_index, index, 1.0, idx_b, val_b);
		/*
#ifdef PRINT_TREE
		std::cout << std::bitset<6>(idxv_type(trans_index)) << std::endl;
		std::flush(std::cout);
#endif
		idx_tree[0] = idx_size_t(trans_index) + 1;
		val_tree[0] = 1.0;

		idx_it_type idx_b = idx_tree.begin();
		idx_it_type idx_e = idx_tree.begin() + 1;
		val_it_type val_b = val_tree.begin();
		val_it_type val_e = val_tree.begin() + 1;

		for(int i_op = ops_.size() - 1; i_op >= 0; --i_op)
		{
			auto & op = ops_[i_op];
			auto idx_w = idx_e;
			auto val_w = val_e;

			for(; idx_b != idx_e; ++idx_b) //convert total index to "label index"
			{
#ifdef PRINT_TREE
				std::cout << std::bitset<6>(*idx_b - 1) << ": "; std::flush(std::cout);
#endif
				if(*idx_b != 0)
				{
					trans_index[*idx_b - 1];
					*idx_b = op -> get_idx() + 1;
				}

				auto idx_ww = idx_w;
				op -> branch(idx_b, val_b, idx_w, val_w);

				for(auto it = idx_ww; it != idx_w; ++it) //convert "label index" to the total index
				{
					if(*it != 0)
					{
						op -> feed_idx(*it - 1);
						(*it) = trans_index + 1;
					}
#ifdef PRINT_TREE
				 std::cout << " " << std::bitset<6>((*it) - 1);std::flush(std::cout);
#endif

				}
#ifdef PRINT_TREE
				std::cout << std::endl << "->" << std::endl;
				std::flush(std::cout);
#endif
				++val_b;
			}

			idx_e = idx_w;
			val_e = val_w;

		}

		for(size_t j = 0; j < n_leaves; ++j)
		{
			//translate back to index
			if(*(idx_b + j) != 0)
			{
				trans_index[*(idx_b + j) - 1];
				//if(index.exists())
				if(index < index.range())
					*(idx_b + j) = index + 1;
				else
					*(idx_b + j) = 0;
			}
#ifdef PRINT_TREE
			std::cout << " " << std::bitset<6>((*(idx_b + j)) - 1); std::flush(std::cout);
#endif
			if(*(idx_b + j) != 0 && std::abs(*(val_b + j)) > eps)
				triplets.emplace_back((*(idx_b + j)) - 1, i, *(val_b + j));

		}
	*/
		for(size_t j = 0; j < n_leaves; ++j)
		{
				if(*(idx_b + j) != 0 && std::abs(*(val_b + j)) > eps)
					triplets.emplace_back((*(idx_b + j)) - 1, i, *(val_b + j));
		}

#ifdef PRINT_TREE
		std::cout << std::endl;
#endif
	}

	mat.setFromTriplets(triplets.begin(), triplets.end());

	return mat;
}

inline size_t count_bits(const idxv_type & bits)
{
        return bits_type(bits).count();
}

} //namespace

#include "operators/bit.hpp"
#include "operators/mode.hpp"
#include "operators/hopbit.hpp"
