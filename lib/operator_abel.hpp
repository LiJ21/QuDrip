
#pragma once
//define addition and subtraction between operators
namespace qudrip
{

template<typename OP, typename STATE>
class OPsiType;

class OperatorSum
{
	using ops_type = std::vector<Operator>;
	ops_type ops_;
	std::vector<value_type> coeffs_;

	public:

	int size() const {return ops_.size();}

	OperatorSum(int i = 2)
	{
		ops_.reserve(i);	
		coeffs_.reserve(i);	
	}

	OperatorSum(const Operator & op, value_type coeff = 1.0):
		ops_(1, op), coeffs_(1, coeff)
	{}

	template<typename IDX, typename = typename std::remove_reference_t<IDX>::is_index>
	SpMatrix operator>>(IDX && index)
	{
		SpMatrix result(index.range(), index.range());
		for(int i = 0; i < ops_.size(); ++i)
		{
			result += coeffs_[i] * (ops_[i] >> std::forward<IDX>(index));
		}

		return result.eval();
	}

	template<typename index_type>
	void map_acc(State<index_type> & to_psi, int t1, State<index_type> & from_psi, int t2, value_type coeff = 1.0)
	{
		for(int i = 0; i < ops_.size(); ++i)
		{
			ops_[i].map_acc(to_psi, t1, from_psi, t2, coeffs_[i] * coeff);
		}
	}

	template<typename index_type>
	void map(State<index_type> & to_psi, int t1, State<index_type> & from_psi, int t2, value_type coeff = 1.0)
	{
		to_psi(t2).mat().setZero();
		map_acc(to_psi, t1, from_psi, t2, coeff);
	}

	OperatorSum operator*(const Operator & op) const
	{
		OperatorSum newOS(*this);
		for(auto & op_ : newOS.ops_)
		{
			op_ = op_ * op;
		}

		return newOS;
	}

	template<typename index_type>
	auto operator*(State<index_type> & psi)
	{
		return OPsiType<OperatorSum, State<index_type>>(*this, psi);
	}

	OperatorSum operator*(value_type mult) const
	{
		OperatorSum newOS(*this);
		for(auto & c : newOS.coeffs_)
		{
			c *= mult;
		}
		return newOS;
	}

	OperatorSum operator*(const OperatorSum & ops) const
	{
		OperatorSum newOS(ops_.size() * ops.ops_.size());
		for(auto i : range(ops_.size()))
		{
			for(auto j : range(ops.ops_.size()))
			{
				newOS.append(ops_[i] * ops.ops_[j], coeffs_[i] * ops.coeffs_[j]);
			}
		}
		return newOS;
	}

	OperatorSum operator+(const Operator & op) const
	{
		OperatorSum newOS(*this);
		newOS.append(op, 1);
		return newOS;
	}

	OperatorSum operator+(const OperatorSum & ops) const
	{
		OperatorSum newOS(*this);
		newOS.ops_.insert(newOS.ops_.end(), ops.ops_.begin(), ops.ops_.end());
		newOS.coeffs_.insert(newOS.coeffs_.end(), ops.coeffs_.begin(), ops.coeffs_.end());

		return newOS;
	}

	OperatorSum operator-(const OperatorSum & ops) const
	{
		OperatorSum newOS(*this);
		newOS.ops_.insert(newOS.ops_.end(), ops.ops_.begin(), ops.ops_.end());

		for(auto i : range(ops.size()))
		{
			newOS.coeffs_.push_back(-ops.coeffs_[i]);
		}

		return newOS;
	}

	OperatorSum operator-(const Operator & op) const
	{
		OperatorSum newOS(*this);
		newOS.append(op, -1);
		return newOS;
	}

	OperatorSum & append(const Operator & op, value_type coeff)
	{
		bool test_eq = false;
		for(auto i : range(ops_.size()))
		{
			if(ops_[i] == op)
			{
				coeffs_[i] += coeff;
				return *this;
			}
		}

		ops_.push_back(op);
		coeffs_.push_back(coeff);
		return *this;
	}
	
};

//=============================================================
OperatorSum operator*(value_type coeff, const Operator & op)
{
	return OperatorSum(op, coeff);
}

OperatorSum operator*(const Operator & op, value_type coeff)
{
	return OperatorSum(op, coeff);
}

OperatorSum operator*(value_type coeff, const OperatorSum & ops)
{
	auto newOS = ops * coeff;

	return newOS;
}

OperatorSum operator+(const Operator & op1, const Operator & op2)
{
	OperatorSum newOS(op1, 1.0);
	newOS.append(op2, 1.0);
	return newOS;
}

OperatorSum operator-(const Operator & op1, const Operator & op2)
{
	OperatorSum newOS(op1, 1.0);
	newOS.append(op2, -1.0);
	return newOS;
}

//=============================================================
OperatorSum getEmptyOperator()
{
	return OperatorSum(2);
}

}//namespace
