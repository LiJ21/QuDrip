
#pragma once
// define addition and subtraction between operators
namespace qudrip {

template <typename OP, typename STATE>
class OPsiType;

class OperatorSum {
  using ops_type = std::vector<Operator>;
  ops_type ops_;
  std::vector<value_type> coeffs_;

  struct ColumnEntry {
    idxv_type row;
    value_type value;
    size_t ordinal;
  };

  template <typename IDX>
  SpMatrix assemble_fused_triplets(IDX &index, idxv_type ndim,
                                   Eigen::Index eigen_dim) const {
    SpMatrix result(eigen_dim, eigen_dim);
    triplets_type triplets;

    const auto dimension = static_cast<size_t>(ndim);
    if (dimension != 0 &&
        ops_.size() <= triplets.max_size() / dimension) {
      triplets.reserve(ops_.size() * dimension);
    }
    for (size_t term = 0; term < ops_.size(); ++term)
      ops_[term].append_triplets(index, coeffs_[term], triplets);

    if (triplets.size() > detail::eigen_sparse_storage_limit())
      throw std::length_error(
          "operator-sum triplet count exceeds Eigen sparse storage range");
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
  }

  template <typename IDX>
  SpMatrix assemble_direct_csc(
      IDX &index, idxv_type ndim, Eigen::Index eigen_dim,
      const std::vector<detail::BoundOperatorPlan> &plans) const {
    SpMatrix result(eigen_dim, eigen_dim);
    const auto storage_limit = detail::eigen_sparse_storage_limit();

    if (ndim != 0 && ops_.size() <= storage_limit &&
        static_cast<idxv_type>(ops_.size()) <= storage_limit / ndim) {
      const auto estimate =
          static_cast<idxv_type>(ops_.size()) * ndim;
      result.reserve(detail::checked_eigen_sparse_index(
          estimate, "sparse matrix reserve exceeds Eigen index range"));
    }

    std::vector<detail::OperatorWorkspace> workspaces;
    workspaces.reserve(plans.size());
    for (const auto &plan : plans) workspaces.emplace_back(plan);

    std::vector<ColumnEntry> entries;
    if (ops_.size() <= entries.max_size()) entries.reserve(ops_.size());

    idxv_type cumulative_nnz = 0;
    for (idxv_type column = 0; column < ndim; ++column) {
      entries.clear();
      const auto source_raw = index.raw_at(column);

      for (size_t term = 0; term < plans.size(); ++term) {
        detail::evaluate_source(
            plans[term], index, column, source_raw, value_type(1.0),
            workspaces[term],
            [&](idxv_type row, value_type amplitude) {
              if (entries.size() == entries.max_size())
                throw std::length_error(
                    "sparse column candidate count exceeds vector capacity");
              entries.push_back(
                  {row, coeffs_[term] * amplitude, entries.size()});
            });
      }

      std::sort(entries.begin(), entries.end(),
                [](const ColumnEntry &lhs, const ColumnEntry &rhs) {
                  if (lhs.row != rhs.row) return lhs.row < rhs.row;
                  return lhs.ordinal < rhs.ordinal;
                });

      const auto eigen_column = detail::checked_eigen_sparse_index(
          column, "sparse matrix column exceeds Eigen index range");
      result.startVec(eigen_column);

      size_t cursor = 0;
      while (cursor < entries.size()) {
        const auto row = entries[cursor].row;
        value_type reduced = entries[cursor].value;
        ++cursor;
        while (cursor < entries.size() && entries[cursor].row == row) {
          reduced += entries[cursor].value;
          ++cursor;
        }

        if (cumulative_nnz == storage_limit)
          throw std::length_error(
              "sparse matrix nonzero count exceeds Eigen storage range");
        const auto eigen_row = detail::checked_eigen_sparse_index(
            row, "sparse matrix row exceeds Eigen index range");
        result.insertBack(eigen_row, eigen_column) = reduced;
        ++cumulative_nnz;
      }
    }

    result.finalize();
    return result;
  }

 public:
  size_t size() const noexcept { return ops_.size(); }

  explicit OperatorSum(size_t reserve_hint = 2) {
    ops_.reserve(reserve_hint);
    coeffs_.reserve(reserve_hint);
  }

  OperatorSum(const Operator& op, value_type coeff = 1.0)
      : ops_(1, op), coeffs_(1, coeff) {}

  template <typename IDX,
            typename = typename std::remove_reference_t<IDX>::is_index>
  SpMatrix operator>>(IDX&& index) const {
    auto& index_ref = index;
    const auto ndim = static_cast<idxv_type>(index_ref.range());
    const auto eigen_dim = detail::checked_eigen_sparse_index(
        ndim, "sparse matrix dimension exceeds Eigen index range");

    std::vector<detail::BoundOperatorPlan> plans;
    plans.reserve(ops_.size());
    bool has_legacy = false;
    for (const auto &op : ops_) {
      plans.push_back(detail::bind_operator_plan(op.ops_, index_ref));
      has_legacy |= plans.back().kind == detail::KernelKind::Legacy;
    }

    if (has_legacy)
      return assemble_fused_triplets(index_ref, ndim, eigen_dim);
    return assemble_direct_csc(index_ref, ndim, eigen_dim, plans);
  }

#ifdef QUDRIP_TESTING
  template <typename IDX,
            typename = typename std::remove_reference_t<IDX>::is_index>
  SpMatrix fused_triplets_for_test(IDX&& index) const {
    auto& index_ref = index;
    const auto ndim = static_cast<idxv_type>(index_ref.range());
    const auto eigen_dim = detail::checked_eigen_sparse_index(
        ndim, "sparse matrix dimension exceeds Eigen index range");
    return assemble_fused_triplets(index_ref, ndim, eigen_dim);
  }
#endif

  template <typename index_type>
  void map_acc(State<index_type>& to_psi, int t1,
               const State<index_type>& from_psi, int t2,
               value_type coeff = 1.0) const {
    for (size_t i = 0; i < ops_.size(); ++i) {
      ops_[i].map_acc(to_psi, t1, from_psi, t2, coeffs_[i] * coeff);
    }
  }

  template <typename index_type>
  void map(State<index_type>& to_psi, int t1,
           const State<index_type>& from_psi, int t2,
           value_type coeff = 1.0) const {
    to_psi(t1).mat().setZero();
    map_acc(to_psi, t1, from_psi, t2, coeff);
  }

  OperatorSum operator*(const Operator& op) const {
    OperatorSum newOS(*this);
    for (auto& op_ : newOS.ops_) {
      op_ = op_ * op;
    }

    return newOS;
  }

  template <typename index_type>
  auto operator*(const State<index_type>& psi) const {
    return OPsiType<OperatorSum, State<index_type>>(*this, psi);
  }

  OperatorSum operator*(value_type mult) const {
    OperatorSum newOS(*this);
    for (auto& c : newOS.coeffs_) {
      c *= mult;
    }
    return newOS;
  }

  OperatorSum operator*(const OperatorSum& ops) const {
    if (ops_.size() != 0 &&
        ops.ops_.size() >
            std::numeric_limits<size_t>::max() / ops_.size())
      throw std::length_error("operator-sum product size overflow");
    OperatorSum newOS(ops_.size() * ops.ops_.size());
    for (auto i : range(ops_.size())) {
      for (auto j : range(ops.ops_.size())) {
        newOS.append(ops_[i] * ops.ops_[j], coeffs_[i] * ops.coeffs_[j]);
      }
    }
    return newOS;
  }

  OperatorSum operator+(const Operator& op) const {
    OperatorSum newOS(*this);
    newOS.append(op, 1);
    return newOS;
  }

  OperatorSum operator+(const OperatorSum& ops) const {
    OperatorSum newOS(*this);
    newOS.ops_.insert(newOS.ops_.end(), ops.ops_.begin(), ops.ops_.end());
    newOS.coeffs_.insert(newOS.coeffs_.end(), ops.coeffs_.begin(),
                         ops.coeffs_.end());

    return newOS;
  }

  OperatorSum operator-(const OperatorSum& ops) const {
    OperatorSum newOS(*this);
    newOS.ops_.insert(newOS.ops_.end(), ops.ops_.begin(), ops.ops_.end());

    for (auto i : range(ops.size())) {
      newOS.coeffs_.push_back(-ops.coeffs_[i]);
    }

    return newOS;
  }

  OperatorSum operator-(const Operator& op) const {
    OperatorSum newOS(*this);
    newOS.append(op, -1);
    return newOS;
  }

  OperatorSum& append(const Operator& op, value_type coeff) {
    for (auto i : range(ops_.size())) {
      if (ops_[i] == op) {
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
inline OperatorSum operator*(value_type coeff, const Operator& op) {
  return OperatorSum(op, coeff);
}

inline OperatorSum operator*(const Operator& op, value_type coeff) {
  return OperatorSum(op, coeff);
}

inline OperatorSum operator*(value_type coeff, const OperatorSum& ops) {
  auto newOS = ops * coeff;

  return newOS;
}

inline OperatorSum operator+(const Operator& op1, const Operator& op2) {
  OperatorSum newOS(op1, 1.0);
  newOS.append(op2, 1.0);
  return newOS;
}

inline OperatorSum operator-(const Operator& op1, const Operator& op2) {
  OperatorSum newOS(op1, 1.0);
  newOS.append(op2, -1.0);
  return newOS;
}

//=============================================================
inline OperatorSum getEmptyOperator() { return OperatorSum(2); }

}  // namespace qudrip
