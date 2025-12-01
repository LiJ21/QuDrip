
#pragma once
#include <iostream>
namespace qudrip {
template <typename STATE, typename OP>
class OPsiType;

template <class T>
class State {
  template <typename TS>
  friend State<TS> operator*(const Matrix& mat, const State<TS>& psi);
  template <typename TS>
  friend State<TS> operator*(const SpMatrix& mat, const State<TS>& psi);

  using index_type = T;
  index_type& idx_;
  Matrix psi_;
  size_t nt_, ndim_;
  value_type null_position;
  size_t t_;

 public:
  Matrix& data() { return psi_; }
  const Matrix& data() const { return psi_; }

  State(T& idx, size_t nt)
      : idx_(idx),
        nt_(nt),
        psi_(Matrix::Zero(idx_.range(), nt)),
        ndim_(idx_.range()),
        null_position(0.0) {}

  template <typename MATRIX>
  value_type eval(int tstp, MATRIX&& mat) const {
    if (mat.cols() != idx_.range()) {
      std::cout
          << "State: evaluation failure due to dimension mismatch. Range = "
          << idx_.range() << " and Matrix col = " << mat.cols() << std::endl;
    }
    value_type result;
    result = (psi_.block(0, tstp, idx_.range(), 1).adjoint() * mat *
              psi_.block(0, tstp, idx_.range(), 1))(0, 0);
    return result;
  }

  State& operator=(const State& s) {
    assert(&idx_ == &s.idx_);
    psi_.col(t_) = s.psi_.col(s.t());
    return *this;
  }

  value_type& operator[](int tstp) {
    t_ = tstp;
    if (idx_ < idx_.range())
      return psi_(idxv_type(idx_), tstp);
    else {
      null_position = 0.0;
      return null_position;
    }
  }

  const value_type& operator[](int tstp) const {
    return const_cast<const value_type&>(
        const_cast<const State*>(this)->operator[](tstp));
  }

  State<T>& operator()(int tstp) {
    t_ = tstp;
    return *this;
  }

  auto mat() { return psi_.block(0, t_, idx_.range(), 1); }

  const auto mat() const { return psi_.block(0, t_, idx_.range(), 1); }

  template <typename OPsi_type, typename = typename OPsi_type::isOPsi>
  State<T>& operator=(OPsi_type&& OPsi) {
    OPsi.op.map(*this, t_, OPsi.psi, OPsi.psi.t());
    return *this;
  }

  index_type& get_index() const { return idx_; }

  int ndim() const { return idx_.range(); }
  int nt() const { return nt_; }
  int t() const { return t_; }

  State<T>& operator*=(value_type alpha) {
    psi_.col(t_) *= alpha;
    return *this;
  }
  State<T>& operator/=(value_type alpha) {
    psi_.col(t_) /= alpha;
    return *this;
  }
};
//===========================================================================

template <typename T>
State<T> getState(T& index, int nt) {
  return State<T>(index, nt);
}

template <typename T>
State<T> operator*(const Matrix& mat, const State<T>& psi) {
  State<T> new_psi(psi);
  new_psi.psi_ = mat * psi.psi_;
  return new_psi;
}

template <typename T>
State<T> operator*(const SpMatrix& mat, const State<T>& psi) {
  State<T> new_psi(psi);
  new_psi.psi_ = mat * psi.psi_;
  return new_psi;
}

//===========================================================================
template <typename IDX_T, typename IDX,
          typename = typename std::remove_reference_t<IDX>::is_index>
value_type Trace(State<IDX_T>& rho, int tstp, IDX&& idxL, IDX&& idxR) {
  value_type trace = 0.0;
  for (auto i : range(idxL.range())) {
    idxL[i];
    idxR[i];
    trace += rho[tstp];
  }
  return trace;
}

//===========================================================================
template <typename T>
std::ostream& operator<<(std::ostream& os, const State<T>& psi) {
  os << psi.mat();
  return os;
}
}  // namespace qudrip
