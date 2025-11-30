//------------------------------------------------------------------
#pragma once
namespace qudrip {
//=========================================================================

template <typename T>
class modeOp : public elOpBase {
  using mode_type = SingleModeIndex<T>;
  mode_type& index_;
  int ndim_;

 public:
  modeOp(mode_type& index)
      : elOpBase(index.range(), index.range()),
        index_(index),
        ndim_(index.range()) {}

  Operator operator()(size_t site) {
    return Operator(
        std::static_pointer_cast<elOpBase>(std::make_shared<modeOp>(*this)));
  }
  virtual void branch(idx_it_type&, val_it_type&, idx_it_type&,
                      val_it_type&) override;

  virtual void feed_idx(idx_size_t idx) override { index_[idx]; }

  virtual idx_size_t get_idx() override { return index_; }
};

//==========================================================================
template <typename T = idx_size_t>
auto getModeOp(SingleModeIndex<T>& index) {
  return modeOp<T>(index);
}
//==========================================================================

enum class Mode { Upper = 1, Lower = -1 };
//==========================================================================================
Matrix LadderMatrix(size_t ndim, Mode mode = Mode::Upper) {
  Matrix mat(ndim, ndim);
  mat.setZero();
  if (mode == Mode::Upper) {
    for (auto i : range(ndim - 1)) {
      mat(i + 1, i) = std::sqrt(i + 1);
    }
  } else if (mode == Mode::Lower) {
    for (auto i : range(ndim - 1)) {
      mat(i, i + 1) = std::sqrt(i + 1);
    }
  }

  return mat;
}

//==========================================================================================
Matrix CoherentProjector(size_t ndim, value_type amp) {
  Matrix astate(ndim, 1);

  astate(0, 0) = 1.;
  for (int i = 1; i < ndim; ++i) {
    astate(i, 0) = astate(i - 1, 0) * amp / std::sqrt((double)i);
  }
  astate *= std::exp(-std::abs(amp) * std::abs(amp) / 2.0);

  Matrix mat = astate * astate.adjoint();

  return mat;
}

//==========================================================================================
Matrix FockProjector(size_t ndim, size_t idx) {
  Matrix mat(ndim, ndim);
  mat.setZero();
  mat(idx, idx) = 1.;

  return mat;
}

//==========================================================================================
template <typename T>
void modeOp<T>::branch(idx_it_type& idx_b, val_it_type& val_b,
                       idx_it_type& idx_w, val_it_type& val_w) {
  std::vector<value_type> val(ndim_);
  std::vector<idxv_type> idx(ndim_);
  if (*idx_b != 0) {
    // auto & U = U_[*idx_b - 1];
    for (auto i = 0; i < ndim_; ++i) {
      // val[i] = U(i,0) * (*val_b);
      val[i] = U_[(*idx_b - 1) * elOpBase::n_output + i] * (*val_b);
      idx[i] = i;
    }
  } else {
    for (auto i = 0; i < ndim_; ++i) {
      val[i] = 0.0;
      idx[i] = -1;
    }
  }

  for (auto i = 0; i < ndim_; ++i) {
    *idx_w = idx[i] + 1;
    *val_w = val[i];
    ++val_w;
    ++idx_w;
  }
}
}  // namespace qudrip
