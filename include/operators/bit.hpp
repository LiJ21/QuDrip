//------------------------------------------------------------------

#pragma once
namespace qudrip {
//=========================================================================
template <typename T>
class bitOp : public elOpBase {
  using qbits_type = QbitsIndex<T>;
  qbits_type& index_;
  size_t site_;
  const int string_;

 public:
  bitOp(qbits_type& index, int string)
      : elOpBase(2, 2), index_(index), string_(string) {}

  Operator operator()(size_t site) {
    (void)detail::checked_site_mask(site, index_.range() - 1);
    site_ = site;
    return Operator(
        std::static_pointer_cast<elOpBase>(std::make_shared<bitOp>(*this)));
  }
  virtual void branch(idx_it_type&, val_it_type&, idx_it_type&,
                      val_it_type&) override;
  virtual void feed_idx(idx_size_t idx) override { index_[idx]; }

  virtual idx_size_t get_idx() override { return index_; }

  bool describe_factor(detail::FactorDescriptor& factor) const noexcept override {
    factor.kind = detail::FactorKind::BitMatrix;
    factor.primitive_identity = std::addressof(index_);
    factor.input_range = elOpBase::n_input;
    factor.output_range = elOpBase::n_output;
    factor.matrix = std::span<const value_type>(elOpBase::U_);
    factor.site = site_;
    factor.site2 = site_;
    factor.string_value = string_;
    return true;
  }
};

//==========================================================================
template <typename T = idx_size_t>
auto getFermiGate(QbitsIndex<T>& index) {
  const short fermion = -1;
  return bitOp<T>(index, fermion);
}
template <typename T = idx_size_t>
auto getBoseGate(QbitsIndex<T>& index) {
  const short boson = 1;
  return bitOp<T>(index, boson);
}

//==========================================================================
inline const Matrix pauli_x = (Matrix(2, 2) << 0, 1, 1, 0).finished();
inline const Matrix pauli_y = (Matrix(2, 2) << 0, -II, II, 0).finished();
inline const Matrix pauli_z = (Matrix(2, 2) << 1, 0, 0, -1).finished();
inline const Matrix id = (Matrix(2, 2) << 1, 0, 0, 1).finished();

//==========================================================================
template <typename T>
void bitOp<T>::branch(idx_it_type& idx_b, val_it_type& val_b,
                      idx_it_type& idx_w, val_it_type& val_w) {
  value_type val0, val1;
  idxv_type idx0, idx1;
  if (*idx_b != 0) {
    bits_type bits(*idx_b - 1);  // get the bitset representation
    // auto & U = U_[bits[site_]];
    auto prev_set =
        bits_type(bit_convert<idxv_type>(bits) >> (site_ + 1)).count();
    value_type jw = prev_set % 2 ? string_ : 1;
    // val0 = jw * U(0,0) * (*val_b);
    // val1 = jw * U(1,0) * (*val_b);
    val0 = jw * U_[bits[site_] * 2] * (*val_b);
    val1 = jw * U_[bits[site_] * 2 + 1] * (*val_b);
    bits[site_] = 0;
    idx0 = bit_convert<idxv_type>(bits) + 1;
    bits[site_] = 1;
    idx1 = bit_convert<idxv_type>(bits) + 1;
  } else {
    idx0 = 0;
    idx1 = 0;
    val0 = 0.0;
    val1 = 0.0;
  }

  *idx_w = idx0;
  *val_w = val0;
  ++val_w;
  ++idx_w;

  *idx_w = idx1;
  *val_w = val1;
  ++val_w;
  ++idx_w;
}
}  // namespace qudrip
