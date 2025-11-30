//------------------------------------------------------------------

#pragma once
namespace qudrip {
//=========================================================================
template <typename T>
class hopBitOp : public elOpBase {
  using qbits_type = QbitsIndex<T>;
  qbits_type& index_;
  size_t site1_, site2_;
  const int string_;

 public:
  hopBitOp(qbits_type& index, int string)
      : elOpBase(2, 2), index_(index), string_(string) {}

  Operator operator()(size_t site2, size_t site1) {
    site1_ = site1;
    site2_ = site2;
    return Operator(
        std::static_pointer_cast<elOpBase>(std::make_shared<hopBitOp>(*this)));
  }
  virtual void branch(idx_it_type&, val_it_type&, idx_it_type&,
                      val_it_type&) override;
  virtual void feed_idx(idx_size_t idx) override { index_[idx]; }

  virtual idx_size_t get_idx() override { return index_; }
};

//==========================================================================
template <typename T = idx_size_t>
auto getHopFermiGate(QbitsIndex<T>& index) {
  const short fermion = -1;
  return hopBitOp<T>(index, fermion);
}
template <typename T = idx_size_t>
auto getHopBoseGate(QbitsIndex<T>& index) {
  const short boson = 1;
  return hopBitOp<T>(index, boson);
}

//==========================================================================
template <typename T>
void hopBitOp<T>::branch(idx_it_type& idx_b, val_it_type& val_b,
                         idx_it_type& idx_w, val_it_type& val_w) {
  int i = 0;
  value_type val0;
  idxv_type idx0;
  if (*idx_b != 0) {
    bits_type bits(*idx_b - 1);  // get the bitset representation

    if (site1_ == site2_) {
      idx0 = bit_convert<idxv_type>(bits) + 1;
      val0 = static_cast<double>(bits[site1_]) * (*val_b);
    } else if (bits[site1_] == 1 && bits[site2_] == 0) {
      // std::cout << " Hop bit: " << site1_ << "," << site2_ << " for " <<
      // std::bitset<6>(*idx_b - 1) << " yield counts "; std::cout <<
      // bits[site1_] << "," << bits[site2_] << std::endl;
      auto jw_set1 =
          bits_type(bit_convert<idxv_type>(bits) >> site1_ + 1).count();
      bits[site1_] = 0;

      auto jw_set2 =
          bits_type(bit_convert<idxv_type>(bits) >> site2_ + 1).count();
      bits[site2_] = 1;
      //	std::cout << jw_set1 << "," << jw_set2 << std::endl;

      idx0 = bit_convert<idxv_type>(bits) + 1;
      value_type jw = (jw_set1 % 2 ? string_ : 1) * (jw_set2 % 2 ? string_ : 1);
      val0 = jw * (*val_b);
    } else {
      idx0 = 0;
      val0 = 0.0;
    }
  } else {
    idx0 = 0;
    val0 = 0.0;
  }

  *idx_w = idx0;
  *val_w = val0;
  ++val_w;
  ++idx_w;
}
}  // namespace qudrip
