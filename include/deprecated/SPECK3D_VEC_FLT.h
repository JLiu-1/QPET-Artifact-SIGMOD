#ifndef SPECK3D_VEC_FLT_H
#define SPECK3D_VEC_FLT_H

#include "SPECK_VEC_FLT.h"

namespace sperr {

class SPECK3D_VEC_FLT : public SPECK_VEC_FLT {
 protected:
  void m_instantiate_encoder() override;
  void m_instantiate_decoder() override;

  void m_wavelet_xform() override;
  void m_inverse_wavelet_xform(bool) override;
};

};  // namespace sperr

#endif
