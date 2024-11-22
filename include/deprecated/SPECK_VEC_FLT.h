#ifndef SPECK_VEC_FLT_H
#define SPECK_VEC_FLT_H

//
// This class serves as the base class of 1D, 2D, and 3D SPECK algorithm on floats.
//

#include "CDF97.h"
#include "Conditioner.h"
#include "Outlier_Coder.h"
#include "SPECK_INT.h"
#include "qoi/QoIInfo.hpp"
#include "zstd_lossless.h"

#include <variant>

namespace sperr {

class SPECK_VEC_FLT {
 public:
  //
  // Virtual Destructor
  //
  virtual ~SPECK_VEC_FLT() = default;

  //
  // Input
  //
  // Accept incoming data: copy from a raw memory block.
  // Note: `len` is the number of values.
  //template <typename T>
  //void copy_data(const T* p1,const T* p2,const T* p3, size_t len);

  // Accept incoming data: take ownership of a memory block
  void take_data(std::array<std::vector<double>,3>&&);

  // Use an encoded bitstream
  // Note: `len` is the number of bytes.
 // virtual auto use_bitstream(const void* p1, size_t len1,const void* p2, size_t len2,const void* p3, size_t len3) -> RTNType;

  //
  // Output
  //
  void append_encoded_bitstream(std::array<vec8_type,3>& buf) const;
  //auto view_decoded_data() const -> const std::array<vecd_type,3>&;
  //auto view_hierarchy() const -> const std::array<std::vector<vecd_type>,3>&;
  //auto release_decoded_data() -> std::array<vecd_type,3>&&;
  //auto release_hierarchy() -> std::array<std::vector<vecd_type>,3>&&;

  //
  // General configuration and info.
  //
  void set_psnr(std::vector<double,3> psnr);
  void set_tolerance(std::vector<double,3> tol);
  void set_bitrate(std::vector<double,3> bpp);
  void set_dims(dims_type);
  void set_qoi(std::shared_ptr<QoZ::concepts::QoIInterface<double> >);
  void set_qoi_tol(double);
  void set_qoi_block_size(int);
  std::shared_ptr<QoZ::concepts::QoIInterface<double> > get_qoi();
  auto integer_len() const -> size_t;//??

  void block_qoi_outlier_correction();

#ifdef EXPERIMENTING
  //void set_direct_q(double q);
#endif

  //
  // Actions
  //
  auto compress() -> RTNType;
  //auto decompress(bool multi_res = false) -> RTNType;

 protected:
  std::vector<UINTType,3> m_uint_flag = UINTType::UINT64;
  //std::vector<bool,3> m_has_outlier = {false,false,false};           // encoding (PWE mode) and decoding
  //std::vector<bool,3> m_has_lossless = {false,false,false};           // encoding (PWE mode) and decoding
  CompMode m_mode = CompMode::Unknown;  // encoding only
  std::vector<double,3> m_q = {0.0,0.0,0.0};                     // encoding and decoding
  std::vector<double,3> m_quality = {0.0,0.0,0.0};               // encoding only, represent either PSNR, PWE, or BPP.
  std::array<vecd_type,3> m_vals_orig;                // encoding only (PWE mode)
  dims_type m_dims = {0, 0, 0};
  std::array<vecd_type,3> m_vals_d;
  std::array<condi_type,3> m_condi_bitstream;
  std::array<Bitmask,3> m_sign_array;
  std::array<std::vector<vecd_type>,3> m_hierarchy;  // multi-resolution decoding

  std::array<CDF97,3> m_cdf;
  Conditioner m_conditioner;
  Outlier_Coder m_out_coder;
  std::shared_ptr<QoZ::concepts::QoIInterface<double> > qoi = nullptr;
  double qoi_tol = 0.0;
  int qoi_block_size = 1;
  Lossless_zstd zstd_encoder;



  std::array<std::variant<std::vector<uint8_t>,
               std::vector<uint16_t>,
               std::vector<uint32_t>,
               std::vector<uint64_t>>,3 >
      m_vals_ui;

  std::variant<std::unique_ptr<SPECK_INT<uint8_t>>,
               std::unique_ptr<SPECK_INT<uint16_t>>,
               std::unique_ptr<SPECK_INT<uint32_t>>,
               std::unique_ptr<SPECK_INT<uint64_t>>>
      m_encoder, m_decoder;

  // Instantiate `m_vals_ui` based on the chosen integer length.
  void m_instantiate_int_vec();

  // Derived classes instantiate the correct `m_encoder` and `m_decoder` depending on
  // 3D/2D/1D classes, and on the integer length in use.
  virtual void m_instantiate_encoder() = 0;
  virtual void m_instantiate_decoder() = 0;

  // Both wavelet transforms operate on `m_vals_d`.
  virtual void m_wavelet_xform() = 0;
  virtual void m_inverse_wavelet_xform(bool multi_res) = 0;

  // This base class provides two midtread quantization implementations.
  //    Quantization reads from `m_vals_d`, and writes to `m_vals_ui` and `m_sign_array`.
  //    Inverse quantization reads from `m_vals_ui` and `m_sign_array`, and writes to `m_vals_d`.
  auto m_midtread_quantize(int i) -> RTNType;
  void m_midtread_inv_quantize(int i);

  // Estimate MSE assuming midtread quantization strategy.
  auto m_estimate_mse_midtread(double q, int i) const -> double;

  // The meaning of inputs `param` and `high_prec` differ depending on the compression mode:
  //    - PWE:  no input is used; they can be anything;
  //    - PSNR: `param` must be the data range of the original input; `high_prec` is not used;
  //    - Rate: `param` must be the biggest magnitude of transformed wavelet coefficients;
  //            `high_prec` should be false at first, and true if not enough bits are produced.
  auto m_estimate_q(double param, bool high_prec) const -> double;
};

};  // namespace sperr

#endif
