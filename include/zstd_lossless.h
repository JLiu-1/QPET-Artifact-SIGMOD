#ifndef ZSTD_LOSSLESS
#define ZSTD_LOSSLESS

//
// This is the base class of 1D, 2D, and 3D integer SPECK implementations.
//

#include "sperr_helper.h"
#include "Bitmask.h"
#include "Bitstream.h"
#include "zstd.h"
#include <vector>
namespace sperr {

class Lossless_zstd {

 public:
  // Constructor and destructor
  Lossless_zstd(){}
  
  template<class T>
  void encode(const std::vector<T> & data){

    size_t dataLength = data.size()*sizeof(T);

    size_t estimatedCompressedSize = (dataLength < 100 ? 200 : size_t(dataLength * 1.2));
    uint8_t *compressBytes = new uint8_t[estimatedCompressedSize];
    uint8_t *compressBytesPos = compressBytes;
    memcpy(compressBytesPos,&dataLength,sizeof(dataLength));
    compressBytesPos+=sizeof(dataLength);

    byte_size = ZSTD_compress(compressBytesPos, estimatedCompressedSize, data.data(), dataLength,
                              compression_level);
    byte_size+= sizeof(dataLength);
    //std::cout<<"bs c: "<<byte_size<<std::endl;
    m_bit_buffer.parse_bitstream(compressBytes, byte_size);
    delete []compressBytes;


  }
  uint8_t* decode(){
    std::cout<<"start decoding: "<<std::endl;
    uint8_t *dataPos = new uint8_t[byte_size+10];
    m_bit_buffer.write_bitstream(dataPos, byte_size*8); 
    std::cout<<"written bitstream "<<std::endl;
    size_t dataLength = 0;
    size_t compressedSize = byte_size;
    memcpy(&dataLength,dataPos,sizeof(dataLength));
    std::cout<<"dl d: "<<dataLength<<std::endl;
    dataPos += sizeof (dataLength);
    compressedSize -= sizeof (dataLength);
    


    uint8_t* oriData = new uint8_t [dataLength];
    ZSTD_decompress(oriData, dataLength, dataPos, compressedSize);
    delete[]dataPos;
    return oriData;
  }

  // Input
  void use_bitstream(const void* pp,size_t &pos){

    const auto* const p = static_cast<const uint8_t*>(pp);
    std::memcpy(&byte_size, p, sizeof(byte_size));
    std::cout<<"bs d: "<<byte_size<<std::endl;
    m_bit_buffer.parse_bitstream(p +sizeof(byte_size), byte_size);
    pos += byte_size+sizeof(byte_size);

  }
  size_t encoded_bitstream_len() const{
    return byte_size+sizeof(byte_size);
  }
  // Output
  //auto encoded_bitstream_len() const -> size_t;
  void append_encoded_bitstream(vec8_type& buffer) const{
    const auto app_size = encoded_bitstream_len();
    const auto orig_size = buffer.size();
    buffer.resize(orig_size + app_size);
    auto* const ptr = buffer.data() + orig_size;

    // Step 2: fill header
    size_t pos = 0;
    std::memcpy(ptr + pos, &byte_size, sizeof(byte_size));
    pos += sizeof(byte_size);

    // Step 3: assemble the right amount of bits into bytes.
    // See discussion on the number of bits to pack in function `encoded_bitstream_len()`.
    m_bit_buffer.write_bitstream(ptr + pos, byte_size);
  }


 protected:

  size_t byte_size = 0;
  Bitstream m_bit_buffer;
  int compression_level = 3;
};

};  // namespace sperr

#endif
