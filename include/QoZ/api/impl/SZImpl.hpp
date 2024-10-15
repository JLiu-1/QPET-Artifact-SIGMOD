#ifndef SZ3_IMPL_SZ_HPP
#define SZ3_IMPL_SZ_HPP

#include "QoZ/def.hpp"
#include "QoZ/api/impl/SZDispatcher.hpp"
//#include "QoZ/api/impl/SZImplOMP.hpp"
#include <cmath>


template<class T, QoZ::uint N>
std::array<char *,3>SZ_compress_impl(std::array<QoZ::Config,3> &confs, std::array<T *,3> &data, std::array<size_t,3> &outSizes) {
    /*
#ifndef _OPENMP
    conf.openmp=false;
#endif
    if (conf.openmp) {
        return SZ_compress_OMP<T, N>(confs, data, outSizes);
    } else {

       */
        auto output=SZ_compress_dispatcher<T, N>(confs, data, outSizes);
     
       
        return output;
    //}
}


template<class T, QoZ::uint N>
void SZ_decompress_impl(std::array<QoZ::Config,3> &conf, std::array<char *,3> &cmpData, std::array<size_t,3> &cmpSizes, std::array<T *,3>&decData) {
    /*
#ifndef _OPENMP
    conf.openmp=false;
#endif
   
    if (conf.openmp) {
        SZ_decompress_OMP<T, N>(conf, cmpData, cmpSize, decData);
    } else {*/
        SZ_decompress_dispatcher<T, N>(confs, cmpData, cmpSizes, decData);
    //}
}

#endif