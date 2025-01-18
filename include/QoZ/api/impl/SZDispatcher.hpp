#ifndef SZ3_IMPL_SZDISPATCHER_HPP
#define SZ3_IMPL_SZDISPATCHER_HPP

#include "QoZ/utils/MemoryUtil.hpp"
#include "QoZ/utils/Statistic.hpp"
#include "QoZ/utils/Config.hpp"
#include "QoZ/api/impl/SZInterp.hpp"
#include "QoZ/api/impl/SZLorenzoReg.hpp"
#include <cmath>


template<class T, QoZ::uint N>
std::array<char *,3>SZ_compress_dispatcher(std::array<QoZ::Config,3> &confs, std::array<T *,3>&data, std::array<size_t,3> &outSizes) {

    assert(N == conf.N);
    for(auto i:{0,1,2})
        QoZ::calAbsErrorBound(confs[i], data[i]);

    std::array<char *,3> cmpData;
    
    if (conf.cmprAlgo == QoZ::ALGO_LORENZO_REG) {
        cmpData = (char *) SZ_compress_LorenzoReg_Vec<T, N>(confs, data, outSizes);
    } else
     if (confs[0].cmprAlgo == QoZ::ALGO_INTERP) {
        cmpData = SZ_compress_Interp<T, N>(confs, data, outSizes);
    } else if (confs[0].cmprAlgo == QoZ::ALGO_INTERP_LORENZO) {
        cmpData = SZ_compress_Interp_lorenzo<T, N>(confs, data, outSizes);
    }
    
    /*
    else if (conf.cmprAlgo == QoZ::ALGO_NEWINTERP) {
        cmpData = (char *) SZ_compress_NewInterp<T, N>(conf, data, outSize);
    }
    */
    /*
    else if (conf.cmprAlgo == QoZ::ALGO_INTERP_BLOCKED) {

        cmpData = (char *) SZ_compress_Interp_blocked<T, N>(conf, data, outSize);
    }
    */
  
    return cmpData;
}


template<class T, QoZ::uint N>
void SZ_decompress_dispatcher(std::array<QoZ::Config,3> &confs, std::array<char *,3> &cmpData, std::array<size_t,3> &cmpSizes, std::array<T *,3> &decData) {
    
    if (conf.cmprAlgo == QoZ::ALGO_LORENZO_REG) {
        SZ_decompress_LorenzoReg_Vec<T, N>(confs, cmpData, cmpSizes, decData);
    } else if (confs[0].cmprAlgo == QoZ::ALGO_INTERP) {
        SZ_decompress_Interp<T, N>(confs, cmpData, cmpSizes, decData);
    } else {
        printf("SZ_decompress_dispatcher, Method not supported\n");
        exit(0);
    }
}

#endif