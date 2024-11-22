#ifndef SZ3_QOI_INFO
#define SZ3_QOI_INFO

#include "QoI.hpp"

#include "FXYZ.hpp"

#include <vector>

namespace QoZ {

    template<class T >
    std::shared_ptr<concepts::QoIInterface<T>> GetQOI(int idx, double qoiEB, double absErrorBound, std::string qoi_string = "x^2+y^2+z^2"){
        switch(idx){
           
            case 1:
                return std::make_shared<QoZ::QoI_FXYZ<T>>(qoiEB,absErrorBound,qoi_string);
            
        }
        return NULL;
    }

}
#endif