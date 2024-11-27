#ifndef SZ3_QOI_INFO
#define SZ3_QOI_INFO

#include "QoI.hpp"

#include "FX.hpp"
#include "FX_abs.hpp"
#include "FX_P.hpp"
#include "XCubic.hpp"
#include <vector>

namespace QoZ {

    template<class T>
    std::shared_ptr<concepts::QoIInterface<T> > GetQOI(int idx, double qoiEB, double absErrorBound, std::string qoi_string = "x^2", bool isolated = false, double threshold = 0.0, std::string qoi_string_2 = "x^2" ){
        switch(idx){
            case 1:
                return std::make_shared<QoZ::QoI_FX<T>>(qoiEB, absErrorBound, qoi_string, isolated, threshold);
            case 2:
                return std::make_shared<QoZ::QoI_FX_P<T>>(qoiEB, absErrorBound, qoi_string, qoi_string_2, threshold, isolated);
            case 3:
                return std::make_shared<QoZ::QoI_FX_ABS<T>>(qoiEB, absErrorBound, qoi_string, isolated, threshold);
            case 3:
                return std::make_shared<QoZ::QoI_FX_ABS<T>>(qoiEB, absErrorBound);
        }
        return NULL;
    }

}
#endif