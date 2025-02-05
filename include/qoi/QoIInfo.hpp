#ifndef SZ3_QOI_INFO
#define SZ3_QOI_INFO

#include "QoI.hpp"
#include "XSquare.hpp"
#include "XCubic.hpp"
#include "XPower.hpp"
#include "XSin.hpp"
#include "XExp.hpp"
#include "XLin.hpp"
#include "XSqrt.hpp"
#include "XAbs.hpp"
#include "XTanh.hpp"
#include "XReciprocal.hpp"
#include "XComposite.hpp"
#include "LogX.hpp"
#include "FX.hpp"
#include "FX_abs.hpp"
#include "FX_P.hpp"
#include <vector>

namespace QoZ {

    template<class T>
    std::shared_ptr<concepts::QoIInterface<T> > GetQOI(const QoIMeta &meta, double qoiEB, double absErrorBound){
        auto analytical = meta.analytical;
        switch(meta.qoi_id){
            case 1:
                return std::make_shared<QoZ::QoI_X_Square<T>>(qoiEB, absErrorBound);
            case 2:{
                if(analytical)
                    return std::make_shared<QoZ::QoI_Log_X<T>>(qoiEB, absErrorBound, meta.qoi_base);
                else
                    return std::make_shared<QoZ::QoI_Log_X_Approx<T>>(qoiEB, absErrorBound, meta.qoi_base);
            }
            case 9:{
                if(analytical)
                    return std::make_shared<QoZ::QoI_X_Cubic<T>>(qoiEB, absErrorBound);
                else
                    return std::make_shared<QoZ::QoI_X_Cubic_Approx<T>>(qoiEB, absErrorBound);
            }       
            case 10:{
                //return std::make_shared<QoZ::QoI_X_Sin<T>>(qoiEB, absErrorBound);
                if(analytical)
                    return std::make_shared<QoZ::QoI_X_Sqrt<T>>(qoiEB, absErrorBound);
                else
                    return std::make_shared<QoZ::QoI_X_Sqrt_Approx<T>>(qoiEB, absErrorBound);
            }
            case 11:
                return std::make_shared<QoZ::QoI_X_Lin<T>>(qoiEB, absErrorBound, meta.lin_A, meta.lin_B);
            case 12:{
                if(analytical)
                    return std::make_shared<QoZ::QoI_X_Exp<T>>(qoiEB, absErrorBound,meta.qoi_base);
                else
                    return std::make_shared<QoZ::QoI_X_Exp_Approx<T>>(qoiEB, absErrorBound,meta.qoi_base);
            }
            case 13:{
                //return std::make_shared<QoZ::QoI_XLog_X<T>>(qoiEB, absErrorBound);
                if(analytical)
                    return std::make_shared<QoZ::QoI_X_Recip<T>>(qoiEB, absErrorBound);
                else
                    return std::make_shared<QoZ::QoI_X_Recip_Approx<T>>(qoiEB, absErrorBound);
            }
            case 14:
                return std::make_shared<QoZ::QoI_FX<T>>(qoiEB, absErrorBound, meta.qoi_string, meta.isolated, meta.threshold);
            case 15:
                return std::make_shared<QoZ::QoI_FX_P<T>>(qoiEB, absErrorBound, meta.qoi_string, meta.qoi_string_2, meta.threshold, meta.isolated);

            case 17:
                return std::make_shared<QoZ::QoI_FX_ABS<T>>(qoiEB, absErrorBound, meta.qoi_string, meta.isolated, meta.threshold);
            case 18:{
                if(analytical)
                    return std::make_shared<QoZ::QoI_X_Power<T>>(qoiEB, absErrorBound,meta.qoi_base);
                else
                    return std::make_shared<QoZ::QoI_X_Power_Approx<T>>(qoiEB, absErrorBound,meta.qoi_base);

            }
            case 19:
                return std::make_shared<QoZ::QoI_X_Abs<T>>(qoiEB, absErrorBound);
            case 20:
                return std::make_shared<QoZ::QoI_X_Composite<T>>(qoiEB, absErrorBound, meta.qoi_string, analytical);
            case 22:{
                if(analytical)
                    return std::make_shared<QoZ::QoI_X_Sin<T>>(qoiEB, absErrorBound);
                else
                    return std::make_shared<QoZ::QoI_X_Sin_Approx<T>>(qoiEB, absErrorBound);
            }
            case 23:{
                if(analytical)
                    return std::make_shared<QoZ::QoI_X_Tanh<T>>(qoiEB, absErrorBound);
                else
                    return std::make_shared<QoZ::QoI_X_Tanh_Approx<T>>(qoiEB, absErrorBound);
            }

        }
        return NULL;
    }

}
#endif