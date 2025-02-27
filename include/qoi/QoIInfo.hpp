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
    std::shared_ptr<concepts::QoIInterface<T> > GetQOI(QoIMeta &meta, double qoiEB, double absErrorBound){
        auto analytical = meta.analytical;
        auto qoi_string = meta.qoi_string;
        auto qoi_string_2 = meta.qoi_string_2;
        if(meta.qoi_id == 14){
            
            if(qoi_string == "x^2" or qoi_string == "x**2")
                meta.qoi_id = 1;
            else if(qoi_string == "logx" or qoi_string == "log(x)" or qoi_string == "Logx" or qoi_string == "Log(x)" or qoi_string == "lnx" or qoi_string == "ln(x)" or qoi_string == "Lnx" or qoi_string == "Ln(x)"){
                meta.qoi_id = 2;
                meta.qoi_base = std::exp(1.0);
            }
            else if(qoi_string == "log2(x)"  or qoi_string == "Log2(x)"  or qoi_string == "Log(x,2)"  or qoi_string == "log(x,2)"){
                meta.qoi_id = 2;
                meta.qoi_base = 2.0;
            }

            else if (qoi_string == "x^3" or qoi_string == "x**3")
                meta.qoi_id = 9;
            else if (qoi_string == "x^0.5" or qoi_string == "x**0.5" or qoi_string == "x**1/2" or qoi_string == "x^1/2" or qoi_string == "sqrt(x)" or qoi_string == "Sqrt(x)")
                meta.qoi_id = 10;
            //temp do not dispath lin beacuse it is preprocessed in compressor. Todo: dispatch lin.
            else if (qoi_string == "2^x" or qoi_string == "2**x"){
                meta.qoi_id = 12;
                meta.qoi_base= 2.0;
            }
            else if (qoi_string == "e^x" or qoi_string == "e**x" or qoi_string == "E**x" or qoi_string == "E^x"){
                meta.qoi_id = 12;
                meta.qoi_base= std::exp(1.0);
            }
            else if (qoi_string == "1/x" or qoi_string == "x**-1" or qoi_string == "x^-1")
                meta.qoi_id = 13;
            //todo: dispatch x^n
            else if (qoi_string == "|x|")
                meta.qoi_id = 19;
            else if (qoi_string == "sinx" or qoi_string == "sin(x)" or qoi_string == "Sinx" or qoi_string == "Sin(x)")
                meta.qoi_id = 22;
            else if (qoi_string == "sin10x" or qoi_string == "sin(10x)" or qoi_string == "Sin10x" or qoi_string == "Sin(10x)"){
                meta.qoi_id = 20;
                qoi_string = "11 1 0 22";
            }
            else if (qoi_string == "tanhx" or qoi_string == "tanh(x)" or qoi_string == "Tanhx" or qoi_string == "Tanh(x)")
                meta.qoi_id = 23;

        }
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
                return std::make_shared<QoZ::QoI_FX<T>>(qoiEB, absErrorBound, qoi_string, meta.isolated, meta.threshold);
            case 15:
                return std::make_shared<QoZ::QoI_FX_P<T>>(qoiEB, absErrorBound, qoi_string, qoi_string_2, meta.threshold, meta.isolated);

            case 17:
                return std::make_shared<QoZ::QoI_FX_ABS<T>>(qoiEB, absErrorBound, qoi_string, meta.isolated, meta.threshold);
            case 18:{
                if(analytical)
                    return std::make_shared<QoZ::QoI_X_Power<T>>(qoiEB, absErrorBound,meta.qoi_base);
                else
                    return std::make_shared<QoZ::QoI_X_Power_Approx<T>>(qoiEB, absErrorBound,meta.qoi_base);

            }
            case 19:
                return std::make_shared<QoZ::QoI_X_Abs<T>>(qoiEB, absErrorBound);
            case 20:
                return std::make_shared<QoZ::QoI_X_Composite<T>>(qoiEB, absErrorBound, qoi_string, analytical);
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