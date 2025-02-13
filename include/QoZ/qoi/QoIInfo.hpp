#ifndef SZ3_QOI_INFO
#define SZ3_QOI_INFO

#include "QoI.hpp"
#include "XSquare.hpp"
#include "LogX.hpp"
#include "RegionalAverage.hpp"
#include "RegionalAverageOfSquare.hpp"
#include "Isoline.hpp"
#include "MultiQoIs.hpp"
#include <vector>
#include "XCubic.hpp"
#include "XSqrt.hpp"
#include "XExp.hpp"
#include "XLin.hpp"

namespace QoZ {

    template<class T, QoZ::uint N>
    std::shared_ptr<concepts::QoIInterface<T, N>> GetQOI(const Config &conf){
        switch(conf.qoi){
            case 1:
                return std::make_shared<QoZ::QoI_X_Square<T, N>>(conf.qoiEB, conf.absErrorBound);
            case 2:
                return std::make_shared<QoZ::QoI_Log_X<T, N>>(conf.qoiEB, conf.absErrorBound);
            case 3:{
                if(!conf.lorenzo && !conf.lorenzo2) return std::make_shared<QoZ::QoI_RegionalAverageOfSquareInterp<T, N>>(conf.qoiEB, conf.absErrorBound, conf.qoiRegionSize, conf.dims);
                else return std::make_shared<QoZ::QoI_RegionalAverageOfSquare<T, N>>(conf.qoiEB, conf.absErrorBound);
                // return std::make_shared<SZ::QoI_RegionalAverage<T, N>>(conf.qoiEB, conf.absErrorBound);
            }
            case 4:{
            	std::vector<T> values;
            	for(int i=0; i<conf.isovalues.size(); i++){
            		values.push_back(conf.isovalues[i]);
            	}
                return std::make_shared<QoZ::QoI_Isoline<T, N>>(conf.dims, values, conf.absErrorBound);            	
            }
            case 5:{
            	// x^2 + log x
            	std::vector<std::shared_ptr<concepts::QoIInterface< T, N>>> qois;
            	qois.push_back(std::make_shared<QoZ::QoI_X_Square<T, N>>(conf.qoiEBs[0], conf.absErrorBound));
            	qois.push_back(std::make_shared<QoZ::QoI_Log_X<T, N>>(conf.qoiEBs[1], conf.absErrorBound));
                return std::make_shared<QoZ::QoI_MultiQoIs<T, N>>(qois);            	
            }
            case 6:{
            	// x^2 + isoline
            	std::vector<std::shared_ptr<concepts::QoIInterface< T, N>>> qois;
            	qois.push_back(std::make_shared<QoZ::QoI_X_Square<T, N>>(conf.qoiEBs[0], conf.absErrorBound));
                std::vector<T> values;
                for(int i=0; i<conf.isovalues.size(); i++){
                    values.push_back(conf.isovalues[i]);
                }
                qois.push_back(std::make_shared<QoZ::QoI_Isoline<T, N>>(conf.dims, values, conf.absErrorBound));
                return std::make_shared<QoZ::QoI_MultiQoIs<T, N>>(qois);            	
            }
            case 7:{
                // log x + isoline
                std::vector<std::shared_ptr<concepts::QoIInterface< T, N>>> qois;
                qois.push_back(std::make_shared<QoZ::QoI_Log_X<T, N>>(conf.qoiEBs[0], conf.absErrorBound));
                std::vector<T> values;
                for(int i=0; i<conf.isovalues.size(); i++){
                    values.push_back(conf.isovalues[i]);
                }
                qois.push_back(std::make_shared<QoZ::QoI_Isoline<T, N>>(conf.dims, values, conf.absErrorBound));
                return std::make_shared<QoZ::QoI_MultiQoIs<T, N>>(qois);             
            }
            case 8:{
            	// x^2 + log x + isoline
            	std::vector<std::shared_ptr<concepts::QoIInterface< T, N>>> qois;
            	qois.push_back(std::make_shared<QoZ::QoI_X_Square<T, N>>(conf.qoiEBs[0], conf.absErrorBound));
                qois.push_back(std::make_shared<QoZ::QoI_Log_X<T, N>>(conf.qoiEBs[1], conf.absErrorBound));
            	std::vector<T> values;
            	for(int i=0; i<conf.isovalues.size(); i++){
            		values.push_back(conf.isovalues[i]);
            	}
                qois.push_back(std::make_shared<QoZ::QoI_Isoline<T, N>>(conf.dims, values, conf.absErrorBound));
                return std::make_shared<QoZ::QoI_MultiQoIs<T, N>>(qois);            	
            }
            case 13:{
                return std::make_shared<QoZ::QoI_X_Cubic<T, N>>(conf.qoiEB, conf.absErrorBound);

            }
            case 10:{
                return std::make_shared<QoZ::QoI_X_Sqrt<T, N>>(conf.qoiEB, conf.absErrorBound);
            }
            case 11:
                return std::make_shared<QoZ::QoI_X_Lin<T, N>>(conf.qoiEB, conf.absErrorBound);
            case 12:
                return std::make_shared<QoZ::QoI_X_Exp<T, N>>(conf.qoiEB, conf.absErrorBound,conf.qoiLogBase);
            case 9:{
                if(!conf.lorenzo && !conf.lorenzo2) return std::make_shared<QoZ::QoI_RegionalAverageInterp<T, N>>(conf.qoiEB, conf.absErrorBound, conf.qoiRegionSize, conf.dims);
                else return std::make_shared<QoZ::QoI_RegionalAverage<T, N>>(conf.qoiEB, conf.absErrorBound);
            }
        }
        return NULL;
    }

}
#endif