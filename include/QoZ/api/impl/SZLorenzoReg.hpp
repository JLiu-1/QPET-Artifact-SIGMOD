#ifndef SZ3_SZ_LORENZO_REG_HPP
#define SZ3_SZ_LORENZO_REG_HPP

#include "QoZ/compressor/SZGeneralCompressor.hpp"
#include "QoZ/frontend/SZFastFrontend.hpp"
#include "QoZ/frontend/SZGeneralFrontend.hpp"
#include "QoZ/frontend/SZQoIFrontend.hpp"
#include "QoZ/quantizer/IntegerQuantizer.hpp"
#include "QoZ/quantizer/QoIIntegerQuantizer.hpp"
#include "QoZ/predictor/ComposedPredictor.hpp"
#include "QoZ/predictor/LorenzoPredictor.hpp"
#include "QoZ/predictor/RegressionPredictor.hpp"
#include "QoZ/predictor/PolyRegressionPredictor.hpp"
#include "QoZ/predictor/ZeroPredictor.hpp"
#include "QoZ/encoder/QoIEncoder.hpp"
#include "QoZ/lossless/Lossless_zstd.hpp"
//#include "QoZ/qoi/XSquare.hpp"
#include "QoZ/qoi/QoIInfo.hpp"
#include "QoZ/utils/Iterator.hpp"
#include "QoZ/utils/Statistic.hpp"
#include "QoZ/utils/Extraction.hpp"
#include "QoZ/utils/QuantOptimization.hpp"
#include "QoZ/utils/Config.hpp"
#include "QoZ/def.hpp"

#include <cmath>
#include <cstdlib>
#include <memory>


template<class T, QoZ::uint N, class Quantizer, class Encoder, class Lossless>
//std::shared_ptr<QoZ::concepts::CompressorInterface<T>>
QoZ::concepts::CompressorInterface<T>*
make_lorenzo_regression_compressor(const QoZ::Config &conf, Quantizer quantizer, Encoder encoder, Lossless lossless) {
    std::vector<std::shared_ptr<QoZ::concepts::PredictorInterface<T, N>>> predictors;

    int methodCnt = (conf.lorenzo + conf.lorenzo2 + conf.regression + conf.regression2);
    int use_single_predictor = (methodCnt == 1);
    if (methodCnt == 0) {
        printf("All lorenzo and regression methods are disabled.\n");
        exit(0);
    }
    if (conf.lorenzo) {
        
        if (use_single_predictor) {
            return QoZ::make_sz_general_compressor<T, N>(
                    QoZ::make_sz_general_frontend<T, N>(conf, QoZ::LorenzoPredictor<T, N, 1>(conf.absErrorBound), quantizer),
                    encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<QoZ::LorenzoPredictor<T, N, 1>>(conf.absErrorBound));
        }
    }
    if (conf.lorenzo2) {
       
        if (use_single_predictor) {
            return QoZ::make_sz_general_compressor<T, N>(
                    QoZ::make_sz_general_frontend<T, N>(conf, QoZ::LorenzoPredictor<T, N, 2>(conf.absErrorBound), quantizer),
                    encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<QoZ::LorenzoPredictor<T, N, 2>>(conf.absErrorBound));
        }
    }
    if (conf.regression) {
        if (use_single_predictor) {
            return QoZ::make_sz_general_compressor<T, N>(
                    QoZ::make_sz_general_frontend<T, N>(conf, QoZ::RegressionPredictor<T, N>(conf.blockSize, conf.absErrorBound),
                                                       quantizer), encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<QoZ::RegressionPredictor<T, N>>(conf.blockSize, conf.absErrorBound));
        }
    }

    if (conf.regression2) {
        if (use_single_predictor) {
            return QoZ::make_sz_general_compressor<T, N>(
                    QoZ::make_sz_general_frontend<T, N>(conf, QoZ::PolyRegressionPredictor<T, N>(conf.blockSize, conf.absErrorBound),
                                                       quantizer), encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<QoZ::PolyRegressionPredictor<T, N>>(conf.blockSize, conf.absErrorBound));
        }
    }
    return QoZ::make_sz_general_compressor<T, N>(
            QoZ::make_sz_general_frontend<T, N>(conf, QoZ::ComposedPredictor<T, N>(predictors),
                                               quantizer), encoder, lossless);
}


template<class T, QoZ::uint N, class Quantizer, class Quantizer_EB>
QoZ::concepts::CompressorInterface<T>*
make_qoi_lorenzo_compressor(const QoZ::Config &conf, std::shared_ptr<QoZ::concepts::QoIInterface<T, N>> qoi, Quantizer quantizer, Quantizer_EB quantizer_eb) {

    quantizer.clear();
    quantizer_eb.clear();
    //std::shared_ptr<QoZ::concepts::CompressorInterface<T>> sz;

    int methodCnt = (conf.lorenzo + conf.lorenzo2);
    int use_single_predictor = (methodCnt == 1);

    if(use_single_predictor){
        if(conf.lorenzo){
            return QoZ::make_sz_general_compressor<T, N>(QoZ::make_sz_qoi_frontend<T, N>(conf, QoZ::LorenzoPredictor<T, N, 1>(conf.absErrorBound), quantizer, quantizer_eb, qoi),
                                                    QoZ::QoIEncoder<int>(), QoZ::Lossless_zstd());
        }
        else if(conf.lorenzo2){
            return QoZ::make_sz_general_compressor<T, N>(QoZ::make_sz_qoi_frontend<T, N>(conf, QoZ::LorenzoPredictor<T, N, 2>(conf.absErrorBound), quantizer, quantizer_eb, qoi),
                                                    QoZ::QoIEncoder<int>(), QoZ::Lossless_zstd());
        }
    }
    //else{
        std::vector<std::shared_ptr<QoZ::concepts::PredictorInterface<T, N>>> predictors;
        predictors.push_back(std::make_shared<QoZ::LorenzoPredictor<T, N, 1>>(conf.absErrorBound));
        predictors.push_back(std::make_shared<QoZ::LorenzoPredictor<T, N, 2>>(conf.absErrorBound));
        return QoZ::make_sz_general_compressor<T, N>(QoZ::make_sz_qoi_frontend<T, N>(conf, QoZ::ComposedPredictor<T, N>(predictors), quantizer, quantizer_eb, qoi),
                                                QoZ::QoIEncoder<int>(), QoZ::Lossless_zstd());
    //}
    //return sz;
    //return NULL;
}

template<class T, QoZ::uint N>
char *SZ_compress_LorenzoReg(QoZ::Config &conf, T *data, size_t &outSize) {

    assert(N == conf.N);
    assert(conf.cmprAlgo == QoZ::ALGO_LORENZO_REG);
    //QoZ::calAbsErrorBound(conf, data);

    char *cmpData;

    if(conf.qoi > 0){
        //std::cout << "absErrorBound = " << conf.absErrorBound << std::endl;
        //std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << " " << conf.qoiRegionSize << std::endl;
        auto quantizer = QoZ::VariableEBLinearQuantizer<T, T>(conf.quantbinCnt / 2);
        auto quantizer_eb = QoZ::EBLogQuantizer<T>(conf.qoiEBBase, conf.qoiEBLogBase, conf.qoiQuantbinCnt / 2);
        auto qoi = QoZ::GetQOI<T, N>(conf);
        double max_abs_eb = 0;
        if(conf.qoi == 3){
            conf.blockSize = conf.qoiRegionSize;
        }
        // use sampling to determine abs bound
        {
            auto dims = conf.dims;
            auto tmp_abs_eb = conf.absErrorBound;
            //double max_abs_eb = 0;

            size_t sampling_num, sampling_block;
            std::vector<size_t> sample_dims(N);
            std::vector<T> samples = QoZ::sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);
            conf.setDims(sample_dims.begin(), sample_dims.end());

            auto sz = make_qoi_lorenzo_compressor(conf, qoi, quantizer, quantizer_eb);
            T * sampling_data = (T *) malloc(sampling_num * sizeof(T));
            // get current ratio
            double ratio = 0;
            {
                size_t sampleOutSize;
                memcpy(sampling_data, samples.data(), sampling_num * sizeof(T));
                auto cmprData = sz->compress(conf, sampling_data, sampleOutSize);
                //max_abs_eb = sz.get_max_eb();
                delete[]cmprData;
                ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;                
                std::cout << "current_eb = " << conf.absErrorBound << ", current_ratio = " << ratio << std::endl;
                max_abs_eb = fabs(sampling_data[0] - samples[0]);
                for(size_t i=1; i<sampling_num; i++){
                    max_abs_eb = std::max(max_abs_eb, fabs(sampling_data[i] - samples[i]));
                }
            }
            double prev_ratio = 1;
            double current_ratio = ratio;
            double best_abs_eb = std::min(conf.absErrorBound, max_abs_eb);
            double best_ratio = current_ratio;
            // check smaller bounds
            int max_iter = 100; 
            int iter = 0;
            while(iter++ < max_iter){
                auto prev_eb = conf.absErrorBound;
                prev_ratio = current_ratio;
                conf.absErrorBound /= 2;
                qoi->set_global_eb(conf.absErrorBound);
                size_t sampleOutSize;
                memcpy(sampling_data, samples.data(), sampling_num * sizeof(T));
                // reset variables for average of square
                if(conf.qoi == 3) qoi->init();
                auto cmprData = sz->compress(conf, sampling_data, sampleOutSize);
                delete[]cmprData;
                current_ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;                
                std::cout << "current_eb = " << conf.absErrorBound << ", current_ratio = " << current_ratio << std::endl;
                if(current_ratio < prev_ratio * 0.99){
                    if(prev_ratio > best_ratio){
                        best_abs_eb = prev_eb;
                        best_ratio = prev_ratio;
                    }
                    break;
                }
            }
            // set error bound
            free(sampling_data);
            //std::cout << "Best abs eb / pre-set eb: " << best_abs_eb / tmp_abs_eb << std::endl; 
            //std::cout << best_abs_eb << " " << tmp_abs_eb << std::endl;
            conf.absErrorBound = best_abs_eb;
            qoi->set_global_eb(best_abs_eb);
            conf.setDims(dims.begin(), dims.end());
        }
        auto sz = make_qoi_lorenzo_compressor(conf, qoi, quantizer, quantizer_eb);
        cmpData = (char *) sz->compress(conf, data, outSize);
        return cmpData;
    }


    auto quantizer = QoZ::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2);

    if (N == 3 and !conf.regression2 ) {
        // use fast version for 3D
        auto sz = QoZ::make_sz_general_compressor<T, N>(QoZ::make_sz_fast_frontend<T, N>(conf, quantizer), QoZ::HuffmanEncoder<int>(),
                                                       QoZ::Lossless_zstd());
        cmpData = (char *) sz->compress(conf, data, outSize);
    } else {
        auto sz = make_lorenzo_regression_compressor<T, N>(conf, quantizer, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
        //std::cout<<"lor1"<<std::endl;
        cmpData = (char *) sz->compress(conf, data, outSize);
    }
    return cmpData;
}


template<class T, QoZ::uint N>
void SZ_decompress_LorenzoReg(const QoZ::Config &theconf, char *cmpData, size_t cmpSize, T *decData) {
    QoZ::Config conf(theconf);
    assert(conf.cmprAlgo == QoZ::ALGO_LORENZO_REG);
    QoZ::uchar const *cmpDataPos = (QoZ::uchar *) cmpData;

    if(conf.qoi > 0){
        //std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << " " << conf.qoiRegionSize << std::endl;
        auto quantizer = QoZ::VariableEBLinearQuantizer<T, T>(conf.quantbinCnt / 2);
        auto quantizer_eb = QoZ::EBLogQuantizer<T>(conf.qoiEBBase, conf.qoiEBLogBase, conf.qoiQuantbinCnt / 2);
        auto qoi = QoZ::GetQOI<T, N>(conf);
        auto sz = make_qoi_lorenzo_compressor(conf, qoi, quantizer, quantizer_eb);
        sz->decompress(cmpDataPos, cmpSize, decData);
        return;
    }  


    QoZ::LinearQuantizer<T> quantizer;
  
        
    if (N == 3 and !conf.regression2) {
        // use fast version for 3D
        auto sz = QoZ::make_sz_general_compressor<T, N>(QoZ::make_sz_fast_frontend<T, N>(conf, quantizer),
                                                       QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
        sz->decompress(cmpDataPos, cmpSize, decData);
       
    } else {
        auto sz = make_lorenzo_regression_compressor<T, N>(conf, quantizer, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
        sz->decompress(cmpDataPos, cmpSize, decData);
       
    }
    
    




}

#endif