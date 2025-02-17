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
#include "QoZ/qoi/QoIInfo.hpp"
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
std::shared_ptr<QoZ::concepts::CompressorInterface<T>>
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
std::shared_ptr<QoZ::concepts::CompressorInterface<T>>
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
char *SZ_compress_LorenzoReg(QoZ::Config &conf, T *data, size_t &outSize, bool tuning=false) {

    assert(N == conf.N);
    assert(conf.cmprAlgo == QoZ::ALGO_LORENZO_REG);
    //QoZ::calAbsErrorBound(conf, data);
    //std::cout<<"ABSEB "<<conf.absErrorBound<<std::endl;
    char *cmpData;

    if(conf.qoi > 0 and !conf.use_global_eb){
        //std::cout << "absErrorBound = " << conf.absErrorBound << std::endl;
        //std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << " " << conf.qoiRegionSize << std::endl;
        auto quantizer = QoZ::VariableEBLinearQuantizer<T, T>(conf.quantbinCnt / 2);
        auto quantizer_eb = QoZ::EBLogQuantizer<T>(conf.qoiEBBase, conf.qoiEBLogBase, conf.qoiQuantbinCnt / 2, conf.absErrorBound);
        //auto qoi = QoZ::GetQOI<T, N>({conf,conf,conf});
        std::shared_ptr<QoZ::concepts::QoIInterface<T, N>> qoi = nullptr;
        if(conf.qoi == 3){
            conf.blockSize = conf.qoiRegionSize;
        }
        //std::cout<<"1"<<std::e
        auto sz = make_qoi_lorenzo_compressor(conf, qoi, quantizer, quantizer_eb);
        cmpData = (char *) sz->compress(conf, data, outSize);
        return cmpData;
    }


    auto quantizer = QoZ::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2);

    if (N == 3 and !conf.regression2) {
        // use fast version for 3D
        auto sz = QoZ::make_sz_general_compressor<T, N>(QoZ::make_sz_fast_frontend<T, N>(conf, quantizer), QoZ::HuffmanEncoder<int>(),
                                                       QoZ::Lossless_zstd());
        cmpData = (char *) sz->compress(conf, data, outSize);
    } else {
        auto sz = make_lorenzo_regression_compressor<T, N>(conf, quantizer, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
        //std::cout<<"lor1"<<std::endl;
        cmpData = (char *) sz->compress(conf, data, outSize);
        
    }
    if(conf.qoi>0 and !tuning)
        conf.qoi = 0;
    return cmpData;
}

template<class T, QoZ::uint N>
std::array<char *,3> SZ_compress_LorenzoReg_Vec(std::array<QoZ::Config,3> &confs, std::array<T *,3>&data, std::array<size_t,3> &outSizes, bool tuning=false) {

    assert(N == confs[0].N and N == confs[1].N and N == confs[2].N);
    assert(confs[0].cmprAlgo == QoZ::ALGO_LORENZO_REG and confs[1].cmprAlgo == QoZ::ALGO_LORENZO_REG and confs[2].cmprAlgo == QoZ::ALGO_LORENZO_REG);
    //QoZ::calAbsErrorBound(conf, data);
    //std::cout<<"ABSEB "<<conf.absErrorBound<<std::endl;
    std::array<char *,3>cmpData;
    std::array<std::vector<T>,3> ori_data;
    int ori_qoi = 0;
    if(confs[0].qoi>0){
        ori_qoi = confs[0].qoi;

        if(confs[0].qoiRegionMode==0){
            for (auto i:{0,1,2})
                ori_data[i]=std::vector(data[i],data[i]+confs[i].num);
        }
        //if(!confs[0].qoi_tuned)
        
        //   QoI_tuning<T,N>(confs, data);
    }

    std::array<bool,3>qoi_used = {false,false,false};
    for (auto i:{0,1,2}){
        //std::cout<<confs[i].qoi<<" "<<confs[i].use_global_eb<<std::endl;
        if (confs[i].qoi>0)
            confs[i].qoi = 99;//empty qoi;
        if ( confs[i].qoi>0 and !confs[i].use_global_eb)
            //std::cout<<"Compress Data "<<i<<" with qoi interpolator"<<std::endl;
            qoi_used[i]=true;

        cmpData[i] = SZ_compress_LorenzoReg<T,N>(confs[i], data[i], outSizes[i]);

        confs[i].ebs.clear();
        confs[i].ebs.shrink_to_fit();
    }


    if(ori_qoi>0 and confs[0].qoiRegionMode==0){
        int conf_ori_qoi = confs[0].qoi;
        confs[0].qoi = ori_qoi;
        auto qoi = QoZ::GetQOI<T, N>(confs);//todo: avoid duplicated initialization.
        confs[0].qoi = conf_ori_qoi;
        
        
        for(size_t i=0;i<confs[0].num;i++){

            if (qoi->check_compliance(ori_data[0][i],ori_data[1][i],ori_data[2][i],data[0][i],data[1][i],data[2][i])){
                for (auto j:{0,1,2})
                    ori_data[j][i]=0;
            }
            else{
                for (auto j:{0,1,2}){
                    T offset = ori_data[j][i]-data[j][i];
                    data[j][i] = ori_data[j][i];
                    ori_data[j][i]=offset;
                }
            }
        }

        auto zstd = QoZ::Lossless_zstd();
        
        for (auto i:{0,1,2}){
            size_t offset_size;
            QoZ::uchar *lossless_data = zstd.compress(reinterpret_cast< QoZ::uchar *>(ori_data[i].data()),
                                                         confs[i].num*sizeof(T),
                                                         offset_size);
            ori_data[i].clear();
            size_t newSize = outSizes[i] + offset_size + QoZ::Config::size_est()  + 100;
            char * newcmpData = new char[newSize];
            memcpy(newcmpData,cmpData[i],outSizes[i]);
            delete [] cmpData[i];
            memcpy(newcmpData+outSizes[i],lossless_data,offset_size);
            
            outSizes[i]+=offset_size;
            //std::cout<<offset_size<<" "<<outSizes[i]<<std::endl;
            delete []lossless_data;
            //lossless_data = NULL;
            memcpy(newcmpData+outSizes[i],&offset_size,sizeof(size_t));
            //
            outSizes[i]+=sizeof(size_t);

            cmpData[i] = newcmpData;





          //  std::cout<<offset_size<<" "<<outSizes[i]<<std::endl;



            
        }
    }
    else{
        size_t offset_size=0;
        for(auto i:{0,1,2}){
            memcpy(cmpData[i]+outSizes[i],&offset_size,sizeof(size_t));
            outSizes[i]+=sizeof(size_t);
        }

    }

    for (auto i:{0,1,2}){
        if (!qoi_used[i])
            confs[i].qoi = 0;
    }

    /*


    if(conf.qoi > 0 and !conf.use_global_eb){
        //std::cout << "absErrorBound = " << conf.absErrorBound << std::endl;
        //std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << " " << conf.qoiRegionSize << std::endl;
        auto quantizer = QoZ::VariableEBLinearQuantizer<T, T>(conf.quantbinCnt / 2);
        auto quantizer_eb = QoZ::EBLogQuantizer<T>(conf.qoiEBBase, conf.qoiEBLogBase, conf.qoiQuantbinCnt / 2, conf.absErrorBound);
        auto qoi = QoZ::GetQOI<T, N>(conf);
        if(conf.qoi == 3){
            conf.blockSize = conf.qoiRegionSize;
        }
        //std::cout<<"1"<<std::e
        auto sz = make_qoi_lorenzo_compressor(conf, qoi, quantizer, quantizer_eb);
        cmpData = (char *) sz->compress(conf, data, outSize);
        return cmpData;
    }

    confs[i].ebs.clear();
    confs[i].ebs.shrink_to_fit();


    auto quantizer = QoZ::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2);

    if (N == 3 and !conf.regression2 and conf.qoi==0) {
        // use fast version for 3D
        auto sz = QoZ::make_sz_general_compressor<T, N>(QoZ::make_sz_fast_frontend<T, N>(conf, quantizer), QoZ::HuffmanEncoder<int>(),
                                                       QoZ::Lossless_zstd());
        cmpData = (char *) sz->compress(conf, data, outSize);
    } else {
        auto sz = make_lorenzo_regression_compressor<T, N>(conf, quantizer, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
        //std::cout<<"lor1"<<std::endl;
        cmpData = (char *) sz->compress(conf, data, outSize);
        if(conf.qoi>0 and !tuning)
            conf.qoi = 99;
    }*/
    
    return cmpData;
}



template<class T, QoZ::uint N>
void SZ_decompress_LorenzoReg(const QoZ::Config &theconf, char *cmpData, size_t cmpSize, T *decData) {

    assert(theconf.cmprAlgo == QoZ::ALGO_LORENZO_REG);

    QoZ::Config conf(theconf);
   // std::cout<<"ABSEB "<<conf.absErrorBound<<std::endl;
    assert(conf.cmprAlgo == QoZ::ALGO_LORENZO_REG);
    QoZ::uchar const *cmpDataPos = (QoZ::uchar *) cmpData;

    if(conf.qoi > 0){
        //std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << " " << conf.qoiRegionSize << std::endl;
        auto quantizer = QoZ::VariableEBLinearQuantizer<T, T>(conf.quantbinCnt / 2);
        auto quantizer_eb = QoZ::EBLogQuantizer<T>(conf.qoiEBBase, conf.qoiEBLogBase, conf.qoiQuantbinCnt / 2, conf.absErrorBound);
        //auto qoi = QoZ::GetQOI<T, N>(conf);
        std::shared_ptr<QoZ::concepts::QoIInterface<T, N>> qoi = nullptr;
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

template<class T, QoZ::uint N>
void SZ_decompress_LorenzoReg_Vec(const std::array<QoZ::Config ,3>&confs, std::array<char *,3> &cmpData, std::array<size_t,3> &cmpSizes, std::array<T *,3>&decData) {
    assert(confs[0].cmprAlgo == QoZ::ALGO_LORENZO_REG and confs[1].cmprAlgo == QoZ::ALGO_LORENZO_REG and confs[2].cmprAlgo == QoZ::ALGO_LORENZO_REG);

    for (auto i:{0,1,2}){

       

        size_t offset_size=0;
        T* offset_data;
        size_t cmpSize = cmpSizes[i];
        //QoZ::read<size_t>(offset_size,reinterpret_cast< QoZ::uchar const *>(cmpData[i])+cmpSize-sizeof(size_t));
        memcpy(&offset_size,cmpData[i]+cmpSize-sizeof(size_t),sizeof(size_t));
        cmpSize-=sizeof(size_t);     
        if (offset_size!=0){
            cmpSize-=offset_size;
            //outlier_data.resize(confs[i].num);
            auto zstd = QoZ::Lossless_zstd();
            offset_data = reinterpret_cast<T *> ( zstd.decompress(reinterpret_cast<QoZ::uchar *>(cmpData[i])+cmpSize, offset_size) );
            

        }   
        SZ_decompress_LorenzoReg<T,N>(confs[i], cmpData[i], cmpSize, decData[i]);
        if (offset_size!=0){
            for(size_t j=0;j<confs[i].num;j++)
                decData[i][j]+=offset_data[j];
            delete []offset_data;
        }
    }


}

#endif