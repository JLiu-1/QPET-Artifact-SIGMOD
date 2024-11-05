#ifndef SZ3_SZINTERP_HPP
#define SZ3_SZINTERP_HPP

#include "QoZ/compressor/SZInterpolationCompressor.hpp"
#include "QoZ/compressor/SZQoIInterpolationCompressor.hpp"
#include "QoZ/compressor/deprecated/SZBlockInterpolationCompressor.hpp"


#include "QoZ/encoder/QoIEncoder.hpp"
#include "QoZ/quantizer/IntegerQuantizer.hpp"
#include "QoZ/quantizer/QoIIntegerQuantizer.hpp"
#include "QoZ/lossless/Lossless_zstd.hpp"
#include "QoZ/utils/Iterator.hpp"
#include "QoZ/utils/Sample.hpp"
#include "QoZ/utils/Statistic.hpp"
#include "QoZ/utils/Extraction.hpp"
#include "QoZ/utils/QuantOptimization.hpp"
#include "QoZ/utils/Config.hpp"
#include "QoZ/utils/Metrics.hpp"
//#include "QoZ/api/impl/SZLorenzoReg.hpp"


#include "QoZ/qoi/QoIInfo.hpp"

//#include <cunistd>
#include <cmath>
#include <algorithm>
#include <memory>
#include <limits>
#include <cstring>
#include <cstdlib>
#include <queue>

double estimate_rate_Hoeffdin(size_t n, size_t N, double q, double k = 2.0){//n: element_per_block N: num_blocks q: confidence
    //no var information
    //if gaussian, just multiply k 
   
    /*
    if (q>=0.95 and N >= 1000){
        return sqrt( -n / ( 2 * ( log(1-q) - log(2*N) ) ) );
    }
    else{
        return sqrt( -n / ( 2 * ( log(1- pow(q,1.0/N) ) - log(2) ) ) );
    }
    */

    double p;
    if (q>=0.95 and N >= 1000){
        p = (1-q)/N;
    }
    else{
        p = 1- pow(q,1.0/N);
    }


    return k*sqrt(0.5*n/log(2.0/p));// / k*sqrt(0.5/sumai^2*log(2.0/p)) 
    
}

double estimate_rate_Bernstein(size_t n, size_t N, double q, double k = 3.0){//n: element_per_block N: num_blocks q: confidence
    //estimate var as gaussian
    double p;
    if (q>=0.95 and N >= 1000){
        p = (1-q)/N;
    }
    else{
        p = 1- pow(q,1.0/N);
    }


    return (k*k/6)*(sqrt(1+18*n/(k*k*log(2/p)))-1);
    
}

/*
double estimate_rate_Gaussian(size_t n, size_t N, double q, double k = 3.0){//n: element_per_block N: num_blocks q: confidence
    //assume gaussian error
    //inaccurate
    double p;
    if (q>=0.95 and N >= 1000){
        p = (1-q)/N;
    }
    else{
        p = 1- pow(q,1.0/N);
    }


    return k*sqrt(0.5*n/log(2.0/p));
    
}
*/

template<class T, QoZ::uint N>
void QoI_tuning(std::array<QoZ::Config,3> &conf, std::array<T *,3> &data);




template<class T, QoZ::uint N>
std::array<char *,3>SZ_compress_Interp(std::array<QoZ::Config,3> &confs, std::array<T *,3> data, std::array<size_t,3> &outSizes) {

//    std::cout << "****************** Interp Compression ****************" << std::endl;
//    std::cout << "Interp Op          = " << interpAlgo << std::endl
//              << "Direction          = " << direction << std::endl
//              << "SZ block size      = " << blockSize << std::endl
//              << "Interp block size  = " << interpBlockSize << std::endl;

    assert(N == confs[0].N&&N == confs[1].N&&N == confs[2].N);
    //assert(confs[0].cmprAlgo == QoZ::ALGO_INTERP);
    for (auto i:{0,1,2})
        QoZ::calAbsErrorBound(confs[i], data[i]);
    std::array<char*,3> cmpData;
    std::array<std::vector<T>,3> ori_data;
    size_t offset_size=0;


    int ori_qoi = 0;
    if(confs[0].qoi>0){
        ori_qoi = confs[0].qoi;

        if(confs[0].qoiRegionMode==0){
            for (auto i:{0,1,2})
                ori_data[i]=std::vector(data[i],data[i]+confs[i].num);
        }

        if(!confs[0].qoi_tuned)
        
           QoI_tuning<T,N>(confs, data);
    
    }

    //conf.print();
    std::array<bool,3>qoi_used = {false,false,false};
    for (auto i:{0,1,2}){
        //std::cout<<confs[i].qoi<<" "<<confs[i].use_global_eb<<std::endl;
        if (confs[i].qoi>0)
            confs[i].qoi = 10;//empty qoi;
        if ( confs[i].qoi>0 and !confs[i].use_global_eb){
            //std::cout<<"Compress Data "<<i<<" with qoi interpolator"<<std::endl;
            qoi_used[i]=true;
            
            auto qoi = QoZ::GetQOI<T, N>(confs);//todo: bring qoi to conf to avoid duplicated initialization.
            
            
            auto quantizer= QoZ::VariableEBLinearQuantizer<T, T>(confs[i].quantbinCnt / 2);

            //auto quantizer = QoZ::VariableEBLinearQuantizer<T, T>(confs[0].quantbinCnt / 2);

            auto quantizer_eb = QoZ::EBLogQuantizer<T>(confs[i].qoiEBBase, confs[i].qoiEBLogBase, confs[i].qoiQuantbinCnt / 2, confs[i].absErrorBound);
                                                                  
            //auto quantizer_eb = QoZ::EBLogQuantizer<T>(confs[0].qoiEBBase, conf.qoiEBLogBase, conf.qoiQuantbinCnt / 2, conf.absErrorBound);
            auto sz = QoZ::SZQoIInterpolationCompressor<T, N, QoZ::VariableEBLinearQuantizer<T, T>, QoZ::EBLogQuantizer<T>, QoZ::QoIEncoder<int>, QoZ::Lossless_zstd>(
                    quantizer, quantizer_eb, qoi, QoZ::QoIEncoder<int>(), QoZ::Lossless_zstd());//maybe we need multiple encoder and zstd

            cmpData[i] = (char *) sz.compress(confs[i], data[i], outSizes[i]);
            
             //double incall_time = timer.stop();
            //std::cout << "incall time = " << incall_time << "s" << std::endl;

        }
        else{
           
                auto sz = QoZ::SZInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
                        QoZ::LinearQuantizer<T>(confs[i].absErrorBound, confs[i].quantbinCnt / 2),
                        QoZ::HuffmanEncoder<int>(),
                        QoZ::Lossless_zstd());

               
                //QoZ::Timer timer;

                //timer.start();
                cmpData[i] = (char *) sz.compress(confs[i], data[i], outSizes[i]);
                //confs[i].qoi = 0;
                 //double incall_time = timer.stop();
                //std::cout << "incall time = " << incall_time << "s" << std::endl;
            
            
        }
        
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
                for (auto j:{0,1,2})
                    ori_data[j][i]-=data[j][i];
            }
        }

        auto zstd = QoZ::Lossless_zstd();
        
        for (auto i:{0,1,2}){
            QoZ::uchar *lossless_data = zstd.compress(reinterpret_cast< QoZ::uchar *>(ori_data[i].data()),
                                                         confs[i].num*sizeof(T),
                                                         offset_size);
            ori_data[i].clear();
            memcpy(cmpData[i]+outSizes[i],lossless_data,offset_size);
            outSizes[i]+=offset_size;
            delete []lossless_data;

            memcpy(cmpData[i]+outSizes[i],&offset_size,sizeof(size_t));
            outSizes[i]+=sizeof(size_t);

            
        }
    }
    else{
        offset_size=0;
        for(auto i:{0,1,2}){
            memcpy(cmpData[i]+outSizes[i],&offset_size,sizeof(size_t));
            outSizes[i]+=sizeof(size_t);
        }

    }

    for (auto i:{0,1,2}){
        if (!qoi_used[i])
            confs[i].qoi = 0;
    }


    return cmpData;
}

template<class T, QoZ::uint N>
void SZ_decompress_Interp(std::array<QoZ::Config ,3>&confs, std::array<char *,3> &cmpData, std::array<size_t,3> &cmpSizes, std::array<T *,3>&decData) {
    assert(confs[0].cmprAlgo == QoZ::ALGO_INTERP&&confs[1].cmprAlgo == QoZ::ALGO_INTERP&&confs[2].cmprAlgo == QoZ::ALGO_INTERP);



    
   //std::cout<<"5 "<<confs[0].qoi<<std::endl;
   for (auto i:{0,1,2}){

        size_t offset_size=0;
        T* offset_data;
        size_t cmpSize = cmpSizes[i];
        //QoZ::read<size_t>(offset_size,reinterpret_cast< QoZ::uchar const *>(cmpData[i])+cmpSize-sizeof(size_t));
        memcpy(&offset_size,cmpData[i]+cmpSize-sizeof(size_t),sizeof(size_t));
        cmpSize-=sizeof(size_t);     
        if (offset_size!=0){
            //outlier_data.resize(confs[i].num);
            auto zstd = QoZ::Lossless_zstd();
            offset_data = reinterpret_cast<T *> ( zstd.decompress(reinterpret_cast<QoZ::uchar *>(cmpData[i])+cmpSize-offset_size, offset_size) );
            cmpSize-=offset_size;

        }   
       if(confs[i].qoi > 0){
            //QoZ::uchar const *cmpDataPos = (QoZ::uchar *) cmpData;
            //std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << std::endl;
            //std::array<QoZ::uchar *,3>cmpDataPos = {(QoZ::uchar *) cmpData[0],(QoZ::uchar *) cmpData[1],(QoZ::uchar *) cmpData[2]};
            confs[i].qoi = 10; //empty qoi
            auto qoi = QoZ::GetQOI<T, N>(confs);
           

                QoZ::uchar const *cmpDataPos = (QoZ::uchar *) cmpData[i];
                auto quantizer=QoZ::VariableEBLinearQuantizer<T, T>(confs[i].quantbinCnt / 2);

                //auto quantizer = QoZ::VariableEBLinearQuantizer<T, T>(confs[0].quantbinCnt / 2);

                auto quantizer_eb =QoZ::EBLogQuantizer<T>(confs[i].qoiEBBase, confs[i].qoiEBLogBase, confs[i].qoiQuantbinCnt / 2, confs[i].absErrorBound);
                
                auto sz = QoZ::SZQoIInterpolationCompressor<T, N, QoZ::VariableEBLinearQuantizer<T, T>, QoZ::EBLogQuantizer<T>, QoZ::QoIEncoder<int>, QoZ::Lossless_zstd>(
                        quantizer, quantizer_eb, qoi, QoZ::QoIEncoder<int>(), QoZ::Lossless_zstd());//may need 3 encoders
                sz.decompress(cmpDataPos, cmpSize, decData[i]);
            
            
        }
        else{   
        
            QoZ::uchar const *cmpDataPos = (QoZ::uchar *) cmpData[i];
            auto sz = QoZ::SZInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
                    QoZ::LinearQuantizer<T>(),
                    QoZ::HuffmanEncoder<int>(),
                    QoZ::Lossless_zstd());
            sz.decompress(cmpDataPos, cmpSize, decData[i]);
        }
        if (offset_size!=0){
            for(size_t j=0;j<confs[i].num;j++)
                decData[i][j]+=offset_data[j];
            delete []offset_data;
        }
    }
        
}


template<class T, QoZ::uint N>
double do_not_use_this_interp_compress_block_test(T *data, std::vector<size_t> dims, size_t num,
                                                  double eb, int interp_op, int direction_op, int block_size) {
    std::vector<T> data1(data, data + num);
    size_t outSize = 0;
    QoZ::Config conf;
    conf.absErrorBound = eb;
    conf.setDims(dims.begin(), dims.end());
    conf.blockSize = block_size;
    conf.interpMeta.interpAlgo = interp_op;
    conf.interpMeta.interpDirection = direction_op;
    auto sz = QoZ::SZBlockInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
            QoZ::LinearQuantizer<T>(eb),
            QoZ::HuffmanEncoder<int>(),
            QoZ::Lossless_zstd());
    char *cmpData = (char *) sz.compress(conf, data1.data(), outSize);
    delete[]cmpData;
    auto compression_ratio = num * sizeof(T) * 1.0 / outSize;
    return compression_ratio;
}



template<class T, QoZ::uint N>
inline void init_alphalist(std::vector<double> &alpha_list,const double &rel_bound, QoZ::Config &conf){

    
    /*
    if(conf.linearReduce){
        alpha_list={0,0.1,0.2,0.3,0.4,0.5};

    }
    */
    //else{
        if (conf.tuningTarget!=QoZ::TUNING_TARGET_CR){
           
            if(conf.abList==0)            
                alpha_list={1,1.25,1.5,1.75,2};
            else if(conf.abList==1)
                alpha_list={1,1.25,1.5,1.75,2,2.25,2.5};
            else
                alpha_list={1,1.25,1.5,1.75,2,2.25,2.5,2.75,3};
           
        }
        else{
            
            alpha_list={-1,1,1.25,1.5,1.75,2};
           
        }
    //}
}
template<class T, QoZ::uint N>
inline void init_betalist(std::vector<double> &beta_list,const double &rel_bound, QoZ::Config &conf){
   
    /*
    if(conf.linearReduce){
        beta_list={1,0.75,0.5,0.33,0.25};
    }
    */
    //else{
        if (conf.tuningTarget!=QoZ::TUNING_TARGET_CR){    
            
            beta_list={1.5,2,3,4};//may remove 1.5
           
        }
        else {
          
            beta_list={-1,1.5,2,3};
           
        }
    //}
}




template<class T, QoZ::uint N>
std::pair<double,double> CompressTest(const QoZ::Config &conf,std::vector< std::vector<T> > & sampled_blocks,QoZ::ALGO algo = QoZ::ALGO_INTERP,
                    QoZ::TUNING_TARGET tuningTarget=QoZ::TUNING_TARGET_RD,bool useFast=true,double profiling_coeff=1,const std::vector<double> &orig_means=std::vector<double>(),
                    const std::vector<double> &orig_sigma2s=std::vector<double>(),const std::vector<double> &orig_ranges=std::vector<double>(),const std::vector<T> &flattened_sampled_data=std::vector<T>()){
    QoZ::Config testConfig(conf);
    size_t ssim_size=conf.SSIMBlockSize;
    std::vector< std::vector<T> > ori_sampled_blocks;
    if(tuningTarget!=QoZ::TUNING_TARGET_CR and (algo != QoZ::ALGO_INTERP or tuningTarget!=QoZ::TUNING_TARGET_RD)) 
        ori_sampled_blocks=sampled_blocks;
    /*    
    if(algo == QoZ::ALGO_LORENZO_REG){
        testConfig.cmprAlgo = QoZ::ALGO_LORENZO_REG;
        testConfig.dims=conf.dims;
        testConfig.num=conf.num;
        testConfig.lorenzo = true;
        testConfig.lorenzo2 = true;
        testConfig.regression = false;
        testConfig.regression2 = false;
        testConfig.openmp = false;
        testConfig.blockSize = 5;//why?
        testConfig.quantbinCnt = 65536 * 2;
    }*/
    double square_error=0.0;
    double bitrate=0.0;
    double metric=0.0;
    size_t sampleBlockSize=testConfig.sampleBlockSize;
    size_t num_sampled_blocks=sampled_blocks.size();
    size_t per_block_ele_num=pow(sampleBlockSize+1,N);
    size_t ele_num=num_sampled_blocks*per_block_ele_num;
    std::vector<T> cur_block(testConfig.num,0);
    std::vector<int> q_bins;
    std::vector<std::vector<int> > block_q_bins;
    std::vector<size_t> q_bin_counts;
    std::vector<T> flattened_cur_blocks;
    size_t idx=0;   
    QoZ::concepts::CompressorInterface<T> *sz;
    size_t totalOutSize=0;
    /*
    if(algo == QoZ::ALGO_LORENZO_REG){
        auto quantizer = QoZ::LinearQuantizer<T>(testConfig.absErrorBound, testConfig.quantbinCnt / 2);
        if (useFast &&N == 3 && !testConfig.regression2) {
            sz = QoZ::make_sz_general_compressor<T, N>(QoZ::make_sz_fast_frontend<T, N>(testConfig, quantizer), QoZ::HuffmanEncoder<int>(),
                                                                   QoZ::Lossless_zstd());
        }
        else{
            sz = make_lorenzo_regression_compressor<T, N>(testConfig, quantizer, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());

        }
    }*/
   // else 
    if(algo == QoZ::ALGO_INTERP){

        sz =  new QoZ::SZInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
                        QoZ::LinearQuantizer<T>(testConfig.absErrorBound),
                        QoZ::HuffmanEncoder<int>(),
                        QoZ::Lossless_zstd());

    }
    else{
        std::cout<<"algo type error!"<<std::endl;
        return std::pair<double,double>(0,0);
    }
                           
    for (int k=0;k<num_sampled_blocks;k++){
        size_t sampleOutSize;
        //std::vector<T> cur_block(testConfig.num);
        //std::copy(sampled_blocks[k].begin(),sampled_blocks[k].end(),cur_block.begin());
        
        char *cmprData;
         
        
        cmprData = (char*)sz->compress(testConfig, sampled_blocks[k].data(), sampleOutSize,1);

        delete[]cmprData;
        
        

        
        if(algo==QoZ::ALGO_INTERP){
            block_q_bins.push_back(testConfig.quant_bins);
        }

        if(tuningTarget==QoZ::TUNING_TARGET_RD){
            if(algo==QoZ::ALGO_INTERP)
                square_error+=testConfig.decomp_square_error;
            else{
               
                for(size_t j=0;j<per_block_ele_num;j++){
                    T value=sampled_blocks[k][j]-ori_sampled_blocks[k][j];
                    square_error+=value*value;
                
                }
            }
        }
        else if (tuningTarget==QoZ::TUNING_TARGET_SSIM){
            size_t ssim_block_num=orig_means.size();                       
            double mean=0,sigma2=0,cov=0,range=0;
            double orig_mean=0,orig_sigma2=0,orig_range=0;  
            std::vector<size_t>block_dims(N,sampleBlockSize+1);                      
            if(N==2){
                for (size_t i=0;i+ssim_size<sampleBlockSize+1;i+=ssim_size){
                    for (size_t j=0;j+ssim_size<sampleBlockSize+1;j+=ssim_size){
                        orig_mean=orig_means[idx];
                        orig_sigma2=orig_sigma2s[idx];
                        orig_range=orig_ranges[idx];
                        std::vector<size_t> starts{i,j};
                        QoZ::blockwise_profiling<T>(sampled_blocks[k].data(),block_dims,starts,ssim_size,mean,sigma2,range);
                        cov=QoZ::blockwise_cov<T>(ori_sampled_blocks[k].data(),sampled_blocks[k].data(),block_dims,starts,ssim_size,orig_mean,mean);
                        metric+=QoZ::SSIM(orig_range,orig_mean,orig_sigma2,mean,sigma2,cov)/ssim_block_num;
                        idx++;


                    }
                }
            }
            else if(N==3){
                for (size_t i=0;i+ssim_size<sampleBlockSize+1;i+=ssim_size){
                    for (size_t j=0;j+ssim_size<sampleBlockSize+1;j+=ssim_size){
                        for (size_t kk=0;kk+ssim_size<sampleBlockSize+1;kk+=ssim_size){
                            orig_mean=orig_means[idx];
                            orig_sigma2=orig_sigma2s[idx];
                            orig_range=orig_ranges[idx];
                            std::vector<size_t> starts{i,j,kk};
                            QoZ::blockwise_profiling<T>(sampled_blocks[k].data(),block_dims,starts,ssim_size,mean,sigma2,range);
                            cov=QoZ::blockwise_cov<T>(ori_sampled_blocks[k].data(),sampled_blocks[k].data(),block_dims,starts,ssim_size,orig_mean,mean);
                            //printf("%.8f %.8f %.8f %.8f %.8f %.8f %.8f\n",orig_range,orig_sigma2,orig_mean,range,sigma2,mean,cov);
                            metric+=QoZ::SSIM(orig_range,orig_mean,orig_sigma2,mean,sigma2,cov)/ssim_block_num;
                     
                            idx++;
                        }
                    }
                }
            }
        }
        else if (tuningTarget==QoZ::TUNING_TARGET_AC){
            flattened_cur_blocks.insert(flattened_cur_blocks.end(),sampled_blocks[k].begin(),sampled_blocks[k].end());
        }                      
    }
    if(algo==QoZ::ALGO_INTERP ){
        q_bin_counts=testConfig.quant_bin_counts;
        size_t level_num=q_bin_counts.size();
        size_t last_pos=0;
        for(int k=level_num-1;k>=0;k--){
            for (size_t l =0;l<num_sampled_blocks;l++){
                for (size_t m=last_pos;m<q_bin_counts[k];m++){
                    q_bins.push_back(block_q_bins[l][m]);
                }
            }
            last_pos=q_bin_counts[k];
        }      
    }
    size_t sampleOutSize;
    
    auto cmprData=sz->encoding_lossless(totalOutSize,q_bins);             
    delete[]cmprData;
  
    
    bitrate=8*double(totalOutSize)/ele_num;
    
    bitrate*=profiling_coeff;
    if(tuningTarget==QoZ::TUNING_TARGET_RD){
        double mse=square_error/ele_num;
        mse*=profiling_coeff;      
        
        metric=QoZ::PSNR(testConfig.rng,mse);
    }
    else if (tuningTarget==QoZ::TUNING_TARGET_AC){                       
        metric=1.0-QoZ::autocorrelation<T>(flattened_sampled_data.data(),flattened_cur_blocks.data(),ele_num);                        
    }                    
   

    if(algo==QoZ::ALGO_LORENZO_REG)    {
        bitrate*=testConfig.lorenzoBrFix;
    }
    delete sz;
    return std::pair(bitrate,metric);
}

template<class T, QoZ::uint N>
double CompressTest_QoI(const QoZ::Config &conf,const std::vector< std::vector<T> > & sampled_blocks, const std::vector<std::vector<double> > & eb_blocks,std::shared_ptr<QoZ::concepts::QoIInterface<T, N>> qoi=nullptr){
    QoZ::Config testConfig(conf);

    if(qoi==nullptr){
        qoi = QoZ::GetQOI<T, N>({testConfig,testConfig,testConfig});
    }
    
    double bitrate=0.0;
    size_t sampleBlockSize=testConfig.sampleBlockSize;
    size_t num_sampled_blocks=sampled_blocks.size();
    size_t per_block_ele_num=pow(sampleBlockSize+1,N);
    size_t ele_num=num_sampled_blocks*per_block_ele_num;

    std::vector<T> cur_block(testConfig.num,0);
    std::vector<int> q_bins;
    std::vector<int> q_bins_eb;
    std::vector<std::vector<int> > block_q_bins;
    std::vector<std::vector<int> > block_q_bins_eb;
    std::vector<size_t> q_bin_counts;
    std::vector<T> flattened_cur_blocks;
    size_t idx=0;   
    size_t totalOutSize=0;

    testConfig.qoiEBBase = testConfig.absErrorBound / 1030;
    auto quantizer = QoZ::VariableEBLinearQuantizer<T, T>(testConfig.quantbinCnt / 2);
    auto quantizer_eb = QoZ::EBLogQuantizer<T>(testConfig.qoiEBBase, testConfig.qoiEBLogBase, testConfig.qoiQuantbinCnt / 2, testConfig.absErrorBound);
            
    auto sz = QoZ::SZQoIInterpolationCompressor<T, N, QoZ::VariableEBLinearQuantizer<T, T>, QoZ::EBLogQuantizer<T>, QoZ::QoIEncoder<int>, QoZ::Lossless_zstd>(
                    quantizer, quantizer_eb, qoi, QoZ::QoIEncoder<int>(), QoZ::Lossless_zstd());
                             
    for (int k=0;k<num_sampled_blocks;k++){
        size_t sampleOutSize;
        std::vector<T> cur_block(testConfig.num);
        std::copy(sampled_blocks[k].begin(),sampled_blocks[k].end(),cur_block.begin());
        
        char *cmprData;

        testConfig.ebs=eb_blocks[k];
         
        
        cmprData = (char*)sz.compress(testConfig, cur_block.data(), sampleOutSize,1);

        delete[]cmprData;
        
        

        
        
        block_q_bins.push_back(testConfig.quant_bins);
        block_q_bins_eb.push_back(testConfig.quant_bins_eb);
        


        
    }
    
    q_bin_counts=testConfig.quant_bin_counts;
    size_t level_num=q_bin_counts.size();
    size_t last_pos=0;
    for(int k=level_num-1;k>=0;k--){
        //std::cout<<q_bin_counts[k]<<std::endl;
        for (size_t l =0;l<num_sampled_blocks;l++){
            for (size_t m=last_pos;m<q_bin_counts[k];m++){
                q_bins.push_back(block_q_bins[l][m]);
                q_bins_eb.push_back(block_q_bins_eb[l][m]);
            }
        }
        last_pos=q_bin_counts[k];
    }      
    //std::cout<<q_bins.size()<<" "<<q_bins_eb.size()<<std::endl;
    
    size_t sampleOutSize;

    q_bins_eb.insert(q_bins_eb.end(),q_bins.begin(),q_bins.end());
   // std::cout<<q_bins_eb.size()<<std::endl;
    
    auto cmprData=sz.encoding_lossless(totalOutSize,q_bins_eb);             
    delete[]cmprData;
  
    
    bitrate=8*double(totalOutSize)/ele_num;
    
    //bitrate*=profiling_coeff;

    return bitrate;
}




template<class T, QoZ::uint N>
void QoI_tuning(std::array<QoZ::Config,3> &confs, std::array<T *,3> &data){
    //std::cout<<confs[0].qoi<<" "<<confs[1].qoi<<" "<<confs[2].qoi<<" "<<std::endl;
    auto qoi_id = confs[0].qoi;
    if(qoi_id){
       
        for (auto i:{0,1,2}){
            if(confs[i].qoiEBBase == 0) 
                confs[i].qoiEBBase = std::numeric_limits<T>::epsilon();
                //confs.qoiEBBase = qoi_rel_eb / 1030
            if(confs[i].qoiEBLogBase == 0)
                confs[i].qoiEBLogBase = 2;   
        }     
        
    }
    else{
        return;
    }
        
        
    auto qoi = QoZ::GetQOI<T, N>(confs);


    if(confs[0].qoiEBMode !=QoZ::EB_ABS){//rel
        
        double max_qoi = -std::numeric_limits<double>::max();
        double min_qoi = std::numeric_limits<double>::max();
       
        for(size_t i=0; i<confs[0].num; i++){
            
            double q = qoi->eval(data[0][i],data[1][i],data[2][i]);
            if(std::isinf(q) or std::isnan(q))
                continue;

            if (max_qoi < q) max_qoi = q;
            if (min_qoi > q) min_qoi = q;
        }

        if (max_qoi == min_qoi){
            max_qoi = 1.0;
            min_qoi = 0.0;
        }
        //std::cout<<max_qoi << " "<<min_qoi<<" "<<conf.qoiEB << std::endl;
        for(auto i:{0,1,2}){
            confs[i].qoiEB *= (max_qoi-min_qoi);
            confs[i].qoiEBMode = QoZ::EB_ABS;
        }
    }
    std::cout<<"ABS QoI eb: " << confs[0].qoiEB << std::endl;
    qoi->set_qoi_tolerance(confs[0].qoiEB);


    for(auto i:{0,1,2})
        confs[i].ebs = std::vector<double>(confs[i].num);

    for (size_t i = 0; i < confs[0].num; i++){
        std::array<T,3>cur_ebs = qoi->interpret_eb(data[0][i],data[1][i],data[2][i]);
        for(auto j:{0,1,2})
           confs[j].ebs[i]=cur_ebs[j];
    }

    //for(auto j:{0,1,2})

    //   std::cout<<j<<" "<<confs[j].ebs.front()<<" "<<confs[j].ebs.back()<<std::endl;

    std::array<double,3> best_abs_ebs;

    double quantile_rate = confs[0].quantile<= 0? confs[0].max_quantile_rate : confs[0].quantile  ;//conf.quantile;//quantile
        //std::cout<<quantile<<std::endl;
    size_t k = std::ceil(quantile_rate * confs[0].num);
    k = std::max((size_t)1, std::min(confs[0].num-1, k)); 

    std::array<std::vector<std::vector<double> >,3>sample_block_ebs;
    std::array<std::vector<std::vector<T> >,3>sampled_blocks;
    std::vector<std::vector<size_t>>starts;


    std::array < std::vector<T>,3 > sampled_regions;
    std::array<std::vector<double>,3>sample_region_ebs;
    size_t sampling_num, sampling_block;
    std::vector<size_t> sample_dims(N);

    QoZ::Config testConf = confs[0];

    std::vector<double>ebs;



    if (confs[0].quantile<=0){
        double best_overall_br = 0.0;


        if(N==2 or N==3){

            if (testConf.sampleBlockSize<=0){
            
                testConf.sampleBlockSize = (N<=2?64:32);
            }

            ebs = std::vector<double>(confs[0].ebs.begin(),confs[0].ebs.end());


            
            testConf.profStride=std::max(1,testConf.sampleBlockSize/4);//todo: bugfix for others
            testConf.profiling = 1;

            size_t totalblock_num=1; 


            std::nth_element(ebs.begin(),ebs.begin()+k, ebs.end()); 

            double prof_abs_threshold = ebs[k];
            double sample_ratio = 5e-3;
            //std::cout<<"p1"<<std::endl;
            for(int i=0;i<N;i++){                      
                totalblock_num*=(size_t)((testConf.dims[i]-1)/testConf.sampleBlockSize);
            }
            if(N==2){
                QoZ::profiling_block_2d<T,N>(data[0],testConf.dims,starts,testConf.sampleBlockSize, prof_abs_threshold,testConf.profStride);
            }
            else if (N==3){
                QoZ::profiling_block_3d<T,N>(data[0],testConf.dims,starts,testConf.sampleBlockSize, prof_abs_threshold,testConf.profStride);
            } 
            //std::cout<<"p2"<<std::endl;
            size_t num_filtered_blocks=starts.size();
            if(num_filtered_blocks<=(int)(0.3*sample_ratio*totalblock_num))//todo: bugfix for others 
                testConf.profiling=0;

            for(auto j:{0,1,2}){

                QoZ::sampleBlocks<T,N>(data[j],testConf.dims,testConf.sampleBlockSize,sampled_blocks[j],sample_ratio,testConf.profiling,starts,false);
            }
            //std::cout<<"p3"<<std::endl;
            //std::cout<<sampled_blocks.size()<<std::endl;
            //std::cout<<sampled_blocks[0].size()<<std::endl;

            testConf.dims=std::vector<size_t>(N,testConf.sampleBlockSize+1);
            testConf.num=pow(testConf.sampleBlockSize+1,N);
            //std::cout<<sampled_blocks[0].size()<<" "<<sampled_blocks[0][0].size()<<" "<<testConf.num<<std::endl;
            //sample_block_ebs.resize(3);
            for(auto j:{0,1,2})
                sample_block_ebs[j].resize(sampled_blocks[0].size());

            for(int i=0;i<sampled_blocks[0].size();i++){
                for(auto j:{0,1,2}){
                    sample_block_ebs[j][i].resize(testConf.num);

                }   

                for (size_t k = 0; k < testConf.num; k++){
                    //std::cout<<" "<<i<<" "<<k<<std::endl;
                    std::array<T,3>cur_ebs = qoi->interpret_eb(sampled_blocks[0][i][k],sampled_blocks[1][i][k],sampled_blocks[2][i][k]);
                    for(auto j:{0,1,2})
                       sample_block_ebs[j][i][k]=cur_ebs[j];
                }

            }
            //std::cout<<"p4"<<std::endl;

        }
        else{
            
            for(auto j:{0,1,2})
                sampled_regions[j] = QoZ::sampling<T, N>(data[j], testConf.dims, sampling_num, sample_dims, sampling_block);
            testConf.setDims(sample_dims.begin(), sample_dims.end());

            for(auto j:{0,1,2})
                sample_region_ebs[j].resize(testConf.num);

            for (size_t i = 0; i < testConf.num; i++){
                std::array<T,3>cur_ebs = qoi->interpret_eb( sampled_regions[0][i],sampled_regions[1][i],sampled_regions[2][i]);
                for(auto j:{0,1,2})
                   sample_region_ebs[j][i]=cur_ebs[j];
            }
            

        }
    

        for(auto j:{0,1,2}){

            auto tmp_abs_eb = confs[j].absErrorBound;
            auto min_abs_eb = confs[j].absErrorBound;
            //QoZ::Config testConf = confs[i];

            

            //T *ebs = new T[conf.num];


            double best_abs_eb;

            ebs=std::vector<double>(confs[j].ebs.begin(),confs[j].ebs.end());
            std::cout<<j<<" "<<confs[j].ebs.front()<<" "<<confs[j].ebs.back()<<std::endl;

       
            std::vector<size_t> quantiles;
            //std::array<double,4> fixrate = {1.0,1.05,1.10,1.15};//or{1.0,1.1,1.2,1.3}
            
            //double quantile_split=0.1;
            for(auto i:{1.0,0.5,0.25,0.10,0.05,0.025,0.01})
                quantiles.push_back((size_t)(i*k));
            int quantile_num = quantiles.size();

            

            
            //std::sort(ebs.begin(),ebs.begin()+k+1);

            if(testConf.QoZ>0){
                if (testConf.maxStep==0){
                    std::array<size_t,4> anchor_strides={256,64,32,16};
                    testConf.maxStep = anchor_strides[N-1];
                }
                testConf.alpha = 1.5;
                testConf.beta = 2.0;
                
            }
            QoZ::Interp_Meta def_meta;
            testConf.interpMeta = def_meta;

            size_t best_quantile = 0;
            size_t last_quantile = confs[j].num;


            if(N==2 or N==3){



                double best_br = 9999;
                best_abs_eb = testConf.absErrorBound;
                                
                int idx = 0;
                for(auto quantile:quantiles)
                {   
                   
                    std::nth_element(ebs.begin(),ebs.begin()+quantile, ebs.begin()+last_quantile);

                    
                    testConf.absErrorBound = ebs[quantile];
                    //qoi->set_global_eb(testConf.absErrorBound);
                    // reset variables for average of square
                    //std::cout <<"quantile = "<<quantile<< " current_eb = " << testConf.absErrorBound<<std::endl;
                    //auto sampled_blocks_copy=;
                    double cur_br = CompressTest_QoI<T,N>(testConf,sampled_blocks[j],sample_block_ebs[j], qoi);        
                    std::cout << "current_eb = " << testConf.absErrorBound << ", current_br = " << cur_br << std::endl;
                    if(cur_br < best_br * 1.02){//todo: optimize
                        best_br = cur_br;
                        best_abs_eb = testConf.absErrorBound;
                        best_quantile = quantile;
                    }
                    else if(cur_br>1.1*best_br and testConf.early_termination){
                        break;
                    }

                    last_quantile = quantile+1;
                    idx++;
                    
                }

               

                if(best_quantile == quantiles.back()){
                    std::sort(ebs.begin(),ebs.begin()+best_quantile);
                    min_abs_eb = ebs[0];
                



                    size_t cur_quantile = best_quantile-1;
                    double init_best_abs_eb = best_abs_eb;
                    double min_ratio = 0.9;
                    while(cur_quantile>0){
                        
                        double temp_best_eb = ebs[cur_quantile];
                        double cur_ratio = min_ratio + (1.0-min_ratio) * ((double)cur_quantile/(double)best_quantile);
                        if (temp_best_eb/init_best_abs_eb<cur_ratio)
                            break;
                        else
                            best_abs_eb = temp_best_eb;
                        cur_quantile--;
                    }
                    best_quantile = cur_quantile + 1;
                }
                else{
                    min_abs_eb = *std::min_element(ebs.begin(),ebs.begin()+last_quantile);
                }


                if( (confs[j].qoiRegionMode == 0 or (confs[j].qoiRegionMode == 1 and confs[j].qoiRegionSize >= 3)) and  min_abs_eb >= 0.95 * best_abs_eb ){//may fix
                    confs[j].use_global_eb = true;
                    //conf.qoiPtr = qoi;
                }

                best_overall_br+=best_br/3;

            }


            else{
                std::array<double,7> fixrate = {1.06,1.04,1.02,1.00,1.00,1.00,1.00};


                testConf.ebs=sample_region_ebs[j];
                T * sampling_data = (T *) malloc(testConf.num * sizeof(T));




                testConf.qoiEBBase = testConf.absErrorBound / 1030;
                auto quantizer = QoZ::VariableEBLinearQuantizer<T, T>(testConf.quantbinCnt / 2);
                auto quantizer_eb = QoZ::EBLogQuantizer<T>(testConf.qoiEBBase, testConf.qoiEBLogBase, testConf.qoiQuantbinCnt / 2, testConf.absErrorBound);
                
                auto sz = QoZ::SZQoIInterpolationCompressor<T, N, QoZ::VariableEBLinearQuantizer<T, T>, QoZ::EBLogQuantizer<T>, QoZ::QoIEncoder<int>, QoZ::Lossless_zstd>(
                        quantizer, quantizer_eb, qoi, QoZ::QoIEncoder<int>(), QoZ::Lossless_zstd());



                double best_ratio = 0;
                best_abs_eb = testConf.absErrorBound;
                int idx = 0;
                for(auto quantile:quantiles)
                {   
                    
                    std::nth_element(ebs.begin(),ebs.begin()+quantile, ebs.begin()+last_quantile);

                    
                    testConf.absErrorBound = ebs[quantile];
                    //qoi->set_global_eb(testConf.absErrorBound);
                    size_t sampleOutSize;
                    memcpy(sampling_data, sampled_regions[j].data(), testConf.num * sizeof(T));
                    // reset variables for average of square
                    auto cmprData = sz.compress(testConf, sampling_data, sampleOutSize,1);
                    sz.clear();
                    delete[]cmprData;
                    double cur_ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;                
                    std::cout << "current_eb = " << testConf.absErrorBound << ", current_ratio = " << cur_ratio << std::endl;
                    double fr = fixrate[idx];
                    if(cur_ratio*fr>best_ratio){
                        best_ratio = cur_ratio*fr;
                        best_abs_eb = testConf.absErrorBound;
                        best_quantile = quantile;
                    }
                    else if(cur_ratio*fr<0.9*best_ratio and testConf.early_termination){
                        break;
                    }
                    last_quantile = quantile+1;

                    idx ++;
                }
                // set error bound
                
                free(sampling_data);

                if(best_quantile == quantiles.back()){
                    std::sort(ebs.begin(),ebs.begin()+best_quantile);
                    min_abs_eb=ebs[0];
                



                    size_t cur_quantile = best_quantile-1;
                    double init_best_abs_eb = best_abs_eb;
                    double min_ratio = 0.9;
                    while(cur_quantile>0){
                        
                        double temp_best_eb = ebs[cur_quantile];
                        double cur_ratio = min_ratio + (1.0-min_ratio) * ((double)cur_quantile/(double)best_quantile);
                        if (temp_best_eb/init_best_abs_eb<cur_ratio)
                            break;
                        else
                            best_abs_eb = temp_best_eb;
                        cur_quantile--;
                    }
                    best_quantile = cur_quantile + 1;
                }
                else{
                    min_abs_eb = (*std::min_element(ebs.begin(),ebs.begin()+last_quantile));
                }



                if( (confs[j].qoiRegionMode == 0 or (confs[j].qoiRegionMode == 1 and confs[j].qoiRegionSize >= 3)) and min_abs_eb >= 0.95 * best_abs_eb ){//may fix
                    confs[j].use_global_eb = true;
                    //conf.qoiPtr = qoi;
                }
                best_overall_br+=8*sizeof(T)/(3*best_ratio);


            }
                
                
            std::cout<<"Selected quantile: "<<(double)best_quantile/(double)confs[0].num<<" best overall bitrate: "<<best_overall_br<<std::endl;
            best_abs_ebs[j]=best_abs_eb;
            
        }
        if((confs[0].qoiRegionMode == 0 or (confs[0].qoiRegionMode == 1 and confs[0].qoiRegionSize >= 3)) and not (confs[0].use_global_eb and confs[1].use_global_eb and confs[2].use_global_eb) ){//global tuning
            double cur_overall_br = 0.0;
            if(N==2 or N==3){

                auto sampled_blocks_dec = sampled_blocks;
                
                for(auto j:{0,1,2}){

                    testConf.absErrorBound = best_abs_ebs[j];
                    testConf.use_global_eb = true;
                    //qoi->set_global_eb(testConf.absErrorBound);
                    // reset variables for average of square
                    //std::cout <<"quantile = "<<quantile<< " current_eb = " << testConf.absErrorBound<<std::endl;
                    std::pair<double,double> results = CompressTest<T,N>(testConf,sampled_blocks_dec[j],QoZ::ALGO_INTERP,QoZ::TUNING_TARGET_CR); 
                    double cur_br = results.first; 


                    cur_overall_br += cur_br/3;    

                    
                }
                std::array<std::vector<T>,3>offsets;
                for(size_t i=0;i<sampled_blocks[0].size();i++){
                    for(size_t j=0;j<testConf.num;j++){
                        if ( qoi->check_compliance(sampled_blocks[0][i][j],sampled_blocks[1][i][j],sampled_blocks[2][i][j],
                                                   sampled_blocks_dec[0][i][j],sampled_blocks_dec[1][i][j],sampled_blocks_dec[2][i][j])){
                            for (auto l:{0,1,2})
                                offsets[l].push_back(0);
                        }
                        else{
                            std::cout<<"incomp "<<i<<" "<<j<<std::endl;
                            for (auto l:{0,1,2})
                                offsets[l].push_back(sampled_blocks[l][i][j]-sampled_blocks_dec[l][i][j]);
                        }
                    }
                }

                auto zstd = QoZ::Lossless_zstd();
                size_t overall_offset_size;
                
                for (auto i:{0,1,2}){
                    size_t offset_size;
                    QoZ::uchar *lossless_data = zstd.compress(reinterpret_cast< QoZ::uchar *>(offsets[i].data()),
                                                                 offsets[i].size()*sizeof(T),
                                                                 offset_size);
                    offsets[i].clear();
                    overall_offset_size+=offset_size;
                    delete []lossless_data;

                }
                cur_overall_br += (double)(overall_offset_size*8)/(double)(3*sampled_blocks[0].size()*testConf.num);
                
            }
            else{
                std::array<T*,3>sampling_data_dec;
                    
                for(auto j:{0,1,2}){
                    sampling_data_dec[j] = (T *) malloc(testConf.num * sizeof(T));
                    testConf.absErrorBound = best_abs_ebs[j];
                    testConf.use_global_eb = true;
                    //testConf.qoiPtr = qoi;

                    size_t sampleOutSize;
                    memcpy(sampling_data_dec[j], sampled_regions[j].data(), testConf.num * sizeof(T));
                    // reset variables for average of square
                    auto sz =  QoZ::SZInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
                        QoZ::LinearQuantizer<T>(testConf.absErrorBound),
                        QoZ::HuffmanEncoder<int>(),
                        QoZ::Lossless_zstd());
                    auto cmprData = sz.compress(testConf, sampling_data_dec[j], sampleOutSize,0);
                    delete[]cmprData;
                    double cur_br = (double)sampleOutSize*8/ (double)sampling_num; 
                    cur_overall_br += cur_br/3;              
                    //double fr = fixrate[idx];
                   
                }

                std::array<std::vector<T>,3>offsets;

                for(size_t i=0;i<testConf.num;i++){
                   
                    if (qoi->check_compliance(sampled_regions[0][i],sampled_regions[1][i],sampled_regions[2][i],
                                                   sampling_data_dec[0][i],sampling_data_dec[1][i],sampling_data_dec[2][i])){
                        for (auto j:{0,1,2})
                            offsets[j].push_back(0);
                    }
                    else{
                        for (auto j:{0,1,2})
                            offsets[j].push_back(sampled_regions[j][i]-sampling_data_dec[j][i]);
                    }
                    
                }
                for(auto j:{0,1,2})
                    free(sampling_data_dec[j]);

                auto zstd = QoZ::Lossless_zstd();
                size_t overall_offset_size;
                
                for (auto i:{0,1,2}){
                    size_t offset_size;
                    QoZ::uchar *lossless_data = zstd.compress(reinterpret_cast< QoZ::uchar *>(offsets[i].data()),
                                                                 offsets[i].size()*sizeof(T),
                                                                 offset_size);
                    offsets[i].clear();
                    overall_offset_size+=offset_size;
                    delete []lossless_data;

                }
                cur_overall_br += (double)(overall_offset_size*8)/(double)(3*testConf.num);


                



            }

            std::cout << "Global test, current_br = " << cur_overall_br << std::endl;
            if(cur_overall_br < best_overall_br * 1.02){//todo: optimize
                best_overall_br = cur_overall_br;
                for(auto j:{0,1,2})
                    confs[j].use_global_eb=true;

            }
        }
    }
    else{
        for(auto j:{0,1,2}){
            ebs =std::vector<double>(confs[j].ebs.begin(),confs[j].ebs.end());
            std::nth_element(ebs.begin(),ebs.begin()+k, ebs.end());
            best_abs_ebs[j] = ebs[k];
            
        }

    }

    //qoi->set_global_eb(best_abs_eb);

    
    //std::cout<<"Smaller ebs: "<<smaller_ebs_ratio<<std::endl;

    for(auto j:{0,1,2}){
        std::cout<<"File "<<j+1<<":"<<std::endl;
        if(confs[j].use_global_eb)
            std::cout<<"Use global eb."<<std::endl; 

        
        std::cout << "Best abs eb : " << best_abs_ebs[j] << std::endl; 
        //std::cout << best_abs_eb << " " << tmp_abs_eb << std::endl;
        confs[j].absErrorBound = best_abs_ebs[j];
        //qoi->set_global_eb(best_abs_eb);
       // conf.setDims(dims.begin(), dims.end());
        // reset dimensions and variables for average of square
        //if(conf.qoi == 3){
         //   qoi->set_dims(dims);
        //    qoi->init();
        //}
        for (size_t i = 0; i < confs[j].num; i++){
            if(confs[j].ebs[i] > best_abs_ebs[j])
                confs[j].ebs[i] = best_abs_ebs[j];
        }

        confs[j].qoiEBBase = confs[j].absErrorBound / 1030;
        //std::cout << conf.qoi << " " << conf.qoiEB << " " << conf.qoiEBBase << " " << conf.qoiEBLogBase << " " << conf.qoiQuantbinCnt << std::endl;
        confs[j].qoi_tuned = true;
    }
    //std::cout<<confs[0].qoi<<" "<<confs[1].qoi<<" "<<confs[2].qoi<<" "<<std::endl;
    return;
        
}
    


std::pair <double,double> setABwithRelBound(double rel_bound,int configuration=0){

    double cur_alpha=-1,cur_beta=-1;
    if(configuration==0){              
        if (rel_bound>=0.01){
            cur_alpha=2;
            cur_beta=2;
        }
        else if (rel_bound>=0.007){
            cur_alpha=1.75;
            cur_beta=2;
        }                 
        else if (rel_bound>=0.004){
            cur_alpha=1.5;
            cur_beta=2;
        }                   
        else if (rel_bound>0.001){
            cur_alpha=1.25;
            cur_beta=1.5;
        }
        else {
            cur_alpha=1;
            cur_beta=1;
        }
    }
    else if(configuration==1){                
        if (rel_bound>=0.01){
            cur_alpha=2;
            cur_beta=4;
        }
        else if (rel_bound>=0.007){
            cur_alpha=1.75;
            cur_beta=3;
        }
        else if (rel_bound>=0.004){
            cur_alpha=1.5;
            cur_beta=2;
        }           
        else if (rel_bound>0.001){
            cur_alpha=1.25;
            cur_beta=1.5;
        }
        else {
                cur_alpha=1;
                cur_beta=1;
            }
    }
    else if(configuration==2){                
        if (rel_bound>=0.01){
            cur_alpha=2;
            cur_beta=4;
        }
        else if (rel_bound>=0.007){
            cur_alpha=1.75;
            cur_beta=3;
        }                    
        else if (rel_bound>=0.004){
            cur_alpha=1.5;
            cur_beta=2;
        }
        else if (rel_bound>0.001){
            cur_alpha=1.5;
            cur_beta=1.5;
        }
        else if (rel_bound>0.0005){
            cur_alpha=1.25;
            cur_beta=1.5;
        }
        else {
            cur_alpha=1;
            cur_beta=1;
        }
    }
    return std::pair<double,double>(cur_alpha,cur_beta);
}

void setLorenzoFixRates(QoZ::Config &conf,double rel_bound){
    double e1=1e-5;
    double e2=1e-4;
    double e3=1e-3;
    //double e4=1e-1;
    /*
    double f1=conf.sampleBlockSize>=64?2: 1;
    // double f2=1.1;old
    double f2=conf.sampleBlockSize>=64?3:1.2; 
    // double f3=1.2;//old need to raise 
    //double f3=1.3;
    double f3=conf.sampleBlockSize>=64?4:1.3;
    */
    double f1=conf.sampleBlockSize>=64?2: 1.1;//updated from 1.0, not tested
    // double f2=1.1;old
    double f2=conf.sampleBlockSize>=64?3:1.2; 
    // double f3=1.2;//old need to raise 
    //double f3=1.3;
    double f3=conf.sampleBlockSize>=64?4:1.3;
    if(rel_bound<=e1)
        conf.lorenzoBrFix=f1;
    else if(rel_bound<=e2)
        conf.lorenzoBrFix=f1-(f1-f2)*(rel_bound-e1)/(e2-e1);
    else if (rel_bound<=e3)
        conf.lorenzoBrFix=f2-(f2-f3)*(rel_bound-e2)/(e3-e2);
    else 
        conf.lorenzoBrFix=f3;
}

template<class T, QoZ::uint N>
double Tuning(QoZ::Config &conf, T *data){
   
    T rng=conf.rng;
    double rel_bound = conf.relErrorBound>0?conf.relErrorBound:conf.absErrorBound/rng;
    //if(rel_bound>1e-3 or conf.tuningTarget==QoZ::TUNING_TARGET_SSIM)//rencently changed, need to fix later
    //    conf.testLorenzo=0;
   // QoZ::Timer timer(true);
    //timer.stop("")
    if(conf.QoZ>0){
        
        //testLorenzo?
        

        //activate
        //conf.testLorenzo=0;//temp deactivated Lorenzo. need further revision
        conf.profiling=1;
        if(conf.autoTuningRate<=0)
            conf.autoTuningRate = (N<=2?0.01:0.005);
        if(conf.predictorTuningRate<=0)
            conf.predictorTuningRate = (N<=2?0.01:0.005);
        if (conf.maxStep<=0){
            std::array<size_t,4> anchor_strides={256,64,32,16};
            conf.maxStep = anchor_strides[N-1];
        }
        if (conf.levelwisePredictionSelection<=0)
            conf.levelwisePredictionSelection = (N<=2?5:4);
        if (conf.sampleBlockSize<=0){
            
            conf.sampleBlockSize = (N<=2?64:32);
        }

        if(conf.QoZ>=2){
            //conf.testLorenzo=1;
            conf.multiDimInterp=1;
            conf.naturalSpline=1;
            conf.fullAdjacentInterp=1;
            conf.freezeDimTest=1;
            
        }
        if(conf.QoZ>=3){
            conf.dynamicDimCoeff=1;
        }
        if(conf.QoZ>=4){
            conf.crossBlock=1;
            conf.blockwiseTuning=1;
            if(conf.blockwiseSampleRate<1.0)
                conf.blockwiseSampleRate=3.0;
        }
    }   
    //Deactivate several features when dim not fit
    if(N!=2 and N!=3){
       // conf.QoZ=0; comment it for level-wise eb
        conf.autoTuningRate=0;
        conf.predictorTuningRate=0;
        conf.levelwisePredictionSelection=0;
        conf.multiDimInterp=0;
        conf.naturalSpline=0;
        conf.fullAdjacentInterp=0;
        conf.freezeDimTest=0;
        conf.dynamicDimCoeff=0;
        conf.blockwiseTuning=0;
    }
    

    if(conf.multiDimInterp==0)
        conf.dynamicDimCoeff=0;
    size_t sampling_num, sampling_block;
    double best_interp_cr=0.0;
    double best_lorenzo_ratio=0.0;
    bool useInterp=true;        
    std::vector<size_t> sample_dims(N);
    std::vector<T> sampling_data;
    double anchor_rate=0;
    int max_interp_level = -1;

    for (size_t i = 0; i < N; i++) {
        if ( max_interp_level < ceil(log2(conf.dims[i]))) {
             max_interp_level = (uint) ceil(log2(conf.dims[i]));
        }
                
    }
    
    if (conf.maxStep>0){
        anchor_rate=1/(pow(conf.maxStep,N));   
        int temp_max_interp_level=(uint)log2(conf.maxStep);//to be catious: the max_interp_level is different from the ones in szinterpcompressor, which includes the level of anchor grid.
        if (temp_max_interp_level<=max_interp_level){                  
            max_interp_level=temp_max_interp_level;
        }
        if (conf.levelwisePredictionSelection>max_interp_level)
            conf.levelwisePredictionSelection=max_interp_level;
    }
            
    

    std::vector<QoZ::Interp_Meta> bestInterpMeta_list(conf.levelwisePredictionSelection);
    QoZ::Interp_Meta bestInterpMeta;


    size_t shortest_edge=conf.dims[0];
    for (size_t i=0;i<N;i++){
        shortest_edge=conf.dims[i]<shortest_edge?conf.dims[i]:shortest_edge;
    }
 

    
    if (conf.sampleBlockSize<=0){
        conf.sampleBlockSize = (N==2?64:32);
            
    }

    
    
    size_t minimum_sbs=8;
    if (conf.sampleBlockSize<minimum_sbs)
        conf.sampleBlockSize=minimum_sbs;


    while(conf.sampleBlockSize>=shortest_edge)
        conf.sampleBlockSize/=2;


    

    while(conf.autoTuningRate>0 and conf.sampleBlockSize>=2*minimum_sbs and (pow(conf.sampleBlockSize+1,N)/(double)conf.num)>1.5*conf.autoTuningRate)
        conf.sampleBlockSize/=2;

    if (conf.sampleBlockSize<minimum_sbs){
        conf.predictorTuningRate=0.0;
        conf.autoTuningRate=0.0;
    }
    else{
        int max_lps_level=(uint)log2(conf.sampleBlockSize);//to be catious: the max_interp_level is different from the ones in szinterpcompressor, which includes the level of anchor grid.

        if (conf.levelwisePredictionSelection>max_lps_level)
            conf.levelwisePredictionSelection=max_lps_level;
    }

    std::vector< std::vector<T> > sampled_blocks;
    size_t sampleBlockSize=conf.sampleBlockSize;
    size_t num_sampled_blocks;
    size_t per_block_ele_num;
    size_t ele_num;

    
           
    size_t totalblock_num=1;  
    for(int i=0;i<N;i++){                      
        totalblock_num*=(size_t)((conf.dims[i]-1)/sampleBlockSize);
    }

    std::vector<std::vector<size_t> >starts;
    if((conf.autoTuningRate>0 or conf.predictorTuningRate>0) and conf.profiling){      
        conf.profStride=conf.sampleBlockSize/4;
        if(N==2){
            QoZ::profiling_block_2d<T,N>(data,conf.dims,starts,sampleBlockSize,conf.absErrorBound,conf.profStride);
        }
        else if (N==3){
            QoZ::profiling_block_3d<T,N>(data,conf.dims,starts,sampleBlockSize,conf.absErrorBound,conf.profStride);
        }
       
    }


    size_t num_filtered_blocks=starts.size();
    if(num_filtered_blocks<=(int)(0.3*conf.predictorTuningRate*totalblock_num))//temp. to refine
        conf.profiling=0;
    double profiling_coeff=1;//It seems that this coefficent is useless. Need further test
  
    if(conf.profiling){
        profiling_coeff=((double)num_filtered_blocks)/(totalblock_num);
    }
    std::vector<size_t> global_dims=conf.dims;
    size_t global_num=conf.num;
    if(conf.autoTuningRate>0){
        if (conf.testLorenzo>0)
            setLorenzoFixRates(conf,rel_bound);
    }
    

    bool blockwiseTuning=conf.blockwiseTuning;
    conf.blockwiseTuning=false;

    if (conf.predictorTuningRate>0 and conf.predictorTuningRate<1){
        if (conf.verbose)
            std::cout<<"Predictor tuning started."<<std::endl;
        double o_alpha=conf.alpha;
        double o_beta=conf.beta;
                    
        
        QoZ::sampleBlocks<T,N>(data,conf.dims,sampleBlockSize,sampled_blocks,conf.predictorTuningRate,conf.profiling,starts,conf.var_first);
              
        num_sampled_blocks=sampled_blocks.size();
        per_block_ele_num=pow(sampleBlockSize+1,N);
        ele_num=num_sampled_blocks*per_block_ele_num;
        conf.dims=std::vector<size_t>(N,sampleBlockSize+1);
        conf.num=per_block_ele_num;
        std::vector<T> cur_block(per_block_ele_num,0);
        
        
        
            
        double ori_eb=conf.absErrorBound;
        std::vector<size_t> coeffs_size;
        /*
        if(conf.testLorenzo and conf.autoTuningRate==0 ){

            std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks,QoZ::ALGO_LORENZO_REG,QoZ::TUNING_TARGET_CR,false);
            best_lorenzo_ratio=sizeof(T)*8.0/results.first;
            
            if(conf.verbose)
                std::cout << "lorenzo best cr = " << best_lorenzo_ratio << std::endl;

        }*/

        if (conf.autoTuningRate>0){

           
            std::pair<double,double> ab=setABwithRelBound(rel_bound,0);//ori pdtuningqabconf
            conf.alpha=ab.first;
            conf.beta=ab.second;
       
        }
        std::vector<uint8_t> interpAlgo_Candidates={QoZ::INTERP_ALGO_LINEAR, QoZ::INTERP_ALGO_CUBIC};
        if(conf.quadInterp){//deprecated
            interpAlgo_Candidates.push_back(QoZ::INTERP_ALGO_QUAD);
        }

        
        std::vector<uint8_t> interpParadigm_Candidates={0};//
        std::vector<uint8_t> cubicSplineType_Candidates={0};
        std::vector<uint8_t> interpDirection_Candidates={0, (uint8_t)(QoZ::factorial(N) -1)};
       
        std::vector<uint8_t> adjInterp_Candidates={0};//

        
        if(conf.multiDimInterp>0){
            for(size_t i=1;i<=conf.multiDimInterp;i++)
                interpParadigm_Candidates.push_back(i);
        }

        if (conf.naturalSpline){
            cubicSplineType_Candidates.push_back(1);
        }


    
        
        if(conf.fullAdjacentInterp){
            adjInterp_Candidates.push_back(1);
        }
        
        if(conf.levelwisePredictionSelection>0){
            std::vector<QoZ::Interp_Meta> interpMeta_list(conf.levelwisePredictionSelection);
            auto sz = QoZ::SZInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
                                    QoZ::LinearQuantizer<T>(conf.absErrorBound),
                                    QoZ::HuffmanEncoder<int>(),
                                    QoZ::Lossless_zstd());   
            double best_accumulated_interp_loss_1=0;
            double best_accumulated_interp_loss_2=0;
            std::vector<std::vector<double> > linear_interp_vars(conf.levelwisePredictionSelection),cubic_noknot_vars(conf.levelwisePredictionSelection),cubic_nat_vars(conf.levelwisePredictionSelection);
            
            for(int level=conf.levelwisePredictionSelection;level>0;level--){
              
                int start_level=(level==conf.levelwisePredictionSelection?9999:level);
                int end_level=level-1;
                
                if((conf.multiDimInterp>0 and conf.dynamicDimCoeff>0) or (conf.freezeDimTest>0 and level==1 and N>=3) ){
                    
                    
                    size_t interp_stride=pow(2,level-1);
                    size_t stride;
                    double cur_eb=level>=3?conf.absErrorBound/2:conf.absErrorBound;
                    if(level>=3)
                        stride=2;
                    else if(level>=2)
                        stride= conf.adaptiveMultiDimStride<=4?2:conf.adaptiveMultiDimStride/2;
                    else
                        stride= conf.adaptiveMultiDimStride;
                    if(conf.multiDimInterp>0 and conf.dynamicDimCoeff){
                        QoZ::calculate_interp_error_vars<T,N>(data, global_dims,linear_interp_vars[level-1],0,0,stride,interp_stride,cur_eb);
                        QoZ::preprocess_vars<N>(linear_interp_vars[level-1]);
                       
                        if (conf.naturalSpline){
                            QoZ::calculate_interp_error_vars<T,N>(data, global_dims,cubic_nat_vars[level-1],1,1,stride,interp_stride,cur_eb);
                            QoZ::preprocess_vars<N>(cubic_nat_vars[level-1]);
                         
                        }
                    }
                    QoZ::calculate_interp_error_vars<T,N>(data, global_dims,cubic_noknot_vars[level-1],1,0,stride,interp_stride,cur_eb);
                    QoZ::preprocess_vars<N>(cubic_noknot_vars[level-1]);
                   
                }


                double best_interp_absloss=std::numeric_limits<double>::max();
                //conf.cmprAlgo = QoZ::ALGO_INTERP;    
                QoZ::Interp_Meta cur_meta;
                QoZ::Interp_Meta best_meta;

                for (auto &interp_op: interpAlgo_Candidates) {
                    cur_meta.interpAlgo=interp_op;
                    for (auto &interp_pd: interpParadigm_Candidates) {
                        cur_meta.interpParadigm=interp_pd;

                        

                        for (auto &interp_direction: interpDirection_Candidates) {
                            if ((interp_pd==1 or  (interp_pd==2 and N<=2)) and interp_direction!=0)
                                continue;
                            cur_meta.interpDirection=interp_direction;
                            for(auto &cubic_spline_type:cubicSplineType_Candidates){
                                if (interp_op!=QoZ::INTERP_ALGO_CUBIC and cubic_spline_type!=0)
                                    break;
                                cur_meta.cubicSplineType=cubic_spline_type;
                                for(auto adj_interp:adjInterp_Candidates){
                                    if (interp_op!=QoZ::INTERP_ALGO_CUBIC and adj_interp!=0)
                                        break;
                                    
                                    cur_meta.adjInterp=adj_interp;
                                    
                                    if(conf.dynamicDimCoeff>0 and interp_pd>0){
                                        if(interp_op==0){
                                            for(size_t i=0;i<N;i++)
                                                cur_meta.dimCoeffs[i]=linear_interp_vars[level-1][i];
                                        }
                                        else if (cubic_spline_type==0){
                                            for(size_t i=0;i<N;i++)
                                                cur_meta.dimCoeffs[i]=cubic_noknot_vars[level-1][i];
                                        }
                                        else{
                                            for(size_t i=0;i<N;i++)
                                                cur_meta.dimCoeffs[i]=cubic_nat_vars[level-1][i];
                                        }

                                    }
                                    
                                    conf.interpMeta=cur_meta;


                                    double cur_absloss=0;
                                    for (int i=0;i<num_sampled_blocks;i++){
                                        cur_block=sampled_blocks[i];  //not so efficient              
                                        size_t outSize=0;                              
                                        auto cmprData =sz.compress(conf, cur_block.data(), outSize,2,start_level,end_level);
                                        delete []cmprData;                              
                                        cur_absloss+=conf.decomp_square_error;
                                    }
                                    
                                    if (cur_absloss<best_interp_absloss){
                                        best_meta=cur_meta;
                                        best_interp_absloss=cur_absloss;
                                    }
                                    cur_meta.dimCoeffs={1.0/3.0,1.0/3.0,1.0/3.0};
                    
                                }
                            }
                        }   
                    }
                }
                best_accumulated_interp_loss_1+=best_interp_absloss;
      

                interpMeta_list[level-1]=best_meta;
                    
                if(conf.pdTuningRealComp){
                    //place to add real compression,need to deal the problem that the sampled_blocks are changed. 
                    
                    conf.interpMeta=best_meta;
                    for (int i=0;i<num_sampled_blocks;i++){

                        size_t outSize=0;
                                   
                        auto cmprData =sz.compress(conf, sampled_blocks[i].data(), outSize,2,start_level,end_level);
                        delete []cmprData;
                    }
                    
                } 
                

            }
            
            bestInterpMeta_list=interpMeta_list;

           

            
            if(conf.pdTuningRealComp and ((conf.autoTuningRate>0 and conf.autoTuningRate==conf.predictorTuningRate) or (conf.adaptiveMultiDimStride>0 and N==3))){
                    //recover sample if real compression used                  
                QoZ::sampleBlocks<T,N>(data,global_dims,sampleBlockSize,sampled_blocks,conf.predictorTuningRate,conf.profiling,starts,conf.var_first);
            }
            
                

            //frozendim
            
            if(conf.freezeDimTest and N>=3 ){

                std::vector<QoZ::Interp_Meta> tempmeta_list=conf.interpMeta_list;
                conf.interpMeta_list=interpMeta_list;  
                auto sampled_blocks_copy=sampled_blocks;
                std::pair<double,double> results=CompressTest<T,N>(conf,sampled_blocks_copy,QoZ::ALGO_INTERP,QoZ::TUNING_TARGET_CR,false);
                double best_interp_cr_1=sizeof(T)*8.0/results.first;     
                conf.interpMeta_list=tempmeta_list;
                




                size_t frozen_dim=0;
                double cur_weight=cubic_noknot_vars[0][0];
                for(size_t i=1;i<N;i++){
                    if(cubic_noknot_vars[0][i]<cur_weight){
                        frozen_dim=i;
                        cur_weight=cubic_noknot_vars[0][i];
                    }
                }
                if(frozen_dim==0)
                    interpDirection_Candidates={6,7};
                else if (frozen_dim==1)
                    interpDirection_Candidates={8,9};
                else
                    interpDirection_Candidates={10,11};


                for(int level=conf.levelwisePredictionSelection;level>0;level--){
                    int start_level=(level==conf.levelwisePredictionSelection?9999:level);
                    int end_level=level-1;
                  
                    double best_interp_absloss=std::numeric_limits<double>::max();
        
                    QoZ::Interp_Meta cur_meta;
                    QoZ::Interp_Meta best_meta;

                    for (auto &interp_op: interpAlgo_Candidates) {
                        cur_meta.interpAlgo=interp_op;
                        for (auto &interp_pd: interpParadigm_Candidates) {
                            if (interp_pd>1)
                                continue;
                            cur_meta.interpParadigm=interp_pd;

                            

                            for (auto &interp_direction: interpDirection_Candidates) {
                                if (interp_pd>=1  and interp_direction%2!=0)//only dims[0] matters.
                                    continue;
                                cur_meta.interpDirection=interp_direction;
                                for(auto &cubic_spline_type:cubicSplineType_Candidates){
                                    if (interp_op!=QoZ::INTERP_ALGO_CUBIC and cubic_spline_type!=0)
                                        break;
                                    cur_meta.cubicSplineType=cubic_spline_type;
                                    for(auto adj_interp:adjInterp_Candidates){
                                        if (interp_op!=QoZ::INTERP_ALGO_CUBIC and adj_interp!=0)
                                            break;
                                       
                                        cur_meta.adjInterp=adj_interp;

                                        if(conf.dynamicDimCoeff>0 and interp_pd>0){
                                            if(interp_op==0){
                                                for(size_t i=0;i<N;i++)
                                                    cur_meta.dimCoeffs[i]=linear_interp_vars[level-1][i];
                                            }
                                            else if (cubic_spline_type==0){
                                                for(size_t i=0;i<N;i++)
                                                    cur_meta.dimCoeffs[i]=cubic_noknot_vars[level-1][i];
                                            }
                                            else{
                                                for(size_t i=0;i<N;i++)
                                                    cur_meta.dimCoeffs[i]=cubic_nat_vars[level-1][i];
                                            }

                                        }
                                        
                                        conf.interpMeta=cur_meta;

                                        
                                        double cur_absloss=0;
                                        for (int i=0;i<num_sampled_blocks;i++){
                                            cur_block=sampled_blocks[i];  //not so efficient              
                                            size_t outSize=0;                              
                                            auto cmprData =sz.compress(conf, cur_block.data(), outSize,2,start_level,end_level);
                                            delete []cmprData;                              
                                            cur_absloss+=conf.decomp_square_error;
                                        }
                                        if (cur_absloss<best_interp_absloss){
                                            best_meta=cur_meta;
                                            best_interp_absloss=cur_absloss;
                                        }
                                        cur_meta.dimCoeffs={1.0/3.0,1.0/3.0,1.0/3.0};
                                        
                                    }
                                }
                            }   
                        }
                    }
                    best_accumulated_interp_loss_2+=best_interp_absloss;
          
                 

                    interpMeta_list[level-1]=best_meta;
                        
                    if(conf.pdTuningRealComp){
                        //place to add real compression,need to deal the problem that the sampled_blocks are changed. 
                  
                        conf.interpMeta=best_meta;
                        for (int i=0;i<num_sampled_blocks;i++){

                            size_t outSize=0;
                                       
                            auto cmprData =sz.compress(conf, sampled_blocks[i].data(), outSize,2,start_level,end_level);
                            delete []cmprData;
                        }
                        
                    } 
                }
                

                tempmeta_list=conf.interpMeta_list;
                conf.interpMeta_list=interpMeta_list;      
                sampled_blocks_copy=sampled_blocks;
                results=CompressTest<T,N>(conf,sampled_blocks_copy,QoZ::ALGO_INTERP,QoZ::TUNING_TARGET_CR,false);
                double best_interp_cr_2=sizeof(T)*8.0/results.first;     
                conf.interpMeta_list=tempmeta_list;

                if(best_interp_cr_2>best_interp_cr_1*1.05){
                    conf.frozen_dim=frozen_dim;
                    bestInterpMeta_list=interpMeta_list;
                    std::cout<<"Dim "<<frozen_dim<<" frozen"<<std::endl;
                }
            

                if(conf.pdTuningRealComp and conf.autoTuningRate>0 and conf.autoTuningRate==conf.predictorTuningRate){
                        //recover sample if real compression used                  
                    QoZ::sampleBlocks<T,N>(data,global_dims,sampleBlockSize,sampled_blocks,conf.predictorTuningRate,conf.profiling,starts,conf.var_first);
                }



            }
            
            



            
            conf.interpMeta_list=bestInterpMeta_list;
            if(conf.autoTuningRate==0){ //Qustionable part.  //when adaptivemdtride>0 there's a duplication of work. To fix.
              
                auto sampled_blocks_copy=sampled_blocks;
                std::pair<double,double> results=CompressTest<T,N>(conf,sampled_blocks_copy,QoZ::ALGO_INTERP,QoZ::TUNING_TARGET_CR,false);
                double cur_best_interp_cr=sizeof(T)*8.0/results.first;     
               
                best_interp_cr=cur_best_interp_cr;
                   
                   
            }
           
        }

        else{
            
            QoZ::Interp_Meta best_meta,cur_meta;
            double cur_best_interp_cr=0.0;
            for (auto &interp_op: interpAlgo_Candidates) {
                cur_meta.interpAlgo=interp_op;
                for (auto &interp_pd: interpParadigm_Candidates) {
                    cur_meta.interpParadigm=interp_pd;
                    for (auto &interp_direction: interpDirection_Candidates) {
                        if (interp_pd==1 or  (interp_pd==2 and N<=2) and interp_direction!=0)
                            continue;
                        cur_meta.interpDirection=interp_direction;
                        for(auto &cubic_spline_type:cubicSplineType_Candidates){
                            if (interp_op!=QoZ::INTERP_ALGO_CUBIC and cubic_spline_type!=0)
                                break;
                            cur_meta.cubicSplineType=cubic_spline_type;
                            for(auto adj_interp:adjInterp_Candidates){
                                if (interp_op!=QoZ::INTERP_ALGO_CUBIC and adj_interp!=0)
                                    break;
                                cur_meta.adjInterp=adj_interp;       
                                conf.interpMeta=cur_meta;
                                double cur_ratio=0;
                                auto sampled_blocks_copy=sampled_blocks;
                                std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks_copy,QoZ::ALGO_INTERP,QoZ::TUNING_TARGET_CR,false);
                                cur_ratio=sizeof(T)*8.0/results.first;
                                
                                if (cur_ratio>cur_best_interp_cr){
                                    cur_best_interp_cr=cur_ratio;
                                    
                                    best_meta=cur_meta;

                                }
                            }
                        }
                    }
                }
            }
           
            bestInterpMeta=best_meta;
            conf.interpMeta=best_meta;
            if(conf.autoTuningRate==0){
                if(cur_best_interp_cr>best_interp_cr){
                    
                    best_interp_cr=cur_best_interp_cr;
                    conf.interpMeta=best_meta;
                
                }
            }
        }
        conf.absErrorBound=ori_eb;


        

        if(conf.verbose)           
            printf("Predictor tuning finished.\n");           
        conf.alpha=o_alpha;
        conf.beta=o_beta;
        conf.dims=global_dims;
        conf.num=global_num;
        useInterp= (best_interp_cr>=best_lorenzo_ratio) or best_lorenzo_ratio>=80 or best_interp_cr>=80;//orig 0.95*lorenzo_ratio
        if(conf.verbose){
            if (conf.levelwisePredictionSelection<=0){
                std::cout << "interp best interpAlgo = " << (bestInterpMeta.interpAlgo == 0 ? "LINEAR" : (bestInterpMeta.interpAlgo == 1?"CUBIC":"QUAD")) << std::endl;
                std::cout << "interp best interpParadigm = " << (bestInterpMeta.interpParadigm == 0 ? "1D" : (bestInterpMeta.interpParadigm == 1 ? "MD" : "HD") ) << std::endl;
                if(bestInterpMeta.interpParadigm!=1)
                    std::cout << "interp best direction = " << (unsigned) bestInterpMeta.interpDirection << std::endl;
                if(bestInterpMeta.interpAlgo!=0){
                    std::cout << "interp best cubic spline = " << (unsigned) bestInterpMeta.cubicSplineType << std::endl;
                    std::cout << "interp best adj = " << (unsigned) bestInterpMeta.adjInterp << std::endl;

                }
            }
            else{
                for(int level=conf.levelwisePredictionSelection;level>0;level--){
                    std::cout << "Level: " << (unsigned) level<<std::endl;
                    std::cout << "\tinterp best interpAlgo = " << (bestInterpMeta_list[level-1].interpAlgo == 0 ? "LINEAR" : (bestInterpMeta_list[level-1].interpAlgo == 1 ? "CUBIC" : "QUAD")) << std::endl;
                    std::cout << "\tinterp best interpParadigm = " << (bestInterpMeta_list[level-1].interpParadigm == 0 ? "1D" : (bestInterpMeta_list[level-1].interpParadigm == 1 ? "MD" : "HD") ) << std::endl;
                    if(bestInterpMeta_list[level-1].interpParadigm!=1)
                        std::cout << "\tinterp best direction = " << (unsigned) bestInterpMeta_list[level-1].interpDirection << std::endl;
                    if(bestInterpMeta_list[level-1].interpAlgo!=0){
                        std::cout << "\tinterp best cubic spline = " << (unsigned) bestInterpMeta_list[level-1].cubicSplineType << std::endl;
                        std::cout << "\tinterp best adj = " << (unsigned) bestInterpMeta_list[level-1].adjInterp << std::endl;

                    }
                }
            }
            if(conf.autoTuningRate==0){
                std::cout << "interp best cr = " << best_interp_cr << std::endl;
                printf("choose %s\n", useInterp ? "interp" : "Lorenzo");
            }
        }

    }
    
    else{// if (!conf.blockwiseTuning){ //recently modified. not sure.

        QoZ::Timer timer(true);
    
        sampling_data = QoZ::sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);
        if (sampling_num == conf.num) {
            conf.cmprAlgo = QoZ::ALGO_INTERP;
       
        }

        else{
            /*
            QoZ::Config lorenzo_config = conf;
            lorenzo_config.cmprAlgo = QoZ::ALGO_LORENZO_REG;
            lorenzo_config.setDims(sample_dims.begin(), sample_dims.end());
            lorenzo_config.lorenzo = true;
            lorenzo_config.lorenzo2 = true;
            lorenzo_config.regression = false;
            lorenzo_config.regression2 = false;
            lorenzo_config.openmp = false;
            lorenzo_config.blockSize = 5;//why?
            size_t sampleOutSize;
            std::vector<T> cur_sampling_data=sampling_data;
            auto cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, cur_sampling_data.data(), sampleOutSize);
                
            delete[]cmprData;
            double ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
            if(conf.verbose)
                printf("Lorenzo ratio = %.4f\n", ratio);

            best_lorenzo_ratio = ratio;
            */
            double best_interp_ratio = 0;


            for (auto &interp_op: {QoZ::INTERP_ALGO_LINEAR, QoZ::INTERP_ALGO_CUBIC}) {
                double ratio = do_not_use_this_interp_compress_block_test<T, N>(sampling_data.data(), sample_dims, sampling_num, conf.absErrorBound,
                                                                             interp_op, conf.interpMeta.interpDirection, sampling_block);
                if (ratio > best_interp_ratio) {
                    best_interp_ratio = ratio;
                    conf.interpMeta.interpAlgo = interp_op;
                }
            }
            if(conf.verbose)
                std::cout << "interp best interpAlgo = " << (conf.interpMeta.interpAlgo == 0 ? "LINEAR" : "CUBIC") << std::endl;
                
            int direction_op = QoZ::factorial(N) - 1;
            double ratio = do_not_use_this_interp_compress_block_test<T, N>(sampling_data.data(), sample_dims, sampling_num, conf.absErrorBound,
                                                                         conf.interpMeta.interpAlgo, direction_op, sampling_block);
            if (ratio > best_interp_ratio * 1.02) {
                best_interp_ratio = ratio;
                conf.interpMeta.interpDirection = direction_op;
            }
            useInterp=!(best_lorenzo_ratio > best_interp_ratio && best_lorenzo_ratio < 80 && best_interp_ratio < 80);
            if(conf.verbose){
                std::cout << "interp best direction = " << (unsigned) conf.interpMeta.interpDirection << std::endl;
                
                printf("Interp ratio = %.4f\n", best_interp_ratio);
                    
                printf("choose %s\n", useInterp ? "interp" : "Lorenzo");
            }
            if (useInterp){
                conf.cmprAlgo=QoZ::ALGO_INTERP;
            }
            else{
                conf.cmprAlgo=QoZ::ALGO_LORENZO_REG;
            }
        }
        if(conf.verbose)
            timer.stop("sz3 tuning");
    }
    /*
    if(conf.verbose){
        timer.stop("PredTuning");
        timer.start();
    }
    */
    if (useInterp and conf.autoTuningRate>0){
            
        if(conf.verbose)
            std::cout<<"B-M tuning started."<<std::endl;
       
        if (conf.autoTuningRate!=conf.predictorTuningRate){//} and (conf.predictorTuningRate!=0 or conf.autoTuningRate!=conf.waveletTuningRate)){
              
            QoZ::sampleBlocks<T,N>(data,conf.dims,sampleBlockSize,sampled_blocks,conf.autoTuningRate,conf.profiling,starts,conf.var_first);
        }

        
        double bestalpha=1;
        double bestbeta=1;
        
        double bestb=9999;

        double bestm=0;
        size_t num_sampled_blocks=sampled_blocks.size();
        size_t per_block_ele_num=pow(sampleBlockSize+1,N);
        size_t ele_num=num_sampled_blocks*per_block_ele_num;
        std::vector<double> orig_means;//(num_sampled_blocks,0);
        std::vector<double> orig_sigma2s;//(num_sampled_blocks,0);
        std::vector<double> orig_ranges;//(num_sampled_blocks,0);
        conf.dims=std::vector<size_t>(N,sampleBlockSize+1);
        conf.num=per_block_ele_num;

        if(conf.tuningTarget==QoZ::TUNING_TARGET_SSIM){
            size_t ssim_size=conf.SSIMBlockSize;
            for (size_t k =0;k<num_sampled_blocks;k++){

                double orig_mean=0,orig_sigma2=0,orig_range=0;       
                if(N==2){
                    for (size_t i=0;i+ssim_size<sampleBlockSize+1;i+=ssim_size){
                        for (size_t j=0;j+ssim_size<sampleBlockSize+1;j+=ssim_size){
                            std::vector<size_t> starts{i,j};
                            QoZ::blockwise_profiling<T>(sampled_blocks[k].data(),conf.dims,starts,ssim_size,orig_mean,orig_sigma2,orig_range);
                            orig_means.push_back(orig_mean);
                            orig_sigma2s.push_back(orig_sigma2);
                            orig_ranges.push_back(orig_range);


                        }
                    }
                }
                else if(N==3){
                    for (size_t i=0;i+ssim_size<sampleBlockSize+1;i+=ssim_size){
                        for (size_t j=0;j+ssim_size<sampleBlockSize+1;j+=ssim_size){
                            for (size_t kk=0;kk+ssim_size<sampleBlockSize+1;kk+=ssim_size){
                                std::vector<size_t> starts{i,j,kk};
                                QoZ::blockwise_profiling<T>(sampled_blocks[k].data(),conf.dims,starts,ssim_size,orig_mean,orig_sigma2,orig_range);
                                orig_means.push_back(orig_mean);
                                orig_sigma2s.push_back(orig_sigma2);
                                orig_ranges.push_back(orig_range);
                            }
                        }
                    }
                }                      
            }

        }
        std::vector<T> flattened_sampled_data;
           
        if (conf.tuningTarget==QoZ::TUNING_TARGET_AC){
            for(int i=0;i<num_sampled_blocks;i++)
                flattened_sampled_data.insert(flattened_sampled_data.end(),sampled_blocks[i].begin(),sampled_blocks[i].end());

        }
        double oriabseb=conf.absErrorBound;
        

        
        
       
        //std::vector<double> flattened_cur_blocks;

        
        conf.dims=std::vector<size_t>(N,sampleBlockSize+1);
        conf.num=per_block_ele_num;

        std::vector<double>alpha_list;
        init_alphalist<T,N>(alpha_list,rel_bound,conf);
        size_t alpha_nums=alpha_list.size();
        std::vector<double>beta_list;
        init_betalist<T,N>(beta_list,rel_bound,conf);
        size_t beta_nums=beta_list.size();  
        
        for (size_t i=0;i<alpha_nums;i++){
            for (size_t j=0;j<beta_nums;j++){
                conf.absErrorBound=oriabseb;
                double alpha=alpha_list[i];
                double beta=beta_list[j];
                if (( (alpha>=1 and alpha>beta) or (alpha<0 and beta!=-1) ) )
                    continue;
                conf.alpha=alpha;
                conf.beta=beta; 
                auto sampled_blocks_copy=sampled_blocks;
                std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks_copy,QoZ::ALGO_INTERP,(QoZ::TUNING_TARGET)conf.tuningTarget,false,profiling_coeff,orig_means,
                                                                    orig_sigma2s,orig_ranges,flattened_sampled_data);
                double bitrate=results.first;
                double metric=results.second;
                if ( (conf.tuningTarget!=QoZ::TUNING_TARGET_CR and metric>=bestm and bitrate<=bestb) or (conf.tuningTarget==QoZ::TUNING_TARGET_CR and bitrate<=bestb ) ){
                    bestalpha=alpha;
                    bestbeta=beta;
                    bestb=bitrate;
                    bestm=metric;
                    useInterp=true;
                }
                else if ( (conf.tuningTarget!=QoZ::TUNING_TARGET_CR and metric<=bestm and bitrate>=bestb) or (conf.tuningTarget==QoZ::TUNING_TARGET_CR and bitrate>bestb) ){
                    if ( ((alpha>=1 and pow(alpha,max_interp_level-1)<=beta) or (alpha<1 and alpha*(max_interp_level-1)<=beta)))
                        break;

                    continue;
                }
                else{
                    double eb_fixrate;

                    eb_fixrate=bitrate/bestb;
                    double orieb=conf.absErrorBound;
                    conf.absErrorBound*=eb_fixrate;
                    auto sampled_blocks_copy=sampled_blocks;
                    std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks_copy,QoZ::ALGO_INTERP,(QoZ::TUNING_TARGET)conf.tuningTarget,false,profiling_coeff,orig_means,
                                                                        orig_sigma2s,orig_ranges,flattened_sampled_data);
                    conf.absErrorBound=orieb;

                    double bitrate_r=results.first;
                    double metric_r=results.second;
                    double a=(metric-metric_r)/(bitrate-bitrate_r);
                    double b=metric-a*bitrate;
                    double reg=a*bestb+b;

                    if (reg>bestm){
                        bestalpha=alpha;
                        bestbeta=beta;
                        bestb=bitrate;
                        bestm=metric;
                        useInterp=true;
                    }
                }
                if ( ( (alpha>=1 and pow(alpha,max_interp_level-1)<=beta) or (alpha<1 and alpha*(max_interp_level-1)<=beta)) )
                    break;

            }
        }

        //add lorenzo
        conf.absErrorBound=oriabseb;
        /*
        if(conf.testLorenzo){    


            std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks,QoZ::ALGO_LORENZO_REG,(QoZ::TUNING_TARGET)conf.tuningTarget,false,profiling_coeff,orig_means,
                    orig_sigma2s,orig_ranges,flattened_sampled_data);

            double bitrate=results.first;
            double metric=results.second;

            
            if ( (conf.tuningTarget!=QoZ::TUNING_TARGET_CR and metric>=bestm and bitrate<=bestb) or (conf.tuningTarget==QoZ::TUNING_TARGET_CR and bitrate<=bestb ) ){
                
                bestb=bitrate;
                bestm=metric;
                bestalpha=-1;
                bestbeta=-1;

                useInterp=false;
                   
            }
            else if ( (conf.tuningTarget!=QoZ::TUNING_TARGET_CR and metric<=bestm and bitrate>=bestb) or (conf.tuningTarget==QoZ::TUNING_TARGET_CR and bitrate>bestb) ){
                useInterp=true;
            }
            else{
                double eb_fixrate;

                eb_fixrate=bitrate/bestb;
                double orieb=conf.absErrorBound;
                conf.absErrorBound*=eb_fixrate;                        
                std::pair<double,double> results=CompressTest<T,N>(conf, sampled_blocks,QoZ::ALGO_LORENZO_REG,(QoZ::TUNING_TARGET)conf.tuningTarget,false,profiling_coeff,orig_means,
                                                                    orig_sigma2s,orig_ranges,flattened_sampled_data);
                conf.absErrorBound=orieb;
                double bitrate_r=results.first;
                double metric_r=results.second;
                double a=(metric-metric_r)/(bitrate-bitrate_r);
                double b=metric-a*bitrate;
                double reg=a*bestb+b;

                if (a>0 and reg>bestm){
                    bestb=bitrate;
                    bestm=metric;
                    bestalpha=-1;
                    bestbeta=-1;
                    useInterp=false;
                }
            }          
        }

        conf.absErrorBound=oriabseb;
        */
        /*
        if(conf.verbose){
            timer.stop("B-M step");
            timer.start();
        }
        */
        
        if(conf.tuningTarget==QoZ::TUNING_TARGET_AC){
            bestm=1-bestm;
        }
        std::string metric_name="Quality";
        if (conf.tuningTarget==QoZ::TUNING_TARGET_RD ){
            metric_name="PSNR";
        }
        else if (conf.tuningTarget==QoZ::TUNING_TARGET_SSIM ){
            metric_name="SSIM";
        }
        else if (conf.tuningTarget==QoZ::TUNING_TARGET_AC ){
            metric_name="AutoCorrelation";
        }
        if(conf.verbose){
            printf("Autotuning finished.\n");
            if (useInterp)
                printf("Interp selected. Selected alpha: %f. Selected beta: %f. Best bitrate: %f. Best %s: %f.\n",bestalpha,bestbeta,bestb, const_cast<char*>(metric_name.c_str()),bestm);
            else
                printf("Lorenzo selected. Best bitrate: %f. Best %s: %f.\n",bestb, const_cast<char*>(metric_name.c_str()),bestm);

        }
        conf.alpha=bestalpha;
        conf.beta=bestbeta;
        conf.dims=global_dims;
        conf.num=global_num;  

    }
    else if(useInterp and conf.QoZ){
        std::pair<double,double> ab=setABwithRelBound(rel_bound,2);
        conf.alpha=ab.first;
        conf.beta=ab.second;


    }
    conf.blockwiseTuning=blockwiseTuning;
    if (useInterp){
        conf.cmprAlgo=QoZ::ALGO_INTERP;
    }
    else{
         conf.cmprAlgo=QoZ::ALGO_LORENZO_REG;
    } 
        
    for(int i=0;i<sampled_blocks.size();i++){
        std::vector< T >().swap(sampled_blocks[i]);              
    }
    std::vector< std::vector<T> >().swap(sampled_blocks);  
    conf.quant_bins.clear();   
    return best_lorenzo_ratio;    
}



template<class T, QoZ::uint N>
std::array<char *,3> SZ_compress_Interp_lorenzo(std::array<QoZ::Config,3> &confs, std::array<T *,3>&data, std::array<size_t,3> &outSizes) {
    assert(confs[0].cmprAlgo == QoZ::ALGO_INTERP_LORENZO);
    for(auto i:{0,1,2})
       QoZ::calAbsErrorBound(confs[i], data[i]);
    for(auto i:{0,1,2}){
        if(N!=2&&N!=3){
            
                confs[i].autoTuningRate=0;
                confs[i].predictorTuningRate=0;
                confs[i].maxStep=0;
                confs[i].levelwisePredictionSelection=0;
                confs[i].multiDimInterp=0;
                confs[i].QoZ=0;
        }//merge this part to other branches
        

    
   
   
    
        if (confs[i].rng<0)
            confs[i].rng=QoZ::data_range<T>(data[i],confs[i].num);

        
        if (confs[i].relErrorBound<=0)
            confs[i].relErrorBound=confs[i].absErrorBound/confs[i].rng;
  
    
    }

    if(confs[0].verbose)
        std::cout << "====================================== BEGIN TUNING ================================" << std::endl;
    QoZ::Timer timer(true);
    if (confs[0].qoi>0)
        QoI_tuning<T,N>(confs,data);
    for(auto i:{0,1,2})
        Tuning<T,N>(confs[i],data[i]);
    
    //std::cout<<confs[0].qoi<<" "<<confs[1].qoi<<" "<<confs[2].qoi<<" "<<std::endl;


    if (1) {
    
        //std::vector<int>().swap(conf.quant_bins);
        double tuning_time = timer.stop();
        if(confs[0].verbose){
            std::cout << "Tuning time = " << tuning_time << "s" << std::endl;
            std::cout << "====================================== END TUNING ======================================" << std::endl;
        }

        //std::cout<<confs[0].qoi<<" "<<confs[1].qoi<<" "<<confs[2].qoi<<" "<<std::endl;
        return SZ_compress_Interp<T, N>(confs, data, outSizes);        

    } 

    else {
        /*
        QoZ::Config lorenzo_config = conf;
        size_t sampling_num, sampling_block;        
        std::vector<size_t> sample_dims(N);
        std::vector<T> sampling_data;

        size_t sampleOutSize;
        double ratio;  
            
        sampling_data = QoZ::sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);    
        lorenzo_config.cmprAlgo = QoZ::ALGO_LORENZO_REG;
        
        lorenzo_config.lorenzo = true;
        lorenzo_config.lorenzo2 = true;
        lorenzo_config.regression = false;
        lorenzo_config.regression2 = false;
        lorenzo_config.openmp = false;
        lorenzo_config.blockSize = 5;//why?
        if (sampling_num != conf.num) {
            lorenzo_config.setDims(sample_dims.begin(), sample_dims.end());
       
        //lorenzo_config.quantbinCnt = 65536 * 2;
                    
            if(conf.autoTuningRate>0 or conf.predictorTuningRate>0){
                auto cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), sampleOutSize);
                delete[]cmprData;
                ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
                //printf("Lorenzo ratio = %.2f\n", ratio);
                best_lorenzo_ratio = ratio;
            }
          
            //further tune lorenzo
            if (N == 3 ) {
                lorenzo_config.quantbinCnt = QoZ::optimize_quant_invl_3d<T>(data, conf.dims[0], conf.dims[1], conf.dims[2], conf.absErrorBound);
                lorenzo_config.pred_dim = 2;
                auto cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), sampleOutSize);
                delete[]cmprData;
                ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
                //printf("Lorenzo, pred_dim=2, ratio = %.4f\n", ratio);
                if (ratio > best_lorenzo_ratio * 1.02) {
                    best_lorenzo_ratio = ratio;
                } else {
                    lorenzo_config.pred_dim = 3;
                }
            }
            if (conf.relErrorBound < 1.01e-6 && best_lorenzo_ratio > 5 && lorenzo_config.quantbinCnt != 16384) {
                auto quant_num = lorenzo_config.quantbinCnt;
                lorenzo_config.quantbinCnt = 16384;
                auto cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), sampleOutSize);
                delete[]cmprData;
                ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
    //            printf("Lorenzo, quant_bin=8192, ratio = %.2f\n", ratio);
                if (ratio > best_lorenzo_ratio * 1.02) {
                    best_lorenzo_ratio = ratio;
                } else {
                    lorenzo_config.quantbinCnt = quant_num;
                }
            }
            lorenzo_config.setDims(conf.dims.begin(), conf.dims.end());
        }
     
        
        conf = lorenzo_config;
        double tuning_time = timer.stop();
        if(conf.verbose){
            std::cout << "Tuning time = " << tuning_time << "s" << std::endl;
            std::cout << "====================================== END TUNING ======================================" << std::endl;
        }
        return SZ_compress_LorenzoReg<T, N>(conf, data, outSize);
        */
        return std::array<char *,3>{};
    }
  
}


#endif