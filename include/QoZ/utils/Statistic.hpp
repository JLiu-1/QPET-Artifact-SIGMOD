
#ifndef SZ_STATISTIC_HPP
#define SZ_STATISTIC_HPP

#include "Config.hpp"

namespace QoZ {
    template<class T>
    T data_range(const T *data, size_t num) {

        T max = data[0];
        T min = data[0];
        for (size_t i = 1; i < num; i++) {
            if (max < data[i]) max = data[i];
            if (min > data[i]) min = data[i];
        }
        return max - min;
    }

    uint8_t factorial(int n) {
        return (n == 0) || (n == 1) ? 1 : n * factorial(n - 1);
    }

    double computeABSErrBoundFromPSNR(double psnr, double threshold, double value_range) {
        double v1 = psnr + 10 * log10(1 - 2.0 / 3.0 * threshold);
        double v2 = v1 / (-20);
        double v3 = pow(10, v2);
        return value_range * v3;
    }

    template<class T>
    void calAbsErrorBound(QoZ::Config &conf, const T *data,T range = 0 ) {
        if (conf.errorBoundMode != EB_ABS) {
            if (conf.errorBoundMode == EB_REL) {
                conf.errorBoundMode = EB_ABS;
                double rng= (range > 0) ? range : QoZ::data_range(data, conf.num);
                conf.rng=rng;
                conf.absErrorBound = conf.relErrorBound * rng;
            } else if (conf.errorBoundMode == EB_PSNR) {
                conf.errorBoundMode = EB_ABS;
                double rng=(range > 0) ? range : QoZ::data_range(data, conf.num);
                conf.rng=rng;
                conf.absErrorBound = computeABSErrBoundFromPSNR(conf.psnrErrorBound, 0.99, rng);
                conf.relErrorBound=conf.absErrorBound/rng;
            } else if (conf.errorBoundMode == EB_L2NORM) {
                conf.errorBoundMode = EB_ABS;
                conf.absErrorBound = sqrt(3.0 / conf.num) * conf.l2normErrorBound;
            } else if (conf.errorBoundMode == EB_ABS_AND_REL) {
                conf.errorBoundMode = EB_ABS;
                double rng=(range > 0) ? range : QoZ::data_range(data, conf.num);
                conf.rng=rng;
                conf.absErrorBound = std::min(conf.absErrorBound, conf.relErrorBound * rng);
            } else if (conf.errorBoundMode == EB_ABS_OR_REL) {
                conf.errorBoundMode = EB_ABS;
                double rng=(range > 0) ? range : QoZ::data_range(data, conf.num);
                conf.rng=rng;
                conf.absErrorBound = std::max(conf.absErrorBound, conf.relErrorBound *rng);
            } else {
                printf("Error, error bound mode not supported\n");
                exit(0);
            }
        }
    }

    template<typename Type>
    double autocorrelation1DLag1(const Type *data, size_t numOfElem, Type avg) {
        double cov = 0;
        for (size_t i = 0; i < numOfElem; i++) {
            cov += (data[i] - avg) * (data[i] - avg);
        }
        cov = cov / numOfElem;

        if (cov == 0) {
            return 0;
        } else {
            int delta = 1;
            double sum = 0;

            for (size_t i = 0; i < numOfElem - delta; i++) {
                sum += (data[i] - avg) * (data[i + delta] - avg);
            }
            return sum / (numOfElem - delta) / cov;
        }
    }

    template<typename Type>
    double calcNormedVariance(const Type *data, size_t num){
        double max=data[0],min=data[0];

        double average=0;

        for(size_t i=0;i<num;i++){
            double val=data[i];
            max=val>max?val:max;
            min=val<min?val:min;
            average+=val;
        }
        average/=num;
        average=(average-min)/(max-min);
        double variance=0;
        for(size_t i=0;i<num;i++){
            double val=data[i];
            val=(val-min)/(max-min);

            variance+=(val-average)*(val-average);

        }
        variance/=num;
        return variance;
        


    }

    template<typename Type>
    void verify(Type *ori_data, Type *data, size_t num_elements, double &psnr, double &nrmse) {
        size_t i = 0;
        double Max = ori_data[0];
        double Min = ori_data[0];
        double diffMax = fabs(data[0] - ori_data[0]);
        double diff_sum = 0;
        double maxpw_relerr = 0;
        double sum1 = 0, sum2 = 0, l2sum = 0;
        for (i = 0; i < num_elements; i++) {
            sum1 += ori_data[i];
            sum2 += data[i];
            l2sum += data[i] * data[i];
        }
        double mean1 = sum1 / num_elements;
        double mean2 = sum2 / num_elements;

        double sum3 = 0, sum4 = 0;
        double sum = 0, prodSum = 0, relerr = 0;

        double *diff = (double *) malloc(num_elements * sizeof(double));

        for (i = 0; i < num_elements; i++) {
            diff[i] = data[i] - ori_data[i];
            diff_sum += data[i] - ori_data[i];
            if (Max < ori_data[i]) Max = ori_data[i];
            if (Min > ori_data[i]) Min = ori_data[i];
            double err = fabs(data[i] - ori_data[i]);
            if (ori_data[i] != 0) {
                relerr = err / fabs(ori_data[i]);
                if (maxpw_relerr < relerr)
                    maxpw_relerr = relerr;
            }

            if (diffMax < err)
                diffMax = err;
            prodSum += (ori_data[i] - mean1) * (data[i] - mean2);
            sum3 += (ori_data[i] - mean1) * (ori_data[i] - mean1);
            sum4 += (data[i] - mean2) * (data[i] - mean2);
            sum += err * err;
        }
        double std1 = sqrt(sum3 / num_elements);
        double std2 = sqrt(sum4 / num_elements);
        double ee = prodSum / num_elements;
        double acEff = ee / std1 / std2;

        double mse = sum / num_elements;
        double range = Max - Min;
        psnr = 20 * log10(range) - 10 * log10(mse);
        nrmse = sqrt(mse) / range;

        double normErr = sqrt(sum);
        double normErr_norm = normErr / sqrt(l2sum);

        printf("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
        printf("Max absolute error = %.2G\n", diffMax);
        printf("Max relative error = %.2G\n", diffMax / (Max - Min));
        printf("Max pw relative error = %.2G\n", maxpw_relerr);
        printf("PSNR = %f, NRMSE= %.10G\n", psnr, nrmse);
        printf("normError = %f, normErr_norm = %f\n", normErr, normErr_norm);
        printf("acEff=%f\n", acEff);
//        printf("errAutoCorr=%.10f\n", autocorrelation1DLag1<double>(diff, num_elements, diff_sum / num_elements));
        free(diff);
    }

    template<typename Type>
    void verify(Type *ori_data, Type *data, size_t num_elements) {
        double psnr, nrmse;
        verify(ori_data, data, num_elements, psnr, nrmse);
    }


    template<typename Type>
    void verifyQoI(Type *ori_data, Type *data, std::vector<size_t> dims, int blockSize=1) {
        size_t num_elements = 1;
        for(const auto d:dims){
            num_elements *= d;
        }
        double psnr = 0;
        double nrmse = 0;
        verify(ori_data, data, num_elements, psnr, nrmse);
        Type max = ori_data[0];
        Type min = ori_data[0];
        for (size_t i = 1; i < num_elements; i++) {
            if (max < ori_data[i]) max = ori_data[i];
            if (min > ori_data[i]) min = ori_data[i];
        }

        double max_abs_val = std::max(fabs(max), fabs(min));
        double max_abs_val_sq = max_abs_val * max_abs_val;
        double max_abs_val_pow = pow(2,max);
        double max_abs_val_cu = max_abs_val * max_abs_val * max_abs_val;
        double max_x_square_diff = 0;
        double max_x_cubic_diff = 0;
        double max_log_diff = 0;
        double max_sin_diff = 0;
        double max_pow_diff = 0;
        double max_xlogx_diff = 0;
        double max_xlogx_val = ori_data[0]!=0 ? ori_data[0]*log(fabs(ori_data[0]))/log(2) : 0;
        double min_xlogx_val = max_xlogx_val;
        double max_tanh_diff = 0;
        double max_relu_diff = 0;

        for(int i=0; i<num_elements; i++){
            double x_square_diff = fabs(ori_data[i] * ori_data[i] - data[i] * data[i]);
            if(x_square_diff > max_x_square_diff) max_x_square_diff = x_square_diff;
            double x_cubic_diff = fabs(ori_data[i] * ori_data[i] * ori_data[i] - data[i] * data[i] * data[i]);
            if(x_cubic_diff > max_x_cubic_diff) max_x_cubic_diff = x_cubic_diff;
            double x_sin_diff = fabs(sin(ori_data[i]) - sin(data[i]));
            if(x_sin_diff > max_sin_diff) max_sin_diff = x_sin_diff;
            double x_pow_diff = fabs(pow(2,ori_data[i]) -pow(2,data[i]));
            if(x_pow_diff > max_pow_diff) max_pow_diff = x_pow_diff;
            // if(x_square_diff / max_abs_val_sq > 1e-5){
            //     std::cout << i << ": ori = " << ori_data[i] << ", dec = " << data[i] << ", err = " << x_square_diff / max_abs_val_sq << std::endl;
            // }
            double x_tanh_diff = fabs(std::tanh(ori_data[i])-std::tanh(data[i]));
            if (x_tanh_diff>max_tanh_diff) max_tanh_diff=x_tanh_diff;
            double x_relu_diff = fabs(ori_data[i]*(double)(ori_data[i]>0)-data[i]*(double)(data[i]>0));
            if (x_relu_diff>max_relu_diff) max_relu_diff=x_relu_diff;
            if(ori_data[i] == 0){
                // for log x only
                // if(data[i] != 0){
                //     std::cout << i << ": dec_data does not equal to 0\n";
                // }
            }
            else{
                if(ori_data[i] == data[i]) continue;
                double log2 = log(2);
                double log_diff = fabs(log(fabs(ori_data[i]))/log2 - log(fabs(data[i]))/log2);
                if(log_diff > max_log_diff) max_log_diff = log_diff;  
                // if(log_diff > 1){
                //     std::cout << i << ": " << ori_data[i] << " " << data[i] << std::endl;
                // }              
                double cur_xlogx = ori_data[i]!=0 ? ori_data[i]*log(fabs(ori_data[i]))/log(2) : 0;
                if (max_xlogx_val < cur_xlogx) max_xlogx_val = cur_xlogx;
                if (min_xlogx_val > cur_xlogx) min_xlogx_val = cur_xlogx;
                double dec_xlogx = data[i]!=0 ? data[i]*log(fabs(data[i]))/log(2) : 0;
                double xlogx_diff = fabs(cur_xlogx-dec_xlogx);
                if(xlogx_diff > max_xlogx_diff) max_xlogx_diff = xlogx_diff;  

            }


        }

        printf("QoI error info:\n");
        printf("Max x^2 error = %.6G, relative x^2 error = %.6G\n, Max x^3 error = %.6G, relative x^3 error = %.6G\n", max_x_square_diff, max_x_square_diff / max_abs_val_sq, max_x_cubic_diff, max_x_cubic_diff / max_abs_val_cu);
        printf("Max log error = %.6G, max sin error = %.6G, max 2^x error = %.6G, relative 2^x error = %.6G\n", max_log_diff, max_sin_diff, max_pow_diff, max_pow_diff / max_abs_val_pow);
        //std::cout<<max_xlogx_val<<" "<<min_xlogx_val<<std::endl;
        printf("Max xlogx error = %.6G, relative xlogx error = %.6G\n", max_xlogx_diff, max_xlogx_diff/(max_xlogx_val-min_xlogx_val));
        printf("Max tanh error = %.6G, max relu error = %.6G\n", max_tanh_diff, max_relu_diff);

        for(int i=0; i<num_elements; i++){
            ori_data[i] = ori_data[i] * ori_data[i];
            data[i] = data[i] * data[i];
        }
        max = ori_data[0];
        min = ori_data[0];
        for (size_t i = 1; i < num_elements; i++) {
            if (max < ori_data[i]) max = ori_data[i];
            if (min > ori_data[i]) min = ori_data[i];
        }
        if(dims.size() == 2) evaluate_average(ori_data, data, max - min, 1, dims[0], dims[1], blockSize);
        else if(dims.size() == 3) evaluate_average(ori_data, data, max - min, dims[0], dims[1], dims[2], blockSize);

    }


};


#endif
