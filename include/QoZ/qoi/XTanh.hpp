
#ifndef SZ_QOI_X_TANH_HPP
#define SZ_QOI_X_TANH_HPP

#include <algorithm>
#include <cmath>
#include "QoZ/def.hpp"
#include "QoZ/qoi/QoI.hpp"
#include "QoZ/utils/Iterator.hpp"

namespace QoZ {

    template<class T, uint N>
    class QoI_X_Tanh : public concepts::QoIInterface<T, N> {

    public:
        QoI_X_Tanh(double tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 23;
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            
            double t = std::tanh(data);
            T eb;
            if(t>=0){
                //if(tolerance*tolerance<=1+t*t)
                    eb = (t-tolerance)>-1 ? data-std::atanh(t-tolerance):global_eb;
                //else
                //    eb = (t+tolerance)<1 ? std::atanh(t+tolerance)-data:global_eb;
            }
            else{
                //if(tolerance*tolerance>=1+t*t)
                //    eb = (t-tolerance)>-1 ? data-std::atanh(t-tolerance):global_eb;
                //else
                    eb = (t+tolerance)<1 ? std::atanh(t+tolerance)-data:global_eb;
            }

            //double low_bound = (t-tolerance)>-1? data-std::atanh(t-tolerance):global_eb;
            //double high_bound = (t+tolerance)<-1? std::atanh(t+tolerance)-data:global_eb;
            //T eb = std::min(low_bound,high_bound);
            return std::min(eb,global_eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return (fabs(std::tanh(data) - std::tanh(dec_data)) < tolerance);
        }

        void update_tolerance(T data, T dec_data){}

        void precompress_block(const std::shared_ptr<Range> &range){}

        void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

        double eval(T val) const{
            
            return std::tanh(val);//todo

        } 
        std::string get_expression(const std::string var="x") const{
            return "tanh("+var+")";
        }

        void pre_compute(const T * data){}

        void set_qoi_tolerance(double tol) {tolerance = tol;}


    private:
        double tolerance;
        T global_eb;
     
    };


    template<class T, uint N>
    class QoI_X_Tanh_Approx : public concepts::QoIInterface<T, N> {

    public:
        QoI_X_Tanh_Approx(double tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            //printf("tolerance = %.4f\n", (double) tolerance );
            concepts::QoIInterface<T, N>::id = 23;
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            
            double t = std::tanh(data);
            double a = 1-t*t;
            double b = fabs(2*a*t);
            T eb = 0;
            if(!std::isnan(a) and !std::isnan(b) and !std::isinf(a) and !std::isinf(b) and b >=1e-10)
                eb = (sqrt(a*a+2*b*tolerance)-a)/b;
            if(eb==0){
                if (!std::isnan(a) and !std::isinf(a) and a!=0 )
                    eb = tolerance/a;
                else 
                    eb = global_eb;
            }
           
            return std::min(eb, global_eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return (fabs(std::tanh(data) - std::tanh(dec_data)) < tolerance);
        }

        void update_tolerance(T data, T dec_data){}

        void precompress_block(const std::shared_ptr<Range> &range){}

        void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

        double eval(T val) const{
            
            return std::tanh(val);//todo

        } 
        std::string get_expression(const std::string var="x") const{
            return "Tanh("+var+")";
        }

        void pre_compute(const T * data){}

        void set_qoi_tolerance(double tol) {tolerance = tol;}


    private:
        double tolerance;
        T global_eb;
     
    };
}
#endif 
