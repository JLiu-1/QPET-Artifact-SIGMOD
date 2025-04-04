

#ifndef SZ_QOI_X_CUBIC_HPP
#define SZ_QOI_X_CUBIC_HPP

#include <algorithm>
#include <cmath>
#include "QoI.hpp"


namespace QoZ {
    template<class T>
    class QoI_X_Cubic : public concepts::QoIInterface<T> {

    public:
        QoI_X_Cubic(double tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T>::id = 9;
           // std::cout<<"init 1 "<< std::endl;
            
        }

        T interpret_eb(T data) const {
            

            T eb = data >= 0 ? std::cbrt(data*data*data+tolerance)-data : data - std::cbrt(data*data*data-tolerance);
            return std::min(eb, global_eb);
        }

        T interpret_eb(const T * data, size_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            //if(isolated and (data-thresold)*(dec_data-thresold)<0)//maybe can remove
            //    return false;
            
            return (fabs(data*data*data - dec_data*dec_data*dec_data) <= tolerance);
        }

        void update_tolerance(T data, T dec_data){}

        void precompress_block(){}

        void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

        double eval(T val) const{
            
            return val*val*val; 

        } 

        std::string get_expression(const std::string var="x") const{
            return var+"^3";
        }

        void pre_compute(const T * data){}

        void set_qoi_tolerance(double tol) {tolerance = tol;}
        double get_qoi_tolerance() {return tolerance;}
        
    private:



        double tolerance;
        T global_eb;
     
    };

    template<class T>
    class QoI_X_Cubic_Approx : public concepts::QoIInterface<T> {

    public:
        QoI_X_Cubic_Approx(double tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T>::id = 9;
           // std::cout<<"init 1 "<< std::endl;
            
        }

        T interpret_eb(T data) const {
            

            data = fabs(data);
            double b = 6*data;
            double a = 0.5*b*data;//datatype may be T

            T eb;
            if(!std::isnan(a) and !std::isnan(b) and !std::isinf(a) and !std::isinf(b) and b >= 1e-10 )
                eb = (sqrt(a*a+2*b*tolerance)-a)/b;
            else if (!std::isnan(a) and !std::isinf(a) and a!=0 )
                eb = tolerance/a;
            else 
                eb = global_eb;
          
            //T eb = data >= 0 ? std::cbrt(data*data*data+tolerance)-data : data - std::cbrt(data*data*data-tolerance);
            return std::min(eb, global_eb);
        }

        T interpret_eb(const T * data, size_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            //if(isolated and (data-thresold)*(dec_data-thresold)<0)//maybe can remove
            //    return false;
            
            return (fabs(data*data*data - dec_data*dec_data*dec_data) <= tolerance);
        }

        void update_tolerance(T data, T dec_data){}

        void precompress_block(){}

        void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

        double eval(T val) const{
            
            return val*val*val; 

        } 

        std::string get_expression(const std::string var="x") const{
            return var+"^3";
        }

        void pre_compute(const T * data){}

        void set_qoi_tolerance(double tol) {tolerance = tol;}
        double get_qoi_tolerance() {return tolerance;}
        
    private:



        double tolerance;
        T global_eb;
     
    };

}
#endif 
