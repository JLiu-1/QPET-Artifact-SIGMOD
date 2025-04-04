

#ifndef SZ_QOI_X_RECIP_HPP
#define SZ_QOI_X_RECIP_HPP

#include <algorithm>

#include "QoI.hpp"

namespace QoZ {
    template<class T>
    class QoI_X_Recip : public concepts::QoIInterface<T> {

    public:
        QoI_X_Recip(double tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T>::id = 13;
        }

 
        T interpret_eb(T data) const {
            
            //if (data == 0)
            //    return global_eb;

            //double b = fabs(6*data);
            //double a = fabs(0.5*b*data);//datatype may be T

            
            //T eb = (sqrt(a*a+2*b*tolerance)-a)/b;
          
            T eb = data >= 0 ?  tolerance*data*data/(1.0+tolerance*data): tolerance*data*data/(1.0-tolerance*data);
            return std::min(eb, global_eb);
        }


        T interpret_eb(const T * data, size_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            if (data == 0) return (dec_data == 0);
            if (dec_data == 0) return false; 
            return (fabs(1.0/data - 1.0/dec_data) < tolerance);
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
            
            return 1.0/val;//todo

        } 

        std::string get_expression(const std::string var="x") const{
            return "1/"+var;
        }

        void pre_compute(const T * data){}

        void set_qoi_tolerance(double tol) {tolerance = tol;}
        double get_qoi_tolerance() {return tolerance;}


    private:
        double tolerance;
        T global_eb;
        
    };

    template<class T>
    class QoI_X_Recip_Approx : public concepts::QoIInterface<T> {

    public:
        QoI_X_Recip_Approx(double tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T>::id = 13;
        }

        T interpret_eb(T data) const {
            
            if (data == 0)
                return global_eb;
            double recip = fabs(1.0/data);
            double a = recip*recip;
            double b = 2*a*recip;
            //double a = -1/(data*data);//datatype may be T

            T eb;
            if(!std::isnan(a) and !std::isnan(b) and !std::isinf(a) and !std::isinf(b)and b >= 1e-10 )
                eb = (sqrt(a*a+2*b*tolerance)-a)/b;
            else if (!std::isnan(a) and !std::isinf(a) and a!=0 )
                eb = tolerance/a;
            else 
                eb = global_eb;



            
            //T eb = (sqrt(a*a+2*b*tolerance)-a)/b;
          
            //T eb = data >= 0 ?  tolerance*data*data/(1.0+tolerance*data): tolerance*data*data/(1.0-tolerance*data);
            return std::min(eb, global_eb);
        }

        T interpret_eb(const T * data, size_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            if (data == 0) return (dec_data == 0);
            if (dec_data == 0) return false; 
            return (fabs(1.0/data - 1.0/dec_data) < tolerance);
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
            
            return 1.0/val;//todo

        } 

        std::string get_expression(const std::string var="x") const{
            return "1/"+var;
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
