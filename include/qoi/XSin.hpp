

#ifndef SZ_QOI_X_SIN_HPP
#define SZ_QOI_X_SIN_HPP

#include <algorithm>
#include <cmath>
#include "QoI.hpp"

namespace QoZ {

    template<class T>
    class QoI_X_Sin : public concepts::QoIInterface<T> {

    public:
        QoI_X_Sin(double tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T>::id = 22;
        }


        T interpret_eb(T data) const {
            
            if (tolerance>=2)
                return global_eb;
            T eb = tolerance;
            return std::min(eb,global_eb);
        }

        T interpret_eb(const T * data, size_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return (fabs(sin(data) - sin(dec_data)) < tolerance);
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
            
            return sin(val);//todo

        } 
        std::string get_expression(const std::string var="x") const{
            return "sin("+var+")";
        }

        void pre_compute(const T * data){}

        void set_qoi_tolerance(double tol) {tolerance = tol;}
        double get_qoi_tolerance() {return tolerance;}



    private:
        double tolerance;
        T global_eb;
     
    };


    template<class T>
    class QoI_X_Sin_Approx : public concepts::QoIInterface<T> {

    public:
        QoI_X_Sin_Approx(double tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T>::id = 22;
        }


        T interpret_eb(T data) const {
            

            double a = fabs(cos(data));//datatype may be T
            double b = fabs(sin(data));
            T eb;
            if(!std::isnan(a) and !std::isnan(b) and !std::isinf(a) and !std::isinf(b)and b >= 1e-10 )
                eb = (sqrt(a*a+2*b*tolerance)-a)/b;
            else if (!std::isnan(a) and !std::isinf(a) and a!=0 )
                eb = tolerance/a;
            else 
                eb = global_eb;
           
            return std::min(eb, global_eb);
        }


        T interpret_eb(const T * data, size_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return (fabs(sin(data) - sin(dec_data)) < tolerance);
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
            
            return sin(val);//todo

        } 
        std::string get_expression(const std::string var="x") const{
            return "sin("+var+")";
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
