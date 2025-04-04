
#ifndef SZ_QOI_X_SQUARE_HPP
#define SZ_QOI_X_SQUARE_HPP

#include <algorithm>
#include <cmath>
#include "QoI.hpp"

namespace QoZ {
    template<class T>
    class QoI_X_Square : public concepts::QoIInterface<T> {

    public:
        QoI_X_Square(double tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T>::id = 1;
           // std::cout<<"init 1 "<< std::endl;
            
        }

        T interpret_eb(T data) const {
            

            T eb = - fabs(data) + sqrt(data * data + tolerance);
            return std::min(eb, global_eb);
        }

        T interpret_eb(const T * data, size_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            //if(isolated and (data-thresold)*(dec_data-thresold)<0)//maybe can remove
            //    return false;
            
            return (fabs(data*data- dec_data*dec_data) <= tolerance);
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
            
            return val*val; 

        } 

        std::string get_expression(const std::string var="x") const{
             return var+"^2";
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
