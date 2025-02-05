//
// Created by Xin Liang on 12/06/2021.
//

#ifndef SZ_QOI_X_ABS_HPP
#define SZ_QOI_X_ABS_HPP

#include <algorithm>
#include <cmath>
#include "QoI.hpp"

namespace QoZ {
    template<class T>
    class QoI_X_Abs: public concepts::QoIInterface<T> {

    public:
        QoI_X_Abs(double tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb)
                 {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T>::id = 19;

        }


        T interpret_eb(T data) const {
            T eb = tolerance;
            return std::min(global_eb,eb);
        }



        T interpret_eb(const T * data, size_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return (fabs(std::abs(data) - std::abs(dec_data)) <= tolerance);
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
            
            return std::abs(val);

        } 

        std::string get_expression(const std::string var="x") const{
            return "|"+var+"|";
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
