//
// Created by Xin Liang on 12/06/2021.
//

#ifndef SZ_QOI_X_EXP_HPP
#define SZ_QOI_X_EXP_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include "QoI.hpp"
#include "SymEngine.hpp"
#include <symengine/expression.h>
#include <symengine/parser.h>
#include <symengine/symbol.h>
#include <symengine/derivative.h>
#include <symengine/eval.h> 
#include <symengine/solve.h>
#include <symengine/functions.h>
#include <set>

using SymEngine::Expression;
using SymEngine::Symbol;
using SymEngine::symbol;
using SymEngine::parse;
using SymEngine::diff;
using SymEngine::RealDouble;
using SymEngine::Integer;
using SymEngine::evalf;
using SymEngine::map_basic_basic;
using SymEngine::down_cast;
using SymEngine::RCP;
using SymEngine::Basic;
using SymEngine::real_double;
using SymEngine::eval_double;


using SymEngine::solve;
using SymEngine::rcp_static_cast;
using SymEngine::Mul;
using SymEngine::Pow;
using SymEngine::Log;
using SymEngine::sqrt;
using SymEngine::is_a;
using SymEngine::FiniteSet;


namespace QoZ {
    template<class T>
    class QoI_X_Exp : public concepts::QoIInterface<T> {

    public:
        QoI_X_Exp(double tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T>::id = 5;
           // std::cout<<"init 1 "<< std::endl;
            
        }

        T interpret_eb(T data) const {
            
            auto p = pow(base,-data);
            double bound = log(tolerance*p + 1.0)/ log_base;
            if (p*tolerance < 1.0)
                bound = std::min(bound, -log(1-tolerance*p)/ log_base);
            T eb = bound;
            return std::min(eb, global_eb);
        }

        T interpret_eb(const T * data, size_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            //if(isolated and (data-thresold)*(dec_data-thresold)<0)//maybe can remove
            //    return false;
            
            return (fabs(pow(base,data) - pow(base,dec_data)) <= tolerance);
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
            
            return pow(base,val); 

        } 

        std::string get_expression() const{
            return std::to_string(base)+"^x";
        }

        void pre_compute(const T * data){}

        void set_qoi_tolerance(double tol) {tolerance = tol;}
        double get_qoi_tolerance() {return tolerance;}
        
    private:



        double tolerance;
        T global_eb;
        double base = std::exp(1.0);
        double log_base = 1.0;
     
    };
}
#endif 
