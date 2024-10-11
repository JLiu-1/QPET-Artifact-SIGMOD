//
// Created by Xin Liang on 03/09/2021.
//

#ifndef SZ_QOI_R_FX_HPP
#define SZ_QOI_R_FX_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include "QoZ/def.hpp"
#include "QoZ/qoi/QoI.hpp"
#include "QoZ/utils/Iterator.hpp"
#include <symengine/expression.h>
#include <symengine/parser.h>
#include <symengine/symbol.h>
#include <symengine/derivative.h>
#include <symengine/eval.h> 
#include <symengine/solve.h>
#include <symengine/functions.h>
#include<set>

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
    template<class T, uint N>
    class QoI_RegionalFX_lorenzo : public concepts::QoIInterface<T, N> {

    public:
        QoI_RegionalFX_lorenzo(T tolerance, T global_eb) : 
                tolerance(tolerance),
                global_eb(global_eb) {
            printf("QoI_RegionalAverageOfSquare\n");
            printf("tolerance = %.4e\n", (double) tolerance);
            printf("global_eb = %.4e\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 16;
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            // eb for data^2
            double eb_x2 = (aggregated_tolerance - fabs(error)) / rest_elements;
            // compute eb based on formula of x^2
            T eb = - fabs(data) + sqrt(data * data + eb_x2);
            return std::min(eb, global_eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        void update_tolerance(T data, T dec_data){
            error += data*data - dec_data*dec_data;
            rest_elements --;
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return true;
        }

        void precompress_block(const std::shared_ptr<Range> &range){
            // compute number of elements
            auto dims = range->get_dimensions();
            size_t num_elements = 1;
            for (const auto &dim: dims) {
                num_elements *= dim;
            }
            // assignment
            rest_elements = num_elements;
            block_elements = num_elements;
            aggregated_tolerance = tolerance * num_elements;
        }

        void postcompress_block(){
            error = 0;
        }

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

        double eval(T val) const{
            
            return 0;//todo

        } 

        std::string get_expression() const{
            return "Regional average";
        }

        void pre_compute(const T * data){}

    private:
        T tolerance;
        T global_eb;
        double error = 0;
        int rest_elements;
        int block_elements;
        double aggregated_tolerance;
    };

//This is sum_{aif(xi)}
//now, ai is just 1/n (average)
    template<class T, uint N>
    class QoI_RegionalFX : public concepts::QoIInterface<T, N> {

    public:
        QoI_RegionalFX(T tolerance, T global_eb, int block_size, std::vector<size_t> dims, std::string ff = "x^2", bool isolated = false, double threshold = 0.0) : 
                tolerance(tolerance),
                global_eb(global_eb),
                dims(dims),
                block_size(block_size),
                func_string (ff) , isolated (isolated), threshold (threshold){
            printf("tolerance = %.4e\n", (double) tolerance);
            printf("global_eb = %.4e\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 16;

            Expression f;
            Expression df;
            Expression ddf;
            x = symbol("x");
    
            f = Expression(ff);

            singularities = find_singularities(f,x);
            std::cout << "Singularities:" << std::endl;
            for (const auto& singularity : singularities) {
                std::cout << singularity << std::endl;
            }
            // std::cout<<"init 2"<< std::endl;
            //df = diff(f,x);
            df = f.diff(x);
            // std::cout<<"init 3 "<< std::endl;
            //ddf = diff(df,x);
            ddf = df.diff(x);
            std::cout<<"f: "<< f<<std::endl;
            std::cout<<"df: "<< df<<std::endl;
            std::cout<<"ddf: "<< ddf<<std::endl;
  
            func = convert_expression_to_function(f, x);
            deri_1 = convert_expression_to_function(df, x);
            deri_2 = convert_expression_to_function(ddf, x);


            init();
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            std::cerr << "Not implemented\n";
            exit(-1);
            return 0;
        }

        T interpret_eb(const iterator &iter) const {
            std::cerr << "Not implemented\n";
            exit(-1);
            return 0;
        }

        T interpret_eb(const T *data, ptrdiff_t offset) {
            block_id = compute_block_id(offset);
            double Li = L_i[offset];
            double ai = 1.0 / block_elements[block_id];
            if(ai==0)
                return global_eb;
            double T_estimation_1,T_estimation_2;
            if(Li==0){
                T_estimation_1 = sum_aiti_square_tolerance_sqrt*sqrt(block_elements[block_id]); //all even T
                T_estimation_2 = tolerance; // tolerance/(sum{a_i})
            }
            else{
                T_estimation_1 = (1/(ai*ai*Li))*sum_aiti_square_tolerance_sqrt/block_sum_aiLi_square_reciprocal_sqrt[block_id];//todo: nan issue
                T_estimation_2 = tolerance*Li/block_sum_aiLi[block_id];//todo: nan issue
            }

            //if(offset%64==0 and block_id%10000==0){
            //    std::cout<<block_id<<" "<<*data<<" "<<Li<<" "<<block_sum_aiLi[block_id]<<" "<<block_sum_aiLi_square_reciprocal_sqrt[block_id]<<" "<<T_estimation_1<<" "<<T_estimation_2<<std::endl;
           // }

            double T_estimation = std::max(T_estimation_1,T_estimation_2);

            double a = Li;//datatype may be T
            double b = fabs(deri_2(*data));
           // 
            T eb;
            if(!std::isnan(a) and !std::isnan(b) and b !=0 )
                eb = (sqrt(a*a+2*b*T_estimation)-a)/b;
            else if (!std::isnan(a) and a!=0 )
                eb = T_estimation/a;
            else 
                eb = global_eb;



            //T eb = L !=0 ? (aggregated_tolerance[block_id] - fabs(accumulated_error[block_id])) / (rest_elements[block_id] * L) : global_eb;


            for (auto sg : singularities){
                T diff = fabs(*data-sg);
                eb = std::min(diff,eb);
             }

            return std::min(eb, global_eb);


        }

        void update_tolerance(T data, T dec_data){
            //if (tuning)
            //    return;
            return;
           /*
            accumulated_error[block_id] += func(data) - func(dec_data);
            rest_elements[block_id] --;
            if(rest_elements[block_id] == 0 and fabs(accumulated_error[block_id]) > aggregated_tolerance[block_id]){
                printf("%d: %.4e / %.4e\n", block_id, accumulated_error[block_id], aggregated_tolerance[block_id]);
                printf("%d / %d\n", rest_elements[block_id], block_elements[block_id]);
                exit(-1);
            }
            */
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return true;//todo: a real check, and rework on a block if it is false 
           
        }

        void precompress_block(const std::shared_ptr<Range> &range){}

        void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){
            block_dims = std::vector<size_t>(dims.size());
            size_t num_blocks = 1;
            num_elements = 1;
            for(int i=0; i<dims.size(); i++){
                num_elements *= dims[i];
                block_dims[i] = (dims[i] - 1) / block_size + 1;
                num_blocks *= block_dims[i];
                std::cout << block_dims[i] << " ";
            }
            std::cout << std::endl;
            //aggregated_tolerance = std::vector<double>(num_blocks);
            block_elements = std::vector<size_t >(num_blocks, 0);
            //rest_elements = std::vector<size_t>(num_blocks, 0);
            L_i = std::vector<double>(num_elements,0);
            block_sum_aiLi=std::vector<double>(num_blocks, 0);
            block_sum_aiLi_square_reciprocal_sqrt=std::vector<double>(num_blocks, 0);


            if(dims.size() == 2){
                for(int i=0; i<block_dims[0]; i++){
                    int size_x = (i < block_dims[0] - 1) ? block_size : dims[0] - i * block_size;
                    for(int j=0; j<block_dims[1]; j++){
                        int size_y = (j < block_dims[1] - 1) ? block_size : dims[1] - j * block_size;
                        size_t num_block_elements = size_x * size_y;
                        //aggregated_tolerance[i * block_dims[1] + j] = num_block_elements * tolerance;
                        block_elements[i * block_dims[1] + j] = num_block_elements;
                        //rest_elements[i * block_dims[1] + j] = num_block_elements;
                        //double a_i = 1.0 / num_block_elements;

                    }
                }
            }
            else if(dims.size() == 3){
                for(int i=0; i<block_dims[0]; i++){
                    int size_x = (i < block_dims[0] - 1) ? block_size : dims[0] - i * block_size;
                    for(int j=0; j<block_dims[1]; j++){
                        int size_y = (j < block_dims[1] - 1) ? block_size : dims[1] - j * block_size;
                        for(int k=0; k<block_dims[2]; k++){
                            int size_z = (k < block_dims[2] - 1) ? block_size : dims[2] - k * block_size;
                            size_t num_block_elements = size_x * size_y * size_z;
                            // printf("%d, %d, %d: %d * %d * %d = %d\n", i, j, k, size_x, size_y, size_z, num_block_elements);
                            //aggregated_tolerance[i * block_dims[1] * block_dims[2] + j * block_dims[2] + k] = num_block_elements * tolerance;
                            block_elements[i * block_dims[1] * block_dims[2] + j * block_dims[2] + k] = num_block_elements;
                            //rest_elements[i * block_dims[1] * block_dims[2] + j * block_dims[2] + k] = num_block_elements;
                        }
                    }
                }
            }
            else{
                std::cerr << "dims other than 2 or 3 are not implemented" << std::endl;
                exit(-1);
            }

            //accumulated_error = std::vector<double>(num_blocks, 0);

            if (q>=0.95 and num_blocks >= 1000){
                sum_aiti_square_tolerance_sqrt = sqrt( -1 / ( 2 * ( log(1-q) - log(2*num_blocks) ) ) ) *tolerance;
            }
            else{
                sum_aiti_square_tolerance_sqrt = sqrt( -1 / ( 2 * ( log(1- pow(q,1.0/num_blocks) ) - log(2) ) ) ) *tolerance;

            }


           //std::cout<<sum_aiti_square_tolerance_sqrt<<std::endl;


    

            std::cout << "end of init\n";            
        }

        void set_dims(const std::vector<size_t>& new_dims){
            dims = new_dims;
        }

        double eval(T val) const{
            
            return func(val); //point_wise
        } 

        std::string get_expression() const{
            return "Regional average of " + func_string;
        }

        void pre_compute(const T * data){
            for(size_t i = 0; i < num_elements ; i ++){
                block_id = compute_block_id(i);
                L_i[i] = fabs(deri_1(data[i]));
                double ai = 1.0 / block_elements[block_id];
                double aiLi = ai*L_i[i]; 
                block_sum_aiLi[block_id]+=aiLi;
                if(aiLi!=0)
                    block_sum_aiLi_square_reciprocal_sqrt[block_id] += 1.0/(aiLi*aiLi);
            }
            for(auto &x:block_sum_aiLi_square_reciprocal_sqrt){
                x = sqrt(x);
            }
        }

    private:
        template<uint NN = N>
        inline typename std::enable_if<NN == 1, size_t>::type compute_block_id(ptrdiff_t offset) const noexcept {
            // 1D data
            return offset / block_size;
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 2, size_t>::type compute_block_id(ptrdiff_t offset) const noexcept {
            // 3D data
            int i = offset / dims[1];
            int j = offset % dims[1];
            return (i / block_size) * block_dims[1] + (j / block_size);
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 3, size_t>::type compute_block_id(ptrdiff_t offset) const noexcept {
            // 3D data
            int i = offset / (dims[1] * dims[2]);
            offset = offset % (dims[1] * dims[2]);
            int j = offset / dims[2];
            int k = offset % dims[2];
            return (i / block_size) * block_dims[1] * block_dims[2] + (j / block_size) * block_dims[2] + (k / block_size);
        }

        template<uint NN = N>
        inline typename std::enable_if<NN == 4, size_t>::type compute_block_id(ptrdiff_t offset) const noexcept {
            // 4D data
            std::cerr << "Not implemented!\n";
            exit(-1);
            return 0;
        }

        inline double evaluate(const Expression & func, T val) const{
            
            return (double)func.subs({{x,real_double(val)}}); 

        } 
        std::function<double(T)> convert_expression_to_function(const Basic &expr, const RCP<const Symbol> &x) {
            //std::cout<<SymEngine::type_code_name(expr.get_type_code())<<std::endl;
            // x
            if (is_a<const SymEngine::Symbol>(expr)) {
                return [](T x_value) { return x_value; };
            }
            // c
            else if (is_a<const RealDouble>(expr) or SymEngine::is_a<const Integer>(expr)) {
                double constant_value = eval_double(expr);
                return [constant_value](T) { return constant_value; };
            }
            // +
            else if (is_a<SymEngine::Add>(expr)) {
                auto args = expr.get_args();
                auto left = convert_expression_to_function(Expression(args[0]), x);
                auto right = convert_expression_to_function(Expression(args[1]), x);
                return [left, right](T x_value) {
                    return left(x_value) + right(x_value);
                };
            }
            // -
            /*
            else if (SymEngine::is_a<SymEngine::sub>(expr)) {
                auto args = expr.get_args();
                auto left = convert_expression_to_function(Expression(args[0]), x);
                auto right = convert_expression_to_function(Expression(args[1]), x);
                return [left, right](T x_value) {
                    return left(x_value) - right(x_value);
                };
            }*/
            // *
            else if (is_a<SymEngine::Mul>(expr)) {
                auto args = expr.get_args();
                auto left = convert_expression_to_function(Expression(args[0]), x);
                auto right = convert_expression_to_function(Expression(args[1]), x);
                return [left, right](T x_value) {
                    return left(x_value) * right(x_value);
                };
            }
            // /
            /*
            else if (SymEngine::is_a<SymEngine::div>(expr)) {
                auto args = expr.get_args();
                auto left = convert_expression_to_function(Expression(args[0]), x);
                auto right = convert_expression_to_function(Expression(args[1]), x);
                return [left, right](T x_value) {
                    return left(x_value) / right(x_value);
                };
            }*/
            // pow
            else if (is_a<SymEngine::Pow>(expr)) {
                auto args = expr.get_args();
                auto base = convert_expression_to_function(Expression(args[0]), x);
                auto exponent = convert_expression_to_function(Expression(args[1]), x);
                return [base, exponent](T x_value) {
                    return std::pow(base(x_value), exponent(x_value));
                };
            }
            // sin
            else if (is_a<SymEngine::Sin>(expr)) {
                auto arg = convert_expression_to_function(Expression(expr.get_args()[0]), x);
                return [arg](T x_value) {
                    return std::sin(arg(x_value));
                };
            }
            // cos
            else if (is_a<SymEngine::Cos>(expr)) {
                auto arg = convert_expression_to_function(Expression(expr.get_args()[0]), x);
                return [arg](T x_value) {
                    return std::cos(arg(x_value));
                };
            }

            else if (is_a<SymEngine::Tan>(expr)) {
                auto arg = convert_expression_to_function(Expression(expr.get_args()[0]), x);
                return [arg](T x_value) {
                    return std::tan(arg(x_value));
                };
            }

            else if (is_a<SymEngine::Sinh>(expr)) {
                auto arg = convert_expression_to_function(Expression(expr.get_args()[0]), x);
                return [arg](T x_value) {
                    return std::sinh(arg(x_value));
                };
            }
            // cos
            else if (is_a<SymEngine::Cosh>(expr)) {
                auto arg = convert_expression_to_function(Expression(expr.get_args()[0]), x);
                return [arg](T x_value) {
                    return std::cosh(arg(x_value));
                };
            }

            else if (is_a<SymEngine::Tanh>(expr)) {
                auto arg = convert_expression_to_function(Expression(expr.get_args()[0]), x);
                return [arg](T x_value) {
                    return std::tanh(arg(x_value));
                };
            }

            else if (is_a<SymEngine::Sign>(expr)) {
                auto arg = convert_expression_to_function(Expression(expr.get_args()[0]), x);
                return [arg](T x_value) {
                    return (x_value > 0) - (0 > x_value);
                };
            }
            //  log
            else if (is_a<SymEngine::Log>(expr)) {
                auto args = expr.get_args();
                auto arg = convert_expression_to_function(Expression(args[0]), x);

                if (args.size() == 2) { // base log
                    auto base = convert_expression_to_function(Expression(args[1]), x);
                    return [arg, base](T x_value) {
                        return std::log(arg(x_value)) / std::log(base(x_value));
                    };
                } else { // ln
                    return [arg](T x_value) {
                        return std::log(arg(x_value));
                    };
                }
            }

            throw std::runtime_error("Unsupported expression type");
        }

        std::set<double> find_singularities(const Expression& expr, const RCP<const Symbol> &x) {
            std::set<double> singularities;

            if (is_a<Mul>(expr)) {
                auto mul_expr = rcp_static_cast<const Mul>(expr.get_basic());
                for (auto arg : mul_expr->get_args()) {
                    if (is_a<Pow>(*arg)) {
                        auto pow_expr = rcp_static_cast<const Pow>(arg);
                        if (pow_expr->get_exp()->__eq__(*SymEngine::minus_one)) {
                            auto denominator = pow_expr->get_base();
                            auto solutions = solve(denominator, x);
                            if (is_a<FiniteSet>(*solutions)) {
                                auto finite_set_casted = rcp_static_cast<const FiniteSet>(solutions);
                                auto elements = finite_set_casted->get_container();
                                for (auto sol : elements) {

                                    if (is_a<const Integer>(*sol))
                                        singularities.insert(SymEngine::rcp_static_cast<const SymEngine::Integer>(sol)->as_int());
                                    else
                                        singularities.insert(SymEngine::rcp_static_cast<const SymEngine::RealDouble>(sol)->as_double());
                                }
                            }
                        }
                    } else {
                        auto sub_singularities = find_singularities(Expression(arg), x);
                        singularities.insert(sub_singularities.begin(), sub_singularities.end());
                    }
                }
            }
            
            if (is_a<Log>(expr)) {
                auto log_expr = rcp_static_cast<const Log>(expr.get_basic());
                auto log_argument = log_expr->get_args()[0];
                auto solutions = solve(log_argument, x);
                if (is_a<FiniteSet>(*solutions)) {
                    auto finite_set_casted = rcp_static_cast<const FiniteSet>(solutions);
                    auto elements = finite_set_casted->get_container();
                    for (auto sol : elements) {
                        if (is_a<const Integer>(*sol))
                            singularities.insert(SymEngine::rcp_static_cast<const SymEngine::Integer>(sol)->as_int());
                        else
                            singularities.insert(SymEngine::rcp_static_cast<const SymEngine::RealDouble>(sol)->as_double());
                    }
                }
            }

           
            if (is_a<Pow>(expr)) {
                auto pow_expr = rcp_static_cast<const Pow>(expr.get_basic());
                auto exponent = pow_expr->get_exp();
                //std::cout<<*exponent<<std::endl;

                if (is_a<const RealDouble>(*exponent) or is_a<const Integer>(*exponent)) {
                    double exp_val;
                    if (is_a<const Integer>(*exponent))
                        exp_val=SymEngine::rcp_static_cast<const SymEngine::Integer>(exponent)->as_int();
                    else
                        exp_val=SymEngine::rcp_static_cast<const SymEngine::RealDouble>(exponent)->as_double();
            
                    if (exp_val < 0 || (exp_val < 2 && std::floor(exp_val) != exp_val)) {
                        auto base = pow_expr->get_base();
                        auto solutions = solve(base, x);
                        if (is_a<FiniteSet>(*solutions)) {
                            auto finite_set_casted = rcp_static_cast<const FiniteSet>(solutions);
                            auto elements = finite_set_casted->get_container();
                            for (auto sol : elements) {
                                if (is_a<const Integer>(*sol))
                                    singularities.insert(SymEngine::rcp_static_cast<const SymEngine::Integer>(sol)->as_int());
                                else
                                    singularities.insert(SymEngine::rcp_static_cast<const SymEngine::RealDouble>(sol)->as_double());
                            }
                        }
                    }
                }
            }

            if (is_a<SymEngine::Add>(expr)) {
                auto add_expr = rcp_static_cast<const SymEngine::Add>(expr.get_basic());
                for (auto arg : add_expr->get_args()) {
                    auto sub_singularities = find_singularities(Expression(arg), x);
                    singularities.insert(sub_singularities.begin(), sub_singularities.end());
                }
            }

            if (is_a<SymEngine::Mul>(expr) && !is_a<Pow>(expr)) {
                auto mul_expr = rcp_static_cast<const Mul>(expr.get_basic());
                for (auto arg : mul_expr->get_args()) {
                    auto sub_singularities = find_singularities(Expression(arg), x);
                    singularities.insert(sub_singularities.begin(), sub_singularities.end());
                }
            }

            return singularities;
        }



        T tolerance;
        T global_eb;
        size_t block_id;
        size_t block_size;
        //std::vector<double> aggregated_tolerance;
        //std::vector<double> accumulated_error;
        std::vector<double> L_i;
        std::vector<double> block_sum_aiLi;
        std::vector<double> block_sum_aiLi_square_reciprocal_sqrt;
        std::vector<size_t> block_elements;
        //std::vector<int> rest_elements;
        std::vector<size_t> dims;
        std::vector<size_t> block_dims;
        size_t num_elements = 1;


        RCP<const Symbol>  x;
        std::function<double(T)> func;
        std::function<double(T)> deri_1;
        std::function<double(T)> deri_2;
        std::set<double>singularities;
        //bool tuning = false;

        std::string func_string;

        double threshold;
        bool isolated;

        double sum_aiti_square_tolerance_sqrt;
        double sum_aiti_tolerance;
        double q = 0.999999;

    };


}
#endif 