//
// Created by Xin Liang on 12/06/2021.
//

#ifndef SZ_QOI_FX_HPP
#define SZ_QOI_FX_HPP

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
using SymEngine::Sqrt;
using SymEngine::is_a;

namespace QoZ {
    template<class T, uint N>
    class QoI_FX : public concepts::QoIInterface<T, N> {

    public:
        QoI_FX(T tolerance, T global_eb, std::string ff = "x^2") : 
                tolerance(tolerance),
                global_eb(global_eb) {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 14;
           // std::cout<<"init 1 "<< std::endl;
            
            Expression f;
            Expression df;
            Expression ddf;
             x = symbol("x");
    
            f = Expression(ff).simplify();;//may also expand or remove simplify

            std::vector<Expression>singularities = find_singularities(f);
            for (const auto& singularity : singularities) {
                std::cout << singularity << std::endl;
            }
            // std::cout<<"init 2"<< std::endl;
            //df = diff(f,x);
            df = f.diff(x).simplify();
            // std::cout<<"init 3 "<< std::endl;
            //ddf = diff(df,x);
            ddf = df.diff(x).simplify();
            std::cout<<"f: "<< f<<std::endl;
            std::cout<<"df: "<< df<<std::endl;
            std::cout<<"ddf: "<< ddf<<std::endl;
  
            func = convert_expression_to_function(f, x);
            deri_1 = convert_expression_to_function(df, x);
            deri_2 = convert_expression_to_function(ddf, x);
            // std::cout<<"init 4 "<< std::endl;
              
           // RCP<const Basic> result = evalf(df.subs(map_basic_basic({{x,RealDouble(2).rcp_from_this()}})),53, SymEngine::EvalfDomain::Real);
           // RCP<const Symbol> value = symbol("2");
           // map_basic_basic mbb=  {{x,value}};
            //std::cout<<"init 5 "<< std::endl;
             //double result = (double)df.subs({{x,real_double(2)}}); 
           
           // std::cout<<"Eval res: "<<result<<std::endl;
            //SymEngine::RCP<const Basic> result = evalf(df,53, SymEngine::EvalfDomain::Real);
            //std::cout<< (down_cast<const RealDouble &>(*result)).as_double()<<std::endl;
        }

        using Range = multi_dimensional_range<T, N>;
        using iterator = typename multi_dimensional_range<T, N>::iterator;

        T interpret_eb(T data) const {
            

            double a = fabs(deri_1(data));//datatype may be T
            double b = fabs(deri_2(data));
           // 
            T eb;
            if(!std::isnan(a) and !std::isnan(b) and b !=0 )
                eb = (sqrt(a*a+2*b*tolerance)-a)/b;
            else if (!std::isnan(a) and a!=0 )
                eb = tolerance/a;
            else 
                eb = global_eb;
           // std::cout<<data<<" "<<a<<" "<<b<<" "<<eb<<" "<<global_eb<<std::endl; 
            return std::min(eb, global_eb);
        }

        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }

        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            return (fabs(func(data) - func(dec_data)) < tolerance);
        }

        void update_tolerance(T data, T dec_data){}

        void precompress_block(const std::shared_ptr<Range> &range){}

        void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

    private:

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

        std::vector<Expression> find_singularities(const Expression& expr, const Expression& x) {
            std::vector<Expression> singularities;

            if (is_a<Mul>(*expr)) {
                auto mul_expr = rcp_static_cast<const Mul>(expr.get_basic());
                for (auto arg : mul_expr->get_args()) {
                    if (is_a<Pow>(*arg)) {
                        auto pow_expr = rcp_static_cast<const Pow>(arg);
                        if (pow_expr->get_exp()->__eq__(*SymEngine::minus_one)) {
                            auto denominator = pow_expr->get_base();
                            auto solutions = solve(denominator, x);
                            for (auto sol : solutions) {
                                singularities.push_back(sol);
                            }
                        }
                    } else {
                        auto sub_singularities = find_singularities(Expression(arg), x);
                        singularities.insert(singularities.end(), sub_singularities.begin(), sub_singularities.end());
                    }
                }
            }
            
            if (is_a<Log>(*expr)) {
                auto log_expr = rcp_static_cast<const Log>(expr.get_basic());
                auto log_argument = log_expr->get_args()[0];
                auto log_solutions = solve(log_argument, x);
                for (auto sol : log_solutions) {
                    singularities.push_back(sol);
                }
            }

            if (is_a<Sqrt>(*expr)) {
                auto sqrt_expr = rcp_static_cast<const Sqrt>(expr.get_basic());
                auto sqrt_argument = sqrt_expr->get_args()[0];
                auto sqrt_solutions = solve(sqrt_argument, x);
                for (auto sol : sqrt_solutions) {
                    singularities.push_back(sol);
                }
            }

            if (is_a<Pow>(*expr)) {
                auto pow_expr = rcp_static_cast<const Pow>(expr.get_basic());
                auto exponent = pow_expr->get_exp();

                if (is_a<SymEngine::Number>(*exponent)) {
                    double exp_val = evalf(exponent).as_double();
                    if (exp_val > 0 && exp_val < 2 && std::floor(exp_val) != exp_val) {
                        auto base = pow_expr->get_base();
                        auto base_solutions = solve(base, x);
                        for (auto sol : base_solutions) {
                            singularities.push_back(sol);
                        }
                    }
                }
            }

            if (is_a<SymEngine::Add>(*expr)) {
                auto add_expr = rcp_static_cast<const SymEngine::Add>(expr.get_basic());
                for (auto arg : add_expr->get_args()) {
                    auto sub_singularities = find_singularities(Expression(arg), x);
                    singularities.insert(singularities.end(), sub_singularities.begin(), sub_singularities.end());
                }
            }

            if (is_a<SymEngine::Mul>(*expr) && !is_a<Pow>(*expr)) {
                auto mul_expr = rcp_static_cast<const Mul>(expr.get_basic());
                for (auto arg : mul_expr->get_args()) {
                    auto sub_singularities = find_singularities(Expression(arg), x);
                    singularities.insert(singularities.end(), sub_singularities.begin(), sub_singularities.end());
                }
            }

            return singularities;
        }



        RCP<const Symbol>  x;
        T tolerance;
        T global_eb;
        std::function<double(T)> func;
        std::function<double(T)> deri_1;
        std::function<double(T)> deri_2;
     
    };
}
#endif 
