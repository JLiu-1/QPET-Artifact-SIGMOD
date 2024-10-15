//
// Created by Xin Liang on 12/06/2021.
//

#ifndef SZ_QOI_FXYZ_HPP
#define SZ_QOI_FXYZ_HPP

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
    template<class T, uint N>
    class QoI_FXYZ : public concepts::QoIInterface<T, N> {

    public:
        QoI_FXYZ(const std::array<Config,3> confs)  {
            // TODO: adjust type for int data
            //printf("global_eb = %.4f\n", (double) global_eb);
            concepts::QoIInterface<T, N>::id = 1;


            tolerance = confs[0].qoiEB;
            func_string = confs[0].qoi_string;
            confidence = confs[0].confidence;
            global_eb = confs[0].absErrorBound;
            for(auto i:{0,1,2})
                global_ebs[i]=confs[i].absErrorBound;

           // std::cout<<"init 1 "<< std::endl;
            
            Expression f;
            Expression dfdx,dfdy,dfdz;
            x = symbol("x");
            y = symbol("y");
            z = symbol("z");
    
            f = Expression(func_string);
            /*
            singularities = find_singularities(f,x);
            std::cout << "Singularities:" << std::endl;
            for (const auto& singularity : singularities) {
                std::cout << singularity << std::endl;
            }
            */
            // std::cout<<"init 2"<< std::endl;
            //df = diff(f,x);
            dfdx = f.diff(x);
            // std::cout<<"init 3 "<< std::endl;
            //ddf = diff(df,x);
            dfdy = f.diff(y);
            dfdz = f.diff(z);
            std::cout<<"f: "<< f<<std::endl;
            std::cout<<"dfdx: "<< dfdx<<std::endl;
            std::cout<<"dfdy: "<< dfdy<<std::endl;
            std::cout<<"dfdz: "<< dfdz<<std::endl;
  
            func = convert_expression_to_function(f, x,y,z);
            dx = convert_expression_to_function(dfdx, x,y,z);
            dy = convert_expression_to_function(dfdy, x,y,z);
            dz = convert_expression_to_function(dfdz, x,y,z);

            std::cout<<dx(1,1,1)<<" "<<dy(1,1,1)<<" "<<dz(1,1,1)<<std::endl;
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

        //using Range = multi_dimensional_range<T, N>;
        //using iterator = typename multi_dimensional_range<T, N>::iterator;

        std::array<T,3> interpret_eb(T x, T y, T z) const {
            

            double alpha = fabs(dx(x,y,z));//datatype may be T
            double beta = fabs(dy(x,y,z));
            double gamma = fabs(dz(x,y,z));
            double sum= alpha+beta+gamma;
            double square_sum= alpha*alpha+beta*beta+gamma*gamma;


           // 
            double k = 1.732;
            T estimation_1 = square_sum!=0?k*sqrt(0.5/(square_sum*log(2.0/(1-confidence))))*tolerance:0;
            T estimation_2 = (sum!=0)?tolerance/sum:0;

            //std::cout<<
            
            T eb = std::max(estimation_1,estimation_2);
            if (eb == 0)
                eb=global_ebs[0]+global_ebs[1]+global_ebs[2];
            std::array<T,3> res;
            for (auto i:{0,1,2})
                res[i]=std::min(eb,global_ebs[i]);
        /*
             for (auto sg : singularities){
                T diff = fabs(data-sg);
                eb = std::min(diff,eb);
             }
             */
           // std::cout<<data<<" "<<a<<" "<<b<<" "<<eb<<" "<<global_eb<<std::endl; 
            return res;
        }
        /*
        T interpret_eb(const iterator &iter) const {
            return interpret_eb(*iter);
        }

        T interpret_eb(const T * data, ptrdiff_t offset) {
            return interpret_eb(*data);
        }
        */
        /*
        bool check_compliance(T data, T dec_data, bool verbose=false) const {
            //if(isolated and (data-thresold)*(dec_data-thresold)<0)//maybe can remove
            //    return false;
            return (fabs(func(data) - func(dec_data)) < tolerance);
        }
        */
        //void update_tolerance(T data, T dec_data){}

        //void precompress_block(const std::shared_ptr<Range> &range){}

        //void postcompress_block(){}

        void print(){}

        T get_global_eb() const { return global_eb; }

        void set_global_eb(T eb) {global_eb = eb;}

        void init(){}

        void set_dims(const std::vector<size_t>& new_dims){}

        double eval(T x, T y, T z) const{
            
            return func(x,y,z); 

        } 

        std::string get_expression() const{
            return func_string;
        }

        //void pre_compute(const T * data){}
        
    private:

        
        std::function<double(T, T, T)> convert_expression_to_function(const Basic &expr, 
                                                              const RCP<const Symbol> &x, 
                                                              const RCP<const Symbol> &y, 
                                                              const RCP<const Symbol> &z) {
            
            if (is_a<const SymEngine::Symbol>(expr) && expr.__eq__(*x)) {
                return [](T x_value, T, T) { /*std::cout<<"x="<<x_value<<std::endl;*/return x_value; };
            }
          
            else if (is_a<const SymEngine::Symbol>(expr) && expr.__eq__(*y)) {
                return [](T, T y_value, T) { /*std::cout<<"y="<<y_value<<std::endl;*/return y_value; };
            }
           
            else if (is_a<const SymEngine::Symbol>(expr) && expr.__eq__(*z)) {
                return [](T, T, T z_value) { /*std::cout<<"z="<<z_value<<std::endl;*/return z_value; };
            }
          
            else if (is_a<const RealDouble>(expr) or SymEngine::is_a<const Integer>(expr)) {
                double constant_value = eval_double(expr);
                return [constant_value](T, T, T) { /*std::cout<<"c="<<constant_value<<std::endl;*/return constant_value; };
            }
           
            else if (is_a<SymEngine::Add>(expr)) {
                auto args = expr.get_args();
                std::cout<<"add "<<args.size()<<std::endl;
                std::vector<std::function<double(T, T, T)> > fs;
                for (size_t i = 0; i < args.size(); ++i) {

                    fs.push_back(convert_expression_to_function(Expression(args[i]), x, y, z));
                }

               // auto first = convert_expression_to_function(Expression(args[0]), x, y, z);

                return [fs](T x_value, T y_value, T z_value) {
                    double result = 0;
                    for (auto &fnc:fs) {
                        result += fnc(x_value, y_value, z_value);
                    }
                    return result;
                };
            }
            else if (is_a<SymEngine::Mul>(expr)) {
                auto args = expr.get_args();
                std::vector<std::function<double(T, T, T)> > fs(args.size());
                for (size_t i = 0; i < args.size(); ++i) 
                    fs.push_back(convert_expression_to_function(Expression(args[i]), x, y, z));

                return [ fs](T x_value, T y_value, T z_value) {
                    double result = 1.0;
                    for (auto &fnc:fs) {
                        result *= fnc(x_value, y_value, z_value);
                    }
                    return result;
                };
            }
           
            else if (is_a<SymEngine::Pow>(expr)) {
                auto args = expr.get_args();
                auto base = convert_expression_to_function(Expression(args[0]), x, y, z);
                auto exponent = convert_expression_to_function(Expression(args[1]), x, y, z);
                return [base, exponent](T x_value, T y_value, T z_value) {
                    //std::cout<<"pow"<<std::endl;
                    return std::pow(base(x_value, y_value, z_value), exponent(x_value, y_value, z_value));
                };
            }
          
            else if (is_a<SymEngine::Sin>(expr)) {
                auto arg = convert_expression_to_function(Expression(expr.get_args()[0]), x, y, z);
                return [arg](T x_value, T y_value, T z_value) {
                    return std::sin(arg(x_value, y_value, z_value));
                };
            }
            
            else if (is_a<SymEngine::Cos>(expr)) {
                auto arg = convert_expression_to_function(Expression(expr.get_args()[0]), x, y, z);
                return [arg](T x_value, T y_value, T z_value) {
                    return std::cos(arg(x_value, y_value, z_value));
                };
            }
            
            else if (is_a<SymEngine::Log>(expr)) {
                auto args = expr.get_args();
                auto arg = convert_expression_to_function(Expression(args[0]), x, y, z);
                if (args.size() == 2) { // base log
                    auto base = convert_expression_to_function(Expression(args[1]), x, y, z);
                    return [arg, base](T x_value, T y_value, T z_value) {
                        return std::log(arg(x_value, y_value, z_value)) / std::log(base(x_value, y_value, z_value));
                    };
                } else { // ln
                    return [arg](T x_value, T y_value, T z_value) {
                        return std::log(arg(x_value, y_value, z_value));
                    };
                }
            }

            throw std::runtime_error("Unsupported expression type");
        }
        /*
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
        }*/



        RCP<const Symbol>  x,y,z;
        T tolerance;
        std::array<T,3> global_ebs;
        T global_eb;
        double confidence = 0.999999;
        std::function<double(T,T,T)> func;
        std::function<double(T,T,T)> dx;
        std::function<double(T,T,T)> dy;
        std::function<double(T,T,T)> dz;
        
        //std::set<double>singularities;
        std::string func_string;

        //double threshold;
       // bool isolated;
     
    };
}
#endif 
