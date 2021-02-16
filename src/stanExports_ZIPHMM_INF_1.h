// Generated by rstantools.  Do not edit by hand.

/*
    Bayeshmmcts is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Bayeshmmcts is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Bayeshmmcts.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_ZIPHMM_INF_1_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_ZIPHMM_INF_1");
    reader.add_event(52, 50, "end", "model_ZIPHMM_INF_1");
    return reader;
}
#include <stan_meta_header.hpp>
class model_ZIPHMM_INF_1
  : public stan::model::model_base_crtp<model_ZIPHMM_INF_1> {
private:
        int N;
        std::vector<int> y;
        int m;
public:
    model_ZIPHMM_INF_1(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_ZIPHMM_INF_1(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_ZIPHMM_INF_1_namespace::model_ZIPHMM_INF_1";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            check_greater_or_equal(function__, "N", N, 0);
            current_statement_begin__ = 6;
            validate_non_negative_index("y", "N", N);
            context__.validate_dims("data initialization", "y", "int", context__.to_vec(N));
            y = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("y");
            pos__ = 0;
            size_t y_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < y_k_0_max__; ++k_0__) {
                y[k_0__] = vals_i__[pos__++];
            }
            size_t y_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < y_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "y[i_0__]", y[i_0__], 0);
            }
            current_statement_begin__ = 7;
            context__.validate_dims("data initialization", "m", "int", context__.to_vec());
            m = int(0);
            vals_i__ = context__.vals_i("m");
            pos__ = 0;
            m = vals_i__[pos__++];
            check_greater_or_equal(function__, "m", m, 1);
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 11;
            num_params_r__ += 1;
            current_statement_begin__ = 12;
            validate_non_negative_index("lambda", "m", m);
            num_params_r__ += m;
            current_statement_begin__ = 13;
            validate_non_negative_index("A", "m", m);
            validate_non_negative_index("A", "m", m);
            num_params_r__ += ((m - 1) * m);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_ZIPHMM_INF_1() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 11;
        if (!(context__.contains_r("theta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable theta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("theta");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "theta", "double", context__.to_vec());
        double theta(0);
        theta = vals_r__[pos__++];
        try {
            writer__.scalar_lub_unconstrain(0, 1, theta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable theta: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 12;
        if (!(context__.contains_r("lambda")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable lambda missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("lambda");
        pos__ = 0U;
        validate_non_negative_index("lambda", "m", m);
        context__.validate_dims("parameter initialization", "lambda", "vector_d", context__.to_vec(m));
        Eigen::Matrix<double, Eigen::Dynamic, 1> lambda(m);
        size_t lambda_j_1_max__ = m;
        for (size_t j_1__ = 0; j_1__ < lambda_j_1_max__; ++j_1__) {
            lambda(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.positive_ordered_unconstrain(lambda);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable lambda: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 13;
        if (!(context__.contains_r("A")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable A missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("A");
        pos__ = 0U;
        validate_non_negative_index("A", "m", m);
        validate_non_negative_index("A", "m", m);
        context__.validate_dims("parameter initialization", "A", "vector_d", context__.to_vec(m,m));
        std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> > A(m, Eigen::Matrix<double, Eigen::Dynamic, 1>(m));
        size_t A_j_1_max__ = m;
        size_t A_k_0_max__ = m;
        for (size_t j_1__ = 0; j_1__ < A_j_1_max__; ++j_1__) {
            for (size_t k_0__ = 0; k_0__ < A_k_0_max__; ++k_0__) {
                A[k_0__](j_1__) = vals_r__[pos__++];
            }
        }
        size_t A_i_0_max__ = m;
        for (size_t i_0__ = 0; i_0__ < A_i_0_max__; ++i_0__) {
            try {
                writer__.simplex_unconstrain(A[i_0__]);
            } catch (const std::exception& e) {
                stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable A: ") + e.what()), current_statement_begin__, prog_reader__());
            }
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 11;
            local_scalar_t__ theta;
            (void) theta;  // dummy to suppress unused var warning
            if (jacobian__)
                theta = in__.scalar_lub_constrain(0, 1, lp__);
            else
                theta = in__.scalar_lub_constrain(0, 1);
            current_statement_begin__ = 12;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> lambda;
            (void) lambda;  // dummy to suppress unused var warning
            if (jacobian__)
                lambda = in__.positive_ordered_constrain(m, lp__);
            else
                lambda = in__.positive_ordered_constrain(m);
            current_statement_begin__ = 13;
            std::vector<Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> > A;
            size_t A_d_0_max__ = m;
            A.reserve(A_d_0_max__);
            for (size_t d_0__ = 0; d_0__ < A_d_0_max__; ++d_0__) {
                if (jacobian__)
                    A.push_back(in__.simplex_constrain(m, lp__));
                else
                    A.push_back(in__.simplex_constrain(m));
            }
            // model body
            {
            current_statement_begin__ = 17;
            validate_non_negative_index("log_A_tr", "m", m);
            validate_non_negative_index("log_A_tr", "m", m);
            std::vector<Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1>  > log_A_tr(m, Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1>(m));
            stan::math::initialize(log_A_tr, DUMMY_VAR__);
            stan::math::fill(log_A_tr, DUMMY_VAR__);
            current_statement_begin__ = 18;
            validate_non_negative_index("lp", "m", m);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> lp(m);
            stan::math::initialize(lp, DUMMY_VAR__);
            stan::math::fill(lp, DUMMY_VAR__);
            current_statement_begin__ = 19;
            validate_non_negative_index("lp_p1", "m", m);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> lp_p1(m);
            stan::math::initialize(lp_p1, DUMMY_VAR__);
            stan::math::fill(lp_p1, DUMMY_VAR__);
            current_statement_begin__ = 22;
            for (int i = 1; i <= m; ++i) {
                current_statement_begin__ = 23;
                for (int j = 1; j <= m; ++j) {
                    current_statement_begin__ = 24;
                    stan::model::assign(log_A_tr, 
                                stan::model::cons_list(stan::model::index_uni(j), stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list())), 
                                stan::math::log(get_base1(get_base1(A, i, "A", 1), j, "A", 2)), 
                                "assigning variable log_A_tr");
                }
            }
            current_statement_begin__ = 29;
            stan::math::assign(lp, rep_vector(stan::math::log(m), m));
            current_statement_begin__ = 31;
            for (int n = 1; n <= N; ++n) {
                current_statement_begin__ = 32;
                for (int j = 1; j <= m; ++j) {
                    current_statement_begin__ = 34;
                    stan::model::assign(lp_p1, 
                                stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list()), 
                                log_sum_exp(add(get_base1(log_A_tr, j, "log_A_tr", 1), lp)), 
                                "assigning variable lp_p1");
                    current_statement_begin__ = 37;
                    if (as_bool(logical_eq(j, 1))) {
                        current_statement_begin__ = 38;
                        if (as_bool(logical_eq(get_base1(y, n, "y", 1), 0))) {
                            current_statement_begin__ = 39;
                            stan::model::assign(lp_p1, 
                                        stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list()), 
                                        (stan::model::rvalue(lp_p1, stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list()), "lp_p1") + log_mix(theta, 0, poisson_log(0, get_base1(lambda, j, "lambda", 1)))), 
                                        "assigning variable lp_p1");
                        } else {
                            current_statement_begin__ = 41;
                            stan::model::assign(lp_p1, 
                                        stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list()), 
                                        (stan::model::rvalue(lp_p1, stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list()), "lp_p1") + (log1m(theta) + poisson_log(get_base1(y, n, "y", 1), get_base1(lambda, j, "lambda", 1)))), 
                                        "assigning variable lp_p1");
                        }
                    } else {
                        current_statement_begin__ = 44;
                        stan::model::assign(lp_p1, 
                                    stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list()), 
                                    (stan::model::rvalue(lp_p1, stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list()), "lp_p1") + poisson_log(get_base1(y, n, "y", 1), get_base1(lambda, j, "lambda", 1))), 
                                    "assigning variable lp_p1");
                    }
                }
                current_statement_begin__ = 47;
                stan::math::assign(lp, lp_p1);
            }
            current_statement_begin__ = 49;
            lp_accum__.add(log_sum_exp(lp));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("theta");
        names__.push_back("lambda");
        names__.push_back("A");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(m);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(m);
        dims__.push_back(m);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_ZIPHMM_INF_1_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double theta = in__.scalar_lub_constrain(0, 1);
        vars__.push_back(theta);
        Eigen::Matrix<double, Eigen::Dynamic, 1> lambda = in__.positive_ordered_constrain(m);
        size_t lambda_j_1_max__ = m;
        for (size_t j_1__ = 0; j_1__ < lambda_j_1_max__; ++j_1__) {
            vars__.push_back(lambda(j_1__));
        }
        std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> > A;
        size_t A_d_0_max__ = m;
        A.reserve(A_d_0_max__);
        for (size_t d_0__ = 0; d_0__ < A_d_0_max__; ++d_0__) {
            A.push_back(in__.simplex_constrain(m));
        }
        size_t A_j_1_max__ = m;
        size_t A_k_0_max__ = m;
        for (size_t j_1__ = 0; j_1__ < A_j_1_max__; ++j_1__) {
            for (size_t k_0__ = 0; k_0__ < A_k_0_max__; ++k_0__) {
                vars__.push_back(A[k_0__](j_1__));
            }
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_ZIPHMM_INF_1";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "theta";
        param_names__.push_back(param_name_stream__.str());
        size_t lambda_j_1_max__ = m;
        for (size_t j_1__ = 0; j_1__ < lambda_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lambda" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t A_j_1_max__ = m;
        size_t A_k_0_max__ = m;
        for (size_t j_1__ = 0; j_1__ < A_j_1_max__; ++j_1__) {
            for (size_t k_0__ = 0; k_0__ < A_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "A" << '.' << k_0__ + 1 << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "theta";
        param_names__.push_back(param_name_stream__.str());
        size_t lambda_j_1_max__ = m;
        for (size_t j_1__ = 0; j_1__ < lambda_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lambda" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t A_j_1_max__ = (m - 1);
        size_t A_k_0_max__ = m;
        for (size_t j_1__ = 0; j_1__ < A_j_1_max__; ++j_1__) {
            for (size_t k_0__ = 0; k_0__ < A_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "A" << '.' << k_0__ + 1 << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_ZIPHMM_INF_1_namespace::model_ZIPHMM_INF_1 stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
