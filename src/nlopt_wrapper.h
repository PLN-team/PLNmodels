#include <RcppArmadillo.h>
#include <nlopt.h>

#include <functional> // lambda wrapping
#include <string>
#include <utility> // move

// Retrieve the algorithm enum value associated to 'name', or throw an error
nlopt_algorithm algorithm_from_name(const std::string & name);

// Required configuration values for using an optimizer
struct OptimizerConfiguration {
    nlopt_algorithm algorithm; // Must be from the supported algorithm list.

    arma::vec xtol_abs; // of size packer.size
    double xtol_rel;

    double ftol_abs;
    double ftol_rel;

    int maxeval;
    double maxtime;

    // Build configuration from R list (with named elements).
    //
    // xtol_abs has special handling, due to having values for each parameter element.
    // 2 modes are supported, depending on the input type:
    // - single float value: use single value for all parameter elements
    // - named list with values for each parameter: use parameter-specific values
    // The parameter-specific mode delegates extraction of the xtol_abs list and packing to xtol_abs array to
    // a user sub function of prototype: void pack_xtol_abs(arma::vec & packed, Rcpp::List xtol_abs_list).
    // These subfunctions are generally specific to the algorithm variant (parameter list types).
    // They use PackedInfo::pack_double_or_arma() to accept as input, for each parameter:
    // - a single double value: use this value for all elements of the parameter
    // - an arma mat/vec with the parameter dimensions: use element-specific values
    template <typename F>
    static OptimizerConfiguration from_r_list(const Rcpp::List & list, arma::uword packer_size, F pack_xtol_abs) {
        // Special handling for xtol_abs
        auto xtol_abs = arma::vec(packer_size);
        SEXP xtol_abs_r_value = list["xtol_abs"];
        if(Rcpp::is<double>(xtol_abs_r_value)) {
            xtol_abs.fill(Rcpp::as<double>(xtol_abs_r_value));
        } else if(Rcpp::is<Rcpp::List>(xtol_abs_r_value)) {
            pack_xtol_abs(xtol_abs, Rcpp::as<Rcpp::List>(xtol_abs_r_value));
        } else {
            throw Rcpp::exception("unsupported config[xtol_abs] type: must be double or list of by-parameter values");
        }
        // All others
        return {
            algorithm_from_name(Rcpp::as<std::string>(list["algorithm"])),

            std::move(xtol_abs),
            Rcpp::as<double>(list["xtol_rel"]),

            Rcpp::as<double>(list["ftol_abs"]),
            Rcpp::as<double>(list["ftol_rel"]),

            Rcpp::as<int>(list["maxeval"]),
            Rcpp::as<double>(list["maxtime"]),
        };
    }
};

// Return value of an optimizer call.
struct OptimizerResult {
    nlopt_result status;
    double objective;
    int nb_iterations;
};

// Find parameters minimizing the given objective function, under the given configuration.
OptimizerResult minimize_objective_on_parameters(
    // Parameters are modified in place
    arma::vec & parameters,
    const OptimizerConfiguration & config,
    // Computation step function (usually initialised with a stateful lambda / closure).
    // It should compute and return the objective value for the given parameters, and store computed gradients.
    // Both vectors are of size nb_parameters.
    std::function<double(const arma::vec & parameters, arma::vec & gradients)> objective_and_grad_fn);