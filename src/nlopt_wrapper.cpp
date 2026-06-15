#include "nlopt_wrapper.h"

// ---------------------------------------------------------------------------------------
// Algorithm naming

// List of all supported algorithms and associated names
struct AlgorithmNameAssociation {
    const char * name;
    nlopt_algorithm enum_value;
};
static const AlgorithmNameAssociation supported_algorithms[] = {
    {"LBFGS", NLOPT_LD_LBFGS},
    {"VAR1", NLOPT_LD_VAR1},
    {"VAR2", NLOPT_LD_VAR2},
    {"TNEWTON", NLOPT_LD_TNEWTON},
    {"TNEWTON_RESTART", NLOPT_LD_TNEWTON_RESTART},
    {"TNEWTON_PRECOND", NLOPT_LD_TNEWTON_PRECOND},
    {"TNEWTON_PRECOND_RESTART", NLOPT_LD_TNEWTON_PRECOND_RESTART},
    {"MMA", NLOPT_LD_MMA},
    {"CCSAQ", NLOPT_LD_CCSAQ},
};

// Retrieve the algorithm enum value associated to 'name', or throw an error
static nlopt_algorithm algorithm_from_name(const std::string & name) {
    for(const AlgorithmNameAssociation & association : supported_algorithms) {
        if(name == association.name) {
            return association.enum_value;
        }
    }
    // None matched
    std::string msg;
    msg += "Unsupported algorithm name: \"";
    msg += name;
    msg += "\"\nSupported:";
    for(const AlgorithmNameAssociation & association : supported_algorithms) {
        msg += " ";
        msg += association.name;
    }
    throw Rcpp::exception(msg.c_str());
}

// ---------------------------------------------------------------------------------------
// nlopt wrapper non-template helpers

std::unique_ptr<NloptStruct, NloptDeleter> new_nlopt_optimizer(const Rcpp::List & config, std::size_t size) {
    auto algorithm = algorithm_from_name(Rcpp::as<std::string>(config["algorithm"]));
    auto opt = std::unique_ptr<NloptStruct, NloptDeleter>(nlopt_create(algorithm, size));
    if(!opt) {
        throw Rcpp::exception("nlopt_create");
    }

    if(config.containsElementNamed("xtol_abs")) {
      if(nlopt_set_xtol_abs1(opt.get(), Rcpp::as<double>(config["xtol_abs"])) != NLOPT_SUCCESS) {
        throw Rcpp::exception("nlopt_set_xtol_abs1");
      }
    }

    if(config.containsElementNamed("xtol_rel")) {
        if(nlopt_set_xtol_rel(opt.get(), Rcpp::as<double>(config["xtol_rel"])) != NLOPT_SUCCESS) {
            throw Rcpp::exception("nlopt_set_xtol_rel");
        }
    }

    if(config.containsElementNamed("ftol_abs")) {
        if(nlopt_set_ftol_abs(opt.get(), Rcpp::as<double>(config["ftol_abs"])) != NLOPT_SUCCESS) {
            throw Rcpp::exception("nlopt_set_ftol_abs");
        }
    }
    if(config.containsElementNamed("ftol_rel")) {
        if(nlopt_set_ftol_rel(opt.get(), Rcpp::as<double>(config["ftol_rel"])) != NLOPT_SUCCESS) {
            throw Rcpp::exception("nlopt_set_ftol_rel");
        }
    }
    if(config.containsElementNamed("maxeval")) {
        if(nlopt_set_maxeval(opt.get(), Rcpp::as<int>(config["maxeval"])) != NLOPT_SUCCESS) {
            throw Rcpp::exception("nlopt_set_maxeval");
        }
    }
    if(config.containsElementNamed("maxtime")) {
        if(nlopt_set_maxtime(opt.get(), Rcpp::as<double>(config["maxtime"])) != NLOPT_SUCCESS) {
            throw Rcpp::exception("nlopt_set_maxtime");
        }
    }

    return opt;
}

// ---------------------------------------------------------------------------------------
// sanity test and example

// [[Rcpp::export]]
bool cpp_test_nlopt() {
    bool success = true;
    auto check = [&success](bool check_value, const char * context) {
        if(!check_value) {
            REprintf("Cpp internals failed: %s\n", context);
            success = false;
        }
    };

    // min_x x^2 -> should be 0. Does not uses packer due to only 1 variable.
    auto config = Rcpp::List::create(
        Rcpp::Named("algorithm", "LBFGS"),
        Rcpp::Named("xtol_rel", 1e-12),
        Rcpp::Named("ftol_abs", 0.0),
        Rcpp::Named("ftol_rel", 0.0),
        Rcpp::Named("xtol_abs", 0.0),
        Rcpp::Named("maxeval",  200),
        Rcpp::Named("maxtime",  100.));

    check(config.containsElementNamed("xtol_rel"), "config parsing using containsElementNamed");

    auto x = std::vector<double>{1.5, -2};

    auto optimizer = new_nlopt_optimizer(config, x.size());

    check(nlopt_get_algorithm(optimizer.get()) == NLOPT_LD_LBFGS, "optim algorithm");
    check(nlopt_get_ftol_abs(optimizer.get()) == 0.0, "optim ftol_abs");
    check(nlopt_get_ftol_rel(optimizer.get()) == 0.0, "optim ftol_rel");
    check(nlopt_get_xtol_rel(optimizer.get()) == 1e-12, "optim xtol_rel");

    auto f_and_grad = [check](const double * x, double * grad) -> double {
        double x1sq = x[0] * x[0] ;
        double obj = 100*std::pow(x[1] - x1sq,2) + std::pow(1-x[0],2);

        grad[0] = -400*(x[1] - x1sq)*x[0] - 2*(1-x[0]);
        grad[1] = 200*(x[1] - x1sq);
        return obj ;

    };
    OptimizerResult r = minimize_objective_on_parameters(optimizer.get(), f_and_grad, x);

    check(std::abs(x[0]-1) < 1e-8, "optim convergence");
    check(r.status != NLOPT_FAILURE, "optim status");

    x = std::vector<double>{1.5, -2};
    r = minimize_objective_on_parameters(optimizer.get(), f_and_grad, x);

    return success;
}
