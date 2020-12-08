#include "nlopt_wrapper.h"

#include <memory>      // unique_ptr
#include <type_traits> // remove_pointer

// This header DEFINES non inline functions that follow the declarations of nlopt.h
// It must be only included once in a project, or it will generate multiple definitions.
#include "nloptrAPI.h"

// ---------------------------------------------------------------------------------------
// Algorithm naming

// List of all supported algorithms and associated names
struct AlgorithmNameAssociation {
    const char * name;
    nlopt_algorithm enum_value;
};
static const AlgorithmNameAssociation supported_algorithms[] = {
    {"LBFGS_NOCEDAL", NLOPT_LD_LBFGS_NOCEDAL},
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

nlopt_algorithm algorithm_from_name(const std::string & name) {
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
// nlopt wrapper

OptimizerResult minimize_objective_on_parameters(
    std::vector<double> & parameters,
    const OptimizerConfiguration & config,
    std::function<double(const double * parameters, double * gradients)> objective_and_grad_fn //
) {
    if(config.xtol_abs.size() != parameters.size()) {
        throw Rcpp::exception("config.xtol_abs size");
    }

    // Create optimizer, stored in a unique_ptr to ensure automatic destruction.
    using Optimizer = std::remove_pointer<nlopt_opt>::type; // Retrieve struct type hidden by nlopt_opt typedef
    struct Deleter {
        void operator()(Optimizer * optimizer) const { nlopt_destroy(optimizer); }
    };
    auto optimizer = std::unique_ptr<Optimizer, Deleter>(nlopt_create(config.algorithm, parameters.size()));
    if(!optimizer) {
        throw Rcpp::exception("nlopt_create");
    }

    // Set optimizer configuration, with error checking
    auto check = [](nlopt_result r, const char * reason) {
        if(r != NLOPT_SUCCESS) {
            throw Rcpp::exception(reason);
        }
    };
    check(nlopt_set_xtol_abs(optimizer.get(), config.xtol_abs.data()), "nlopt_set_xtol_abs");
    check(nlopt_set_xtol_rel(optimizer.get(), config.xtol_rel), "nlopt_set_xtol_rel");
    check(nlopt_set_ftol_abs(optimizer.get(), config.ftol_abs), "nlopt_set_ftol_abs");
    check(nlopt_set_ftol_rel(optimizer.get(), config.ftol_rel), "nlopt_set_ftol_rel");
    check(nlopt_set_maxeval(optimizer.get(), config.maxeval), "nlopt_set_maxeval");
    check(nlopt_set_maxtime(optimizer.get(), config.maxtime), "nlopt_set_maxtime");

    // Optimize.
    // nlopt requires a pair of function pointer and custom data pointer for the objective function.
    // The OptimData struct stores iteration count and the step function ; it is used as void* data.
    // optim_fn is a stateless lambda, and can be cast to a function pointer as required by nlopt.
    // It is an adapter: convert nlopt raw arrays to arma values, and call the closure.
    struct OptimData {
        int nb_iterations;
        std::function<double(const double *, double *)> objective_and_grad_fn;
    };
    OptimData optim_data = {0, std::move(objective_and_grad_fn)};

    auto optim_fn = [](unsigned n, const double * x, double * grad, void * data) -> double {
        OptimData & optim_data = *static_cast<OptimData *>(data);
        optim_data.nb_iterations += 1;
        return optim_data.objective_and_grad_fn(x, grad);
    };
    if(nlopt_set_min_objective(optimizer.get(), optim_fn, &optim_data) != NLOPT_SUCCESS) {
        throw Rcpp::exception("nlopt_set_min_objective");
    }

    double objective = 0.;
    nlopt_result status = nlopt_optimize(optimizer.get(), parameters.data(), &objective);
    return OptimizerResult{status, objective, optim_data.nb_iterations};
}

// ---------------------------------------------------------------------------------------
// sanity test and example

// [[Rcpp::export]]
bool cpp_test_nlopt() {
    bool success = true;
    auto check = [&success](bool check_value, const char * context) {
        if(!check_value) {
            REprintf("Cpp internals failed: %s", context);
            success = false;
        }
    };
    const double epsilon = 1e-6;

    // min_x x^2 -> should be 0. Does not uses packer due to only 1 variable.
    auto config = OptimizerConfiguration{
        algorithm_from_name("LBFGS"),
        std::vector<double>{epsilon}, // xtol_abs
        epsilon,                      // xtol_rel
        epsilon,                      // ftol_abs
        epsilon,                      // ftol_rel
        100,                          // maxeval
        100.,                         // maxtime
    };
    auto x = std::vector<double>{42.};
    auto f_and_grad = [check](const double * x, double * grad) -> double {
        double v = x[0];
        grad[0] = 2. * v;
        return v * v;
    };
    OptimizerResult r = minimize_objective_on_parameters(x, config, f_and_grad);
    check(std::abs(x[0]) < epsilon, "optim convergence");
    check(r.status != NLOPT_FAILURE, "optim status");

    return success;
}