// Provide a nice wrapper for nloptr optimizer.
// It accepts a C++11 lambda anonymous closure for the computation of objective & gradients.
// This is easier than generating by hand the void* data and C style callback.
#pragma once

#include <Rcpp.h>
#include <nloptrAPI.h>

#include <memory>      // unique_ptr
#include <type_traits> // remove_pointer
#include <vector>

// Get the struct type in the "opaque pointer" nlopt_opt. Required by unique_ptr.
using NloptStruct = typename std::remove_pointer<nlopt_opt>::type;

struct NloptDeleter {
    void operator()(NloptStruct * optimizer) const { nlopt_destroy(optimizer); }
};

// Create an optimizer with the config in the List.
// Requires key 'algorithm' with algorithm name (string).
// Sets the following configuration scalar values if the keys are found in the list:
// xtol_rel, ftol_abs, ftol_rel, maxeval, maxtime.
std::unique_ptr<NloptStruct, NloptDeleter> new_nlopt_optimizer(const Rcpp::List & config, std::size_t size);

// Helpers to set xtol_abs (uniform or per-parameter packed array).
// This is not done by new_nlopt_optimizer as it may require packing values, which must be user specified.
void set_uniform_xtol_abs(NloptStruct * opt, double value);
void set_per_value_xtol_abs(NloptStruct * opt, const std::vector<double> & xtol_abs);

// Helpers to set x_weights (uniform or per-parameter packed array).
// This is not done by new_nlopt_optimizer as it may require packing values, which must be user specified.
// void set_uniform_x_weights(NloptStruct * opt, double value);
// void set_per_value_x_weights(NloptStruct * opt, const std::vector<double> & x_weigths);

struct OptimizerResult {
    nlopt_result status;
    double objective;
    int nb_iterations;
};

// Using a configured optimizer, find the minimum of the objective function by modifying the parameters.
template <typename ObjectiveAndGradFunc> OptimizerResult minimize_objective_on_parameters(
    // Configured optimizer
    NloptStruct * optimizer,
    // Function / closure, computes and return objective, stores gradients.
    // (const double * params, double * gradients) -> double.
    // params and gradients are packed.
    const ObjectiveAndGradFunc & objective_and_grad,
    // Packed parameters, modified in place.
    std::vector<double> & parameters //
) {
    if(parameters.size() != nlopt_get_dimension(optimizer)) {
        throw Rcpp::exception("minimize_objective: parameter size mismatch");
    }

    // nlopt requires a pair of function pointer and custom data pointer for the objective function.
    // The OptimData struct stores iteration count and the step function ; it is used as void* data.
    // optim_fn is a stateless lambda, and can be cast to a function pointer as required by nlopt.
    struct OptimData {
        int nb_iterations;
        const ObjectiveAndGradFunc & objective_and_grad;
    };
    OptimData optim_data = {0, objective_and_grad};

    auto optim_fn = [](unsigned n, const double * p, double * grad, void * data) -> double {
        OptimData & optim_data = *static_cast<OptimData *>(data);
        optim_data.nb_iterations += 1;
        return optim_data.objective_and_grad(p, grad);
    };
    if(nlopt_set_min_objective(optimizer, optim_fn, &optim_data) != NLOPT_SUCCESS) {
        throw Rcpp::exception("nlopt_set_min_objective");
    }

    double objective = 0.;
    nlopt_result status = nlopt_optimize(optimizer, parameters.data(), &objective);
    return {status, objective, optim_data.nb_iterations};
}
