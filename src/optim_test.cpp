#include <RcppArmadillo.h> // Includes R stuff
#include <nlopt.h>

#include <cstddef>
#include <functional> // lambda wrapping
#include <memory>
#include <string>
#include <tuple>       // packer system
#include <type_traits> // remove_pointer
#include <utility>     // move, forward

// ---------------------------------------------------------------------------------------
// Misc

inline arma::vec logfact(arma::mat y) {
    y.replace(0., 1.);
    return sum(y % log(y) - y + log(8 * pow(y, 3) + 4 * pow(y, 2) + y + 1. / 30.) / 6. + std::log(M_PI) / 2., 1);
}

inline arma::vec ki(arma::mat y) {
    arma::uword p = y.n_cols;
    return -logfact(std::move(y)) + 0.5 * (1. + (1. - double(p)) * std::log(2. * M_PI));
}

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

// Retrieve the algorithm enum value associated to 'name', or throw an error
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

// Required configuration values for using an optimizer
struct OptimizerConfiguration {
    nlopt_algorithm algorithm; // Must be from the supported algorithm list.

    arma::vec xtol_abs; // of size packer.size
    double xtol_rel;

    double ftol_abs;
    double ftol_rel;

    int maxeval;
    double maxtime;

    // Build from R list (with named elements).
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
    std::function<double(const arma::vec & parameters, arma::vec & gradients)> objective_and_grad_fn //
) {
    if(!(config.xtol_abs.n_elem == parameters.n_elem)) {
        throw Rcpp::exception("config.xtol_abs size");
    }

    // Create optimizer, stored in a unique_ptr to ensure automatic destruction.
    using Optimizer = std::remove_pointer<nlopt_opt>::type; // Retrieve struct type hidden by nlopt_opt typedef
    struct Deleter {
        void operator()(Optimizer * optimizer) const { nlopt_destroy(optimizer); }
    };
    auto optimizer = std::unique_ptr<Optimizer, Deleter>(nlopt_create(config.algorithm, parameters.n_elem));
    if(!optimizer) {
        throw Rcpp::exception("nlopt_create");
    }

    // Set optimizer configuration, with error checking
    auto check = [](nlopt_result r, const char * reason) {
        if(r != NLOPT_SUCCESS) {
            throw Rcpp::exception(reason);
        }
    };
    check(nlopt_set_xtol_abs(optimizer.get(), config.xtol_abs.memptr()), "nlopt_set_xtol_abs");
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
        std::function<double(const arma::vec &, arma::vec &)> objective_and_grad_fn;
    };
    OptimData optim_data = {0, std::move(objective_and_grad_fn)};

    auto optim_fn = [](unsigned n, const double * x, double * grad, void * data) -> double {
        // Wrap raw C arrays from nlopt into arma::vec (no copy)
        const auto parameters = arma::vec(const_cast<double *>(x), n, false, true);
        auto grad_storage = arma::vec(grad, n, false, true);
        // Restore optim_data and use it to perform computation step
        OptimData & optim_data = *static_cast<OptimData *>(data);
        optim_data.nb_iterations += 1;
        return optim_data.objective_and_grad_fn(parameters, grad_storage);
    };
    if(nlopt_set_min_objective(optimizer.get(), optim_fn, &optim_data) != NLOPT_SUCCESS) {
        throw Rcpp::exception("nlopt_set_min_objective");
    }

    double objective = 0.;
    nlopt_result status = nlopt_optimize(optimizer.get(), parameters.memptr(), &objective);
    return OptimizerResult{status, objective, optim_data.nb_iterations};
}

// ---------------------------------------------------------------------------------------
// Packing / unpacking utils

// Stores type, dimensions and offset for a single T object
// Must be specialised ; see specialisations for arma::vec/arma::mat below
//
// Required API when specialised:
// Constructor(const T & reference_object, arma::uword & current_offset);
// T unpack(const arma::vec & packed_storage);
// void pack(arma::vec & packed_storage, "T-like arma expression type" expr);
// void pack_double_or_arma(arma::vec & packed_storage, SEXP r_value); (see OptimizerConfiguration)
template <typename T> struct PackedInfo;

// All following implementation use vec.subvec() to access slices of the packed vector.
// The "prefered" way to give indexces is using span(start, end).
// However end is inclusive and unsigned, which causes underflow for 0-sized span of offset 0.
// Thus I use the less intuitive (but correct) form: subvec(offset, arma::size(size, 1))

template <> struct PackedInfo<arma::vec> {
    arma::uword offset;
    arma::uword size;

    PackedInfo(const arma::vec & v, arma::uword & current_offset) {
        offset = current_offset;
        size = v.n_elem;
        current_offset += size;
    }

    arma::vec unpack(const arma::vec & packed) const { return packed.subvec(offset, arma::size(size, 1)); }

    template <typename Expr> void pack(arma::vec & packed, Expr && expr) const {
        packed.subvec(offset, arma::size(size, 1)) = std::forward<Expr>(expr);
    }

    void pack_double_or_arma(arma::vec & packed, SEXP r_value) const {
        if(Rcpp::is<double>(r_value)) {
            packed.subvec(offset, arma::size(size, 1)).fill(Rcpp::as<double>(r_value));
        } else {
            pack(packed, Rcpp::as<arma::vec>(r_value));
        }
    }
};

template <> struct PackedInfo<arma::mat> {
    arma::uword offset;
    arma::uword rows;
    arma::uword cols;

    PackedInfo(const arma::mat & m, arma::uword & current_offset) {
        offset = current_offset;
        rows = m.n_rows;
        cols = m.n_cols;
        current_offset += rows * cols;
    }

    arma::mat unpack(const arma::vec & packed) const {
        return arma::reshape(packed.subvec(offset, arma::size(rows * cols, 1)), arma::size(rows, cols));
    }

    template <typename Expr> void pack(arma::vec & packed, Expr && expr) const {
        // Handles: mat expressions, vec expressions
        packed.subvec(offset, arma::size(rows * cols, 1)) = arma::vectorise(std::forward<Expr>(expr));
    }

    void pack_double_or_arma(arma::vec & packed, SEXP r_value) const {
        if(Rcpp::is<double>(r_value)) {
            packed.subvec(offset, arma::size(rows * cols, 1)).fill(Rcpp::as<double>(r_value));
        } else {
            pack(packed, Rcpp::as<arma::mat>(r_value));
        }
    }
};

// Stores packing information for multiple objects of types T1,T2,...,TN.
// Created (with type deduction) using make_packer(T1, ..., TN) below.
template <typename... Types> struct Packer {
    std::tuple<PackedInfo<Types>...> elements; // Packing info for each element (offset, type, dimensions)
    arma::uword size;                          // Total number of packed elements

    // Allows nicer syntax for unpack/pack/pack_double_or_arma
    template <std::size_t Index> auto unpack(const arma::vec & packed) const
        -> decltype(std::get<Index>(elements).unpack(packed)) {
        return std::get<Index>(elements).unpack(packed);
    }
    template <std::size_t Index, typename Expr> void pack(arma::vec & packed, Expr && expr) const {
        std::get<Index>(elements).pack(packed, std::forward<Expr>(expr));
    }
    template <std::size_t Index> void pack_double_or_arma(arma::vec & packed, SEXP r_value) const {
        std::get<Index>(elements).pack_double_or_arma(packed, r_value);
    }
};

template <typename... Types> Packer<Types...> make_packer(const Types &... values) {
    // Initialize Packer<Types...> using brace init, which guarantees evaluation order (required here !).
    // Will call each PackedInfo<T> constructor in order, increasing offset.
    // Then the final offset value will be copied into the size field.
    arma::uword current_offset = 0;
    return {
        {PackedInfo<Types>(values, current_offset)...}, // Increases offset sequentially for each value
        current_offset,                                 // Final value
    };
}

// ---------------------------------------------------------------------------------------
// Fully parametrized covariance

// [[Rcpp::export]]
Rcpp::List cpp_optimize_full(
    const Rcpp::List & init_parameters, // List(Theta, M, S)
    const arma::mat & y,                // responses (n,p)
    const arma::mat & x,                // covariates (n,d)
    const arma::mat & o,                // offsets (n,p)
    const arma::vec & w,                // weights (n)
    const Rcpp::List & configuration    // OptimizerConfiguration
) {
    // Conversion from R, prepare optimization
    const auto init_theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_m = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,p)
    const auto init_s = Rcpp::as<arma::mat>(init_parameters["S"]);         // (n,p)

    const auto packer = make_packer(init_theta, init_m, init_s);
    enum { THETA, M, S }; // Nice names for packer indexes

    auto parameters = arma::vec(packer.size);
    packer.pack<THETA>(parameters, init_theta);
    packer.pack<M>(parameters, init_m);
    packer.pack<S>(parameters, init_s);

    auto pack_xtol_abs = [&packer](arma::vec & packed, Rcpp::List list) {
        packer.pack_double_or_arma<THETA>(packed, list["Theta"]);
        packer.pack_double_or_arma<M>(packed, list["M"]);
        packer.pack_double_or_arma<S>(packed, list["S"]);
    };
    const auto config = OptimizerConfiguration::from_r_list(configuration, packer.size, pack_xtol_abs);

    const double w_bar = accu(w);

    // Optimize
    auto objective_and_grad =
        [&packer, &y, &x, &o, &w, &w_bar](const arma::vec & parameters, arma::vec & grad_storage) -> double {
        arma::mat theta = packer.unpack<THETA>(parameters);
        arma::mat m = packer.unpack<M>(parameters);
        arma::mat s = packer.unpack<S>(parameters);

        arma::mat s2 = s % s;
        arma::mat z = o + x * theta.t() + m;
        arma::mat a = exp(z + 0.5 * s2);
        arma::mat omega = w_bar * inv_sympd(m.t() * (m.each_col() % w) + diagmat(w.t() * s2));
        double objective = accu(w.t() * (a - y % z - 0.5 * log(s2))) - 0.5 * w_bar * real(log_det(omega));

        packer.pack<THETA>(grad_storage, (a - y).t() * (x.each_col() % w));
        packer.pack<M>(grad_storage, diagmat(w) * (m * omega + a - y));
        packer.pack<S>(grad_storage, diagmat(w) * (s.each_row() % diagvec(omega).t() + s % a - pow(s, -1)));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(parameters, config, objective_and_grad);

    // Variational parameters
    arma::mat m = packer.unpack<M>(parameters);
    arma::mat s = packer.unpack<S>(parameters);
    arma::mat s2 = s % s;
    // Regression parameters
    arma::mat theta = packer.unpack<THETA>(parameters);
    // Variance parameters
    arma::mat sigma = (1. / w_bar) * (m.t() * (m.each_col() % w) + diagmat(sum(s2.each_col() % w, 0)));
    arma::mat omega = inv_sympd(sigma);
    // Element-wise log-likehood
    arma::mat z = o + x * theta.t() + m;
    arma::mat a = exp(z + 0.5 * s2);
    arma::vec loglik = sum(y % z - a + 0.5 * log(s2) - 0.5 * ((m * omega) % m + s2 * diagmat(omega)), 1) +
                       0.5 * real(log_det(omega)) + ki(y);

    return Rcpp::List::create(
        Rcpp::Named("status", static_cast<int>(result.status)),
        Rcpp::Named("iterations", result.nb_iterations),
        Rcpp::Named("Theta", theta),
        Rcpp::Named("M", m),
        Rcpp::Named("S", s),
        Rcpp::Named("Z", z),
        Rcpp::Named("A", a),
        Rcpp::Named("Sigma", sigma),
        Rcpp::Named("Omega", omega),
        Rcpp::Named("loglik", loglik));
}

// ---------------------------------------------------------------------------------------
// Spherical covariance

// [[Rcpp::export]]
Rcpp::List cpp_optimize_spherical(
    const Rcpp::List & init_parameters, // List(Theta, M, S)
    const arma::mat & y,                // responses (n,p)
    const arma::mat & x,                // covariates (n,d)
    const arma::mat & o,                // offsets (n,p)
    const arma::vec & w,                // weights (n)
    const Rcpp::List & configuration    // OptimizerConfiguration
) {
    // Conversion from R, prepare optimization
    const auto init_theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_m = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,p)
    const auto init_s = Rcpp::as<arma::vec>(init_parameters["S"]);         // (n)

    const auto packer = make_packer(init_theta, init_m, init_s);
    enum { THETA, M, S }; // Nice names for packer indexes

    auto parameters = arma::vec(packer.size);
    packer.pack<THETA>(parameters, init_theta);
    packer.pack<M>(parameters, init_m);
    packer.pack<S>(parameters, init_s);

    auto pack_xtol_abs = [&packer](arma::vec & packed, Rcpp::List list) {
        packer.pack_double_or_arma<THETA>(packed, list["Theta"]);
        packer.pack_double_or_arma<M>(packed, list["M"]);
        packer.pack_double_or_arma<S>(packed, list["S"]);
    };
    const auto config = OptimizerConfiguration::from_r_list(configuration, packer.size, pack_xtol_abs);

    const double w_bar = accu(w);

    // Optimize
    auto objective_and_grad =
        [&packer, &o, &x, &y, &w, &w_bar](const arma::vec & parameters, arma::vec & grad_storage) -> double {
        arma::mat theta = packer.unpack<THETA>(parameters);
        arma::mat m = packer.unpack<M>(parameters);
        arma::vec s = packer.unpack<S>(parameters);

        arma::vec s2 = s % s;
        const arma::uword p = y.n_cols;
        arma::mat z = o + x * theta.t() + m;
        arma::mat a = exp(z.each_col() + 0.5 * s2);
        double sigma2 = arma::as_scalar(accu(m % (m.each_col() % w)) / (w_bar * double(p)) + accu(w % s2) / w_bar);
        double objective = accu(diagmat(w) * (a - y % z)) - 0.5 * double(p) * accu(w % log(s2)) +
                           0.5 * w_bar * double(p) * log(sigma2);

        packer.pack<THETA>(grad_storage, (a - y).t() * (x.each_col() % w));
        packer.pack<M>(grad_storage, diagmat(w) * (m / sigma2 + a - y));
        packer.pack<S>(grad_storage, w % (s % sum(a, 1) - double(p) * pow(s, -1) - double(p) * s / sigma2));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(parameters, config, objective_and_grad);

    // Variational parameters
    arma::mat m = packer.unpack<M>(parameters);
    arma::mat s = packer.unpack<S>(parameters); // vec(n) -> mat(n, 1)
    arma::vec s2 = s % s;
    // Regression parameters
    arma::mat theta = packer.unpack<THETA>(parameters);
    // Variance parameters
    const arma::uword p = y.n_cols;
    const double n_sigma2 = arma::as_scalar(dot(w, sum(pow(m, 2), 1) + double(p) * s2));
    const double sigma2 = n_sigma2 / (double(p) * w_bar);
    arma::mat sigma = arma::eye(p, p) * sigma2;
    arma::mat omega = arma::eye(p, p) * pow(sigma2, -1);
    // Element-wise log-likelihood
    arma::mat z = o + x * theta.t() + m;
    arma::mat a = exp(z.each_col() + 0.5 * s2);
    arma::mat loglik = sum(y % z - a - 0.5 * pow(m, 2) / sigma2, 1) - double(p) * s / sigma2 +
                       0.5 * double(p) * log(s2 / sigma2) + ki(y);

    return Rcpp::List::create(
        Rcpp::Named("status", static_cast<int>(result.status)),
        Rcpp::Named("iterations", result.nb_iterations),
        Rcpp::Named("Theta", theta),
        Rcpp::Named("M", m),
        Rcpp::Named("S", s),
        Rcpp::Named("Z", z),
        Rcpp::Named("A", a),
        Rcpp::Named("Sigma", sigma),
        Rcpp::Named("Omega", omega),
        Rcpp::Named("loglik", loglik));
}

// ---------------------------------------------------------------------------------------
// Diagonal covariance

// [[Rcpp::export]]
Rcpp::List cpp_optimize_diagonal(
    const Rcpp::List & init_parameters, // List(Theta, M, S)
    const arma::mat & y,                // responses (n,p)
    const arma::mat & x,                // covariates (n,d)
    const arma::mat & o,                // offsets (n,p)
    const arma::vec & w,                // weights (n)
    const Rcpp::List & configuration    // OptimizerConfiguration
) {
    // Conversion from R, prepare optimization
    const auto init_theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_m = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,p)
    const auto init_s = Rcpp::as<arma::mat>(init_parameters["S"]);         // (n,p)

    const auto packer = make_packer(init_theta, init_m, init_s);
    enum { THETA, M, S }; // Nice names for packer indexes

    auto parameters = arma::vec(packer.size);
    packer.pack<THETA>(parameters, init_theta);
    packer.pack<M>(parameters, init_m);
    packer.pack<S>(parameters, init_s);

    auto pack_xtol_abs = [&packer](arma::vec & packed, Rcpp::List list) {
        packer.pack_double_or_arma<THETA>(packed, list["Theta"]);
        packer.pack_double_or_arma<M>(packed, list["M"]);
        packer.pack_double_or_arma<S>(packed, list["S"]);
    };
    const auto config = OptimizerConfiguration::from_r_list(configuration, packer.size, pack_xtol_abs);

    const double w_bar = accu(w);

    // Optimize
    auto objective_and_grad =
        [&packer, &o, &x, &y, &w, &w_bar](const arma::vec & parameters, arma::vec & grad_storage) -> double {
        arma::mat theta = packer.unpack<THETA>(parameters);
        arma::mat m = packer.unpack<M>(parameters);
        arma::mat s = packer.unpack<S>(parameters);

        arma::mat s2 = s % s;
        arma::mat z = o + x * theta.t() + m;
        arma::mat a = exp(z + 0.5 * s2);
        arma::rowvec diag_sigma = sum(m % (m.each_col() % w) + (s2.each_col() % w), 0) / w_bar;
        double objective = accu(diagmat(w) * (a - y % z - 0.5 * log(s2))) + 0.5 * w_bar * accu(log(diag_sigma));

        packer.pack<THETA>(grad_storage, (a - y).t() * (x.each_col() % w));
        packer.pack<M>(grad_storage, diagmat(w) * ((m.each_row() / diag_sigma) + a - y));
        packer.pack<S>(grad_storage, diagmat(w) * (s.each_row() % pow(diag_sigma, -1) + s % a - pow(s, -1)));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(parameters, config, objective_and_grad);

    // Variational parameters
    arma::mat m = packer.unpack<M>(parameters);
    arma::mat s = packer.unpack<S>(parameters);
    arma::mat s2 = s % s;
    // Regression parameters
    arma::mat theta = packer.unpack<THETA>(parameters);
    // Variance parameters
    arma::rowvec sigma2 = w.t() * (pow(m, 2) + s2) / w_bar;
    arma::vec omega2 = pow(sigma2.t(), -1);
    arma::mat sigma = diagmat(sigma2);
    arma::mat omega = diagmat(omega2);
    // Element-wise log-likelihood
    arma::mat z = o + x * theta.t() + m;
    arma::mat a = exp(z + 0.5 * s2);
    arma::mat loglik =
        sum(y % z - a + 0.5 * log(s2), 1) - 0.5 * (pow(m, 2) + s2) * omega2 + 0.5 * sum(log(omega2)) + ki(y);

    return Rcpp::List::create(
        Rcpp::Named("status", static_cast<int>(result.status)),
        Rcpp::Named("iterations", result.nb_iterations),
        Rcpp::Named("Theta", theta),
        Rcpp::Named("M", m),
        Rcpp::Named("S", s),
        Rcpp::Named("Z", z),
        Rcpp::Named("A", a),
        Rcpp::Named("Sigma", sigma),
        Rcpp::Named("Omega", omega),
        Rcpp::Named("loglik", loglik));
}

// ---------------------------------------------------------------------------------------
// Rank-constrained covariance

// Rank (q) is already determined by param dimensions ; not passed anywhere

// [[Rcpp::export]]
Rcpp::List cpp_optimize_rank(
    const Rcpp::List & init_parameters, // List(Theta, B, M, S)
    const arma::mat & y,                // responses (n,p)
    const arma::mat & x,                // covariates (n,d)
    const arma::mat & o,                // offsets (n,p)
    const arma::vec & w,                // weights (n)
    const Rcpp::List & configuration    // OptimizerConfiguration
) {
    // Conversion from R, prepare optimization
    const auto init_theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_b = Rcpp::as<arma::mat>(init_parameters["B"]);         // (p,q)
    const auto init_m = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,q)
    const auto init_s = Rcpp::as<arma::mat>(init_parameters["S"]);         // (n,q)

    const auto packer = make_packer(init_theta, init_b, init_m, init_s);
    enum { THETA, B, M, S }; // Nice names for packer indexes

    auto parameters = arma::vec(packer.size);
    packer.pack<THETA>(parameters, init_theta);
    packer.pack<B>(parameters, init_m);
    packer.pack<M>(parameters, init_m);
    packer.pack<S>(parameters, init_s);

    auto pack_xtol_abs = [&packer](arma::vec & packed, Rcpp::List list) {
        packer.pack_double_or_arma<THETA>(packed, list["Theta"]);
        packer.pack_double_or_arma<B>(packed, list["B"]);
        packer.pack_double_or_arma<M>(packed, list["M"]);
        packer.pack_double_or_arma<S>(packed, list["S"]);
    };
    const auto config = OptimizerConfiguration::from_r_list(configuration, packer.size, pack_xtol_abs);

    // Optimize
    auto objective_and_grad =
        [&packer, &o, &x, &y, &w](const arma::vec & parameters, arma::vec & grad_storage) -> double {
        arma::mat theta = packer.unpack<THETA>(parameters);
        arma::mat b = packer.unpack<B>(parameters);
        arma::mat m = packer.unpack<M>(parameters);
        arma::mat s = packer.unpack<S>(parameters);

        arma::mat s2 = s % s;
        arma::mat z = o + x * theta.t() + m * b.t();
        arma::mat a = exp(z + 0.5 * s2 * (b % b).t());
        double objective = accu(diagmat(w) * (a - y % z)) + 0.5 * accu(diagmat(w) * (m % m + s2 - log(s2) - 1.));

        packer.pack<THETA>(grad_storage, (a - y).t() * (x.each_col() % w));
        packer.pack<B>(grad_storage, (diagmat(w) * (a - y)).t() * m + (a.t() * (s2.each_col() % w)) % b);
        packer.pack<M>(grad_storage, diagmat(w) * ((a - y) * b + m));
        packer.pack<S>(grad_storage, diagmat(w) * (s - pow(s, -1) + a * (b % b) % s));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(parameters, config, objective_and_grad);

    // Model and variational parameters
    arma::mat theta = packer.unpack<THETA>(parameters);
    arma::mat b = packer.unpack<B>(parameters);
    arma::mat m = packer.unpack<M>(parameters);
    arma::mat s = packer.unpack<S>(parameters);
    arma::mat s2 = s % s;
    arma::mat sigma = b * (m.t() * (m.each_col() % w) + diagmat(sum(s2.each_col() % w, 0))) * b.t() / accu(w);
    // Element-wise log-likelihood
    arma::mat z = o + x * theta.t() + m * b.t();
    arma::mat a = exp(z + 0.5 * s2 * (b % b).t());
    arma::mat loglik = arma::sum(y % z - a, 1) - 0.5 * sum(m % m + s2 - log(s2) - 1., 1) + ki(y);

    return Rcpp::List::create(
        Rcpp::Named("status", static_cast<int>(result.status)),
        Rcpp::Named("iterations", result.nb_iterations),
        Rcpp::Named("Theta", theta),
        Rcpp::Named("B", b),
        Rcpp::Named("M", m),
        Rcpp::Named("S", s),
        Rcpp::Named("Z", z),
        Rcpp::Named("A", a),
        Rcpp::Named("Sigma", sigma),
        Rcpp::Named("loglik", loglik));
}

// ---------------------------------------------------------------------------------------
// Sparse inverse covariance

// [[Rcpp::export]]
Rcpp::List cpp_optimize_sparse(
    const Rcpp::List & init_parameters, // List(Theta, M, S)
    const arma::mat & y,                // responses (n,p)
    const arma::mat & x,                // covariates (n,d)
    const arma::mat & o,                // offsets (n,p)
    const arma::vec & w,                // weights (n)
    const arma::mat & omega,            // covinv (p,p)
    const Rcpp::List & configuration    // OptimizerConfiguration
) {
    // Conversion from R, prepare optimization
    const auto init_theta = Rcpp::as<arma::mat>(init_parameters["Theta"]); // (p,d)
    const auto init_m = Rcpp::as<arma::mat>(init_parameters["M"]);         // (n,p)
    const auto init_s = Rcpp::as<arma::mat>(init_parameters["S"]);         // (n,p)

    const auto packer = make_packer(init_theta, init_m, init_s);
    enum { THETA, M, S }; // Nice names for packer indexes

    auto parameters = arma::vec(packer.size);
    packer.pack<THETA>(parameters, init_theta);
    packer.pack<M>(parameters, init_m);
    packer.pack<S>(parameters, init_s);

    auto pack_xtol_abs = [&packer](arma::vec & packed, Rcpp::List list) {
        packer.pack_double_or_arma<THETA>(packed, list["Theta"]);
        packer.pack_double_or_arma<M>(packed, list["M"]);
        packer.pack_double_or_arma<S>(packed, list["S"]);
    };
    const auto config = OptimizerConfiguration::from_r_list(configuration, packer.size, pack_xtol_abs);

    // Optimize
    auto objective_and_grad =
        [&packer, &o, &x, &y, &w, &omega](const arma::vec & parameters, arma::vec & grad_storage) -> double {
        arma::mat theta = packer.unpack<THETA>(parameters);
        arma::mat m = packer.unpack<M>(parameters);
        arma::mat s = packer.unpack<S>(parameters);

        arma::mat s2 = s % s;
        arma::mat z = o + x * theta.t() + m;
        arma::mat a = exp(z + 0.5 * s);
        arma::mat nSigma = m.t() * (m.each_col() % w) + diagmat(w.t() * s2);
        double objective = accu(w.t() * (a - y % z - 0.5 * log(s2))) - trace(omega * nSigma);

        packer.pack<THETA>(grad_storage, (a - y).t() * (x.each_col() % w));
        packer.pack<M>(grad_storage, diagmat(w) * (m * omega + a - y));
        packer.pack<S>(grad_storage, diagmat(w) * (s.each_row() % diagvec(omega).t() + s % a - pow(s, -1)));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(parameters, config, objective_and_grad);

    // Model and variational parameters
    arma::mat theta = packer.unpack<THETA>(parameters);
    arma::mat m = packer.unpack<M>(parameters);
    arma::mat s = packer.unpack<S>(parameters);
    arma::mat s2 = s % s;
    arma::mat sigma = (m.t() * (m.each_col() % w) + diagmat(w.t() * s2)) / accu(w);
    // Element-wise log-likelihood
    arma::mat z = o + x * theta.t() + m;
    arma::mat a = exp(z + 0.5 * s2);
    arma::mat loglik = sum(y % z - a - 0.5 * ((m * omega) % m - log(s2) + s2 * diagmat(omega)), 1) +
                       0.5 * real(log_det(omega)) + ki(y);

    return Rcpp::List::create(
        Rcpp::Named("status", static_cast<int>(result.status)),
        Rcpp::Named("iterations", result.nb_iterations),
        Rcpp::Named("Theta", theta),
        Rcpp::Named("M", m),
        Rcpp::Named("S", s),
        Rcpp::Named("Z", z),
        Rcpp::Named("A", a),
        Rcpp::Named("Sigma", sigma),
        Rcpp::Named("loglik", loglik));
}

// ---------------------------------------------------------------------------------------
// Single VE step
// TODO adapt all variants

// [[Rcpp::export]]
Rcpp::List cpp_optimize_ve(
    const Rcpp::List & init_parameters, // List(M, S)
    const arma::mat & y,                // responses (n,p)
    const arma::mat & x,                // covariates (n,d)
    const arma::mat & o,                // offsets (n,p)
    const arma::mat & theta,            // regression_parameters (p,d)
    const arma::mat & sigma,            // (p,p)
    const Rcpp::List & configuration    // OptimizerConfiguration
) {
    // Conversion from R, prepare optimization
    const auto init_m = Rcpp::as<arma::mat>(init_parameters["M"]); // (n,p)
    const auto init_s = Rcpp::as<arma::mat>(init_parameters["S"]); // (n,p)

    const auto packer = make_packer(init_m, init_s);
    enum { M, S }; // Nice names for packer indexes

    auto parameters = arma::vec(packer.size);
    packer.pack<M>(parameters, init_m);
    packer.pack<S>(parameters, init_s);

    auto pack_xtol_abs = [&packer](arma::vec & packed, Rcpp::List list) {
        packer.pack_double_or_arma<M>(packed, list["M"]);
        packer.pack_double_or_arma<S>(packed, list["S"]);
    };
    const auto config = OptimizerConfiguration::from_r_list(configuration, packer.size, pack_xtol_abs);

    const arma::mat omega = inv_sympd(sigma); // (p,p)
    const double log_det_omega = real(log_det(omega));

    // Optimize
    auto objective_and_grad = [&packer, &o, &x, &y, &theta, &omega, &log_det_omega](
                                  const arma::vec & parameters, arma::vec & grad_storage) -> double {
        arma::mat m = packer.unpack<M>(parameters);
        arma::mat s = packer.unpack<S>(parameters);

        const arma::uword n = y.n_rows;
        arma::mat z = o + x * theta.t() + m;
        arma::mat a = exp(z + 0.5 * s);
        // 0.5 tr(\Omega M'M) + 0.5 tr(\bar{S} \Omega)
        double prior = 0.5 * accu(omega % (m.t() * m)) + 0.5 * dot(arma::ones(n).t() * s, diagvec(omega));
        // J(M, S, \Theta, \Omega, Y, X, O)
        double objective = accu(a - y % z - 0.5 * log(s)) + prior - 0.5 * double(n) * log_det_omega;

        packer.pack<M>(grad_storage, m * omega + a - y);
        packer.pack<S>(grad_storage, 0.5 * (arma::ones(n) * diagvec(omega).t() + a - 1. / s));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(parameters, config, objective_and_grad);

    // Post process
    arma::mat m = packer.unpack<M>(parameters);
    arma::mat s = packer.unpack<S>(parameters);

    const arma::uword p = y.n_cols;
    arma::mat z = o + x * theta.t() + m;
    arma::mat a = exp(z + 0.5 * s);
    arma::vec loglik = arma::sum(y % z - a + 0.5 * log(s) - 0.5 * ((m * omega) % m + s * diagmat(omega)), 1) +
                       0.5 * log_det_omega - logfact(y) + 0.5 * double(p);

    return Rcpp::List::create(
        Rcpp::Named("status", static_cast<int>(result.status)),
        Rcpp::Named("iterations", result.nb_iterations),
        Rcpp::Named("objective", result.objective + accu(logfact(y))),
        Rcpp::Named("M", m),
        Rcpp::Named("S", s),
        Rcpp::Named("loglik", loglik));
}

// TODO : New model

// ---------------------------------------------------------------------------------------
// Internals tests

// [[Rcpp::export]]
bool cpp_internal_tests() {
    bool success = true;
    auto check = [&success](bool check_value, const char * context) {
        if(!check_value) {
            REprintf("Cpp internals failed: %s", context);
            success = false;
        }
    };
    const double epsilon = 1e-6;

    // Test packing / unpacking
    auto z = arma::vec{}; // 0-sized
    auto a = arma::mat(4, 10, arma::fill::randu);
    auto b = arma::vec(7, arma::fill::randu);

    const auto packer = make_packer(z, a, b, b);
    check(packer.size == 4 * 10 + 7 + 7, "packer size computation");
    check(std::get<0>(packer.elements).offset == 0, "packer offset 0");
    check(std::get<1>(packer.elements).offset == 0, "packer offset 1");
    check(std::get<2>(packer.elements).offset == 4 * 10, "packer offset 2");
    check(std::get<3>(packer.elements).offset == 4 * 10 + 7, "packer offset 3");

    auto packed = arma::vec(packer.size);
    packer.pack<0>(packed, z);
    packer.pack<1>(packed, a);
    packer.pack<2>(packed, b);
    packer.pack<3>(packed, b);
    check(packer.unpack<0>(packed).n_elem == 0, "unpack 0");
    check(arma::approx_equal(a, packer.unpack<1>(packed), "absdiff", epsilon), "unpack 1");
    check(arma::approx_equal(b, packer.unpack<2>(packed), "absdiff", epsilon), "unpack 2");
    check(arma::approx_equal(b, packer.unpack<3>(packed), "absdiff", epsilon), "unpack 3");

    packer.pack_double_or_arma<1>(packed, Rcpp::wrap(0.));
    check(packer.unpack<1>(packed).is_zero(), "pack_double_or_arma double(0.) in mat");
    packer.pack_double_or_arma<1>(packed, Rcpp::wrap(a));
    check(arma::approx_equal(a, packer.unpack<1>(packed), "absdiff", epsilon), "pack_double_or_arma mat");

    packer.pack_double_or_arma<2>(packed, Rcpp::wrap(0.));
    check(packer.unpack<2>(packed).is_zero(), "pack_double_or_arma double(0.) in vec");
    packer.pack_double_or_arma<2>(packed, Rcpp::wrap(b));
    check(arma::approx_equal(b, packer.unpack<2>(packed), "absdiff", epsilon), "pack_double_or_arma vec");

    // Test nlopt wrapper
    // min_x x^2 -> should be 0. Does not uses packer due to only 1 variable.
    auto config = OptimizerConfiguration{
        algorithm_from_name("LBFGS"),
        arma::vec{epsilon}, // xtol_abs
        epsilon,            // xtol_rel
        epsilon,            // ftol_abs
        epsilon,            // ftol_rel
        100,                // maxeval
        100.,               // maxtime
    };
    auto x = arma::vec{42.};
    auto f_and_grad = [check](const arma::vec & x, arma::vec & grad) -> double {
        check(x.n_elem == 1, "opt x size");
        check(grad.n_elem == 1, "opt grad size");
        double v = x[0];
        grad[0] = 2. * v;
        return v * v;
    };
    OptimizerResult r = minimize_objective_on_parameters(x, config, f_and_grad);
    check(arma::approx_equal(x, arma::vec{0.}, "absdiff", 10. * epsilon), "optim convergence");
    check(r.status != NLOPT_FAILURE, "optim status");

    return success;
}