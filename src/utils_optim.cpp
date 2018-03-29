#include "utils_optim.h"

using namespace Rcpp;

double K(arma::mat Y) {
  arma::colvec v = arma::nonzeros(Y);
  return accu(v % log(v) - v + log(8*pow(v,3) + 4*pow(v, 2) + v + 1/30)/6 + log(M_PI)/2);
}

// Convert string to nlopt_alogirthm
//
// restrict the choices to algorithms meaningful for PLN optimization
nlopt::algorithm getAlgorithmCode( const std::string algorithm_str) {

    nlopt::algorithm algorithm;

    if ( algorithm_str.compare("LBFGS_NOCEDAL") == 0 ) {
        algorithm = nlopt::LD_LBFGS_NOCEDAL;
    }
    else if ( algorithm_str.compare("LBFGS" ) == 0 ) {
        algorithm = nlopt::LD_LBFGS;
    }
    else if ( algorithm_str.compare("VAR1" ) == 0 ) {
        algorithm = nlopt::LD_VAR1;
    }
    else if ( algorithm_str.compare("VAR2" ) == 0 ) {
        algorithm = nlopt::LD_VAR2;
    }
    else if ( algorithm_str.compare("TNEWTON" ) == 0 ) {
        algorithm = nlopt::LD_TNEWTON;
    }
    else if ( algorithm_str.compare("TNEWTON_RESTART" ) == 0 ) {
        algorithm = nlopt::LD_TNEWTON_RESTART;
    }
    else if ( algorithm_str.compare("TNEWTON_PRECOND" ) == 0 ) {
        algorithm = nlopt::LD_TNEWTON_PRECOND;
    }
    else if ( algorithm_str.compare("TNEWTON_PRECOND_RESTART" ) == 0 ) {
        algorithm = nlopt::LD_TNEWTON_PRECOND_RESTART;
    }
    else if ( algorithm_str.compare("MMA" ) == 0 ) {
        algorithm = nlopt::LD_MMA;
    }
    else if ( algorithm_str.compare("CCSAQ" ) == 0 ) {
        algorithm = nlopt::LD_CCSAQ;
    }
    else {
        // unknown algorithm code
        algorithm = nlopt::NUM_ALGORITHMS;       // Not an algorithm, so this should result in a runtime error.
    }

    return algorithm;
}

nlopt::opt initNLOPT(int n_param, List options) {

  // Prepare optimization by setting nlopt options
  nlopt::algorithm algo = getAlgorithmCode(as<std::string>(options["algorithm"])) ;
  nlopt::opt opt(algo, n_param);
  opt.set_xtol_rel(as<double>(options["xtol_rel"]));
  opt.set_ftol_abs(as<double>(options["ftol_abs"]));
  opt.set_ftol_rel(as<double>(options["ftol_rel"]));
  opt.set_maxeval (as<int>   (options["maxeval" ]));
  opt.set_xtol_abs(as<stdvec>(options["xtol_abs"]));
  opt.set_lower_bounds(as<stdvec>(options["lower_bound"]));

  return opt;
}
