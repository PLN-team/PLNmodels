#include "utils.h"

using namespace Rcpp;

arma::mat logfact(arma::mat Y) {
  arma::mat v = Y.replace(0, 1);
  return sum(v % arma::log(v) - v + arma::log(8*pow(v,3) + 4*pow(v, 2) + v + 1/30)/6 + std::log(M_PI)/2, 1);
}

// Convert string to nlopt_algorithm
//
// restrict the choices to algorithms meaningful for PLN optimization
nlopt_algorithm getAlgorithmCode( const std::string & algorithm_str) {

    nlopt_algorithm algorithm;

    if ( algorithm_str.compare("LBFGS_NOCEDAL") == 0 ) {
        algorithm = NLOPT_LD_LBFGS_NOCEDAL;
    }
    else if ( algorithm_str.compare("LBFGS" ) == 0 ) {
        algorithm = NLOPT_LD_LBFGS;
    }
    else if ( algorithm_str.compare("VAR1" ) == 0 ) {
        algorithm = NLOPT_LD_VAR1;
    }
    else if ( algorithm_str.compare("VAR2" ) == 0 ) {
        algorithm = NLOPT_LD_VAR2;
    }
    else if ( algorithm_str.compare("TNEWTON" ) == 0 ) {
        algorithm = NLOPT_LD_TNEWTON;
    }
    else if ( algorithm_str.compare("TNEWTON_RESTART" ) == 0 ) {
        algorithm = NLOPT_LD_TNEWTON_RESTART;
    }
    else if ( algorithm_str.compare("TNEWTON_PRECOND" ) == 0 ) {
        algorithm = NLOPT_LD_TNEWTON_PRECOND;
    }
    else if ( algorithm_str.compare("TNEWTON_PRECOND_RESTART" ) == 0 ) {
        algorithm = NLOPT_LD_TNEWTON_PRECOND_RESTART;
    }
    else if ( algorithm_str.compare("MMA" ) == 0 ) {
        algorithm = NLOPT_LD_MMA;
    }
    else if ( algorithm_str.compare("CCSAQ" ) == 0 ) {
        algorithm = NLOPT_LD_CCSAQ;
    }
    else {
        // unknown algorithm code
        algorithm = NLOPT_NUM_ALGORITHMS;       // Not an algorithm, so this should result in a runtime error.
    }

    return algorithm;
}

nlopt_opt initNLOPT(int n_param, List options) {

  // Prepare optimization by setting nlopt options

  nlopt_algorithm algo = getAlgorithmCode(as<std::string>(options["algorithm"])) ;
  nlopt_opt opt = nlopt_create(algo, n_param);

  stdvec xtol_abs    = as<stdvec>(options["xtol_abs"   ]) ;
  stdvec lower_bound = as<stdvec>(options["lower_bound"]) ;

  nlopt_set_xtol_rel    (opt, as<double>(options["xtol_rel"]));
  nlopt_set_ftol_abs    (opt, as<double>(options["ftol_abs"]));
  nlopt_set_ftol_rel    (opt, as<double>(options["ftol_rel"]));
  nlopt_set_maxeval     (opt, as<int>   (options["maxeval" ]));
  nlopt_set_maxtime     (opt, as<double>(options["maxtime" ]));
  nlopt_set_xtol_abs    (opt, &xtol_abs[0]   );
  nlopt_set_lower_bounds(opt, &lower_bound[0]) ;

  return opt;
}
