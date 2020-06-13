#include "nlopt_utils.h"

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

nlopt_opt initNLOPT(int n_param, Rcpp::List options) {

  // Prepare optimization by setting nlopt options

  nlopt_algorithm algo = getAlgorithmCode(Rcpp::as<std::string>(options["algorithm"])) ;
  nlopt_opt opt = nlopt_create(algo, n_param);

  stdvec xtol_abs    = Rcpp::as<stdvec>(options["xtol_abs"   ]) ;

  nlopt_set_xtol_rel    (opt, Rcpp::as<double>(options["xtol_rel"]));
  nlopt_set_ftol_abs    (opt, Rcpp::as<double>(options["ftol_abs"]));
  nlopt_set_ftol_rel    (opt, Rcpp::as<double>(options["ftol_rel"]));
  nlopt_set_maxeval     (opt, Rcpp::as<int>   (options["maxeval" ]));
  nlopt_set_maxtime     (opt, Rcpp::as<double>(options["maxtime" ]));
  nlopt_set_xtol_abs    (opt, &xtol_abs[0]   );

  return opt;
}
