{ rPackages, ... }:
let
  r-deps = with rPackages; [
    corrplot
    dplyr
    future
    future_apply
    ggplot2
    glassoFast
    gridExtra
    igraph
    magrittr
    MASS
    Matrix
    nloptr
    purrr
    R6
    Rcpp
    RcppArmadillo
    rlang
    tidyr
    torch
  ];
in
rPackages.buildRPackage {
  name = "PLNmodels";
  src = ./.;
  propagatedBuildInputs = r-deps;
}

