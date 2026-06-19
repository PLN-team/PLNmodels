# Helper function for ZIPLN initialization.

Fast LM-based starting point for ZIPLN: one multivariate `lm.fit` for
the PLN component and empirical zero rates / binomial GLMs for the ZI
component. Replaces the previous per-species
[`pscl::zeroinfl`](https://rdrr.io/pkg/pscl/man/zeroinfl.html) loop.

## Usage

``` r
compute_ZIPLN_starting_point(Y, X, X0, O, w = NULL)
```

## Arguments

- Y:

  Response count matrix (n × p)

- X:

  Design matrix for the PLN component (n × d)

- X0:

  Design matrix for the ZI component (n × d0, empty `matrix(NA,0,0)`
  when unused)

- O:

  Offset matrix in log-scale (n × p)

- w:

  Weight vector of length n (defaults to uniform weights)

## Value

Named list: `B` (d × p), `M` (n × p), `S2` (n × p), `R` (n × p), `B0`
(d0 × p)
