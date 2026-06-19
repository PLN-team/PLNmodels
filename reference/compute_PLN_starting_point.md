# Helper function for PLN initialization.

Barebone function to compute starting points for B, M and S2 when
fitting a PLN. Mostly intended for internal use.

## Usage

``` r
compute_PLN_starting_point(Y, X, O, w, method = c("LM", "GLM"))
```

## Arguments

- Y:

  Response count matrix

- X:

  Covariate matrix. Note that initialization will fail if the model
  matrix is singular.

- O:

  Offset matrix (in log-scale)

- w:

  Weight vector (defaults to 1)

- method:

  character: strategy used to initialize B. Either `"LM"` (default, fast
  weighted log-linear regression) or `"GLM"` (p independent Poisson
  GLMs, more accurate for complex or unbalanced designs but slower).

## Value

a named list of starting values for model parameter B and variational
parameters M and S2 used in the iterative optimization algorithm of
[`PLN()`](https://pln-team.github.io/PLNmodels/reference/PLN.md)

## Details

- **B**: estimated by weighted LM (`method = "LM"`, default) or p
  independent Poisson GLMs (`method = "GLM"`). The GLM option gives
  better B estimates for factorial or unbalanced designs at the cost of
  p IRLS fits.

- **M**: initialized to `log((1 + Y) / exp(O))` (M_full in the X\*B +
  M_res parameterization).

- **S**: initialized element-wise to `1 / sqrt(2 + Y)`, the approximate
  VE-step optimum at Omega = I. This adapts automatically to count
  levels: high S for zero counts (high uncertainty), low S for large
  counts.

## Examples

``` r
if (FALSE) { # \dontrun{
data(barents)
Y <- barents$Abundance
X <- model.matrix(Abundance ~ Latitude + Longitude + Depth + Temperature, data = barents)
O <- log(barents$Offset)
w <- rep(1, nrow(Y))
compute_PLN_starting_point(Y, X, O, w)
compute_PLN_starting_point(Y, X, O, w, method = "GLM")
} # }
```
