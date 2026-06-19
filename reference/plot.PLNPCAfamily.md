# Display the criteria associated with a collection of PLNPCA fits (a PLNPCAfamily)

Display the criteria associated with a collection of PLNPCA fits (a
PLNPCAfamily)

## Usage

``` r
# S3 method for class 'PLNPCAfamily'
plot(x, criteria = c("loglik", "BIC", "ICL"), reverse = FALSE, ...)
```

## Arguments

- x:

  an R6 object with class
  [`PLNPCAfamily`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfamily.md)

- criteria:

  vector of characters. The criteria to plot in c("loglik", "BIC",
  "ICL"). Default is c("loglik", "BIC", "ICL").

- reverse:

  A logical indicating whether to plot the value of the criteria in the
  "natural" direction (loglik - 0.5 penalty) or in the "reverse"
  direction (-2 loglik + penalty). Default to FALSE, i.e use the natural
  direction, on the same scale as the log-likelihood.

- ...:

  additional parameters for S3 compatibility. Not used

## Value

Produces a plot representing the evolution of the criteria of the
different models considered, highlighting the best model in terms of BIC
and ICL (see details).

## Details

The BIC and ICL criteria have the form 'loglik - 1/2 \* penalty' so that
they are on the same scale as the model log-likelihood. You can change
this direction and use the alternate form '-2\*loglik + penalty', as
some authors do, by setting `reverse = TRUE`.

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPCAs <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, ranks = 1:5)
#> 
#>  Initialization...
#> 
#>  Adjusting 5 PLN models for PCA analysis.
#>   Rank approximation = 1      Rank approximation = 2      Rank approximation = 3      Rank approximation = 4      Rank approximation = 5 
#>  Post-treatments
#>  DONE!
if (FALSE) { # \dontrun{
plot(myPCAs)
} # }
```
