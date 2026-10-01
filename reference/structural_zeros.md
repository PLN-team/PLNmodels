# Species absent from a whole level of a factor

Finds the species that are absent from all the samples of a level of a
factor covariate while present in the other samples, beyond what chance
alone would explain. Such structural zeros are what a PLN model without
the corresponding covariate cannot represent: it fits them by sending
the latent means to minus infinity, the latent variance of the species
blows up, and in a network the species ends up connected to most of the
others (see the field `degenerate_species` of a
[`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)).
[`prepare_data()`](https://pln-team.github.io/PLNmodels/reference/prepare_data.md)
runs this check and reports its result in a message.

## Usage

``` r
structural_zeros(counts, covariates, alpha = 0.05)
```

## Arguments

- counts:

  An abundance count table, with species as columns.

- covariates:

  A covariates data frame, with the samples of `counts` as rows. Only
  its factor, character and logical columns are used.

- alpha:

  Level of the test, after a Bonferroni correction over all the pairs of
  a species and a level that are tested. Default is `0.05`.

## Value

A data frame with one row per reported pair, sorted by p-value, with
columns `species`, `covariate`, `level`, `n_samples` (the number of
samples of the level), `prevalence_elsewhere` (the proportion of the
other samples where the species is present) and `p_value`
(Bonferroni-adjusted). It has no row when nothing is reported.

## Details

For a species present in `K` of the `N` samples, the probability that
none of the `n` samples of a level is among them, were the presences
spread at random, is hypergeometric: `phyper(0, K, N - K, n)`. A rare
species is easily absent from a level by chance, and is not reported; a
species present in most of the other samples is.

## See also

[`prepare_data()`](https://pln-team.github.io/PLNmodels/reference/prepare_data.md)

## Examples

``` r
data(oaks)
## species that are absent from a whole tree, or from a whole type of branch
structural_zeros(oaks$Abundance, oaks[c("tree", "branch")])
#>       species covariate        level n_samples prevalence_elsewhere
#> 1    f_OTU_63      tree  susceptible        39            1.0000000
#> 2    f_OTU_30      tree  susceptible        39            0.9740260
#> 3    f_OTU_65      tree    resistant        39            0.7142857
#> 4  f_OTU_1011      tree  susceptible        39            0.4935065
#> 5  f_OTU_1011      tree    resistant        39            0.4935065
#> 6  f_OTU_1090      tree  susceptible        39            0.4935065
#> 7    f_OTU_46      tree  susceptible        39            0.4935065
#> 8    f_OTU_46      tree    resistant        39            0.4935065
#> 9  f_OTU_1090      tree intermediate        38            0.4871795
#> 10  f_OTU_579      tree    resistant        39            0.4805195
#> 11  f_OTU_579      tree intermediate        38            0.4743590
#> 12   f_OTU_32      tree    resistant        39            0.4155844
#>         p_value
#> 1  2.984898e-29
#> 2  2.447616e-26
#> 3  6.608189e-13
#> 4  8.124001e-07
#> 5  8.124001e-07
#> 6  8.124001e-07
#> 7  8.124001e-07
#> 8  8.124001e-07
#> 9  1.584180e-06
#> 10 1.604490e-06
#> 11 3.052445e-06
#> 12 4.054306e-05
```
