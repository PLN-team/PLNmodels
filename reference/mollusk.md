# Mollusk data set

This data set gives the abundance of 32 mollusk species in 163 samples.
For each sample, 4 additional covariates are known.

## Usage

``` r
mollusk
```

## Format

A list with 2 two data frames:

- Abundance:

  a 163 x 32 data frame of abundancies/counts (163 samples and 32
  mollusk species)

- Covariate:

  a 163 x 4 data frame of covariates:

  site

  :   a factor with 8 levels indicating the sampling site

  season

  :   a factor with 4 levels indicating the season

  method

  :   a factor with 2 levels for the method of sampling - wood or string

  duration

  :   a numeric with 3 levels for the time of exposure in week

In order to prepare the data for using formula in multivariate analysis
(multiple outputs and inputs), use
[`prepare_data()`](https://pln-team.github.io/PLNmodels/reference/prepare_data.md).
Original data set has been extracted from ade4.

## Source

Data from Richardot-Coulet, Chessel and Bournaud.

## References

Richardot-Coulet, M., Chessel D. and Bournaud M. (1986) Typological
value of the benthos of old beds of a large river. Methodological
approach. Archiv fùr Hydrobiologie, 107, 363–383.

## See also

[`prepare_data()`](https://pln-team.github.io/PLNmodels/reference/prepare_data.md)

## Examples

``` r
data(mollusk)
mollusc <- prepare_data(mollusk$Abundance, mollusk$Covariate)
#> Warning: ! There is at least one empty sample in `counts`.
#> ℹ <4> samples (<134/137/145/146>) in `counts` have been dropped for lack of
#>   positive counts.
```
