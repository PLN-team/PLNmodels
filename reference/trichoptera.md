# Trichoptera data set

Data gathered between 1959 and 1960 during 49 insect trapping nights.
For each trapping night, the abundance of 17 Trichoptera species is
recorded as well as 6 meteorological variables which may influence the
abundance of each species. Finally, the observations (that is to say,
the trapping nights), have been classified into 12 groups corresponding
to contiguous nights between summer 1959 and summer 1960.

## Usage

``` r
trichoptera
```

## Format

A list with 2 two data frames:

- Abundance:

  a 49 x 17 matrix of abundancies/counts (49 trapping nights and 17
  trichoptera species)

- Covariate:

  a 49 x 7 data frame of covariates:

  Temperature

  :   Evening Temperature in Celsius

  Wind

  :   Wind in m/s

  Pressure

  :   Pressure in mm Hg

  Humidity

  :   relative to evening humidity in percent

  Cloudiness

  :   proportion of sky coverage at 9pm

  Precipitation

  :   Nighttime precipitation in mm

  Group

  :   a factor of 12 levels for the definition of the consecutive night
      groups

In order to prepare the data for using formula in multivariate analysis
(multiple outputs and inputs), use
[`prepare_data()`](https://pln-team.github.io/PLNmodels/reference/prepare_data.md).
We only kept a subset of the original meteorological covariates for
illustration purposes.

## Source

Data from P. Usseglio-Polatera.

## References

Usseglio-Polatera, P. and Auda, Y. (1987) Influence des facteurs
météorologiques sur les résultats de piégeage lumineux. Annales de
Limnologie, 23, 65–79. (code des espèces p. 76) See a data description
at <http://pbil.univ-lyon1.fr/R/pdf/pps034.pdf> (in French)

## See also

[`prepare_data()`](https://pln-team.github.io/PLNmodels/reference/prepare_data.md)

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
```
