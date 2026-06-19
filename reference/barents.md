# Barents fish data set

This data set gives the abundance of 30 fish species observed in 89
sites in the Barents sea. For each site, 4 additional covariates are
known. Subsample of the original datasets studied by Fossheim et al,
2006.

## Usage

``` r
barents
```

## Format

A data frame with 6 variables:

- Abundance: A 30 fish species by 89 sites count matrix

- Offset: A 30 fish species by 89 samples offset matrix, measuring the
  sampling effort in each site

- 4 covariates for latitude, longitude, depth (in meters), temperature
  (in Celsius degrees).

## Source

Data from M. Fossheim and coauthors.

## References

Fossheim, Maria, Einar M. Nilssen, and Michaela Aschan. "Fish
assemblages in the Barents Sea." Marine Biology Research 2.4 (2006).
[doi:10.1080/17451000600815698](https://doi.org/10.1080/17451000600815698)

## Examples

``` r
data(barents)
```
