# Prepare data for use in PLN models

Prepare data in proper format for use in PLN model and its variants. The
function (i) merges a count table and a covariate data frame in the most
comprehensive way and (ii) computes offsets from the count table using
one of several normalization schemes (TSS, CSS, RLE, GMPR, Wrench, etc).
The function fails with informative messages when the heuristics used
for sample matching fail.

## Usage

``` r
prepare_data(
  counts,
  covariates,
  offset = "TSS",
  call = rlang::caller_env(),
  ...
)
```

## Arguments

- counts:

  Required. An abundance count table, preferably with dimensions names
  and species as columns.

- covariates:

  Required. A covariates data frame, preferably with row names.

- offset:

  Optional. Normalization scheme used to compute scaling factors used as
  offset during PLN inference. Available schemes are "TSS" (Total Sum
  Scaling, default), "CSS" (Cumulative Sum Scaling, used in
  metagenomeSeq), "RLE" (Relative Log Expression, used in DESeq2),
  "GMPR" (Geometric Mean of Pairwise Ratio, introduced in Chen et al.,
  2018), Wrench (introduced in Kumar et al., 2018) or "none".
  Alternatively the user can supply its own vector or matrix of offsets
  (see note for specification of the user-supplied offsets).

- call:

  Optional. The execution environment in which to set the local error
  call.

- ...:

  Additional parameters passed on to
  [`compute_offset()`](https://pln-team.github.io/PLNmodels/reference/compute_offset.md)

## Value

A data.frame suited for use in
[`PLN()`](https://pln-team.github.io/PLNmodels/reference/PLN.md) and its
variants with two specials components: an abundance count matrix (in
component "Abundance") and an offset vector/matrix (in component
"Offset", only if offset is not set to "none")

## Note

User supplied offsets should be either vectors/column-matrices or have
the same number of column as the original count matrix and either (i)
dimension names or (ii) the same dimensions as the count matrix. Samples
are trimmed in exactly the same way to remove empty samples.

## References

Chen, L., Reeve, J., Zhang, L., Huang, S., Wang, X. and Chen, J. (2018)
GMPR: A robust normalization method for zero-inflated count data with
application to microbiome sequencing data. PeerJ, 6, e4600
[doi:10.7717/peerj.4600](https://doi.org/10.7717/peerj.4600)

Paulson, J. N., Colin Stine, O., Bravo, H. C. and Pop, M. (2013)
Differential abundance analysis for microbial marker-gene surveys.
Nature Methods, 10, 1200-1202
[doi:10.1038/nmeth.2658](https://doi.org/10.1038/nmeth.2658)

Anders, S. and Huber, W. (2010) Differential expression analysis for
sequence count data. Genome Biology, 11, R106
[doi:10.1186/gb-2010-11-10-r106](https://doi.org/10.1186/gb-2010-11-10-r106)

Kumar, M., Slud, E., Okrah, K. et al. (2018) Analysis and correction of
compositional bias in sparse sequencing count data. BMC Genomics 19, 799
[doi:10.1186/s12864-018-5160-5](https://doi.org/10.1186/s12864-018-5160-5)

Robinson, M.D., Oshlack, A. (2010) A scaling normalization method for
differential expression analysis of RNA-seq data. Genome Biol 11, R25
[doi:10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25)

## See also

[`compute_offset()`](https://pln-team.github.io/PLNmodels/reference/compute_offset.md)
for details on the different normalization schemes

## Examples

``` r
data(trichoptera)
proper_data <- prepare_data(
 counts     = trichoptera$Abundance,
 covariates = trichoptera$Covariate,
 offset     = "GMPR",
 scale      = "count"
)
proper_data$Abundance
#>    Che Hyc Hym Hys  Psy Aga Glo Ath Cea Ced Set All Han Hfo Hsp Hve Sta
#> 1    0   0   5   0   17   0   0   0   0   2   0   1   0   1   2   0   1
#> 2    0   0   3   0    8   0   0   0   0   0   0   2   0   0   0   0   0
#> 3    0   0   1   0   32   0   0   0   0   0   0   4   0   0   1   0   0
#> 4    0   0   3   0  176   4   0   0   0   1   0   3   0   0   3   0   2
#> 5    0   0   4   0   69   2   0   0   0   0   0   1   0   0   1   0   2
#> 6    0   0   2   0   14   1   0   0   0   0   0   1   0   0   0   0   0
#> 7    0   0   2   0    4   0   0   0   0   2   0   0   0   0   0   0   0
#> 8    0   0   1   0   29   1   0   0   0   0   0   0   0   0   0   0   3
#> 9    0   0   4   0    8   0   0   0   0   0   0   0   0   0   0   0   0
#> 10   0   0   2   0    2   0   0   0   0   0   0   0   0   0   0   0   0
#> 11   0   0   1   0    3   0   0   0   0   0   0   0   0   0   0   0   0
#> 12   0   0   0   0    3   0   0   0   0   0   0   0   0   0   0   0   0
#> 13   1   0   1   0   40   2   0   1   0   2   0   1   0   0   0   1   0
#> 14   0   0   4   0   20   2   1   0   0   1   0   0   0   0   1   0   4
#> 15   0   0   3   0  500  19   0   2   0  17   5   2   0   0   2   1  49
#> 16   0   0   2   0  142   9   0   0   0   3   3   2   0   0   0   0  11
#> 17   0   0   2   0   44   1   0   0   1   4   3   0   0   0   0   0   3
#> 18   0   0   3   1   31   4   0   0   0   1   7   0   0   0   0   0   4
#> 19   0   0   6   0   30   6   0   0   0   2   5   0   0   0   0   0   7
#> 20   0   1   7   0   71   5   1   1   0   5  18   0   0   1   5   0  12
#> 21   0   0   5   0   20   0   0   0   0   5   1   0   0   0   0   0   4
#> 22   0   0   3   0    5   0   0   0   0   2   0   0   0   0   2   0   1
#> 23   0   0   0   1    8   0   0   0   0   1   0   0   0   1   4   0   2
#> 24   0   0   1   0    2   0   0   0   0   0   0   0   0   0   0   0   0
#> 25   0   0   3   0   10   1   1   0   0   0   2   0   0   2   7   0   1
#> 26   0   0   1   0    9   0   0   0   0   1   2   1   3   3  14   0   6
#> 27   0   0   1   0   14   0   0   0   0   5   3   0   7   0  10   0   4
#> 28   0   0   0   0    2   0   0   1   0   0   0   0   0   1   3   0   1
#> 29   0   0   0   0    4   0   0   0   0   0   1   0   0   0   3   0   1
#> 30   0   1   0   0 1173   2   3   1   2   9  71   0  27  69 227   1  13
#> 31   0   0  12   0 2671   1   4   1   3   7  49   1  15  52 161   1   2
#> 32   0   0   7   0   33   4   0   0   0   0   0   3  33   1   0   0   7
#> 33   0   0  12   0   62   8   0   0   0   3   1  15  23   0   5   0   6
#> 34   0   0  15   0  220  20   0   0   0   2   0   5  27   0   5   0  33
#> 35   0   0   6   0   29   6   0   0   0   3   0   2   3   0   0   0  17
#> 36   0   0   4   0   60   3   0   0   0   0   0   0   2   0   0   0  21
#> 37   0   0  18   1   21   2   0   0   0   2   1   2   5   0   0   0  11
#> 38   0   0   4   0    3   0   0   0   0   4   0   0   1   0   0   0   3
#> 39   0   0   5   0    5   0   0   0   0   1   0   0   1   0   0   0   2
#> 40   0   0   1   0   18   0   0   0   0   0   0   1   0   0   0   0   0
#> 41   0   0   4   0   49   0   0   0   0   3   0   1   0   0   0   2  11
#> 42   0   0   3   1   33   0   0   2   0   1   0   0   2   0   0   0  11
#> 43   0   0   3   0   71   1   0   0   0  11   1   0   3   0   0   0   5
#> 44   0   0   6   1   28   1   2   1   0   3   0   1   0   0   0   0   0
#> 45   0   0   5   1   37   3   1   2   0   7   1   1   0   0   0   1   3
#> 46   1   0   4   0  103   1   1   0   0   2  10   2   0   2   3   2  18
#> 47   0   0   2   0   11   0   0   1   0   1   1   0   0   0   0   0   0
#> 48   1   0   0   0   17   0   0   1   1   2   0   0   3   0   0   0   6
#> 49   0   1   2   0   27   0   0   0   0   1   4   0  36   0  11   0   4
proper_data$Offset
#>  [1]  30.084856  29.463845  37.495950  57.505735  49.888839  29.513765
#>  [7]  19.949548  31.864374  32.075396  12.733865   9.104921   4.550832
#> [13]  41.915322  36.609494 249.546262 100.457779  46.401773  47.994943
#> [19]  69.219635 106.791751  48.521274  19.577977  18.898231   7.790568
#> [25]  23.345452  25.239568  34.703381   9.006205  10.116829 542.364607
#> [31] 243.692871  82.500856 108.703715 241.444947  77.274416  95.984236
#> [37]  62.017749  29.133040  22.133947  25.942671  75.746825  51.798296
#> [43]  67.640234  52.596146  57.108005  83.702490  22.685087  43.059821
#> [49]  41.384158
```
