# Compute offsets from a count data using one of several normalization schemes

Computes offsets from the count table using one of several normalization
schemes (TSS, CSS, RLE, GMPR, Wrench, TMM, etc) described in the
literature.

## Usage

``` r
compute_offset(
  counts,
  offset = c("TSS", "GMPR", "RLE", "CSS", "Wrench", "TMM", "none"),
  scale = c("none", "count"),
  ...
)
```

## Arguments

- counts:

  Required. An abundance count table, preferably with dimensions names
  and species as columns.

- offset:

  Optional. Normalization scheme used to compute scaling factors used as
  offset during PLN inference. Available schemes are "TSS" (Total Sum
  Scaling, default), "CSS" (Cumulative Sum Scaling, used in
  metagenomeSeq), "RLE" (Relative Log Expression, used in DESeq2),
  "GMPR" (Geometric Mean of Pairwise Ratio, introduced in Chen et al.,
  2018), Wrench (introduced in Kumar et al., 2018) or "none".
  Alternatively the user can supply its own vector or matrix of offsets
  (see note for specification of the user-supplied offsets).

- scale:

  Either `"none"` (default) or `"count"`. Should the offset be
  normalized to be on the same scale as the counts ?

- ...:

  Additional parameters passed on to specific methods (for now CSS and
  RLE)

## Value

If `offset = "none"`, `NULL` else a vector of length `nrow(counts)` with
one offset per sample.

## Details

RLE has additional `pseudocounts` and `type` arguments to add
pseudocounts to the observed counts (defaults to 0L) and to compute
offsets using only positive counts (if `type == "poscounts"`). This
mimics the behavior of `DESeq2::DESeq()` when using
`sfType == "poscounts"`. CSS has an additional `reference` argument to
choose the location function used to compute the reference quantiles
(defaults to `median` as in the Nature publication but can be set to
`mean` to reproduce behavior of functions cumNormStat\* from
metagenomeSeq). Wrench has two additional parameters: `groups` to
specify sample groups and `type` to either reproduce exactly the default
`Wrench::wrench()` behavior (`type = "wrench"`, default) or to use
simpler heuristics (`type = "simple"`). Note that (i) CSS normalization
fails when the median absolute deviation around quantiles does not
become instable for high quantiles (limited count variations both within
and across samples) and/or one sample has less than two positive counts,
(ii) RLE fails when there are no common species across all samples
(unless `type == "poscounts"` has been specified) and (iii) GMPR fails
if a sample does not share any species with all other samples. TMM code
between two libraries is simplified and adapted from M. Robinson
(edgeR:::.calcFactorTMM). The final output is however different from the
one produced by edgeR:::.calcFactorTMM as they are intended to be used
as such in the model (whereas they need to be multiplied by sequencing
depths in edgeR)

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

## Examples

``` r
data(trichoptera)
counts <- trichoptera$Abundance
compute_offset(counts)
#>    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
#>   29   13   38  192   79   18    8   34   12    4    4    3   49   33  600  172 
#>   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
#>   58   51   56  127   35   13   17    3   27   40   44    8    9 1599 2980   88 
#>   33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48 
#>  135  327   66   90   63   15   14   20   70   53   95   43   62  149   16   31 
#>   49 
#>   86 
## Other normalization schemes
compute_offset(counts, offset = "RLE", pseudocounts = 1)
#>         1         2         3         4         5         6         7         8 
#> 0.9186270 0.8349121 0.8570257 0.9186270 0.9186270 0.8349121 0.8192245 0.8570257 
#>         9        10        11        12        13        14        15        16 
#> 0.7322797 0.7322797 0.6321923 0.6321923 0.9240361 0.9186270 1.3788037 0.9240361 
#>        17        18        19        20        21        22        23        24 
#> 0.9186270 0.9240361 0.9240361 1.7140514 0.8908577 0.8570257 0.8349121 0.6321923 
#>        25        26        27        28        29        30        31        32 
#> 0.9240361 0.9240361 0.9186270 0.8570257 0.8349121 2.7721084 3.2934218 0.9584503 
#>        33        34        35        36        37        38        39        40 
#> 1.0406547 0.9584503 0.9584503 0.9186270 0.9584503 0.8908577 0.8349121 0.7322797 
#>        41        42        43        44        45        46        47        48 
#> 0.9240361 0.9240361 0.9584503 0.9584503 1.2643846 1.7140514 0.8233555 0.8908577 
#>        49 
#> 0.9186270 
compute_offset(counts, offset = "Wrench", groups = trichoptera$Covariate$Group)
#>  [1]  0.69617112  0.17772439  0.31918810  0.58503048  0.37114188  0.19846015
#>  [7]  0.05871729  0.17732815  0.06477092  0.03038695  0.02636997  0.01335544
#> [13] 15.97878650  1.25308073  2.59166652  0.73153479  8.31955593  3.29726449
#> [19]  0.68605117 35.74486487  0.29922507  0.20183174  3.21536182  0.02964492
#> [25]  2.42898846  2.11066579  0.95438702  1.17336576  0.22739391 70.13928757
#> [31] 31.14789183  2.24096704  2.14256759  1.74908073  0.84165887  0.55533980
#> [37]  3.59125649  0.22105352  0.20919575  0.16386981  2.70312687  3.73203456
#> [43]  0.53690902  5.31410880  6.53316557 23.33458085  0.56217888 21.34169565
#> [49] 34.52859945
compute_offset(counts, offset = "GMPR")
#>  [1]  0.7981233  0.7816485  0.9947327  1.5255738  1.3235046  0.7829728
#>  [7]  0.5292430  0.8453323  0.8509305  0.3378176  0.2415451  0.1207293
#> [13]  1.1119746  0.9712159  6.6202307  2.6650516  1.2309960  1.2732613
#> [19]  1.8363327  2.8330860  1.2872244  0.5193856  0.5013525  0.2066765
#> [25]  0.6193332  0.6695823  0.9206485  0.2389263  0.2683901 14.3884296
#> [31]  6.4649457  2.1886711  2.8838087  6.4053103  2.0500185  2.5463727
#> [37]  1.6452733  0.7728725  0.5871931  0.6882350  2.0094930  1.3741607
#> [43]  1.7944326  1.3953269  1.5150224  2.2205494  0.6018143  1.1423371
#> [49]  1.0978833
compute_offset(counts, offset = "TMM")
#>         1         2         3         4         5         6         7         8 
#> 0.8267930 0.4195206 1.6780823 7.6180131 3.2441013 0.7341610 0.2615234 1.4001773 
#>         9        10        11        12        13        14        15        16 
#> 0.4195206 0.1048801 0.1573202 0.1573202 0.8914812 0.9609351 7.3465136 5.9940572 
#>        17        18        19        20        21        22        23        24 
#> 1.4845609 1.4328717 1.4574326 1.8815971 1.0316007 0.2885066 0.4022524 0.1048801 
#>        25        26        27        28        29        30        31        32 
#> 0.4768149 0.4719606 0.8597281 0.1525071 0.1999304 3.3028204 3.2417791 1.9540165 
#>        33        34        35        36        37        38        39        40 
#> 3.1306172 5.7554868 1.5100080 2.9737473 1.1372035 0.4035513 0.2971604 0.9439213 
#>        41        42        43        44        45        46        47        48 
#> 2.3690225 1.6450422 3.4859969 1.4482832 1.9361593 2.6744436 0.5781787 0.8914812 
#>        49 
#> 1.2731958 
## User supplied offsets
my_offset <- setNames(rep(1, nrow(counts)), rownames(counts))
compute_offset(counts, offset = my_offset)
#>    Che Hyc Hym Hys Psy Aga Glo Ath Cea Ced Set All Han Hfo Hsp Hve Sta
#> 1    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 2    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 3    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 4    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 5    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 6    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 7    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 8    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 9    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 10   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 11   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 12   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 13   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 14   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 15   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 16   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 17   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 18   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 19   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 20   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 21   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 22   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 23   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 24   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 25   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 26   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 27   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 28   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 29   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 30   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 31   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 32   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 33   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 34   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 35   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 36   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 37   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 38   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 39   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 40   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 41   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 42   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 43   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 44   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 45   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 46   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 47   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 48   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
#> 49   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
```
