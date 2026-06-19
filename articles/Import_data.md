# Data importation in PLNmodels

## Preliminaries

This vignette documents the data format used in **PLNmodel** by `PLN`
and its variants. It also shows how to create an object in the proper
format for further analyses from (i) tabular data, (ii) biom-class
objects and (iii) phyloseq-class objects.

## Format description

We illustrate the format using trichoptera data set, a full description
of which can be found in [the corresponding
vignette](https://pln-team.github.io/PLNmodels/articles/Trichoptera.md).

``` r

library(PLNmodels)
data(trichoptera)
```

The trichoptera data set is a list made of two data frames: `Abundance`
(hereafter referred to as the *counts*) and `Covariate` (hereafter the
*covariates*).

``` r

str(trichoptera, max.level = 1)
```

    ## List of 2
    ##  $ Abundance:'data.frame':   49 obs. of  17 variables:
    ##  $ Covariate:'data.frame':   49 obs. of  7 variables:

The covariates include, among others, the wind, pressure and humidity.

``` r

names(trichoptera$Covariate)
```

    ## [1] "Temperature"   "Wind"          "Pressure"      "Humidity"     
    ## [5] "Cloudiness"    "Precipitation" "Group"

In the PLN framework, we model the counts from the covariates, let’s say
wind and pressure, using a Poisson Log-Normal model. Most models in R
use the so-called *formula interface* and it would thus be natural to
write something like

``` r

PLN(Abundance ~ Wind + Pressure, data = trichoptera)
```

Unfortunately and unlike many generalized linear models, the response in
PLN is intrinsically **multivariate**: it has 17 dimensions in our
example. The left hand side (LHS) must encode a multivariate response
across multiple samples, using a 2D-array (e.g. a matrix or a data
frame).

We must therefore prepare a data structure where `Abundance` refers to a
count *matrix* whereas `Wind` and `Pressure` refer to *vectors* before
feeding it to `PLN`. That’s the purpose of `prepare_data`.

``` r

trichoptera2 <- prepare_data(counts     = trichoptera$Abundance, 
                             covariates = trichoptera$Covariate)
str(trichoptera2)
```

    ## 'data.frame':    49 obs. of  9 variables:
    ##  $ Abundance    : num [1:49, 1:17] 0 0 0 0 0 0 0 0 0 0 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:49] "1" "2" "3" "4" ...
    ##   .. ..$ : chr [1:17] "Che" "Hyc" "Hym" "Hys" ...
    ##  $ Temperature  : num  18.7 19.8 22 23 22.5 23.9 15 17.2 15.4 14.1 ...
    ##  $ Wind         : num  -2.3 -2.7 -0.7 2.3 2.3 -2 -4.7 -1 -2.7 -3.7 ...
    ##  $ Pressure     : num  998 1000 997 991 990 ...
    ##  $ Humidity     : num  60 63 73 71 62 64 93 84 88 75 ...
    ##  $ Cloudiness   : num  19 0 6 81 50 50 100 19 69 6 ...
    ##  $ Precipitation: num  0 0 0 0 0 0 1.6 0 1.6 0 ...
    ##  $ Group        : Factor w/ 12 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Offset       : num  29 13 38 192 79 18 8 34 12 4 ...

If you look carefully, you can notice a few difference between
`trichoptera` and `trichoptera2`:

- the first is a `list` whereas the second is a `data.frame`[^1];
- `Abundance` is a matrix-column of `trichoptera2` that you can extract
  using the usual functions `[` and `[[` to retrieve the count matrix;
- `trichoptera2` has an additional `Offset` column (more on that later).

## Computing offsets

It is common practice when modeling count data to introduce an offset
term to control for different sampling efforts, exposures, baselines,
etc. The *proper way* to compute sample-specific offsets in still
debated and may vary depending on the field. There are nevertheless a
few popular methods:

- Total Sum Scaling (TSS), where the offset of a sample is the total
  count in that sample
- Cumulative Sum Scaling (CSS), introduced in ([Paulson et al.
  2013](#ref-CSS)), where the offset of a sample if the cumulative sum
  of counts in that sample, up to a quantile determined in a data driven
  way.
- Relative Log-Expression (RLE), implemented in ([Anders and Huber
  2010](#ref-DESeq2)), where all samples are used to compute a reference
  sample, each sample is compared to the reference sample using
  log-ratios and the offset is the median log-ratio.
- Geometric Mean of Pairwise Ratio (GMPR), introduced in ([Chen et al.
  2018](#ref-GMPR)) where each sample is compared to each other to
  compute a median log-ratio and the offset of a sample is the geometric
  means of those pairwise ratios.
- Wrench, introduced in ([Kumar et al. 2018](#ref-Kumar2018)) and fully
  implemented in the [Wrench
  package](https://bioconductor.org/packages/release/bioc/html/Wrench.html),
  where all samples are used to compute reference proportions and each
  sample is compared to the reference using ratios (and **not
  log-ratios**) of proportions to compute compositional correction
  factors. In that case, the offset is the product of (geometrically
  centered) compositional factors and (geometrically centered) depths.

Each of these offset be computed from a counts matrix using the
`compute_offset` function and changing its `offset` argument:

``` r

## same as compute_offset(trichoptera$Abundance, offset = "TSS")
compute_offset(trichoptera$Abundance) 
```

    ##    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
    ##   29   13   38  192   79   18    8   34   12    4    4    3   49   33  600  172 
    ##   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
    ##   58   51   56  127   35   13   17    3   27   40   44    8    9 1599 2980   88 
    ##   33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48 
    ##  135  327   66   90   63   15   14   20   70   53   95   43   62  149   16   31 
    ##   49 
    ##   86

In this particular example, the counts are too sparse and sophisticated
offset methods all fail (numeric output hidden)

``` r

compute_offset(trichoptera$Abundance, "CSS")
```

    ## Warning in offset_function(counts, ...): Some samples only have 1 positive
    ## values. Can't compute quantiles and fall back to TSS normalization

``` r

compute_offset(trichoptera$Abundance, "RLE")
```

    ## Warning in offset_function(counts, ...): Because of high sparsity, some samples
    ## have null or infinite offset.

``` r

compute_offset(trichoptera$Abundance, "GMPR")
```

We can mitigate this problem for the RLE offset by adding pseudocounts
to the counts although doing so has its own drawbacks.

``` r

compute_offset(trichoptera$Abundance, "RLE", pseudocounts = 1)
```

    ##         1         2         3         4         5         6         7         8 
    ## 0.9186270 0.8349121 0.8570257 0.9186270 0.9186270 0.8349121 0.8192245 0.8570257 
    ##         9        10        11        12        13        14        15        16 
    ## 0.7322797 0.7322797 0.6321923 0.6321923 0.9240361 0.9186270 1.3788037 0.9240361 
    ##        17        18        19        20        21        22        23        24 
    ## 0.9186270 0.9240361 0.9240361 1.7140514 0.8908577 0.8570257 0.8349121 0.6321923 
    ##        25        26        27        28        29        30        31        32 
    ## 0.9240361 0.9240361 0.9186270 0.8570257 0.8349121 2.7721084 3.2934218 0.9584503 
    ##        33        34        35        36        37        38        39        40 
    ## 1.0406547 0.9584503 0.9584503 0.9186270 0.9584503 0.8908577 0.8349121 0.7322797 
    ##        41        42        43        44        45        46        47        48 
    ## 0.9240361 0.9240361 0.9584503 0.9584503 1.2643846 1.7140514 0.8233555 0.8908577 
    ##        49 
    ## 0.9186270

A better solution consists in using only positive counts to compute the
offsets:

``` r

compute_offset(trichoptera$Abundance, "RLE", type = "poscounts")
```

    ##         1         2         3         4         5         6         7         8 
    ## 0.5631099 0.9462046 0.8299806 0.9462046 0.6460415 0.5756774 0.6308031 0.4947789 
    ##         9        10        11        12        13        14        15        16 
    ## 0.7988730 0.3574190 0.2207270 0.1260525 0.8051707 0.7289732 1.5591914 1.1850090 
    ##        17        18        19        20        21        22        23        24 
    ## 0.6389431 1.0000000 1.4606879 1.8224330 0.8403498 0.3942114 0.4058586 0.1997183 
    ##        25        26        27        28        29        30        31        32 
    ## 0.6286560 0.5631099 0.8777257 0.3440303 0.2504662 2.7086423 2.5800932 1.4584997 
    ##        33        34        35        36        37        38        39        40 
    ## 2.6050843 4.5898981 1.2185072 1.2616062 0.8823673 0.6250713 0.3950030 0.5631099 
    ##        41        42        43        44        45        46        47        48 
    ## 1.4511384 1.0000000 0.9462046 1.0882448 1.0000000 1.0631099 0.4621924 0.7900060 
    ##        49 
    ## 1.0672361

Finally, we can use wrench to compute the offsets:

``` r

compute_offset(trichoptera$Abundance, "Wrench")
```

    ##  [1]  0.41269451  0.17385897  0.31925785  0.70682391  0.44749382  0.20676142
    ##  [7]  0.10226641  0.21204241  0.09228293  0.03919178  0.03946910  0.02295077
    ## [13] 24.18096281  1.54038084  3.95187144  0.80172207 10.52136472  4.55854738
    ## [19]  0.46968859 51.88528943  0.32803605  0.14764973  4.13635353  0.03152826
    ## [25]  1.63626578  0.53769290  0.38482355  0.62689594  0.08030685 83.87692476
    ## [31] 33.40810265  0.75508731  1.05809625  1.18385013  0.58657751  0.43005272
    ## [37]  4.79925224  0.18381986  0.16078594  0.18490773  4.25999137  5.46689591
    ## [43]  0.60139850  7.56594658  9.89766712 33.73727074  0.74264459 29.13747796
    ## [49] 47.12017975

**Note** TSS is the only methods that produces offset on the same scale
as the counts, all others produces offsets that are (hopefully)
*proportional* to library sizes but on a different scale. To force the
offsets to be on the same scale as the counts for all methods, you can
use the option `scale = "count"`.

``` r

compute_offset(trichoptera$Abundance, "Wrench", scale = "count")
```

    ##            1            2            3            4            5            6 
    ##   17.2788462    7.2791915   13.3668057   29.5935647   18.7358369    8.6567635 
    ##            7            8            9           10           11           12 
    ##    4.2817277    8.8778701    3.8637358    1.6408958    1.6525069    0.9609112 
    ##           13           14           15           16           17           18 
    ## 1012.4174950   64.4932346  165.4584151   33.5667961  440.5123899  190.8589478 
    ##           19           20           21           22           23           24 
    ##   19.6650957 2172.3524890   13.7343345    6.1818537  173.1823794    1.3200369 
    ##           25           26           27           28           29           30 
    ##   68.5077809   22.5123254   16.1119346   26.2471113    3.3623170 3511.7901098 
    ##           31           32           33           34           35           36 
    ## 1398.7427987   31.6142749   44.3007651   49.5658750   24.5590443   18.0056062 
    ##           37           38           39           40           41           42 
    ##  200.9368679    7.6962380    6.7318451    7.7417851  178.3588942  228.8900199 
    ##           43           44           45           46           47           48 
    ##   25.1795750  316.7738500  414.3991879 1412.5245302   31.0933184 1219.9387043 
    ##           49 
    ## 1972.8451146

## Building data frame using `prepare_data`

We’ll already learned that `prepare_data` can join counts and covariates
into a single data.frame. It can also compute offset through
`compute_offset` and does so by default with `offset = "TSS"`, hence the
`Offset` column in `trichoptera2`. You can change the offset method and
provide additional arguments that will passed on to `compute_offset`.

``` r

str(prepare_data(trichoptera$Abundance, 
             trichoptera$Covariate, 
             offset = "RLE", pseudocounts = 1))
```

    ## 'data.frame':    49 obs. of  9 variables:
    ##  $ Abundance    : num [1:49, 1:17] 0 0 0 0 0 0 0 0 0 0 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:49] "1" "2" "3" "4" ...
    ##   .. ..$ : chr [1:17] "Che" "Hyc" "Hym" "Hys" ...
    ##  $ Temperature  : num  18.7 19.8 22 23 22.5 23.9 15 17.2 15.4 14.1 ...
    ##  $ Wind         : num  -2.3 -2.7 -0.7 2.3 2.3 -2 -4.7 -1 -2.7 -3.7 ...
    ##  $ Pressure     : num  998 1000 997 991 990 ...
    ##  $ Humidity     : num  60 63 73 71 62 64 93 84 88 75 ...
    ##  $ Cloudiness   : num  19 0 6 81 50 50 100 19 69 6 ...
    ##  $ Precipitation: num  0 0 0 0 0 0 1.6 0 1.6 0 ...
    ##  $ Group        : Factor w/ 12 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Offset       : num  0.919 0.835 0.857 0.919 0.919 ...

Different communities use different standard for the count data where
samples are either or columns of the counts matrix. `prepare_data` uses
heuristics to guess the direction of the counts matrix (or fail
informatively doing so) and automatically transpose it if needed.

Finally, `prepare_data` enforces sample-consistency between the counts
and the covariates and automatically trims away: - samples for which
only covariates or only counts are available; - samples with no positive
counts

For example, if we remove the first sample from the counts and the last
one from the covariates, we end up with 49 - 2 = 47 samples left, as
expected.

``` r

nrow(prepare_data(trichoptera$Abundance[-1, ], ## remove first sample
                  trichoptera$Covariate[-49,]  ## remove last sample
                  ))
```

    ## [1] 47

## Importing data from biom and phyloseq objects using `prepare_data_from_[phyloseq|biom]`

Community composition data are quite popular in microbial ecology and
usually stored in flat files using the [biom
format](http://biom-format.org/) and/or imported in R as phyloseq-class
objects ([McMurdie 2013](#ref-phyloseq)) using the Bioconductor
[phyloseq](https://joey711.github.io/phyloseq/) package.

We show here how to import data from a biom file (or biom-class object)
and form a phyloseq-class object.

### Reading from a biom file

Reading from a biom file requires the bioconductor package
[biomformat](https://www.bioconductor.org/packages/release/bioc/html/biomformat.html).
This package is **not** a standard dependency of PLNmodels and needs to
be installed separately.

You can easily prepare your data from a biom file using the following
steps:

- read your biom file with `biomformat::read_biom()`
- extract the count table with `biomformat::biom_data()`
- extract the covariates with `biomformat::sample_metadata()` (or build
  your own)
- feed them to `prepare_data`

as illustrated below:

``` r

## If biomformat is not installed, uncomment the following lines
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("biomformat")
library(biomformat)
biomfile <- system.file("extdata", "rich_dense_otu_table.biom", package = "biomformat")
biom <- biomformat::read_biom(biomfile)
## extract counts
counts <- as(biomformat::biom_data(biom), "matrix")
## extract covariates (or prepare your own)
covariates <- biomformat::sample_metadata(biom)
## prepare data
my_data <- prepare_data(counts = counts, covariates = covariates)
str(my_data)
```

### Reading from a phyloseq-class object

Likewise, preparing data from a phyloseq-class object requires the
bioconductor package
[phyloseq](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html).
This package is **not** a standard dependency of PLNmodels and needs to
be installed separately.

You can easily prepare your data from a phyloseq object using the
following steps:

- extract the count table with `phyloseq::otu_table()`
- extract the covariates with `phyloseq::sample_data()` (or build your
  own)
- feed them to `prepare_data`

as illustrated below:

``` r

## If biomformat is not installed, uncomment the following lines
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
library(phyloseq)
data("enterotype")
## extract counts
counts <- as(phyloseq::otu_table(enterotype), "matrix")
## extract covariates (or prepare your own)
covariates <- phyloseq::sample_data(enterotype)
## prepare data
my_data <- prepare_data(counts = counts, covariates = covariates)
str(my_data)
```

## Mathematical details about the offsets

We detail here the mathematical background behind the various offsets
and the way they are computed. Note $`\mathbf{Y} = (Y_{ij})`$ the counts
matrix where $`Y_{ij}`$ is the count of species $`j`$ in sample $`i`$.
Assume that there are $`p`$ species and $`n`$ samples in total. The
offset of sample $`i`$ is noted $`O_i`$ and computed in the following
way.

### Total Sum Scaling

Offsets are simply the total counts of a sample (frequently called
depths in the metabarcoding literature):
``` math
O_i = \sum_{j=1}^p Y_{ij}
```

### Cumulative Sum Scaling

Positive counts are used to compute sample-specific quantiles $`q_i^l`$
and cumulative sums $`s_i^l`$ defined as
``` math
q_i^l = \min \{q \text{ such that } \sum_j 1_{Y_{ij} \leq q} \geq l \sum_j 1_{Y_{ij} > 0} \} \qquad s_i^l = \sum_{j: Y_{ij} \leq q_i^l} Y_{ij}
```
The sample-specific quantiles are then used to compute reference
quantiles defined as $`q^l = \text{median} \{q^i_l\}`$ and median
average deviation around the quantile $`q^l`$ as
$`d^l = \text{median} |q_i^l - q^l|`$. The method then searches for the
smallest quantile $`l`$ for which it detects instability, defined as
large relative increase in the $`d^l`$. Formally, $`\hat{l}`$ is the
smallest $`l`$ satisfying $`\frac{d^{l+1} - d^l}{d^l} \geq 0.1`$. The
scaling sample-specific offset are then chosen as:
``` math
O_i = s_i^{\hat{l}} / \text{median}_i \{ s_i^{\hat{l}} \}
```
Dividing by the median of the $`s_i^{\hat{l}}`$ ensures that offsets are
centered around $`1`$ and compare sizes differences with respect to the
reference sample. Note also that the reference quantiles $`q^l`$ can be
computed using either the median (default, as in the original Paulson et
al. ([2013](#ref-CSS)) paper) or the mean, by specifying
`reference = mean`, as implemented in `metagenomeseq`.

### Relative Log Expression

A reference sample $`(q_j)_j`$ is first built by computing the geometric
means of each species count:
``` math
q_j = \exp \left( \frac{1}{n} \sum_{i} \log(Y_{ij})\right)
```
Each sample is then compared to the reference sample to compute one
ratio per species and the final offset $`O_i`$ is the median of those
ratios:
``` math
O_i = \text{median}_j \frac{Y_{ij}}{q_j}
```
The method fails when no species is shared across all sample (as all
$`q_j`$ are then $`0`$) or when a sample shares less than 50% of species
with the reference (in which case the median of the ratios may be null
or infinite). The problem can be alleviated by adding pseudocounts to
the $`c_{ij}`$ with `pseudocounts = 1` or using positive counts in the
computations (`type = "poscounts"`)

### Geometric Mean of Pairwise Ratio

This method is similar to RLE but does create a reference sample.
Instead, each sample is compared to each other to compute a median ratio
(similar to RLE)
``` math
r_{ii'} = {\text{median}}_{j: Y_{ij}.Y_{i'j} > 0} \frac{Y_{ij}}{Y_{i'j}}
```
The offset is then taken as the median of all the $`r_{ii'}`$:
``` math
O_i = \text{median}_{i' != i} r_{ii'}
```
The method fails when there is only one sample in the data set or when a
sample shares no species with any other.

### Wrench normalisation

This method is fully detailed in Kumar et al. ([2018](#ref-Kumar2018))
and we only provide a barebone implementation corresponding to the
defaults parameters of `Wrench::wrench()`. Assume that samples belong to
$`K`$ discrete groups and note $`g_i`$ the group of sample $`i`$. Wrench
is based on the following (simplified) log-normal model for counts:
``` math
Y_{ij} \sim \pi_{ij} \delta_0 + (1 - \pi_{ij})\log\mathcal{N}(\mu_{ij}, \sigma^2_j)
```
where the $`Y_{ij}`$ are independent and the mean $`\mu_{ij}`$ is
decomposed as:
``` math
\mu_{ij} = \underbrace{\log{p_{0j}}}_{\text{log-ref. prop.}} +  \underbrace{\log{d_i}}_{\text{log-depth}} + \underbrace{\log{\zeta_{0g_i}}}_{\text{log effect of group } g_i} + \underbrace{a_{i}}_{\text{(f|m)ixed effect}} + \underbrace{b_{ij}}_{\text{mixed effects}}
```
where the random effects are independents centered gaussian and the
depths is the total sum of counts:
``` math
\begin{align*}
d_i & = \sum_{j=1}^p c_{ij} \\
b_{ij} & \sim \mathcal{N}(0, \eta^2_{g_i}) \\
\end{align*}
```

The **net** log fold change $`\theta_{ij}`$ of the **proportion ratio**
$`r_{ij} = c_{ij} / d_i p_{0j}`$ of species $`j`$ relative to the
reference is
$`\log(\theta_{ij}) \overset{\Delta}{=} \mathbb{E}[\log(r_{ij}) | a_i, b_{ij}] = \log{\zeta_{0g_i}} + a_i + b_{ij}`$.
We can decompose it as $`\theta_{ij} = \Lambda_i^{-1} v_{ij}`$ where
$`\Lambda_i^{-1}`$ is the *compositional correction factor* and
$`v_{ij}`$ is the fold change of **true abundances**.

With the above notations, the net fold change compounds both the fold
change of true abundances and the compositional correction factors. With
the assumption that the $`b_{ij}`$ are centered,
$`\log(\hat{\Lambda}_i)`$ can be estimated through a robust average of
the $`\hat{\theta}_{ij}`$, which can themselves be computed from the
log-ratio of proportions.

We detail here how the different parameters and/or effects are
estimated.

- The reference proportions $`p_{0j}`$ are constructed as averages of
  the sample proportions $`p_{ij}`$ and the ratio are derived from both
  quantities
  ``` math
  p_{ij} = \frac{Y_{ij}}{\sum_{j=1}^p Y_{ij}} \qquad p_{0j} = \frac{1}{n} \sum_{i=1}^n p_{ij} \qquad r_{ij} = \frac{p_{ij}}{p_{0j}}
  ```
- The probabilities of absence $`\pi_{ij}`$ are estimated by fitting the
  following Bernoulli models:
  ``` math
  1_{\{Y_{ij} = 0\}} \sim \mathcal{B}(\pi_{j}^{d_i})
  ```
  and setting $`\hat{\pi}_{ij} = \hat{\pi}^{d_i}`$
- The species variances $`\sigma^2_j`$ are estimated by fitting the
  following linear model (with no zero-inflation component)
  ``` math
  \log Y_{ij} \sim \log(d_i) + \mu_{g_i} + \mathcal{N}(0, \sigma^2_j)
  ```
  Note that in the original `Wrench::wrench()`, the log depth
  $`\log(d_i)`$ is used as predictor but I believe it makes more sense
  to use it an offset.
- set the group proportions $`p_{gj}`$ and group ratios $`r_{gj}`$ to:
  ``` math
  p_{gj} = \frac{\sum_{i : g_i = g} Y_{ij}}{\sum_{j, i : g_i = g} Y_{ij}} \qquad r_{gj} = \frac{p_{gj}}{p_{0j}}
  ```
- Estimate the location and dispersion parameters as:
  ``` math
  \hat{\zeta}_{0g} = \frac{\sum_{j=1}^p r_{gj}}{p} 
  \qquad
  \log{r_{g.}} = \frac{\sum_{j: r_{gj} \neq 1} \log{r_{gj}}}{\sum_{j: r_{gj} \neq 0} 1}   
  \qquad
  \hat{\eta}_{g}^2 = \frac{\sum_{j: r_{gj} \neq 1} (\log{r_{gj}} - \log{r_{g.}})^2}{\sum_{j: r_{gj} \neq 0} 1}
  ```
- Estimate the mixed effects as shrunken (and scaled) averages of the
  ratios
  ``` math
  \hat{a}_i = \frac{\sum_{j = 1}^p \frac{1}{\hat{\eta}^2_{g_i} + \hat{\sigma}^2_j} (\log{r_{ij}} - \log{\hat{\zeta}_{0g_i}})}{\sum_{j = 1}^p \frac{1}{\hat{\eta}^2_{g_i} + \hat{\sigma}^2_j}}
  \qquad 
  \hat{b}_{ij} = \frac{\hat{\eta}^2_{g_i}}{\hat{\eta}^2_{g_i} + \hat{\sigma}^2_j} \left( \log{r_{ij}} - \log\hat{\zeta}_{0g_i} - \hat{a}_i\right)
  ```
- Estimate the regularized ratios as:
  ``` math
  \hat{\theta}_{ij} = \exp\left( \log\hat{\zeta}_{0g_i} + \hat{a}_i + \hat{b}_{ij} \right)
  ```
- Estimate the compositional correction factors as (weighted) means of
  the regularized (and possibly corrected) ratios:
  ``` math
  \hat{\Lambda}_i = 
  \begin{cases}
  \sum_{j = 1}^p \hat{\theta}_{ij} \bigg/ p& \text{ if type = "simple"} \\
  \sum_{j = 1}^p \hat{\theta}_{ij} e^{-\hat{\sigma}_j^2 / 2} /  w_{ij} \bigg/ \sum_{j=1}^p 1/w_{ij}  & \text{ if type = "wrench"} \\
  \end{cases}
  ```
  where
  $`w_{ij} = (1 - \hat{\pi}_{ij})(\hat{\pi}_{ij} + e^{\hat{\sigma}_j^2 + \hat{\eta}_i^2} - 1)`$.
  The correction term $`e^{\hat{\sigma}_j^2 / 2}`$ arises from the
  relation
  $`\mathbb{E}[r_{ij} | r_{ij} > 0] = \theta_{ij} e^{\sigma_j^2/2}`$ and
  the weight $`w_{ij}`$ are marginal variances:
  $`\mathbb{V}[r_{ij}] = w_{ij}`$.

The offsets are then the product of compositional correction factors and
depths:
``` math
O_i = \frac{\hat{\Lambda}_i}{(\prod_{i = 1}^n \hat{\Lambda}_i)^{1/n}} \times \frac{d_i}{(\prod_{i = 1}^n d_i)^{1/n}}
```

## References

Anders, Simon, and Wolfgang Huber. 2010. “Differential Expression
Analysis for Sequence Count Data.” *Genome Biology* 11 (10): R106.
<https://doi.org/10.1186/gb-2010-11-10-r106>.

Chen, Li, James Reeve, Lujun Zhang, Shengbing Huang, Xuefeng Wang, and
Jun Chen. 2018. “GMPR: A Robust Normalization Method for Zero-Inflated
Count Data with Application to Microbiome Sequencing Data.” *PeerJ* 6
(April): e4600. <https://doi.org/10.7717/peerj.4600>.

Kumar, M. Senthil, Eric V. Slud, Kwame Okrah, Stephanie C. Hicks,
Sridhar Hannenhalli, and Héctor Corrada Bravo. 2018. “Analysis and
Correction of Compositional Bias in Sparse Sequencing Count Data.” *BMC
Genomics* 19 (1). <https://doi.org/10.1186/s12864-018-5160-5>.

McMurdie, Paul J. AND Holmes. 2013. “Phyloseq: An r Package for
Reproducible Interactive Analysis and Graphics of Microbiome Census
Data.” *PLoS ONE* 8 (4): e61217.
<https://doi.org/10.1371/journal.pone.0061217>.

Paulson, Joseph N, O. Colin Stine, Héctor Corrada Bravo, and Mihai Pop.
2013. “Differential Abundance Analysis for Microbial Marker-Gene
Surveys.” *Nat Methods* 10 (September): 1200–1202.
<https://doi.org/10.1038/nmeth.2658>.

[^1]: although a `data.frame` is technically a `list`
