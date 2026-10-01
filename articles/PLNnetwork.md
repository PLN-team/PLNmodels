# Sparse structure estimation for multivariate count data with PLN-network

## Preliminaries

This vignette illustrates the standard use of the `PLNnetwork` function
and the methods accompanying the R6 Classes `PLNnetworkfamily` and
`PLNnetworkfit`.

### Requirements

The packages required for the analysis are **PLNmodels** plus some
others for data manipulation and representation:

``` r

library(PLNmodels)
library(ggplot2)
```

### Data set

We illustrate our point with the trichoptera data set, a full
description of which can be found in [the corresponding
vignette](https://pln-team.github.io/PLNmodels/articles/Trichoptera.md).
Data preparation is also detailed in [the specific
vignette](https://pln-team.github.io/PLNmodels/articles/Import_data.md).

``` r

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
```

The `trichoptera` data frame stores a matrix of counts
(`trichoptera$Abundance`), a matrix of offsets (`trichoptera$Offset`)
and some vectors of covariates (`trichoptera$Wind`,
`trichoptera$Temperature`, etc.)

### Mathematical background

The network model for multivariate count data that we introduce in
Chiquet et al. ([2019](#ref-PLNnetwork)) is a variant of the Poisson
Lognormal model of Aitchison and Ho ([1989](#ref-AiH89)), see [the PLN
vignette](https://pln-team.github.io/PLNmodels/articles/PLN.md) as a
reminder. Compare to the standard PLN model we add a sparsity constraint
on the inverse covariance matrix
$`{\boldsymbol\Sigma}^{-1}\triangleq \boldsymbol\Omega`$ by means of the
$`\ell_1`$-norm, such that $`\|\boldsymbol\Omega\|_1 < c`$. PLN-network
is the equivalent of the sparse multivariate Gaussian model ([Banerjee
et al. 2008](#ref-banerjee2008)) in the PLN framework. It relates some
$`p`$-dimensional observation vectors $`\mathbf{Y}_i`$ to some
$`p`$-dimensional vectors of Gaussian latent variables $`\mathbf{Z}_i`$
as follows
``` math
\begin{equation}
  \begin{array}{rcl}
  \text{latent space } &   \mathbf{Z}_i \sim \mathcal{N}\left({\boldsymbol\mu},\boldsymbol\Omega^{-1}\right) &  \|\boldsymbol\Omega\|_1 < c \\
  \text{observation space } &  Y_{ij} | Z_{ij} \quad \text{indep.} & Y_{ij} | Z_{ij} \sim \mathcal{P}\left(\exp\{Z_{ij}\}\right)
  \end{array}
\end{equation}
```

The parameter $`{\boldsymbol\mu}`$ corresponds to the main effects and
the latent covariance matrix $`\boldsymbol\Sigma`$ describes the
underlying structure of dependence between the $`p`$ variables.

The $`\ell_1`$-penalty on $`\boldsymbol\Omega`$ induces sparsity and
selection of important direct relationships between entities. Hence, the
support of $`\boldsymbol\Omega`$ correspond to a network of underlying
interactions. The sparsity level ($`c`$ in the above mathematical
model), which corresponds to the number of edges in the network, is
controlled by a penalty parameter in the optimization process sometimes
referred to as $`\lambda`$. All mathematical details can be found in
Chiquet et al. ([2019](#ref-PLNnetwork)).

#### Covariates and offsets

Just like PLN, PLN-network generalizes to a formulation close to a
multivariate generalized linear model where the main effect is due to a
linear combination of $`d`$ covariates $`\mathbf{x}_i`$ and to a vector
$`\mathbf{o}_i`$ of $`p`$ offsets in sample $`i`$. The latent layer then
reads
``` math
\begin{equation}
  \mathbf{Z}_i \sim \mathcal{N}\left({\mathbf{o}_i + \mathbf{x}_i^\top\mathbf{B}},\boldsymbol\Omega^{-1}\right), \qquad \|\boldsymbol\Omega\|_1 < c ,
\end{equation}
```
where $`\mathbf{B}`$ is a $`d\times p`$ matrix of regression parameters.

#### Alternating optimization

Regularization via sparsification of $`\boldsymbol\Omega`$ and
visualization of the consecutive network is the main objective in
PLN-network. To reach this goal, we need to first estimate the model
parameters. Inference in PLN-network focuses on the regression
parameters $`\mathbf{B}`$ and the inverse covariance
$`\boldsymbol\Omega`$. Technically speaking, we adopt a variational
strategy to approximate the $`\ell_1`$-penalized log-likelihood function
and optimize the consecutive sparse variational surrogate with an
optimization scheme that alternates between two step

1.  a gradient-ascent-step, performed with the CCSA algorithm of
    Svanberg ([2002](#ref-Svan02)) implemented in the C++ library
    ([Johnson 2011](#ref-nlopt)), which we link to the package.
2.  a penalized log-likelihood step, performed with the graphical-Lasso
    of Friedman et al. ([2008](#ref-FHT08)), with an internal C++ port
    of the GLASSOFAST implementation ([Sustik and Calderhead
    2012](#ref-glassofast)), also available to users as
    [`graphical_lasso()`](https://pln-team.github.io/PLNmodels/reference/graphical_lasso.md).

More technical details can be found in Chiquet et al.
([2019](#ref-PLNnetwork))

## Analysis of trichoptera data with a PLNnetwork model

In the package, the sparse PLN-network model is adjusted with the
function `PLNnetwork`, which we review in this section. This function
adjusts the model for a series of value of the penalty parameter
controlling the number of edges in the network. It then provides a
collection of objects with class `PLNnetworkfit`, corresponding to
networks with different levels of density, all stored in an object with
class `PLNnetworkfamily`.

### Adjusting a collection of network - a.k.a. a regularization path

`PLNnetwork` finds an hopefully appropriate set of penalties on its own.
This set can be controlled by the user, but use it with care and check
details in
[`?PLNnetwork`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md).
The collection of models is fitted as follows:

``` r

network_models <- PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data = trichoptera)
```

    ## 
    ##  Initialization...
    ##  Adjusting 30 PLN with sparse inverse covariance estimation
    ##  Joint optimization alternating gradient descent and graphical-lasso
    ##  sparsifying penalty = 0.8883632     sparsifying penalty = 0.8205552     sparsifying penalty = 0.7579229     sparsifying penalty = 0.7000713     sparsifying penalty = 0.6466354     sparsifying penalty = 0.5972783     sparsifying penalty = 0.5516886     sparsifying penalty = 0.5095787     sparsifying penalty = 0.470683  sparsifying penalty = 0.4347561     sparsifying penalty = 0.4015716     sparsifying penalty = 0.37092   sparsifying penalty = 0.342608  sparsifying penalty = 0.316457  sparsifying penalty = 0.2923021     sparsifying penalty = 0.2699909     sparsifying penalty = 0.2493827     sparsifying penalty = 0.2303476     sparsifying penalty = 0.2127653     sparsifying penalty = 0.1965251     sparsifying penalty = 0.1815246     sparsifying penalty = 0.1676689     sparsifying penalty = 0.1548709     sparsifying penalty = 0.1430497     sparsifying penalty = 0.1321309     sparsifying penalty = 0.1220454     sparsifying penalty = 0.1127298     sparsifying penalty = 0.1041253     sparsifying penalty = 0.09617746    sparsifying penalty = 0.08883632 
    ##  Post-treatments
    ##  DONE!

Note the use of the `formula` object to specify the model, similar to
the one used in the function `PLN`.

### Structure of `PLNnetworkfamily`

The `network_models` variable is an `R6` object with class
`PLNnetworkfamily`, which comes with a couple of methods. The most basic
is the `show/print` method, which sends a very basic summary of the
estimation process:

``` r

network_models
```

    ## --------------------------------------------------------
    ## COLLECTION OF 30 POISSON LOGNORMAL MODELS
    ## --------------------------------------------------------
    ##  Task: Network Inference 
    ## ========================================================
    ##  - 30 penalties considered: from 0.08883632 to 0.8883632 
    ##  - Best model (greater BIC): lambda = 0.552 
    ##  - Best model (greater EBIC): lambda = 0.552

One can also easily access the successive values of the criteria in the
collection

``` r

network_models$criteria %>% head() %>% knitr::kable()
```

| param | nb_param | loglik | BIC | AIC | ICL | n_edges | EBIC | pen_loglik | density | stability |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|
| 0.8883632 | 35 | -1109.600 | -1177.707 | -1144.600 | -2187.804 | 1 | -1180.540 | -1131.228 | 0.0073529 | NA |
| 0.8205552 | 35 | -1108.511 | -1176.617 | -1143.511 | -2167.318 | 1 | -1179.451 | -1129.466 | 0.0073529 | NA |
| 0.7579229 | 35 | -1107.916 | -1176.023 | -1142.916 | -2153.575 | 1 | -1178.856 | -1128.113 | 0.0073529 | NA |
| 0.7000713 | 35 | -1107.575 | -1175.682 | -1142.575 | -2148.803 | 1 | -1178.515 | -1126.551 | 0.0073529 | NA |
| 0.6466354 | 35 | -1107.298 | -1175.405 | -1142.298 | -2145.526 | 1 | -1178.238 | -1125.043 | 0.0073529 | NA |
| 0.5972783 | 35 | -1107.055 | -1175.162 | -1142.055 | -2142.413 | 1 | -1177.995 | -1123.650 | 0.0073529 | NA |

A diagnostic of the optimization process is available via the
`convergence` field:

``` r

network_models$convergence %>% head() %>% knitr::kable()
```

|  | param | nb_param | status | backend | objective | iterations | convergence | glasso_nonconverged | glasso_stalled | glasso_indefinite | n_floor | floored |
|:---|---:|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| out | 0.8883632 | 35 | 3 | newton | 1109.6 | 20 | 4.766354e-05 | 0 | 0 | 0 | 0 | FALSE |
| elt | 0.8205552 | 35 | 3 | newton | 1108.511 | 20 | 1.710477e-05 | 0 | 0 | 0 | 0 | FALSE |
| elt.1 | 0.7579229 | 35 | 3 | newton | 1107.916 | 18 | 9.943289e-06 | 0 | 0 | 0 | 0 | FALSE |
| elt.2 | 0.7000713 | 35 | 3 | newton | 1107.575 | 7 | 9.48024e-06 | 0 | 0 | 0 | 0 | FALSE |
| elt.3 | 0.6466354 | 35 | 3 | newton | 1107.298 | 5 | 9.947209e-06 | 0 | 0 | 0 | 0 | FALSE |
| elt.4 | 0.5972783 | 35 | 3 | newton | 1107.055 | 5 | 8.646742e-06 | 0 | 0 | 0 | 0 | FALSE |

An nicer view of this output comes with the option “diagnostic” in the
`plot` method:

``` r

plot(network_models, "diagnostic")
```

![](PLNnetwork_files/figure-html/diagnostic-1.png)

### Exploring the path of networks

By default, the `plot` method of `PLNnetworkfamily` displays evolution
of the criteria mentioned above, and is a good starting point for model
selection:

``` r

plot(network_models)
```

![](PLNnetwork_files/figure-html/plot-1.png)

Note that we use the original definition of the BIC/ICL criterion
($`\texttt{loglik} - \frac{1}{2}\texttt{pen}`$), which is on the same
scale as the log-likelihood. A [popular
alternative](https://en.wikipedia.org/wiki/Bayesian_information_criterion)
consists in using $`-2\texttt{loglik} + \texttt{pen}`$ instead. You can
do so by specifying `reverse = TRUE`:

``` r

plot(network_models, reverse = TRUE)
```

![](PLNnetwork_files/figure-html/plot-reverse-1.png)

In this case, the variational lower bound of the log-likelihood is
hopefully strictly increasing (or rather decreasing if using
`reverse = TRUE`) with a lower level of penalty (meaning more edges in
the network). The same holds true for the penalized counterpart of the
variational surrogate. Generally, smoothness of these criteria is a good
sanity check of optimization process. BIC and its extended-version
high-dimensional version EBIC are classically used for selecting the
correct amount of penalization with sparse estimator like the one used
by PLN-network. However, we will consider later a more robust albeit
more computationally intensive strategy to chose the appropriate number
of edges in the network.

To pursue the analysis, we can represent the coefficient path (i.e.,
value of the edges in the network according to the penalty level) to see
if some edges clearly come off. An alternative and more intuitive view
consists in plotting the values of the partial correlations along the
path, which can be obtained with the options `corr = TRUE`. To this end,
we provide the S3 function `coefficient_path`:

``` r

coefficient_path(network_models, corr = TRUE) %>%
  ggplot(aes(x = Penalty, y = Coeff, group = Edge, colour = Edge)) +
    geom_line(show.legend = FALSE) +  coord_transform(x="log10") + theme_bw()
```

![](PLNnetwork_files/figure-html/path_coeff-1.png)

### Model selection issue: choosing a network

To select a network with a specific level of penalty, one uses the
`getModel(lambda)` S3 method. We can also extract the best model
according to the BIC or EBIC with the method
[`getBestModel()`](https://pln-team.github.io/PLNmodels/reference/getBestModel.md).

``` r

model_pen <- getModel(network_models, network_models$penalties[20]) # give some sparsity
model_BIC <- getBestModel(network_models, "BIC")   # if no criteria is specified, the best BIC is used
```

An alternative strategy is to use StARS ([Liu et al. 2010](#ref-stars)),
which performs resampling to evaluate the robustness of the network
along the path of solutions in a similar fashion as the stability
selection approach of Meinshausen and Bühlmann
([2010](#ref-stabilitySelection)), but in a network inference context.

Resampling can be computationally demanding but is easily parallelized:
the function `stability_selection` relies on
[`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) to
perform parallel computing. Set the number of workers with the
`mc.cores` option (forking-based, so only effective on Unix-like
systems; ignored on Windows):

``` r

options(mc.cores = 2)
```

We first invoke `stability_selection` explicitly for pedagogical
purpose. In this case, we need to build our sub-samples manually:

``` r

n <- nrow(trichoptera)
subs <- replicate(10, sample.int(n, size = n/2), simplify = FALSE)
stability_selection(network_models, subsamples = subs)
```

    ## 
    ## Stability Selection for PLNnetwork: 
    ## subsampling: ++++++++++

Requesting ‘StARS’ in `gestBestmodel` automatically invokes
`stability_selection` with 20 sub-samples, if it has not yet been run.

``` r

model_StARS <- getBestModel(network_models, "StARS")
```

When “StARS” is requested for the first time, `getBestModel`
automatically calls the method `stability_selection` with the default
parameters. After the first call, the stability path is available from
the `plot` function:

``` r

plot(network_models, "stability")
```

![](PLNnetwork_files/figure-html/plot%20stability-1.png)

When you are done, do not forget to get back to the default (sequential)
behavior.

``` r

options(mc.cores = 1)
```

### Structure of a `PLNnetworkfit`

The variables `model_BIC`, `model_StARS` and `model_pen` are other
`R6Class` objects with class `PLNnetworkfit`. They all inherits from the
class `PLNfit` and thus own all its methods, with a couple of specific
one, mostly for network visualization purposes. Most fields and methods
are recalled when such an object is printed:

``` r

model_StARS
```

    ## Poisson Lognormal with sparse inverse covariance (penalty = 0.597)
    ## ==================================================================
    ##  nb_param    loglik       BIC       AIC       ICL n_edges      EBIC pen_loglik
    ##        35 -1107.055 -1175.162 -1142.055 -2142.413       1 -1177.995   -1123.65
    ##  density
    ##    0.007
    ## ==================================================================
    ## * Useful fields
    ##     $model_par, $latent, $latent_pos, $var_par, $optim_par
    ##     $loglik, $BIC, $ICL, $loglik_vec, $nb_param, $criteria
    ## * Useful S3 methods
    ##     print(), coef(), sigma(), vcov(), fitted()
    ##     predict(), predict_cond(), standard_error()
    ## * Additional fields for sparse network
    ##     $EBIC, $density, $penalty 
    ## * Additional S3 methods for network
    ##     plot.PLNnetworkfit()

The `plot` method provides a quick representation of the inferred
network, with various options (either as a matrix, a graph, and always
send back the plotted object invisibly if users needs to perform
additional analyses).

``` r

my_graph <- plot(model_StARS, plot = FALSE)
my_graph
```

    ## IGRAPH e5e6f9c UNW- 17 1 -- 
    ## + attr: name (v/c), label (v/c), label.cex (v/n), size (v/n),
    ## | label.color (v/c), weight (e/n), width (e/n), color (e/c)
    ## + edge from e5e6f9c (vertex names):
    ## [1] Hfo--Hsp

``` r

plot(model_StARS)
```

![](PLNnetwork_files/figure-html/stars_network-1.png)

``` r

plot(model_StARS, type = "support", output = "corrplot")
```

![](PLNnetwork_files/figure-html/stars_network-2.png)

We can finally check that the fitted value of the counts – even with
sparse regularization of the covariance matrix – are close to the
observed ones:

``` r

data.frame(
  fitted   = as.vector(fitted(model_StARS)),
  observed = as.vector(trichoptera$Abundance)
) %>%
  ggplot(aes(x = observed, y = fitted)) +
    geom_point(size = .5, alpha =.25 ) +
    scale_x_log10(limits = c(1,1000)) +
    scale_y_log10(limits = c(1,1000)) +
    theme_bw() + annotation_logticks()
```

![fitted value vs.
observation](PLNnetwork_files/figure-html/fitted-1.png)

fitted value vs. observation

## References

Aitchison, J., and C. H. Ho. 1989. “The Multivariate Poisson-Log Normal
Distribution.” *Biometrika* 76 (4): 643–53.

Banerjee, Onureena, Laurent El Ghaoui, and Alexandre d’Aspremont. 2008.
“Model Selection Through Sparse Maximum Likelihood Estimation for
Multivariate Gaussian or Binary Data.” *Journal of Machine Learning
Research* 9 (Mar): 485–516.

Chiquet, Julien, Stephane Robin, and Mahendra Mariadassou. 2019.
“Variational Inference for Sparse Network Reconstruction from Count
Data.” In *Proceedings of the 36th International Conference on Machine
Learning*, edited by Kamalika Chaudhuri and Ruslan Salakhutdinov, vol.
97. Proceedings of Machine Learning Research. PMLR.
[http://proceedings.mlr.press/v97/chiquet19a.html](http://proceedings.mlr.press/v97/chiquet19a.md).

Friedman, J., T. Hastie, and R. Tibshirani. 2008. “Sparse Inverse
Covariance Estimation with the Graphical Lasso.” *Biostatistics* 9 (3):
432–41.

Johnson, Steven G. 2011. *The NLopt Nonlinear-Optimization Package*.
<https://nlopt.readthedocs.io/en/latest/>.

Liu, Han, Kathryn Roeder, and Larry Wasserman. 2010. “Stability Approach
to Regularization Selection (StARS) for High Dimensional Graphical
Models.” *Proceedings of the 23rd International Conference on Neural
Information Processing Systems - Volume 2* (USA), 1432–40.

Meinshausen, Nicolai, and Peter Bühlmann. 2010. “Stability Selection.”
*Journal of the Royal Statistical Society: Series B (Statistical
Methodology)* 72 (4): 417–73.

Sustik, Mátyás A, and Ben Calderhead. 2012. “GLASSOFAST: An Efficient
GLASSO Implementation.” *UTCS Technical Report TR-12-29 2012*.

Svanberg, Krister. 2002. “A Class of Globally Convergent Optimization
Methods Based on Conservative Convex Separable Approximations.” *SIAM
Journal on Optimization* 12 (2): 555–73.
