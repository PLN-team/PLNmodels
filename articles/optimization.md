# Optimization methods in PLNmodels

All the models in **PLNmodels** are fitted by *variational inference*:
each one maximises a tractable lower bound of an otherwise intractable
log-likelihood. This document describes, model by model, **which
optimizer is used** and, for the home-made `"builtin"` backend, **the
mathematical principles** behind it.

## 1 The common variational framework

Every model shares the same Poisson log-normal core. For a count matrix
$`\mathbf Y \in \mathbb N^{n\times p}`$, offsets $`\mathbf O`$ and
covariates $`\mathbf X`$, the latent Gaussian vectors
$`\mathbf Z_i \sim \mathcal N(\mathbf x_i^\top\mathbf B,\, \boldsymbol\Sigma)`$
drive Poisson emissions
$`Y_{ij}\,|\,Z_{ij}\sim\mathcal P(\exp(O_{ij}+Z_{ij}))`$.

Because the marginal likelihood integrates out $`\mathbf Z`$, inference
maximises a **variational lower bound** (the ELBO) obtained with a
factorised Gaussian approximation
$`q(\mathbf Z_i)=\mathcal N(\mathbf x_i^\top\mathbf B+\mathbf m_i,\ \mathrm{diag}(\mathbf s_i^2))`$.
Writing $`\boldsymbol\Omega=\boldsymbol\Sigma^{-1}`$,
$`\mathbf A=\exp(\mathbf O+\mathbf X\mathbf B+\mathbf M+\tfrac12\mathbf S^2)`$
and $`\mathbf Z=\mathbf O+\mathbf X\mathbf B+\mathbf M`$, the bound
reads

``` math
J(\mathbf B,\boldsymbol\Omega,\mathbf M,\mathbf S^2)=
\sum_{ij}\Big(Y_{ij}Z_{ij}-A_{ij}+\tfrac12\log S^2_{ij}\Big)
-\tfrac12\sum_i\Big(\mathbf m_i^\top\boldsymbol\Omega\,\mathbf m_i+\mathrm{tr}(\mathrm{diag}(\mathbf s_i^2)\boldsymbol\Omega)\Big)
+\tfrac{n}{2}\log\det\boldsymbol\Omega + \text{cst}.
```

The **variational parameters** are the means $`\mathbf M`$ and variances
$`\mathbf S^2`$ (one pair per observation and per latent dimension); the
**model parameters** are the regression coefficients $`\mathbf B`$ and
the covariance structure $`\boldsymbol\Sigma`$ (or
$`\boldsymbol\Omega`$). Optimization alternates a **VE-step** (update
$`\mathbf M,\mathbf S^2`$ at fixed model parameters) with an **M-step**
(update $`\mathbf B,\boldsymbol\Sigma`$), until the ELBO stabilises. The
different models change *what structure is imposed* on
$`\boldsymbol\Sigma`$ (or on the emission), not this general scheme.
Throughout, the variances are optimised in the unconstrained
parametrisation $`\boldsymbol\psi=\log\mathbf S^2`$.

## 2 Which optimizer for which model

Three backends are available, selected through the `backend` argument of
each `*_param()` helper:

- **`nlopt`** — a first-order solver from the NLOPT library (the *CCSAQ*
  conservative convex separable approximation). Robust and cheap per
  iteration but needs many iterations.
- **`builtin`** — the home-made *second-order* optimizer implemented
  directly in the package (the subject of §3). It uses the analytic
  Hessian of the bound.
- **`torch`** — automatic differentiation with first-order optimizers
  (Adam, RPROP, …). **Experimental**; typically reaches a lower bound
  than the other two.

| Model | Imposed structure | Backends | Default | builtin optimizer |
|----|----|----|----|----|
| `PLN` | $`\boldsymbol\Sigma`$ full / diagonal / spherical / fixed / genetic | nlopt, builtin, torch | **nlopt** | joint coordinate Newton (§3.1) |
| `PLNLDA` | PLN with per-group means | (via `PLN`) | **nlopt** | as `PLN` |
| `PLNPCA` | rank-$`q`$$`\boldsymbol\Sigma=\mathbf C\mathbf C^\top`$ | nlopt, builtin, torch | **nlopt** | profiled trust-region Newton (§3.4) |
| `PLNnetwork` | sparse $`\boldsymbol\Omega`$ ($`\ell_1`$) | builtin, nlopt | **builtin** | Newton VE + graphical-lasso M-step (§3.2) |
| `PLNmixture` | mixture of $`K`$ PLN | builtin, nlopt, torch | **builtin** | outer mixture-EM + per-component builtin (§3.5) |
| `ZIPLN` | zero-inflation + $`\boldsymbol\Sigma`$ structure | builtin, nlopt | **builtin** | joint Newton on $`(\mathbf M,\boldsymbol\psi,\mathbf R)`$ (§3.3) |
| `ZIPLNnetwork` | zero-inflation + sparse $`\boldsymbol\Omega`$ | builtin, nlopt | **builtin** | as `ZIPLN` + graphical-lasso |

The default is `nlopt` for `PLN`/`PLNPCA` (conservative, well tested)
and `builtin` for the penalised / structured models, where the
second-order steps pay off most.

## 3 The `builtin` backend: mathematical principles

The common thread is that the **VE-step is concave and separable**,
which makes it an ideal target for Newton’s method with a cheaply
invertible Hessian; the M-steps are then either closed-form or a small
convex problem.

### 3.1 PLN and its covariance variants — joint coordinate Newton

At fixed $`(\mathbf B,\boldsymbol\Omega)`$ the bound decouples across
observations, and for each observation the pair
$`(\mathbf m_i,\boldsymbol\psi_i)`$ is updated by a **Newton step that
is joint in $`(m,\psi)`$ but diagonal across the $`p`$ coordinates**.
Per entry $`(i,j)`$ the gradient is

``` math
\partial_{m}J = (\mathbf M\boldsymbol\Omega)_{ij}+A_{ij}-Y_{ij},
\qquad
\partial_{\psi}J = \tfrac12\big(A_{ij}S^2_{ij}+\omega_{jj}S^2_{ij}-1\big),
```

and the $`2\times2`$ Hessian block, using $`\omega_{jj}=\Omega_{jj}`$,
is

``` math
H=\begin{pmatrix} A_{ij}+\omega_{jj} & \tfrac12 A_{ij}S^2_{ij}\\[2pt]
\tfrac12 A_{ij}S^2_{ij} & \tfrac12 S^2_{ij}\big(A_{ij}(1+\tfrac12 S^2_{ij})+\omega_{jj}\big)\end{pmatrix},
```

so the joint step is obtained by inverting a $`2\times2`$ matrix in
closed form (`compute_joint_step_MS` in `src/covariance_pln.h`). Solving
$`(m,\psi)`$*together* captures the coupling between the mean and its
variance and converges in a handful of passes. The cross term
$`\tfrac12 A_{ij}S^2_{ij}`$ is exactly this coupling; a separate fixed
point for $`\psi`$ would ignore it.

The **M-step** is analytic. The regression coefficients are profiled
out,
$`\widehat{\mathbf B}=(\mathbf X^\top\mathbf W\mathbf X)^{-1}\mathbf X^\top\mathbf W\,\mathbf M_{\text{full}}`$,
and the covariance is the empirical variational covariance,
$`\widehat{\boldsymbol\Sigma}=\tfrac1n\big(\mathbf M^\top\mathbf M+\mathrm{diag}(\textstyle\sum_i\mathbf s_i^2)\big)`$,
projected onto the requested structure (full, diagonal, spherical,
fixed, or the genetic $`\sigma^2(\rho\mathbf C+(1-\rho)\mathbf I)`$
solved by 1-D bisection on $`\rho`$). All variants share one templated
implementation (`CovTraitsBase`), each supplying only its
covariance-specific primitives.

### 3.2 PLNnetwork — Newton VE-step + graphical-lasso M-step

`PLNnetwork` adds an $`\ell_1`$ penalty on the precision
$`\boldsymbol\Omega`$ to recover a sparse conditional-dependence graph.
The VE-step is unchanged (§3.1 Newton), and the covariance M-step
becomes a **graphical lasso**:

``` math
\widehat{\boldsymbol\Omega}=\arg\max_{\boldsymbol\Omega\succ0}\
\tfrac{n}{2}\log\det\boldsymbol\Omega-\tfrac{n}{2}\mathrm{tr}(\widehat{\boldsymbol\Sigma}\boldsymbol\Omega)-\lambda\lVert\boldsymbol\Omega\rVert_{1},
```

solved efficiently with `glassoFast` along the penalty path $`\lambda`$.
Successive penalties are warm-started, and the inception model is only
*partially* converged (`maxit_ve = 1`, `inception_niter = 5`) so that
the latent means do not over-fit the unpenalised optimum before the
sparse grid begins.

### 3.3 ZIPLN — joint Newton on $`(\mathbf M,\boldsymbol\psi,\mathbf R)`$

Zero-inflated PLN augments the model with Bernoulli “excess-zero”
indicators; the variational approximation adds posterior probabilities
$`\mathbf R`$ of being an excess zero. The `builtin` VE-step is again a
**Newton step, now joint in
$`(\mathbf m_i,\boldsymbol\psi_i,\mathbf r_i)`$**, with $`\mathbf R`$
updated *inside* the C++ solver (a closed-form logistic update given
$`\mathbf M,\mathbf S^2`$). The effective Poisson rate becomes
$`(1-\mathbf R)\odot\mathbf A`$, which reuses exactly the `PLN`
covariance machinery with the substitution
$`\mathbf A\!\leftarrow\!(1-\mathbf R)\odot\mathbf A`$. The M-step
updates the zero-inflation coefficients and the covariance
(full/diagonal/spherical/fixed, or sparse via graphical lasso for
`ZIPLNnetwork`).

### 3.4 PLNPCA — profiled trust-region Newton

`PLNPCA` constrains $`\boldsymbol\Sigma=\mathbf C\mathbf C^\top`$ to
rank $`q`$, so the latent is
$`\mathbf Z=\mathbf O+\mathbf X\mathbf B+\mathbf M\mathbf C^\top`$ with
$`q`$-dimensional scores $`\mathbf M\in\mathbb R^{n\times q}`$ and
loadings $`\mathbf C\in\mathbb R^{p\times q}`$. Naïve block-coordinate
ascent stalls on the **saddle points** of this non-convex landscape; the
`builtin` optimizer avoids them with a profiled, saddle-aware
trust-region Newton.

**Profiling.** The score block $`(\mathbf M,\boldsymbol\psi)`$ is
concave and is solved by a per-observation Newton VE-step (a
$`q\times q`$ solve for each $`\mathbf m_i`$, plus a fixed point for
$`\mathbf s_i^2`$). This defines the profiled objective on the loadings
alone,
``` math
g(\mathbf B,\mathbf C)=\max_{\mathbf M,\mathbf S^2}J(\mathbf B,\mathbf C,\mathbf M,\mathbf S^2).
```
By the **envelope theorem** its gradient is free — it is simply the
partial gradient $`\partial J/\partial(\mathbf B,\mathbf C)`$ evaluated
at the optimal $`(\mathbf M,\mathbf S^2)`$.

**Reduced (Schur) Hessian.** Writing $`\theta=(\mathbf B,\mathbf C)`$
(the loadings) and $`\varphi=(\mathbf M,\boldsymbol\psi)`$ (the profiled
scores), the Hessian of $`g`$ is the Schur complement
``` math
H_{\text{red}}=J_{\theta\theta}-J_{\theta\varphi}\,J_{\varphi\varphi}^{-1}\,J_{\varphi\theta}.
```
It is never formed. The inner block $`J_{\varphi\varphi}`$ is
**block-diagonal per observation** ($`2q\times2q`$, inverted
analytically); the cross/outer terms are applied through analytic
directional derivatives of the gradient, so a Hessian–vector product
$`H_{\text{red}}\mathbf v`$ costs only a couple of gradient evaluations
plus $`n`$ small solves.

**Saddle-aware trust region.** The reduced Hessian is *indefinite* (the
landscape is saddle-rich), so $`g`$ is maximised with a **trust-region**
method whose subproblem is solved by a **preconditioned Steihaug
conjugate gradient**. The preconditioner is the analytic diagonal of
$`J_{\theta\theta}`$ (a Jacobi metric that absorbs the very different
scales of $`\mathbf B`$ and $`\mathbf C`$), and the conjugate gradient
follows negative-curvature directions to the trust-region boundary —
precisely the directions that let the iterate *escape saddles* rather
than stall at them. This reaches a variational bound at least as high as
`nlopt` across datasets, in one to two orders of magnitude fewer outer
iterations.

> **Large data**
>
> On large problems the reduced Hessian is strongly indefinite and $`g`$
> is highly non-quadratic for large steps, so the trust region — not the
> linear solver — bounds progress: the outer iteration keeps improving
> the bound but converges slowly and is capped by `maxit_out` (default
> 150) rather than by the gradient tolerance. For the best bound on
> large data, raise it, e.g.
> `PLNPCA_param(backend = "builtin", config_optim = list(maxit_out = 300))`.

### 3.5 PLNmixture — mixture-EM with builtin components

A PLN mixture with $`K`$ components is fitted by a genuine EM loop at
the R level: an **E-step** updates the posterior class probabilities
$`\boldsymbol\tau`$ (a soft-max of the component variational bounds plus
the log mixing proportions), and an **M-step** re-fits each component as
a *weighted* PLN — using the `builtin` (or `nlopt`) optimizer of §3.1
with observation weights $`\tau_{ik}`$ — and updates the shared
covariate effects. The mixing proportions are the column means of
$`\boldsymbol\tau`$.

## 4 Summary

- All models maximise the same kind of Poisson log-normal variational
  bound; the optimizers differ in how they exploit its structure.
- The `builtin` backend is a **second-order** method throughout: a joint
  $`(\text{mean},\text{variance})`$ Newton VE-step everywhere,
  closed-form or graphical-lasso covariance M-steps, and — for `PLNPCA`
  — a profiled, saddle-aware **trust-region Newton** on the loadings.
- `nlopt` remains the default for `PLN`/`PLNPCA` for its robustness;
  `builtin` is the default for the penalised and structured variants,
  where its second-order steps and better optima matter most. `torch` is
  experimental.
