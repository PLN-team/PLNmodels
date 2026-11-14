## Rscript toy_ve_cg.R [out.png]  (figure toy_ve_cg.png in this directory)
## Toy illustration: why the cell-wise Newton VE step of ZIPLN is slow on zeros, and
## why a Newton step by preconditioned conjugate gradients fixes it.
##
## One observation, two species with strongly correlated latent variables. The VE
## objective in the variational means m = (m1, m2) (variances and R held fixed) is
##   f(m) = sum_j A_eff_j(m_j) - Y_j m_j + 1/2 (m - mu)' Omega (m - mu),
##   A_eff_j = (1 - R_j) exp(o_j + m_j + s2_j / 2),
## with Hessian diag(A_eff) + Omega. The cell-wise step uses only the diagonal of the
## Hessian, A_eff_j + Omega_jj: fine when A_eff dominates (large counts), but when the
## zeros are explained by the zero inflation, A_eff ~ 0 and the Hessian is Omega
## itself, whose off-diagonal coupling the cell-wise step ignores.
library(ggplot2); library(patchwork)

sigma2 <- 4; rho <- 0.95
Sigma  <- sigma2 * matrix(c(1, rho, rho, 1), 2)
Omega  <- solve(Sigma)
mu     <- c(1, 1); s2 <- c(0.1, 0.1)

make_case <- function(Y, R) {
  aeff <- function(m) (1 - R) * exp(m + s2 / 2)
  list(Y = Y,
       f    = function(m) sum(aeff(m) - Y * m) + 0.5 * drop(t(m - mu) %*% Omega %*% (m - mu)),
       grad = function(m) aeff(m) - Y + drop(Omega %*% (m - mu)),
       hess = function(m) diag(aeff(m), 2) + Omega)
}

## Newton-type descent with the Armijo backtracking of the package (c1 = 1e-4, halving)
descent <- function(cs, m0, direction, maxit = 2000, tol = 1e-10) {
  m <- m0; path <- matrix(m0, 1)
  for (it in seq_len(maxit)) {
    g <- cs$grad(m); d <- direction(m, g)
    f0 <- cs$f(m); alpha <- 1
    while (cs$f(m - alpha * d) > f0 - 1e-4 * alpha * sum(g * d) && alpha > 1e-10) alpha <- alpha / 2
    m <- m - alpha * d; path <- rbind(path, m)
    if (abs(f0 - cs$f(m)) < tol * (1 + abs(f0))) break
  }
  path
}
## cell-wise step: diagonal of the Hessian only (Jacobi-like)
cellwise <- function(cs) function(m, g) g / diag(cs$hess(m))
## truncated Newton step: at most cg_max PCG iterations on H d = g, preconditioned by the diagonal
newton_cg <- function(cs, cg_max = 5, cg_tol = 0.1) function(m, g) {
  H <- cs$hess(m); P <- 1 / diag(H)
  d <- 0 * g; r <- g; z <- P * r; p <- z; rz <- sum(r * z)
  for (k in seq_len(cg_max)) {
    Hp <- drop(H %*% p); a <- rz / sum(p * Hp)
    d <- d + a * p; r <- r - a * Hp
    if (sqrt(sum(r^2)) < cg_tol^2 * sqrt(sum(g^2))) break
    z <- P * r; rz_new <- sum(r * z); p <- z + rz_new / rz * p; rz <- rz_new
  }
  d
}

cases <- list(
  "Zeros explained by the zero inflation\n(Y = 0, R = 0.99: A_eff ≈ 0)" = make_case(Y = c(0, 0),   R = c(0.99, 0.99)),
  "Positive counts\n(Y = (20, 30), R = 0: A_eff large)"                = make_case(Y = c(20, 30), R = c(0, 0)))
m0 <- c(1.5, -2.5)

paths <- list(); contours <- list(); curves <- list()
for (cn in names(cases)) {
  cs <- cases[[cn]]
  ## optimum by full Newton steps
  m_star <- descent(cs, m0, function(m, g) solve(cs$hess(m), g), tol = 1e-15); f_star <- cs$f(m_star[nrow(m_star), ])
  for (meth in c("cell-wise Newton", "Newton-CG")) {
    P <- descent(cs, m0, if (meth == "cell-wise Newton") cellwise(cs) else newton_cg(cs))
    paths[[length(paths) + 1]] <- data.frame(case = cn, method = meth, it = seq_len(nrow(P)) - 1, m1 = P[, 1], m2 = P[, 2])
    curves[[length(curves) + 1]] <- data.frame(case = cn, method = meth, it = seq_len(nrow(P)) - 1,
                                               gap = pmax(apply(P, 1, cs$f) - f_star, 1e-16))
  }
  lim <- c(-3.5, 4.5)
  g <- expand.grid(m1 = seq(lim[1], lim[2], length.out = 200), m2 = seq(lim[1], lim[2], length.out = 200))
  g$z <- log10(apply(g, 1, cs$f) - f_star + 1e-3)
  contours[[cn]] <- data.frame(case = cn, g)
}
paths <- do.call(rbind, paths); curves <- do.call(rbind, curves); contours <- do.call(rbind, contours)
for (v in c("paths", "curves", "contours")) assign(v, transform(get(v), case = factor(case, levels = names(cases))))
n_it <- aggregate(it ~ case + method, paths, max)
cols <- c("cell-wise Newton" = "#D55E00", "Newton-CG" = "#0072B2")

p1 <- ggplot(contours, aes(m1, m2)) +
  geom_contour(aes(z = z), bins = 18, colour = "grey60", linewidth = 0.3) +
  geom_path(data = paths, aes(colour = method, group = method), linewidth = 0.5) +
  geom_point(data = paths, aes(colour = method), size = 0.8) +
  geom_text(data = transform(n_it, m1 = -Inf, m2 = Inf, lab = paste0(method, ": ", it, " steps")),
            aes(label = lab, colour = method), hjust = -0.05, vjust = ifelse(n_it$method == "Newton-CG", 3.2, 1.5),
            size = 3.2, show.legend = FALSE) +
  facet_wrap(~ case) + coord_equal() +
  scale_colour_manual(values = cols, name = NULL) +
  labs(x = expression(m[1]), y = expression(m[2]),
       title = "VE step of ZIPLN at fixed Omega: two species with correlated latent variables (rho = 0.95)",
       subtitle = "Level lines of the objective in the variational means; same start, Armijo line search") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
p2 <- ggplot(curves, aes(it + 1, gap, colour = method)) + geom_line() +
  facet_wrap(~ case) + scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = cols, name = NULL) +
  labs(x = "Newton steps (+1, log scale)", y = "objective − optimum") +
  theme_minimal(base_size = 11) + theme(legend.position = "none")
out <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "toy_ve_cg.png"  # run from this directory
ggsave(out, p1 / p2 + plot_layout(heights = c(2.2, 1)), width = 10, height = 8.5, dpi = 130, bg = "white")
print(n_it)
