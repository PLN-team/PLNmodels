## Which parameters drift in the slow mode of the outer loop? Plain iterations on the slowest
## case; decomposition of the late steps by block and by species.
src <- readLines("inst/simus_ZIPLN/convergence/squarem.R")
eval(parse(text = src[1:(grep("^plain <- function", src) - 1)]))
args <- commandArgs(TRUE); case_name <- if (length(args)) args[1] else "microcosm site col"; init <- if (length(args) > 1) args[2] else "lm"
RhpcBLASctl::blas_set_num_threads(1)
cs <- Filter(function(x) x$name == case_name, cases)[[1]]
mp <- make_map(cs, init)
s <- mp$s0; hist <- list(); obj <- numeric(0)
for (k in 1:1500) {
  s <- mp$F(s); obj[k] <- s$objective
  if (k %in% c(300, 301, 1499, 1500)) hist[[as.character(k)]] <- s
}
Y <- cs$data$Abundance; zero <- as.matrix(Y) == 0
blocks <- function(a, b) {
  dM <- b$M - a$M; dS <- log(b$S2) - log(a$S2); dR <- b$R - a$R
  c(M_zeros = sqrt(sum(dM[zero]^2)), M_positive = sqrt(sum(dM[!zero]^2)),
    logS2_zeros = sqrt(sum(dS[zero]^2)), logS2_positive = sqrt(sum(dS[!zero]^2)), R = sqrt(sum(dR^2)))
}
cat("objective decrease per iteration at 300 and 1500:", signif(-diff(obj)[c(300, 1499)], 3), "\n")
cat("norm of the step by block, at 300 and at 1500:\n"); print(rbind(it300 = blocks(hist[["300"]], hist[["301"]]), it1500 = blocks(hist[["1499"]], hist[["1500"]])), digits = 3)
s1 <- hist[["1499"]]; s2 <- hist[["1500"]]
step_sp <- sqrt(colSums((s2$M - s1$M)^2) + colSums((log(s2$S2) - log(s1$S2))^2))
o <- order(step_sp, decreasing = TRUE)
cat("\nshare of the step carried by the top species:", round(cumsum(step_sp[o]^2)[c(1, 3, 10, 30)] / sum(step_sp^2), 3), "\n")
top <- o[1:8]
sv <- function(st, j) c(mean_M_zeros = mean(st$M[zero[, j], j]), mean_S2_zeros = mean(st$S2[zero[, j], j]), mean_R_zeros = mean(st$R[zero[, j], j]))
print(data.frame(species = top, zeros = round(colMeans(zero)[top], 2), step = signif(step_sp[top], 3),
                 t(sapply(top, function(j) sv(hist[["300"]], j))), t(sapply(top, function(j) sv(s2, j)))), digits = 3, row.names = FALSE)
