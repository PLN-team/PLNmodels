.logfactorial <- function (n) { # Ramanujan's formula
  n <- n[n!=0]
  return(n*log(n) - n + log(8*n^3 + 4*n^2 + n + 1/30)/6 + log(pi)/2)
}

.trunc <- function(x, M = 300) {
  return(pmin(x, M))
}

logLikPoisson <- function(responses, lambda) {
  loglik <- sum(responses * lambda, na.rm=TRUE) - sum(exp(lambda)) - sum(.logfactorial(responses))
  loglik
}

nullModelPoisson <- function(responses, covariates, offsets) {
  Theta <- do.call(rbind, lapply(1:ncol(responses), function(j)
    coefficients(glm.fit(covariates, responses[, j], offset = offsets[,j], family = stats::poisson()))))
  lambda <- offsets + tcrossprod(covariates, Theta)
  lambda
}

fullModelPoisson <- function(responses) {
  lambda <- log(responses)
  lambda
}

circle <- function(center = c(0, 0), radius = 1, npoints = 100) {
  r = radius
  tt = seq(0, 2 * pi, length = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

GeomCircle <- ggplot2::ggproto("GeomCircle", ggplot2::Geom,
  required_aes = c("x", "y", "radius"),
  default_aes = ggplot2::aes(
    colour = "grey30", fill=NA, alpha=NA, linewidth=1, linetype="solid"),

  draw_key = function (data, params, size)
  {
    grid::circleGrob(
      0.5, 0.5,
      r=0.35,
      gp = grid::gpar(
        col = scales::alpha(data$colour, data$alpha),
        fill = scales::alpha(data$fill, data$alpha),
        lty = data$linetype,
        lwd = data$linewidth
      )
    )
  },

  draw_panel = function(data, panel_scales, coord,  na.rm = TRUE) {
    coords <- coord$transform(data, panel_scales)
    grid::circleGrob(
      x=coords$x, y=coords$y,
      r=coords$radius,
      gp = grid::gpar(
        col = alpha(coords$colour, coords$alpha),
        fill = alpha(coords$fill, coords$alpha),
        lty = coords$linetype,
        lwd = coords$linewidth
      )
    )
  }
)

geom_circle <- function(mapping = NULL, data = NULL, stat = "identity",
                        position = "identity", na.rm = FALSE, show.legend = NA,
                        inherit.aes = TRUE, ...) {
  ggplot2::layer(
    geom = GeomCircle, mapping = mapping,  data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

## Pareil avec un grid.arrange?
#' @importFrom grid grid.newpage pushViewport viewport grid.layout grid.draw
multiplot <- function(..., legend=FALSE, plotlist=NULL, cols) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # Make the panel
  plotCols = cols                       # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols

  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)

  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  if (legend) {
    grid.draw(g_legend(plots[[i]]))
  }

}

PLN_param <- function(control, n, p) {
  ctrl <- list(
    ftol_rel  = ifelse(n < 1.5*p, 1e-6, 1e-8),
    ftol_abs  = 0,
    xtol_rel  = 1e-4,
    xtol_abs  = 1e-4,
    maxeval   = 10000,
    method    = "CCSAQ",
    lbvar     = 1e-4,
    trace     = 1,
    inception = "LM"
  )
  ctrl[names(control)] <- control
  ctrl
}

PLNPCA_param <- function(control, n, p, type = c("init", "main")) {
  type <- match.arg(type)

  ctrl <- switch(match.arg(type),
    "init" = list(
      inception = ifelse(n >= 1.5*p, "PLN", "LM"),
      ftol_rel  = ifelse(n < 1.5*p, 1e-6, 1e-8),
      ftol_abs = 0,
      xtol_rel = 1e-4,
      xtol_abs = 1e-4,
      maxeval  = 10000,
      method   = "CCSAQ",
      lbvar    = 1e-4,
      trace    = 0
    ),
    "main" = list(
      ftol_rel = 1e-6,
      ftol_abs = 0,
      xtol_rel = 1e-4,
      xtol_abs = 1e-4,
      maxeval  = 10000,
      method   = "CCSAQ",
      lbvar    = 1e-4,
      trace    = 1,
      cores    = 1
    )
  )
  ctrl[names(control)] <- control
  ctrl
}

PLNnetwork_param <- function(control, n, p, type = c("init", "main")) {
  type <- match.arg(type)

  ctrl <- switch(match.arg(type),
    "init" = list(
      inception = ifelse(n >= 1.5*p, "PLN", "LM"),
      ftol_rel  = ifelse(n < 1.5*p, 1e-6, 1e-8),
      ftol_abs = 0,
      xtol_rel = 1e-4,
      xtol_abs = 1e-4,
      maxeval  = 10000,
      method   = "CCSAQ",
      lbvar    = 1e-4,
      nPenalties = 20,
      min.ratio = ifelse(n >= 1.5*p, 0.1, 0.05),
      trace = 0),
    "main" = list(
      ftol_out  = 1e-5,
      maxit_out = 50,
      penalize.diagonal = FALSE,
      warm      = FALSE,
      ftol_abs  = 0,
      ftol_rel  = 1e-8,
      xtol_rel  = 1e-4,
      xtol_abs  = 1e-4,
      maxeval   = 10000,
      method    = "CCSAQ",
      lbvar     = 1e-4,
      trace = 1)
  )
  ctrl[names(control)] <- control
  ctrl
}

statusToMessage <- function(status) {
    message <- switch( status,
        "1"  = "NLOPT_SUCCESS: Generic success return value.",
        "2"  = "NLOPT_STOPVAL_REACHED: Optimization stopped because stopval (above) was reached.",
        "3"  = "NLOPT_FTOL_REACHED: Optimization stopped because ftol_rel or ftol_abs (above) was reached.",
        "4"  = "NLOPT_XTOL_REACHED: Optimization stopped because xtol_rel or xtol_abs (above) was reached.",
        "5"  = "NLOPT_MAXEVAL_REACHED: Optimization stopped because maxeval (above) was reached.",
        "6"  = "NLOPT_MAXTIME_REACHED: Optimization stopped because maxtime (above) was reached.",
        "-1" = "NLOPT_FAILURE: Generic failure code.",
        "-2" = "NLOPT_INVALID_ARGS: Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera).",
        "-3" = "NLOPT_OUT_OF_MEMORY: Ran out of memory.",
        "-4" = "NLOPT_ROUNDOFF_LIMITED: Roundoff errors led to a breakdown of the optimization algorithm. In this case, the returned minimum may still be useful.",
        "-5" = "Halted because of a forced termination: the user called nlopt_force_stop(opt) on the optimization's nlopt_opt object opt from the user's objective function.",
        "Return status not recognized."
    )
    message
}
