get_ggplot_ind_map <- function(scores, axes_label, main) {

  if (length(axes_label) > 1 ) {
    if (length(axes_label) >= 3) {
      message("Three axes or more provided. Using only the first two")
    }
    p <- ggplot(scores, aes_(x = ~a1, y = ~a2, label = ~names, colour = ~labels)) +
      geom_hline(yintercept = 0, colour = "gray65") +
      geom_vline(xintercept = 0, colour = "gray65") +
      geom_text(alpha = 0.8, size = 4) +
      ggtitle(main) +
      theme_bw() +
      labs(x = axes_label[1], y = axes_label[2])
  } else {
    p <- ggplot(scores, aes_(x = ~a1, group = ~labels, fill = ~labels, colour = ~labels)) +
      geom_density(alpha = .4) +
      geom_rug() +
      ggtitle(main) +
      theme_bw() +
      labs(x = axes_label[1])
  }
  p
}

get_ggplot_corr <- function(correlations, axes_label, main) {
  if (length(axes_label) > 1 ) {
    if (length(axes_label) >= 3) {
      message("Three axes or more provided. Using only the first two")
    }
    ## data frame with arrows coordinates
    arrows <- data.frame(x1 = rep(0, nrow(correlations)), y1 = rep(0, nrow(correlations)),
                         x2 = correlations$axe1,  y2 = correlations$axe2,
                         names = correlations$names, labels = correlations$labels)
    ## geom_path will do open circles
    p <- ggplot() + xlim(-1.1, 1.1) + ylim(-1.1, 1.1)  +
      geom_hline(yintercept = 0, colour = "gray65") +
      geom_vline(xintercept = 0, colour = "gray65") +
      geom_segment(data = arrows, aes_(x = ~x1, y = ~y1, xend = ~x2, yend = ~y2, colour = ~labels)) +
      geom_text(data = correlations, aes_(x = ~axe1, y = ~axe2, label = ~names, colour = ~labels), size=3) +
      theme_bw() +  theme(legend.position = "none") + ggtitle(main) + labs(x = axes_label[1], y = axes_label[2])
  } else {
    sign  <- rep(c(-1,1), each = ceiling(nrow(correlations)/2))[1:nrow(correlations)]
    value <- stats::runif(nrow(correlations), 0.25, 1)
    arrows <- data.frame(x1 = correlations$axe1, y1 = rep(0, nrow(correlations)),
                         x2 = correlations$axe1, y2 = sign * value,
                         names = correlations$names, labels = correlations$labels)
    p <- ggplot(arrows) + xlim(-1.1, 1.1) + ylim(-1.1, 1.1) +
      geom_hline(yintercept = 0, colour = "gray65") +
      geom_segment(aes_(x = ~x1, y = ~y1, xend = ~x2, yend = ~y2, colour = ~labels)) +
      geom_text(aes_(x = ~x2, y = ~y2, label = rownames(correlations), colour = ~labels), vjust = -.5, angle = 90, size=5) +
      geom_point(aes_(x = ~x2, y = 0)) +
      theme_bw() +  theme(axis.title.y = element_blank(),
                          axis.text.y  = element_blank(),
                          axis.ticks.y = element_blank(), legend.position = "none") + ggtitle(main) +
      labs(x = paste("correlation with", axes_label[1]))
  }
  p
}

get_ggplot_corr_circle <- function(correlations, axes_label, main) {
  p <- get_ggplot_corr(correlations, axes_label, main)
  if (length(axes_label) > 1 ) {
    ## add correlation circle
    corcir <- circle(c(0, 0), npoints = 100)
    p <- p + geom_path(data = corcir, aes_(x = ~x,y = ~y), colour = "gray65")
  }
  p
}

get_ggplot_corr_square <- function(correlations, axes_label, main) {
  p <- get_ggplot_corr(correlations, axes_label, main)
  if (length(axes_label) > 1 ) {
    ## add square of size 2 centered at the origin
    square <- data.frame(x = c(-1, -1, 1, 1, -1), y = c(-1, 1, 1, -1, -1))
    p <- p + geom_path(data = square, aes_(x = ~x,y = ~y), colour = "gray65")
  }
  p
}

circle <- function(center = c(0, 0), radius = 1, npoints = 100) {
  r = radius
  tt = seq(0, 2 * pi, length = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

GeomCircle <- ggplot2::ggproto("GeomCircle",
    ggplot2::Geom,
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

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

