.logfactorial <- function (n) { # Ramanujan's formula
  n <- n[n!=0]
  return(n*log(n) - n + log(8*n^3 + 4*n^2 + n + 1/30)/6 + log(pi)/2)
}

.trunc <- function(x, M = 300) {
  return(pmin(x, M))
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
multiplot <- function(..., legend=FALSE, plotlist=NULL, cols) {
  require(grid)

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
