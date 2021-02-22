#' @importFrom rlang .data
get_ggplot_ind_map <- function(scores, axes_label, main) {

  if (length(axes_label) > 1 ) {
    if (length(axes_label) >= 3) {
      message("Three axes or more provided. Using only the first two")
    }
    p <- ggplot(scores, aes(x = .data$a1, y = .data$a2, label = .data$names, colour = .data$labels)) +
      geom_hline(yintercept = 0, colour = "gray65") +
      geom_vline(xintercept = 0, colour = "gray65") +
      geom_text(alpha = 0.8, size = 4) +
      ggtitle(main) +
      theme_bw() +
      labs(x = axes_label[1], y = axes_label[2])
  } else {
    p <- ggplot(scores, aes(x = .data$a1, group = .data$labels, fill = .data$labels, colour = .data$labels)) +
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
      geom_segment(data = arrows, aes(x = .data$x1, y = .data$y1, xend = .data$x2, yend = .data$y2, colour = .data$labels)) +
      geom_text(data = correlations, aes(x = .data$axe1, y = .data$axe2, label = .data$names, colour = .data$labels), size=3) +
      theme_bw() +  theme(legend.position = "none") + ggtitle(main) + labs(x = axes_label[1], y = axes_label[2])
  } else {
    sign  <- rep(c(-1,1), each = ceiling(nrow(correlations)/2))[1:nrow(correlations)]
    value <- stats::runif(nrow(correlations), 0.25, 1)
    arrows <- data.frame(x1 = correlations$axe1, y1 = rep(0, nrow(correlations)),
                         x2 = correlations$axe1, y2 = sign * value,
                         names = correlations$names, labels = correlations$labels)
    p <- ggplot(arrows) + xlim(-1.1, 1.1) + ylim(-1.1, 1.1) +
      geom_hline(yintercept = 0, colour = "gray65") +
      geom_segment(aes(x = .data$x1, y = .data$y1, xend = .data$x2, yend = .data$y2, colour = .data$labels)) +
      geom_text(aes(x = .data$x2, y = .data$y2, label = rownames(correlations), colour = .data$labels), vjust = -.5, angle = 90, size=5) +
      geom_point(aes(x = .data$x2, y = 0)) +
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
    p <- p + geom_path(data = corcir, aes(x = .data$x,y = .data$y), colour = "gray65")
  }
  p
}

get_ggplot_corr_square <- function(correlations, axes_label, main) {
  p <- get_ggplot_corr(correlations, axes_label, main)
  if (length(axes_label) > 1 ) {
    ## add square of size 2 centered at the origin
    square <- data.frame(x = c(-1, -1, 1, 1, -1), y = c(-1, 1, 1, -1, -1))
    p <- p + geom_path(data = square, aes(x = .data$x,y = .data$y), colour = "gray65")
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

# courtesy of S. Donnet
plot_matrix = function(Mat, rowFG = "sample", colFG = "variable", clustering = NULL, log_scale = TRUE){

  n1 <- dim(Mat)[1]
  n2 <- dim(Mat)[2]
  u <- range(c(Mat))

  binary = FALSE
  val <- sort(unique(c(Mat)))
  if (setequal(val ,c(0,1))) {binary = TRUE}

  if (!is.null(clustering)) {
    oRow <- oCol <- order(clustering)
    uRow <- cumsum(table(clustering)) + 0.5
    uRow <- uRow[-length(uRow)]
    sepRow <- as.data.frame(uRow)
    Mat <- Mat[oRow, , drop = FALSE]
    names(sepRow) <- 'sep'
    sepRow = n1 - sepRow
  }

  index_row = rep(1:dim(Mat)[1],each = dim(Mat)[2])
  index_col = rep(1:dim(Mat)[2],dim(Mat)[1])

  melted_Mat =  data.frame(n1 - index_row , index_col)
  link = rep(-10,dim(Mat)[2]*dim(Mat)[1])
  for (k in 1:(dim(Mat)[2] * dim(Mat)[1])) {link[k] = Mat[index_row[k],index_col[k]]}
  melted_Mat$count = link
  if (binary){
    melted_Mat$count <- as.factor(melted_Mat$count)
  }
  colnames(melted_Mat) <- c('index_row', 'index_col', 'count')

  g <- ggplot(data = melted_Mat, aes(y = index_row, x = index_col, fill = count))
  g <- g + geom_tile()
  if (!binary & log_scale) {g <-  g +  scale_fill_viridis_b(limits = u, na.value = "transparent", trans = "log10")}
  if (!binary & !log_scale) {g <-  g +  scale_fill_viridis_b(limits = u, na.value = "transparent")}
  if (binary) {g <- g + scale_fill_manual(breaks = c("0", "1"),values = c("white", "black"),na.value = "transparent")}

  g <- g  +  scale_x_discrete(drop = FALSE, expand = expansion(0, 0)) + scale_y_discrete(drop = FALSE, expand = expansion(0, 0))
  g <- g + theme(axis.text.x = element_text(angle = 270, hjust = 0))
  g <- g +  labs(x = colFG, y = rowFG) +  theme(aspect.ratio = n1/n2)

  if (!is.null(clustering)){
    g <- g + geom_hline(data = sepRow, mapping = aes(yintercept = .data$sep),col = 'grey')
  }
  g
}
