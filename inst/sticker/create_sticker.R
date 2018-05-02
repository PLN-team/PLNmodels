# library(ggplot2)
#
#
# p1 <- ggplot(data = data.frame(x = c(0.1, 3)), aes(x)) +
#   stat_function(fun = dlnorm, n = 101, geom = "density", args = list(meanlog = 0, sdlog = 1), color = "#ece7f2", fill = "#2b8cbe", alpha=0.5) + theme_void()
# ggsave(filename = "wave.png", plot = p1, bg = 'transparent', width = 10, height = 8)
# p1


library(hexSticker)
imgurl <- "poissonnet_vague.png"
p <- sticker(imgurl, package="PLNmodels", p_size = 17, p_y = .5, s_x = 1, s_y = 1.25, s_width = .68,
        h_color = "#2b8cbe", h_fill = "white", h_size = 1.5, p_color = "#2b8cbe", filename = "PLNmodels.png")
