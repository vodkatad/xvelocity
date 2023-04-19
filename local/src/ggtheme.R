library(ggplot2)

textSize <- 6
largerSize <- textSize + 1

gtheme <- theme_bw() +
  theme(
    text = element_text(size = textSize, family='sans'),
    axis.title = element_text(size = largerSize),
    axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
    axis.text.y = element_text(size = textSize, color="black"),
    plot.title = element_text(size = largerSize, hjust = 0.5),
    legend.title = element_text(size=largerSize),
    legend.text = element_text(size=textSize)
  )