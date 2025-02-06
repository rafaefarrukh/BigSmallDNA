###################################################################################################

# BigSmallDNA - Time Complexity Comparison

###################################################################################################

# libraries ----
library(ggplot2) # visualize
library(gridExtra) # place in grid

# Data ----

# import data from example
nd <- read.csv(paste0(wd, "/data/nd.csv"))

# functions
normalize <- function(x) {(x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))}
m1 <- function(n,l) {(l * n^2 * 2^n)} # ideal approach
m2 <- function(n,l) {(l * n^4 * 4)} # BigSmallDNA approach

# sample size for m1 and m2 calculations
n <- seq(2,10000,100); l <- seq(10,10000,100)

# data
m1 <- data.frame(type = "Ideal", n = rep(n,length(l)), l = rep(l, each=length(n)), t = normalize(as.numeric(gsub(Inf, NaN, m1(rep(n,length(l)), rep(l, each=length(n)))))))
m2 <- data.frame(type = "BigSmallDNA", n = rep(n,length(l)), l = rep(l, each=length(n)), t = normalize(m2(rep(n,length(l)), rep(l, each=length(n)))))
m3 <- data.frame(type = "Example Data", n = nd$old.seq, l = nd$length, t = normalize(nd$time))

data <- data.frame(rbind(m1, m2, m3))

# Graphs ----

graphs <- grid.arrange(
    ggplot(data, aes(x=n, y=l, color=t)) +
        geom_point() +
        scale_color_viridis_c() +
        guides(fill = guide_colourbar(position = "left", barwidth = 0.5, label = F)) +
        facet_wrap(vars(type)) +
        labs(x = "n", y = "L", color = "time")
    ,
    
    ggplot(data, aes(x=n, y=t, color=l)) +
        geom_point() +
        scale_color_viridis_c() +
        guides(fill = guide_colourbar(position = "left", barwidth = 0.5, label = F)) +
        facet_wrap(vars(type)) +
        labs(x = "n", y = "time", color = "L")
    ,
    
    ggplot(data, aes(x=l, y=t, color=n)) +
        geom_point() +
        scale_color_viridis_c() +
        guides(fill = guide_colourbar(position = "left", barwidth = 0.5, label = F)) +
        facet_wrap(vars(type)) +
        labs(x = "L", y = "time", color = "n")
)

ggsave(plot = graphs, filename = "TimeComplexity/Figure3.png",
       units = "cm", width = 40, height = 20, dpi=300)

###################################################################################################