library(devtools)
library(plyr)
library(tidyverse)


# importFrom("grDevices", "Type1Font", "col2rgb", "colorRampPalette", "dev.off", "pdfFonts", "postscriptFonts", "rgb", "svg")
# importFrom("graphics", "Axis", "hist", "legend", "lines", "par", "plot", "plot.new")
# importFrom("stats", "coef", "complete.cases", "frequency", "is.ts", "lm", "mad", "na.omit", "predict", "quantile", "runif", "sd", "start", "time", "ts", "tsp<-", "var")
# importFrom("utils", "write.table")

# Add packages to NAMESPACE ----
l_ply(sort(c("dplyr","fractal","ggplot2","gridExtra","ifultools","igraph","qgraph","lattice","latticeExtra","rio","Matrix","nonlinearTseries","plyr","pracma","proxy","scales","tidyr","xts","zoo","cowplot")), function(p) devtools::use_package(p))

# Add dependencies
l_ply(sort(c("grDevices","graphics","stats","utils")), function(p) devtools::use_package(p,"Suggests"))


# Install locally

install(build_vignettes = TRUE)

# BUILD tarballs / zip ---

path = "~/Documents/GitHub/casnet/pkg/"

# Source
build(binary=FALSE, vignettes = TRUE, manual = TRUE, path = path)

# Binary
build(binary=TRUE,  vignettes = TRUE, manual = TRUE, path = path)

# Windows
build_win()


# Build DATA sets ---

# Oomens et al. 2015 ---
df <- openxlsx::read.xlsx("~/Documents/GitHub/casnet/build/package_data/Oomens2015.xlsx",sheet=2,colNames=F)
df <- df[df[,3]==1,]
df <- df[,-c(2,3)]
colnames(df) <- c("ID",paste0("i_",1:100))
RNG <- df
RNG <-gather(RNG, key=time, value=number,-ID)
RNG$ID <- factor(RNG$ID)
RNG$time <- as.numeric(gsub("(i[_])+","",RNG$time))
devtools::use_data(RNG, overwrite = TRUE)




