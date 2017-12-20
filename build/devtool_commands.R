library(tidyverse)

l_ply(sort(c("dplyr","fractal","ggplot2","gridExtra","ifultools","igraph","qgraph","lattice","latticeExtra","Matrix","nonlinearTseries","plyr","pracma","proxy","scales","tidyr","xts","zoo","cowplot")), function(p) devtools::use_package(p))

devtools::use_vignette("commandline_RQA")

df <- openxlsx::read.xlsx("/Users/Fred/Documents/Projects/RNGproject/2014 - Oomens/Ruw_totaal.xlsx",sheet=2,colNames=F)
df <- df[df[,3]==1,]
df <- df[,-c(2,3)]
colnames(df) <- c("ID",paste0("i_",1:100))
RNG <- df
RNG <-gather(RNG, key=time, value=number,-ID)
RNG$ID <- factor(RNG$ID)
RNG$time <- as.numeric(gsub("(i[_])+","",RNG$time))
devtools::use_data(RNG, overwrite = TRUE)


tmpd  <- tempdir()
tmpf1 <- tempfile(tmpdir = tmpd, fileext = ".csv")
write.table(as.data.frame(y1), tmpf1, col.names = FALSE, row.names = FALSE)

devtools::install(build_vignettes = TRUE)
vignette("cl_RQA")


build_win(args = c("--no-examples"))



devtools::RCMD(paste0(getOption("casnet.rp_prefix"),"rp"), options = opts, path = normalizePath(path.expand(path_to_rp), mustWork = FALSE), quiet = FALSE)

opts ='-i /var/folders/zr/n8mgv2nj5rz1qq04xsj_x_c80000gp/T//RtmpMTSUw8/file51172e578ad.csv -r /var/folders/zr/n8mgv2nj5rz1qq04xsj_x_c80000gp/T//RtmpMTSUw8/RQAplot_1.txt -o /var/folders/zr/n8mgv2nj5rz1qq04xsj_x_c80000gp/T//RtmpMTSUw8/RQAmeasures_1.txt -p /var/folders/zr/n8mgv2nj5rz1qq04xsj_x_c80000gp/T//RtmpMTSUw8/RQAhist_1.txt -m 1 -t 1 -e 1 -l 2 -v 2 -w 0 -n EUCLIDEAN -s'


y1 <- rnorm(100)
rio::export(data_frame(y1),'/Users/Fred/y1.csv')
rng_out <- crqa_fast(y1 = y1, eDim  = 1, eLag = 1, eRad= 1)

rp<-recmat(y1=y1,emDim = 1,emLag = 1)

# rp[rp >  0] <- NA
# rp[rp == 0] <- 1
rp<-di2bi(rp,radius = 1)
recmat_plot(rp)

attributes(rp)





surr <- surrogate(beamchaos, method="dh")

## print the result
print(surr)

## plot and compare various statistics of the
## surrogate and original time series
plot(surr, type="time")
plot(surr, type="sdf")
plot(surr, type="lag")
plot(surr, type="pdf")

## create comparison time history
plot(surr, show="both", type="time")



library(RandomFields)
# 1d time series
n <- 256
rf <- RFsimulate(x = c(0,1, 1/n), model = "stable",
              grid = TRUE, gridtriple = TRUE,
              param = c(mean=0, variance=1, nugget=0, scale=1, kappa=1))
par(mfrow=c(4,2))




## first let us look at the list of implemented models
RFgetModelNames(type="positive definite", domain="single variable",
                iso="isotropic")

## our choice is the exponential model;
## the model includes nugget effect and the mean:
model <- RMexp(var=5, scale=10) + # with variance 4 and scale 10
  RMnugget(var=1) + # nugget
  RMtrend(mean=0.5) # and mean

## define the locations:
from <- 0
to <- 20
x.seq <- seq(from, to, length=200)
y.seq <- seq(from, to, length=200)

simu <- RFsimulate(model, x=x.seq, y=y.seq)
plot(simu)
