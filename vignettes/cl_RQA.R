## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = FALSE,
  comment = ">",
  fig.align = "center",
  fig.width = 7 
)

## ----RNG1, echo=TRUE, message=FALSE, warning=FALSE, paged.print=FALSE----
library(casnet)
library(ggplot2)

# Load the random number sequence data from Oomens et al. (2015)
data(RNG)

# Select a subject
IDs <- RNG$ID%in%c(163,291)

# Look at the sequence
ggplot(RNG[IDs,],aes(x=time,y=number,group=ID)) +
  geom_line(aes(colour=ID))+
  facet_grid(~ID) +
  scale_y_continuous(breaks = 1:9) +
  ggtitle("Which of these number sequences is more 'random'?") +
  theme_bw()

## ----RNG2, message=FALSE, warning=FALSE, cache=TRUE, collapse=FALSE, include=FALSE, echo=TRUE, paged.print=FALSE----
# Run the RQA analysis
y_1  <- RNG$number[RNG$ID==163]
y_2 <-  RNG$number[RNG$ID==291]
crqa_1 <- crqa_cl(y1 = y_1, emDim  = 1, emLag = 1, emRad= 1)
crqa_2 <- crqa_cl(y1 = y_2, emDim  = 1, emLag = 1, emRad= 1)

## Plot the recurrence matrix
# Get the matrix
rp_1 <- recmat(y1=y_1, y2=y_1, emDim = 1, emLag = 1)
rp_2 <- recmat(y1=y_2, y2=y_2, emDim = 1, emLag = 1)

# Turn it into a recurrence matrix
rp_1 <- di2bi(rp_1,radius = 1)
rp_2 <- di2bi(rp_2,radius = 1)

# Get the plots
g_1 <- recmat_plot(rp_1, doPlot=FALSE)
g_2 <- recmat_plot(rp_2, doPlot=FALSE)

## ----RNG3, echo=TRUE, message=FALSE, warning=FALSE, cache=TRUE, paged.print=FALSE----
library(cowplot)
# The recurrence plots
plot_grid(g_1, g_2, labels = c("ID 163", "ID 291"), align = "h")

# The RQA measures
cbind.data.frame(subj163=t(crqa_1), subj291=t(crqa_2))

## ----RNG4, echo=TRUE, message=FALSE, warning=FALSE, cache=TRUE, paged.print=FALSE----
# Reproduce the same randomisation
set.seed(123456789)

# Randomise the number sequences
y_1rnd <- y_1[sample(1:NROW(y_1),size = NROW(y_1))]
y_2rnd <- y_2[sample(1:NROW(y_2),size = NROW(y_2))]

# Calculate RQA measures
crqa_1rnd <- crqa_cl(y1 = y_1rnd, emDim  = 1, emLag = 1, emRad= 1)
crqa_2rnd <- crqa_cl(y1 = y_2rnd, emDim  = 1, emLag = 1, emRad= 1)

# Get the RPs in one statement
g_1rnd <- recmat_plot(di2bi(recmat(y1=y_1rnd, y2=y_1rnd, emDim = 1, emLag = 1), radius=1), doPlot=FALSE)
g_2rnd <- recmat_plot(di2bi(recmat(y1=y_2rnd, y2=y_2rnd, emDim = 1, emLag = 1), radius=1), doPlot=FALSE)

# Display recurrence plots
cowplot::plot_grid(g_1rnd, g_2rnd, labels = c("ID 163 shuffled", "ID 291 shuffled"), align = "h")

# Display the RQA measures
cbind.data.frame(subj163=t(crqa_1), subj163rnd=t(crqa_1rnd), subj291=t(crqa_2),  subj291rnd=t(crqa_2rnd))

## ----RNG5, echo=TRUE, message=FALSE, warning=FALSE, cache=TRUE, paged.print=FALSE----
library(plyr)
library(dplyr)

set.seed(123456789)

y_1rnd_sur <- ldply(1:99, function(s) y_1[sample(1:NROW(y_1),size = NROW(y_1))])
y_2rnd_sur <- ldply(1:99, function(s) y_2[sample(1:NROW(y_2),size = NROW(y_2))])

crqa_1rnd_sur <- ldply(seq_along(y_1rnd_sur$V1), function(r) crqa_cl(y1 = as.numeric(y_1rnd_sur[r,]), emDim  = 1, emLag = 1, emRad= 1))
crqa_1rnd_sur[NROW(crqa_1rnd_sur)+1,] <- crqa_1
crqa_2rnd_sur <- ldply(seq_along(y_2rnd_sur$V1), function(r) crqa_cl(y1 = as.numeric(y_2rnd_sur[r,]), emDim  = 1, emLag = 1, emRad= 1))
crqa_2rnd_sur[NROW(crqa_2rnd_sur)+1,] <- crqa_2

## ----RNG5a, echo=TRUE, message=FALSE, warning=FALSE, cache=TRUE, paged.print=FALSE----
# Get point estimates for p-values based on rank of observation (discrete distribution)
#   99 = (1 / alpha) - 1
# 99+1 = (1 / alpha)
alpha = 1/100

p_1 <- point_p(surrogateValues = crqa_1rnd_sur$DET, observedValue = crqa_1$DET, alpha = alpha,measureName = "DET")
p_2 <-point_p(surrogateValues = crqa_2rnd_sur$DET, observedValue = crqa_2$DET, alpha = alpha, measureName = "DET")

cowplot::plot_grid(p_1$plot, p_2$plot, labels = c("ID 163","ID 291"), align = "h")

## ----RNG6, echo=TRUE, message=FALSE, warning=FALSE, cache=TRUE, paged.print=FALSE----
# Get point estimates for p-values based on rank of observation (discrete distribution)
#   99 = (1 / alpha) - 1
# 99+1 = (1 / alpha)
alpha = 1/100

p_1 <- point_p(surrogateValues = crqa_1rnd_sur$LAM, observedValue = crqa_1$LAM, alpha = alpha)
p_2 <- point_p(surrogateValues = crqa_2rnd_sur$LAM, observedValue = crqa_2$LAM, alpha = alpha)

cowplot::plot_grid(p_1$plot, p_2$plot, labels = c("ID 163","ID 291"), align = "h")

