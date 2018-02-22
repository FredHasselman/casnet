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

## ----RNG5a, echo=FALSE---------------------------------------------------
cowplot::plot_grid(g_1sur, g_2sur, labels = c(paste0("ID 163 - P(X=",round(crqa_1$DET,2),") = ",p_1), paste0("ID 291 - P(X=",round(crqa_2$DET,2),") = ",p_2), align = "h"))

