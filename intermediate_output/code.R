# code for plotting PR curves
library(PRROC)
load("TSG.rdata")
for (i in 1:length(TSG.predictions)) {
  prc <- pr.curve(scores.class0=TSG.predictions[[i]], weights.class0=y_TSG, curve=T)
  plot(prc, main=names(TSG.predictions)[i])
}

load("OG.rdata")
for (i in 1:length(OG.predictions)) {
  prc <- pr.curve(scores.class0=OG.predictions[[i]], weights.class0=y_OG, curve=T)
  plot(prc, main=names(OG.predictions)[i])
}