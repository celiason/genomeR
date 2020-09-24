# 

library(ape)
trees <- lapply(list.files("output", pattern=".phy", full=T), read.tree)

plot(trees[[1]])
plot(trees[[2]])


library(phytools)

trees[[1]]$tip.label <- c('ceyArg','ceyMar','corVin','corMad','todChl')
trees[[2]]$tip.label <- c('ceyArg','ceyMar','todChl','corVin','corMad')

# Fish
# SRR9614116 = argentatus
# SRR9644504 = vintsiodes

# No fish
# SRR9642614 = margarethae
# SRR9644342 = madagascariensis

assoc <- cbind(trees[[1]]$tip, trees[[1]]$tip)

obj<-cophylo(trees[[1]], trees[[2]], rotate=F)
plot(obj, assoc=assoc, link.type="curved")

cophyloplot(trees[[1]], trees[[2]], assoc=assoc, link='curved')
