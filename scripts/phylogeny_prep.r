library("MCMCtreeR")
library(ape)
library(laser)
library(geiger)
library(TreeSim)
library(phytools)
library(ggplot2)
library(reshape2)
library(gridExtra)

phy_mcmc=readMCMCtree("/Users/anton/Dropbox/droso_phylogeny/schemeA.tre")
phy=phy_mcmc$apePhy
tt=read.table("/Users/Anton/Downloads/schemeA_mcmc.txt",header=T)
tt=tt[,1:159]


phy.small = extract.clade(phy, 164)
tt=tt[,c(1,6:159)]
MCMC.tree.plot(phy.small,analysis.type = "MCMCtree",MCMC.chain =tt,plot.type = "distributions",density.col = adjustcolor( "navy", alpha.f = 0.5),density.border.col = "navy", lwd.bar = 3,scale.res = c("Period","Epoch"), node.method = "bar",col.age = adjustcolor( "navy", alpha.f = 0.5), no.margin = T, cex.labels = 0.01,cex.tips = 0.6,ladderize.tree = F,pos.age=-7,abs.age.lwd.ticks=0,relative.height=0.05,cex.age = 0.6)