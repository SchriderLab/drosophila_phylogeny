library("MCMCtreeR")
<<<<<<< HEAD
library(ape)
library(laser)
library(geiger)
library(TreeSim)
library(phytools)
library(ggplot2)
library(reshape2)
library(gridExtra)

phy_mcmc=readMCMCtree("/Users/anton/Dropbox/droso_phylogeny/schemeA.tre")
=======
library("ape")


phy_mcmc=readMCMCtree("schemeA.tre")
>>>>>>> 70d273d427c8591e90fbc13fa458fd33c999251c
phy=phy_mcmc$apePhy
tt=read.table("schemeA_mcmc.txt",header=T)

quartz(width=8,height=4)
quartz(width=8,height=4)
MCMC.tree.plot(rotate(keep.tip(phy,c(1:2,80,157:159)),9),analysis.type = "MCMCtree",MCMC.chain =tt[,1:6],plot.type = "distributions",density.col = adjustcolor( "navy", alpha.f = 0.5),density.border.col = "navy", lwd.bar = 3,scale.res = c("Period"), node.method = "bar",col.age = adjustcolor( "navy", alpha.f = 0.5), cex.labels = 0.5,cex.tips = 0.7,ladderize.tree = T,pos.age=-0.4,abs.age.lwd.ticks=0,relative.height=0.01,cex.age = 0.6,label.timescale.names = T,burn.in = 0,no.margin =T)
quartz.save("Fig1_droso.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12) 







quartz(width=5.5,height=12) 
MCMC.tree.plot(extract.clade(phy, 164),analysis.type = "MCMCtree",MCMC.chain=tt[,c(1,6:162)],plot.type = "distributions",density.col = adjustcolor( "navy", alpha.f = 0.5),density.border.col = "navy", lwd.bar = 3,scale.res = c("Period"), node.method = "bar",col.age = adjustcolor( "navy", alpha.f = 0.5), cex.labels = 0.5,cex.tips = 0.5,ladderize.tree = T,pos.age=-6.5,abs.age.lwd.ticks=0,relative.height=0.02,cex.age = 0.6,label.timescale.names = T,burn.in = 0,no.margin=T,label.offset = 0.5)
quartz.save("Fig1_drosoAB.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)

















quartz(width=8.5,height=7.7)
MCMC.tree.plot(extract.clade(phy, 165),analysis.type = "MCMCtree",MCMC.chain=tt[,c(1,7:83)],plot.type = "distributions",density.col = adjustcolor( "navy", alpha.f = 0.5),density.border.col = "navy", lwd.bar = 3,scale.res = c("Period","Epoch"), node.method = "bar",col.age = adjustcolor( "navy", alpha.f = 0.5), cex.labels = 0.5,cex.tips = 0.5,ladderize.tree = T,pos.age=-11.5,abs.age.lwd.ticks=0,relative.height=0.05,cex.age = 0.6,label.timescale.names = T,burn.in = 0,no.margin=T)
quartz.save("Fig1_drosoB.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)















phy.small = extract.clade(phy, 164)
tt=tt[,c(1,6:159)]
MCMC.tree.plot(extract.clade(phy, 164),analysis.type = "MCMCtree",MCMC.chain =tt,plot.type = "distributions",density.col = adjustcolor( "navy", alpha.f = 0.5),density.border.col = "navy", lwd.bar = 3,scale.res = c("Period","Epoch"), node.method = "bar",col.age = adjustcolor( "navy", alpha.f = 0.5), no.margin = T, cex.labels = 0.01,cex.tips = 0.6,ladderize.tree = F,pos.age=-7,abs.age.lwd.ticks=0,relative.height=0.05,cex.age = 0.6)