library("MCMCtreeR")
library("ape")
library("phytools")
library("PerformanceAnalytics")

#main tree plot
phy_mcmc=readMCMCtree("schemeA.tre")
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

#time tree comparisons schemes

phy_mcmc=readMCMCtree("schemeA.tre")
phyA=phy_mcmc$apePhy
phy_mcmc=readMCMCtree("schemeB.tre")
phyB=phy_mcmc$apePhy
phy_mcmc=readMCMCtree("schemeC.tre")
phyC=phy_mcmc$apePhy
phy_mcmc=readMCMCtree("schemeD.tre")
phyD=phy_mcmc$apePhy
phy_mcmc=readMCMCtree("schemeRusso.tre")
phyR=phy_mcmc$apePhy

n_a=abs(node.depth.edgelength(phyA)[(phyA$Nnode+2):length(node.depth.edgelength(phyA))]-node.depth.edgelength(phyA)[1])
n_b=abs(node.depth.edgelength(phyB)[(phyB$Nnode+2):length(node.depth.edgelength(phyB))]-node.depth.edgelength(phyB)[1])
n_c=abs(node.depth.edgelength(phyC)[(phyC$Nnode+2):length(node.depth.edgelength(phyC))]-node.depth.edgelength(phyC)[1])
n_d=abs(node.depth.edgelength(phyD)[(phyD$Nnode+2):length(node.depth.edgelength(phyD))]-node.depth.edgelength(phyD)[1])
n_r=abs(node.depth.edgelength(phyR)[(phyR$Nnode+2):length(node.depth.edgelength(phyR))]-node.depth.edgelength(phyR)[1])

n_ages=data.frame(n_a=n_a,n_b=n_b,n_c=n_c,n_d=n_d,n_r=n_r)
n_ages=n_ages[order(n_ages$n_a),]
names(n_ages)=c("schemeA","schemeB","schemeC","schemeD","schemeRusso")


plot(log(n_ages[,2]),type="o",pch=16,col="red",xlab="node",ylab="log(age)")
lines(log(n_ages[,3]),type="o",pch=16,col="blue")
lines(log(n_ages[,4]),type="o",pch=16,col="purple")
lines(log(n_ages[,5]),type="o",pch=16,col="orange")
lines(log(n_ages[,1]),type="o",pch=16)
legend("bottomright",legend=c("scheme A","scheme B","scheme C","scheme D","scheme Russo"),col=c("black","red","blue","purple","orange"),lty=1,pch=16)
chart.Correlation(n_ages, histogram=F, pch="+",col="red")


#time tree comparisons between different gene samplings 

phy_mcmc=readMCMCtree("schemeA.tre")
phyA=phy_mcmc$apePhy
phy_mcmc=readMCMCtree("schemeA_10loci.tre")
phyB=phy_mcmc$apePhy
phy_mcmc=readMCMCtree("schemeA_100loci.tre")
phyC=phy_mcmc$apePhy

n_a=abs(node.depth.edgelength(phyA)[(phyA$Nnode+2):length(node.depth.edgelength(phyA))]-node.depth.edgelength(phyA)[1])
n_b=abs(node.depth.edgelength(phyB)[(phyB$Nnode+2):length(node.depth.edgelength(phyB))]-node.depth.edgelength(phyB)[1])
n_c=abs(node.depth.edgelength(phyC)[(phyC$Nnode+2):length(node.depth.edgelength(phyC))]-node.depth.edgelength(phyC)[1])

n_ages=data.frame(n_a=n_a,n_b=n_b,n_c=n_c)
n_ages=n_ages[order(n_ages$n_a),]
names(n_ages)=c("schemeA_1000","schemeA_10","schemeA_100")

plot(log(n_ages[,2]),type="o",pch=16,col="red",xlab="node",ylab="log(age)")
lines(log(n_ages[,3]),type="o",pch=16,col="blue")
lines(log(n_ages[,1]),type="o",pch=16)
chart.Correlation(n_ages, histogram=F, pch="+",col="blue")
legend("bottomright",legend=c("scheme A: 10 loci","scheme A: 100 loci","scheme A: 1000 loci"),col=c("red","blue","black"),lty=1,pch=16)
chart.Correlation(n_ages, histogram=F, pch="+",col="red")