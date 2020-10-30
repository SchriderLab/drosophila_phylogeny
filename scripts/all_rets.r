library("ape")

tt=read.tree("MLrooted.tre")
c2=extract.clade(tt,245)
c2$edge.length=NULL
c2n=evonet(c2,20,13)
write.evonet(c2n)
nodelabels()