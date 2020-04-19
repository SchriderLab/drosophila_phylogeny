library("ape")
library('ggplot2')
library('pals')
library('reshape2')
library('dplyr')
library('gridExtra')

#Identify introgressing taxa pair from QuiBL result
get_intropair=function(m)
{
    l=m[,c("P1","P2","P3")]!=m[,c("outgroup","outgroup","outgroup")]
    sp_list=m[,c("P1","P2","P3")]
    pair_v=rbind()
    for (i in 1:nrow(l))        
    {
        pair_intro=as.character(sp_list[i,][l[i,]])
        pair_v=rbind(pair_v,pair_intro)
    }
    pair_v=data.frame(pair_v)
    names(pair_v)=c("i1","i2")
    rownames(pair_v) = c()
    return(pair_v)
}



tt=read.tree("MLrooted.tre")


C1=extract.clade(tt,168)$tip.label
C2=extract.clade(tt,185)$tip.label
C3=extract.clade(tt,204)$tip.label
C4=extract.clade(tt,226)$tip.label
C5=extract.clade(tt,245)$tip.label
C6=extract.clade(tt,256)$tip.label
C7=extract.clade(tt,265)$tip.label
C8=extract.clade(tt,289)$tip.label
C9=extract.clade(tt,307)$tip.label




total_q=read.csv("/Users/Anton/Downloads/droso_quibl_results.txt",header=F,stringsAsFactors=F)
names(total_q)=c("clade","id","P1","P2","P3","outgroup","Com1","Com2","mixprop1","mixprop2","lambda2Dist","lambda1Dist","BIC2Dist","BIC1Dist","count")
total_q$clade=unlist(lapply(strsplit(as.character(total_q$clade), "_"),"[",1))
total_q=total_q[complete.cases(total_q), ] 
total_q$BICdiff = total_q$BIC2-total_q$BIC1
total_q$triplet=as.character(apply(total_q[,c("P1","P2","P3")],1,paste,collapse="_"))
total_q=cbind(total_q,get_intropair(total_q))


#total_q=total_q %>% distinct(triplet,C2 ,  mixprop1,  mixprop2, lambda2Dist, BIC2Dist ,   BIC1Dist, count ,.keep_all = T)
total_q_min=data.frame(total_q %>%  group_by(triplet) %>% filter(count!=max(count)))
total_q_max=data.frame(total_q %>%  group_by(triplet) %>% filter(count==max(count)))
total_q_min$common=FALSE
total_q_max$common=TRUE
total_q=rbind(total_q_min,total_q_max)
total_q$sig=total_q$BICdiff < -10
total_q$genus="Drosophila"
total_q$type=ifelse(total_q$common==TRUE,"Concordant",ifelse(total_q$sig==FALSE & total_q$common==FALSE,"ILS","Introgression"))

total_mixprop=melt(total_q[total_q$sig==TRUE & total_q$common==FALSE,c("mixprop2","clade","genus")],value.name="taxon",id=c("mixprop2"))

g1=ggplot(total_mixprop, aes(x=taxon, y=mixprop2))+geom_violin(fill='salmon')+facet_grid(~variable,scales = "free", space = "free")+stat_summary(fun.y=median, geom="point", size=2, color="black")+geom_boxplot(width=0.01,outlier.size=-1)+ylab(expression(pi[2]))+xlab("")

total_p=melt(total_q[ ,c("type","clade","genus")],id="type",value.name="taxon")


g2=ggplot(total_p, aes(x=taxon, y=..count../sum(..count..),fill=type))+geom_bar(position="fill")+facet_grid(~variable,scales = "free", space = "free")+geom_text(aes(label=..count..),stat="count",position=position_fill(vjust=0.5))+theme(legend.position=c("top") ,legend.direction="horizontal")+ylab("Proportion")+xlab("")+scale_fill_manual(values=c("gray48","gold", "salmon"),name="")

grid.arrange(g1,g2,nrow=2)







get_intomatrix=function(taxa,dat)
{
    m=matrix(0, nrow=length(taxa),ncol=length(taxa)) 
    m=data.frame(m)
    colnames(m)=taxa
    rownames(m)=taxa
    dat_int=dat[dat$sig==TRUE & dat$common==FALSE & (dat$i1 %in% taxa | dat$i2 %in% taxa ),c("i1","i2")]
    for (i in 1:nrow(dat_int))
    {
        a=dat_int[i,"i1"]
        b=dat_int[i,"i2"]
        m[rownames(m)==a,colnames(m)==b]=m[rownames(m)==a,colnames(m)==b]+1
    }    
    m=m+t(m)
    m[upper.tri(m)]=NA
    melted_m=melt(as.matrix(m),na.rm = TRUE)
    melted_m$value=melted_m$value/sum(melted_m$value)
    g1=ggplot(melted_m, aes(Var2, Var1, fill = value))+geom_tile(color = "white")+scale_fill_gradientn(colors=viridis(100),space = "Lab", name="Introgression frequency") +theme_minimal()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+coord_fixed()+labs(x="",y="")
    return(g1)
}    




         
        