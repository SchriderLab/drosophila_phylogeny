library("ape")
library('ggplot2')
library('pals')
library('reshape2')
library('dplyr')
library('gridExtra')
library("svMisc")

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

# Plot matrix from QuiBL results 
get_intomatrix_quibl=function(taxa,dat)
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
    #Normalazing by the triplets that contain an introgressing pair (= Ntaxa-2)
    melted_m$value=melted_m$value/(length(taxa)-2)
    g1=ggplot(melted_m, aes(Var2, Var1, fill = value))+geom_tile(color = "white")+scale_fill_gradientn(colors=jet(100),space = "Lab", name="N significant triplets/\nN triplets cointaining introgressing taxa") +theme_minimal()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+coord_fixed()+labs(x="",y="")
    return(g1)
}    



# Plot matrix from QuiBL counts results 
get_intomatrix_quibl_c=function(taxa,dat)
{
    m=matrix(0, nrow=length(taxa),ncol=length(taxa)) 
    m=data.frame(m)
    colnames(m)=taxa
    rownames(m)=taxa
    dat_int=dat[dat$sig==TRUE & dat$common==FALSE & (dat$i1 %in% taxa | dat$i2 %in% taxa ),c("i1","i2","count")]
    for (i in 1:nrow(dat_int))
    {
        a=dat_int[i,"i1"]
        b=dat_int[i,"i2"]
        m[rownames(m)==a,colnames(m)==b]=m[rownames(m)==a,colnames(m)==b]+dat_int[i,"count"]
    }    
    m=m+t(m)
    m[upper.tri(m)]=NA
    melted_m=melt(as.matrix(m),na.rm = TRUE)
    #Normalazing by the triplets that contain an introgressing pair (= Ntaxa-2)
    melted_m$value=melted_m$value/(length(taxa)-2)
    g1=ggplot(melted_m, aes(Var2, Var1, fill = value))+geom_tile(color = "white")+scale_fill_gradientn(colors=jet(100),space = "Lab", name="Introgression discordant counts") +theme_minimal()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+coord_fixed()+labs(x="",y="")
    return(g1)
}    

# Plot matrix from QuiBL chi-square results 
get_intomatrix_quibl_s=function(taxa,dat)
{
    m=matrix(0, nrow=length(taxa),ncol=length(taxa)) 
    m=data.frame(m)
    colnames(m)=taxa
    rownames(m)=taxa
    dat_int=dat[dat$Qsig==TRUE & dat$common==FALSE & (dat$i1 %in% taxa | dat$i2 %in% taxa ),c("i1","i2")]
    for (i in 1:nrow(dat_int))
    {
        a=dat_int[i,"i1"]
        b=dat_int[i,"i2"]
        m[rownames(m)==a,colnames(m)==b]=m[rownames(m)==a,colnames(m)==b]+1
    }    
    m=m+t(m)
    m[upper.tri(m)]=NA
    melted_m=melt(as.matrix(m),na.rm = TRUE)
    #Normalazing by the triplets that contain an introgressing pair (= Ntaxa-2)
    melted_m$value=melted_m$value/(length(taxa)-2)
    g1=ggplot(melted_m, aes(Var2, Var1, fill = value))+geom_tile(color = "white")+scale_fill_gradientn(colors=jet(100),space = "Lab", name="N significant triplets/\nN triplets cointaining introgresing taxa") +theme_minimal()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+coord_fixed()+labs(x="",y="")
    return(g1)
}    



#Identify introgressing taxa pair from Dfoil result
get_intropair_dfoil=function(m)
{
    pair_v=rbind()
    for (i in 1:nrow(m))
    {
        progress(i,nrow(m))
        if(m[i,"introgression"]!="none" & m[i,"introgression"]!="123" & m[i,"introgression"]!="124")
        {
            pair=as.character(unlist(m[i,c("P1","P2","P3","P4")][sort(as.numeric(unlist(strsplit(as.character(m[i,"introgression"]),split=""))))]))
            pair_v=rbind(pair_v,pair)
        }else{
            pair=c("none","none")
            pair_v=rbind(pair_v,pair)
        }    
    }
    pair_v=data.frame(pair_v)
    names(pair_v)=c("i1","i2")
    rownames(pair_v) = c()
    return(pair_v)
}

# Plot matrix from Dfoil results 
get_intomatrix_dfoil=function(taxa,dat)
{
    m=matrix(0, nrow=length(taxa),ncol=length(taxa)) 
    m=data.frame(m)
    colnames(m)=taxa
    rownames(m)=taxa
    total_pairs=rbind(as.matrix(dat[,c("P1","P4")]),as.matrix(dat[,c("P1","P3")]),as.matrix(dat[,c("P2","P3")]),as.matrix(dat[,c("P2","P4")]))
    total_pairs=data.frame(total_pairs)
    names(total_pairs)=c("i1","i2")
    total_pairs=apply(total_pairs,2,as.character)
    dat_int=dat[dat$introgression!="none" & (dat$i1 %in% taxa | dat$i2 %in% taxa),c("i1","i2")]
    if (nrow(dat_int)!=0)
    {
        for (i in 1:nrow(dat_int))
        {
            a=dat_int[i,"i1"]
            b=dat_int[i,"i2"]
            m[rownames(m)==a,colnames(m)==b]=m[rownames(m)==a,colnames(m)==b]+1
        }    
    }    
    
    m=m+t(m)
    m[upper.tri(m)]=NA
    melted_m=melt(as.matrix(m),na.rm = TRUE)
    for (n in 1:nrow(melted_m))
    {
        if (melted_m[n,"value"]!=0)
        {
            f=nrow(total_pairs[(total_pairs[,"i1"]==as.character(melted_m[n,"Var1"]) &  total_pairs[,"i2"]==as.character(melted_m[n,"Var2"])) | (total_pairs[,"i2"]==as.character(melted_m[n,"Var1"]) &  total_pairs[,"i1"]==as.character(melted_m[n,"Var2"])),])
            
            melted_m[n,"value"]=melted_m[n,"value"]/f
        }    
        
    }
    
    g1=ggplot(melted_m, aes(Var2, Var1, fill = value))+geom_tile(color = "white")+scale_fill_gradientn(colors=jet(100),space = "Lab", name="N significant quintets/\nN quintets cointaining introgressing taxa") +theme_minimal()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+coord_fixed()+labs(x="",y="")
    return(g1)
}    

#Identify introgressing taxa pair from BLT result
get_intropair_blt=function(m)
{
    pair_v=rbind()
    for (i in 1:nrow(m))
    {
        progress(i,nrow(m))
        if(m[i,"Pvalue"]<0.05)
        {
            pair=as.character(m[i,c("P1out","P2out","P3out")][as.character(m[i,c("P1out","P2out","P3out")])!=as.character(m[i,c("P1out","P2out","P3out")][1:3!=which.max(m[i,c("CountP1","CountP2","CountP3")])][which.min(m[i,c("meanT_discord1","meanT_discord2")])])])
            pair_v=rbind(pair_v,pair)
        }else{
            pair=c("none","none")
            pair_v=rbind(pair_v,pair)
        }    
    }
    pair_v=data.frame(pair_v)
    names(pair_v)=c("i1","i2")
    rownames(pair_v) = c()
    return(pair_v)
}

# Plot matrix from BLT results 
get_intomatrix_blt=function(taxa,dat)
{
    m=matrix(0, nrow=length(taxa),ncol=length(taxa)) 
    m=data.frame(m)
    colnames(m)=taxa
    rownames(m)=taxa
    dat_int=dat[dat$i1 %in% taxa | dat$i2 %in% taxa ,c("i1","i2")]
    for (i in 1:nrow(dat_int))
    {
        a=dat_int[i,"i1"]
        b=dat_int[i,"i2"]
        m[rownames(m)==a,colnames(m)==b]=m[rownames(m)==a,colnames(m)==b]+1
    }    
    m=m+t(m)
    m[upper.tri(m)]=NA
    melted_m=melt(as.matrix(m),na.rm = TRUE)
    #Normalazing by the triplets that contain an introgressing pair (= Ntaxa-2)
    melted_m$value=melted_m$value/(length(taxa)-2)
    g1=ggplot(melted_m, aes(Var2, Var1, fill = value))+geom_tile(color = "white")+scale_fill_gradientn(colors=jet(100),space = "Lab", name="N significant triplets/\nN triplets cointaining introgressing taxa") +theme_minimal()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+coord_fixed()+labs(x="",y="")
    return(g1)
}    




#################################################################### QuiBL ##############################################################

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
#total_q=total_q[complete.cases(total_q), ] 
total_q$BICdiff = total_q$BIC2-total_q$BIC1
total_q$triplet=as.character(apply(total_q[,c("P1","P2","P3")],1,paste,collapse="_"))
total_q=cbind(total_q,get_intropair(total_q))
total_q$common=FALSE

chiPtri=c()
for (trip in unique(as.character(total_q$triplet)))
{
    p_v=chisq.test(subset(total_q,triplet==trip)$count)$p.value
    chiPtri=c(chiPtri,p_v)
}    
chiPtri=p.adjust(chiPtri,method="fdr")
total_q$Qtri=rep(chiPtri,each=3)







for (trip in unique(as.character(total_q$triplet)))
{
  maxVal=max(subset(total_q,triplet==trip)$count)
  # Handle ties
  if(nrow(total_q[which(total_q$triplet==trip & total_q$count==maxVal),])>1)
  {
     total_q[which(total_q$triplet==trip & total_q$count==maxVal),][1,]$common=TRUE 
  }else{    
    total_q[which(total_q$triplet==trip & total_q$count==maxVal),]$common=TRUE
  }    
}

total_min=total_q[total_q$common==F,]
total_max=total_q[total_q$common==T,]
chiP=c()
for (trip in unique(as.character(total_min$triplet)))
{
    co=subset(total_min,triplet==trip)$count
    # Handle cases when both incongruent topologies have 0 counts (very few cases)
    if (all(co==0))
    {
        p_v=chisq.test(co+1)$p.value
        chiP=c(chiP,p_v)
    }else{
        
        p_v=chisq.test(subset(total_min,triplet==trip)$count)$p.value
        chiP=c(chiP,p_v)
    }
}
       
chiP=p.adjust(chiP,method="fdr")
total_min$Q=rep(chiP,each=2)
total_max$Q=0
total_q=rbind(total_min,total_max)
total_q=total_q[complete.cases(total_q), ] 
total_q$BICdiff = total_q$BIC2-total_q$BIC1
total_q$Qsig = total_q$Q<0.05
total_q$Qtrisig = total_q$Qtri<0.05
total_q$sig=total_q$BICdiff < -10    





#total_q=total_q %>% distinct(triplet,C2 ,  mixprop1,  mixprop2, lambda2Dist, BIC2Dist ,   BIC1Dist, count ,.keep_all = T)
#total_q_min=data.frame(total_q %>%  group_by(triplet) %>% filter(count!=max(count)))
#total_q_max=data.frame(total_q %>%  group_by(triplet) %>% filter(count==max(count)))

total_q$genus="Drosophila"
total_q$type=ifelse(total_q$common==TRUE & total_q$sig==TRUE,"Concordant",ifelse(total_q$common==TRUE & total_q$sig==FALSE,"Extreme ILS",ifelse(total_q$sig==FALSE & total_q$common==FALSE,"ILS","Introgression+ILS")))
total_q$Qtype=ifelse(total_q$common==TRUE & total_q$Qsig==TRUE,"Concordant",ifelse(total_q$common==TRUE & total_q$Qtrisig==FALSE,"Extreme ILS",ifelse(total_q$Qsig==FALSE & total_q$common==FALSE ,"ILS","Introgression+ILS")))


####Plotting
total_mixprop=melt(total_q[total_q$sig==TRUE & total_q$common==FALSE,c("mixprop2","clade","genus")],value.name="taxon",id=c("mixprop2"))
g1=ggplot(total_mixprop, aes(x=taxon, y=mixprop2))+geom_violin(fill='salmon')+facet_grid(~variable,scales = "free", space = "free")+stat_summary(fun.y=median, geom="point", size=2, color="black")+geom_boxplot(width=0.01,outlier.size=-1)+ylab(expression(pi[2]))+xlab("")

total_mixprop=melt(total_q[total_q$Qsig==TRUE & total_q$common==FALSE,c("mixprop2","clade","genus")],value.name="taxon",id=c("mixprop2"))
gq1=ggplot(total_mixprop, aes(x=taxon, y=mixprop2))+geom_violin(fill='salmon')+facet_grid(~variable,scales = "free", space = "free")+stat_summary(fun.y=median, geom="point", size=2, color="black")+geom_boxplot(width=0.01,outlier.size=-1)+ylab(expression(pi[2]))+xlab("")


total_p=melt(total_q[ ,c("type","clade","genus")],id="type",value.name="taxon")
g2=ggplot(total_p, aes(x=taxon, y=..count../sum(..count..),fill=type))+geom_bar(position="fill")+facet_grid(~variable,scales = "free", space = "free")+geom_text(aes(label=..count..),stat="count",position=position_fill(vjust=0.5))+theme(legend.position=c("top") ,legend.direction="horizontal")+ylab("Proportion")+xlab("")+scale_fill_manual(values=c("gray48","orange","gold", "salmon"),name="")


total_p=melt(total_q[ ,c("Qtype","clade","genus")],id="Qtype",value.name="taxon")
gq2=ggplot(total_p, aes(x=taxon, y=..count../sum(..count..),fill=Qtype))+geom_bar(position="fill")+facet_grid(~variable,scales = "free", space = "free")+geom_text(aes(label=..count..),stat="count",position=position_fill(vjust=0.5))+theme(legend.position=c("top") ,legend.direction="horizontal")+ylab("Proportion")+xlab("")+scale_fill_manual(values=c("gray48","orange","gold", "salmon"),name="")

quartz(width=15.4,height=10.6)
grid.arrange(g1,gq1,g2,gq2,nrow=2,ncol=2)

#################################################################### Dfoil ##############################################################



names_v=c("clade","P1","P2","P3","P4","Out","chrom1","position", "AAAAA" , "AAABA" , "AABAA" , "AABBA" , "ABAAA" , "ABABA" , "ABBAA" , "ABBBA" , "BAAAA" , "BAABA" , "BABAA" , "BABBA" , "BBAAA" ,"BBABA","BBBAA","BBBBA",'chromdup','coord','total','dtotal','T12','T34','T1234','DFO_left','DFO_right','DFO_total','DFO_stat','DFO_chisq','DFO_Pvalue','DIL_left','DIL_right','DIL_total','DIL_stat','DIL_chisq','DIL_Pvalue','DFI_left','DFI_right','DFI_total','DFI_stat','DFI_chisq','DFI_Pvalue','DOL_left','DOL_right','DOL_total','DOL_stat','DOL_chisq','DOL_Pvalue','introgression','introgna','intrognone','introg13','introg14','introg23','introg24','introg31','introg41','introg32','introg42','introg123','introg124')

total_d=read.csv("/Users/Anton/Downloads/droso_dfoil_results.txt",stringsAsFactors=FALSE)
names(total_d)=names_v
total_d=total_d[total_d$T12<total_d$T34,]

total_d$Genus="Drosophila"
total_d$introgressionid=ifelse(total_d$introgression=="none","None",
                             ifelse(total_d$introgression=="123" | total_d$introgression=="124","Ancestral","Inter-group"))
         
total_d=cbind(total_d,get_intropair_dfoil(total_d))        



#################################################################### Branch Length Test ################################################

names_vb=c("clade","P1out","P2out","P3out","CountP1","CountP2","CountP3","PvalueChi","meanT_concord","meanT_discord1","meanT_discord2","PvalueWCOMC1","PvalueWCOMC2","PvalueWC1C2")
total_b=read.csv("/Users/Anton/Downloads/droso_blt_results.txt",stringsAsFactors=FALSE)
names(total_b)=names_vb
total_b$PvalueChi=p.adjust(total_b$PvalueChi,method="fdr")
total_b$PvalueWCOMC1=p.adjust(total_b$PvalueWCOMC1,method="fdr") 
total_b$PvalueWCOMC2=p.adjust(total_b$PvalueWCOMC2,method="fdr") 
total_b$PvalueWC1C2=p.adjust(total_b$PvalueWC1C2,method="fdr") 



total_b=cbind(total_b,get_intropair_blt(total_b))  



names(counts),counts,chi,mean(ccom),mean(c1),mean(c2),w_testc1,w_testc2,w_test


#####################################Plotting Compare######################
d1=get_intomatrix_dfoil(C1,total_d)
d2=get_intomatrix_dfoil(C2,total_d)
d3=get_intomatrix_dfoil(C3,total_d)
d4=get_intomatrix_dfoil(C4,total_d)
d5=get_intomatrix_dfoil(C5,total_d)
d6=get_intomatrix_dfoil(C6,total_d)
d7=get_intomatrix_dfoil(C7,total_d)
d8=get_intomatrix_dfoil(C8,total_d)
d9=get_intomatrix_dfoil(C9,total_d)

q1=get_intomatrix_quibl(C1,total_q)
q2=get_intomatrix_quibl(C2,total_q)
q3=get_intomatrix_quibl(C3,total_q)
q4=get_intomatrix_quibl(C4,total_q)
q5=get_intomatrix_quibl(C5,total_q)
q6=get_intomatrix_quibl(C6,total_q)
q7=get_intomatrix_quibl(C7,total_q)
q8=get_intomatrix_quibl(C8,total_q)
q9=get_intomatrix_quibl(C9,total_q)

qc1=get_intomatrix_quibl_s(C1,total_q)
qc2=get_intomatrix_quibl_s(C2,total_q)
qc3=get_intomatrix_quibl_s(C3,total_q)
qc4=get_intomatrix_quibl_s(C4,total_q)
qc5=get_intomatrix_quibl_s(C5,total_q)
qc6=get_intomatrix_quibl_s(C6,total_q)
qc7=get_intomatrix_quibl_s(C7,total_q)
qc8=get_intomatrix_quibl_s(C8,total_q)
qc9=get_intomatrix_quibl_s(C9,total_q)

b1=get_intomatrix_blt(C1,total_b)
b2=get_intomatrix_blt(C2,total_b)
b3=get_intomatrix_blt(C3,total_b)
b4=get_intomatrix_blt(C4,total_b)
b5=get_intomatrix_blt(C5,total_b)
b6=get_intomatrix_blt(C6,total_b)
b7=get_intomatrix_blt(C7,total_b)
b8=get_intomatrix_blt(C8,total_b)
b9=get_intomatrix_blt(C9,total_b)



grid.arrange(b1,b2,b3,b4,ncol=1,nrow=4)
quartz.save("C1_C4_blt.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)
grid.arrange(b5,b6,b7,b8,b9,ncol=1,nrow=5)
quartz.save("C5_C9_blt.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)



grid.arrange(d1,q1,qc1,b1,d2,q2,qc2,b2,d3,q3,qc3,b3,d4,q4,qc4,b4,ncol=4,nrow=4)
grid.arrange(d5,q5,qc5,b5,d6,q6,qc6,b6,d7,q7,qc7,b7,d8,q8,qc8,b8,d9,q9,qc9,b9,ncol=4,nrow=5)

quartz.save("C1_C4_dfoil_quibl.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)
quartz.save("C5_C9_dfoil_quibl.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)