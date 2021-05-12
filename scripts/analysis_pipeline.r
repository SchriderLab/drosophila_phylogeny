library("ape")
library("ggplot2")
library("pals")
library("reshape2")
library("dplyr")
library("gridExtra")
library("svMisc")
library("MCMCtreeR")
library("phytools")
options(repr.matrix.max.cols=100, repr.matrix.max.rows=100)

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


#Identify introgressing taxa pair from DCT/BLT result 

get_intropair_dctblt=function(m)
{
    pair_v=rbind()
    for (i in 1:nrow(m))
    {
       pair_dct=as.character(m[i,c("P1out","P2out","P3out")][which(rank(m[i,c("CountP1","CountP2","CountP3")],ties.method = "random")!=2)])
       if (m[i,"meanT_discord2"] < m[i,"meanT_discord1"])       
       {
           pair_blt=pair_dct
           
       }else{
           
           pair_blt=as.character(m[i,c("P1out","P2out","P3out")][which(rank(m[i,c("CountP1","CountP2","CountP3")],ties.method = "random")!=1)])
       } 
       pair_v=rbind(pair_v,c(sort(pair_dct),m[i,"PvalueChi"]<0.05,sort(pair_blt),m[i,"PvalueWC1C2"]<0.05,m[i,"PvalueChi"]<0.05 & m[i,"PvalueWC1C2"]<0.05,all(pair_dct %in% pair_blt))) 
    }
    pair_v=data.frame(pair_v)
    names(pair_v)=c("i1_dct","i2_dct","pass_dct","i1_blt","i2_blt","pass_blt","pass_dctblt","i1i2_overlap")
    return(pair_v)
}





#Identify introgressing taxa pair from HyDe result 
get_intropair_hyde=function(m)
{
    pair_v=rbind()
    for (i in 1:nrow(m))
    {
        if(m[i,"Gamma"]>0.5)
        {
            pair=as.character(m[i,c("Hybrid","P2")]) 
        }else{
            pair=as.character(m[i,c("Hybrid","P1")]) 
        }    
        
        if(m[i,"Pvalue"]<0.05)
        {
            pair_v=rbind(pair_v,c(pair,"TRUE"))
        }else{
            pair_v=rbind(pair_v,c(pair,"FALSE"))
        }    
    }
    pair_v=data.frame(pair_v)
    names(pair_v)=c("i1","i2","pass")
    rownames(pair_v) = c()
    pair_v[,c("i1","i2")]=t(apply(pair_v[,c("i1","i2")],1,sort))
    return(pair_v)
}




# Plot matrix from BLT/Chisq/QuiBl results 
get_intromatrix_blt=function(taxa,dat,sig)
{
    m=matrix(0, nrow=length(taxa),ncol=length(taxa)) 
    m=data.frame(m)
    colnames(m)=taxa
    rownames(m)=taxa
    n=m
    dat_int=dat[dat$i1 %in% taxa | dat$i2 %in% taxa ,c("i1","i2","pass")]
    for (i in 1:nrow(dat_int))
    {
        a=dat_int[i,"i1"]
        b=dat_int[i,"i2"]
        if (dat_int[i,"pass"]=="TRUE")
        {    
            m[rownames(m)==a,colnames(m)==b]=m[rownames(m)==a,colnames(m)==b]+1
            n[rownames(n)==a,colnames(n)==b]=m[rownames(n)==a,colnames(n)==b]+1
        }else{
            n[rownames(n)==a,colnames(n)==b]=m[rownames(n)==a,colnames(n)==b]+1  
        }   
        
    }    
    m=m+t(m)
    n=n+t(n)
    m[upper.tri(m)]=NA
    n[upper.tri(n)]=NA
    #Normalazing by the triplets that contain an introgressing pair
    norm_m=as.matrix(m/n)
    norm_m[is.nan(norm_m)]=0
    melted_m=melt(as.matrix(norm_m),na.rm = TRUE)
    g1=ggplot(melted_m, aes(Var2, Var1, fill = value))+geom_tile(color = "white")+scale_fill_gradientn(colors=viridis(100),space = "Lab",name="") +theme_minimal()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+coord_fixed()+labs(x="",y="")
    return(g1)
}    

#################################################################### Define Clades #######################################################

tt=read.tree("MLrooted.tre")

C7=extract.clade(tt,168)$tip.label
C6=extract.clade(tt,185)$tip.label
C9=extract.clade(tt,204)$tip.label
C8=extract.clade(tt,226)$tip.label
C2=extract.clade(tt,245)$tip.label
C3=extract.clade(tt,256)$tip.label
C5=extract.clade(tt,265)$tip.label
C4=extract.clade(tt,289)$tip.label
C1=extract.clade(tt,307)$tip.label

sp_space=list(C1,C2,C3,C4,C5,C6,C7,C8,C9)
#################################################################### QuiBL ##############################################################

total_q=read.csv("droso_quibl_results.txt",header=F,stringsAsFactors=F)
names(total_q)=c("clade","id","P1","P2","P3","outgroup","Com1","Com2","mixprop1","mixprop2","lambda2Dist","lambda1Dist","BIC2Dist","BIC1Dist","count")
total_q$clade=unlist(lapply(strsplit(as.character(total_q$clade), "_"),"[",1))
#total_q=total_q[complete.cases(total_q), ] 
total_q$BICdiff = total_q$BIC2-total_q$BIC1
total_q$triplet=as.character(apply(total_q[,c("P1","P2","P3")],1,paste,collapse="_"))
total_q=cbind(total_q,get_intropair(total_q))
total_q$BICdiff = total_q$BIC2-total_q$BIC1
total_q$pass=total_q$BICdiff < -30 

total_q$clade=ifelse(total_q$clade=="C1","C7",
              ifelse(total_q$clade=="C2","C6",
              ifelse(total_q$clade=="C3","C9",
              ifelse(total_q$clade=="C4","C8",
              ifelse(total_q$clade=="C5","C2",
              ifelse(total_q$clade=="C6","C3",
              ifelse(total_q$clade=="C7","C5",
              ifelse(total_q$clade=="C8","C4","C1"))))))))


#identify triplets with the most counts
common=c()
i=0
for (trip in unique(as.character(total_q$triplet)))
{
    i=i+1
    progress(i,length(unique(as.character(total_q$triplet))))
    gcount=total_q[total_q$triplet==trip,"count"]+1
    pos=rank(gcount,ties.method ="random")
    cl=c("discord2","discord1","common")
    common=c(common,cl[pos])
}

total_q$common=common
total_q=total_q[complete.cases(total_q$pass),]
total_q=total_q[total_q$common!="common",]
#total_q=total_q[total_q$common!="discord1",]

  

#################################################################### Branch Length Test ################################################
names_vb=c("clade","P1out","P2out","P3out","CountP1","CountP2","CountP3","PvalueChi","meanT_concord","meanT_discord1","meanT_discord2","PvalueWCOMC1","PvalueWCOMC2","PvalueWC1C2")
total_b=read.csv("droso_blt_results.txt",stringsAsFactors=FALSE,header=F)
#total_b=read.csv("droso_blt_results_treeshrink.txt",stringsAsFactors=FALSE,header=F)
names(total_b)=names_vb
total_b$PvalueChi=p.adjust(total_b$PvalueChi,method="fdr")
total_b$PvalueWCOMC1=p.adjust(total_b$PvalueWCOMC1,method="fdr") 
total_b$PvalueWCOMC2=p.adjust(total_b$PvalueWCOMC2,method="fdr") 
total_b$PvalueWC1C2=p.adjust(total_b$PvalueWC1C2,method="fdr") 


bltdct=get_intropair_dctblt(total_b)
total_b=cbind(total_b,bltdct)
#write.csv(total_b,"droso_blt_results_analyzed.txt",quote = F,row.names = F)
bltdct_overlap=total_b[total_b$i1i2_overlap == TRUE,c("i1_dct","i2_dct","pass_dctblt")]
names(bltdct_overlap)=c("i1","i2","pass")

blt_alone=total_b[,c("i1_blt","i2_blt","pass_blt")]
names(blt_alone)=c("i1","i2","pass")

dct_alone=total_b[,c("i1_dct","i2_dct","pass_dct")]
names(dct_alone)=c("i1","i2","pass")



#################################################################### Hyde ################################################
total_h=read.table("all_hyde.txt",stringsAsFactors=FALSE,header=T)
total_h=total_h[total_h$Gamma <= 1 & total_h$Gamma >= 0,]
total_h$Pvalue=p.adjust(total_h$Pvalue, method = "bonferroni")
h=get_intropair_hyde(total_h)
names(h)=c("i1","i2","pass")

total_h=cbind(total_h,h)



################################################################### Overlap Hypergeometric test ######################################## 
pair_name=apply(total_b[,c("i1_dct","i2_dct")],1,paste,collapse="_")
total_b$pair_name=pair_name
total_b_i1i2=total_b



#Stringent unique introgression (main figure)
for (cl in c("C1","C2","C3","C4","C5","C6","C7","C8","C9"))
{    
    overl=0
    sig_dct=0
    sig_blt=0
    no_sig=0
    a=total_b_i1i2[total_b_i1i2$clade==cl,]
    for (p in unique(a$pair_name))
    {
        if (any(a$pair_name==p & a$pass_dctblt=="TRUE" & a$i1i2_overlap=="TRUE"))
        {
            overl=overl+1
        } else if (all(c(any(a$pair_name==p & a$pass_dct!="FALSE"),any(a$pair_name==p & a$pass_blt!="FALSE")))) {
            sig_dct=sig_dct+1
            sig_blt=sig_blt+1
        } else if (any(a$pair_name==p & a$pass_dct!="FALSE")) {
            sig_dct=sig_dct+1
        } else if (any(a$pair_name==p & a$pass_blt!="FALSE")) {
            sig_blt=sig_blt+1
        } else {
            no_sig=no_sig+1  
        }    
    }
    tot=sig_dct+sig_blt+overl+no_sig
    tot_blt=sig_blt+overl
    tot_dct=sig_dct+overl
    print(cl)
    hp=phyper(overl, tot_blt, tot - tot_blt, tot_dct, lower.tail = FALSE)
    print(c(tot_dct-overl,overl,tot_blt-overl,no_sig))
    print(tot)
    quartz(width=4, height=3.2)
    plot(c(0.8,1.2),c(1,1),cex=15,xlim=c(0,2),main=paste("P = ",hp),xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
    text(1,1,overl)
    text(0.55,1,tot_dct-overl)
    text(1.45,1,tot_blt-overl)
    text(1.9,0.65,no_sig)
    quartz.save(paste(cl,"_vennstringentoverlap_uniq.pdf",sep=""), type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
    dev.off()
}


#Stringent unique introgression power test
hyper_subsample_uniq=function(tabl,size=8,cl)
{
    sub_tabl=tabl[tabl$clade==cl,]
    sps=unique(c(sub_tabl$P1out,sub_tabl$P2out,sub_tabl$P3out))
    p_vals=c()
    for(i in 1:10000)
    {
        overl=0
        sig_dct=0
        sig_blt=0
        no_sig=0
        sps_sub=sample(sps,size=size)
        a=sub_tabl[sub_tabl$P1out %in% sps_sub & sub_tabl$P2out %in% sps_sub & sub_tabl$P3out %in% sps_sub,]
        for (p in unique(a$pair_name))
        {
            if (any(a$pair_name==p & a$pass_dctblt=="TRUE" & a$i1i2_overlap=="TRUE"))
            {
                overl=overl+1 
            } else if (all(c(any(a$pair_name==p & a$pass_dct!="FALSE"),any(a$pair_name==p & a$pass_blt!="FALSE")))) {
                sig_dct=sig_dct+1
                sig_blt=sig_blt+1
            } else if (any(a$pair_name==p & a$pass_dct!="FALSE")) {
                sig_dct=sig_dct+1
            } else if (any(a$pair_name==p & a$pass_blt!="FALSE")) {
                sig_blt=sig_blt+1
            } else {
                no_sig=no_sig+1  
            }    
        }
        tot=sig_dct+sig_blt+overl+no_sig
        tot_blt=sig_blt+overl
        tot_dct=sig_dct+overl
        hp=phyper(overl, tot_blt, tot - tot_blt, tot_dct, lower.tail = FALSE)
        p_vals=c(p_vals,hp)
        
    }    
    return(p_vals)
    
}    

quartz(width=10, height=10)
par(mfrow=c(3,3))
i=1
for (cl in c("C1","C2","C3","C4","C5","C6","C7","C8","C9"))
{ 
    pval=hyper_subsample_uniq(total_b_i1i2,8,cl)
    n=round(sum(pval< 0.05 & pval!=0)/sum(pval!=0),digits=3)
   
    hist(log(pval),col="grey",xlab="",ylab="Frequency",prob=TRUE,main=paste("Clade ",i," (",n,")"),nclass=20)
    abline(v=log(0.05),lty=3,col="red",lwd=2)
      
    
    i=i+1
}

quartz.save("Power_uniqstringent.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
dev.off()


#Relaxed unique introgression 
for (cl in c("C1","C2","C3","C4","C5","C6","C7","C8","C9"))
{    
    overl=0
    sig_dct=0
    sig_blt=0
    no_sig=0
    a=total_b_i1i2[total_b_i1i2$clade==cl,]
    for (p in unique(a$pair_name))
    {
        if (any(a$pair_name==p & a$pass_dctblt=="TRUE" & a$i1i2_overlap=="TRUE"))
        {
            overl=overl+1 
        } else if (all(c(any(a$pair_name==p & a$pass_dct!="FALSE"),any(a$pair_name==p & a$pass_blt!="FALSE")))) {
            overl=overl+1
        } else if (any(a$pair_name==p & a$pass_dct!="FALSE")) {
            sig_dct=sig_dct+1
        } else if (any(a$pair_name==p & a$pass_blt!="FALSE")) {
            sig_blt=sig_blt+1
        } else {
            no_sig=no_sig+1  
        }    
    }
    tot=sig_dct+sig_blt+overl+no_sig
    tot_blt=sig_blt+overl
    tot_dct=sig_dct+overl
    print(cl)
    hp=phyper(overl, tot_blt, tot - tot_blt, tot_dct, lower.tail = FALSE)
    print(c(tot_dct-overl,overl,tot_blt-overl,no_sig))
    print(tot)
    quartz(width=4, height=3.2)
    plot(c(0.8,1.2),c(1,1),cex=15,xlim=c(0,2),main=paste("P = ",hp),xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
    text(1,1,overl)
    text(0.55,1,tot_dct-overl)
    text(1.45,1,tot_blt-overl)
    text(1.9,0.65,no_sig)
    quartz.save(paste(cl,"_vennrelaxedoverlap_uniq.pdf",sep=""), type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
    dev.off()
}


#Relaxed unique introgression power test
hyper_subsample_uniqrelaxed=function(tabl,size=8,cl)
{
    sub_tabl=tabl[tabl$clade==cl,]
    sps=unique(c(sub_tabl$P1out,sub_tabl$P2out,sub_tabl$P3out))
    p_vals=c()
    for(i in 1:10000)
    {
        overl=0
        sig_dct=0
        sig_blt=0
        no_sig=0
        sps_sub=sample(sps,size=size)
        a=sub_tabl[sub_tabl$P1out %in% sps_sub & sub_tabl$P2out %in% sps_sub & sub_tabl$P3out %in% sps_sub,]
        for (p in unique(a$pair_name))
        {
            if (any(a$pair_name==p & a$pass_dctblt=="TRUE" & a$i1i2_overlap=="TRUE"))
            {
                overl=overl+1 
            } else if (all(c(any(a$pair_name==p & a$pass_dct!="FALSE"),any(a$pair_name==p & a$pass_blt!="FALSE")))) {
                overl=overl+1
            } else if (any(a$pair_name==p & a$pass_dct!="FALSE")) {
                sig_dct=sig_dct+1
            } else if (any(a$pair_name==p & a$pass_blt!="FALSE")) {
                sig_blt=sig_blt+1
            } else {
                no_sig=no_sig+1  
            }    
        }
        tot=sig_dct+sig_blt+2*overl+no_sig
        tot_blt=sig_blt+overl
        tot_dct=sig_dct+overl
        hp=phyper(overl, tot_blt, tot - tot_blt, tot_dct, lower.tail = FALSE)
        p_vals=c(p_vals,hp)
        
    }    
    return(p_vals)
    
}    

quartz(width=10, height=10)
par(mfrow=c(3,3))
i=1
for (cl in c("C1","C2","C3","C4","C5","C6","C7","C8","C9"))
{ 
    pval=hyper_subsample_uniqrelaxed(total_b,8,cl)
    n=round(sum(pval< 0.05 & pval!=0)/sum(pval!=0),digits=3)
    
    hist(log(pval),col="grey",xlab="",ylab="Frequency",prob=TRUE,main=paste("Clade ",i," (",n,")"),nclass=20)
    abline(v=log(0.05),lty=3,col="red",lwd=2)
    i=i+1
}

quartz.save("Power_uniqrelaxed.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
dev.off()


#Introgression all triplets
for (cl in c("C1","C2","C3","C4","C5","C6","C7","C8","C9"))
{    
    a=total_b_i1i2[total_b_i1i2$clade==cl,]
    tot=nrow(a)
    overl=sum(as.numeric(a$pass_dctblt=="TRUE" & a$i1i2_overlap=="TRUE"))
    tot_dct=sum(as.numeric(a$pass_dct=="TRUE"))
    tot_blt=sum(as.numeric(a$pass_blt=="TRUE"))
    print(cl)
    hp=phyper(overl, tot_blt, tot - tot_blt, tot_dct, lower.tail = FALSE)
    print(c(tot_dct-overl,overl,tot_blt-overl,sum(as.numeric(a$pass_dct=="FALSE" & a$pass_blt=="FALSE"))))
    print(tot)
    quartz(width=4, height=3.2)
    plot(c(0.8,1.2),c(1,1),cex=15,xlim=c(0,2),main=paste("P = ",hp),xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
    text(1,1,overl)
    text(0.55,1,tot_dct-overl)
    text(1.45,1,tot_blt-overl)
    text(1.9,0.65,sum(as.numeric(a$pass_dct=="FALSE" & a$pass_blt=="FALSE")))
    quartz.save(paste(cl,"_venn_alltrips.pdf",sep=""), type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
    dev.off()
}

hyper_subsample=function(tabl,size=8,cl)
{
    sub_tabl=tabl[tabl$clade==cl,]
    sps=unique(c(sub_tabl$P1out,sub_tabl$P2out,sub_tabl$P3out))
    p_vals=c()
    for(i in 1:10000)
    {
        
        sps_sub=sample(sps,size=size)
        a=sub_tabl[sub_tabl$P1out %in% sps_sub & sub_tabl$P2out %in% sps_sub & sub_tabl$P3out %in% sps_sub,]
        tot=nrow(a)
        overl=sum(as.numeric(a$pass_dctblt=="TRUE" & a$i1i2_overlap=="TRUE"))
        tot_dct=sum(as.numeric(a$pass_dct=="TRUE"))
        tot_blt=sum(as.numeric(a$pass_blt=="TRUE"))
        hyper_p=phyper(overl, tot_blt, tot - tot_blt, tot_dct, lower.tail = FALSE)
        p_vals=c(p_vals,hyper_p)
        
    }    
    return(p_vals)
    
}    

#Introgression power test all triplets
quartz(width=10, height=10)
par(mfrow=c(3,3))
i=1
for (cl in c("C1","C2","C3","C4","C5","C6","C7","C8","C9"))
{ 
    pval=hyper_subsample(total_b,8,cl)
    n=round(sum(pval< 0.05 & pval!=0)/sum(pval!=0),digits=3)
    
    hist(log(pval),col="grey",xlab="",ylab="Frequency",prob=TRUE,main=paste("Clade ",i," (",n,")"),nclass=20)
    abline(v=log(0.05),lty=3,col="red",lwd=2)
    i=i+1
}
quartz.save("Power_alltriplets.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
dev.off()


#####################################Printing pairwise introgression matrices#####################################
print_save_matrix=function(clades,data,name)
{
    ids=1
    for (cl in clades)
    {
       pair_m=get_intomatrix_blt(cl,data)
       print(pair_m)
       quartz.save(paste("C",ids,name,".pdf",sep=""), type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
       ids=ids+1  
    }    
}    

#Agreement between BLT and DCT
print_save_matrix(sp_space,bltdct_overlap,"_blt_chi")
#DCT
print_save_matrix(sp_space,dct_alone,"_chi")
#BLT
print_save_matrix(sp_space,blt_alone,"_blt")
#QuibL
print_save_matrix(sp_space,total_q,"_quibl")
#HyDe
print_save_matrix(sp_space,total_h,"_hyde")



#####################################Time-introgression plot######################################################
phy_mcmc=readMCMCtree("schemeA.tre")
phy=phy_mcmc$apePhy
tree_h=nodeheight(phy, node=159)
age_v=c()
node_v=c()
for (i in 1:nrow(total_b))
{
    age_pair=tree_h-findMRCA(phy,c(total_b[i,"i1_chi"],total_b[i,"i2_chi"]),type="height")
    node_pair=findMRCA(phy,c(total_b[i,"i1_chi"],total_b[i,"i2_chi"]),type="node")
    age_v=c(age_v,age_pair)
    node_v=c(node_v,node_pair)
}    
total_b$tmrca=age_v
total_b$mrca=node_v


#Time-introgression plot for the entire tree 
all_a=total_b[,"tmrca"]
sig_a=total_b[total_b$pass_chi==TRUE & total_b$pass_wilx==TRUE ,"tmrca"]
d_a=data.frame(Age=c(all_a,sig_a),Distribution=c(rep("All",length(all_a)),rep("Sig.",length(sig_a))))
quartz(width=6.5, height=4.3)
ggplot(d_a, aes(x=Age, fill=Distribution))+geom_density(alpha=1,position = "stack")+scale_fill_manual(values=c("dodgerblue4", "gold"))+ geom_vline(xintercept = unique(sig_a),linetype="dashed", color = "red", size=1))



for (cl in c("C1","C2","C3","C4","C5","C6","C7","C8","C9"))
{    
    all_a=total_b[total_b$clade==cl,"tmrca"]
    sig_a=total_b[total_b$clade==cl & total_b$pass_chi==TRUE & total_b$pass_wilx==TRUE ,"tmrca"]
    d_a=data.frame(Age=c(all_a,sig_a),Distribution=c(rep("All",length(all_a)),rep("Sig.",length(sig_a))))
    quartz(width=6.5, height=4.3)
    p=ggplot(d_a, aes(x=Age, fill=Distribution))+geom_density(alpha=1,position = "stack")+scale_fill_manual(values=c("dodgerblue4", "gold"))
    print(p)
    quartz.save(paste(cl,"_agedistr.png",sep=""), type = "png",antialias=F,bg="white",dpi=400,pointsize=12)
    dev.off()
}

#####################################Number of introgression events############################################
for (cl in c("C1","C2","C3","C4","C5","C6","C7","C8","C9"))
{   
    print(cl)
    sig_mrca=total_b[total_b$clade==cl & total_b$pass_chi==TRUE & total_b$pass_wilx==TRUE ,"mrca"]
    print(table(sig_mrca))
    tot_mrca=total_b[total_b$clade==cl & total_b$mrca %in% sig_mrca ,"mrca"]
    print(table(tot_mrca))
}    


cl=1
quartz(width=6.5, height=5.5)
for (cl_node in c(307,245,256,289,265,185,168,226,204))
{   
    plot(ladderize(extract.clade(phy,cl_node)),cex=0.8)
    quartz.save(paste("C",cl,"_clade.pdf",sep=""), type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
    #dev.off()
    cl=cl+1
}    

#####################################HyDe loci introgression############################################
h=read.table("hyde_loci_wintrogression.txt")
names(h)=c("id","P1","P2","P3","Gamma","len")
h$clade=ifelse(apply(apply(h[,c("P1","P2","P3")],2,"%in%",C1),1,all),"C1",
        ifelse(apply(apply(h[,c("P1","P2","P3")],2,"%in%",C2),1,all),"C2",
        ifelse(apply(apply(h[,c("P1","P2","P3")],2,"%in%",C3),1,all),"C3",
        ifelse(apply(apply(h[,c("P1","P2","P3")],2,"%in%",C4),1,all),"C4",
        ifelse(apply(apply(h[,c("P1","P2","P3")],2,"%in%",C5),1,all),"C5",
        ifelse(apply(apply(h[,c("P1","P2","P3")],2,"%in%",C6),1,all),"C6",
        ifelse(apply(apply(h[,c("P1","P2","P3")],2,"%in%",C7),1,all),"C7",
        ifelse(apply(apply(h[,c("P1","P2","P3")],2,"%in%",C8),1,all),"C8",       
        ifelse(apply(apply(h[,c("P1","P2","P3")],2,"%in%",C9),1,all),"C9","Noclade")))))))))
h=h[h$clade!="Noclade",]

m_occ=c()
for (l in unique(h$id))
{
    v=c(0,0,0,0,0,0,0,0,0)
    names(v)=c("C1","C2","C3","C4","C5","C6","C7","C8","C9")
    occ=names(table(h[h$id==l,"clade"]))
    v[occ]=1
    m_occ=rbind(m_occ,c(id=l,v))
    
}
m_occ=data.frame(m_occ)    
heatmap.2(apply(m_occ[,2:10],2,as.numeric))

zz=melt(m_occ,id=c("id"))    

