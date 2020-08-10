library("ape")
library('ggplot2')
library('pals')
library('reshape2')
library('dplyr')
library('gridExtra')
library("svMisc")
library("MCMCtreeR")
library('phytools')

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
get_intropair_blt_wilcox=function(m)
{
    pair_v=rbind()
    for (i in 1:nrow(m))
    {
        progress(i,nrow(m))
        pair=as.character(m[i,c("P1out","P2out","P3out")][which(rank(m[i,c("CountP1","CountP2","CountP3")],ties.method = "random")!=2)])
        if(m[i,"PvalueWC1C2"]<0.05)
        {
            
            if(m[i,c("CountP1","CountP2","CountP3")][which(rank(m[i,c("CountP1","CountP2","CountP3")])!=3)][1] > m[i,c("CountP1","CountP2","CountP3")][which(rank(m[i,c("CountP1","CountP2","CountP3")])!=3)][2] & m[i,"meanT_discord1"] < m[i,"meanT_discord2"] | m[i,c("CountP1","CountP2","CountP3")][which(rank(m[i,c("CountP1","CountP2","CountP3")])!=3)][1] < m[i,c("CountP1","CountP2","CountP3")][which(rank(m[i,c("CountP1","CountP2","CountP3")])!=3)][2] & m[i,"meanT_discord1"] > m[i,"meanT_discord2"])
            {
               pair=c(pair,"TRUE","TRUE") 
            }else{
               pair=c(pair,"TRUE","FALSE")  
            }
            pair_v=rbind(pair_v,pair)
        }else{
            pair=c(pair,"FALSE","FALSE")
            pair_v=rbind(pair_v,pair)
        }    
    }
    pair_v=data.frame(pair_v)
    names(pair_v)=c("i1","i2","pass","totalbl_pass")
    rownames(pair_v) = c()
    pair_v[,c("i1","i2")]=t(apply(pair_v[,c("i1","i2")],1,sort))
    return(pair_v)
}


get_intropair_blt_chisq=function(m)
{
    pair_v=rbind()
    for (i in 1:nrow(m))
    {
        progress(i,nrow(m))
        pair=as.character(m[i,c("P1out","P2out","P3out")][which(rank(m[i,c("CountP1","CountP2","CountP3")],ties.method = "random")!=2)])
        if(m[i,"PvalueChi"]<0.05)
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
get_intomatrix_blt=function(taxa,dat,sig)
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




#################################################################### QuiBL ##############################################################

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





#################################################################### Dfoil ##############################################################



names_v=c("clade","P1","P2","P3","P4","Out","chrom1","position", "AAAAA" , "AAABA" , "AABAA" , "AABBA" , "ABAAA" , "ABABA" , "ABBAA" , "ABBBA" , "BAAAA" , "BAABA" , "BABAA" , "BABBA" , "BBAAA" ,"BBABA","BBBAA","BBBBA",'chromdup','coord','total','dtotal','T12','T34','T1234','DFO_left','DFO_right','DFO_total','DFO_stat','DFO_chisq','DFO_Pvalue','DIL_left','DIL_right','DIL_total','DIL_stat','DIL_chisq','DIL_Pvalue','DFI_left','DFI_right','DFI_total','DFI_stat','DFI_chisq','DFI_Pvalue','DOL_left','DOL_right','DOL_total','DOL_stat','DOL_chisq','DOL_Pvalue','introgression','introgna','intrognone','introg13','introg14','introg23','introg24','introg31','introg41','introg32','introg42','introg123','introg124')

total_d=read.csv("droso_dfoil_results.txt",stringsAsFactors=FALSE,header=F)
names(total_d)=names_v
total_d=total_d[total_d$T12<total_d$T34,]

total_d$Genus="Drosophila"
total_d$introgressionid=ifelse(total_d$introgression=="none","None",
                             ifelse(total_d$introgression=="123" | total_d$introgression=="124","Ancestral","Inter-group"))
         
total_d=cbind(total_d,get_intropair_dfoil(total_d))        



#################################################################### Branch Length Test ################################################

names_vb=c("clade","P1out","P2out","P3out","CountP1","CountP2","CountP3","PvalueChi","meanT_concord","meanT_discord1","meanT_discord2","PvalueWCOMC1","PvalueWCOMC2","PvalueWC1C2")
total_b=read.csv("droso_blt_results.txt",stringsAsFactors=FALSE,header=F)
names(total_b)=names_vb
total_b$PvalueChi=p.adjust(total_b$PvalueChi,method="fdr")
total_b$PvalueWCOMC1=p.adjust(total_b$PvalueWCOMC1,method="fdr") 
total_b$PvalueWCOMC2=p.adjust(total_b$PvalueWCOMC2,method="fdr") 
total_b$PvalueWC1C2=p.adjust(total_b$PvalueWC1C2,method="fdr") 

b_wilcox=get_intropair_blt_wilcox(total_b)
b_chisq=get_intropair_blt_chisq(total_b)

m_ch=b_chisq
m_wilx=b_wilcox
m_overlap=cbind(b_chisq[,1:2],b_chisq$pass=="TRUE" & b_wilcox$pass=="TRUE")
names(b_wilcox)=c("i1_wilx","i2_wilx","pass_wilx","totalbl_pass")
names(b_chisq)=c("i1_chi","i2_chi","pass_chi")
names(m_overlap)=c("i1","i2","pass")



total_b=cbind(total_b,b_chisq,b_wilcox)
total_b$clade=ifelse(total_b$clade=="C1","C7",
              ifelse(total_b$clade=="C2","C6",
              ifelse(total_b$clade=="C3","C9",
              ifelse(total_b$clade=="C4","C8",
              ifelse(total_b$clade=="C5","C2",
              ifelse(total_b$clade=="C6","C3",
              ifelse(total_b$clade=="C7","C5",
              ifelse(total_b$clade=="C8","C4","C1"))))))))

######################################Overlap Hypergeometric test##################################### 

pair_name=apply(total_b[,c("i1_chi","i2_chi")],1,paste,collapse="_")
total_b$pair_name=pair_name
for (cl in c("C1","C2","C3","C4","C5","C6","C7","C8","C9"))
{    
    overl=0
    sig_chi=0
    sig_wilx=0
    no_sig=0
    a=total_b[total_b$clade==cl,]
    for (p in unique(a$pair_name))
    {
        if (any(a$pair_name==p & a$pass_chi!="FALSE" & a$pass_wilx!="FALSE"))
        {
            overl=overl+1 
        } else if (all(c(any(a$pair_name==p & a$pass_chi!="FALSE"),any(a$pair_name==p & a$pass_wilx!="FALSE")))) {
            sig_chi=sig_chi+1
            sig_wilx=sig_wilx+1
        } else if (any(a$pair_name==p & a$pass_chi!="FALSE")) {
            sig_chi=sig_chi+1
        } else if (any(a$pair_name==p & a$pass_wilx!="FALSE")) {
            sig_wilx=sig_wilx+1
        } else {
            no_sig=no_sig+1  
        }    
    }
    tot=sig_chi+sig_wilx+2*overl+no_sig
    tot_wilx=sig_wilx+overl
    tot_chi=sig_chi+overl
    print(cl)
    hp=phyper(overl, tot_wilx, tot - tot_wilx, tot_chi, lower.tail = FALSE)
    print(c(tot_chi-overl,overl,tot_wilx-overl,no_sig))
    print(tot)
    quartz(width=4, height=3.2)
    plot(c(0.8,1.2),c(1,1),cex=15,xlim=c(0,2),main=paste("P = ",hp),xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
    text(1,1,overl)
    text(0.55,1,tot_chi-overl)
    text(1.45,1,tot_wilx-overl)
    text(1.9,0.65,no_sig)
    quartz.save(paste(cl,"_venn.pdf",sep=""), type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
    dev.off()
}


for (cl in c("C1","C2","C3","C4","C5","C6","C7","C8","C9"))
{    
    a=total_b[total_b$clade==cl,]
    tot=nrow(a)
    overl=sum(as.numeric(a$pass_chi!="FALSE" & a$pass_wilx!="FALSE"))
    tot_chi=sum(as.numeric(a$pass_chi=="TRUE"))
    tot_wilx=sum(as.numeric(a$pass_wilx=="TRUE"))
    print(cl)
    hp=phyper(overl, tot_wilx, tot - tot_wilx, tot_chi, lower.tail = FALSE)
    print(c(tot_chi-overl,overl,tot_wilx-overl,sum(as.numeric(a$pass_chi=="FALSE" & a$pass_wilx=="FALSE"))))
    print(tot)
    quartz(width=4, height=3.2)
    plot(c(0.8,1.2),c(1,1),cex=15,xlim=c(0,2),main=paste("P = ",hp),xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
    text(1,1,overl)
    text(0.55,1,tot_chi-overl)
    text(1.45,1,tot_wilx-overl)
    text(1.9,0.65,sum(as.numeric(a$pass_chi=="FALSE" & a$pass_wilx=="FALSE")))
    quartz.save(paste(cl,"_venn_alltrips.pdf",sep=""), type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
    dev.off()
}

hyper_subsample=function(tabl,size=9,cl)
{
    sub_tabl=tabl[tabl$clade==cl,]
    sps=unique(c(sub_tabl$P1out,sub_tabl$P2out,sub_tabl$P3out))
    p_vals=c()
    for(i in 1:10000)
    {
        
        sps_sub=sample(sps,size=size)
        a=sub_tabl[sub_tabl$P1out %in% sps_sub & sub_tabl$P2out %in% sps_sub & sub_tabl$P3out %in% sps_sub,]
        tot=nrow(a)
        overl=sum(as.numeric(a$pass_chi!="FALSE" & a$pass_wilx!="FALSE"))
        tot_chi=sum(as.numeric(a$pass_chi=="TRUE"))
        tot_wilx=sum(as.numeric(a$pass_wilx=="TRUE"))
        hyper_p=phyper(overl, tot_wilx, tot - tot_wilx, tot_chi, lower.tail = FALSE)
        p_vals=c(p_vals,hyper_p)
        
    }    
    return(p_vals)
    
}    

quartz(width=7, height=7)
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
quartz.save("Power.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
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

#Agreement between BLT and Chi-square
print_save_matrix(sp_space,m_overlap,"_blt_chi")
#Chi-square
print_save_matrix(sp_space,m_ch,"_chi")
#BLT
print_save_matrix(sp_space,m_wilx,"_blt")
#QuibL
print_save_matrix(sp_space,total_q,"_quibl")

#####################################Time-introgression plot#########################################


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

#####################################Number of introgression events#########################################


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


bc1=get_intomatrix_blt(C1,b_chisq)
bc2=get_intomatrix_blt(C2,b_chisq)
bc3=get_intomatrix_blt(C3,b_chisq)
bc4=get_intomatrix_blt(C4,b_chisq)
bc5=get_intomatrix_blt(C5,b_chisq)
bc6=get_intomatrix_blt(C6,b_chisq)
bc7=get_intomatrix_blt(C7,b_chisq)
bc8=get_intomatrix_blt(C8,b_chisq)
bc9=get_intomatrix_blt(C9,b_chisq)

bw1=get_intomatrix_blt(C1,b_wilcox)
bw2=get_intomatrix_blt(C2,b_wilcox)
bw3=get_intomatrix_blt(C3,b_wilcox)
bw4=get_intomatrix_blt(C4,b_wilcox)
bw5=get_intomatrix_blt(C5,b_wilcox)
bw6=get_intomatrix_blt(C6,b_wilcox)
bw7=get_intomatrix_blt(C7,b_wilcox)
bw8=get_intomatrix_blt(C8,b_wilcox)
bw9=get_intomatrix_blt(C9,b_wilcox)







f=1
for (cl in sp_space)
{    
    assign(paste("bw",f,sep=""),get_intomatrix_blt(cl,b_wilcox))
    assign(paste("bwtrue",f,sep=""),get_intomatrix_blt(cl,b_wilcox[b_wilcox$totalbl_pass=="TRUE",]))
    assign(paste("bwfalse",f,sep=""),get_intomatrix_blt(cl,b_wilcox[b_wilcox$totalbl_pass=="FALSE",]))
    f=f+1
}    

grid.arrange(bw1,bwtrue1,bwfalse1,bw2,bwtrue2,bwfalse2,bw3,bwtrue3,bwfalse3,bw4,bwtrue4,bwfalse4,ncol=3,nrow=4)
quartz.save("C1_C4_bw_split.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)
grid.arrange(bw5,bwtrue5,bwfalse5,bw6,bwtrue6,bwfalse6,bw7,bwtrue7,bwfalse7,bw8,bwtrue8,bwfalse8,bw9,bwtrue9,bwfalse9,ncol=3,nrow=5)
quartz.save("C5_C9_bw_split.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)





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






quartz(width=12, height=9)
grid.arrange(bc1,bw1,bc2,bw2,bc3,bw3,bc4,bw4,bc5,bw5,ncol=2,nrow=5)
quartz.save("C1_C5_bc_bw.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)
dev.off()
quartz(width=12, height=7.2)
grid.arrange(bc6,bw6,bc7,bw7,bc8,bw8,bc9,bw9,ncol=2,nrow=4)
quartz.save("C6_C9_bc_bw.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)
dev.off()


grid.arrange(b1,b2,b3,b4,ncol=1,nrow=4)
quartz.save("C1_C4_blt.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)
grid.arrange(b5,b6,b7,b8,b9,ncol=1,nrow=5)
quartz.save("C5_C9_blt.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)



grid.arrange(d1,q1,qc1,b1,d2,q2,qc2,b2,d3,q3,qc3,b3,d4,q4,qc4,b4,ncol=4,nrow=4)
grid.arrange(d5,q5,qc5,b5,d6,q6,qc6,b6,d7,q7,qc7,b7,d8,q8,qc8,b8,d9,q9,qc9,b9,ncol=4,nrow=5)

quartz.save("C1_C4_dfoil_quibl.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)
quartz.save("C5_C9_dfoil_quibl.png", type = "png",antialias=F,bg="white",dpi=400,pointsize=12)