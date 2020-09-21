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
#total_b$PvalueChi=p.adjust(total_b$PvalueChi,method="fdr")
#total_b$PvalueWCOMC1=p.adjust(total_b$PvalueWCOMC1,method="fdr") 
#total_b$PvalueWCOMC2=p.adjust(total_b$PvalueWCOMC2,method="fdr") 
#total_b$PvalueWC1C2=p.adjust(total_b$PvalueWC1C2,method="fdr") 

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

#Stringent unique introgression (main figure)
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
        sig_chi=0
        sig_wilx=0
        no_sig=0
        sps_sub=sample(sps,size=size)
        a=sub_tabl[sub_tabl$P1out %in% sps_sub & sub_tabl$P2out %in% sps_sub & sub_tabl$P3out %in% sps_sub,]
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
        hp=phyper(overl, tot_wilx, tot - tot_wilx, tot_chi, lower.tail = FALSE)
        p_vals=c(p_vals,hp)
        
    }    
    return(p_vals)
    
}    

quartz(width=3, height=27)
par(mfrow=c(9,1))
i=1
for (cl in c("C1","C2","C3","C4","C5","C6","C7","C8","C9"))
{ 
    pval=hyper_subsample_uniq(total_b,8,cl)
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
            overl=overl+1
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
        sig_chi=0
        sig_wilx=0
        no_sig=0
        sps_sub=sample(sps,size=size)
        a=sub_tabl[sub_tabl$P1out %in% sps_sub & sub_tabl$P2out %in% sps_sub & sub_tabl$P3out %in% sps_sub,]
        for (p in unique(a$pair_name))
        {
            if (any(a$pair_name==p & a$pass_chi!="FALSE" & a$pass_wilx!="FALSE"))
            {
                overl=overl+1 
            } else if (all(c(any(a$pair_name==p & a$pass_chi!="FALSE"),any(a$pair_name==p & a$pass_wilx!="FALSE")))) {
                overl=overl+1
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
        hp=phyper(overl, tot_wilx, tot - tot_wilx, tot_chi, lower.tail = FALSE)
        p_vals=c(p_vals,hp)
        
    }    
    return(p_vals)
    
}    

quartz(width=3, height=27)
par(mfrow=c(9,1))
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
        overl=sum(as.numeric(a$pass_chi!="FALSE" & a$pass_wilx!="FALSE"))
        tot_chi=sum(as.numeric(a$pass_chi=="TRUE"))
        tot_wilx=sum(as.numeric(a$pass_wilx=="TRUE"))
        hyper_p=phyper(overl, tot_wilx, tot - tot_wilx, tot_chi, lower.tail = FALSE)
        p_vals=c(p_vals,hyper_p)
        
    }    
    return(p_vals)
    
}    


#Introgression power test all triplets
quartz(width=3, height=27)
par(mfrow=c(9,1))
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

#Agreement between BLT and Chi-square
print_save_matrix(sp_space,m_overlap,"_blt_chi")
#Chi-square
print_save_matrix(sp_space,m_ch,"_chi_nofdr")
#BLT
print_save_matrix(sp_space,m_wilx,"_blt_nofdr")
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




