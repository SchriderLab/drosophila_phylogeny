library("ape")
library('ggplot2')
library('pals')
library('reshape2')
library('dplyr')
library('gridExtra')
library("svMisc")
library("svMisc")


#Disable scientific notation for branch length 
.write.tree2 <- function(phy, digits = 10, tree.prefix = "", check_tips)
{
    brl <- !is.null(phy$edge.length)
    nodelab <- !is.null(phy$node.label)
    #if (check_tips) phy$tip.label <- checkLabel(phy$tip.label)
    #if (nodelab) phy$node.label <- checkLabel(phy$node.label)
    f.d <- paste("%.", digits, "f", sep = "")
    cp <- function(x){
        STRING[k] <<- x
        k <<- k + 1
    }
    add.internal <- function(i) {
        cp("(")
        desc <- kids[[i]]
        for (j in desc) {
            if (j > n) add.internal(j)
            else add.terminal(ind[j])
            if (j != desc[length(desc)]) cp(",")
        }
        cp(")")
        if (nodelab && i > n) cp(phy$node.label[i - n]) # fixed by Naim Matasci (2010-12-07)
        if (brl) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[ind[i]]))
        }
    }
    add.terminal <- function(i) {
        cp(phy$tip.label[phy$edge[i, 2]])
        if (brl) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[i]))
        }
    }

    n <- length(phy$tip.label)

    ## borrowed from phangorn:
    parent <- phy$edge[, 1]
    children <- phy$edge[, 2]
    kids <- vector("list", n + phy$Nnode)
    for (i in 1:length(parent))
        kids[[parent[i]]] <- c(kids[[parent[i]]], children[i])

    ind <- match(1:max(phy$edge), phy$edge[, 2])

    LS <- 4*n + 5
    if (brl) LS <- LS + 4*n
    if (nodelab)  LS <- LS + n
    STRING <- character(LS)
    k <- 1
    cp(tree.prefix)
    cp("(")
    getRoot <- function(phy)
        phy$edge[, 1][!match(phy$edge[, 1], phy$edge[, 2], 0)][1]
    root <- getRoot(phy) # replaced n+1 with root - root has not be n+1
    desc <- kids[[root]]
    for (j in desc) {
        if (j > n) add.internal(j)
        else add.terminal(ind[j])
        if (j != desc[length(desc)]) cp(",")
    }

    if (is.null(phy$root.edge)) {
        cp(")")
        if (nodelab) cp(phy$node.label[1])
        cp(";")
    }
    else {
        cp(")")
        if (nodelab) cp(phy$node.label[1])
        cp(":")
        cp(sprintf(f.d, phy$root.edge))
        cp(";")
    }
    paste(STRING, collapse = "")
}

assignInNamespace(".write.tree2", .write.tree2, "ape")


keep.tip.new=function(phy, tip) 
{
    if (!inherits(phy, "phylo")) 
        stop("object \"phy\" is not of class \"phylo\"")
    Ntip <- length(phy$tip.label)
    if (is.character(tip)) 
    {
        idx <- match(tip, phy$tip.label)
        if (!anyNA(idx)) 
        {
            tip <- idx
            toDrop <- setdiff(1:Ntip, tip)
            drop.tip(phy, toDrop)
        }
    }
}



#QuIBL input generator 

all_triplets=function(taxa_list,gene_trees,dir_name)
{
    dir.create(dir_name)
    setwd(dir_name)
    phy=read.tree(gene_trees)
    taxa_combn=combn(taxa_list,m=3)
    for (i in 1:ncol(taxa_combn))
    {
        triplet=taxa_combn[,i]
        sub_triplet=lapply(phy,keep.tip.new, c(triplet,"Anopheles_gambiae"))
        sub_triplet=Filter(Negate(is.null),sub_triplet)
        config="
[Input]
treefile: ./INNAME
numdistributions: 2
likelihoodthresh: 0.01
numsteps: 50
gradascentscalar: 0.5
totaloutgroup: Anopheles_gambiae
multiproc: True
maxcores:1000
[Output]
OutputPath: ./OUTNAME" 
        config=gsub("INNAME",paste("quibltrees_",i,sep=""),config)
        config=gsub("OUTNAME",paste("quibl_out",i,sep=""),config)
        write(config,paste("quiblin_",i,sep=""))
        write(unlist(lapply(sub_triplet,write.tree)),paste("quibltrees_",i,sep=""))
    
    }    
    setwd("..")
}    

#################################### GENERATE INPUT ###############################

tt=read.tree("MLrooted.tre")
node_n=c(168,185,204,226,245,256,265,289,307)

all_triplets(extract.clade(tt,168)$tip.label,"/Users/Anton/Downloads/gene_trees_wboot_dna_mafft",paste("C",1,"_node",168,"_quibl",sep=""))
all_triplets(extract.clade(tt,185)$tip.label,"/Users/Anton/Downloads/gene_trees_wboot_dna_mafft",paste("C",2,"_node",185,"_quibl",sep=""))
all_triplets(extract.clade(tt,204)$tip.label,"/Users/Anton/Downloads/gene_trees_wboot_dna_mafft",paste("C",3,"_node",204,"_quibl",sep=""))
all_triplets(extract.clade(tt,226)$tip.label,"/Users/Anton/Downloads/gene_trees_wboot_dna_mafft",paste("C",4,"_node",226,"_quibl",sep=""))
all_triplets(extract.clade(tt,245)$tip.label,"/Users/Anton/Downloads/gene_trees_wboot_dna_mafft",paste("C",5,"_node",245,"_quibl",sep=""))
all_triplets(extract.clade(tt,256)$tip.label,"/Users/Anton/Downloads/gene_trees_wboot_dna_mafft",paste("C",6,"_node",256,"_quibl",sep=""))
all_triplets(extract.clade(tt,265)$tip.label,"/Users/Anton/Downloads/gene_trees_wboot_dna_mafft",paste("C",7,"_node",265,"_quibl",sep=""))
all_triplets(extract.clade(tt,289)$tip.label,"/Users/Anton/Downloads/gene_trees_wboot_dna_mafft",paste("C",8,"_node",289,"_quibl",sep=""))
all_triplets(extract.clade(tt,307)$tip.label,"/Users/Anton/Downloads/gene_trees_wboot_dna_mafft",paste("C",9,"_node",307,"_quibl",sep=""))



################################## ANALYSES OF THE RESULTS #########################

total_q=read.csv("/Users/Anton/Downloads/droso_quibl_results.txt",header=F,stringsAsFactors=F)
names(total_q)=c("clade","id","P1","P2","P3","outgroup","Com1","Com2","mixprop1","mixprop2","lambda2Dist","lambda1Dist","BIC2Dist","BIC1Dist","count")
total_q$triplid=as.character(apply(total_q[,c("P1","P2","P3")],1,paste,collapse="_"))
total_q$BICdiff = total_q$BIC2-total_q$BIC1

total_b=read.csv("/Users/Anton/Downloads/droso_blt_results.txt",stringsAsFactors=FALSE,header=F)
names(total_b)=c("clade","P1out","P2out","P3out","CountP1","CountP2","CountP3","PvalueChi","meanT_concord","meanT_discord1","meanT_discord2","PvalueWCOMC1","PvalueWCOMC2","PvalueWC1C2")
total_b$PvalueChi=p.adjust(total_b$PvalueChi,method="fdr")
total_b$PvalueWCOMC1=p.adjust(total_b$PvalueWCOMC1,method="fdr") 
total_b$PvalueWCOMC2=p.adjust(total_b$PvalueWCOMC2,method="fdr") 
total_b$PvalueWC1C2=p.adjust(total_b$PvalueWC1C2,method="fdr") 



total_b_f=c()
for (n in 1:nrow(total_q))
{
    
    progress(n,nrow(total_q))
    sp=as.character(unlist(total_q[n,c("P1","P2","P3")]))
    to_select=apply(apply(total_b[,c("P1out","P2out","P3out")],1,"%in%",sp),2,all)
   
    if (any(to_select))
    {    
        intropair=sp[!sp %in% total_q[n,"outgroup"]]
        total_b_v=total_b[to_select,]
        ord_means=total_b_v[,c("meanT_concord","meanT_discord1","meanT_discord2")][order(total_b_v[,c("CountP1","CountP2","CountP3")],decreasing = T)]
        names(ord_means)=c("meanP1","meanP2","meanP3")
        total_b_v=cbind(total_b_v,ord_means)
        pos=which(total_b_v[c("P1out","P2out","P3out")]==total_q[n,"outgroup"])
        Count=as.numeric(total_b_v[,c("CountP1","CountP2","CountP3")][pos])
        MeanT=as.numeric(total_b_v[,c("meanP1","meanP2","meanP3")][pos])
        ord=order(total_b_v[,c("CountP1","CountP2","CountP3")])[pos]
        
        if (ord==3)
        {
            type="common"
            wilcoxP=NA
            chisqP=NA
            chisq_introg=NA
            blt_introg_st=NA
            blt_introg_realx=NA
            

        } else if (ord==1) {

            type="discord1"
            wilcoxP=total_b_v$PvalueWC1C2
            chisqP=NA
            chisq_introg=NA
            blt_introg_st=MeanT < total_b_v[,c("meanP1","meanP2","meanP3")][order(total_b_v[,c("CountP1","CountP2","CountP3")])==3] & MeanT < total_b_v[,c("meanP1","meanP2","meanP3")][order(total_b_v[,c("CountP1","CountP2","CountP3")])==2] & total_b_v$PvalueWC1C2 < 0.05
            blt_introg_realx = MeanT < total_b_v[,c("meanP1","meanP2","meanP3")][order(total_b_v[,c("CountP1","CountP2","CountP3")])==2] & total_b_v$PvalueWC1C2 < 0.05
    
        } else {

            type="discord2"
            wilcoxP=total_b_v$PvalueWC1C2
            chisqP=total_b_v$PvalueChi
            chisq_introg=total_b_v$PvalueChi<0.05
            blt_introg_st=MeanT < total_b_v[,c("meanP1","meanP2","meanP3")][order(total_b_v[,c("CountP1","CountP2","CountP3")])==3]  & MeanT < total_b_v[,c("meanP1","meanP2","meanP3")][order(total_b_v[,c("CountP1","CountP2","CountP3")])==1] & total_b_v$PvalueWC1C2 < 0.05
            blt_introg_realx= MeanT < total_b_v[,c("meanP1","meanP2","meanP3")][order(total_b_v[,c("CountP1","CountP2","CountP3")])==1] & total_b_v$PvalueWC1C2 < 0.05
        }
        total_b_f=rbind(total_b_f,c(total_q[n,],MeanT,wilcoxP,chisqP,Count,type,intropair,chisq_introg,blt_introg_st,blt_introg_realx))
     } else { 
        
        print("No ILS/introgression")
        
    }   
    
   
      
}

total_b_f=as.data.frame(total_b_f)

names(total_b_f)=c('clade','id','P1','P2','P3','outgroup','Com1','Com2','mixprop1','mixprop2','lambda2Dist','lambda1Dist','BIC2Dist','BIC1Dist','count','triplid','BICdiff',"MeanT","wilcoxP","chisqP","Count","type","i1","i2","chisq_introg","blt_introg_st","blt_introg_relax")

total_b_f$quibl_introg=total_b_f$BICdiff< -30

total_b_f$intersect_introg_all=apply(total_b_f[,c("chisq_introg","quibl_introg","blt_introg_relax")],1,all)
total_b_f$intersect_introg_chiblt=apply(total_b_f[,c("chisq_introg","blt_introg_relax")],1,all)
total_b_f$clade_id=unlist(lapply(lapply((lapply(total_b_f[,1],strsplit,"_")),"[[",1),"[",1))


get_intomatrix=function(taxa,dat,sig)
{
    m=matrix(0, nrow=length(taxa),ncol=length(taxa)) 
    m=data.frame(m)
    colnames(m)=taxa
    rownames(m)=taxa
    dat_int=dat[ dat$type!="common" & sig =="TRUE" & (dat$i1 %in% taxa | dat$i2 %in% taxa),c("i1","i2")]
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
    g1=ggplot(melted_m, aes(Var2, Var1, fill = value))+geom_tile(color = "white")+scale_fill_gradientn(colors=jet(100),space = "Lab", limits=c(0,1)) +theme_minimal()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+coord_fixed()+labs(x="",y="")
    return(g1)
}   



hygtest=function(cl,dat)
{
    sset=dat[dat$clade_id==cl & dat$type!="common" ,c("chisq_introg","blt_introg_relax","intersect_introg_chiblt")]
    group_chi=nrow(sset[sset$chisq_introg=="TRUE",])
    group_blt=nrow(sset[sset$blt_introg_relax=="TRUE",])
    over=nrow(sset[sset$intersect_introg_chiblt=="TRUE",])
    total=nrow(sset)
    pval=phyper(over,group_blt,total-group_blt,group_chi,lower.tail= FALSE)
    cat(cl,pval,group_chi,group_blt,total,over,fill=TRUE)
}    



for (cl in c("C1","C2","C3","C4","C5","C6","C7","C8","C9"))
{
   hygtest(cl,total_b_f) 
}    











f=1
for (cl in sp_space)
{    
    assign(paste("c",f,sep=""),get_intomatrix(cl,total_b_f,unlist(total_b_f$chisq_introg)))
    assign(paste("q",f,sep=""),get_intomatrix(cl,total_b_f,unlist(total_b_f$quibl_introg)))
    assign(paste("b",f,sep=""),get_intomatrix(cl,total_b_f,unlist(total_b_f$blt_introg_relax)))
    assign(paste("a",f,sep=""),get_intomatrix(cl,total_b_f,unlist(total_b_f$intersect_introg_chiblt)))
    f=f+1
}    

grid.arrange(c1,q1,b1,a1,c2,q2,b2,a2,c3,q3,b3,a3,c4,q4,b4,a4,ncol=4,nrow=4)
grid.arrange(c5,q5,b5,a5,c6,q6,b6,a6,c7,q7,b7,a7,c8,q8,b8,a8,c9,q9,b9,a9,ncol=4,nrow=5)


library("VennDiagram")