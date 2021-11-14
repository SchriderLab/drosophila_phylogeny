#!/usr/bin/env Rscript
library("optparse")
library("doSNOW")

#Arguments
option_list = list(
  make_option(c("-t", "--trees"), type="character", default=NULL,help="input gene trees in newick", metavar="character"),
  make_option(c("-s", "--species_tree"), type="character", default=NULL, help="rooted species tree in newick", metavar="character"),
  make_option(c("-n", "--node"), type="numeric", default=NULL, help="node number of a tested clade in a species tree", metavar="numeric"),
  make_option(c("-c", "--cores"), type="numeric", default=NULL, help="number of cores for parallel computing", metavar="numeric"),
  make_option(c("-p", "--prefix"), type="character", default=NULL, help="clade prefix", metavar="character"),
  make_option(c("-o", "--outgroup"), type="character", default=NULL, help="outgroup species (only one allowed)", metavar="character")  
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#Create cluster
cl = makeCluster(opt$cores,type = "SOCK") 

library("phangorn")
library("foreach")
library("ape")

args = commandArgs(trailingOnly=TRUE)

d3_stat=function(taxa,species_tree,pw_distance,outg)
{

    tre_trip=keep.tip(species_tree,taxa)
    root_tip=tre_trip$tip.label[!tre_trip$tip.label %in% extract.clade(tre_trip,max(tre_trip$edge))$tip.label]
    ingr=sort(tre_trip$tip.label[tre_trip$tip.label!=root_tip])
    dis1=pw_distance[row.names(pw_distance)==ingr[1],root_tip]
    dis2=pw_distance[row.names(pw_distance)==ingr[2],root_tip]
    d3=(dis1-dis2)/(dis1+dis2)
    tre_trip$edge.length=NULL
    v_out=c(ingr[1],ingr[2],root_tip,d3,write.tree(tre_trip))
    return(as.vector(v_out))
    
}


test_triplet=function(taxa,gene_trees,clade_name,outg)
{
    gene_trees=gene_trees 
    trl_all=c()
    brls_all1=c()
    brls_all2=c()
    internal_all=c()
    out_all=c()
    root_tip_all=c()
    for (tre in gene_trees)
    {
        if(outg %in% tre$tip.label & all(taxa %in% tre$tip.label))
        {
            
            tre=root(tre,outg)
            trl=sum(tre$edge.length)
            tre_trip=keep.tip(tre,taxa)
            brls=extract.clade(tre_trip,max(tre_trip$edge))$edge.length
            root_tip=tre_trip$tip.label[!tre_trip$tip.label %in% extract.clade(tre_trip,max(tre_trip$edge))$tip.label]
            trl_all=c(trl_all,trl)
            outl=tre_trip$edge.length[which(tre_trip$edge[,1]==4 & (tre_trip$edge[,2]==1 | tre_trip$edge[,2]==2 | tre_trip$edge[,2]==3))]
            out_all=c(out_all,outl)
            internall=tre_trip$edge.length[which(tre_trip$edge[,1]==4 & tre_trip$edge[,2]==5)]
            internal_all=c(internal_all,internall)
            brls_all1=c(brls_all1,brls[1])
            brls_all2=c(brls_all2,brls[2])
            root_tip_all=c(root_tip_all,root_tip)
                   
        }    
    }
    #Get counts/ compute common and non-common
    counts=table(root_tip_all)
    com=names(which.max(counts))
    dis=names(which.min(counts))
    m=data.frame(clade_name,P1=tre_trip$tip.label[1],P2=tre_trip$tip.label[2],P3=tre_trip$tip.label[3],brl1=brls_all1,brl2=brls_all2,trl=trl_all,brl_out=out_all,brl_int=internal_all,root_tip=root_tip_all,topo=ifelse(root_tip_all %in% com,"concord",ifelse(root_tip_all %in% dis,"discord1","discord2")))
    #write.table(m,paste(c(taxa,"csv"),collapse="."),quote=F,row.names=F,col.names=F)
    if(!any(table(m$root_tip)==0) & length(table(m$root_tip))==3)
    {
        #Normalize ditances by total gene tree lengths
        m$proxy_t=(m$brl1+m$brl2)/m$trl
        ccom=m[m$topo=="concord","proxy_t"]
        c1=m[m$topo=="discord1","proxy_t"]
        c2=m[m$topo=="discord2","proxy_t"]
        not_com_count=table(as.vector(m[m$topo!="concord","topo"]))
        w_testc1=wilcox.test(ccom,c1)$p.value
        w_testc2=wilcox.test(ccom,c2)$p.value
        w_test=wilcox.test(c1,c2)$p.value
        chi=chisq.test(not_com_count)$p.value
        v_out=c(clade_name,names(counts),counts,chi,mean(ccom),mean(c1),mean(c2),w_testc1,w_testc2,w_test)
        return(as.vector(v_out))
    } 
}    


getstats_triplets=function(taxa_list,gene_trees,clade_name,outg,species_tree,pw_distance)
{
   
    taxa_combn=combn(taxa_list,m=3)
    print(paste("N triplets:",ncol(taxa_combn)))
    pb=txtProgressBar(0,ncol(taxa_combn),style=3)
    progress=function(n){
    setTxtProgressBar(pb,n)
    }
    opts=list(progress=progress)
    out_t=foreach(i=1:ncol(taxa_combn),.combine='rbind',.options.snow=opts) %dopar% 
    {
       
        triplet=taxa_combn[,i]
        stats1=test_triplet(triplet,gene_trees,clade_name,outg)
        stats2=d3_stat(triplet,species_tree,pw_distance,outg)
        return(c(stats1,stats2))
        
    }
    write.table(as.data.frame(out_t),clade_name,quote = F, row.names = F, col.names = F,sep=",")
    
}    

clusterExport(cl, c("read.tree","cophenetic.phylo","root","keep.tip","extract.clade","setTxtProgressBar","test_triplet","drop.tip","write.tree","d3_stat"))
registerDoSNOW(cl)

tt=read.tree(opt$species_tree)
cat("Read species topology. Done.\n")
pw=cophenetic.phylo(tt)
cat("Calculate pairwise distances from species tree. Done.\n")
phy=read.tree(opt$trees)
cat("Read gene trees. Done.\n")
clade=extract.clade(tt,as.numeric(opt$node))
cat("Exatract clade. Done.\n")
cat(paste(clade$tip.label,collapse="\n"))
cat("\nCalculating DCT/BLT.\n")
getstats_triplets(clade$tip.label,phy,opt$prefix,opt$outgroup,clade,pw)
cat("\nDone.\n") 
stopCluster(cl)




