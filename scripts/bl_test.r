library("ape")
library("phangorn")
library("svMisc")

args = commandArgs(trailingOnly=TRUE)

test_triplet=function(taxa,gene_trees,clade_name)
{
    busco_len=gene_trees[,1:2]
    gene_trees=read.tree(text=gene_trees[,"V3"]) 
    trl_all=c()
    brls_all1=c()
    brls_all2=c()
    internal_all=c()
    out_all=c()
    root_tip_all=c()
    busco_id_all=c()
    aln_length_all=c()
    ind=0
    for (tre in gene_trees)
    {
        
        ind=ind+1
        if("M_domestica" %in% tre$tip.label & all(taxa %in% tre$tip.label))
        {
            
            tre=root(tre,"M_domestica")
            trl=sum(tre$edge.length)
            tre_trip=keep.tip(tre,taxa)
            brls=extract.clade(tre_trip,max(tre_trip$edge))$edge.length
            root_tip=tre_trip$tip.label[!tre_trip$tip.label %in% extract.clade(tre_trip,max(tre_trip$edge))$tip.label]
            trl_all=c(trl_all,trl)
            out_all=c(out_all,tre_trip$edge.length[1])
            internal_all=c(internal_all,tre_trip$edge.length[2])
            brls_all1=c(brls_all1,brls[1])
            brls_all2=c(brls_all2,brls[2])
            root_tip_all=c(root_tip_all,root_tip)
            busco_id_all=c(busco_id_all,busco_len[ind,1])
            aln_length_all=c(aln_length_all,busco_len[ind,2])
            
           
        }    
    }
    counts=table(root_tip_all)
    con=names(which.max(counts))
    dis=names(which.min(counts))
    m=data.frame(brl1=brls_all1,brl2=brls_all2,trl=trl_all,brl_out=out_all,brl_int=internal_all,root_tip=root_tip_all,topo=ifelse(root_tip_all %in% con,"concord",ifelse(root_tip_all %in% dis,"discord2","discord1")))
    write.table(m,paste(c(taxa,"csv"),collapse="."),quote=F,row.names=F)
    if(!any(table(m$root_tip)==0) & length(table(m$root_tip))==3)
    {
        counts=table(m$root_tip)
        com=names(which.max(counts))
        not_com=names(counts)[!names(counts) %in% com]
        not_com_c=counts[!names(counts) %in% com]
        m$common=ifelse(m$root_tip==com,"TRUE","FALSE")
        m$proxy_t=(m$brl1+m$brl2)/m$trl
        ccom=m[m$common==TRUE,"proxy_t"]
        c1=m[m$root_tip==not_com[1],"proxy_t"]
        c2=m[m$root_tip==not_com[2],"proxy_t"]
        w_testc1=wilcox.test(ccom,c1)$p.value
        w_testc2=wilcox.test(ccom,c2)$p.value
        w_test=wilcox.test(c1,c2)$p.value
        chi=chisq.test(not_com_c)$p.value
        v_out=c(clade_name,names(counts),counts,chi,mean(ccom),mean(c1),mean(c2),w_testc1,w_testc2,w_test)
        return(as.vector(v_out))
    } 
  

}    


getstats_triplets=function(taxa_list,gene_trees,clade_name)
{
   
    taxa_combn=combn(taxa_list,m=3)
    out_t=c()
    for (i in 1:ncol(taxa_combn))
    {
        progress(i,ncol(taxa_combn))
        triplet=taxa_combn[,i]
        stats=test_triplet(triplet,gene_trees,clade_name)
        out_t=rbind(out_t,stats) 
    }
    write.table(as.data.frame(out_t),clade_name,quote = F, row.names = F, col.names = F,sep=",")
    
}    



tt=read.tree(args[1])
phy=read.table(args[2],stringsAsFactors = F)
#Arg 3 = the node number
clade=extract.clade(tt,as.numeric(args[3]))$tip.label
#Arg 4 = clade name
getstats_triplets(clade,phy,args[4])




