library("ape")
library("phangorn")
library("svMisc")

args = commandArgs(trailingOnly=TRUE)

test_triplet=function(taxa,gene_trees)
{
    trl_all=c()
    brls_all1=c()
    brls_all2=c()
    root_tip_all=c()
    for (tre in gene_trees)
    {
        
        if("Anopheles_gambiae" %in% tre$tip.label & all(taxa %in% tre$tip.label))
        {
            
            trl=sum(tre$edge.length)
            tre_trip=keep.tip(tre,taxa)
            brls=extract.clade(tre_trip,max(tre_trip$edge))$edge.length
            root_tip=tre_trip$tip.label[!tre_trip$tip.label %in% extract.clade(tre_trip,max(tre_trip$edge))$tip.label]
            trl_all=c(trl_all,trl)
            brls_all1=c(brls_all1,brls[1])
            brls_all2=c(brls_all2,brls[2])
            root_tip_all=c(root_tip_all,root_tip)
            
        }    
    }    
    m=data.frame(brl1=brls_all1,brl2=brls_all2,trl=trl_all,root_tip=root_tip_all)
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
        w_testc2=wilcox.test(ccom,c1)$p.value
        w_test=wilcox.test(c1,c2)$p.value
        chi=chisq.test(not_com_c)$p.value
        v_out=c(names(counts),counts,chi,mean(ccom),mean(c1),mean(c2),w_testc1,w_testc2,w_test)
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
        stats=test_triplet(triplet,gene_trees)
        out_t=rbind(out_t,stats) 
    }
    write.table(as.data.frame(out_t),clade_name,quote = F, row.names = F, col.names = F,sep=",")
    
}    



tt=read.tree(args[1])
phy=read.tree(args[2])
clade=extract.clade(tt,as.numeric(args[3]))$tip.label

getstats_triplets(clade,phy,args[4])




