library("ape")

phylonet_full_topo=function(sp_tree,gene_trees,n_retic)
{
    sp_phy=read.tree(sp_tree)
    sp_phy$node.label=NULL
    phy=read.tree(gene_trees)
    rooted_phy=c()
    out_v=c("Anopheles_gambiae","L_varia","S_lebanonensis","C_costata")
    for (gt in phy)
    {
        
        if (any(out_v %in% gt$tip.label))
        {
            out_sp=out_v[which(out_v %in% gt$tip.label)[1]]
            gt_r=root(gt,out_sp,resolve.root=T)
            gt_r$node.label[gt_r$node.label=="Root"]=100
            gt_r$node.label[gt_r$node.label==""]=1
            gt_r$node.label=as.numeric(gt_r$node.label)/100
            rooted_phy=c(rooted_phy,write.tree(gt_r))
            
        } 
        
       
    } 
    phy=rooted_phy
    for (i in 1:n_retic)
    {    
        f_n=paste("phylonet_genes_",i,"ret.nex",sep="")
        write("#NEXUS\n\nBEGIN TREES;",f_n)
        write(paste("Tree fixtr = ",write.tree(sp_phy)),f_n,append=T)
        d=data.frame(rep("Tree",length(phy)),paste("gt",1:length(phy),"=",sep=""),phy)
        write.table(d,f_n,quote = F,row.names = F, col.names=F,append=T)
        write("END;\n\nBEGIN PHYLONET;",f_n,append=T)
        write(paste("InferNetwork_MPL (all)",i,"-s fixtr -fs -di -pl 15 -x 1;","\nEND;"),f_n,append=T)
    }    
}    

phylonet_full_topo("/Users/Anton/Downloads/BUSCO50_dna_pasta_nopart_iqtree_root.tre","/Users/Anton/Downloads/BUSCO50_dna_pasta_iqtree_all",1)