library("ape")


dfoil_select=function(phy,id)
{
    dfoil_out=c()
    taxa_combn=combn(phy$tip.label,m=5)
    test_phy=read.tree(text="((('A','A'),('A','A')),A);")
    for (i in 1:ncol(taxa_combn))
    {
        sub_phy=keep.tip(phy,taxa_combn[,i])
        sub_phy_dup=sub_phy
        sub_phy_dup$tip.label=rep("A",length(sub_phy_dup$tip.label))
        if (all.equal.phylo(sub_phy_dup,test_phy,use.edge.length=F,use.tip.label=F))
        {
            node_number_order=rev(order(node.depth.edgelength(sub_phy)[(length(sub_phy$tip.label)+1):length(node.depth.edgelength(sub_phy))])+length(sub_phy$tip.label))[1:2]
            taxa_order=c()
            for (n in node_number_order)
            {
               node_taxa=extract.clade(sub_phy,n)$tip.label
               taxa_order=c(taxa_order,node_taxa) 
            }
            dfoil_order=c(taxa_order,sub_phy$tip.label[!sub_phy$tip.label %in% taxa_order])
            dfoil_out=c(dfoil_out,paste(dfoil_order,collapse=","))
            
        }    
        
    }
    write(dfoil_out,id)
}    


dfoil_select_all=function(phy,id)
{
    dfoil_out=c()
    taxa_combn=combn(phy$tip.label,m=5)
    taxa_combn=data.frame(t(taxa_combn))
    write(apply(taxa_combn,1,paste,collapse=","),id)
}    



tt=read.tree("MLrooted.tre")
nodelabels()
node_n=c(168,185,204,226,245,256,265,289,307)



cl_id=1
for (i in node_n)
{
   dfoil_select(extract.clade(tt,i),paste("C",cl_id,"_node",i,"_dfoil",sep=""))
   cl_id=cl_id+1 
}    


#All clades 
cl_id=1
for (i in node_n)
{
   dfoil_select_all(extract.clade(tt,i),paste("C",cl_id,"_node",i,"_dfoil_all",sep=""))
   cl_id=cl_id+1 
}    