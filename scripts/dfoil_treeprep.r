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


################################################################ RUN DFOIL ######################################################################## 


names_v=c("clade","P1","P2","P3","P4","Out","chrom1","position", "AAAAA" , "AAABA" , "AABAA" , "AABBA" , "ABAAA" , "ABABA" , "ABBAA" , "ABBBA" , "BAAAA" , "BAABA" , "BABAA" , "BABBA" , "BBAAA" ,"BBABA","BBBAA","BBBBA",'chromdup','coord','total','dtotal','T12','T34','T1234','DFO_left','DFO_right','DFO_total','DFO_stat','DFO_chisq','DFO_Pvalue','DIL_left','DIL_right','DIL_total','DIL_stat','DIL_chisq','DIL_Pvalue','DFI_left','DFI_right','DFI_total','DFI_stat','DFI_chisq','DFI_Pvalue','DOL_left','DOL_right','DOL_total','DOL_stat','DOL_chisq','DOL_Pvalue','introgression','introgna','intrognone','introg13','introg14','introg23','introg24','introg31','introg41','introg32','introg42','introg123','introg124')

total=read.csv("/Users/Anton/Downloads/droso_dfoil_results.txt")
names(total)=names_v
total=total[total$T12<total$T34,]

total$Genus="Drosophila"
total$introgressionid=ifelse(total$introgression=="none","None",
                             ifelse(total$introgression=="123" | total$introgression=="124","Ancestral","Inter-group"))















total_order=melt(total[,c("introgressionid","Genus")])
g1=ggplot(total_order, aes(x=Genus, y=..count../sum(..count..),fill=introgressionid))+geom_bar(position="fill")+labs(x="Genus", y = "")+scale_fill_manual(values=c("#f0ad4e", "#5cb85c","#337ab7"),name="Introgression")
total_suborder=melt(total[,c("introgressionid","clade")])
g2=ggplot(total_suborder, aes(x=clade, y=..count../sum(..count..),fill=introgressionid))+geom_bar(position="fill")+geom_bar(position="fill")+labs(x="Clade", y = "Proportion of Quintets")+scale_fill_manual(values=c("#f0ad4e", "#5cb85c","#337ab7"),name="Introgression")

grid.arrange(g1,g2,nrow=3)