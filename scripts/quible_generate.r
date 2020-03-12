library("ape")

quible_trees=function(triplet,gene_trees)
{
    phy=read.tree(gene_trees)
    trees=c()
    for (gt in phy)
    {
        if (all(triplet %in% gt$tip.label))
        {
            gt$node.label=NULL
            gt_sub=keep.tip(gt,triplet)
            trees=c(trees,write.tree(gt_sub))
        }    
        
    }    
    write.table(trees,"quiblin", quote = FALSE,row.names = FALSE, col.names=FALSE,append=TRUE)
    
}    

C16=c('D_ananassae','D_parapallidosa','D_bipectinata','D_parabipectinata','D_m_pallens','D_malerkotliana_malerkotliana','D_p_nigrens','D_pseudoananassae_pseudoananassae','D_ercepeae')

C16_ret=c("Anopheles_gambiae","D_p_nigrens","D_m_pallens","D_parabipectinata")
quible_trees(C16_ret,"gene_trees_wboot_dna_mafft")