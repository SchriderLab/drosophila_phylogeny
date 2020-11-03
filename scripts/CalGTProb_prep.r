library('ape')

all_retictrees=function(retic_trees,gene_trees,dir_name)
{
    dir.create(dir_name)
    setwd(dir_name)
    retic_phy=read.tree(retic_trees)
    phy=read.tree(gene_trees)
    rooted_phy=c()
    out_v=c("Anopheles_gambiae","L_varia","S_lebanonensis","C_costata")
    sp_list=retic_phy[[1]]$tip.label[!retic_phy[[1]]$tip.label %in% out_v]
    for (gt in phy)
    {
        if (any(sp_list %in% gt$tip.label) & any(out_v %in% gt$tip.label) & sum(sp_list %in% gt$tip.label)>=3)
        {
              
            out_sp=out_v[which(out_v %in% gt$tip.label)[1]]
            gt_r=root(gt,out_sp,resolve.root=T)
            gt_sub=keep.tip(gt_r,c(sp_list[sp_list %in% gt$tip.label]))
            gt_sub$node.label[1]=100
            gt_sub$node.label[gt_sub$node.label=="0"]="1"
            gt_sub$node.label[gt_sub$node.label==""]="1"
            gt_sub$node.label[gt_sub$node.label=="NA"]="1"
            gt_sub$node.label[gt_sub$node.label=="Root"]=100
            gt_sub$node.label=as.numeric(gt_sub$node.label)/100
            rooted_phy=c(rooted_phy,write.tree(gt_sub))
               
        } 
        
    }
    phy=rooted_phy
    for (i in 1:length(retic_phy))
    {    
        f_n=paste("calgtprob_",i,"ret.nex",sep="")
        write("#NEXUS\n\nBEGIN NETWORKS;",f_n)
        write(paste("Network rettree = ",write.tree(retic_phy[[i]])),f_n,append=TRUE)
        write("END;\n\nBEGIN TREES;",f_n,append=TRUE)
        d=data.frame(rep("Tree",length(phy)),paste("gt",1:length(phy),"=",sep=""),phy)
        write.table(d,f_n, quote = FALSE,row.names = FALSE, col.names=FALSE,append=TRUE)
        write("END;\n\nBEGIN PHYLONET;",f_n,append=TRUE)
        write(paste("CalGTProb rettree (all) -o;","\nEND;"),f_n,append=TRUE)
        #return(phy)
    }        
    setwd("..")
}    

all_retictrees("/Users/anton/Desktop/data/droso_project/netsC1_score_nets.txt","/Users/anton/Desktop/data/droso_project/gene_trees_wboot_dna_mafft","clade1_score")
all_retictrees("/Users/anton/Desktop/data/droso_project/netsC2_score_nets.txt","/Users/anton/Desktop/data/droso_project/gene_trees_wboot_dna_mafft","clade2_score")
all_retictrees("/Users/anton/Desktop/data/droso_project/netsC3_score_nets.txt","/Users/anton/Desktop/data/droso_project/gene_trees_wboot_dna_mafft","clade3_score")
all_retictrees("/Users/anton/Desktop/data/droso_project/netsC4_score_nets.txt","/Users/anton/Desktop/data/droso_project/gene_trees_wboot_dna_mafft","clade4_score")
all_retictrees("/Users/anton/Desktop/data/droso_project/netsC5_score_nets.txt","/Users/anton/Desktop/data/droso_project/gene_trees_wboot_dna_mafft","clade5_score")
all_retictrees("/Users/anton/Desktop/data/droso_project/netsC6_score_nets.txt","/Users/anton/Desktop/data/droso_project/gene_trees_wboot_dna_mafft","clade6_score")
all_retictrees("/Users/anton/Desktop/data/droso_project/netsC7_score_nets.txt","/Users/anton/Desktop/data/droso_project/gene_trees_wboot_dna_mafft","clade7_score")
all_retictrees("/Users/anton/Desktop/data/droso_project/netsC8_score_nets.txt","/Users/anton/Desktop/data/droso_project/gene_trees_wboot_dna_mafft","clade8_score")
all_retictrees("/Users/anton/Desktop/data/droso_project/netsC9_score_nets.txt","/Users/anton/Desktop/data/droso_project/gene_trees_wboot_dna_mafft","clade9_score")