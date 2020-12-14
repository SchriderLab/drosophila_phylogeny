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

all_triplets(extract.clade(tt,168)$tip.label,"gene_trees_wboot_dna_mafft",paste("C",1,"_node",168,"_quibl",sep=""))
all_triplets(extract.clade(tt,185)$tip.label,"gene_trees_wboot_dna_mafft",paste("C",2,"_node",185,"_quibl",sep=""))
all_triplets(extract.clade(tt,204)$tip.label,"gene_trees_wboot_dna_mafft",paste("C",3,"_node",204,"_quibl",sep=""))
all_triplets(extract.clade(tt,226)$tip.label,"gene_trees_wboot_dna_mafft",paste("C",4,"_node",226,"_quibl",sep=""))
all_triplets(extract.clade(tt,245)$tip.label,"gene_trees_wboot_dna_mafft",paste("C",5,"_node",245,"_quibl",sep=""))
all_triplets(extract.clade(tt,256)$tip.label,"gene_trees_wboot_dna_mafft",paste("C",6,"_node",256,"_quibl",sep=""))
all_triplets(extract.clade(tt,265)$tip.label,"gene_trees_wboot_dna_mafft",paste("C",7,"_node",265,"_quibl",sep=""))
all_triplets(extract.clade(tt,289)$tip.label,"gene_trees_wboot_dna_mafft",paste("C",8,"_node",289,"_quibl",sep=""))
all_triplets(extract.clade(tt,307)$tip.label,"gene_trees_wboot_dna_mafft",paste("C",9,"_node",307,"_quibl",sep=""))



