#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library("ape")
trees=read.tree(args[1])
index=1
for (tr in trees)
{
    tr$node.label=NULL
    trees[[index]]=tr
    index=index+1
}    
write.tree(trees,"trees_woboot")