library('paleotree')
library('seqinr')
### parameters 
options(scipen=999)
args = commandArgs(trailingOnly=TRUE)

mrbayesnexus_gen=function(fasta,tree,outgroup,nsim,filename="mrbayes.out") 
{
    #tree: path to ROOTED tree file in newick format 
    #outgroup: outgroup taxon 
    #nsim: total number of MCMC samples   
    read_tree=read.tree(tree)
    #Create dummy data
    fasta_f=read.fasta(fasta)
    write.nexus.data(fasta_f,"mrbayes_config.nex",interleaved = F)
    #Create constrains
    constr_partial=read_tree$tip.label[!read_tree$tip.label %in% outgroup]
    constr_partial=paste("\tconstraint ingroup partial =",paste(constr_partial,collapse=" "),":",outgroup,";")
    constr_hard=createMrBayesConstraints(read_tree, partial = F,includeIngroupConstraint = T)
    constr_hard=paste("\t",constr_hard,sep="")
    #Nodes to calibrate
    nodes_calibrate=c()
    for (n in 1:read_tree$Nnode)
    {
        if (read_tree$node.label[1]=="")
        {
            stop("Error: No information for treeagepr")
        } else {
            if (n == 1)
            {
                tr_age_pr=paste("\tprset treeagepr =",read_tree$node.label[1],";\n")
                nodes_calibrate=c(nodes_calibrate,tr_age_pr)
            } else if (read_tree$node.label[n]!="") {
                
                n_age_pr=paste("\tcalibrate node",n-1," = ",read_tree$node.label[n],";\n",sep="")
                nodes_calibrate=c(nodes_calibrate,n_age_pr)
            }   
        }   
    }    
     
    bd_params=
    "\tprset brlenspr=clock:birthdeath;
     \tprset speciationpr = fixed(0.03);
     \tprset extinctionpr = fixed(0.4);
     \tprset sampleprob = 1;
     \tprset nodeagepr = calibrated;"

    mcmc_params=paste(
    "\tmcmc ngen=",nsim*10000," samplefreq=10000 printfreq=10000 nruns=1 nchains=1 data=no savetrees=yes;\n\tsumt burninfrac=0.5 output=",filename,"savebrparams=yes;")
    
     write("begin mrbayes;","mrbayes_config.nex",append=T)
     write(constr_partial,"mrbayes_config.nex",append=T)
     write(constr_hard,"mrbayes_config.nex",append=T)
     write(bd_params,"mrbayes_config.nex",append=T)
     write(nodes_calibrate,"mrbayes_config.nex",append=T)
     write(mcmc_params,"mrbayes_config.nex",append=T)
     write("end;","mrbayes_config.nex",append=T)    
}    