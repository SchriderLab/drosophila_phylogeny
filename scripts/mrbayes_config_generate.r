library('paleotree')
library('seqinr')
### parameters 
options(scipen=999)
args = commandArgs(trailingOnly=TRUE)

mrbayesnexus_gen=function(fasta,tree,partition,outgroup,nsim,filename="mrbayes.out") 
{
    #tree: path to ROOTED tree file in newick format 
    #outgroup: outgroup taxon 
    #nsim: total number of MCMC samples   
    read_tree=read.tree(tree)
    #Create data NEXUS block
    fasta_f=read.fasta(fasta)
    write.nexus.data(fasta_f,"mrbayes_config.nex",interleaved = F)
    #Create partitions
    parts=read.csv(partition,header=F)
    parts=cbind(V0="\tcharset",V1=parts$V1,"=",V2=parts$V2)
    #Create constrains
    constr_partial=read_tree$tip.label[!read_tree$tip.label %in% outgroup]
    constr_partial=paste("\tconstraint ingroup partial =",paste(constr_partial,collapse=" "),":",outgroup,";")
    constr_hard=createMrBayesConstraints(read_tree, partial = F,includeIngroupConstraint = T)
    constr_hard=paste("\t",constr_hard,sep="")
    #Nodes to calibrate
    tr_age_pr=paste("\tprset treeagepr = ;\n")
                
            
    clock_params=
    "\tlset applyto=(all) nst=6 rates = gamma;
     \tprset applyto = (all) ratepr = variable;
     \tunlink statefreq = (all) revmat = (all) shape = (all); 
     \tprset brlenspr=clock:birthdeath;
     \tprset speciationpr = uniform(0,1);
     \tprset extinctionpr = beta(1,1);
     \tprset sampleprob = 1;
     \tprset clockvarpr = igr;
     \tprset igrvarpr = exp(10);
     \tprset nodeagepr = calibrated;"

    mcmc_params=paste(
    "\tmcmc ngen=",nsim*10000," samplefreq=100 printfreq=100 nruns=1 nchains=1 data=yes savetrees=yes;\n\tsumt burninfrac=0.5 output=",filename,"savebrparams=yes;")
    
     write("begin mrbayes;","mrbayes_config.nex",append=T)
     write.table(parts,"mrbayes_config.nex",quote = F,row.names = F,col.names = F,append=T)
     write(constr_partial,"mrbayes_config.nex",append=T)
     write(constr_hard,"mrbayes_config.nex",append=T)
     write(clock_params,"mrbayes_config.nex",append=T)
     write(tr_age_pr,"mrbayes_config.nex",append=T)
     write(mcmc_params,"mrbayes_config.nex",append=T)
     write("end;","mrbayes_config.nex",append=T)    
}    