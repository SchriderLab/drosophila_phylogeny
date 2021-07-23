library("MCMCtreeR")
library("ape")
library("phytools")
library("PerformanceAnalytics")
library("dendextend")
library("BioGeoBEARS")
library("ggplot2")


#main tree plot
phy_mcmc=readMCMCtree("schemeA.tre")
phy=phy_mcmc$apePhy
tt=read.table("schemeA_mcmc.txt",header=T)

quartz(width=8,height=4)
quartz(width=8,height=4)
MCMC.tree.plot(rotate(keep.tip(phy,c(1:2,80,157:159)),9),analysis.type = "MCMCtree",MCMC.chain =tt[,1:6],plot.type = "distributions",density.col = adjustcolor( "navy", alpha.f = 0.5),density.border.col = "navy", lwd.bar = 3,scale.res = c("Period"), node.method = "bar",col.age = adjustcolor( "navy", alpha.f = 0.5), cex.labels = 0.5,cex.tips = 0.7,ladderize.tree = T,pos.age=-0.4,abs.age.lwd.ticks=0,relative.height=0.01,cex.age = 0.6,label.timescale.names = T,burn.in = 0,no.margin =T)
quartz.save("Fig1_droso.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12) 


quartz(width=5.5,height=12) 
MCMC.tree.plot(extract.clade(phy, 164),analysis.type = "MCMCtree",MCMC.chain=tt[,c(1,6:162)],plot.type = "distributions",density.col = adjustcolor( "navy", alpha.f = 0.5),density.border.col = "navy", lwd.bar = 3,scale.res = c("Period"), node.method = "bar",col.age = adjustcolor( "navy", alpha.f = 0.5), cex.labels = 0.5,cex.tips = 0.5,ladderize.tree = T,pos.age=-6.5,abs.age.lwd.ticks=0,relative.height=0.02,cex.age = 0.6,label.timescale.names = T,burn.in = 0,no.margin=T,label.offset = 0.5)
quartz.save("Fig1_drosoAB.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)


quartz(width=8.5,height=7.7)
MCMC.tree.plot(extract.clade(phy, 165),analysis.type = "MCMCtree",MCMC.chain=tt[,c(1,7:83)],plot.type = "distributions",density.col = adjustcolor( "navy", alpha.f = 0.5),density.border.col = "navy", lwd.bar = 3,scale.res = c("Period","Epoch"), node.method = "bar",col.age = adjustcolor( "navy", alpha.f = 0.5), cex.labels = 0.5,cex.tips = 0.5,ladderize.tree = T,pos.age=-11.5,abs.age.lwd.ticks=0,relative.height=0.05,cex.age = 0.6,label.timescale.names = T,burn.in = 0,no.margin=T)
quartz.save("Fig1_drosoB.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)


phy.small = extract.clade(phy, 164)
tt=tt[,c(1,6:159)]
MCMC.tree.plot(extract.clade(phy, 164),analysis.type = "MCMCtree",MCMC.chain =tt,plot.type = "distributions",density.col = adjustcolor( "navy", alpha.f = 0.5),density.border.col = "navy", lwd.bar = 3,scale.res = c("Period","Epoch"), node.method = "bar",col.age = adjustcolor( "navy", alpha.f = 0.5), no.margin = T, cex.labels = 0.01,cex.tips = 0.6,ladderize.tree = F,pos.age=-7,abs.age.lwd.ticks=0,relative.height=0.05,cex.age = 0.6)


#Dencitree plots
trim_trees=function(taxa,gene_trees,outg)
{
    gene_trees=gene_trees 
    trimmed_trees=c()
    
    for (tre in gene_trees)
    {
        if(any(taxa %in% tre$tip.label ) & outg %in% tre$tip.label)
        {
            tre=root(tre,outg)
            tre=drop.tip(tre,tre$tip.label[!tre$tip.label %in% taxa])
            tre$edge.length=NULL
            if(length(trimmed_trees)==0)
            {    
                trimmed_trees=tre
            }else{
                trimmed_trees=c(trimmed_trees,tre)
            }    
        }    
    }
    return(trimmed_trees)
}


#Co-plot tree topologies


outgr_remove1=c("Anopheles_gambiae","L_varia","C_costata","S_lebanonensis")
outgr_remove2=c("M_domestica","L_varia","L_trifolii","C_hians","E_gracilis","C_costata","S_lebanonensis")
outgr1="Anopheles_gambiae"
outgr2="M_domestica"


cophy=function(phy1,phy2,outgr_remove1,outgr1,outgr_remove2,outgr2,mleft,mright)
{
    phy1=root.phylo(phy1,outgroup=outgr1)
    phy2=root.phylo(phy2,outgroup=outgr2)
    phy1=drop.tip(phy1,outgr_remove1)
    phy2=drop.tip(phy2,outgr_remove2)
    phy1$edge.length=rep(1,length(phy1$edge.length))
    phy2$edge.length=rep(1,length(phy2$edge.length))
    phy1=chronos(phy1)
    phy2=chronos(phy2)

    #phy1$edge.length=NULL
    #phy2$edge.length=NULL
    tanglegram(phy1,phy2,sort=T,lab.cex=0.5,lwd=1.5,common_subtrees_color_branches = TRUE,margin_inner=6.5,rank_branches=T,main_left=mleft,main_right=mright)
    
}    
    

ml_tree=read.tree("iqtree.tre")
astral_tree=read.tree("astral.tre")
ml_prot_tree=read.tree("iqtree_protein.tre")
ml_tree_musca=read.tree("iqtree_differentout.tre")
astral_tree_treeshrink=read.tree("astral_treeshrink.tre")


quartz(width=12,height=11)
cophy(ml_tree,astral_tree,outgr_remove1,outgr1,outgr_remove1,outgr1,mleft="ML DNA Anopheles outgr",mright="ASTRAL DNA Anopheles outgr")
quartz.save("ML_ASTRAL.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
dev.off()

quartz(width=12,height=11)
cophy(ml_tree,ml_prot_tree,outgr_remove1,outgr1,outgr_remove1,outgr1,mleft="ML DNA Anopheles outgr",mright="ML AA Anopheles outgr")
quartz.save("ML_MLprot.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
dev.off()

quartz(width=12,height=11)
cophy(ml_tree,ml_tree_musca,outgr_remove1,outgr1,outgr_remove2,outgr2,mleft="ML DNA Anopheles outgr",mright="ML DNA Musca outgr")
quartz.save("ML_MLmusca.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
dev.off()

quartz(width=12,height=11)
cophy(ml_tree,astral_tree_treeshrink,outgr_remove1,outgr1,outgr_remove1,outgr1,mleft="ML DNA Anopheles outgr",mright="ASTRAL DNA TreeShrink Anopheles outgr")
quartz.save("ML_ASTRALtreeshrink.pdf", type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
dev.off()


#BEAST2/MCMCTREE age distributions


get_age_distr=function(trees,set_name)
{
    age_distr=c()
    for (i in 1:length(trees))
    {
        tre=trees[[i]]
        tre=drop.tip(tre,c("Oligophryne","Phytomyzites","Electrophortica","Scaptomyza_dominicana"))
        tab_age=prt(tre,printflag = F)
        rootage=tab_age[tab_age$label=="Anopheles_gambiae",c("ancestor","edge.length")]
        names(rootage)=c("node","time_bp")
        intern=tab_age[tab_age$node.type=="internal" ,c("node","time_bp")]
        tab_tree=rbind(rootage,intern)
        age_distr=rbind(age_distr,tab_tree)
        
    }
    age_distr$setname=set_name
    return(age_distr)   
}    

get_age_distr2=function(tree,set_name)
{

        
    tre=drop.tip(tree,c("Oligophryne","Phytomyzites","Electrophortica","Scaptomyza_dominicana"))
    tab_age=prt(tre,printflag = F)
    rootage=tab_age[tab_age$label=="Anopheles_gambiae",c("ancestor","edge.length")]
    names(rootage)=c("node","time_bp")
    intern=tab_age[tab_age$node.type=="internal" ,c("node","time_bp")]
    tab_tree=rbind(rootage,intern)
    tab_tree$setname=set_name
    return(tab_tree)   
}    


#Treannotator command
#treeannotator -burnin 0 -heights mean -lowMem rep_1.trees output.tre
beast_files=c("prior.trees.summary.tre",paste(paste("rep",1:10,sep="_"),".trees.summary.tre",sep=""))
all_ages=c()
for (f in beast_files)
{
    tre=read.nexus(f)
    age_tab=get_age_distr2(tre,f)
    age_tab$method="BEAST"
    all_ages=rbind(all_ages,age_tab)
}    
mcmctree_files=c("schemeA.tre","schemeB.tre","schemeC.tre","schemeD.tre","schemeRusso.tre")

for (f in mcmctree_files)
{
    phy_mcmc=readMCMCtree(f)
    phy=phy_mcmc$apePhy
    age_tab=get_age_distr2(phy,f)
    age_tab$method="MCMCTREE"
    all_ages=rbind(all_ages,age_tab)
}



ggplot(all_ages,aes(x=as.character(node),y=(time_bp)))+
geom_point(aes(color=setname))+
theme(axis.text.x = element_text(size=5.5, angle=45))+
scale_y_continuous(breaks=seq(0,301,25))+
xlab("Node number")+
ylab("Age in MYA")+
scale_colour_manual(values = c("black",heat.colors(10),viridis(10)[3:7]),labels = c("BEAST FBD Prior",paste("BEAST",1:10),"MCMCTREE schemeA","MCMCTREE schemeB","MCMCTREE schemeC","MCMCTREE schemeD","MCMCTREE schemeRusso"))

# Comparisons of MCMCTREE with diffrent loci  

all_ages_genes=c()
mcmctree_files_genes=c("schemeA.tre","schemeA_10loci.tre","schemeA_100loci.tre")
for (f in mcmctree_files_genes)
{
    phy_mcmc=readMCMCtree(f)
    phy=phy_mcmc$apePhy
    age_tab=get_age_distr2(phy,f)
    age_tab$method="MCMCTREE"
    all_ages_genes=rbind(all_ages_genes,age_tab)
}

ggplot(all_ages_genes,aes(x=reorder(node, time_bp),y=time_bp,group=setname))+
geom_point(aes(color=setname))+
theme_light()+
theme(axis.text.x = element_text(size=5.5, angle=45))+
coord_trans(y="log")+
scale_y_continuous(breaks=c(0,1,2,3,4,5,10,25,50,100,150,200,250))+
xlab("Node number")+
ylab("Age")+
geom_line(aes(color=setname))+
scale_colour_manual(name="",values = c("black","red","blue"),labels = c("MCMCTREE schemeA 100 loci","MCMCTREE schemeA 10 loci","MCMCTREE schemeA 1000 loci"))+
theme(legend.position="top", legend.box = "horizontal")
ggsave("mtcarsee.pdf",width = 491, height = 266,units="mm")

#Plot BEAST2 age distributions 


phy_beast=read.nexus("rep_5.trees.summary.tre")
phy_beast=drop.tip(phy,c("Oligophryne","Phytomyzites","Electrophortica","Scaptomyza_dominicana"))
trees=read.nexus("rep_5.trees.top5")
beast_ages=get_age_distr(trees,"BEAST2")
node_posteriors = vector(mode = "list", length = Nnode(phy_beast))
names(node_posteriors) = unique(beast_ages$node)
pos=1
for (n in unique(beast_ages$node))
{
    node_posteriors[[pos]]=beast_ages[beast_ages$node==n,"time_bp"]
    pos=pos+1
}



mcmc_df = data.frame(t(data.frame(matrix(unlist(node_posteriors), nrow=length(node_posteriors), byrow=TRUE))))
names(mcmc_df)=paste("t_n",unique(beast_ages$node),sep="")
mcmc_df=cbind(Gen=1:nrow(mcmc_df),mcmc_df)


quartz(width=14,height=11)
MCMC.tree.plot(phy_beast,analysis.type = "MCMCtree",MCMC.chain =mcmc_df,plot.type = "cladogram",density.col = adjustcolor( "red", alpha.f = 0.5),density.border.col = "red", lwd.bar = 3,scale.res = c("Period"), node.method = "bar",col.age = adjustcolor( "red", alpha.f = 0.5), cex.labels = 0.5,cex.tips = 0.5,ladderize.tree = T,pos.age=-0.4,abs.age.lwd.ticks=0,relative.height=0.01,cex.age = 0.6,label.timescale.names = T,burn.in = 0,no.margin =T)
quartz.save("beast_ci.psd", type = "psd",antialias=F,bg="white",dpi=400,pointsize=12) 


quartz(width=14,height=11)
MCMC.tree.plot(phy_beast,analysis.type = "MCMCtree",MCMC.chain =tt[,1:162],plot.type = "cladogram",density.col = adjustcolor( "navy", alpha.f = 0.5),density.border.col = "navy", lwd.bar = 3,add.time.scale=F, node.method = "bar",col.age = adjustcolor( "navy", alpha.f = 0.5), cex.labels = 0.5,cex.tips = 0.0000001,ladderize.tree = T,pos.age=-0.4,abs.age.lwd.ticks=0,relative.height=0.01,cex.age = 0.6,label.timescale.names = T,burn.in = 0,no.margin =T,edge.width=0.1)
quartz.save("mcmctree_ci.psd", type = "psd",antialias=F,bg="white",dpi=400,pointsize=12)


#ASTRAL-ML Branch length comparison
ml_tree=read.tree("iqtree.tre")
astral_tree=read.tree("astral.tre")
ml_tree=root.phylo(ml_tree,outgroup="Anopheles_gambiae",resolve.root = T)
astral_tree=root.phylo(astral_tree,outgroup="Anopheles_gambiae",resolve.root = T)
ml_tree_sub=drop.tip(ml_tree,c("D_vulcana","D_simulans","D_villosipedis","Anopheles_gambiae","L_varia"))
astral_tree_sub=drop.tip(astral_tree,c("D_vulcana","D_simulans","D_villosipedis","Anopheles_gambiae","L_varia"))




<init id="SBI" spec="starbeast2.StarBeastInitializer" birthRate="@netDiversificationRate.t:Species" estimate="false" speciesTree="@Tree.t:Species" newick="((((C_costata1:48.321835,S_lebanonensis1:48.321835):5.093982,((((((((((((D_acanthoptera1:3.845434,D_wassermani1:3.845434):0.399792,D_pachea1:4.245226):1.403876,D_nannoptera1:5.649102):8.824571,D_bromeliae1:14.473673):5.300776,(((((D_arizonae1:1.960916,D_mojavensis1:1.960916):1.895776,D_navojoa1:3.856691):4.415384,D_seriema1:8.272075):2.180421,D_repleta1:10.452496):1.498834,D_hydei1:11.95133):7.823118):1.527575,D_pseudotalamancana1:21.302023):3.034649,((((D_americana1:1.637444,D_novamexicana1:1.637444):0.747678,D_lummei1:2.385121):0.836676,D_virilis1:3.221798):2.114214,D_montana1:5.336012):19.00066):2.410028,D_robusta1:26.7467):3.324544,((((D_cyrtoloma1:2.421318,((D_differens1:1.232136,D_planitibia1:1.232136):0.288142,(D_heteroneura1:0.845331,D_silvestris1:0.845331):0.674947):0.90104):1.035397,((D_grimshawi1:2.489819,(((D_limitata1:1.175918,D_ochracea1:1.175918):0.516549,D_villosipedis1:1.692467):0.120043,(D_murphyi1:1.322782,D_sproati1:1.322782):0.489728):0.677309):0.309197,(D_paucipunta1:1.746868,D_prolacticillia1:1.746868):1.052148):0.657698):1.94271,D_primaeva1:5.399425):15.93909,((S_elmoi1:3.512012,S_pallida1:3.512012):11.724135,(S_flava1:2.687229,S_montana1:2.687229):12.548918):6.102368):8.732729):4.043368,((((((((((D_albomicans1:1.872232,D_nasuta1:1.872232):0.490283,D_kepulauana1:2.362516):0.609009,((D_neonasuta1:1.592289,D_sulfurigaster_albostrigata1:1.592289):0.434222,(D_pulaua1:1.526559,(D_sulfurigaster_bilimbata1:1.300024,D_sulfurigaster_sulfurigaster1:1.300024):0.226535):0.499952):0.945013):7.98714,D_neohypocausta1:10.958664):2.007233,(D_immigrans1:1.019262,D_immigrans_kari171:1.019262):11.946635):6.579911,D_pruinosa1:19.545808):5.190667,(((((D_arawakana1:3.70265,D_dunni1:3.70265):3.767982,D_cardini1:7.470632):5.854747,(D_ornatifrons1:2.296925,D_subbadia1:2.296925):11.028454):2.813074,D_pallidipennis1:16.138453):4.231955,(D_funebris1:17.055686,(D_guttifera1:13.462469,(D_innubila1:9.821877,D_mush_saotome1:9.821877):3.640592):3.593216):3.314722):4.366067):3.241803,D_quadrilineata1:27.978278):2.692728,((((Z_africanus1:3.549113,(Z_gabonicus1:1.503132,Z_indianus_16GNV011:1.503132):2.045981):2.269993,((((Z_camerounensis1:1.687499,Z_nigranus1:1.687499):0.59985,Z_lachaisei1:2.287349):0.235784,Z_vittiger1:2.523133):1.317864,(Z_capensis1:3.037572,(Z_davidi1:2.487502,Z_taronus1:2.487502):0.55007):0.803425):1.978109):1.043449,Z_ornatusmayotte1:6.862555):7.337531,(Z_ghesquierei1:12.16053,(Z_inermis1:7.586313,(Z_kolodkinae1:4.266969,(Z_sepsoides_mayotte1:1.841424,(Z_tsacasi1:1.318827,Z_tsacasi_21:1.318827):0.522597):2.425545):3.319344):4.574217):2.039556):16.47092):1.702121,D_repletoides1:32.373127):1.741485):4.957086,D_busckii1:39.071698):7.764337,(((((((D_affinis1:2.664562,D_athabasca1:2.664562):0.980901,D_azteca1:3.645463):3.619446,(D_lowei1:4.799589,(D_miranda1:3.274875,(D_persimilis1:2.370625,D_pseudoobscura1:2.370625):0.90425):1.524713):2.465321):2.839571,((D_bifasciata1:7.234309,D_obscura1:7.234309):1.239389,(D_guanche1:2.719657,D_subobscura1:2.719657):5.754041):1.630782):22.678244,((((D_ananassae1:3.299068,D_parapallidosa1:3.299068):5.434774,(((D_bipectinata1:2.25316,D_parabipectinata1:2.25316):0.380818,(D_m_pallens1:2.160349,D_malerkotliana_malerkotliana1:2.160349):0.473628):1.608527,(D_p_nigrens1:1.84064,D_pseudoananassae_pseudoananassae1:1.84064):2.401865):4.491337):4.333355,D_ercepeae1:13.067197):12.620076,(((((((D_asahinai1:1.708814,D_lacteicornis1:1.708814):0.535947,D_rufa1:2.244761):1.181091,D_tani1:3.425852):3.217834,(D_aurauria1:1.976286,D_triauraria1:1.976286):4.667399):3.672303,D_pectinifera1:10.315988):2.371841,(((((((D_bakoue1:4.941256,(((D_burlai1:1.970973,D_sp_aff_chauvacae1:1.970973):0.54447,D_bocqueti1:2.515443):1.633026,D_nikananu1:4.148469):0.792787):0.374147,(D_seguyi1:4.876809,D_vulcana1:4.876809):0.438595):0.645877,D_jambulina1:5.96128):1.91916,(((D_birchii1:3.165172,D_mayri1:3.165172):2.264285,D_truncata1:5.429457):1.756392,(D_bunnanda1:5.257097,D_serrata1:5.257097):1.928752):0.694591):0.383898,(D_punjabiensis1:3.623767,D_watanabei1:3.623767):4.640571):0.883101,(D_kikkawai1:2.225816,D_leontia1:2.225816):6.921623):2.19085,D_kanapiae1:11.338289):1.34954):9.466695,(((((D_biarmipes1:6.771533,(D_subpulchrella1:3.960115,D_suzukii1:3.960115):2.811419):3.862028,D_takahashii1:10.633561):2.553753,((((D_erecta1:3.30116,D_orena1:3.30116):2.167638,(D_teissieri_273_31:3.479504,D_yakuba1:3.479504):1.989295):1.256779,(((D_mauritiana1:1.945627,D_simulans1:1.945627):0.218613,D_sechellia1:2.16424):1.457463,D_melanogaster1:3.621703):3.103875):4.996078,D_eugracilis1:11.721656):1.465659):1.814617,((((D_carrolli1:2.044694,D_rhopaloa1:2.044694):1.017073,D_kurseongensis1:3.061766):0.316665,D_fuyamai1:3.378431):7.135146,D_elegans1:10.513576):4.488354):1.621014,D_ficusphila1:16.622945):5.531579):3.532749):7.095452):7.52439,((((((D_equinoxialis1:2.600592,(D_paulistorum1:1.232377,D_paulistorum_L121:1.232377):1.368214):1.663272,D_willistoni1:4.263863):0.839402,D_tropicalis1:5.103265):1.518778,D_insularis1:6.622043):10.732717,((D_neocordata1:8.227214,(D_prosaltans1:1.568998,D_saltans1:1.568998):6.658216):2.023527,D_sturtevanti1:10.250741):7.10402):13.919696,(L_magnipectinata1:8.746987,L_stackelbergi1:8.746987):22.52747):9.032659):3.161823,H_duncani1:43.468939):3.367097):6.579782):9.776416,L_varia1:63.192233):140.582809,Anopheles_gambiae1:203.775042);">

Reheight.t:Species
SubtreeSlide.t:Species
Narrow.t:Species
Wide.t:Species
WilsonBalding.t:Species




#FBD Uncertanties

<operator spec='SampledNodeDateRandomWalker' windowSize="1"  tree="@Tree.t:rep_1" weight="10">

        <taxonset spec="TaxonSet">
            <taxon id="Oligophryne" spec="Taxon"/>
            <taxon id="Phytomyzites" spec="Taxon"/>
            <taxon id="Electrophortica" spec="Taxon"/>
            <taxon id="Scaptomyza_dominicana" spec="Taxon"/>
        </taxonset>

        <samplingDates id="samplingDate1" spec="beast.evolution.tree.SamplingDate" taxon="Oligophryne" upper="201.6" lower="189.6"/>
        <samplingDates id="samplingDate2" spec="beast.evolution.tree.SamplingDate" taxon="Phytomyzites" upper="66." lower="61.7"/>
        <samplingDates id="samplingDate3" spec="beast.evolution.tree.SamplingDate" taxon="Electrophortica" upper="48." lower="34."/>
        <samplingDates id="samplingDate4" spec="beast.evolution.tree.SamplingDate" taxon="Scaptomyza_dominicana" upper="20.43" lower="13.65"/>
</operator>