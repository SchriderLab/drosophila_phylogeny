library("ape")
library("phytools")
library("phangorn")

#Gene trees densitree 
genefordensitree=function(alltree,filename)
{
    all_trees=c()
    
    outgr=c("S_lebanonensis", "C_costata","L_varia","Anopheles_gambiae")
    for (tr in alltree)
    {
        
        if (any(outgr %in% tr$tip.label) & length(tr$tip.label) >= 6 )
        {
            
            for (taxon in outgr)
            {
                if (taxon %in% tr$tip.label)
                {
                    tr=root(tr,taxon,resolve.root = T)
                    tr=drop.tip(tr,taxon)
                }      
             
            } 
            tr$node.label=NULL
            tr$edge.length=NULL
            #tr$edge.length=NULL
            all_trees=c(all_trees,write.tree(tr))        
        }
       
    }
    write.table(all_trees,filename,quote = F,row.names = F, col.names=F)
}   

occup=function(alltree)
{
    allsps=c()
    for (tr in alltree)
    {
        allsps=c(allsps,tr$tip.label)
    }
    return(allsps)
    
}    
    

cladefordensitree=function(alltree,sps)
{
    all_trees=c()
    
    for (tr in alltree)
    {
        
        if (sum(sps %in% tr$tip.label) >= 6 )
        {
            tr_clade=keep.tip(tr,tr$tip.label[tr$tip.label %in% sps])
        }
        all_trees=c(all_trees,write.tree(tr_clade))
        
    }
    return(read.tree(text=all_trees))
}    
    



trees=read.tree("gene_trees_wboot_dna_mafft")
genefordensitree(trees,"gene_trees_fordensi")
gg=occup(trees)
barplot(table(gg),las=2,cex.names=0.5)
trees=read.tree("gene_trees_fordensi")
densiTree(trees[1:1000],scaleX = T,consensus=rev(taxa_order),alpha = 1,col=adjustcolor("black", alpha.f = 0.008),label.offset=0.01,cex=0.5,jitter = list(amount = 0.2, random=TRUE))


C1=c('H_duncani','L_stackelbergi','L_magnipectinata','D_sturtevanti','D_neocordata','D_saltans','D_prosaltans','D_insularis','D_tropicalis','D_willistoni','D_equinoxialis','D_paulistorum_L12','D_paulistorum')
C2=c('D_subobscura','D_guanche','D_obscura','D_bifasciata','D_azteca','D_athabasca','D_affinis','D_lowei','D_miranda','D_pseudoobscura','D_persimilis')
C3=c('D_ercepeae','D_parapallidosa','D_ananassae','D_pseudoananassae_pseudoananassae','D_p_nigrens','D_malerkotliana_malerkotliana','D_m_pallens','D_parabipectinata','D_bipectinata')
C4=c('D_ficusphila','D_elegans','D_fuyamai','D_kurseongensis','D_rhopaloa','D_carrolli','D_takahashii','D_biarmipes','D_suzukii','D_subpulchrella','D_eugracilis','D_melanogaster','D_sechellia','D_simulans','D_mauritiana','D_yakuba','D_teissieri_273_3','D_orena','D_erecta')
C5=c('D_pectinifera','D_triauraria','D_aurauria','D_tani','D_rufa','D_lacteicornis','D_asahinai','D_kanapiae','D_leontia','D_kikkawai','D_watanabei','D_punjabiensis','D_serrata','D_bunnanda','D_truncata','D_mayri','D_birchii','D_jambulina','D_vulcana','D_seguyi','D_bakoue','D_nikananu','D_montium_STLow','D_montium_STUp','D_burlai')
'D_busckii'
C6=c('S_montana','S_flava','S_pallida','S_elmoi','D_primaeva','D_cyrtoloma','D_silvestris','D_heteroneura','D_planitibia','D_differens','D_prolacticillia','D_paucipunta','D_grimshawi','D_sproati','D_murphyi','D_villosipedis','D_ochracea','D_limitata')
C7=c('D_robusta','D_montana','D_virilis','D_lummei','D_novamexicana','D_americana','D_pseudotalamancana','D_bromeliae','D_nannoptera','D_pachea','D_wassermani','D_acanthoptera','D_hydei','D_repleta','D_seriema','D_navojoa','D_mojavensis','D_arizonae')
'D_repletoides'
C8=c('Z_ghesquierei','Z_inermis','Z_kolodkinae','Z_sepsoides_mayotte','Z_tuberculatus','Z_tsacasi','Z_ornatusmayotte','Z_africanus','Z_indianus_16GNV01','Z_gabonicus','Z_capensis','Z_taronus','Z_davidi','Z_vittiger','Z_lachaisei','Z_nigranus','Z_camerounensis')
C9=c('D_quadrilineata','D_funebris','D_guttifera','D_mush_saotome','D_innubila','D_pallidipennis','D_subbadia','D_ornatifrons','D_cardini','D_dunni','D_arawakana','D_pruinosa','D_immigrans_kari17','D_immigrans','D_neohypocausta','D_kepulauana','D_nasuta','D_albomicans','D_sulfurigaster_albostrigata','D_neonasuta','D_pulaua','D_sulfurigaster_sulfurigaster','D_sulfurigaster_bilimbata')

0.53
quartz(width=4.8,height=0.53*12)
sp_space=list(C1=C1,C2=C2,C3=C3,C4=C4,C5=C5,C6=C6,C7=C7,C8=C8,C9=C9)
i=1
for (cl in sp_space)
{
    subtr=cladefordensitree(trees,cl)
    quartz(width=4.8,height=0.53*length(cl))
    densiTree(subtr,scaleX = T,consensus=rev(cl),alpha = 1,col=adjustcolor("black", alpha.f = 0.008),label.offset=0.01,cex=0.5,jitter = list(amount = 0, random=TRUE))
    quartz.save(paste("C",i,"densi",".pdf",sep=""), type = "pdf",antialias=F,bg="white",dpi=400,pointsize=12)
    quartz.save(paste("C",i,"densi",".png",sep=""), type = "png",antialias=F,bg="white",dpi=400,pointsize=12)
    dev.off()
    i=i+1
}    

