library(ape)
library(phytools)
library(phangorn)


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
            #tr$edge.length=NULL
            all_trees=c(all_trees,write.tree(tr))        
        }
       
    }
    write.table(all_trees,filename,quote = F,row.names = F, col.names=F)
}    

taxa_order=c("L_magnipectinata","L_stackelbergi","D_sturtevanti","D_neocordata","D_prosaltans","D_saltans","D_insularis","D_tropicalis","D_willistoni","D_equinoxialis","D_paulistorum","D_paulistorum_L12","D_guanche","D_subobscura","D_bifasciata","D_obscura","D_azteca","D_affinis","D_athabasca","D_lowei","D_miranda","D_persimilis","D_pseudoobscura","D_ercepeae","D_ananassae","D_parapallidosa","D_p_nigrens","D_pseudoananassae_pseudoananassae","D_malerkotliana_malerkotliana","D_m_pallens","D_bipectinata","D_parabipectinata","D_ficusphila","D_elegans","D_fuyamai","D_kurseongensis","D_carrolli","D_rhopaloa","D_takahashii","D_biarmipes","D_subpulchrella","D_suzukii","D_eugracilis","D_melanogaster","D_mauritiana","D_sechellia","D_simulans","D_erecta","D_orena","D_teissieri_273_3","D_yakuba","D_pectinifera","D_aurauria","D_triauraria","D_tani","D_rufa","D_asahinai","D_lacteicornis","D_kanapiae","D_kikkawai","D_leontia","D_punjabiensis","D_watanabei","D_bunnanda","D_serrata","D_truncata","D_birchii","D_mayri","D_jambulina","D_seguyi","D_vulcana","D_bakoue","D_nikananu","D_montium_STLow","D_burlai","D_montium_STUp","H_duncani","D_busckii","S_flava","S_montana","S_pallida","S_elmoi","D_primaeva","D_cyrtoloma","D_heteroneura","D_silvestris","D_differens","D_planitibia","D_paucipunta","D_prolacticillia","D_grimshawi","D_villosipedis","D_murphyi","D_sproati","D_limitata","D_ochracea","D_robusta","D_montana","D_virilis","D_lummei","D_americana","D_novamexicana","D_pseudotalamancana","D_bromeliae","D_nannoptera","D_pachea","D_acanthoptera","D_wassermani","D_hydei","D_repleta","D_seriema","D_navojoa","D_arizonae","D_mojavensis","D_repletoides","Z_ghesquierei","Z_inermis","Z_kolodkinae","Z_sepsoides_mayotte","Z_tsacasi","Z_tuberculatus","Z_ornatusmayotte","Z_africanus","Z_gabonicus","Z_indianus_16GNV01","Z_capensis","Z_davidi","Z_taronus","Z_vittiger","Z_lachaisei","Z_camerounensis","Z_nigranus","D_quadrilineata","D_funebris","D_guttifera","D_innubila","D_mush_saotome","D_pallidipennis","D_ornatifrons","D_subbadia","D_cardini","D_arawakana","D_dunni","D_pruinosa","D_immigrans","D_immigrans_kari17","D_neohypocausta","D_kepulauana","D_albomicans","D_nasuta","D_neonasuta","D_sulfurigaster_albostrigata","D_pulaua","D_sulfurigaster_bilimbata","D_sulfurigaster_sulfurigaster")

trees=read.tree("")
densiTree(trees[1:1000],scaleX = T,consensus=rev(taxa_order),alpha = 1,col=adjustcolor("black", alpha.f = 0.008),label.offset=0.01,cex=0.5,jitter = list(amount = 0.2, random=TRUE))