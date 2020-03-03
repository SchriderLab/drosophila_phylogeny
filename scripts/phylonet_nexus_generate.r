library("ape")

phylonet_full_topo=function(sp_tree,gene_trees,n_retic)
{
    sp_phy=read.tree(sp_tree)
    sp_phy$node.label=NULL
    phy=read.tree(gene_trees)
    rooted_phy=c()
    out_v=c("Anopheles_gambiae","L_varia","S_lebanonensis","C_costata")
    for (gt in phy)
    {
        
        if (any(out_v %in% gt$tip.label))
        {
            out_sp=out_v[which(out_v %in% gt$tip.label)[1]]
            gt_r=root(gt,out_sp,resolve.root=T)
            gt_r$node.label[gt_r$node.label=="Root"]=100
            gt_r$node.label[gt_r$node.label==""]=1
            gt_r$node.label=as.numeric(gt_r$node.label)/100
            rooted_phy=c(rooted_phy,write.tree(gt_r))
            
        } 
        
       
    } 
    phy=rooted_phy
    for (i in 1:n_retic)
    {    
        f_n=paste("phylonet_genes_",i,"ret.nex",sep="")
        write("#NEXUS\n\nBEGIN TREES;",f_n)
        write(paste("Tree fixtr = ",write.tree(sp_phy)),f_n,append=T)
        d=data.frame(rep("Tree",length(phy)),paste("gt",1:length(phy),"=",sep=""),phy)
        write.table(d,f_n,quote = F,row.names = F, col.names=F,append=T)
        write("END;\n\nBEGIN PHYLONET;",f_n,append=T)
        write(paste("InferNetwork_MPL (all)",i,"-s fixtr -fs -di -pl 15 -x 1;","\nEND;"),f_n,append=T)
    }    
}    

phylonet_full_topo("/Users/Anton/Downloads/BUSCO50_dna_pasta_nopart_iqtree_root.tre","/Users/Anton/Downloads/BUSCO50_dna_pasta_iqtree_all",1)



busckii=c('C_costata','S_lebanonensis','D_acanthoptera','D_wassermani','D_pachea','D_nannoptera','D_bromeliae','D_arizonae','D_mojavensis','D_navojoa','D_seriema','D_repleta','D_hydei','D_pseudotalamancana','D_americana','D_novamexicana','D_lummei','D_virilis','D_montana','D_robusta','D_cyrtoloma','D_differens','D_planitibia','D_heteroneura','D_silvestris','D_grimshawi','D_limitata','D_ochracea','D_villosipedis','D_murphyi','D_sproati','D_paucipunta','D_prolacticillia','D_primaeva','S_elmoi','S_pallida','S_flava','S_montana','D_albomicans','D_nasuta','D_kepulauana','D_neonasuta','D_sulfurigaster_albostrigata','D_pulaua','D_sulfurigaster_bilimbata','D_sulfurigaster_sulfurigaster','D_neohypocausta','D_immigrans','D_immigrans_kari17','D_pruinosa','D_arawakana','D_dunni','D_cardini','D_ornatifrons','D_subbadia','D_pallidipennis','D_funebris','D_guttifera','D_innubila','D_mush_saotome','D_quadrilineata','Z_africanus','Z_gabonicus','Z_indianus_16GNV01','Z_camerounensis','Z_nigranus','Z_lachaisei','Z_vittiger','Z_capensis','Z_davidi','Z_taronus','Z_ornatusmayotte','Z_ghesquierei','Z_inermis','Z_kolodkinae','Z_sepsoides_mayotte','Z_tsacasi','Z_tuberculatus','D_repletoides','D_busckii','L_varia','Anopheles_gambiae')

duncani=c('D_affinis','D_athabasca','D_azteca','D_lowei','D_miranda','D_persimilis','D_pseudoobscura','D_bifasciata','D_obscura','D_guanche','D_subobscura','D_ananassae','D_parapallidosa','D_bipectinata','D_parabipectinata','D_m_pallens','D_malerkotliana_malerkotliana','D_p_nigrens','D_pseudoananassae_pseudoananassae','D_ercepeae','D_asahinai','D_lacteicornis','D_rufa','D_tani','D_aurauria','D_triauraria','D_pectinifera','D_bakoue','D_burlai','D_montium_STUp','D_montium_STLow','D_nikananu','D_seguyi','D_vulcana','D_jambulina','D_birchii','D_mayri','D_truncata','D_bunnanda','D_serrata','D_punjabiensis','D_watanabei','D_kikkawai','D_leontia','D_kanapiae','D_biarmipes','D_subpulchrella','D_suzukii','D_takahashii','D_erecta','D_orena','D_teissieri_273_3','D_yakuba','D_mauritiana','D_simulans','D_sechellia','D_melanogaster','D_eugracilis','D_carrolli','D_rhopaloa','D_kurseongensis','D_fuyamai','D_elegans','D_ficusphila','D_equinoxialis','D_paulistorum','D_paulistorum_L12','D_willistoni','D_tropicalis','D_insularis','D_neocordata','D_prosaltans','D_saltans','D_sturtevanti','L_magnipectinata','L_stackelbergi','H_duncani','L_varia','Anopheles_gambiae','C_costata','S_lebanonensis')

A=c('D_equinoxialis','D_paulistorum','D_paulistorum_L12','D_willistoni','D_tropicalis','D_insularis','D_neocordata','D_prosaltans','D_saltans','D_sturtevanti','L_magnipectinata','L_stackelbergi')
B=c('D_affinis','D_athabasca','D_azteca','D_lowei','D_miranda','D_persimilis','D_pseudoobscura','D_bifasciata','D_obscura','D_guanche','D_subobscura')
C=c('D_ananassae','D_parapallidosa','D_bipectinata','D_parabipectinata','D_m_pallens','D_malerkotliana_malerkotliana','D_p_nigrens','D_pseudoananassae_pseudoananassae','D_ercepeae')
D=c('D_biarmipes','D_subpulchrella','D_suzukii','D_takahashii','D_erecta','D_orena','D_teissieri_273_3','D_yakuba','D_mauritiana','D_simulans','D_sechellia','D_melanogaster','D_eugracilis','D_carrolli','D_rhopaloa','D_kurseongensis','D_fuyamai','D_elegans','D_ficusphila')
E=c('D_asahinai','D_lacteicornis','D_rufa','D_tani','D_aurauria','D_triauraria','D_pectinifera','D_bakoue','D_burlai','D_montium_STUp','D_montium_STLow','D_nikananu','D_seguyi','D_vulcana','D_jambulina','D_birchii','D_mayri','D_truncata','D_bunnanda','D_serrata','D_punjabiensis','D_watanabei','D_kikkawai','D_leontia','D_kanapiae')
F=c('D_cyrtoloma','D_differens','D_planitibia','D_heteroneura','D_silvestris','D_grimshawi','D_limitata','D_ochracea','D_villosipedis','D_murphyi','D_sproati','D_paucipunta','D_prolacticillia','D_primaeva','S_elmoi','S_pallida','S_flava','S_montana')
G=c('D_acanthoptera','D_wassermani','D_pachea','D_nannoptera','D_bromeliae','D_arizonae','D_mojavensis','D_navojoa','D_seriema','D_repleta','D_hydei','D_pseudotalamancana','D_americana','D_novamexicana','D_lummei','D_virilis','D_montana','D_robusta')
H=c('Z_africanus','Z_gabonicus','Z_indianus_16GNV01','Z_camerounensis','Z_nigranus','Z_lachaisei','Z_vittiger','Z_capensis','Z_davidi','Z_taronus','Z_ornatusmayotte','Z_ghesquierei','Z_inermis','Z_kolodkinae','Z_sepsoides_mayotte','Z_tsacasi','Z_tuberculatus')
J=c('D_arawakana','D_dunni','D_cardini','D_ornatifrons','D_subbadia','D_pallidipennis','D_funebris','D_guttifera','D_innubila','D_mush_saotome')
K=c('D_albomicans','D_nasuta','D_kepulauana','D_neonasuta','D_sulfurigaster_albostrigata','D_pulaua','D_sulfurigaster_bilimbata','D_sulfurigaster_sulfurigaster','D_neohypocausta','D_immigrans','D_immigrans_kari17','D_pruinosa')



phylonet_chunk_topo=function(sp_tree,gene_trees,sp_list,n_retic)
{
    sp_phy=read.tree(sp_tree)
    sp_phy$node.label=NULL
    out_v=c("Anopheles_gambiae")
    sp_phy_sub=keep.tip(sp_phy,c(sp_list[sp_list %in% sp_phy$tip.label],out_v))
    phy=read.tree(gene_trees)
    rooted_phy=c()
    for (gt in phy)
    {
        if (any(sp_list %in% gt$tip.label) & any(out_v %in% gt$tip.label))
        {
            if (is.monophyletic(gt,sp_list[sp_list %in% gt$tip.label])) 
            {
                
                out_sp=out_v[which(out_v %in% gt$tip.label)[1]]
                gt_r=root(gt,out_sp,resolve.root=T)
                gt_sub=keep.tip(gt_r,c(sp_list[sp_list %in% gt$tip.label],out_v))
                gt_sub$node.label[gt_sub$node.label=="Root"]=100
                gt_sub$node.label=as.numeric(gt_sub$node.label)/100
                rooted_phy=c(rooted_phy,write.tree(gt_sub))
            }    
        } 
        
    } 
    phy=rooted_phy
    for (i in 1:n_retic)
    {    
        filename=deparse(substitute(sp_list))
        f_n=paste(filename,"_phylonet_genes_",i,"ret.nex",sep="")
        write("#NEXUS\n\nBEGIN TREES;",f_n)
        write(paste("Tree fixtr = ",write.tree(sp_phy_sub)),f_n,append=TRUE)
        d=data.frame(rep("Tree",length(phy)),paste("gt",1:length(phy),"=",sep=""),phy)
        write.table(d,f_n, quote = FALSE,row.names = FALSE, col.names=FALSE,append=TRUE)
        write("END;\n\nBEGIN PHYLONET;",f_n,append=TRUE)
        write(paste("InferNetwork_ML (all)",i,"-s fixtr -di -pl 10 -x 10 -b 0.9 -po;","\nEND;"),f_n,append=TRUE)
    }    
}    

phylonet_chunk_topo("MLrooted.tre","gene_trees_wboot_dna_mafft",I,1)

for (i in list(A,B,C))
{
    phylonet_chunk_topo("MLrooted.tre","gene_trees_wboot_dna_mafft",i,1)
}    



list_triples=function(taxa)
{
    res=rbind()
    ks = taxa
    for (i in 1:(length(ks)-2))
    {
        for (j in (i+1):(length(ks)-1))
        {
           for (k in (j+1):length(ks))
           {
               res=rbind(res,c(ks[i],ks[j], ks[k]))
               res=rbind(res,c(ks[i],ks[k], ks[j]))
               res=rbind(res,c(ks[j],ks[i], ks[k]))
               
            }
        }    
    }
    write.table(res,"tri.txt",quote=F,row.names=F,col.names=F)
    
}        
    
taxa_map=function(tree)
{
    phy=read.tree(tree)
    d=data.frame(phy$tip.label,paste("P",1:length(phy$tip.label),"_",phy$tip.label,sep=""))
    write.table(d,"taxa_map.txt",quote=F,row.names=F,col.names=F)
}
    
tripls_from_taxa=function(taxa,taxa_map)
{
    d=read.table(taxa_map)
    taxa=as.vector(d[d$V1 %in% taxa,"V2"])
    list_triples(taxa)
}

    
    







