library("ape")
library("stringr")
tt=read.tree("MLrooted.tre")



netout=function(phy,from,to)
{
    phy_txt=write.tree(phy)
    
    if(from > length(phy$tip.label))
    {
        phy_txt_from=gsub(strsplit(write.tree(extract.clade(phy,from)),split=";")[[1]],paste("((",strsplit(write.tree(extract.clade(phy,from)),split=";")[[1]],"),#H1)",sep=""),phy_txt,fixed=T)
        
    }else{
        
        phy_txt_from=gsub(phy$tip.label[from],paste("((",phy$tip.label[from],"),#H1)",sep=""),phy_txt,fixed=T)
    
    }
    
    if(to > length(phy$tip.label))
    {
     
        phy_txt_out=gsub(strsplit(write.tree(extract.clade(phy,to)),split=";")[[1]],paste("((",strsplit(write.tree(extract.clade(phy,to)),split=";")[[1]],")#H1)",sep=""),phy_txt_from,fixed=T)
    
    }else{
        
        phy_txt_out=gsub(phy$tip.label[to],paste("((",phy$tip.label[to],")#H1)",sep=""),phy_txt_from,fixed=T)
    
    }    
    
    if (str_count(phy_txt_out,"#")<2)
    {
        phy_txt_out = netout(phy,to,from)
        return(phy_txt_out)
    }else{
        return(phy_txt_out)
    }    
}    


net2score=function(phy,from,to,clade_name)
{
    phy$edge.length=NULL
    f_n=paste(clade_name,"_score_nets.txt",sep="")
    write(write.tree(phy),f_n)
    for (i in 1:length(from))
    {
        net_store=netout(phy,from=from[i],to=to[i])
        write(net_store,f_n,append=TRUE)
        net_store=netout(phy,from=to[i],to=from[i])
        write(net_store,f_n,append=TRUE)
    }    
}    
    

net_combn=function(phy,clde_name)
{
    n=length(phy$tip.label)+phy$Nnode
    net_nodes=combn(n,2)
    net2score(phy,from=net_nodes[1,],to=net_nodes[2,],clade_name=clde_name)
}    




#Subtree prep

#Clade1 

c1=drop.tip(extract.clade(tt,307),c("D_paulistorum","D_saltans"))
c2=drop.tip(extract.clade(tt,245),c("D_persimilis"))
c3=extract.clade(tt,256)
c4=drop.tip(extract.clade(tt,289),c("D_kurseongensis","D_carrolli","D_biarmipes","D_subpulchrella","D_orena","D_teissieri_273_3","D_sechellia","D_mauritiana","D_yakuba"))
c5=keep.tip(extract.clade(tt,265),c('D_asahinai','D_tani','D_triauraria','D_pectinifera','D_burlai','D_mayri','D_serrata','D_watanabei','D_kikkawai','D_kanapiae'))
c6=keep.tip(extract.clade(tt,185),c('D_cyrtoloma','D_planitibia','D_heteroneura','D_grimshawi','D_limitata','D_prolacticillia','D_primaeva','S_elmoi','S_pallida','S_montana'))
c7=keep.tip(extract.clade(tt,168),c('D_acanthoptera','D_pachea','D_bromeliae','D_repleta','D_hydei','D_pseudotalamancana','D_americana','D_virilis','D_montana','D_robusta'))
c8=keep.tip(extract.clade(tt,226),c('Z_africanus','Z_nigranus','Z_lachaisei','Z_vittiger','Z_capensis','Z_ornatusmayotte','Z_ghesquierei','Z_inermis','Z_kolodkinae','Z_tsacasi'))
c9=keep.tip(extract.clade(tt,204),c('D_kepulauana','D_neohypocausta','D_immigrans','D_pruinosa','D_ornatifrons','D_pallidipennis','D_funebris','D_guttifera','D_innubila','D_quadrilineata'))


#Clade 2
#n1=c(19,13,8,9,16,4,4,4,21,20,9)
#n2=c(13,20,13,13,14,1,3,15,15,18,14)
#net2score(extract.clade(tt,245),from=c(n1,n2),to=c(n2,n1),"C2")
net_combn(c1,"netsC1")
net_combn(c2,"netsC2")
net_combn(c3,"netsC3")
net_combn(c4,"netsC4")
net_combn(c5,"netsC5")
net_combn(c6,"netsC6")
net_combn(c7,"netsC7")
net_combn(c8,"netsC8")
net_combn(c9,"netsC9")

