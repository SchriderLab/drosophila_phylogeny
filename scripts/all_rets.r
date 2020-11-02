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
    }    
}    
    

net_combn=function(phy,clde_name)
{
    n=length(phy$tip.label)+phy$Nnode
    net_nodes=combn(n,2)
    net2score(phy,from=net_nodes[1,],to=net_nodes[2,],clade_name=clde_name)
}    





#Clade 2
#n1=c(19,13,8,9,16,4,4,4,21,20,9)
#n2=c(13,20,13,13,14,1,3,15,15,18,14)
#net2score(extract.clade(tt,245),from=c(n1,n2),to=c(n2,n1),"C2")
net_combn(extract.clade(tt,245),"netsC2")
net_combn(extract.clade(tt,289),"netsC4")
net_combn(extract.clade(tt,185),"netsC6")
net_combn(extract.clade(tt,226),"netsC8")
net_combn(extract.clade(tt,204),"netsC9")

#Clade 4
n1=c(19,19,19,19,34,34,34)
n2=c(21,34,26,23,22,26,23)
net2score(extract.clade(tt,289),from=c(n1,n2),to=c(n2,n1),"C4")