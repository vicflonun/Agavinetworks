##get_hubs
#This code creates correlation networks of OTUs across differentof the agaves and cacti 
#microbiome and identify consistent microbial hubs and higly central and connected taxa 
#Data from (Coleman-Der et al, 2016;Fonseca-García, 2016; Partida-Martínez, 2018)

#output. Netwoks in gexf format, Lits with hubs, central and connected taxa in all the networks
#and dataframes with the consistent microbial hubs, connected and central taxa. 

## Libraries ##
library(Hmisc)
library(igraph)
library(ggplot2)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(tidyr)

tema=theme(axis.text.x = element_text(color="black",size=16, angle = 90, hjust = 0.5, vjust = 0.5),
           axis.text.y = element_text(color="black",size=16),
           axis.title = element_text(color="black",size=16),
           legend.text = element_text(color = "black",size=16),
           panel.border =element_rect(color = "black", fill = NA) , #element_blank(),
           panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "none")


aes4=c("#024AAC","#1DACA8","#10B62F", "#E2E41C","#F48F06","#F2252C","#D140B7", "grey80","grey20")
getPalette=colorRampPalette(aes4) 


## Directories ##
#WINDOWS
setwd("C:/Users/victorfn/Google Drive/network_analysis_2018")
network.dir="C:/Users/victorfn/Google Drive/network_analysis_2018/book_chapter_soil" 
diff.dir="C:/Users/victorfn/Google Drive/taxa_enrichment_2020/book_chapter_210420"
dir.create(result.dir,showWarnings=F)
base.name= "CAM"
plot.dir=paste(network.dir,"draws3",sep = "/") ; dir.create(plot.dir)

## Select samples to work with ##
compartment=c("soil","roo.zone.soil","rhizosphere","root.endosphere","leaf.endosphere","phyllosphere")[c(1)]
plants=c("Agave.salmiana","Agave.tequilana","Agave.deserti","Cacti")[c(1,2,4)]
sites=c("cultivated","wild" )[1]

load(file = paste(getwd(),paste(base.name,"OTU_table",sep = "_"),sep = "/"))
load(file = paste(getwd(),paste(base.name,"metadata_table",sep = "_"),sep = "/"))
load(file = paste(getwd(),paste(base.name,"taxonomy_table",sep = "_"),sep = "/"))


  for (comp in compartment) {

print(sites)

formas=c(21,22,23) ; names(formas)=c("Bacteria","Fungi","Archaea")
colores=c("#3DCCC9","#116DC8","#EF990D") ; names(colores)=c("Archaea","Bacteria","Fungi")


sdiff=rbind(read.delim(file = paste(diff.dir,paste("Bacteria","OTU.ID","logPvalues","season","txt",sep = "."), sep = "/")[1]),
            read.delim(file = paste(diff.dir,paste("Fungi","OTU.ID","logPvalues","season","txt",sep = "."), sep = "/")[1]))

#Network data

load(file = paste(network.dir,paste(sites,comp,"igraph",sep = "."),sep = "/"))
nodes=read.delim(file = paste(network.dir,paste(sites,comp,"R_nodes","txt",sep = "."),sep = "/"))
random=read.delim(file = paste(network.dir,paste(sites,comp,"random_networks","txt",sep = "."),sep = "/"))

#Enrichment by compartment and species

pdiff=rbind(read.delim(paste(diff.dir,paste("Bacteria","OTU.ID","logPvalues",comp,"txt",sep="."),sep = "/")[1]),
            read.delim(paste(diff.dir,paste("Fungi","OTU.ID","logPvalues",comp,"txt",sep="."),sep = "/")[1]))


#network relevance
nodes$network="none"
nodes[nodes$OTU.ID %in%  nodes[nodes$central==1 & nodes$random.central==1, "OTU.ID"],"network"] = "central"
nodes[nodes$OTU.ID %in%  nodes[nodes$grado==1 & nodes$random.grado==1, "OTU.ID"],"network"] = "connected"
nodes[nodes$OTU.ID %in%  nodes[nodes$hub==1 & nodes$random.hub==1, "OTU.ID"],"network"] = "hub"


#Enrichment season data
if (sites == "wild" & comp=="soil"){
  #Enrichment by compartment and species
  sal_cdiff=rbind(read.delim(paste(diff.dir,paste("Bacteria","OTU.ID","logPvalues","Agave.salmiana","epis","txt",sep="."),sep = "/")[1]),
              read.delim(paste(diff.dir,paste("Fungi","OTU.ID","logPvalues","Agave.salmiana","epis","txt",sep="."),sep = "/")[1]))
  
  cac_cdiff=rbind(read.delim(paste(diff.dir,paste("Bacteria","OTU.ID","logPvalues","Cacti","epis","txt",sep="."),sep = "/")[1]),
                  read.delim(paste(diff.dir,paste("Fungi","OTU.ID","logPvalues","Cacti","epis","txt",sep="."),sep = "/")[1]))
  
  
  #OTUs enriched in seasons
  
  r=unique(rownames(sdiff[sdiff[,"soil.Agave.salmiana"] >= -log10(0.05),]),rownames(sdiff[sdiff[,"soil.Cacti"] >= -log10(0.05),]))
  d=unique(rownames(sdiff[sdiff[,"soil.Agave.salmiana"] <= log10(0.05),]), rownames(sdiff[sdiff[,"soil.Cacti"] <= log10(0.05),]) )
  
  nodes$dry=ifelse(nodes$OTU.ID %in%  r,"enriched","not.enriched")
  nodes$rainy=ifelse(nodes$OTU.ID %in%  d,"enriched","not.enriched")
  
  #OTUS enriched in compartments
  
  p=unique(rownames(sal_cdiff[sal_cdiff[,"phyllosphere...soil"] <= log10(0.05),]),rownames(cac_cdiff[cac_cdiff[,"phyllosphere...soil"] <= log10(0.05),]))
  r=unique(rownames(sal_cdiff[sal_cdiff[,"rhizosphere...soil"] <= log10(0.05),]), rownames(cac_cdiff[cac_cdiff[,"rhizosphere...soil"] <= log10(0.05),]))
  
  
  nodes$phyllosphere=ifelse(nodes$OTU.ID %in%  p,"enriched","not.enriched")
  nodes$rhizosphere=ifelse(nodes$OTU.ID %in% r ,"enriched","not.enriched")
  
  e=unique(rownames(pdiff[pdiff[,"Agave.salmiana...Agave.tequilana"] >= -log10(0.05),]),rownames(pdiff[pdiff[,"Cacti...Agave.tequilana"] >= -log10(0.05),]))
  
  #OTUS enriched in plant
  nodes$Agave.tequilana=ifelse(nodes$OTU.ID %in% e,"enriched","not.enriched")
  
  #Asign network atributes
  vertex_attr(g, name="kingdom", index = V(g)) <- taxon[V(g)$name,1]
  vertex_attr(g, name="hub", index = V(g)) <- nodes$network
  vertex_attr(g, name="dry", index = V(g)) <- nodes$dry
  vertex_attr(g, name="rainy", index = V(g)) <- nodes$rainy
  vertex_attr(g, name="leaf", index = V(g)) <- nodes$phyllosphere
  vertex_attr(g, name="root", index = V(g)) <- nodes$rhizosphere
  vertex_attr(g, name="A.tequilana", index = V(g)) <- nodes$Agave.tequilana
  vertex_attr(g, name="Cacti", index = V(g)) <- 0
  edge_attr(g, name="peso",index=E(g)) <- ifelse(E(g)$weight>0, "grey50", "firebrick3")
  vertex_attr(g, name="hub2", index = V(g)) <- ifelse(nodes$network=="hub","hub","no_hub")
  vertex_attr(g, name="random", index = V(g)) <- random$HF/100
  
  comparar=names(vertex.attributes(g))[-c(1:3,9,10,11)]
}

  
#####################################################################################


if (sites == "cultivated" & comp=="soil"){
  #Enrichment by compartment and species
  cdiff=rbind(read.delim(paste(diff.dir,paste("Bacteria","OTU.ID","logPvalues","Agave.tequilana","epis","txt",sep="."),sep = "/")[1]),
              read.delim(paste(diff.dir,paste("Fungi","OTU.ID","logPvalues","Agave.tequilana","epis","txt",sep="."),sep = "/")[1]))
  #OTUs enriched in seasons
  nodes$dry=ifelse(nodes$OTU.ID %in%  rownames(sdiff[sdiff[,"soil.Agave.tequilana"] >= -log10(0.05),]),"enriched","not.enriched")
  nodes$rainy=ifelse(nodes$OTU.ID %in%  rownames(sdiff[sdiff[,"soil.Agave.tequilana"] <= log10(0.05),]),"enriched","not.enriched")
  #OTUS enriched in compartments
  nodes$phyllosphere=ifelse(nodes$OTU.ID %in%  rownames(cdiff[cdiff[,"phyllosphere...soil"] <= log10(0.05),]),"enriched","not.enriched")
  nodes$rhizosphere=ifelse(nodes$OTU.ID %in%  rownames(cdiff[cdiff[,"rhizosphere...soil"] <= log10(0.05),]),"enriched","not.enriched")
  #OTUS enriched in species
  
  e=unique(rownames(pdiff[pdiff[,"Agave.salmiana...Agave.tequilana"] <= log10(0.05),]),rownames(pdiff[pdiff[,"Cacti...Agave.tequilana"] <= log10(0.05),]))
  nodes$Agave.salmiana=ifelse(nodes$OTU.ID %in% e,"enriched","not.enriched")

  #Asign network atributes
  vertex_attr(g, name="kingdom", index = V(g)) <- taxon[V(g)$name,1]
  vertex_attr(g, name="hub", index = V(g)) <- nodes$network
  vertex_attr(g, name="dry", index = V(g)) <- nodes$dry
  vertex_attr(g, name="rainy", index = V(g)) <- nodes$rainy
  vertex_attr(g, name="leaf", index = V(g)) <- nodes$phyllosphere
  vertex_attr(g, name="root", index = V(g)) <- nodes$rhizosphere
  vertex_attr(g, name="A.salmiana", index = V(g)) <- nodes$Agave.salmiana
  vertex_attr(g, name="Cacti", index = V(g)) <- 0
  edge_attr(g, name="peso",index=E(g)) <- ifelse(E(g)$weight>0, "grey50", "red")
  vertex_attr(g, name="hub2", index = V(g)) <- ifelse(nodes$network=="hub","hub","no_hub")
  vertex_attr(g, name="random", index = V(g)) <- random$HF
  
  comparar=names(vertex.attributes(g))[-c(1:3,9,10,11)]
}
          
  result2=data.frame(row.names = c("Nodes","Hub","Bacteria.E","Fungal.E","Dry","Rainy","vs leaf","vs Soil", "vsplant1","vsplant2"))
  
  #nodos del total de nodos    
  p=length(grep("p",nodes$OTU.ID,value = T)) / dim(nodes)[1] * 100
  f=length(grep("f",nodes$OTU.ID,value = T)) / dim(nodes)[1] * 100
  d=0
  result2[1,c("Prokaryotic","Fungal","No_enriched")]=c(p,f,d)
  
  #hubs del total de nodos
  p=length(grep("p",nodes[nodes$network=="hub","OTU.ID"],value = T)) / dim(nodes)[1] * 100
  f=length(grep("f",nodes[nodes$network=="hub","OTU.ID"],value = T)) / dim(nodes)[1] * 100
  d=100-p-f
  result2[2,c("Prokaryotic","Fungal","No_enriched")]=c(p,f,d)
  

  elist=data.frame(get.edgelist(g))
  con=paste(elist$X1,elist$X2,sep = "/")
  #nodos bacteria
  p=length(grep("f_",con[grep("p_",con)], invert = T))/length(grep("p_",con))*100
  f=length(grep("f_",con[grep("p_",con)]))/length(grep("p_",con))*100
  d=0
  result2[3,c("Prokaryotic","Fungal","No_enriched")]=c(p,f,d)
  
  #nodos hongos
  f=length(grep("p_",con[grep("f_",con)], invert = T))/length(grep("f_",con))*100
  p=length(grep("p_",con[grep("f_",con)]))/length(grep("f_",con))*100
  d=0
  result2[4,c("Prokaryotic","Fungal","No_enriched")]=c(p,f,d)
  
  #en secas del total de hubs
  p=length(grep("p",nodes[nodes$dry=="enriched" & nodes$network=="hub","OTU.ID"],value = T)) / length(nodes[nodes$network=="hub" ,"OTU.ID"]) *100
  f=length(grep("f",nodes[nodes$dry=="enriched" & nodes$network=="hub","OTU.ID"],value = T)) / length(nodes[nodes$network=="hub" ,"OTU.ID"]) *100
  d=100-p-f
  result2[5,c("Prokaryotic","Fungal","No_enriched")]=c(p,f,d)
  
  #en lluvias del total de hubs
  p=length(grep("p",nodes[nodes$rainy=="enriched" & nodes$network=="hub","OTU.ID"],value = T)) / length(nodes[nodes$network=="hub" ,"OTU.ID"]) *100
  f=length(grep("f",nodes[nodes$rainy=="enriched" & nodes$network=="hub","OTU.ID"],value = T)) / length(nodes[nodes$network=="hub" ,"OTU.ID"]) *100
  d=100-p-f
  result2[6,c("Prokaryotic","Fungal","No_enriched")]=c(p,f,d)
  
  #en compartimento (root)
  p=length(grep("p",nodes[nodes[,19]=="enriched" & nodes$network=="hub","OTU.ID"],value = T)) / length(nodes[nodes$network=="hub" ,"OTU.ID"]) *100
  f=length(grep("f",nodes[nodes[,19]=="enriched" & nodes$network=="hub","OTU.ID"],value = T)) / length(nodes[nodes$network=="hub" ,"OTU.ID"]) *100
  d=100-p-f
  result2[7,c("Prokaryotic","Fungal","No_enriched")]=c(p,f,d)
  
  #en compartimento (soil)
  p=length(grep("p",nodes[nodes[,20]=="enriched" & nodes$network=="hub","OTU.ID"],value = T)) / length(nodes[nodes$network=="hub" ,"OTU.ID"]) *100
  f=length(grep("f",nodes[nodes[,20]=="enriched" & nodes$network=="hub","OTU.ID"],value = T)) / length(nodes[nodes$network=="hub" ,"OTU.ID"]) *100
  d=100-p-f
  result2[8,c("Prokaryotic","Fungal","No_enriched")]=c(p,f,d)
  
  #en planta (1)
  p=length(grep("p",nodes[nodes[,21]=="enriched" & nodes$network=="hub","OTU.ID"],value = T)) / length(nodes[nodes$network=="hub" ,"OTU.ID"]) *100
  f=length(grep("f",nodes[nodes[,21]=="enriched" & nodes$network=="hub","OTU.ID"],value = T)) / length(nodes[nodes$network=="hub" ,"OTU.ID"]) *100
  d=100-p-f
  result2[9,c("Prokaryotic","Fungal","No_enriched")]=c(p,f,d)
  
  #en planta (2)
  #p=length(grep("p",nodes[nodes[,22]=="enriched" & nodes$network=="hub","OTU.ID"],value = T)) / length(nodes[nodes$network=="hub" ,"OTU.ID"]) *100
  #f=length(grep("f",nodes[nodes[,22]=="enriched" & nodes$network=="hub","OTU.ID"],value = T)) / length(nodes[nodes$network=="hub" ,"OTU.ID"]) *100
  #d=100-p-f
  #result2[10,c("Prokaryotic","Fungal","No_enriched")]=c(p,f,d)
  
  result2$enrichment=rownames(result2)
  result2t=gather(result2, value = "proportion", key="taxa",c("Prokaryotic","Fungal","No_enriched"))
  result2t$taxa=factor(result2t$taxa, levels = c("No_enriched","Prokaryotic","Fungal"))
  result2t$enrichment=factor(result2t$enrichment,levels=c("Nodes","Hub","Bacteria.E","Fungal.E","Dry","Rainy","vs leaf","vs Soil", "vsplant1","vsplant2"))
  
  
  
  pdf(file = paste(plot.dir,paste(comp,sites,"pdf",sep = "."),sep = "/"))
  for (i in comparar){
    
    a=ggnet2(g, size = 0, node.shape = i, node.size = "hub2", node.color="kingdom",
           edge.color = "peso", edge.size =0.7, mode="fruchtermanreingold")+
      geom_point(aes(shape=shape,fill=color, size=size))+
      scale_shape_manual(values = c(24,21))+
      scale_fill_manual(values =  colores[unique(vertex_attr(g, name="kingdom", index = V(g)))])+
      scale_size_manual(values = c(6,3))+ggtitle(paste("enriched in or against" , i, sep = " "))
    
    print(a)
  } 
 
  
  b=ggplot(result2t, aes(x=enrichment, y=proportion, fill=taxa))+
    geom_bar(stat = "identity")+tema+
    scale_fill_manual(values=c("grey80","grey50","grey20"))
  print(b)

  dev.off()
  
  a=ggnet2(g, size = 0, node.shape = "dry", node.size = "hub2", node.color="kingdom",
           edge.color = "peso", edge.size =0.85, mode="fruchtermanreingold",
           layout.par = list(niter=600,cell.jitter=0.75))+
    geom_point(aes(shape=shape,fill=color, size=size))+
    scale_shape_manual(values = c(24,21))+
    scale_fill_manual(values =  colores[unique(vertex_attr(g, name="kingdom", index = V(g)))])+
    scale_size_manual(values = c(6,3))+ggtitle("enriched in dry")
  
  png(paste(plot.dir,paste(comp,sites,"net","png",sep = "."),sep = "/"), 
      width = 1200, height = 700, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
 print(a)
  
  dev.off()
  
  
 a=ggplot(result2t, aes(x=enrichment, y=proportion, fill=taxa))+
    geom_bar(stat = "identity", width = 0.8)+tema+
    scale_fill_manual(values=c("grey80","#116DC8","#EF990D"))+
    #facet_wrap(sep~.,scales = "free_x")+
    scale_y_continuous(breaks = seq(0,100,20))
 
  png(paste(plot.dir,paste(comp,sites,"freq","png",sep = "."),sep = "/"), 
      width = 900, height = 650, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")
  
 print(a)
  
  dev.off()
  
  
  write.table(nodes, file = paste(network.dir,paste(sites,comp,"enriched.nodes","txt",sep = "."),sep = "/"),
              row.names = T, quote = F, sep = "\t", col.names = T)
  
  
  
}


## Select samples to work with ##
compartment=c("soil","roo.zone.soil","rhizosphere","root.endosphere","leaf.endosphere","phyllosphere")[c(1,3,4,5,6)]
plants=c("Agave.salmiana","Agave.tequilana","Agave.deserti","Cacti")[c(1,2,4)]

for (comp in compartment) {
  
  for(plant in plants){
    print("________________________________")
    grupo=paste(comp,plant,sep = ".")
    print(grupo)

    load(file = paste(network.dir,paste(plant,comp,"igraph",sep = "."),sep = "/"))
    
    print(paste("nodes",length(V(g))))
    print(paste("edges",length(E(g))))
    
    print("________________________________")
    
    
  }
}
