#Generates venn diagrams for microbial hub families

## Libraries ##
library(Hmisc)
library(igraph)
library(ggplot2)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(tidyr)
library(VennDiagram)

tema=theme(axis.text.x = element_text(color="black",size=17),# angle = 90, hjust = 0.5, vjust = 0.5),
           axis.text.y = element_text(color="black",size=17),
           axis.title = element_text(color="black",size=17),
           legend.text = element_text(color = "black",size=17),
           panel.border =element_rect(color = "white", fill = NA) , #element_blank(),
           panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "right",
           strip.text.x = element_text(size = 15, color = "white", face = "bold"),
           strip.background = element_rect(color="black", fill="black", size=1.5, linetype="solid")
)

aes4=c("#024AAC","#1DACA8","#10B62F", "#E2E41C","#F48F06","#F2252C","#D140B7", "grey80","grey20")
getPalette=colorRampPalette(aes4) 


## Directories ##
#WINDOWS
setwd("C:/Users/victorfn/Google Drive/network_analysis_2018")
network.dir="C:/Users/victorfn/Google Drive/network_analysis_2018/github" 
diff.dir="C:/Users/victorfn/Google Drive/taxa_enrichment_2020/book_chapter_family" #result enrichment analysis by family
#dir.create(result.dir,showWarnings=F)
base.name= "CAM"
plot.dir=paste(network.dir,"venn",sep = "/") ; dir.create(plot.dir)

## Select samples to work with ##
compartment=c("soil","roo.zone.soil", "root.endosphere","rhizosphere","leaf.endosphere","phyllosphere")[c(3:6)]
plants=c("Agave.salmiana","Agave.tequilana","Cacti")[c(1,2,3)]

load(file = paste(getwd(),paste(base.name,"OTU_table",sep = "_"),sep = "/"))
load(file = paste(getwd(),paste(base.name,"metadata_table",sep = "_"),sep = "/"))
load(file = paste(getwd(),paste(base.name,"taxonomy_table",sep = "_"),sep = "/"))

## Generate venn diagrams and overlaps  ##
  
m="family"
hubs=data.frame()
grado=data.frame()
central=data.frame()

for (comp in compartment){ print(comp)
  
  if (comp=="soil"){
    
    nodes.at=read.delim(file = paste(network.dir,paste("wild",comp,"R_nodes","txt",sep = "."),sep = "/"))
    nodes.as=read.delim(file = paste(network.dir,paste("cultivated",comp,"R_nodes","txt",sep = "."),sep = "/"))
    at=as.character(nodes.at[nodes.at$random.hub==1, m])
    as=as.character(nodes.as[nodes.as$random.hub==1, m])
    at=at[at!=""] ; as=as[as!=""] 
    do=list(at,as)
    h=unique(taxon[taxon[,5] %in% calculate.overlap(do)$a3,c(2,3,4,5)])
    
    if (dim(h)[1]!=0){
      h=data.frame(row.names=1:dim(h)[1],h, status="hub", compartment=comp)
      hubs=rbind(hubs,h)} else { print(paste("No core hubs in", comp))}
   
     at=as.character(nodes.at[nodes.at$random.grado==1, m])
    as=as.character(nodes.as[nodes.as$random.grado==1, m])
    at=at[at!=""] ; as=as[as!=""] 
    do=list(at,as)
    g=unique(taxon[taxon[,5] %in% calculate.overlap(do)$a3,c(2,3,4,5)])
    if (dim(g)[1]!=0){
      g=data.frame(row.names=1:dim(g)[1],g, status="grado", compartment=comp)
      grado=rbind(grado,g) } else { print(paste("No core connected in", comp))}
    
    at=as.character(nodes.at[nodes.at$random.central==1, m])
    as=as.character(nodes.as[nodes.as$random.central==1, m])
    at=at[at!=""] ; as=as[as!=""] 
    do=list(at,as)
    c=unique(taxon[taxon[,5] %in% calculate.overlap(do)$a3,c(2,3,4,5)])
    if (dim(c)[1]!=0){
      c=data.frame(row.names=1:dim(c)[1],c, status="central", compartment=comp)
      central=rbind(central,c) } else { print(paste("No core central in", comp))}
    
    
  } else {
  
  nodes.at=read.delim(file = paste(network.dir,paste("Agave.tequilana",comp,"R_nodes","txt",sep = "."),sep = "/"))
  nodes.as=read.delim(file = paste(network.dir,paste("Agave.salmiana",comp,"R_nodes","txt",sep = "."),sep = "/"))
  nodes.ca=read.delim(file = paste(network.dir,paste("Cacti",comp,"R_nodes","txt",sep = "."),sep = "/"))
  
  at=as.character(nodes.at[nodes.at$random.hub==1, m])
  as=as.character(nodes.as[nodes.as$random.hub==1, m])
  ca=as.character(nodes.ca[nodes.ca$random.hub==1, m])
  
  at=at[at!=""] ; as=as[as!=""] ; ca=ca[ca!=""]
  
  do=list(at,as,ca)
  
  h=unique(taxon[taxon[,5] %in% calculate.overlap(do)$a5,c(2,3,4,5)])
  
  
  if (dim(h)[1]!=0){
  h=data.frame(row.names=1:dim(h)[1],h, status="hub", compartment=comp)
  hubs=rbind(hubs,h)} else { print(paste("No core hubs in", comp))}
  
  h=calculate.overlap(do)$a5
  
  venn.diagram(do, filename = paste(plot.dir,paste("hub","venn_taxa",comp,m,"png",sep = "."),sep = "/"), 
               fill=getPalette(5)[c(1,2,4)], alpha=0.50, cex=4.5,cat.cex = rep(0.5,3),
               category.names=c("","",""))
  
  
  at=as.character(nodes.at[nodes.at$random.grado==1, m])
  as=as.character(nodes.as[nodes.as$random.grado==1, m])
  ca=as.character(nodes.ca[nodes.ca$random.grado==1, m])
  at=at[at!=""] ; as=as[as!=""] ; ca=ca[ca!=""]
  do=list(at,as,ca)
  
  g=unique(taxon[taxon[,5] %in% calculate.overlap(do)$a5,c(2,3,4,5)])
  
  if (dim(g)[1]!=0){
  g=data.frame(row.names=1:dim(g)[1],g, status="grado", compartment=comp)
  grado=rbind(grado,g) } else { print(paste("No core connected in", comp))}
  
  venn.diagram(do, filename = paste(plot.dir,paste("grado","venn_taxa",comp,m,"png",sep = "."),sep = "/"), 
               fill=getPalette(5)[c(1,2,4)], alpha=0.50, cex=4.5,cat.cex = rep(0.5,3),
               category.names=c("","",""))
  
  
  at=as.character(nodes.at[nodes.at$random.central==1, m])
  as=as.character(nodes.as[nodes.as$random.central==1, m])
  ca=as.character(nodes.ca[nodes.ca$random.central==1, m])
  at=at[at!=""] ; as=as[as!=""] ; ca=ca[ca!=""]
  do=list(at,as,ca)
  
  c=unique(taxon[taxon[,5] %in% calculate.overlap(do)$a5,c(2,3,4,5)])
  
  if (dim(c)[1]!=0){
  c=data.frame(row.names=1:dim(c)[1],c, status="central", compartment=comp)
  central=rbind(central,c) } else { print(paste("No core central in", comp))}
  
  venn.diagram(do, filename = paste(plot.dir,paste("central","venn_taxa",comp,m,"png",sep = "."),sep = "/"), 
               fill=getPalette(5)[c(1,2,4)], alpha=0.50, cex=4.5,cat.cex = rep(0.5,3),
               category.names=c("","",""))
 }              
  
}


## Generate heatmaps of core hub, central and connected taxa. ##
library(pheatmap)

## Select samples to work with ##
compartment=c("soil","roo.zone.soil", "root.endosphere","rhizosphere","leaf.endosphere","phyllosphere")[c(3:6)]
plant=c("Agave.salmiana","Agave.tequilana","Cacti")[c(1,2,3)]


for(comp in compartment){

data=rbind(hubs[hubs$compartment==comp,],grado[grado$compartment==comp,],central[central$compartment==comp,] )

grupo=paste(comp,plant,sep = ".")
sdiff=rbind(read.delim(file = paste(diff.dir,paste("Bacteria","family","logPvalues","season","txt",sep = "."), sep = "/")[1]),
            read.delim(file = paste(diff.dir,paste("Fungi","family","logPvalues","season","txt",sep = "."), sep = "/")[1]))

#Enrichment by compartment and species
pdiff=rbind(read.delim(paste(diff.dir,paste("Bacteria","family","logPvalues",comp,"txt",sep="."),sep = "/")[1]),
            read.delim(paste(diff.dir,paste("Fungi","family","logPvalues",comp,"txt",sep="."),sep = "/")[1]))


if (comp=="phyllosphere"){
  
  #Enrichment by compartment and species
  cdiff.at=rbind(read.delim(paste(diff.dir,paste("Bacteria","family","logPvalues","Agave.tequilana","epis","txt",sep="."),sep = "/")[1]),
              read.delim(paste(diff.dir,paste("Fungi","family","logPvalues","Agave.tequilana","epis","txt",sep="."),sep = "/")[1]))
  
  cdiff.as=rbind(read.delim(paste(diff.dir,paste("Bacteria","family","logPvalues","Agave.salmiana","epis","txt",sep="."),sep = "/")[1]),
              read.delim(paste(diff.dir,paste("Fungi","family","logPvalues","Agave.salmiana","epis","txt",sep="."),sep = "/")[1]))
  
  cdiff.ca=rbind(read.delim(paste(diff.dir,paste("Bacteria","family","logPvalues","Cacti","epis","txt",sep="."),sep = "/")[1]),
              read.delim(paste(diff.dir,paste("Fungi","family","logPvalues","Cacti","epis","txt",sep="."),sep = "/")[1]))
  
  
  
  #OTUs enriched in seasons
  data$dry.at=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[2]] >= -log10(0.05),]),1,0)
  data$dry.as=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[1]] >= -log10(0.05),]),1,0)
  data$dry.ca=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[3]] >= -log10(0.05),]),1,0)
  
  data$rainy.at=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[2]] <= log10(0.05),]),1,0)
  data$rainy.as=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[1]] <= log10(0.05),]),1,0)
  data$rainy.ca=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[3]] <= log10(0.05),]),1,0)
  
  #OTUS enriched in compartments
  data$root.at=ifelse(data$family %in%  rownames(cdiff.at[cdiff.at[,"phyllosphere...rhizosphere"] >= -log10(0.05),]),1,0)
  data$root.as=ifelse(data$family %in%  rownames(cdiff.as[cdiff.as[,"phyllosphere...rhizosphere"] >= -log10(0.05),]),1,0)
  data$root.ca=ifelse(data$family %in%  rownames(cdiff.ca[cdiff.ca[,"phyllosphere...rhizosphere"] >= -log10(0.05),]),1,0)
  
  data$soil.at=ifelse(data$family %in%  rownames(cdiff.at[cdiff.at[,"phyllosphere...soil"] >= -log10(0.05),]),1,0)
  data$soil.as=ifelse(data$family %in%  rownames(cdiff.as[cdiff.as[,"phyllosphere...soil"] >= -log10(0.05),]),1,0)
  data$soil.ca=ifelse(data$family %in%  rownames(cdiff.ca[cdiff.ca[,"phyllosphere...soil"] >= -log10(0.05),]),1,0)
  
  #OTUS enriched in spcies
  data$cultivated.as=ifelse(data$family %in% rownames(pdiff[pdiff[,"Agave.salmiana...Agave.tequilana"] <= log10(0.05),]),1,0)
  data$cultivated.ca=ifelse(data$family %in% rownames(pdiff[pdiff[,"Cacti...Agave.tequilana"] <= log10(0.05),]),1,0)
  
  data$native.as=ifelse(data$family %in% rownames(pdiff[pdiff[,"Agave.salmiana...Agave.tequilana"] >= -log10(0.05),]),1,0)
  data$native.ca=ifelse(data$family %in% rownames(pdiff[pdiff[,"Cacti...Agave.tequilana"] >= -log10(0.05),]),1,0)
  
}


if (comp=="rhizosphere"){
  
  #Enrichment by compartment and species
  cdiff.at=rbind(read.delim(paste(diff.dir,paste("Bacteria","family","logPvalues","Agave.tequilana","epis","txt",sep="."),sep = "/")[1]),
                 read.delim(paste(diff.dir,paste("Fungi","family","logPvalues","Agave.tequilana","epis","txt",sep="."),sep = "/")[1]))
  
  cdiff.as=rbind(read.delim(paste(diff.dir,paste("Bacteria","family","logPvalues","Agave.salmiana","epis","txt",sep="."),sep = "/")[1]),
                 read.delim(paste(diff.dir,paste("Fungi","family","logPvalues","Agave.salmiana","epis","txt",sep="."),sep = "/")[1]))
  
  cdiff.ca=rbind(read.delim(paste(diff.dir,paste("Bacteria","family","logPvalues","Cacti","epis","txt",sep="."),sep = "/")[1]),
                 read.delim(paste(diff.dir,paste("Fungi","family","logPvalues","Cacti","epis","txt",sep="."),sep = "/")[1]))
  
  
  
  #OTUs enriched in seasons
  data$dry.at=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[2]] >= -log10(0.05),]),1,0)
  data$dry.as=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[1]] >= -log10(0.05),]),1,0)
  data$dry.ca=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[3]] >= -log10(0.05),]),1,0)
  
  data$rainy.at=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[2]] <= log10(0.05),]),1,0)
  data$rainy.as=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[1]] <= log10(0.05),]),1,0)
  data$rainy.ca=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[3]] <= log10(0.05),]),1,0)
  
  #OTUS enriched in compartments
  data$leaf.at=ifelse(data$family %in%  rownames(cdiff.at[cdiff.at[,"phyllosphere...rhizosphere"] <= log10(0.05),]),1,0)
  data$leaf.as=ifelse(data$family %in%  rownames(cdiff.as[cdiff.as[,"phyllosphere...rhizosphere"] <= log10(0.05),]),1,0)
  data$leaf.ca=ifelse(data$family %in%  rownames(cdiff.ca[cdiff.ca[,"phyllosphere...rhizosphere"] <= log10(0.05),]),1,0)
  
  data$soil.at=ifelse(data$family %in%  rownames(cdiff.at[cdiff.at[,"rhizosphere...soil"] >= -log10(0.05),]),1,0)
  data$soil.as=ifelse(data$family %in%  rownames(cdiff.as[cdiff.as[,"rhizosphere...soil"] >= -log10(0.05),]),1,0)
  data$soil.ca=ifelse(data$family %in%  rownames(cdiff.ca[cdiff.ca[,"rhizosphere...soil"] >= -log10(0.05),]),1,0)
  
  #OTUS enriched in spcies
  data$cultivated.as=ifelse(data$family %in% rownames(pdiff[pdiff[,"Agave.salmiana...Agave.tequilana"] <= log10(0.05),]),1,0)
  data$cultivated.ca=ifelse(data$family %in% rownames(pdiff[pdiff[,"Cacti...Agave.tequilana"] <= log10(0.05),]),1,0)
  
  data$native.as=ifelse(data$family %in% rownames(pdiff[pdiff[,"Agave.salmiana...Agave.tequilana"] >= -log10(0.05),]),1,0)
  data$native.ca=ifelse(data$family %in% rownames(pdiff[pdiff[,"Cacti...Agave.tequilana"] >= -log10(0.05),]),1,0)
  
}


if (comp=="leaf.endosphere"){
  
  #Enrichment by compartment and species
  cdiff.at=rbind(read.delim(paste(diff.dir,paste("Bacteria","family","logPvalues","Agave.tequilana","endo","txt",sep="."),sep = "/")[1]),
                 read.delim(paste(diff.dir,paste("Fungi","family","logPvalues","Agave.tequilana","endo","txt",sep="."),sep = "/")[1]))
  
  cdiff.as=rbind(read.delim(paste(diff.dir,paste("Bacteria","family","logPvalues","Agave.salmiana","endo","txt",sep="."),sep = "/")[1]),
                 read.delim(paste(diff.dir,paste("Fungi","family","logPvalues","Agave.salmiana","endo","txt",sep="."),sep = "/")[1]))
  
  cdiff.ca=rbind(read.delim(paste(diff.dir,paste("Bacteria","family","logPvalues","Cacti","endo","txt",sep="."),sep = "/")[1]),
                 read.delim(paste(diff.dir,paste("Fungi","family","logPvalues","Cacti","endo","txt",sep="."),sep = "/")[1]))
  
  
  
  #OTUs enriched in seasons
  data$dry.at=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[2]] >= -log10(0.05),]),1,0)
  data$dry.as=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[1]] >= -log10(0.05),]),1,0)
  data$dry.ca=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[3]] >= -log10(0.05),]),1,0)
  
  data$rainy.at=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[2]] <= log10(0.05),]),1,0)
  data$rainy.as=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[1]] <= log10(0.05),]),1,0)
  data$rainy.ca=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[3]] <= log10(0.05),]),1,0)
  
  #OTUS enriched in compartments
  data$root.at=ifelse(data$family %in%  rownames(cdiff.at[cdiff.at[,"leaf.endosphere...root.endosphere"] >= -log10(0.05),]),1,0)
  data$root.as=ifelse(data$family %in%  rownames(cdiff.as[cdiff.as[,"leaf.endosphere...root.endosphere"] >= -log10(0.05),]),1,0)
  data$root.ca=ifelse(data$family %in%  rownames(cdiff.ca[cdiff.ca[,"leaf.endosphere...root.endosphere"] >= -log10(0.05),]),1,0)
  
  data$soil.at=ifelse(data$family %in%  rownames(cdiff.at[cdiff.at[,"leaf.endosphere...soil"] >= -log10(0.05),]),1,0)
  data$soil.as=ifelse(data$family %in%  rownames(cdiff.as[cdiff.as[,"leaf.endosphere...soil"] >= -log10(0.05),]),1,0)
  data$soil.ca=ifelse(data$family %in%  rownames(cdiff.ca[cdiff.ca[,"leaf.endosphere...soil"] >= -log10(0.05),]),1,0)
  
  #OTUS enriched in spcies
  data$cultivated.as=ifelse(data$family %in% rownames(pdiff[pdiff[,"Agave.salmiana...Agave.tequilana"] <= log10(0.05),]),1,0)
  data$cultivated.ca=ifelse(data$family %in% rownames(pdiff[pdiff[,"Cacti...Agave.tequilana"] <= log10(0.05),]),1,0)
  
  data$native.as=ifelse(data$family %in% rownames(pdiff[pdiff[,"Agave.salmiana...Agave.tequilana"] >= -log10(0.05),]),1,0)
  data$native.ca=ifelse(data$family %in% rownames(pdiff[pdiff[,"Cacti...Agave.tequilana"] >= -log10(0.05),]),1,0)
  
}

if (comp=="root.endosphere"){
  
  #Enrichment by compartment and species
  cdiff.at=rbind(read.delim(paste(diff.dir,paste("Bacteria","family","logPvalues","Agave.tequilana","endo","txt",sep="."),sep = "/")[1]),
                 read.delim(paste(diff.dir,paste("Fungi","family","logPvalues","Agave.tequilana","endo","txt",sep="."),sep = "/")[1]))
  
  cdiff.as=rbind(read.delim(paste(diff.dir,paste("Bacteria","family","logPvalues","Agave.salmiana","endo","txt",sep="."),sep = "/")[1]),
                 read.delim(paste(diff.dir,paste("Fungi","family","logPvalues","Agave.salmiana","endo","txt",sep="."),sep = "/")[1]))
  
  cdiff.ca=rbind(read.delim(paste(diff.dir,paste("Bacteria","family","logPvalues","Cacti","endo","txt",sep="."),sep = "/")[1]),
                 read.delim(paste(diff.dir,paste("Fungi","family","logPvalues","Cacti","endo","txt",sep="."),sep = "/")[1]))
  
  
  
  #OTUs enriched in seasons
  data$dry.at=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[2]] >= -log10(0.05),]),1,0)
  data$dry.as=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[1]] >= -log10(0.05),]),1,0)
  data$dry.ca=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[3]] >= -log10(0.05),]),1,0)
  
  data$rainy.at=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[2]] <= log10(0.05),]),1,0)
  data$rainy.as=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[1]] <= log10(0.05),]),1,0)
  data$rainy.ca=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[3]] <= log10(0.05),]),1,0)
  
  #OTUS enriched in compartments
  data$leaf.at=ifelse(data$family %in%  rownames(cdiff.at[cdiff.at[,"leaf.endosphere...root.endosphere"] <= log10(0.05),]),1,0)
  data$leaf.as=ifelse(data$family %in%  rownames(cdiff.as[cdiff.as[,"leaf.endosphere...root.endosphere"] <= log10(0.05),]),1,0)
  data$leaf.ca=ifelse(data$family %in%  rownames(cdiff.ca[cdiff.ca[,"leaf.endosphere...root.endosphere"] <= log10(0.05),]),1,0)
  
  data$soil.at=ifelse(data$family %in%  rownames(cdiff.at[cdiff.at[,"root.endosphere...soil"] >= -log10(0.05),]),1,0)
  data$soil.as=ifelse(data$family %in%  rownames(cdiff.as[cdiff.as[,"root.endosphere...soil"] >= -log10(0.05),]),1,0)
  data$soil.ca=ifelse(data$family %in%  rownames(cdiff.ca[cdiff.ca[,"root.endosphere...soil"] >= -log10(0.05),]),1,0)
  
  #OTUS enriched in spcies
  data$cultivated.as=ifelse(data$family %in% rownames(pdiff[pdiff[,"Agave.salmiana...Agave.tequilana"] <= log10(0.05),]),1,0)
  data$cultivated.ca=ifelse(data$family %in% rownames(pdiff[pdiff[,"Cacti...Agave.tequilana"] <= log10(0.05),]),1,0)
  
  data$native.as=ifelse(data$family %in% rownames(pdiff[pdiff[,"Agave.salmiana...Agave.tequilana"] >= -log10(0.05),]),1,0)
  data$native.ca=ifelse(data$family %in% rownames(pdiff[pdiff[,"Cacti...Agave.tequilana"] >= -log10(0.05),]),1,0)
  
  

}


if (comp=="soil"){
  
  #Enrichment by compartment and species
  cdiff.at=rbind(read.delim(paste(diff.dir,paste("Bacteria","family","logPvalues","Agave.tequilana","epis","txt",sep="."),sep = "/")[1]),
                 read.delim(paste(diff.dir,paste("Fungi","family","logPvalues","Agave.tequilana","epis","txt",sep="."),sep = "/")[1]))
  
  cdiff.as=rbind(read.delim(paste(diff.dir,paste("Bacteria","family","logPvalues","Agave.salmiana","epis","txt",sep="."),sep = "/")[1]),
                 read.delim(paste(diff.dir,paste("Fungi","family","logPvalues","Agave.salmiana","epis","txt",sep="."),sep = "/")[1]))
  
  cdiff.ca=rbind(read.delim(paste(diff.dir,paste("Bacteria","family","logPvalues","Cacti","epis","txt",sep="."),sep = "/")[1]),
                 read.delim(paste(diff.dir,paste("Fungi","family","logPvalues","Cacti","epis","txt",sep="."),sep = "/")[1]))
  
  #OTUs enriched in seasons
  data$dry.at=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[2]] >= -log10(0.05),]),1,0)
  data$dry.as=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[1]] >= -log10(0.05),]),1,0)
  data$dry.ca=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[3]] >= -log10(0.05),]),1,0)
  
  data$rainy.at=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[2]] <= log10(0.05),]),1,0)
  data$rainy.as=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[1]] <= log10(0.05),]),1,0)
  data$rainy.ca=ifelse(data$family %in%  rownames(sdiff[sdiff[,grupo[3]] <= log10(0.05),]),1,0)
  
  #OTUS enriched in compartments
  data$root.at=ifelse(data$family %in%  rownames(cdiff.at[cdiff.at[,"rhizosphere...soil"] <= log10(0.05),]),1,0)
  data$root.as=ifelse(data$family %in%  rownames(cdiff.as[cdiff.as[,"rhizosphere...soil"] <= log10(0.05),]),1,0)
  data$root.ca=ifelse(data$family %in%  rownames(cdiff.ca[cdiff.ca[,"rhizosphere...soil"] <= log10(0.05),]),1,0)
  
  data$leaf.at=ifelse(data$family %in%  rownames(cdiff.at[cdiff.at[,"phyllosphere...soil"] <= log10(0.05),]),1,0)
  data$leaf.as=ifelse(data$family %in%  rownames(cdiff.as[cdiff.as[,"phyllosphere...soil"] <= log10(0.05),]),1,0)
  data$leaf.ca=ifelse(data$family %in%  rownames(cdiff.ca[cdiff.ca[,"phyllosphere...soil"] <= log10(0.05),]),1,0)
  
  #OTUS enriched in spcies
  data$cultivated.as=ifelse(data$family %in% rownames(pdiff[pdiff[,"Agave.salmiana...Agave.tequilana"] <= log10(0.05),]),1,0)
  data$cultivated.ca=ifelse(data$family %in% rownames(pdiff[pdiff[,"Cacti...Agave.tequilana"] <= log10(0.05),]),1,0)
  
  data$native.as=ifelse(data$family %in% rownames(pdiff[pdiff[,"Agave.salmiana...Agave.tequilana"] >= -log10(0.05),]),1,0)
  data$native.ca=ifelse(data$family %in% rownames(pdiff[pdiff[,"Cacti...Agave.tequilana"] >= -log10(0.05),]),1,0)
  
  
}


data2=unique(data[,-c(1,2,3,5,6)])
rownames(data2)=data2$family
data2=data2[,-1]
data2=data2[rownames(data2)!="g",]


png(paste(plot.dir,paste(comp,"corehubs2enrichment_family","png",sep = "."),sep = "/"), 
    width = 6000, height = 6000, units = "px", pointsize = 15, bg = "white", res=1000, type="cairo")

print(pheatmap(data2, cluster_rows = F, cluster_cols = F, color = c("grey80","black"),border_color = "white",
         gaps_col=c(3,6,9,12,14), cellwidth = 10, cellheight = 10))


dev.off()

}





