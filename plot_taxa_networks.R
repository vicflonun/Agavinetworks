#Generate barplots for the taxonomy (class level) of the hub OTUs

## Libraries ##
library(Hmisc)
library(igraph)
library(rgexf)
library(ggplot2)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(tidyr)
library(VennDiagram)

tema=theme(axis.text.x = element_text(color="black",size=14), #angle=90,hjust=0.95,vjust=0.2),
           axis.text.y = element_text(color="black",size=14),
           axis.title = element_text(color="black",size=14),
           legend.text = element_text(color = "black",size=13),
           strip.text.x = element_text(size=15, color="black"),
           panel.border =element_rect(color = "white", fill = NA) ,
           panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "right")


aes4=c("#024AAC","#1DACA8","#10B62F", "#E2E41C","#F48F06","#F2252C","#D140B7", "grey80","grey20")
getPalette=colorRampPalette(aes4) 


## Directories ##
#WINDOWS
setwd("C:/Users/victorfn/Google Drive/network_analysis_2018")
network.dir="C:/Users/victorfn/Google Drive/network_analysis_2018/github" 
base.name= "CAM"
plot.dir=paste(network.dir,"plots",sep = "/") ; dir.create(plot.dir)

## Select samples to work with ##
compartment=c("soil","roo.zone.soil", "root.endosphere","rhizosphere","leaf.endosphere","phyllosphere")[c(3:6)]
plants=c("Agave.salmiana","Agave.tequilana","Agave.deserti","Cacti")[c(1,2,4)]

load(file = paste(getwd(),paste(base.name,"OTU_table",sep = "_"),sep = "/"))
load(file = paste(getwd(),paste(base.name,"metadata_table",sep = "_"),sep = "/"))
load(file = paste(getwd(),paste(base.name,"taxonomy_table",sep = "_"),sep = "/"))

#What kind of vertices do you wanna plot?
vert=c("all","hub","central","connected")[2]

#Please change the taxonomy rank manually (phylum,class recomended)

#generate plot data

plot=data.frame()
for(plant in plants){
  for (comp in compartment) {
   
    print(paste(plant,comp))
    
    #Network data
    load(file = paste(network.dir,paste(plant,comp,"igraph",sep = "."),sep = "/"))
    nodes=read.delim(file = paste(network.dir,paste(plant,comp,"R_nodes","txt",sep = "."),sep = "/"))
    random=read.delim(file = paste(network.dir,paste(plant,comp,"random_networks","txt",sep = "."),sep = "/"))
    
    #select data
    if (vert=="hub"){
      nodes=nodes[nodes$random.hub==1,] 
    } else if (vert=="central") {
      nodes=nodes[nodes$random.central==1,] 
    } else if (vert=="connected"){
      nodes=nodes[nodes$random.grado==1,] 
    }
    
    #set levels
    nodes$class=factor(nodes$class, levels = unique(nodes$class))
  
    #Calculate relative abundances
    sub=data.frame(row.names=1:length(unique(nodes$class)),
                   count=summary(nodes$class, maxsum = dim(nodes)[1])/sum(summary(nodes$class, maxsum = dim(nodes)[1]))*100, 
                   lineage=names(summary(nodes$class, maxsum = dim(nodes)[1])),
                   compartment=comp,
                   species=plant)
    
    #Merge
    plot=rbind(plot,sub)
    
  }
  
}

#set levels for plot
taxas=sort(unique(as.character(plot$lineage)))[-1]
plot$lineage=factor(plot$lineage, levels = c(taxas,"Uncharacterized_class","low_abundant","Other"))
plot[is.na(plot$lineage),"lineage"]="Uncharacterized_class"
plot$species=factor(plot$species, levels = plants[c(2,1,3)])

#select desired taxa
pro=sort(unique(taxon[taxon[,2]=="Proteobacteria",3]))[-1]
fir=sort(unique(taxon[taxon[,2]=="Firmicutes",3]))[-1]
act=sort(unique(taxon[taxon[,2]=="Actinobacteria",3]))
bac=sort(unique(taxon[taxon[,2]=="Bacteroidetes",3]))[-1]
asc=sort(unique(taxon[taxon[,2]=="Ascomycota",3]))[-1]
plot[!(plot$lineage %in% c(pro,fir,act,bac,asc)),"lineage"]="Other"

#set levels for taxa
print(length(unique(plot$lineage)))
plot$compartment=factor(plot$compartment, levels =compartment )


a=ggplot(plot, aes(x=species ,y=count, fill=lineage))+
  geom_bar(stat = "identity")+
  facet_grid(cols = vars(compartment) )+
  scale_fill_manual(values = c(getPalette(length(unique(plot$lineage))+1),"grey50"))+
  scale_x_discrete(labels=c("At","As","Ca"))+
  tema+ylab("Relative abundance %")

a

png(paste(plot.dir,paste(vert,"network_taxa","png",sep = "."),sep = "/"), 
    width = 1250, height = 600, units = "px", pointsize = 15, bg = "white", res=200, type="cairo")

print(a)

dev.off()


