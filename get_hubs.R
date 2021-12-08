##get_hubs
  #This code creates correlation networks of OTUs across differentof the agaves and cacti 
  #microbiome and identify consistent microbial hubs and higly central and connected taxa 
  #Data from (Coleman-Der et al, 2016;Fonseca-García, 2016; Partida-Martínez, 2018)

  #Extended explanation can be found in Flores-Núñez, et al. 2021. 

## Libraries ##
library(Hmisc)
library(igraph)
library(rgexf)
library(ggplot2)

## Directories ##
#WINDOWS
result_folder="github" 
setwd("C:/Users/victorfn/Google Drive/network_analysis_2018")
main.dir=getwd()
result.dir=paste(main.dir,result_folder,sep = "/")
dir.create(result.dir,showWarnings=F)
base.name= "CAM"

## Load data ##
load(file = paste(main.dir,paste(base.name,"OTU_table",sep = "_"),sep = "/"))
load(file = paste(main.dir,paste(base.name,"metadata_table",sep = "_"),sep = "/"))
load(file = paste(main.dir,paste(base.name,"taxonomy_table",sep = "_"),sep = "/"))

## Select samples to work with ##
##If you need to subsample based on site or season you just need to code the for-loop below

compartimentos=c("soil","roo.zone.soil","rhizosphere","root.endosphere","leaf.endosphere","phyllosphere")[c(3:6)]
plantas=c("Agave.salmiana","Agave.tequilana","Agave.deserti","Cacti")[c(1,2,4)]
#temp=c("dry","rainy")[2] #If you want to work with seasons #not possible for cacti
#site=c()

for (comp in compartimentos){
  for (plant in plantas){

    ##  Parameters ##
    
    #limit  #Interval of counts to consider when filtering the taxa and calculate the robust hubs. 
    #thr    #Abundance threshold. An OTU is considered if it has n counts are greater or equalto thr% of the samples 
    #per    #Proportion of time a hub OTU continuous to be a hub OTU in random networks. 
    if (comp=="leaf.endosphere" | comp=="root.endosphere"){limit=1 ; thr=0.50 ; per=50} else {limit=3 ; thr= 0.75 ; per=50}
    
    
## Construct the network ##

#Select samples
samples=rownames(meta[meta$Sample==comp & meta$Specie %in% plant,])
#samples=rownames(meta[meta$Sample==comp & meta$Specie==plant & meta$Season==temp,])
otus=otu[,samples]

#Filter by abundance and frequency #Rare OTUs generate are non-significant or generate spurious correlations  
otus=otus[rowSums(otus >= limit) >= length(samples)*thr,] 

#Generate correlations
cor=rcorr(t(otus), type = "spearman")
rom=cor$r
pvm=cor$P

#Filter data based on the correlation value (ro) 
rom[rom > -0.6 & rom < 0.6] <- 0

#Remove NaN
rom[is.nan(rom)] <- 0

#Filter correlations based on P-value
for (i in 1:dim(rom)[1]){
  logic=!(pvm[i,]<=0.01)
  rom[i, logic]<-0
}

#generate igraph object
g=graph.adjacency(rom, mode = "undirected", weighted = TRUE, diag = FALSE)

#delete isolated vertices
g=delete.vertices(g, V(g)[degree(g) == 0])

#Save network
save(g, file = paste(result.dir,paste(plant,comp,"igraph",sep = "."),sep = "/"))

#export to gexf format fo Gephi if needed
#g.gexf <- igraph.to.gexf(g)
#nombre=paste(plant,comp,"gexf",sep = ".")
#f <- file(paste(result.dir,nombre,sep="/"))
#writeLines(g.gexf$graph, con = f)
#close(f)
  
#Select the most connected OTUS
grado=degree(g, v = V(g), mode = c("all"),loops = F, normalized = FALSE)
grado=grado[grado>summary(grado)[5]]

#Select the most central OTUS
central=betweenness(g, v = V(g), directed = F, weights = abs(E(g)$weight),nobigint = F, normalized = FALSE)
central=central[central>summary(central)[5]]

#Select the hub OTUS
hub=names(grado[names(grado) %in% names(central)])

#Generate the result data.frames
grafica=data.frame(degree=degree(g, v = V(g), mode = c("all"),loops = F, normalized = FALSE),betweenesscentrality=betweenness(g, v = V(g), directed = F, weights = abs(E(g)$weight),nobigint = F, normalized = FALSE))
grafica=cbind(grafica,taxon[rownames(grafica),]) #add taxonomy
grafica$hub=0 ; grafica$central=0 ; grafica$grado=0 #add extra columns
grafica[hub,"hub"]=1 ; grafica[names(central),"central"]=1 ; grafica[names(grado),"grado"]=1 #indicate the hub, highly conncted and cetral OTUs



## RANDOM NETWORKS ##

#Generate the result data.frames
result=data.frame(row.names = rownames(grafica), 
                  PRS=rep(0, times=dim(grafica)[1]), #vertex present in random sample
                  HRS=rep(0, times=dim(grafica)[1]), #hub present in random sample
                  GRS=rep(0, times=dim(grafica)[1]), #highly connected vertex present in random sample
                  CRS=rep(0, times=dim(grafica)[1]), #highly central vertex present in random sample
                  NRS=rep(0, times=dim(grafica)[1])) #vertex non present in random sample

#Set seeds (every loop will subsample with a different seed)
initial.seed=1587571314
the.seed=71314

for (seed in seq(the.seed, initial.seed, 1588000 )){ 
  
  print((seed))
  random=data.frame(row.names = rownames(grafica), 
                    PRS=rep(0, times=dim(grafica)[1]),
                    HRS=rep(0, times=dim(grafica)[1]),
                    GRS=rep(0, times=dim(grafica)[1]),
                    CRS=rep(0, times=dim(grafica)[1]),
                    NRS=rep(0, times=dim(grafica)[1]))
#Random networks
otus2=otu[,samples]
set.seed(seed)
otus2=otus2[sample(x=rownames(otus2),size= 4316),]
rotus=otus2[rowSums(otus2 >= limit) >= length(samples)*thr,]

#Generate correlations
cor=rcorr(t(rotus), type = "spearman")
rom=cor$r
pvm=cor$P

#Filter data based on ro 
rom[rom > -0.6 & rom < 0.6] <- 0

#Remove NaN
rom[is.nan(rom)] <- 0

#Filter data based on P-value
for (i in 1:dim(rom)[1]){
  logic=!(pvm[i,]<=0.01)
  rom[i, logic]<-0
}

#generate igraph object
r.g=graph.adjacency(rom, mode = "undirected", weighted = TRUE, diag = FALSE)
#delate isolated vertices
r.g=delete.vertices(r.g, V(r.g)[degree(r.g) == 0])

#Select the most connected OTUS
r.grado=degree(r.g, v = V(r.g), mode = c("all"),loops = F, normalized = FALSE)
r.grado=r.grado[r.grado>summary(r.grado)[5]]

#Select the most central OTUS
r.central=betweenness(r.g, v = V(r.g), directed = F, weights = abs(E(r.g)$weight),nobigint = F, normalized = FALSE)
r.central=r.central[r.central>summary(r.central)[5]]

#Select the hub OTUS
r.hub=names(r.grado[names(r.grado) %in% names(r.central)])

#Determine the presence of the vertices
#PRS
random[rownames(random) %in% rownames(otus)[rownames(otus) %in% rownames(otus2)],"PRS"]=1
#HRS
random[rownames(random) %in% rownames(otus2)[rownames(otus2) %in% r.hub],"HRS"]=1
#GRS
random[rownames(random) %in% rownames(otus2)[rownames(otus2) %in% names(r.grado)],"GRS"]=1
#CRS
random[rownames(random) %in% rownames(otus2)[rownames(otus2) %in% names(r.central)],"CRS"]=1
#NRS
random[random$PRS==1 & random$HRS==0, "NRS"]=1

#Generate the counts 
result=result+random

}

#Calculate frequencies
result$HF=round(result$HRS/result$PRS*100,digits=1)
result$GF=round(result$GRS/result$PRS*100,digits=1)
result$CF=round(result$CRS/result$PRS*100,digits=1)

#Merge dataframes
grafica$random.hub=grafica$hub
grafica$random.grado=grafica$grado
grafica$random.central=grafica$central
sum(rownames(grafica)==rownames(result))==dim(grafica)[1]

#Indicate if a vertex is a hub, central or connected based on their frequencies
grafica[result$HF<per,"random.hub"]=0
grafica[result$GF<per,"random.grado"]=0
grafica[result$CF<per,"random.central"]=0

#Save results

#Each vertex with: centrality, degree, taxonomy, hub, central, connected (both in a single network and across random networks)
write.table(grafica, file = paste(result.dir,paste(plant,comp,"R_nodes","txt",sep = "."),sep = "/"),
            row.names = F, quote = F, sep = "\t" ,col.names = T)

#hub, connected and central frequencies in random networks
write.table(result, file = paste(result.dir,paste(plant,comp,"random_networks","txt",sep = "."),sep = "/"),
            row.names = T, quote = F, sep = "\t", col.names = T)

#Generate plots
freq=result
freq$OBH="NO HUB"
freq[hub,"OBH"]="HUB"
freq[,c("HRS","NRS")]=freq[,c("HRS","NRS")]/freq[,c("PRS")]*100
freq.hub=freq[hub,]
freq.nohub=freq[rownames(otus)[!(rownames(otus)%in%hub)],]

a=ggplot(data=freq, aes(x=HRS,  color=OBH))+
  geom_histogram(fill="white",alpha=0,position="identity", size=0.9)+
  scale_color_manual(values=c("orchid","forestgreen"))+
  ggtitle("Hub frequency in random networks (all otus)")+
  xlab("frequency of an OTU being a hub")+ylab("OTU Count")+
  geom_vline(xintercept = 40)


b=ggplot(data=freq.hub, aes(x=HRS))+
  geom_histogram(color="orchid",fill="white",alpha=.90,position="identity", size=0.9)+
  ggtitle("Hub frequency in random networks (only hub otus)")+
  xlab("frequency of a hub OTU being a hub")+ylab("OTU Count")+
  geom_vline(xintercept = 40)


c=ggplot(data=freq.nohub, aes(x=HRS))+
  geom_histogram(color="forestgreen",fill="white",alpha=.90,position="identity", size=0.9)+
  ggtitle("Hub frequency in random networks (only no hub otus)")+
  xlab("frequency of a no hub OTU being a hub")+ylab("OTU Count")+
  geom_vline(xintercept = 40)


pdf(file = paste(result.dir,paste(comp,plant,"hub_frequency.pdf",sep = "."),sep = "/"))
print(a);print(b);print(c)
dev.off()

  }
}


