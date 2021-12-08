#Performs an enrichment analysis between plant species per each compartment
#Please indicate below the Domain you want to analyze and the taxa rank

#Directories
setwd("C:/Users/victorfn/Google Drive/taxa_enrichment_2020/")
result="C:/Users/victorfn/Google Drive/taxa_enrichment_2020/github"
dir.create(result)
raw.data="C:/Users/victorfn/Google Drive/Microbiome_base"

#Packages
library("FSA")
library("Rmisc")
library("tidyverse")
library("phyloseq")
library("pheatmap")
library("colorspace")
library("ggrepel")

# Load data 
load(file = paste(getwd(),"CAM_metadata_table", sep = "/"))
load(file = paste(getwd(),"CAM_OTU_table",sep = "/"))
load(file = paste(getwd(),"CAM_taxonomy_table",sep = "/"))

#Merge cacti if needed
meta[meta$Specie %in% c("O.robusta","M.geometrizans"),"Specie"]=c("Cacti")

#Merge data in PHYLOSEQ object
cam=phyloseq(otu_table(otu, taxa_are_rows = T),tax_table(taxon),sample_data(meta))

# Subset your data by compartment and lineage using logic formulas
microbe=c("Bacteria","Fungi")[c(3)]
#Chose 3 compartments only 
compartment=c("soil","roo.zone.soil","root.endosphere","rhizosphere","leaf.endosphere","phyllosphere")[c(1,3:6)]
plant=c("Agave.tequilana","Agave.deserti","Agave.salmiana","O.robusta","M.geometrizans","Cacti")[c(1,3,6)]
#temp=c("dry","rainy")
taxa.names=rank_names(cam)[c(5)]
#grupo="endo"
#thr=0.05 

#Filter data set and generate relative abundance
cam0=subset_taxa(cam, kingdom %in% microbe) %>%
  transform_sample_counts(function(OTU) OTU/sum(OTU)*100) # ; sample_sums(cam0)
#core(detection = 0, prevalence = 0.75) #for core usage


for (m in taxa.names){
  #Sums or colapse the lineages
  if (m=="OTU.ID"){print(m)
    cam1=cam0
  } else if (m!="OTU.ID"){ 
    print(m)
    cam1=tax_glom(cam0,m, bad_empty = c(NA, "", " ", "\t",NaN, "SJA","PRR","o","f","g","c","p","4CL"))
  }


for (i in compartment){
  
  taxas=get_taxa_unique(cam1, m) #creates a vector for the lineages
  result.frame=data.frame(row.names = taxas) #generates result data frame of pvalues
    
    cam2=subset_samples(cam1, Sample %in% i) %>%
      subset_samples(Specie %in% plant) %>% #Filter the data by plant and compartment
      filter_taxa(function(x) sum(x > 0) > (0.50*length(x)), TRUE) #Prune 0% abundante in 75% of samples
    
    sample.names=sample_names(cam2) #Extract sample names
    
    #Extract phyloseq dataframes
    otu.frame=data.frame(otu_table(cam2), lineage = get_taxa_unique(cam2, m)) 
    
    #Tranform to tidy data
    otu.tidy=gather(otu.frame, key = "library", value = "relative.abundance", sample.names)
    
    for (h in sample.names){ #Add metadata
      otu.tidy[otu.tidy$library==h,c("Sample", "Season", "Specie", "Location")]=meta[h,c("Sample", "Season", "Specie", "Location")]
    } 
    
    
    ### Kruskal TETS ### 
    otu.tidy$lineage=factor(otu.tidy$lineage, levels = taxas) #Change strings to factors 
    otu.tidy$Specie=factor(otu.tidy$Specie) #Change strings to factors 
    kruskal=data.frame(row.names = taxas,p.value=rep(1,times=length(taxas))) #temporal data frame
    for (h in as.character(unique(otu.frame$lineage))){
      sub1=otu.tidy[otu.tidy$lineage==h,]
      kruskal[h,"p.value"]=kruskal.test(relative.abundance ~ Specie, data = sub1)$p.value
    }
    ### 
    
    ### Dunn for multiple comparison
    if (sum(kruskal$p.value<=0.05)>0){
      dunn=data.frame(row.names = taxas)
      for (h in as.character(unique(otu.frame$lineage))){
        sub1=otu.tidy[otu.tidy$lineage==h,]
        test1=dunnTest(relative.abundance ~ Specie, data = sub1, method = "bh")
        p.value=test1$res$P.adj
        dunn[h,1:choose(length(plant),2)]=p.value
      }
      colnames(dunn)=c(as.character(test1$res$Comparison[1:2]), paste(plant[3],plant[1],sep = " - "))
      dunn[is.na(dunn)]=1
      
      
      write.table(dunn,paste(result,paste(microbe,m,"Pvalues",i,"txt",sep="."),sep = "/")[1], 
                  quote = F, row.names = T, col.names = T, sep = "\t")
      
      
    }else{print(paste("No differentially abundant",m,"in", j))}
    ###
    
    
    
    ### Calcultaes the dirrection of the enrichment ###
    log.frame=-log10(dunn)
    otu.se=summarySE(data=otu.tidy, measurevar = "relative.abundance", groupvars = c("Specie","lineage"))
    otu.dif=data.frame(row.names=unique(otu.se$lineage),
                       V1=otu.se[otu.se$Specie==plant[2],"relative.abundance"]-otu.se[otu.se$Specie==plant[1],"relative.abundance"],
                       V2=otu.se[otu.se$Specie==plant[2],"relative.abundance"]-otu.se[otu.se$Specie==plant[3],"relative.abundance"],
                       V1=otu.se[otu.se$Specie==plant[3],"relative.abundance"]-otu.se[otu.se$Specie==plant[1],"relative.abundance"])
    otu.dif[otu.dif>0]=1 
    otu.dif[otu.dif<0]=-1 
    colnames(otu.dif)=c(as.character(test1$res$Comparison[1:2]), paste(plant[3],plant[1],sep = " - "))
    
    for (h in as.character(unique(otu.frame$lineage))){ #merges de direction of the enrichment
      log.frame[h,]= log.frame[h,]*otu.dif[h,] 
      }

    write.table(log.frame,paste(result,paste(microbe,m,"logPvalues",i,"txt",sep="."),sep = "/")[1], 
                quote = F, row.names = T, col.names = T, sep = "\t")
    
    
}

}
