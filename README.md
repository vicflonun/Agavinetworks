# Agavinetworks

This repository contains the code used to perform the downstream statistical analyses in Flores-Núñez, et al. 2021. Chapter 7. The Agave and Cacti microbiome: models for a planet under global warming. Metagenomics. Academic Press, in press
Data is reported in Fonseca-García, et al 2016 and Coleman-Derr, et al 2016. 

get_hubs.R: Generates correlation networks for groups of samples and calculates hub OTUs based on random networks
plot_taxa_networks.R: Generates taxa bar plots of hub OTUs

kruskal_season.R: Performs an enrichment analysis between seasons in each plant compartment and species
kruskal_species.R: Performs an enrichment analysis between plant compartments per each plants species
kruskal_sample.R: Performs an enrichment analysis between plant species per each compartment

venn_networks.R: Generates venn diagrams of shared hub, connected and central OTUs, and generates a heatmap showing their enrichment. 

plot_networks_soil.R: Generate network layouts for soil and frequency barplots
plot_networks_leaf.R: Generate network layouts for leaf compartments and frequency barplots
plot_networks_root.R: Generate network layouts for root compartments and frequency barplots

