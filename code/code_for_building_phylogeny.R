#This code is to obtain a phylogenetic relationship among species for Alonso et al. 201X
#by using available information in phylomatic and phylocom
#User: Oscar Godoy and Alejandro Alonso, Dec 15th 2017

setwd("/Users/oscargodoy/Documents/Estudiantes/Alejandro Alonso US TFM/phylo_interactions_caracoles")

#first load brranching to create the species topology. 
library(brranching)
library(ape)
library(phytools)

#load the speceis list and build the topology
focalspecieslist <- read.csv("data/focalspecies.csv")  # lookup table for updating taxonomy
taxa <- c(as.vector(focalspecieslist$species), "Amborella trichopoda", "Magnolia grandiflora")
#build the tree with available information from APGIII. For more info check http://phylodiversity.net/phylomatic/
topology <- phylomatic(taxa=taxa, get = 'POST')
# is the topology ok?
plot(topology, no.margin=TRUE)
write.tree(topology, file="data/phylo_topology_polytomies.tre") # we have modified by hand the topology to eliminate polytomies in 
#Amaranthaceae


#This file that I am loading comes from downloading phylocom 4.2, and the ages included are according to Wikstrom et al. 2001 PRBS 
ages_df <- read.table("data/wikstrom.ages", sep="\t")
topology2 <- readLines("data/phylo_topology.tre")

#then load phylocomr to use bladj to calibrate the tree in million yr.
library(devtools) #to install phylocomr from git
install_github("ropensci/phylocomr") #this is in case needed. 
library(phylocomr)
phylo <- as.character(ph_bladj(ages = ages_df, phylo = topology2))
# eliminate "\n" after angiosperms 
phylo<- gsub("\n\\b", "", phylo)
#write the phylogeny as a txt
write.table(phylo, file="data/phylo_calibrated.tre")

#prune the two species used to calibrate the tree that we do not longer need
#then store it. 
library(phytools)
phylo2 <- read.tree(file="data/phylo_calibrated.tre")
pruned.phylo <- drop.tip(phylo2, c("amborella_trichopoda", "magnolia_grandiflora"))
plot(pruned.phylo, no.margin=TRUE)
write.tree(pruned.phylo, file="data/phylo_calibrated_without_outgroup.tre")
