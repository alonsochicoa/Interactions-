#This code is to obtain a phylogenetic relationship among species for Petry et al. 2018
#by using available information in phylomatic and phylocom
#User: Oscar Godoy Dec 10th 2017

setwd("/Users/oscargodoy/Documents/Post Doc California/AntFitnessSedgwick")

#first load brranching to create the species topology. 
library(brranching)
library(ape)

#load the speceis list and build the topology
focalspecieslist <- read.csv("data/focalspecies.csv")  # lookup table for updating taxonomy
taxa <- c(as.vector(focalspecieslist$species_name), "Amborella trichopoda", "Magnolia grandiflora")
#build the tree with available information from APGIII. For more info check http://phylodiversity.net/phylomatic/
topology <- phylomatic(taxa=taxa, get = 'POST')
# is the topology ok?
plot(topology, no.margin=TRUE)
write.tree(topology, file="data/phylo_topology.tre")


#This file that I am loading comes from downloading phylocom 4.2, and the ages included are according to Wikstrom et al. 2001 PRBS 
ages_df <- read.table("data/wikstrom.ages", sep="\t")
topology2 <- readLines("data/phylo_topology.tre")

#then load phylocomr to use bladj to calibrate the tree in million yr.
library(devtools) #to install phylocomr from git
install_github("ropensci/phylocomr") #this is in case needed. 
library(phylocomr)
phylo <- as.character(ph_bladj(ages = ages_df, phylo = topology2))
# eliminate "\n" after angiosperms 
gsub("\n\\b", "", phylo)
#write the phylogeny as a txt
write.table(phylo, file="data/phylo_calibrated.tre")

#prune the two species used to calibrate the tree that we do not longer need
#then store it. 
library(phytools)
phylo2 <- read.tree(file="data/phylo_calibrated.tre")
pruned.phylo <- drop.tip(phylo2, c("amborella_trichopoda", "magnolia_grandiflora"))
plot(pruned.phylo, no.margin=TRUE)
write.tree(pruned.phylo, file="data/phylo_calibrated.tre")