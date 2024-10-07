library(tidyverse)
library(Morpho)
library(geomorph)
library(Rvcg)
library(paleomorph)
library(phytools)
library(rjson)
library(data.table)
library(geiger)
library(abind)

setwd("~/Ch_III_Thamno_Morphology/")

source("Analyses_Scripts/Supplementary_Functions/Matched_Local_Superimpositions.R")
meta_lms <- read.csv("Data/Landmarks_Bilat.csv")
Thamno_Tree <- drop.tip(read.tree("Data/Thamno_WGS1_82_ALRT_CCR_PLC.tre"), c("Natrix_natrix","Thamnophis_fulvus","Thamnophis_sirtalis_X_radix"))
Thamno_Tree <- force.ultrametric(Thamno_Tree)
load("Data/slid_coords.RData")

# Use Local Superimposition for Each Partition
lgpa.partition <- meta_lms_full$full_model #Full Skull
test <- local.gpa(coords = slid_coords$dataslide, partition = lgpa.partition) #Full Skull

# Read in Info on specimens and Taxa
scanned <- read.csv("Data/Specimens_Used.csv", header=T)

# Sort dataset by the order of specimens in the GPA
scanned <- scanned[match(dimnames(test)[[3]], scanned$catnum), ]

#--------------------------------------------------------------
# Group GPA coordinates by species and get species mean shapes
#--------------------------------------------------------------
# Create a 2D array that aggregate can work on
two.d.coords <- two.d.array(test)

# Aggregate by species
species.means <- aggregate(two.d.coords ~ scanned$species, FUN = mean)
rownames(species.means) <- species.means$`scanned$species`
species.means <- species.means[, -1]

test <- arrayspecs(species.means, dim(test)[1], dim(test)[2])

#Prune Tree to just taxa in array
All <- unlist(dimnames(test)[3])
All_Tree <- force.ultrametric(keep.tip(Thamno_Tree, All))
plot(ladderize(All_Tree))

# PCAs for all skull modules and bones
# Full Skull
Skull_Test <- test[1:1592,,]
Skull_Physig <- physignal(Skull_Test, All_Tree)
phypca_test <- gm.prcomp(Skull_Test, phy = All_Tree, GLS = T)
paca_test <- gm.prcomp(Skull_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_test$x, "Data/Skull_PhyPCA_PC.csv")
write.csv(paca_test$x, "Data/Skull_PACA_PC.csv")

#Premaxilla
Premax_Test <- test[1:76,,]
Premax_Physig <- physignal(Premax_Test, All_Tree)
phypca_premax <- gm.prcomp(Premax_Test, phy = All_Tree, GLS = T)
paca_premax <- gm.prcomp(Premax_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_premax$x, "Data/Premaxilla_PhyPCA_PC.csv")
write.csv(paca_premax$x, "Data/Premaxilla_PACA_PC.csv")

#Nasal
Nasal_Test <- test[77:152,,]
Nasal_Physig <- physignal(Nasal_Test, All_Tree)
phypca_Nasal <- gm.prcomp(Nasal_Test, phy = All_Tree, GLS = T)
paca_Nasal <- gm.prcomp(Nasal_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Nasal$x, "Data/Nasal_PhyPCA_PC.csv")
write.csv(paca_Nasal$x, "Data/Nasal_PACA_PC.csv")

#Frontal
Frontal_Test <- test[153:228,,]
Frontal_Physig <- physignal(Frontal_Test, All_Tree)
phypca_Frontal <- gm.prcomp(Frontal_Test, phy = All_Tree, GLS = T)
paca_Frontal <- gm.prcomp(Frontal_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Frontal$x, "Data/Frontal_PhyPCA_PC.csv")
write.csv(paca_Frontal$x, "Data/Frontal_PACA_PC.csv")

#Parietal
Parietal_Test <- test[229:304,,]
Parietal_Physig <- physignal(Parietal_Test, All_Tree)
phypca_Parietal <- gm.prcomp(Parietal_Test, phy = All_Tree, GLS = T)
paca_Parietal <- gm.prcomp(Parietal_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Parietal$x, "Data/Parietal_PhyPCA_PC.csv")
write.csv(paca_Parietal$x, "Data/Parietal_PACA_PC.csv")

#Supraoccipital
Supraocc_Test <- test[305:402,,]
Supraocc_Physig <- physignal(Supraocc_Test, All_Tree)
phypca_Supraocc <- gm.prcomp(Supraocc_Test, phy = All_Tree, GLS = T)
paca_Supraocc <- gm.prcomp(Supraocc_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Supraocc$x, "Data/Supraoccipital_PhyPCA_PC.csv")
write.csv(paca_Supraocc$x, "Data/Supraoccipital_PACA_PC.csv")

#Basisphenoid
Basi_Test <- test[403:478,,]
Basi_Physig <- physignal(Basi_Test, All_Tree)
phypca_Basi <- gm.prcomp(Basi_Test, phy = All_Tree, GLS = T)
paca_Basi <- gm.prcomp(Basi_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Basi$x, "Data/Basisphenoid_PhyPCA_PC.csv")
write.csv(paca_Basi$x, "Data/Basisphenoid_PACA_PC.csv")

#Maxilla
Maxilla_Test <- test[c(479:484,491:560),,] #Left Only
Maxilla_Physig <- physignal(Maxilla_Test, All_Tree)
phypca_Maxilla <- gm.prcomp(Maxilla_Test, phy = All_Tree, GLS = T)
paca_Maxilla <- gm.prcomp(Maxilla_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Maxilla$x, "Data/Maxilla_PhyPCA_PC.csv")
write.csv(paca_Maxilla$x, "Data/Maxilla_PACA_PC.csv")

#Palatine
Palatine_Test <- test[c(631:640,651:740),,]
Palatine_Physig <- physignal(Palatine_Test, All_Tree)
phypca_Palatine <- gm.prcomp(Palatine_Test, phy = All_Tree, GLS = T)
paca_Palatine <- gm.prcomp(Palatine_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Palatine$x, "Data/Palatine_PhyPCA_PC.csv")
write.csv(paca_Palatine$x, "Data/Palatine_PACA_PC.csv")

#Ectopterygoid
Ecto_Test <- test[c(831:833,837:866),,]
Ecto_Physig <- physignal(Ecto_Test, All_Tree)
phypca_Ecto <- gm.prcomp(Ecto_Test, phy = All_Tree, GLS = T)
paca_Ecto <- gm.prcomp(Ecto_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Ecto$x, "Data/Ectopterygoid_PhyPCA_PC.csv")
write.csv(paca_Ecto$x, "Data/Ectopterygoid_PACA_PC.csv")

#Pterygoid
Pterygoid_Test <- test[c(897:901,907:976),,]
Pterygoid_Physig <- physignal(Pterygoid_Test, All_Tree)
phypca_Pterygoid <- gm.prcomp(Pterygoid_Test, phy = All_Tree, GLS = T)
paca_Pterygoid <- gm.prcomp(Pterygoid_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Pterygoid$x, "Data/Pterygoid_PhyPCA_PC.csv")
write.csv(paca_Pterygoid$x, "Data/Pterygoid_PACA_PC.csv")

#Dentary
Dentary_Test <- test[c(1047:1053,1061:1120),,]
Dentary_Physig <- physignal(Dentary_Test, All_Tree)
phypca_Dentary <- gm.prcomp(Dentary_Test, phy = All_Tree, GLS = T)
paca_Dentary <- gm.prcomp(Dentary_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Dentary$x, "Data/Dentary_PhyPCA_PC.csv")
write.csv(paca_Dentary$x, "Data/Dentary_PACA_PC.csv")

#Compound
Compound_Test <- test[c(1181:1192,1205:1264),,]
Compound_Physig <- physignal(Compound_Test, All_Tree)
phypca_Compound <- gm.prcomp(Compound_Test, phy = All_Tree, GLS = T)
paca_Compound <- gm.prcomp(Compound_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Compound$x, "Data/Compound_PhyPCA_PC.csv")
write.csv(paca_Compound$x, "Data/Compound_PACA_PC.csv")

#Quadrate
Quadrate_Test <- test[c(1325:1336,1349:1448),,]
Quadrate_Physig <- physignal(Quadrate_Test, All_Tree)
phypca_Quadrate <- gm.prcomp(Quadrate_Test, phy = All_Tree, GLS = T)
paca_Quadrate <- gm.prcomp(Quadrate_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Quadrate$x, "Data/Quadrate_PhyPCA_PC.csv")
write.csv(paca_Quadrate$x, "Data/Quadrate_PACA_PC.csv")

#Supratemporal
Supratemp_Test <- test[c(1549:1550,1553:1572),,]
Supratemp_Physig <- physignal(Supratemp_Test, All_Tree)
phypca_Supratemp <- gm.prcomp(Supratemp_Test, phy = All_Tree, GLS = T)
paca_Supratemp <- gm.prcomp(Supratemp_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Supratemp$x, "Data/Supratemporal_PhyPCA_PC.csv")
write.csv(paca_Supratemp$x, "Data/Supratemporal_PACA_PC.csv")

#Snout
Snout_Test <- test[1:152,,]
Snout_Physig <- physignal(Snout_Test, All_Tree)
phypca_Snout <- gm.prcomp(Snout_Test, phy = All_Tree, GLS = T)
paca_Snout <- gm.prcomp(Snout_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Snout$x, "Data/Snout_PhyPCA_PC.csv")
write.csv(paca_Snout$x, "Data/Snout_PACA_PC.csv")

#Braincase
Braincase_Test <- test[153:478,,]
Braincase_Physig <- physignal(Braincase_Test, All_Tree)
phypca_Braincase <- gm.prcomp(Braincase_Test, phy = All_Tree, GLS = T)
paca_Braincase <- gm.prcomp(Braincase_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Braincase$x, "Data/Braincase_PhyPCA_PC.csv")
write.csv(paca_Braincase$x, "Data/Braincase_PACA_PC.csv")

#Palatopterygoid Arch
Palatoptery_Test <- test[c(631:640,651:740,831:833,837:866,897:901,907:976),,]
Palatoptery_Physig <- physignal(Palatoptery_Test, All_Tree)
phypca_Palatoptery <- gm.prcomp(Palatoptery_Test, phy = All_Tree, GLS = T)
paca_Palatoptery <- gm.prcomp(Palatoptery_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Palatoptery$x, "Data/Palatopterygoid_Arch_PhyPCA_PC.csv")
write.csv(paca_Palatoptery$x, "Data/Palatopterygoid_Arch_PACA_PC.csv")

#Mandible
Mandible_Test <- test[c(1047:1053,1061:1120,1181:1192,1205:1264),,]
Mandible_Physig <- physignal(Mandible_Test, All_Tree)
phypca_Mandible <- gm.prcomp(Mandible_Test, phy = All_Tree, GLS = T)
paca_Mandible <- gm.prcomp(Mandible_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Mandible$x, "Data/Mandible_PhyPCA_PC.csv")
write.csv(paca_Mandible$x, "Data/Mandible_PACA_PC.csv")

#Suspensorium
Suspensorium_Test <- test[c(1325:1336,1349:1448,1549:1550,1553:1572),,]
Suspensorium_Physig <- physignal(Suspensorium_Test, All_Tree)
phypca_Suspensorium <- gm.prcomp(Suspensorium_Test, phy = All_Tree, GLS = T)
paca_Suspensorium <- gm.prcomp(Suspensorium_Test, phy = All_Tree, align.to.phy = T)
write.csv(phypca_Suspensorium$x, "Data/Suspensorium_PhyPCA_PC.csv")
write.csv(paca_Suspensorium$x, "Data/Suspensorium_PACA_PC.csv")
