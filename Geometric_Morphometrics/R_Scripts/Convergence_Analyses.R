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
library(treedata.table)
library(RRphylo)


setwd("~/Ch_III_Thamno_Morphology/")

# Read in tree
Thamno_Tree <- drop.tip(read.tree("Geometric_Morphometrics/Data/Thamno_WGS1_82_ALRT_CCR_PLC.tre"), c("Natrix_natrix","Thamnophis_fulvus","Thamnophis_sirtalis_X_radix"))
Thamno_Tree <- force.ultrametric(Thamno_Tree)

# Read in PC Scores
paca_skull <- as.matrix(read.csv("Geometric_Morphometrics/Results/Skull_PACA_PC.csv", row.names = 1))[,1:4]
paca_snout <- as.matrix(read.csv("Geometric_Morphometrics/Results/Snout_PACA_PC.csv", row.names = 1))[,1:4]
paca_braincase <- as.matrix(read.csv("Geometric_Morphometrics/Results/Braincase_PACA_PC.csv", row.names = 1))[,1:4]
paca_maxilla <- as.matrix(read.csv("Geometric_Morphometrics/Results/Maxilla_PACA_PC.csv", row.names = 1))[,1:4]
paca_palatoptery <- as.matrix(read.csv("Geometric_Morphometrics/Results/Palatopterygoid_Arch_PACA_PC.csv", row.names = 1))[,1:4]
paca_mandible <- as.matrix(read.csv("Geometric_Morphometrics/Results/Mandible_PACA_PC.csv", row.names = 1))[,1:4]
paca_suspensorium <- as.matrix(read.csv("Geometric_Morphometrics/Results/Suspensorium_PACA_PC.csv", row.names = 1))[,1:4]

# Read in Trait Data for Convergence Analyses
traits <- read.csv("Geometric_Morphometrics/Data/All_Thamno_Traits.csv")
Thamno <- as.treedata.table(Thamno_Tree, traits, name_column = "Name_in_Tree")

# Extract Trait Data from treedata.table
diet <- extractVector(Thamno, "Diet")
habitat <- extractVector(Thamno, "Lifestyle")
clade <- extractVector(Thamno, "Clade")
genus <- extractVector(Thamno, "Genus")

# Convergence Tests using RRphylo
Partitions <- c("Skull", "Snout", "Braincase", "Palatomaxillary_Arch", "Palatopterygoid_Arch", "Mandible", "Suspensorium", "Lower_Jaw", "Feeding_Apparatus", "Premaxilla", "Nasal", "Frontal", "Parietal", "Supraoccipital", "Basisphenoid", "Maxilla", "Palatine", "Ectopterygoid", "Pterygoid", "Dentary", "Compound", "Quadrate", "Supratemporal")
Reduced_Partitions <- c("Skull", "Snout", "Braincase", "Maxilla", "Palatopterygoid_Arch", "Mandible", "Suspensorium")

# Habitat - Aquatic
aquatic_vector <- rep("nostate", 51)
names(aquatic_vector)<-rownames(paca_skull)
Aquatic_Taxa <- which(habitat=="Aquatic")
Aquatic_Taxa <- names(Aquatic_Taxa)
aquatic_vector[which(names(aquatic_vector)%in%Aquatic_Taxa)] <- "aquatic"
aquatic_conv_skull.test <- search.conv(tree=Thamno_Tree, y=paca_skull,  state = aquatic_vector, declust = T)
aquatic_conv_snout.test <- search.conv(tree=Thamno_Tree, y=paca_snout,  state = aquatic_vector, declust = T)
aquatic_conv_braincase.test <- search.conv(tree=Thamno_Tree, y=paca_braincase,  state = aquatic_vector, declust = T)
aquatic_conv_maxilla.test <- search.conv(tree=Thamno_Tree, y=paca_maxilla,  state = aquatic_vector, declust = T)
aquatic_conv_palatoptery.test <- search.conv(tree=Thamno_Tree, y=paca_palatoptery,  state = aquatic_vector, declust = T)
aquatic_conv_mandible.test <- search.conv(tree=Thamno_Tree, y=paca_mandible,  state = aquatic_vector, declust = T)
aquatic_conv_suspensorium.test <- search.conv(tree=Thamno_Tree, y=paca_suspensorium,  state = aquatic_vector, declust = T)

# Habitat - Terrestrial
terrestrial_vector <- rep("nostate", 51)
names(terrestrial_vector)<-rownames(paca_skull)
Terrestrial_Taxa <- which(habitat=="Terrestrial")
Terrestrial_Taxa <- names(Terrestrial_Taxa)
terrestrial_vector[which(names(terrestrial_vector)%in%Terrestrial_Taxa)] <- "terrestrial"
terrestrial_conv_skull.test <- search.conv(tree=Thamno_Tree, y=paca_skull,  state = terrestrial_vector, declust = T)
Terrestrial_conv_snout.test <- search.conv(tree=Thamno_Tree, y=paca_snout,  state = terrestrial_vector, declust = T)
terrestrial_conv_braincase.test <- search.conv(tree=Thamno_Tree, y=paca_braincase,  state = terrestrial_vector, declust = T)
terrestrial_conv_maxilla.test <- search.conv(tree=Thamno_Tree, y=paca_maxilla,  state = terrestrial_vector, declust = T)
terrestrial_conv_palatoptery.test <- search.conv(tree=Thamno_Tree, y=paca_palatoptery,  state = terrestrial_vector, declust = T)
terrestrial_conv_mandible.test <- search.conv(tree=Thamno_Tree, y=paca_mandible,  state = terrestrial_vector, declust = T)
terrestrial_conv_suspensorium.test <- search.conv(tree=Thamno_Tree, y=paca_suspensorium,  state = terrestrial_vector, declust = T)

# Habitat - Fossorial
Fossorial_vector <- rep("nostate", 51)
names(Fossorial_vector)<-rownames(paca_skull)
Fossorial_Taxa <- which(habitat=="Fossorial")
Fossorial_Taxa <- names(Fossorial_Taxa)
Fossorial_vector[which(names(Fossorial_vector)%in%Fossorial_Taxa)] <- "Fossorial"
Fossorial_conv_skull.test <- search.conv(tree=Thamno_Tree, y=paca_skull,  state = Fossorial_vector, declust = T)
Fossorial_conv_snout.test <- search.conv(tree=Thamno_Tree, y=paca_snout,  state = Fossorial_vector, declust = T)
Fossorial_conv_braincase.test <- search.conv(tree=Thamno_Tree, y=paca_braincase,  state = Fossorial_vector, declust = T)
Fossorial_conv_maxilla.test <- search.conv(tree=Thamno_Tree, y=paca_maxilla,  state = Fossorial_vector, declust = T)
Fossorial_conv_palatoptery.test <- search.conv(tree=Thamno_Tree, y=paca_palatoptery,  state = Fossorial_vector, declust = T)
Fossorial_conv_mandible.test <- search.conv(tree=Thamno_Tree, y=paca_mandible,  state = Fossorial_vector, declust = T)
Fossorial_conv_suspensorium.test <- search.conv(tree=Thamno_Tree, y=paca_suspensorium,  state = Fossorial_vector, declust = T)

# Diet - Generalist
Generalist_vector <- rep("nostate", 51)
names(Generalist_vector)<-rownames(paca_skull)
Generalist_Taxa <- which(diet=="Generalist")
Generalist_Taxa <- names(Generalist_Taxa)
Generalist_vector[which(names(Generalist_vector)%in%Generalist_Taxa)] <- "Generalist"
Generalist_conv_skull.test <- search.conv(tree=Thamno_Tree, y=paca_skull,  state = Generalist_vector, declust = T)
Generalist_conv_snout.test <- search.conv(tree=Thamno_Tree, y=paca_snout,  state = Generalist_vector, declust = T)
Generalist_conv_braincase.test <- search.conv(tree=Thamno_Tree, y=paca_braincase,  state = Generalist_vector, declust = T)
Generalist_conv_maxilla.test <- search.conv(tree=Thamno_Tree, y=paca_maxilla,  state = Generalist_vector, declust = T)
Generalist_conv_palatoptery.test <- search.conv(tree=Thamno_Tree, y=paca_palatoptery,  state = Generalist_vector, declust = T)
Generalist_conv_mandible.test <- search.conv(tree=Thamno_Tree, y=paca_mandible,  state = Generalist_vector, declust = T)
Generalist_conv_suspensorium.test <- search.conv(tree=Thamno_Tree, y=paca_suspensorium,  state = Generalist_vector, declust = T)

# Diet - Aquatic Generalist
Aquatic_Generalist_vector <- rep("nostate", 51)
names(Aquatic_Generalist_vector)<-rownames(paca_skull)
Aquatic_Generalist_Taxa <- which(diet=="Aquatic_Generalist")
Aquatic_Generalist_Taxa <- names(Aquatic_Generalist_Taxa)
Aquatic_Generalist_vector[which(names(Aquatic_Generalist_vector)%in%Aquatic_Generalist_Taxa)] <- "Aquatic_Generalist"
Aquatic_Generalist_conv_skull.test <- search.conv(tree=Thamno_Tree, y=paca_skull,  state = Aquatic_Generalist_vector, declust = T)
Aquatic_Generalist_conv_snout.test <- search.conv(tree=Thamno_Tree, y=paca_snout,  state = Aquatic_Generalist_vector, declust = T)
Aquatic_Generalist_conv_braincase.test <- search.conv(tree=Thamno_Tree, y=paca_braincase,  state = Aquatic_Generalist_vector, declust = T)
Aquatic_Generalist_conv_maxilla.test <- search.conv(tree=Thamno_Tree, y=paca_maxilla,  state = Aquatic_Generalist_vector, declust = T)
Aquatic_Generalist_conv_palatoptery.test <- search.conv(tree=Thamno_Tree, y=paca_palatoptery,  state = Aquatic_Generalist_vector, declust = T)
Aquatic_Generalist_conv_mandible.test <- search.conv(tree=Thamno_Tree, y=paca_mandible,  state = Aquatic_Generalist_vector, declust = T)
Aquatic_Generalist_conv_suspensorium.test <- search.conv(tree=Thamno_Tree, y=paca_suspensorium,  state = Aquatic_Generalist_vector, declust = T)

# Diet - Vermivorous
Vermivorous_vector <- rep("nostate", 51)
names(Vermivorous_vector)<-rownames(paca_skull)
Vermivorous_Taxa <- which(diet=="Vermivorous")
Vermivorous_Taxa <- names(Vermivorous_Taxa)
Vermivorous_vector[which(names(Vermivorous_vector)%in%Vermivorous_Taxa)] <- "Vermivorous"
Vermivorous_conv_skull.test <- search.conv(tree=Thamno_Tree, y=paca_skull,  state = Vermivorous_vector, declust = T)
Vermivorous_conv_snout.test <- search.conv(tree=Thamno_Tree, y=paca_snout,  state = Vermivorous_vector, declust = T)
Vermivorous_conv_braincase.test <- search.conv(tree=Thamno_Tree, y=paca_braincase,  state = Vermivorous_vector, declust = T)
Vermivorous_conv_maxilla.test <- search.conv(tree=Thamno_Tree, y=paca_maxilla,  state = Vermivorous_vector, declust = T)
Vermivorous_conv_palatoptery.test <- search.conv(tree=Thamno_Tree, y=paca_palatoptery,  state = Vermivorous_vector, declust = T)
Vermivorous_conv_mandible.test <- search.conv(tree=Thamno_Tree, y=paca_mandible,  state = Vermivorous_vector, declust = T)
Vermivorous_conv_suspensorium.test <- search.conv(tree=Thamno_Tree, y=paca_suspensorium,  state = Vermivorous_vector, declust = T)

# Diet - Durophagous
Durophagous_vector <- rep("nostate", 51)
names(Durophagous_vector)<-rownames(paca_skull)
Durophagous_Taxa <- which(diet=="Durophagous")
Durophagous_Taxa <- names(Durophagous_Taxa)
Durophagous_vector[which(names(Durophagous_vector)%in%Durophagous_Taxa)] <- "Durophagous"
Durophagous_conv_skull.test <- search.conv(tree=Thamno_Tree, y=paca_skull,  state = Durophagous_vector, declust = T)
Durophagous_conv_snout.test <- search.conv(tree=Thamno_Tree, y=paca_snout,  state = Durophagous_vector, declust = T)
Durophagous_conv_braincase.test <- search.conv(tree=Thamno_Tree, y=paca_braincase,  state = Durophagous_vector, declust = T)
Durophagous_conv_maxilla.test <- search.conv(tree=Thamno_Tree, y=paca_maxilla,  state = Durophagous_vector, declust = T)
Durophagous_conv_palatoptery.test <- search.conv(tree=Thamno_Tree, y=paca_palatoptery,  state = Durophagous_vector, declust = T)
Durophagous_conv_mandible.test <- search.conv(tree=Thamno_Tree, y=paca_mandible,  state = Durophagous_vector, declust = T)
Durophagous_conv_suspensorium.test <- search.conv(tree=Thamno_Tree, y=paca_suspensorium,  state = Durophagous_vector, declust = T)


# Plot Convergence Test Results
# Habitat - Aquatic
# Skull
aquatic_conv.plot <- plotConv(SC=aquatic_conv_skull.test,y=paca_skull,variable=1,state=aquatic_vector)
aquatic_conv.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#3B9AB2"),lty=3),
                             points.args=list(pch=c(23,21),bg=c("#D9D0D3","#3B9AB2")),
                             legend.args=list(pch=c(23,21),x="top"))
aquatic_conv.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                            polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                            line.args=list(line.col="#3B9AB2",lwd=4))
# Snout
aquatic_conv_Snout.plot <- plotConv(SC=aquatic_conv_snout.test,y=paca_snout,variable=1,state=aquatic_vector)
aquatic_conv_Snout.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#3B9AB2"),lty=3),
                                   points.args=list(pch=c(23,21),bg=c("#D9D0D3","#3B9AB2")),
                                   legend.args=list(pch=c(23,21),x="top"))
aquatic_conv_Snout.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                  polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                  line.args=list(line.col="#3B9AB2",lwd=4))
# Braincase
aquatic_conv_braincase.plot <- plotConv(SC=aquatic_conv_braincase.test,y=paca_braincase,variable=1,state=aquatic_vector)
aquatic_conv_braincase.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#3B9AB2"),lty=3),
                                       points.args=list(pch=c(23,21),bg=c("#D9D0D3","#3B9AB2")),
                                       legend.args=list(pch=c(23,21),x="top"))
aquatic_conv_braincase.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                      polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                      line.args=list(line.col="#3B9AB2",lwd=4))
# Maxilla
aquatic_conv_maxilla.plot <- plotConv(SC=aquatic_conv_maxilla.test,y=paca_maxilla,variable=1,state=aquatic_vector)
aquatic_conv_maxilla.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#3B9AB2"),lty=3),
                                     points.args=list(pch=c(23,21),bg=c("#D9D0D3","#3B9AB2")),
                                     legend.args=list(pch=c(23,21),x="top"))
aquatic_conv_maxilla.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                    polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                    line.args=list(line.col="#3B9AB2",lwd=4))
# Palatopterygoid Arch
aquatic_conv_palatoptery.plot <- plotConv(SC=aquatic_conv_palatoptery.test,y=paca_palatoptery,variable=1,state=aquatic_vector)
aquatic_conv_palatoptery.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#3B9AB2"),lty=3),
                                         points.args=list(pch=c(23,21),bg=c("#D9D0D3","#3B9AB2")),
                                         legend.args=list(pch=c(23,21),x="top"))
aquatic_conv_palatoptery.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                        polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                        line.args=list(line.col="#3B9AB2",lwd=4))
# Mandible
aquatic_conv_mandible.plot <- plotConv(SC=aquatic_conv_mandible.test,y=paca_mandible,variable=1,state=aquatic_vector)
aquatic_conv_mandible.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#3B9AB2"),lty=3),
                                      points.args=list(pch=c(23,21),bg=c("#D9D0D3","#3B9AB2")),
                                      legend.args=list(pch=c(23,21),x="top"))
aquatic_conv_mandible.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                     polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                     line.args=list(line.col="#3B9AB2",lwd=4))
# Suspensorium
aquatic_conv_suspensorium.plot <- plotConv(SC=aquatic_conv_suspensorium.test,y=paca_suspensorium,variable=1,state=aquatic_vector)
aquatic_conv_suspensorium.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#3B9AB2"),lty=3),
                                          points.args=list(pch=c(23,21),bg=c("#D9D0D3","#3B9AB2")),
                                          legend.args=list(pch=c(23,21),x="top"))
aquatic_conv_suspensorium.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                         polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                         line.args=list(line.col="#3B9AB2",lwd=4))

# Habitat - Terrestrial
# Skull
terrestrial_conv.plot <- plotConv(SC=terrestrial_conv_skull.test,y=paca_skull,variable=1,state=terrestrial_vector)
terrestrial_conv.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F21A00"),lty=3),
                                 points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F21A00")),
                                 legend.args=list(pch=c(23,21),x="top"))
terrestrial_conv.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                line.args=list(line.col="#F21A00",lwd=4))
# Snout
terrestrial_conv_Snout.plot <- plotConv(SC=Terrestrial_conv_snout.test,y=paca_snout,variable=1,state=terrestrial_vector)
terrestrial_conv_Snout.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F21A00"),lty=3),
                                       points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F21A00")),
                                       legend.args=list(pch=c(23,21),x="top"))
terrestrial_conv_Snout.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                      polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                      line.args=list(line.col="#F21A00",lwd=4))
# Braincase
terrestrial_conv_braincase.plot <- plotConv(SC=terrestrial_conv_braincase.test,y=paca_braincase,variable=1,state=terrestrial_vector)
terrestrial_conv_braincase.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F21A00"),lty=3),
                                           points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F21A00")),
                                           legend.args=list(pch=c(23,21),x="top"))
terrestrial_conv_braincase.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                          polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                          line.args=list(line.col="#F21A00",lwd=4))
# Maxilla
terrestrial_conv_maxilla.plot <- plotConv(SC=terrestrial_conv_maxilla.test,y=paca_maxilla,variable=1,state=terrestrial_vector)
terrestrial_conv_maxilla.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F21A00"),lty=3),
                                         points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F21A00")),
                                         legend.args=list(pch=c(23,21),x="top"))
terrestrial_conv_maxilla.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                        polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                        line.args=list(line.col="#F21A00",lwd=4))
# Palatopterygoid Arch
terrestrial_conv_palatoptery.plot <- plotConv(SC=terrestrial_conv_palatoptery.test,y=paca_palatoptery,variable=1,state=terrestrial_vector)
terrestrial_conv_palatoptery.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F21A00"),lty=3),
                                             points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F21A00")),
                                             legend.args=list(pch=c(23,21),x="top"))
terrestrial_conv_palatoptery.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                            polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                            line.args=list(line.col="#F21A00",lwd=4))
# Mandible
terrestrial_conv_mandible.plot <- plotConv(SC=terrestrial_conv_mandible.test,y=paca_mandible,variable=1,state=terrestrial_vector)
terrestrial_conv_mandible.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F21A00"),lty=3),
                                          points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F21A00")),
                                          legend.args=list(pch=c(23,21),x="top"))
terrestrial_conv_mandible.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                         polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                         line.args=list(line.col="#F21A00",lwd=4))
# Suspensorium
terrestrial_conv_suspensorium.plot <- plotConv(SC=terrestrial_conv_suspensorium.test,y=paca_suspensorium,variable=1,state=terrestrial_vector)
terrestrial_conv_suspensorium.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F21A00"),lty=3),
                                              points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F21A00")),
                                              legend.args=list(pch=c(23,21),x="top"))
terrestrial_conv_suspensorium.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                             polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                             line.args=list(line.col="#F21A00",lwd=4))
# Habitat - Fossorial
# Skull
Fossorial_conv.plot <- plotConv(SC=Fossorial_conv_skull.test,y=paca_skull,variable=1,state=Fossorial_vector)
Fossorial_conv.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#EBCC2A"),lty=3),
                               points.args=list(pch=c(23,21),bg=c("#D9D0D3","#EBCC2A")),
                               legend.args=list(pch=c(23,21),x="top"))
Fossorial_conv.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                              polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                              line.args=list(line.col="#EBCC2A",lwd=4))
# Snout
Fossorial_conv_Snout.plot <- plotConv(SC=Fossorial_conv_snout.test,y=paca_snout,variable=1,state=Fossorial_vector)
Fossorial_conv_Snout.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#EBCC2A"),lty=3),
                                     points.args=list(pch=c(23,21),bg=c("#D9D0D3","#EBCC2A")),
                                     legend.args=list(pch=c(23,21),x="top"))
Fossorial_conv_Snout.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                    polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                    line.args=list(line.col="#EBCC2A",lwd=4))
# Braincase
Fossorial_conv_braincase.plot <- plotConv(SC=Fossorial_conv_braincase.test,y=paca_braincase,variable=1,state=Fossorial_vector)
Fossorial_conv_braincase.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#EBCC2A"),lty=3),
                                         points.args=list(pch=c(23,21),bg=c("#D9D0D3","#EBCC2A")),
                                         legend.args=list(pch=c(23,21),x="top"))
Fossorial_conv_braincase.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                        polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                        line.args=list(line.col="#EBCC2A",lwd=4))
# Maxilla
Fossorial_conv_maxilla.plot <- plotConv(SC=Fossorial_conv_maxilla.test,y=paca_maxilla,variable=1,state=Fossorial_vector)
Fossorial_conv_maxilla.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#EBCC2A"),lty=3),
                                       points.args=list(pch=c(23,21),bg=c("#D9D0D3","#EBCC2A")),
                                       legend.args=list(pch=c(23,21),x="top"))
Fossorial_conv_maxilla.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                      polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                      line.args=list(line.col="#EBCC2A",lwd=4))
# Palatopterygoid Arch
Fossorial_conv_palatoptery.plot <- plotConv(SC=Fossorial_conv_palatoptery.test,y=paca_palatoptery,variable=1,state=Fossorial_vector)
Fossorial_conv_palatoptery.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#EBCC2A"),lty=3),
                                           points.args=list(pch=c(23,21),bg=c("#D9D0D3","#EBCC2A")),
                                           legend.args=list(pch=c(23,21),x="top"))
Fossorial_conv_palatoptery.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                          polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                          line.args=list(line.col="#EBCC2A",lwd=4))
# Mandible
Fossorial_conv_mandible.plot <- plotConv(SC=Fossorial_conv_mandible.test,y=paca_mandible,variable=1,state=Fossorial_vector)
Fossorial_conv_mandible.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#EBCC2A"),lty=3),
                                        points.args=list(pch=c(23,21),bg=c("#D9D0D3","#EBCC2A")),
                                        legend.args=list(pch=c(23,21),x="top"))
Fossorial_conv_mandible.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                       polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                       line.args=list(line.col="#EBCC2A",lwd=4))
# Suspensorium
Fossorial_conv_suspensorium.plot <- plotConv(SC=Fossorial_conv_suspensorium.test,y=paca_suspensorium,variable=1,state=Fossorial_vector)
Fossorial_conv_suspensorium.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#EBCC2A"),lty=3),
                                            points.args=list(pch=c(23,21),bg=c("#D9D0D3","#EBCC2A")),
                                            legend.args=list(pch=c(23,21),x="top"))
Fossorial_conv_suspensorium.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                           polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                           line.args=list(line.col="#EBCC2A",lwd=4))

# Diet - Generalist
# Skull
Generalist_conv.plot <- plotConv(SC=Generalist_conv_skull.test,y=paca_skull,variable=1,state=Generalist_vector)
Generalist_conv.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#DD3C51"),lty=3),
                                points.args=list(pch=c(23,21),bg=c("#D9D0D3","#DD3C51")),
                                legend.args=list(pch=c(23,21),x="top"))
Generalist_conv.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                               polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                               line.args=list(line.col="#DD3C51",lwd=4))
# Snout
Generalist_conv_Snout.plot <- plotConv(SC=Generalist_conv_snout.test,y=paca_snout,variable=1,state=Generalist_vector)
Generalist_conv_Snout.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#DD3C51"),lty=3),
                                      points.args=list(pch=c(23,21),bg=c("#D9D0D3","#DD3C51")),
                                      legend.args=list(pch=c(23,21),x="top"))
Generalist_conv_Snout.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                     polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                     line.args=list(line.col="#DD3C51",lwd=4))
# Braincase
Generalist_conv_braincase.plot <- plotConv(SC=Generalist_conv_braincase.test,y=paca_braincase,variable=1,state=Generalist_vector)
Generalist_conv_braincase.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#DD3C51"),lty=3),
                                          points.args=list(pch=c(23,21),bg=c("#D9D0D3","#DD3C51")),
                                          legend.args=list(pch=c(23,21),x="top"))
Generalist_conv_braincase.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                         polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                         line.args=list(line.col="#DD3C51",lwd=4))
# Maxilla
Generalist_conv_maxilla.plot <- plotConv(SC=Generalist_conv_maxilla.test,y=paca_maxilla,variable=1,state=Generalist_vector)
Generalist_conv_maxilla.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#DD3C51"),lty=3),
                                        points.args=list(pch=c(23,21),bg=c("#D9D0D3","#DD3C51")),
                                        legend.args=list(pch=c(23,21),x="top"))
Generalist_conv_maxilla.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                       polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                       line.args=list(line.col="#DD3C51",lwd=4))
# Palatopterygoid Arch
Generalist_conv_palatoptery.plot <- plotConv(SC=Generalist_conv_palatoptery.test,y=paca_palatoptery,variable=1,state=Generalist_vector)
Generalist_conv_palatoptery.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#DD3C51"),lty=3),
                                            points.args=list(pch=c(23,21),bg=c("#D9D0D3","#DD3C51")),
                                            legend.args=list(pch=c(23,21),x="top"))
Generalist_conv_palatoptery.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                           polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                           line.args=list(line.col="#DD3C51",lwd=4))
# Mandible
Generalist_conv_mandible.plot <- plotConv(SC=Generalist_conv_mandible.test,y=paca_mandible,variable=1,state=Generalist_vector)
Generalist_conv_mandible.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#DD3C51"),lty=3),
                                         points.args=list(pch=c(23,21),bg=c("#D9D0D3","#DD3C51")),
                                         legend.args=list(pch=c(23,21),x="top"))
Generalist_conv_mandible.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                        polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                        line.args=list(line.col="#DD3C51",lwd=4))
# Suspensorium
Generalist_conv_suspensorium.plot <- plotConv(SC=Generalist_conv_suspensorium.test,y=paca_suspensorium,variable=1,state=Generalist_vector)
Generalist_conv_suspensorium.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#DD3C51"),lty=3),
                                             points.args=list(pch=c(23,21),bg=c("#D9D0D3","#DD3C51")),
                                             legend.args=list(pch=c(23,21),x="top"))
Generalist_conv_suspensorium.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                            polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                            line.args=list(line.col="#DD3C51",lwd=4))
# Diet - Aquatic Generalist
# Skull
Aquatic_Generalist_conv.plot <- plotConv(SC=Aquatic_Generalist_conv_skull.test,y=paca_skull,variable=1,state=Aquatic_Generalist_vector)
Aquatic_Generalist_conv.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#46ACC8"),lty=3),
                                        points.args=list(pch=c(23,21),bg=c("#D9D0D3","#46ACC8")),
                                        legend.args=list(pch=c(23,21),x="top"))
Aquatic_Generalist_conv.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                       polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                       line.args=list(line.col="#46ACC8",lwd=4))
# Snout
Aquatic_Generalist_conv_Snout.plot <- plotConv(SC=Aquatic_Generalist_conv_snout.test,y=paca_snout,variable=1,state=Aquatic_Generalist_vector)
Aquatic_Generalist_conv_Snout.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#46ACC8"),lty=3),
                                              points.args=list(pch=c(23,21),bg=c("#D9D0D3","#46ACC8")),
                                              legend.args=list(pch=c(23,21),x="top"))
Aquatic_Generalist_conv_Snout.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                             polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                             line.args=list(line.col="#46ACC8",lwd=4))
# Braincase
Aquatic_Generalist_conv_braincase.plot <- plotConv(SC=Aquatic_Generalist_conv_braincase.test,y=paca_braincase,variable=1,state=Aquatic_Generalist_vector)
Aquatic_Generalist_conv_braincase.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#46ACC8"),lty=3),
                                                  points.args=list(pch=c(23,21),bg=c("#D9D0D3","#46ACC8")),
                                                  legend.args=list(pch=c(23,21),x="top"))
Aquatic_Generalist_conv_braincase.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                                 polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                                 line.args=list(line.col="#46ACC8",lwd=4))
# Maxilla
Aquatic_Generalist_conv_maxilla.plot <- plotConv(SC=Aquatic_Generalist_conv_maxilla.test,y=paca_maxilla,variable=1,state=Aquatic_Generalist_vector)
Aquatic_Generalist_conv_maxilla.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#46ACC8"),lty=3),
                                                points.args=list(pch=c(23,21),bg=c("#D9D0D3","#46ACC8")),
                                                legend.args=list(pch=c(23,21),x="top"))
Aquatic_Generalist_conv_maxilla.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                               polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                               line.args=list(line.col="#46ACC8",lwd=4))
# Palatopterygoid Arch
Aquatic_Generalist_conv_palatoptery.plot <- plotConv(SC=Aquatic_Generalist_conv_palatoptery.test,y=paca_palatoptery,variable=1,state=Aquatic_Generalist_vector)
Aquatic_Generalist_conv_palatoptery.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#46ACC8"),lty=3),
                                                    points.args=list(pch=c(23,21),bg=c("#D9D0D3","#46ACC8")),
                                                    legend.args=list(pch=c(23,21),x="top"))
Aquatic_Generalist_conv_palatoptery.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                                   polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                                   line.args=list(line.col="#46ACC8",lwd=4))
# Mandible
Aquatic_Generalist_conv_mandible.plot <- plotConv(SC=Aquatic_Generalist_conv_mandible.test,y=paca_mandible,variable=1,state=Aquatic_Generalist_vector)
Aquatic_Generalist_conv_mandible.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#46ACC8"),lty=3),
                                                 points.args=list(pch=c(23,21),bg=c("#D9D0D3","#46ACC8")),
                                                 legend.args=list(pch=c(23,21),x="top"))
Aquatic_Generalist_conv_mandible.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                                polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                                line.args=list(line.col="#46ACC8",lwd=4))
# Suspensorium
Aquatic_Generalist_conv_suspensorium.plot <- plotConv(SC=Aquatic_Generalist_conv_suspensorium.test,y=paca_suspensorium,variable=1,state=Aquatic_Generalist_vector)
Aquatic_Generalist_conv_suspensorium.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#46ACC8"),lty=3),
                                                     points.args=list(pch=c(23,21),bg=c("#D9D0D3","#46ACC8")),
                                                     legend.args=list(pch=c(23,21),x="top"))
Aquatic_Generalist_conv_suspensorium.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                                    polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                                    line.args=list(line.col="#46ACC8",lwd=4))
# Diet - Vermivorous
# Skull
Vermivorous_conv.plot <- plotConv(SC=Vermivorous_conv_skull.test,y=paca_skull,variable=1,state=Vermivorous_vector)
Vermivorous_conv.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F98400"),lty=3),
                                 points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F98400")),
                                 legend.args=list(pch=c(23,21),x="top"))
Vermivorous_conv.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                line.args=list(line.col="#F98400",lwd=4))
# Snout
Vermivorous_conv_Snout.plot <- plotConv(SC=Vermivorous_conv_snout.test,y=paca_snout,variable=1,state=Vermivorous_vector)
Vermivorous_conv_Snout.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F98400"),lty=3),
                                       points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F98400")),
                                       legend.args=list(pch=c(23,21),x="top"))
Vermivorous_conv_Snout.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                      polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                      line.args=list(line.col="#F98400",lwd=4))
# Braincase
Vermivorous_conv_braincase.plot <- plotConv(SC=Vermivorous_conv_braincase.test,y=paca_braincase,variable=1,state=Vermivorous_vector)
Vermivorous_conv_braincase.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F98400"),lty=3),
                                           points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F98400")),
                                           legend.args=list(pch=c(23,21),x="top"))
Vermivorous_conv_braincase.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                          polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                          line.args=list(line.col="#F98400",lwd=4))
# Maxilla
Vermivorous_conv_maxilla.plot <- plotConv(SC=Vermivorous_conv_maxilla.test,y=paca_maxilla,variable=1,state=Vermivorous_vector)
Vermivorous_conv_maxilla.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F98400"),lty=3),
                                         points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F98400")),
                                         legend.args=list(pch=c(23,21),x="top"))
Vermivorous_conv_maxilla.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                        polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                        line.args=list(line.col="#F98400",lwd=4))
# Palatopterygoid Arch
Vermivorous_conv_palatoptery.plot <- plotConv(SC=Vermivorous_conv_palatoptery.test,y=paca_palatoptery,variable=1,state=Vermivorous_vector)
Vermivorous_conv_palatoptery.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F98400"),lty=3),
                                             points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F98400")),
                                             legend.args=list(pch=c(23,21),x="top"))
Vermivorous_conv_palatoptery.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                            polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                            line.args=list(line.col="#F98400",lwd=4))
# Mandible
Vermivorous_conv_mandible.plot <- plotConv(SC=Vermivorous_conv_mandible.test,y=paca_mandible,variable=1,state=Vermivorous_vector)
Vermivorous_conv_mandible.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F98400"),lty=3),
                                          points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F98400")),
                                          legend.args=list(pch=c(23,21),x="top"))
Vermivorous_conv_mandible.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                         polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                         line.args=list(line.col="#F98400",lwd=4))
# Suspensorium
Vermivorous_conv_suspensorium.plot <- plotConv(SC=Vermivorous_conv_suspensorium.test,y=paca_suspensorium,variable=1,state=Vermivorous_vector)
Vermivorous_conv_suspensorium.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F98400"),lty=3),
                                              points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F98400")),
                                              legend.args=list(pch=c(23,21),x="top"))
Vermivorous_conv_suspensorium.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                             polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                             line.args=list(line.col="#F98400",lwd=4))

# Diet - Durophagous
# Skull
Durophagous_conv.plot <- plotConv(SC=Durophagous_conv_skull.test,y=paca_skull,variable=1,state=Durophagous_vector)
Durophagous_conv.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#33A02C"),lty=3),
                                 points.args=list(pch=c(23,21),bg=c("#D9D0D3","#33A02C")),
                                 legend.args=list(pch=c(23,21),x="top"))
Durophagous_conv.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                line.args=list(line.col="#33A02C",lwd=4))
# Snout
Durophagous_conv_Snout.plot <- plotConv(SC=Durophagous_conv_snout.test,y=paca_snout,variable=1,state=Durophagous_vector)
Durophagous_conv_Snout.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#33A02C"),lty=3),
                                       points.args=list(pch=c(23,21),bg=c("#D9D0D3","#33A02C")),
                                       legend.args=list(pch=c(23,21),x="top"))
Durophagous_conv_Snout.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                      polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                      line.args=list(line.col="#33A02C",lwd=4))
# Braincase
Durophagous_conv_braincase.plot <- plotConv(SC=Durophagous_conv_braincase.test,y=paca_braincase,variable=1,state=Durophagous_vector)
Durophagous_conv_braincase.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#33A02C"),lty=3),
                                           points.args=list(pch=c(23,21),bg=c("#D9D0D3","#33A02C")),
                                           legend.args=list(pch=c(23,21),x="top"))
Durophagous_conv_braincase.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                          polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                          line.args=list(line.col="#33A02C",lwd=4))
# Maxilla
Durophagous_conv_maxilla.plot <- plotConv(SC=Durophagous_conv_maxilla.test,y=paca_maxilla,variable=1,state=Durophagous_vector)
Durophagous_conv_maxilla.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#33A02C"),lty=3),
                                         points.args=list(pch=c(23,21),bg=c("#D9D0D3","#33A02C")),
                                         legend.args=list(pch=c(23,21),x="top"))
Durophagous_conv_maxilla.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                        polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                        line.args=list(line.col="#33A02C",lwd=4))
# Palatopterygoid Arch
Durophagous_conv_palatoptery.plot <- plotConv(SC=Durophagous_conv_palatoptery.test,y=paca_palatoptery,variable=1,state=Durophagous_vector)
Durophagous_conv_palatoptery.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#33A02C"),lty=3),
                                             points.args=list(pch=c(23,21),bg=c("#D9D0D3","#33A02C")),
                                             legend.args=list(pch=c(23,21),x="top"))
Durophagous_conv_palatoptery.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                            polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                            line.args=list(line.col="#33A02C",lwd=4))
# Mandible
Durophagous_conv_mandible.plot <- plotConv(SC=Durophagous_conv_mandible.test,y=paca_mandible,variable=1,state=Durophagous_vector)
Durophagous_conv_mandible.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#33A02C"),lty=3),
                                          points.args=list(pch=c(23,21),bg=c("#D9D0D3","#33A02C")),
                                          legend.args=list(pch=c(23,21),x="top"))
Durophagous_conv_mandible.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                         polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                         line.args=list(line.col="#33A02C",lwd=4))
# Suspensorium
Durophagous_conv_suspensorium.plot <- plotConv(SC=Durophagous_conv_suspensorium.test,y=paca_suspensorium,variable=1,state=Durophagous_vector)
Durophagous_conv_suspensorium.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#33A02C"),lty=3),
                                              points.args=list(pch=c(23,21),bg=c("#D9D0D3","#33A02C")),
                                              legend.args=list(pch=c(23,21),x="top"))
Durophagous_conv_suspensorium.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                             polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                             line.args=list(line.col="#33A02C",lwd=4))
