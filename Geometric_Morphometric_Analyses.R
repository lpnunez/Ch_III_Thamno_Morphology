library(tidyverse)
library(Morpho)
library(geomorph)
library(Rvcg)
library(paleomorph)
library(EMMLi)
library(qgraph)
library(ape)
library(geiger)
library(abind)
library(phytools)
library(doParallel)
library(RRphylo)
library(patchwork)
library(gridExtra)
library(treeio)
library(extrafont)
library(deeptime)
library(ggtree)
library(hot.dots)
library(SURGE) 
library(rjson)
library(data.table)
library(treedata.table)
library(phytools)
library(geiger)
library(convevol)
library(dispRity)
library(wesanderson)

setwd("C:/Users/maver/Dropbox (AMNH)/Thamnophiini/Morph_Data/µCT/")
source("C:/Users/maver/Dropbox (AMNH)/Thamnophiini/Scripts/R_Scripts/SlicerMorph_Rexamples-main/SlicerMorph_Rexamples-main/read.markups.json.R")
source("C:/Users/maver/Dropbox (AMNH)/Thamnophiini/Scripts/From_Folks/Danny_Rhoda/Matched_Local_Superimpositions.R")
meta_lms <- read.csv("Landmarks_Bilat.csv")
meta_lms_left <- read.csv("Landmarks.csv")
meta_lms_full <- read.csv("Landmarks_Bilat_Curves.csv")
lms_names <- meta_lms$Name
Thamno_Tree <- drop.tip(read.tree("Trees/Thamno_WGS1_82_ALRT_CCR_PLC.tre"), c("Natrix_natrix","Thamnophis_fulvus","Thamnophis_sirtalis_X_radix"))
Thamno_Tree <- force.ultrametric(Thamno_Tree)

Female_coords <- readland.tps("C:/Users/maver/Dropbox (AMNH)/Thamnophiini/Morph_Data/µCT/Landmarks/All_Taxa_PTS/Female_Thamnophiini_Ana_Curves_UCE_CatNum.tps",specID = c("ID"))
#Female_coords <- readland.tps("C:/Users/maver/Dropbox (AMNH)/Thamnophiini/Morph_Data/µCT/Landmarks/All_Taxa_PTS/Female_Thamnophiini_Ana_Curves_UCE_CatNum_Final.tps",specID = c("ID"))
Male_coords <- readland.tps("C:/Users/maver/Dropbox (AMNH)/Thamnophiini/Morph_Data/µCT/Landmarks/All_Taxa_PTS/Male_Thamnophiini_Ana_Curves_UCE_CatNum.tps",specID = c("ID"))
All_coords <- abind(Female_coords, Male_coords)

# Mirror Data
bilat.landmarks <- cbind(which(meta_lms$Position=="l"),which(meta_lms$Position=="r"))

bilat.landmarks <- bilat.landmarks[c(1:70),]

midline<- which(meta_lms$Position=="m")#a list of which landmarks fall on the midline

all.lm.present <- c(1:922)
ana.lm.present <- c(1:152)

left.curves <- all.lm.present[-midline[13:length(midline)]][-ana.lm.present]#a list of the curves which are going to be mirrored to the other side (excludes midline curves)
leftside<-c(bilat.landmarks[,1],left.curves)
num.missing<-(length(leftside)-length(bilat.landmarks[,2]))
blanks<-c((dim(All_coords)[1]+1):(dim(All_coords)[1]+num.missing))         
rightside<-c(bilat.landmarks[,2],blanks)

module_defs_left <- meta_lms %>% filter(., Position != "r")

add_col_or_row = function(x, n = 1, add_col = T, fill = 0)
{
  m1 = matrix(x, ncol = if(add_col) nrow(x) * ncol(x) else nrow(x), byrow = T)
  m2 = matrix(fill, nrow = if(add_col) dim(x)[3] else prod(dim(x)[-1]), 
              ncol = if(add_col) nrow(x) * n else n)
  array(t(cbind(m1, m2)), 
        c(nrow(x) + ((!add_col) * n), ncol(x) + (add_col * n), dim(x)[3]))
}

coords2<-add_col_or_row(All_coords,n=num.missing,add_col=FALSE,fill=NA)
dimnames(coords2)[3]<-dimnames(All_coords)[3]
bilats<-cbind(leftside, rightside)
newarray<-paleomorph::mirrorfill(coords2,l1=as.integer(midline),l2=bilats)
dimnames(newarray)[3]<-dimnames(All_coords)[3]


# Slide Semilandmarks
ana.lm.present <- c(1:152)
sm_sliders <- list(c(153:162),c(163:172),c(173:182),c(183:192),c(193:202),c(203:212),c(213:222),c(223:232),c(233:242),c(243:252),c(253:262),c(263:272),c(273:282),c(283:292),c(293:302),c(303:312),c(313:322),c(323:332),c(333:342),c(343:352),c(353:362),c(363:372),c(373:382),c(383:392),c(393:402),c(403:412),c(413:422),c(423:432),c(433:442),c(443:452),c(453:462),c(463:472),c(473:482),c(483:492),c(493:502),c(503:512),c(513:522),c(523:532),c(533:542),c(543:552),c(553:562),c(563:572),c(573:582),c(583:592),c(593:602),c(603:612),c(613:622),c(623:632),c(633:642),c(643:652),c(653:662),c(663:672),c(673:682),c(683:692),c(693:702),c(703:712),c(713:722),c(723:732),c(733:742),c(743:752),c(753:762),c(763:772),c(773:782),c(783:792),c(793:802),c(803:812),c(813:822),c(823:832),c(833:842),c(843:852),c(853:862),c(863:872),c(873:882),c(883:892),c(893:902),c(903:912),c(913:922),c(913:922),c(923:932),c(933:942),c(943:952),c(953:962),c(963:972),c(973:982),c(983:992),c(993:1002),c(1003:1012),c(1013:1022),c(1023:1032),c(1033:1042),c(1043:1052),c(1053:1062),c(1063:1072),c(1073:1082),c(1083:1092),c(1093:1102),c(1103:1112),c(1113:1122),c(1123:1132),c(1133:1142),c(1143:1152),c(1153:1162),c(1163:1172),c(1173:1182),c(1183:1192),c(1193:1202),c(1203:1212),c(1213:1222),c(1223:1232),c(1233:1242),c(1243:1252),c(1253:1262),c(1263:1272),c(1273:1282),c(1283:1292),c(1293:1302),c(1303:1312),c(1313:1322),c(1323:1332),c(1333:1342),c(1343:1352),c(1353:1362),c(1363:1372),c(1373:1382),c(1383:1392),c(1393:1402),c(1403:1412),c(1413:1422),c(1423:1432),c(1433:1442),c(1443:1452),c(1453:1462),c(1463:1472),c(1473:1482),c(1483:1492),c(1493:1502),c(1503:1512),c(1513:1522),c(1523:1532),c(1533:1542),c(1543:1552),c(1553:1562),c(1563:1572),c(1573:1582),c(1583:1592),c(1593:1602),c(1603:1612),c(1613:1622),c(1623:1632))
slid_coords <- slider3d(newarray, SMvector=ana.lm.present, deselect=T, outlines=sm_sliders, iterations=0,mc.cores=1)
save(slid_coords, file = "Coords_Data/slid_coords.RData")
load("Coords_Data/slid_coords.RData")

# Use Local Superimposition for Each Partition
lgpa.partition <- meta_lms_full$full_model #Full Skull
test <- local.gpa(coords = slid_coords$dataslide, partition = lgpa.partition) #Full Skull

# Read in Info on specimens and Taxa
scanned <- read.csv("Specimens_Used.csv", header=T)

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

# Full Skull
Skull_Test <- test[1:1592,,]
Skull_Physig <- physignal(Skull_Test, All_Tree)
pc_test <- gm.prcomp(Skull_Test)
phypca_test <- gm.prcomp(Skull_Test, phy = All_Tree, GLS = T)
write.csv(phypca_test$x, "PC_Scores/Skull_PhyPCA_PC.csv")
paca_test <- gm.prcomp(Skull_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_test$x, "PC_Scores/Skull_PACA_PC.csv")
make_ggplot(plot(phypca_test, phylo = T))
make_ggplot(plot(paca_test, phylo = T))
plot(Skull_Physig)

#Premaxilla
Premax_Test <- test[1:76,,]
Premax_Physig <- physignal(Premax_Test, All_Tree)
phypca_premax <- gm.prcomp(Premax_Test, phy = All_Tree, GLS = T)
write.csv(phypca_premax$x, "PC_Scores/Premaxilla_PhyPCA_PC.csv")
paca_premax <- gm.prcomp(Premax_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_premax$x, "PC_Scores/Premaxilla_PACA_PC.csv")
make_ggplot(plot(phypca_premax, phylo = T))
make_ggplot(plot(paca_premax, phylo = T))
plot(Premax_Physig)

#Nasal
Nasal_Test <- test[77:152,,]
Nasal_Physig <- physignal(Nasal_Test, All_Tree)
phypca_Nasal <- gm.prcomp(Nasal_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Nasal$x, "PC_Scores/Nasal_PhyPCA_PC.csv")
paca_Nasal <- gm.prcomp(Nasal_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Nasal$x, "PC_Scores/Nasal_PACA_PC.csv")
make_ggplot(plot(phypca_Nasal, phylo = T))
make_ggplot(plot(paca_Nasal, phylo = T))
plot(Nasal_Physig)

#Frontal
Frontal_Test <- test[153:228,,]
Frontal_Physig <- physignal(Frontal_Test, All_Tree)
phypca_Frontal <- gm.prcomp(Frontal_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Frontal$x, "PC_Scores/Frontal_PhyPCA_PC.csv")
paca_Frontal <- gm.prcomp(Frontal_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Frontal$x, "PC_Scores/Frontal_PACA_PC.csv")
make_ggplot(plot(phypca_Frontal, phylo = T))
make_ggplot(plot(paca_Frontal, phylo = T))
plot(Frontal_Physig)

#Parietal
Parietal_Test <- test[229:304,,]
Parietal_Physig <- physignal(Parietal_Test, All_Tree)
phypca_Parietal <- gm.prcomp(Parietal_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Parietal$x, "PC_Scores/Parietal_PhyPCA_PC.csv")
paca_Parietal <- gm.prcomp(Parietal_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Parietal$x, "PC_Scores/Parietal_PACA_PC.csv")
make_ggplot(plot(phypca_Parietal, phylo = T))
make_ggplot(plot(paca_Parietal, phylo = T))
plot(Parietal_Physig)

#Supraoccipital
Supraocc_Test <- test[305:402,,]
Supraocc_Physig <- physignal(Supraocc_Test, All_Tree)
phypca_Supraocc <- gm.prcomp(Supraocc_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Supraocc$x, "PC_Scores/Supraoccipital_PhyPCA_PC.csv")
paca_Supraocc <- gm.prcomp(Supraocc_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Supraocc$x, "PC_Scores/Supraoccipital_PACA_PC.csv")
make_ggplot(plot(phypca_Supraocc, phylo = T))
make_ggplot(plot(paca_Supraocc, phylo = T))
plot(Supraocc_Physig)

#Basisphenoid
Basi_Test <- test[403:478,,]
Basi_Physig <- physignal(Basi_Test, All_Tree)
phypca_Basi <- gm.prcomp(Basi_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Basi$x, "PC_Scores/Basisphenoid_PhyPCA_PC.csv")
paca_Basi <- gm.prcomp(Basi_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Basi$x, "PC_Scores/Basisphenoid_PACA_PC.csv")
make_ggplot(plot(phypca_Basi, phylo = T))
make_ggplot(plot(paca_Basi, phylo = T))
plot(Basi_Physig)

#Maxilla
Maxilla_Test <- test[c(479:484,491:560),,] #Left Only
Maxilla_Physig <- physignal(Maxilla_Test, All_Tree)
phypca_Maxilla <- gm.prcomp(Maxilla_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Maxilla$x, "PC_Scores/Maxilla_PhyPCA_PC.csv")
paca_Maxilla <- gm.prcomp(Maxilla_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Maxilla$x, "PC_Scores/Maxilla_PACA_PC.csv")
make_ggplot(plot(phypca_Maxilla, phylo = T))
make_ggplot(plot(paca_Maxilla, phylo = T))
plot(Maxilla_Physig)

#Palatine
Palatine_Test <- test[c(631:640,651:740),,]
Palatine_Physig <- physignal(Palatine_Test, All_Tree)
phypca_Palatine <- gm.prcomp(Palatine_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Palatine$x, "PC_Scores/Palatine_PhyPCA_PC.csv")
paca_Palatine <- gm.prcomp(Palatine_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Palatine$x, "PC_Scores/Palatine_PACA_PC.csv")
make_ggplot(plot(phypca_Palatine, phylo = T))
make_ggplot(plot(paca_Palatine, phylo = T))
plot(Palatine_Physig)

#Ectopterygoid
Ecto_Test <- test[c(831:833,837:866),,]
Ecto_Physig <- physignal(Ecto_Test, All_Tree)
phypca_Ecto <- gm.prcomp(Ecto_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Ecto$x, "PC_Scores/Ectopterygoid_PhyPCA_PC.csv")
paca_Ecto <- gm.prcomp(Ecto_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Ecto$x, "PC_Scores/Ectopterygoid_PACA_PC.csv")
make_ggplot(plot(phypca_Ecto, phylo = T))
make_ggplot(plot(paca_Ecto, phylo = T))
plot(Ecto_Physig)

#Pterygoid
Pterygoid_Test <- test[c(897:901,907:976),,]
Pterygoid_Physig <- physignal(Pterygoid_Test, All_Tree)
phypca_Pterygoid <- gm.prcomp(Pterygoid_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Pterygoid$x, "PC_Scores/Pterygoid_PhyPCA_PC.csv")
paca_Pterygoid <- gm.prcomp(Pterygoid_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Pterygoid$x, "PC_Scores/Pterygoid_PACA_PC.csv")
make_ggplot(plot(phypca_Pterygoid, phylo = T))
make_ggplot(plot(paca_Pterygoid, phylo = T))
plot(Pterygoid_Physig)

#Dentary
Dentary_Test <- test[c(1047:1053,1061:1120),,]
Dentary_Physig <- physignal(Dentary_Test, All_Tree)
phypca_Dentary <- gm.prcomp(Dentary_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Dentary$x, "PC_Scores/Dentary_PhyPCA_PC.csv")
paca_Dentary <- gm.prcomp(Dentary_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Dentary$x, "PC_Scores/Dentary_PACA_PC.csv")
make_ggplot(plot(phypca_Dentary, phylo = T))
make_ggplot(plot(paca_Dentary, phylo = T))
plot(Dentary_Physig)

#Compound
Compound_Test <- test[c(1181:1192,1205:1264),,]
Compound_Physig <- physignal(Compound_Test, All_Tree)
phypca_Compound <- gm.prcomp(Compound_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Compound$x, "PC_Scores/Compound_PhyPCA_PC.csv")
paca_Compound <- gm.prcomp(Compound_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Compound$x, "PC_Scores/Compound_PACA_PC.csv")
make_ggplot(plot(phypca_Compound, phylo = T))
make_ggplot(plot(paca_Compound, phylo = T))
plot(Compound_Physig)

#Quadrate
Quadrate_Test <- test[c(1325:1336,1349:1448),,]
Quadrate_Physig <- physignal(Quadrate_Test, All_Tree)
phypca_Quadrate <- gm.prcomp(Quadrate_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Quadrate$x, "PC_Scores/Quadrate_PhyPCA_PC.csv")
paca_Quadrate <- gm.prcomp(Quadrate_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Quadrate$x, "PC_Scores/Quadrate_PACA_PC.csv")
make_ggplot(plot(phypca_Quadrate, phylo = T))
make_ggplot(plot(paca_Quadrate, phylo = T))
plot(Quadrate_Physig)

#Supratemporal
Supratemp_Test <- test[c(1549:1550,1553:1572),,]
Supratemp_Physig <- physignal(Supratemp_Test, All_Tree)
phypca_Supratemp <- gm.prcomp(Supratemp_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Supratemp$x, "PC_Scores/Supratemporal_PhyPCA_PC.csv")
paca_Supratemp <- gm.prcomp(Supratemp_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Supratemp$x, "PC_Scores/Supratemporal_PACA_PC.csv")
make_ggplot(plot(phypca_Supratemp, phylo = T))
make_ggplot(plot(paca_Supratemp, phylo = T))
plot(Supratemp_Physig)

#Snout
Snout_Test <- test[1:152,,]
Snout_Physig <- physignal(Snout_Test, All_Tree)
phypca_Snout <- gm.prcomp(Snout_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Snout$x, "PC_Scores/Snout_PhyPCA_PC.csv")
paca_Snout <- gm.prcomp(Snout_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Snout$x, "PC_Scores/Snout_PACA_PC.csv")
make_ggplot(plot(phypca_Snout, phylo = T))
make_ggplot(plot(paca_Snout, phylo = T))
plot(Snout_Physig)

#Braincase
Braincase_Test <- test[153:478,,]
Braincase_Physig <- physignal(Braincase_Test, All_Tree)
phypca_Braincase <- gm.prcomp(Braincase_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Braincase$x, "PC_Scores/Braincase_PhyPCA_PC.csv")
paca_Braincase <- gm.prcomp(Braincase_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Braincase$x, "PC_Scores/Braincase_PACA_PC.csv")
make_ggplot(plot(phypca_Braincase, phylo = T))
make_ggplot(plot(paca_Braincase, phylo = T))
plot(Braincase_Physig)

#Palatomaxilary_Arch
Palatomax_Test <- test[c(479:484,491:560,631:640,651:740,831:833,837:866,897:901,907:976),,]
Palatomax_Physig <- physignal(Palatomax_Test, All_Tree)
phypca_Palatomax <- gm.prcomp(Palatomax_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Palatomax$x, "PC_Scores/Palatomaxillary_Arch_PhyPCA_PC.csv")
paca_Palatomax <- gm.prcomp(Palatomax_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Palatomax$x, "PC_Scores/Palatomaxillary_Arch_PACA_PC.csv")
make_ggplot(plot(phypca_Palatomax, phylo = T))
make_ggplot(plot(paca_Palatomax, phylo = T))
plot(Palatomax_Physig)

#Palatopterygoid Arch
Palatoptery_Test <- test[c(631:640,651:740,831:833,837:866,897:901,907:976),,]
Palatoptery_Physig <- physignal(Palatoptery_Test, All_Tree)
phypca_Palatoptery <- gm.prcomp(Palatoptery_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Palatoptery$x, "PC_Scores/Palatopterygoid_Arch_PhyPCA_PC.csv")
paca_Palatoptery <- gm.prcomp(Palatoptery_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Palatoptery$x, "PC_Scores/Palatopterygoid_Arch_PACA_PC.csv")
make_ggplot(plot(phypca_Palatoptery, phylo = T))
make_ggplot(plot(paca_Palatoptery, phylo = T))
plot(Palatoptery_Physig)

#Mandible
Mandible_Test <- test[c(1047:1053,1061:1120,1181:1192,1205:1264),,]
Mandible_Physig <- physignal(Mandible_Test, All_Tree)
phypca_Mandible <- gm.prcomp(Mandible_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Mandible$x, "PC_Scores/Mandible_PhyPCA_PC.csv")
paca_Mandible <- gm.prcomp(Mandible_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Mandible$x, "PC_Scores/Mandible_PACA_PC.csv")
make_ggplot(plot(phypca_Mandible, phylo = T))
make_ggplot(plot(paca_Mandible, phylo = T))
plot(Mandible_Physig)

#Suspensorium
Suspensorium_Test <- test[c(1325:1336,1349:1448,1549:1550,1553:1572),,]
Suspensorium_Physig <- physignal(Suspensorium_Test, All_Tree)
phypca_Suspensorium <- gm.prcomp(Suspensorium_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Suspensorium$x, "PC_Scores/Suspensorium_PhyPCA_PC.csv")
paca_Suspensorium <- gm.prcomp(Suspensorium_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Suspensorium$x, "PC_Scores/Suspensorium_PACA_PC.csv")
make_ggplot(plot(phypca_Suspensorium, phylo = T))
make_ggplot(plot(paca_Suspensorium, phylo = T))
plot(Suspensorium_Physig)

#Lower Jaw
Lower_Jaw_Test <- test[c(1047:1053,1061:1120,1181:1192,1205:1264,1325:1336,1349:1448,1549:1550,1553:1572),,]
Lower_Jaw_Physig <- physignal(Lower_Jaw_Test, All_Tree)
phypca_Lower_Jaw <- gm.prcomp(Lower_Jaw_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Lower_Jaw$x, "PC_Scores/Lower_Jaw_PhyPCA_PC.csv")
paca_Lower_Jaw <- gm.prcomp(Lower_Jaw_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Lower_Jaw$x, "PC_Scores/Lower_Jaw_PACA_PC.csv")
make_ggplot(plot(phypca_Lower_Jaw, phylo = T))
make_ggplot(plot(paca_Lower_Jaw, phylo = T))
plot(Lower_Jaw_Physig)

#Feeding Apparatus
Feed_App_Test <- test[c(479:484,491:560,631:640,651:740,831:833,837:866,897:901,907:976,1047:1053,1061:1120,1181:1192,1205:1264,1325:1336,1349:1448,1549:1550,1553:1572),,]
Feed_App_Physig <- physignal(Feed_App_Test, All_Tree)
phypca_Feed_App <- gm.prcomp(Feed_App_Test, phy = All_Tree, GLS = T)
write.csv(phypca_Feed_App$x, "PC_Scores/Feeding_Apparatus_PhyPCA_PC.csv")
paca_Feed_App <- gm.prcomp(Feed_App_Test, phy = All_Tree, align.to.phy = T)
write.csv(paca_Feed_App$x, "PC_Scores/Feeding_Apparatus_PACA_PC.csv")
make_ggplot(plot(phypca_Feed_App, phylo = T))
make_ggplot(plot(paca_Feed_App, phylo = T))
plot(Feed_App_Physig)

# Read in Trait Data for Convergence Analyses
traits <- read.csv("C:/Users/maver/Dropbox (AMNH)/Thamnophiini/All_Thamno_Traits.csv")
Thamno <- as.treedata.table(Thamno_Tree, traits, name_column = "Name_in_Tree")

# Extract Trait Data from treedata.table
diet <- extractVector(Thamno, "Diet")
habitat <- extractVector(Thamno, "Lifestyle")
clade <- extractVector(Thamno, "Clade")
genus <- extractVector(Thamno, "Genus")

# Convergence Tests using RRphylo
Partitions <- c("Skull", "Snout", "Braincase", "Palatomaxillary_Arch", "Palatopterygoid_Arch", "Mandible", "Suspensorium", "Lower_Jaw", "Feeding_Apparatus", "Premaxilla", "Nasal", "Frontal", "Parietal", "Supraoccipital", "Basisphenoid", "Maxilla", "Palatine", "Ectopterygoid", "Pterygoid", "Dentary", "Compound", "Quadrate", "Supratemporal")
Reduced_Partitions <- c("Skull", "Snout", "Braincase", "Maxilla", "Palatopterygoid_Arch", "Mandible", "Suspensorium"
                        )
# Habitat - Aquatic
aquatic_vector <- rep("nostate", 51)
names(aquatic_vector)<-rownames(paca_test$x)
Aquatic_Taxa <- which(habitat=="Aquatic")
Aquatic_Taxa <- names(Aquatic_Taxa)
aquatic_vector[which(names(aquatic_vector)%in%Aquatic_Taxa)] <- "aquatic"
aquatic_conv.test <- search.conv(tree=All_Tree, y=paca_test$x[,1:4],  state = aquatic_vector, declust = T)
aquatic_conv_snout.test <- search.conv(tree=All_Tree, y=paca_Snout$x[,1:4],  state = aquatic_vector, declust = T)
aquatic_conv_Braincase.test <- search.conv(tree=All_Tree, y=paca_Braincase$x[,1:4],  state = aquatic_vector, declust = T)
aquatic_conv_Maxilla.test <- search.conv(tree=All_Tree, y=paca_Maxilla$x[,1:4],  state = aquatic_vector, declust = T)
aquatic_conv_Palatoptery.test <- search.conv(tree=All_Tree, y=paca_Palatoptery$x[,1:4],  state = aquatic_vector, declust = T)
aquatic_conv_Mandible.test <- search.conv(tree=All_Tree, y=paca_Mandible$x[,1:4],  state = aquatic_vector, declust = T)
aquatic_conv_Suspensorium.test <- search.conv(tree=All_Tree, y=paca_Suspensorium$x[,1:4],  state = aquatic_vector, declust = T)
#aquatic_conv_Palatomax.test <- search.conv(tree=All_Tree, y=paca_Palatomax$x[,1:4],  state = aquatic_vector, declust = T)
#aquatic_conv_Lower_Jaw.test <- search.conv(tree=All_Tree, y=paca_Lower_Jaw$x[,1:4],  state = aquatic_vector, declust = T)
#aquatic_conv_Feed_App.test <- search.conv(tree=All_Tree, y=paca_Feed_App$x[,1:4],  state = aquatic_vector, declust = T)
#aquatic_conv_Premax.test <- search.conv(tree=All_Tree, y=paca_premax$x[,1:4],  state = aquatic_vector, declust = T)
#aquatic_conv_Nasal.test <- search.conv(tree=All_Tree, y=paca_Nasal$x[,1:4],  state = aquatic_vector, declust = T)
#aquatic_conv_Frontal.test <- search.conv(tree=All_Tree, y=paca_Frontal$x[,1:4],  state = aquatic_vector, declust = T)
#aquatic_conv_Parietal.test <- search.conv(tree=All_Tree, y=paca_Parietal$x[,1:4],  state = aquatic_vector, declust = T)
#aquatic_conv_Supraocc.test <- search.conv(tree=All_Tree, y=paca_Supraocc$x[,1:4],  state = aquatic_vector, declust = T)
#aquatic_conv_Basi.test <- search.conv(tree=All_Tree, y=paca_Basi$x[,1:4],  state = aquatic_vector, declust = T)
#aquatic_conv_Palatine.test <- search.conv(tree=All_Tree, y=paca_Palatine$x[,1:4],  state = aquatic_vector, declust = T)
#aquatic_conv_Ecto.test <- search.conv(tree=All_Tree, y=paca_Ecto$x[,1:4],  state = aquatic_vector, declust = T)
#aquatic_conv_Pterygoid.test <- search.conv(tree=All_Tree, y=paca_Pterygoid$x[,1:4],  state = aquatic_vector, declust = T)
#aquatic_conv_Dentary.test <- search.conv(tree=All_Tree, y=paca_Dentary$x[,1:4],  state = aquatic_vector, declust = T)
#aquatic_conv_Compound.test <- search.conv(tree=All_Tree, y=paca_Compound$x[,1:4],  state = aquatic_vector, declust = T)
#aquatic_conv_Quadrate.test <- search.conv(tree=All_Tree, y=paca_Quadrate$x[,1:4],  state = aquatic_vector, declust = T)
#aquatic_conv_Supratemp.test <- search.conv(tree=All_Tree, y=paca_Supratemp$x[,1:4],  state = aquatic_vector, declust = T)
#Aquatic_Habitat_Conv_Tests <- rbind(aquatic_conv.test$state.res,aquatic_conv_snout.test$state.res,aquatic_conv_Braincase.test$state.res,aquatic_conv_Palatomax.test$state.res,aquatic_conv_Palatoptery.test$state.res,aquatic_conv_Mandible.test$state.res,aquatic_conv_Suspensorium.test$state.res,aquatic_conv_Lower_Jaw.test$state.res,aquatic_conv_Feed_App.test$state.res,aquatic_conv_Premax.test$state.res,aquatic_conv_Nasal.test$state.res,aquatic_conv_Frontal.test$state.res,aquatic_conv_Parietal.test$state.res,aquatic_conv_Supraocc.test$state.res,aquatic_conv_Basi.test$state.res,aquatic_conv_Maxilla.test$state.res,aquatic_conv_Palatine.test$state.res,aquatic_conv_Ecto.test$state.res,aquatic_conv_Pterygoid.test$state.res,aquatic_conv_Dentary.test$state.res,aquatic_conv_Compound.test$state.res,aquatic_conv_Quadrate.test$state.res,aquatic_conv_Supratemp.test$state.res)
#rownames(Aquatic_Habitat_Conv_Tests) <- Partitions
#write.csv(Aquatic_Habitat_Conv_Tests, file = "Convergence/Tables/Aquatic_Habitat_Conv_Tests.csv")

# Habitat - Terrestrial
terrestrial_vector <- rep("nostate", 51)
names(terrestrial_vector)<-rownames(paca_test$x)
Terrestrial_Taxa <- which(habitat=="Terrestrial")
Terrestrial_Taxa <- names(Terrestrial_Taxa)
terrestrial_vector[which(names(terrestrial_vector)%in%Terrestrial_Taxa)] <- "terrestrial"
terrestrial_conv.test <- search.conv(tree=All_Tree, y=paca_test$x[,1:4],  state = terrestrial_vector, declust = T)
Terrestrial_conv_snout.test <- search.conv(tree=All_Tree, y=paca_Snout$x[,1:4],  state = terrestrial_vector, declust = T)
terrestrial_conv_Braincase.test <- search.conv(tree=All_Tree, y=paca_Braincase$x[,1:4],  state = terrestrial_vector, declust = T)
terrestrial_conv_Maxilla.test <- search.conv(tree=All_Tree, y=paca_Maxilla$x[,1:4],  state = terrestrial_vector, declust = T)
terrestrial_conv_Palatoptery.test <- search.conv(tree=All_Tree, y=paca_Palatoptery$x[,1:4],  state = terrestrial_vector, declust = T)
terrestrial_conv_Mandible.test <- search.conv(tree=All_Tree, y=paca_Mandible$x[,1:4],  state = terrestrial_vector, declust = T)
terrestrial_conv_Suspensorium.test <- search.conv(tree=All_Tree, y=paca_Suspensorium$x[,1:4],  state = terrestrial_vector, declust = T)
#terrestrial_conv_Palatomax.test <- search.conv(tree=All_Tree, y=paca_Palatomax$x[,1:4],  state = terrestrial_vector, declust = T)
#terrestrial_conv_Lower_Jaw.test <- search.conv(tree=All_Tree, y=paca_Lower_Jaw$x[,1:4],  state = terrestrial_vector, declust = T)
#terrestrial_conv_Feed_App.test <- search.conv(tree=All_Tree, y=paca_Feed_App$x[,1:4],  state = terrestrial_vector, declust = T)
#terrestrial_conv_Premax.test <- search.conv(tree=All_Tree, y=paca_premax$x[,1:4],  state = terrestrial_vector, declust = T)
#terrestrial_conv_Nasal.test <- search.conv(tree=All_Tree, y=paca_Nasal$x[,1:4],  state = terrestrial_vector, declust = T)
#terrestrial_conv_Frontal.test <- search.conv(tree=All_Tree, y=paca_Frontal$x[,1:4],  state = terrestrial_vector, declust = T)
#terrestrial_conv_Parietal.test <- search.conv(tree=All_Tree, y=paca_Parietal$x[,1:4],  state = terrestrial_vector, declust = T)
#terrestrial_conv_Supraocc.test <- search.conv(tree=All_Tree, y=paca_Supraocc$x[,1:4],  state = terrestrial_vector, declust = T)
#terrestrial_conv_Basi.test <- search.conv(tree=All_Tree, y=paca_Basi$x[,1:4],  state = terrestrial_vector, declust = T)
#terrestrial_conv_Palatine.test <- search.conv(tree=All_Tree, y=paca_Palatine$x[,1:4],  state = terrestrial_vector, declust = T)
#terrestrial_conv_Ecto.test <- search.conv(tree=All_Tree, y=paca_Ecto$x[,1:4],  state = terrestrial_vector, declust = T)
#terrestrial_conv_Pterygoid.test <- search.conv(tree=All_Tree, y=paca_Pterygoid$x[,1:4],  state = terrestrial_vector, declust = T)
#terrestrial_conv_Dentary.test <- search.conv(tree=All_Tree, y=paca_Dentary$x[,1:4],  state = terrestrial_vector, declust = T)
#terrestrial_conv_Compound.test <- search.conv(tree=All_Tree, y=paca_Compound$x[,1:4],  state = terrestrial_vector, declust = T)
#terrestrial_conv_Quadrate.test <- search.conv(tree=All_Tree, y=paca_Quadrate$x[,1:4],  state = terrestrial_vector, declust = T)
#terrestrial_conv_Supratemp.test <- search.conv(tree=All_Tree, y=paca_Supratemp$x[,1:4],  state = terrestrial_vector, declust = T)
#Terrestrial_Habitat_Conv_Tests <- rbind(terrestrial_conv.test$state.res,Terrestrial_conv_snout.test$state.res,terrestrial_conv_Braincase.test$state.res,terrestrial_conv_Palatomax.test$state.res,terrestrial_conv_Palatoptery.test$state.res,terrestrial_conv_Mandible.test$state.res,terrestrial_conv_Suspensorium.test$state.res,terrestrial_conv_Lower_Jaw.test$state.res,terrestrial_conv_Feed_App.test$state.res,terrestrial_conv_Premax.test$state.res,terrestrial_conv_Nasal.test$state.res,terrestrial_conv_Frontal.test$state.res,terrestrial_conv_Parietal.test$state.res,terrestrial_conv_Supraocc.test$state.res,terrestrial_conv_Basi.test$state.res,terrestrial_conv_Maxilla.test$state.res,terrestrial_conv_Palatine.test$state.res,terrestrial_conv_Ecto.test$state.res,terrestrial_conv_Pterygoid.test$state.res,terrestrial_conv_Dentary.test$state.res,terrestrial_conv_Compound.test$state.res,terrestrial_conv_Quadrate.test$state.res,terrestrial_conv_Supratemp.test$state.res)
#rownames(Terrestrial_Habitat_Conv_Tests) <- Partitions
#write.csv(Terrestrial_Habitat_Conv_Tests, file = "Convergence/Tables/Terrestrial_Habitat_Conv_Tests.csv")

# Habitat - Fossorial
Fossorial_vector <- rep("nostate", 51)
names(Fossorial_vector)<-rownames(paca_test$x)
Fossorial_Taxa <- which(habitat=="Fossorial")
Fossorial_Taxa <- names(Fossorial_Taxa)
Fossorial_vector[which(names(Fossorial_vector)%in%Fossorial_Taxa)] <- "Fossorial"
Fossorial_conv.test <- search.conv(tree=All_Tree, y=paca_test$x[,1:4],  state = Fossorial_vector, declust = T)
Fossorial_conv_snout.test <- search.conv(tree=All_Tree, y=paca_Snout$x[,1:4],  state = Fossorial_vector, declust = T)
Fossorial_conv_Braincase.test <- search.conv(tree=All_Tree, y=paca_Braincase$x[,1:4],  state = Fossorial_vector, declust = T)
Fossorial_conv_Maxilla.test <- search.conv(tree=All_Tree, y=paca_Maxilla$x[,1:4],  state = Fossorial_vector, declust = T)
Fossorial_conv_Palatoptery.test <- search.conv(tree=All_Tree, y=paca_Palatoptery$x[,1:4],  state = Fossorial_vector, declust = T)
Fossorial_conv_Mandible.test <- search.conv(tree=All_Tree, y=paca_Mandible$x[,1:4],  state = Fossorial_vector, declust = T)
Fossorial_conv_Suspensorium.test <- search.conv(tree=All_Tree, y=paca_Suspensorium$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_conv_Palatomax.test <- search.conv(tree=All_Tree, y=paca_Palatomax$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_conv_Lower_Jaw.test <- search.conv(tree=All_Tree, y=paca_Lower_Jaw$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_conv_Feed_App.test <- search.conv(tree=All_Tree, y=paca_Feed_App$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_conv_Premax.test <- search.conv(tree=All_Tree, y=paca_premax$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_conv_Nasal.test <- search.conv(tree=All_Tree, y=paca_Nasal$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_conv_Frontal.test <- search.conv(tree=All_Tree, y=paca_Frontal$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_conv_Parietal.test <- search.conv(tree=All_Tree, y=paca_Parietal$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_conv_Supraocc.test <- search.conv(tree=All_Tree, y=paca_Supraocc$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_conv_Basi.test <- search.conv(tree=All_Tree, y=paca_Basi$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_conv_Palatine.test <- search.conv(tree=All_Tree, y=paca_Palatine$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_conv_Ecto.test <- search.conv(tree=All_Tree, y=paca_Ecto$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_conv_Pterygoid.test <- search.conv(tree=All_Tree, y=paca_Pterygoid$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_conv_Dentary.test <- search.conv(tree=All_Tree, y=paca_Dentary$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_conv_Compound.test <- search.conv(tree=All_Tree, y=paca_Compound$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_conv_Quadrate.test <- search.conv(tree=All_Tree, y=paca_Quadrate$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_conv_Supratemp.test <- search.conv(tree=All_Tree, y=paca_Supratemp$x[,1:4],  state = Fossorial_vector, declust = T)
#Fossorial_Habitat_Conv_Tests <- rbind(Fossorial_conv.test$state.res,Fossorial_conv_snout.test$state.res,Fossorial_conv_Braincase.test$state.res,Fossorial_conv_Palatomax.test$state.res,Fossorial_conv_Palatoptery.test$state.res,Fossorial_conv_Mandible.test$state.res,Fossorial_conv_Suspensorium.test$state.res,Fossorial_conv_Lower_Jaw.test$state.res,Fossorial_conv_Feed_App.test$state.res,Fossorial_conv_Premax.test$state.res,Fossorial_conv_Nasal.test$state.res,Fossorial_conv_Frontal.test$state.res,Fossorial_conv_Parietal.test$state.res,Fossorial_conv_Supraocc.test$state.res,Fossorial_conv_Basi.test$state.res,Fossorial_conv_Maxilla.test$state.res,Fossorial_conv_Palatine.test$state.res,Fossorial_conv_Ecto.test$state.res,Fossorial_conv_Pterygoid.test$state.res,Fossorial_conv_Dentary.test$state.res,Fossorial_conv_Compound.test$state.res,Fossorial_conv_Quadrate.test$state.res,Fossorial_conv_Supratemp.test$state.res)
#rownames(Fossorial_Habitat_Conv_Tests) <- Partitions
#write.csv(Fossorial_Habitat_Conv_Tests, file = "Convergence/Tables/Fossorial_Habitat_Conv_Tests.csv")

# Diet - Generalist
Generalist_vector <- rep("nostate", 51)
names(Generalist_vector)<-rownames(paca_test$x)
Generalist_Taxa <- which(diet=="Generalist")
Generalist_Taxa <- names(Generalist_Taxa)
Generalist_vector[which(names(Generalist_vector)%in%Generalist_Taxa)] <- "Generalist"
Generalist_conv.test <- search.conv(tree=All_Tree, y=paca_test$x[,1:4],  state = Generalist_vector, declust = T)
Generalist_conv_snout.test <- search.conv(tree=All_Tree, y=paca_Snout$x[,1:4],  state = Generalist_vector, declust = T)
Generalist_conv_Braincase.test <- search.conv(tree=All_Tree, y=paca_Braincase$x[,1:4],  state = Generalist_vector, declust = T)
Generalist_conv_Maxilla.test <- search.conv(tree=All_Tree, y=paca_Maxilla$x[,1:4],  state = Generalist_vector, declust = T)
Generalist_conv_Palatoptery.test <- search.conv(tree=All_Tree, y=paca_Palatoptery$x[,1:4],  state = Generalist_vector, declust = T)
Generalist_conv_Mandible.test <- search.conv(tree=All_Tree, y=paca_Mandible$x[,1:4],  state = Generalist_vector, declust = T)
Generalist_conv_Suspensorium.test <- search.conv(tree=All_Tree, y=paca_Suspensorium$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_conv_Palatomax.test <- search.conv(tree=All_Tree, y=paca_Palatomax$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_conv_Lower_Jaw.test <- search.conv(tree=All_Tree, y=paca_Lower_Jaw$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_conv_Feed_App.test <- search.conv(tree=All_Tree, y=paca_Feed_App$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_conv_Premax.test <- search.conv(tree=All_Tree, y=paca_premax$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_conv_Nasal.test <- search.conv(tree=All_Tree, y=paca_Nasal$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_conv_Frontal.test <- search.conv(tree=All_Tree, y=paca_Frontal$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_conv_Parietal.test <- search.conv(tree=All_Tree, y=paca_Parietal$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_conv_Supraocc.test <- search.conv(tree=All_Tree, y=paca_Supraocc$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_conv_Basi.test <- search.conv(tree=All_Tree, y=paca_Basi$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_conv_Palatine.test <- search.conv(tree=All_Tree, y=paca_Palatine$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_conv_Ecto.test <- search.conv(tree=All_Tree, y=paca_Ecto$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_conv_Pterygoid.test <- search.conv(tree=All_Tree, y=paca_Pterygoid$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_conv_Dentary.test <- search.conv(tree=All_Tree, y=paca_Dentary$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_conv_Compound.test <- search.conv(tree=All_Tree, y=paca_Compound$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_conv_Quadrate.test <- search.conv(tree=All_Tree, y=paca_Quadrate$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_conv_Supratemp.test <- search.conv(tree=All_Tree, y=paca_Supratemp$x[,1:4],  state = Generalist_vector, declust = T)
#Generalist_Diet_Conv_Tests <- rbind(Generalist_conv.test$state,Generalist_conv_snout.test$state,Generalist_conv_Braincase.test$state,Generalist_conv_Palatomax.test$state,Generalist_conv_Palatoptery.test$state,Generalist_conv_Mandible.test$state,Generalist_conv_Suspensorium.test$state,Generalist_conv_Lower_Jaw.test$state,Generalist_conv_Feed_App.test$state,Generalist_conv_Premax.test$state,Generalist_conv_Nasal.test$state,Generalist_conv_Frontal.test$state,Generalist_conv_Parietal.test$state,Generalist_conv_Supraocc.test$state,Generalist_conv_Basi.test$state,Generalist_conv_Maxilla.test$state,Generalist_conv_Palatine.test$state,Generalist_conv_Ecto.test$state,Generalist_conv_Pterygoid.test$state,Generalist_conv_Dentary.test$state,Generalist_conv_Compound.test$state,Generalist_conv_Quadrate.test$state,Generalist_conv_Supratemp.test$state)
#rownames(Generalist_Diet_Conv_Tests) <- Partitions
#write.csv(Generalist_Diet_Conv_Tests, file = "Convergence/Tables/Generalist_Diet_Conv_Tests.csv")

# Diet - Aquatic Generalist
Aquatic_Generalist_vector <- rep("nostate", 51)
names(Aquatic_Generalist_vector)<-rownames(paca_test$x)
Aquatic_Generalist_Taxa <- which(diet=="Aquatic_Generalist")
Aquatic_Generalist_Taxa <- names(Aquatic_Generalist_Taxa)
Aquatic_Generalist_vector[which(names(Aquatic_Generalist_vector)%in%Aquatic_Generalist_Taxa)] <- "Aquatic_Generalist"
Aquatic_Generalist_conv.test <- search.conv(tree=All_Tree, y=paca_test$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
Aquatic_Generalist_conv_snout.test <- search.conv(tree=All_Tree, y=paca_Snout$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
Aquatic_Generalist_conv_Braincase.test <- search.conv(tree=All_Tree, y=paca_Braincase$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
Aquatic_Generalist_conv_Maxilla.test <- search.conv(tree=All_Tree, y=paca_Maxilla$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
Aquatic_Generalist_conv_Palatoptery.test <- search.conv(tree=All_Tree, y=paca_Palatoptery$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
Aquatic_Generalist_conv_Mandible.test <- search.conv(tree=All_Tree, y=paca_Mandible$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
Aquatic_Generalist_conv_Suspensorium.test <- search.conv(tree=All_Tree, y=paca_Suspensorium$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_conv_Palatomax.test <- search.conv(tree=All_Tree, y=paca_Palatomax$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_conv_Lower_Jaw.test <- search.conv(tree=All_Tree, y=paca_Lower_Jaw$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_conv_Feed_App.test <- search.conv(tree=All_Tree, y=paca_Feed_App$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_conv_Premax.test <- search.conv(tree=All_Tree, y=paca_premax$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_conv_Nasal.test <- search.conv(tree=All_Tree, y=paca_Nasal$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_conv_Frontal.test <- search.conv(tree=All_Tree, y=paca_Frontal$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_conv_Parietal.test <- search.conv(tree=All_Tree, y=paca_Parietal$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_conv_Supraocc.test <- search.conv(tree=All_Tree, y=paca_Supraocc$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_conv_Basi.test <- search.conv(tree=All_Tree, y=paca_Basi$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_conv_Palatine.test <- search.conv(tree=All_Tree, y=paca_Palatine$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_conv_Ecto.test <- search.conv(tree=All_Tree, y=paca_Ecto$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_conv_Pterygoid.test <- search.conv(tree=All_Tree, y=paca_Pterygoid$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_conv_Dentary.test <- search.conv(tree=All_Tree, y=paca_Dentary$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_conv_Compound.test <- search.conv(tree=All_Tree, y=paca_Compound$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_conv_Quadrate.test <- search.conv(tree=All_Tree, y=paca_Quadrate$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_conv_Supratemp.test <- search.conv(tree=All_Tree, y=paca_Supratemp$x[,1:4],  state = Aquatic_Generalist_vector, declust = T)
#Aquatic_Generalist_Diet_Conv_Tests <- rbind(Aquatic_Generalist_conv.test$state,Aquatic_Generalist_conv_snout.test$state,Aquatic_Generalist_conv_Braincase.test$state,Aquatic_Generalist_conv_Palatomax.test$state,Aquatic_Generalist_conv_Palatoptery.test$state,Aquatic_Generalist_conv_Mandible.test$state,Aquatic_Generalist_conv_Suspensorium.test$state,Aquatic_Generalist_conv_Lower_Jaw.test$state,Aquatic_Generalist_conv_Feed_App.test$state,Aquatic_Generalist_conv_Premax.test$state,Aquatic_Generalist_conv_Nasal.test$state,Aquatic_Generalist_conv_Frontal.test$state,Aquatic_Generalist_conv_Parietal.test$state,Aquatic_Generalist_conv_Supraocc.test$state,Aquatic_Generalist_conv_Basi.test$state,Aquatic_Generalist_conv_Maxilla.test$state,Aquatic_Generalist_conv_Palatine.test$state,Aquatic_Generalist_conv_Ecto.test$state,Aquatic_Generalist_conv_Pterygoid.test$state,Aquatic_Generalist_conv_Dentary.test$state,Aquatic_Generalist_conv_Compound.test$state,Aquatic_Generalist_conv_Quadrate.test$state,Aquatic_Generalist_conv_Supratemp.test$state)
#rownames(Aquatic_Generalist_Diet_Conv_Tests) <- Partitions
#write.csv(Aquatic_Generalist_Diet_Conv_Tests, file = "Convergence/Tables/Aquatic_Generalist_Diet_Conv_Tests.csv")

# Diet - Vermivorous
Vermivorous_vector <- rep("nostate", 51)
names(Vermivorous_vector)<-rownames(paca_test$x)
Vermivorous_Taxa <- which(diet=="Vermivorous")
Vermivorous_Taxa <- names(Vermivorous_Taxa)
Vermivorous_vector[which(names(Vermivorous_vector)%in%Vermivorous_Taxa)] <- "Vermivorous"
Vermivorous_conv.test <- search.conv(tree=All_Tree, y=paca_test$x[,1:4],  state = Vermivorous_vector, declust = T)
Vermivorous_conv_snout.test <- search.conv(tree=All_Tree, y=paca_Snout$x[,1:4],  state = Vermivorous_vector, declust = T)
Vermivorous_conv_Braincase.test <- search.conv(tree=All_Tree, y=paca_Braincase$x[,1:4],  state = Vermivorous_vector, declust = T)
Vermivorous_conv_Maxilla.test <- search.conv(tree=All_Tree, y=paca_Maxilla$x[,1:4],  state = Vermivorous_vector, declust = T)
Vermivorous_conv_Palatoptery.test <- search.conv(tree=All_Tree, y=paca_Palatoptery$x[,1:4],  state = Vermivorous_vector, declust = T)
Vermivorous_conv_Mandible.test <- search.conv(tree=All_Tree, y=paca_Mandible$x[,1:4],  state = Vermivorous_vector, declust = T)
Vermivorous_conv_Suspensorium.test <- search.conv(tree=All_Tree, y=paca_Suspensorium$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_conv_Palatomax.test <- search.conv(tree=All_Tree, y=paca_Palatomax$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_conv_Lower_Jaw.test <- search.conv(tree=All_Tree, y=paca_Palatomax$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_conv_Feed_App.test <- search.conv(tree=All_Tree, y=paca_Feed_App$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_conv_Premax.test <- search.conv(tree=All_Tree, y=paca_premax$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_conv_Nasal.test <- search.conv(tree=All_Tree, y=paca_Nasal$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_conv_Frontal.test <- search.conv(tree=All_Tree, y=paca_Frontal$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_conv_Parietal.test <- search.conv(tree=All_Tree, y=paca_Parietal$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_conv_Supraocc.test <- search.conv(tree=All_Tree, y=paca_Supraocc$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_conv_Basi.test <- search.conv(tree=All_Tree, y=paca_Basi$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_conv_Palatine.test <- search.conv(tree=All_Tree, y=paca_Palatine$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_conv_Ecto.test <- search.conv(tree=All_Tree, y=paca_Ecto$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_conv_Pterygoid.test <- search.conv(tree=All_Tree, y=paca_Pterygoid$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_conv_Dentary.test <- search.conv(tree=All_Tree, y=paca_Dentary$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_conv_Compound.test <- search.conv(tree=All_Tree, y=paca_Compound$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_conv_Quadrate.test <- search.conv(tree=All_Tree, y=paca_Quadrate$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_conv_Supratemp.test <- search.conv(tree=All_Tree, y=paca_Supratemp$x[,1:4],  state = Vermivorous_vector, declust = T)
#Vermivorous_Diet_Conv_Tests <- rbind(Vermivorous_conv.test$state,Vermivorous_conv_snout.test$state,Vermivorous_conv_Braincase.test$state,Vermivorous_conv_Palatomax.test$state,Vermivorous_conv_Palatoptery.test$state,Vermivorous_conv_Mandible.test$state,Vermivorous_conv_Suspensorium.test$state,Vermivorous_conv_Lower_Jaw.test$state,Vermivorous_conv_Feed_App.test$state,Vermivorous_conv_Premax.test$state,Vermivorous_conv_Nasal.test$state,Vermivorous_conv_Frontal.test$state,Vermivorous_conv_Parietal.test$state,Vermivorous_conv_Supraocc.test$state,Vermivorous_conv_Basi.test$state,Vermivorous_conv_Maxilla.test$state,Vermivorous_conv_Palatine.test$state,Vermivorous_conv_Ecto.test$state,Vermivorous_conv_Pterygoid.test$state,Vermivorous_conv_Dentary.test$state,Vermivorous_conv_Compound.test$state,Vermivorous_conv_Quadrate.test$state,Vermivorous_conv_Supratemp.test$state)
#rownames(Vermivorous_Diet_Conv_Tests) <- Partitions
#write.csv(Vermivorous_Diet_Conv_Tests, file = "Convergence/Tables/Vermivorous_Diet_Conv_Tests.csv")

# Diet - Durophagous
Durophagous_vector <- rep("nostate", 51)
names(Durophagous_vector)<-rownames(paca_test$x)
Durophagous_Taxa <- which(diet=="Durophagous")
Durophagous_Taxa <- names(Durophagous_Taxa)
Durophagous_vector[which(names(Durophagous_vector)%in%Durophagous_Taxa)] <- "Durophagous"
Durophagous_conv.test <- search.conv(tree=All_Tree, y=paca_test$x[,1:4],  state = Durophagous_vector, declust = T)
Durophagous_conv_snout.test <- search.conv(tree=All_Tree, y=paca_Snout$x[,1:4],  state = Durophagous_vector, declust = T)
Durophagous_conv_Braincase.test <- search.conv(tree=All_Tree, y=paca_Braincase$x[,1:4],  state = Durophagous_vector, declust = T)
Durophagous_conv_Maxilla.test <- search.conv(tree=All_Tree, y=paca_Maxilla$x[,1:4],  state = Durophagous_vector, declust = T)
Durophagous_conv_Palatoptery.test <- search.conv(tree=All_Tree, y=paca_Palatoptery$x[,1:4],  state = Durophagous_vector, declust = T)
Durophagous_conv_Mandible.test <- search.conv(tree=All_Tree, y=paca_Mandible$x[,1:4],  state = Durophagous_vector, declust = T)
Durophagous_conv_Suspensorium.test <- search.conv(tree=All_Tree, y=paca_Suspensorium$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_conv_Palatomax.test <- search.conv(tree=All_Tree, y=paca_Palatomax$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_conv_Lower_Jaw.test <- search.conv(tree=All_Tree, y=paca_Lower_Jaw$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_conv_Feed_App.test <- search.conv(tree=All_Tree, y=paca_Feed_App$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_conv_Premax.test <- search.conv(tree=All_Tree, y=paca_premax$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_conv_Nasal.test <- search.conv(tree=All_Tree, y=paca_Nasal$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_conv_Frontal.test <- search.conv(tree=All_Tree, y=paca_Frontal$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_conv_Parietal.test <- search.conv(tree=All_Tree, y=paca_Parietal$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_conv_Supraocc.test <- search.conv(tree=All_Tree, y=paca_Supraocc$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_conv_Basi.test <- search.conv(tree=All_Tree, y=paca_Basi$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_conv_Palatine.test <- search.conv(tree=All_Tree, y=paca_Palatine$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_conv_Ecto.test <- search.conv(tree=All_Tree, y=paca_Ecto$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_conv_Pterygoid.test <- search.conv(tree=All_Tree, y=paca_Pterygoid$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_conv_Dentary.test <- search.conv(tree=All_Tree, y=paca_Dentary$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_conv_Compound.test <- search.conv(tree=All_Tree, y=paca_Compound$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_conv_Quadrate.test <- search.conv(tree=All_Tree, y=paca_Quadrate$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_conv_Supratemp.test <- search.conv(tree=All_Tree, y=paca_Supratemp$x[,1:4],  state = Durophagous_vector, declust = T)
#Durophagous_Diet_Conv_Tests <- rbind(Durophagous_conv.test$state,Durophagous_conv_snout.test$state,Durophagous_conv_Braincase.test$state,Durophagous_conv_Palatomax.test$state,Durophagous_conv_Palatoptery.test$state,Durophagous_conv_Mandible.test$state,Durophagous_conv_Suspensorium.test$state,Durophagous_conv_Lower_Jaw.test$state,Durophagous_conv_Feed_App.test$state,Durophagous_conv_Premax.test$state,Durophagous_conv_Nasal.test$state,Durophagous_conv_Frontal.test$state,Durophagous_conv_Parietal.test$state,Durophagous_conv_Supraocc.test$state,Durophagous_conv_Basi.test$state,Durophagous_conv_Maxilla.test$state,Durophagous_conv_Palatine.test$state,Durophagous_conv_Ecto.test$state,Durophagous_conv_Pterygoid.test$state,Durophagous_conv_Dentary.test$state,Durophagous_conv_Compound.test$state,Durophagous_conv_Quadrate.test$state,Durophagous_conv_Supratemp.test$state)
#rownames(Durophagous_Diet_Conv_Tests) <- Partitions
#write.csv(Durophagous_Diet_Conv_Tests, file = "Convergence/Tables/Durophagous_Diet_Conv_Tests.csv")


# Plot Convergence Test Results
# Habitat - Aquatic
# Skull
aquatic_conv.plot <- plotConv(SC=aquatic_conv.test,y=paca_test$x[,1:4],variable=1,state=aquatic_vector)
aquatic_conv.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#3B9AB2"),lty=3),
                             points.args=list(pch=c(23,21),bg=c("#D9D0D3","#3B9AB2")),
                             legend.args=list(pch=c(23,21),x="top"))
aquatic_conv.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                            polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                            line.args=list(line.col="#3B9AB2",lwd=4))
# Snout
aquatic_conv_Snout.plot <- plotConv(SC=aquatic_conv_snout.test,y=paca_Snout$x[,1:4],variable=1,state=aquatic_vector)
aquatic_conv_Snout.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#3B9AB2"),lty=3),
                                         points.args=list(pch=c(23,21),bg=c("#D9D0D3","#3B9AB2")),
                                         legend.args=list(pch=c(23,21),x="top"))
aquatic_conv_Snout.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                        polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                        line.args=list(line.col="#3B9AB2",lwd=4))
# Braincase
aquatic_conv_Braincase.plot <- plotConv(SC=aquatic_conv_Braincase.test,y=paca_Braincase$x[,1:4],variable=1,state=aquatic_vector)
aquatic_conv_Braincase.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#3B9AB2"),lty=3),
                             points.args=list(pch=c(23,21),bg=c("#D9D0D3","#3B9AB2")),
                             legend.args=list(pch=c(23,21),x="top"))
aquatic_conv_Braincase.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                            polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                            line.args=list(line.col="#3B9AB2",lwd=4))
# Maxilla
aquatic_conv_Maxilla.plot <- plotConv(SC=aquatic_conv_Maxilla.test,y=paca_Maxilla$x[,1:4],variable=1,state=aquatic_vector)
aquatic_conv_Maxilla.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#3B9AB2"),lty=3),
                                         points.args=list(pch=c(23,21),bg=c("#D9D0D3","#3B9AB2")),
                                         legend.args=list(pch=c(23,21),x="top"))
aquatic_conv_Maxilla.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                        polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                        line.args=list(line.col="#3B9AB2",lwd=4))
# Palatopterygoid Arch
aquatic_conv_Palatoptery.plot <- plotConv(SC=aquatic_conv_Palatoptery.test,y=paca_Palatoptery$x[,1:4],variable=1,state=aquatic_vector)
aquatic_conv_Palatoptery.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#3B9AB2"),lty=3),
                                         points.args=list(pch=c(23,21),bg=c("#D9D0D3","#3B9AB2")),
                                         legend.args=list(pch=c(23,21),x="top"))
aquatic_conv_Palatoptery.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                        polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                        line.args=list(line.col="#3B9AB2",lwd=4))
# Mandible
aquatic_conv_Mandible.plot <- plotConv(SC=aquatic_conv_Mandible.test,y=paca_Mandible$x[,1:4],variable=1,state=aquatic_vector)
aquatic_conv_Mandible.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#3B9AB2"),lty=3),
                                         points.args=list(pch=c(23,21),bg=c("#D9D0D3","#3B9AB2")),
                                         legend.args=list(pch=c(23,21),x="top"))
aquatic_conv_Mandible.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                        polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                        line.args=list(line.col="#3B9AB2",lwd=4))
# Suspensorium
aquatic_conv_Suspensorium.plot <- plotConv(SC=aquatic_conv_Suspensorium.test,y=paca_Suspensorium$x[,1:4],variable=1,state=aquatic_vector)
aquatic_conv_Suspensorium.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#3B9AB2"),lty=3),
                                         points.args=list(pch=c(23,21),bg=c("#D9D0D3","#3B9AB2")),
                                         legend.args=list(pch=c(23,21),x="top"))
aquatic_conv_Suspensorium.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                        polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                        line.args=list(line.col="#3B9AB2",lwd=4))

# Habitat - Terrestrial
# Skull
terrestrial_conv.plot <- plotConv(SC=terrestrial_conv.test,y=paca_test$x[,1:4],variable=1,state=terrestrial_vector)
terrestrial_conv.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F21A00"),lty=3),
                             points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F21A00")),
                             legend.args=list(pch=c(23,21),x="top"))
terrestrial_conv.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                            polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                            line.args=list(line.col="#F21A00",lwd=4))
# Snout
terrestrial_conv_Snout.plot <- plotConv(SC=Terrestrial_conv_snout.test,y=paca_Snout$x[,1:4],variable=1,state=terrestrial_vector)
terrestrial_conv_Snout.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F21A00"),lty=3),
                                   points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F21A00")),
                                   legend.args=list(pch=c(23,21),x="top"))
terrestrial_conv_Snout.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                  polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                  line.args=list(line.col="#F21A00",lwd=4))
# Braincase
terrestrial_conv_Braincase.plot <- plotConv(SC=terrestrial_conv_Braincase.test,y=paca_Braincase$x[,1:4],variable=1,state=terrestrial_vector)
terrestrial_conv_Braincase.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F21A00"),lty=3),
                                       points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F21A00")),
                                       legend.args=list(pch=c(23,21),x="top"))
terrestrial_conv_Braincase.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                      polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                      line.args=list(line.col="#F21A00",lwd=4))
# Maxilla
terrestrial_conv_Maxilla.plot <- plotConv(SC=terrestrial_conv_Maxilla.test,y=paca_Maxilla$x[,1:4],variable=1,state=terrestrial_vector)
terrestrial_conv_Maxilla.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F21A00"),lty=3),
                                     points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F21A00")),
                                     legend.args=list(pch=c(23,21),x="top"))
terrestrial_conv_Maxilla.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                    polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                    line.args=list(line.col="#F21A00",lwd=4))
# Palatopterygoid Arch
terrestrial_conv_Palatoptery.plot <- plotConv(SC=terrestrial_conv_Palatoptery.test,y=paca_Palatoptery$x[,1:4],variable=1,state=terrestrial_vector)
terrestrial_conv_Palatoptery.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F21A00"),lty=3),
                                         points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F21A00")),
                                         legend.args=list(pch=c(23,21),x="top"))
terrestrial_conv_Palatoptery.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                        polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                        line.args=list(line.col="#F21A00",lwd=4))
# Mandible
terrestrial_conv_Mandible.plot <- plotConv(SC=terrestrial_conv_Mandible.test,y=paca_Mandible$x[,1:4],variable=1,state=terrestrial_vector)
terrestrial_conv_Mandible.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F21A00"),lty=3),
                                      points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F21A00")),
                                      legend.args=list(pch=c(23,21),x="top"))
terrestrial_conv_Mandible.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                     polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                     line.args=list(line.col="#F21A00",lwd=4))
# Suspensorium
terrestrial_conv_Suspensorium.plot <- plotConv(SC=terrestrial_conv_Suspensorium.test,y=paca_Suspensorium$x[,1:4],variable=1,state=terrestrial_vector)
terrestrial_conv_Suspensorium.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F21A00"),lty=3),
                                          points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F21A00")),
                                          legend.args=list(pch=c(23,21),x="top"))
terrestrial_conv_Suspensorium.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                         polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                         line.args=list(line.col="#F21A00",lwd=4))
# Habitat - Fossorial
# Skull
Fossorial_conv.plot <- plotConv(SC=Fossorial_conv.test,y=paca_test$x[,1:4],variable=1,state=Fossorial_vector)
Fossorial_conv.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#EBCC2A"),lty=3),
                                 points.args=list(pch=c(23,21),bg=c("#D9D0D3","#EBCC2A")),
                                 legend.args=list(pch=c(23,21),x="top"))
Fossorial_conv.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                line.args=list(line.col="#EBCC2A",lwd=4))
# Snout
Fossorial_conv_Snout.plot <- plotConv(SC=Fossorial_conv_snout.test,y=paca_Snout$x[,1:4],variable=1,state=Fossorial_vector)
Fossorial_conv_Snout.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#EBCC2A"),lty=3),
                                       points.args=list(pch=c(23,21),bg=c("#D9D0D3","#EBCC2A")),
                                       legend.args=list(pch=c(23,21),x="top"))
Fossorial_conv_Snout.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                      polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                      line.args=list(line.col="#EBCC2A",lwd=4))
# Braincase
Fossorial_conv_Braincase.plot <- plotConv(SC=Fossorial_conv_Braincase.test,y=paca_Braincase$x[,1:4],variable=1,state=Fossorial_vector)
Fossorial_conv_Braincase.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#EBCC2A"),lty=3),
                                           points.args=list(pch=c(23,21),bg=c("#D9D0D3","#EBCC2A")),
                                           legend.args=list(pch=c(23,21),x="top"))
Fossorial_conv_Braincase.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                          polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                          line.args=list(line.col="#EBCC2A",lwd=4))
# Maxilla
Fossorial_conv_Maxilla.plot <- plotConv(SC=Fossorial_conv_Maxilla.test,y=paca_Maxilla$x[,1:4],variable=1,state=Fossorial_vector)
Fossorial_conv_Maxilla.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#EBCC2A"),lty=3),
                                         points.args=list(pch=c(23,21),bg=c("#D9D0D3","#EBCC2A")),
                                         legend.args=list(pch=c(23,21),x="top"))
Fossorial_conv_Maxilla.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                        polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                        line.args=list(line.col="#EBCC2A",lwd=4))
# Palatopterygoid Arch
Fossorial_conv_Palatoptery.plot <- plotConv(SC=Fossorial_conv_Palatoptery.test,y=paca_Palatoptery$x[,1:4],variable=1,state=Fossorial_vector)
Fossorial_conv_Palatoptery.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#EBCC2A"),lty=3),
                                             points.args=list(pch=c(23,21),bg=c("#D9D0D3","#EBCC2A")),
                                             legend.args=list(pch=c(23,21),x="top"))
Fossorial_conv_Palatoptery.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                            polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                            line.args=list(line.col="#EBCC2A",lwd=4))
# Mandible
Fossorial_conv_Mandible.plot <- plotConv(SC=Fossorial_conv_Mandible.test,y=paca_Mandible$x[,1:4],variable=1,state=Fossorial_vector)
Fossorial_conv_Mandible.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#EBCC2A"),lty=3),
                                          points.args=list(pch=c(23,21),bg=c("#D9D0D3","#EBCC2A")),
                                          legend.args=list(pch=c(23,21),x="top"))
Fossorial_conv_Mandible.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                         polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                         line.args=list(line.col="#EBCC2A",lwd=4))
# Suspensorium
Fossorial_conv_Suspensorium.plot <- plotConv(SC=Fossorial_conv_Suspensorium.test,y=paca_Suspensorium$x[,1:4],variable=1,state=Fossorial_vector)
Fossorial_conv_Suspensorium.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#EBCC2A"),lty=3),
                                              points.args=list(pch=c(23,21),bg=c("#D9D0D3","#EBCC2A")),
                                              legend.args=list(pch=c(23,21),x="top"))
Fossorial_conv_Suspensorium.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                             polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                             line.args=list(line.col="#EBCC2A",lwd=4))

# Diet - Generalist
# Skull
Generalist_conv.plot <- plotConv(SC=Generalist_conv.test,y=paca_test$x[,1:4],variable=1,state=Generalist_vector)
Generalist_conv.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#DD3C51"),lty=3),
                             points.args=list(pch=c(23,21),bg=c("#D9D0D3","#DD3C51")),
                             legend.args=list(pch=c(23,21),x="top"))
Generalist_conv.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                            polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                            line.args=list(line.col="#DD3C51",lwd=4))
# Snout
Generalist_conv_Snout.plot <- plotConv(SC=Generalist_conv_snout.test,y=paca_Snout$x[,1:4],variable=1,state=Generalist_vector)
Generalist_conv_Snout.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#DD3C51"),lty=3),
                                     points.args=list(pch=c(23,21),bg=c("#D9D0D3","#DD3C51")),
                                     legend.args=list(pch=c(23,21),x="top"))
Generalist_conv_Snout.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                    polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                    line.args=list(line.col="#DD3C51",lwd=4))
# Braincase
Generalist_conv_Braincase.plot <- plotConv(SC=Generalist_conv_Braincase.test,y=paca_Braincase$x[,1:4],variable=1,state=Generalist_vector)
Generalist_conv_Braincase.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#DD3C51"),lty=3),
                                         points.args=list(pch=c(23,21),bg=c("#D9D0D3","#DD3C51")),
                                         legend.args=list(pch=c(23,21),x="top"))
Generalist_conv_Braincase.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                        polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                        line.args=list(line.col="#DD3C51",lwd=4))
# Maxilla
Generalist_conv_Maxilla.plot <- plotConv(SC=Generalist_conv_Maxilla.test,y=paca_Maxilla$x[,1:4],variable=1,state=Generalist_vector)
Generalist_conv_Maxilla.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#DD3C51"),lty=3),
                                       points.args=list(pch=c(23,21),bg=c("#D9D0D3","#DD3C51")),
                                       legend.args=list(pch=c(23,21),x="top"))
Generalist_conv_Maxilla.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                      polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                      line.args=list(line.col="#DD3C51",lwd=4))
# Palatopterygoid Arch
Generalist_conv_Palatoptery.plot <- plotConv(SC=Generalist_conv_Palatoptery.test,y=paca_Palatoptery$x[,1:4],variable=1,state=Generalist_vector)
Generalist_conv_Palatoptery.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#DD3C51"),lty=3),
                                           points.args=list(pch=c(23,21),bg=c("#D9D0D3","#DD3C51")),
                                           legend.args=list(pch=c(23,21),x="top"))
Generalist_conv_Palatoptery.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                          polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                          line.args=list(line.col="#DD3C51",lwd=4))
# Mandible
Generalist_conv_Mandible.plot <- plotConv(SC=Generalist_conv_Mandible.test,y=paca_Mandible$x[,1:4],variable=1,state=Generalist_vector)
Generalist_conv_Mandible.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#DD3C51"),lty=3),
                                        points.args=list(pch=c(23,21),bg=c("#D9D0D3","#DD3C51")),
                                        legend.args=list(pch=c(23,21),x="top"))
Generalist_conv_Mandible.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                       polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                       line.args=list(line.col="#DD3C51",lwd=4))
# Suspensorium
Generalist_conv_Suspensorium.plot <- plotConv(SC=Generalist_conv_Suspensorium.test,y=paca_Suspensorium$x[,1:4],variable=1,state=Generalist_vector)
Generalist_conv_Suspensorium.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#DD3C51"),lty=3),
                                            points.args=list(pch=c(23,21),bg=c("#D9D0D3","#DD3C51")),
                                            legend.args=list(pch=c(23,21),x="top"))
Generalist_conv_Suspensorium.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                           polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                           line.args=list(line.col="#DD3C51",lwd=4))
# Diet - Aquatic Generalist
# Skull
Aquatic_Generalist_conv.plot <- plotConv(SC=Aquatic_Generalist_conv.test,y=paca_test$x[,1:4],variable=1,state=Aquatic_Generalist_vector)
Aquatic_Generalist_conv.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#46ACC8"),lty=3),
                             points.args=list(pch=c(23,21),bg=c("#D9D0D3","#46ACC8")),
                             legend.args=list(pch=c(23,21),x="top"))
Aquatic_Generalist_conv.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                            polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                            line.args=list(line.col="#46ACC8",lwd=4))
# Snout
Aquatic_Generalist_conv_Snout.plot <- plotConv(SC=Aquatic_Generalist_conv_snout.test,y=paca_Snout$x[,1:4],variable=1,state=Aquatic_Generalist_vector)
Aquatic_Generalist_conv_Snout.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#46ACC8"),lty=3),
                                      points.args=list(pch=c(23,21),bg=c("#D9D0D3","#46ACC8")),
                                      legend.args=list(pch=c(23,21),x="top"))
Aquatic_Generalist_conv_Snout.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                     polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                     line.args=list(line.col="#46ACC8",lwd=4))
# Braincase
Aquatic_Generalist_conv_Braincase.plot <- plotConv(SC=Aquatic_Generalist_conv_Braincase.test,y=paca_Braincase$x[,1:4],variable=1,state=Aquatic_Generalist_vector)
Aquatic_Generalist_conv_Braincase.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#46ACC8"),lty=3),
                                          points.args=list(pch=c(23,21),bg=c("#D9D0D3","#46ACC8")),
                                          legend.args=list(pch=c(23,21),x="top"))
Aquatic_Generalist_conv_Braincase.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                         polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                         line.args=list(line.col="#46ACC8",lwd=4))
# Maxilla
Aquatic_Generalist_conv_Maxilla.plot <- plotConv(SC=Aquatic_Generalist_conv_Maxilla.test,y=paca_Maxilla$x[,1:4],variable=1,state=Aquatic_Generalist_vector)
Aquatic_Generalist_conv_Maxilla.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#46ACC8"),lty=3),
                                        points.args=list(pch=c(23,21),bg=c("#D9D0D3","#46ACC8")),
                                        legend.args=list(pch=c(23,21),x="top"))
Aquatic_Generalist_conv_Maxilla.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                       polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                       line.args=list(line.col="#46ACC8",lwd=4))
# Palatopterygoid Arch
Aquatic_Generalist_conv_Palatoptery.plot <- plotConv(SC=Aquatic_Generalist_conv_Palatoptery.test,y=paca_Palatoptery$x[,1:4],variable=1,state=Aquatic_Generalist_vector)
Aquatic_Generalist_conv_Palatoptery.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#46ACC8"),lty=3),
                                            points.args=list(pch=c(23,21),bg=c("#D9D0D3","#46ACC8")),
                                            legend.args=list(pch=c(23,21),x="top"))
Aquatic_Generalist_conv_Palatoptery.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                           polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                           line.args=list(line.col="#46ACC8",lwd=4))
# Mandible
Aquatic_Generalist_conv_Mandible.plot <- plotConv(SC=Aquatic_Generalist_conv_Mandible.test,y=paca_Mandible$x[,1:4],variable=1,state=Aquatic_Generalist_vector)
Aquatic_Generalist_conv_Mandible.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#46ACC8"),lty=3),
                                         points.args=list(pch=c(23,21),bg=c("#D9D0D3","#46ACC8")),
                                         legend.args=list(pch=c(23,21),x="top"))
Aquatic_Generalist_conv_Mandible.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                        polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                        line.args=list(line.col="#46ACC8",lwd=4))
# Suspensorium
Aquatic_Generalist_conv_Suspensorium.plot <- plotConv(SC=Aquatic_Generalist_conv_Suspensorium.test,y=paca_Suspensorium$x[,1:4],variable=1,state=Aquatic_Generalist_vector)
Aquatic_Generalist_conv_Suspensorium.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#46ACC8"),lty=3),
                                             points.args=list(pch=c(23,21),bg=c("#D9D0D3","#46ACC8")),
                                             legend.args=list(pch=c(23,21),x="top"))
Aquatic_Generalist_conv_Suspensorium.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                            polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                            line.args=list(line.col="#46ACC8",lwd=4))
# Diet - Vermivorous
# Skull
Vermivorous_conv.plot <- plotConv(SC=Vermivorous_conv.test,y=paca_test$x[,1:4],variable=1,state=Vermivorous_vector)
Vermivorous_conv.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F98400"),lty=3),
                             points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F98400")),
                             legend.args=list(pch=c(23,21),x="top"))
Vermivorous_conv.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                            polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                            line.args=list(line.col="#F98400",lwd=4))
# Snout
Vermivorous_conv_Snout.plot <- plotConv(SC=Vermivorous_conv_snout.test,y=paca_Snout$x[,1:4],variable=1,state=Vermivorous_vector)
Vermivorous_conv_Snout.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F98400"),lty=3),
                                              points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F98400")),
                                              legend.args=list(pch=c(23,21),x="top"))
Vermivorous_conv_Snout.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                             polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                             line.args=list(line.col="#F98400",lwd=4))
# Braincase
Vermivorous_conv_Braincase.plot <- plotConv(SC=Vermivorous_conv_Braincase.test,y=paca_Braincase$x[,1:4],variable=1,state=Vermivorous_vector)
Vermivorous_conv_Braincase.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F98400"),lty=3),
                                                  points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F98400")),
                                                  legend.args=list(pch=c(23,21),x="top"))
Vermivorous_conv_Braincase.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                                 polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                                 line.args=list(line.col="#F98400",lwd=4))
# Maxilla
Vermivorous_conv_Maxilla.plot <- plotConv(SC=Vermivorous_conv_Maxilla.test,y=paca_Maxilla$x[,1:4],variable=1,state=Vermivorous_vector)
Vermivorous_conv_Maxilla.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F98400"),lty=3),
                                                points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F98400")),
                                                legend.args=list(pch=c(23,21),x="top"))
Vermivorous_conv_Maxilla.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                               polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                               line.args=list(line.col="#F98400",lwd=4))
# Palatopterygoid Arch
Vermivorous_conv_Palatoptery.plot <- plotConv(SC=Vermivorous_conv_Palatoptery.test,y=paca_Palatoptery$x[,1:4],variable=1,state=Vermivorous_vector)
Vermivorous_conv_Palatoptery.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F98400"),lty=3),
                                                    points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F98400")),
                                                    legend.args=list(pch=c(23,21),x="top"))
Vermivorous_conv_Palatoptery.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                                   polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                                   line.args=list(line.col="#F98400",lwd=4))
# Mandible
Vermivorous_conv_Mandible.plot <- plotConv(SC=Vermivorous_conv_Mandible.test,y=paca_Mandible$x[,1:4],variable=1,state=Vermivorous_vector)
Vermivorous_conv_Mandible.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F98400"),lty=3),
                                                 points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F98400")),
                                                 legend.args=list(pch=c(23,21),x="top"))
Vermivorous_conv_Mandible.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                                polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                                line.args=list(line.col="#F98400",lwd=4))
# Suspensorium
Vermivorous_conv_Suspensorium.plot <- plotConv(SC=Vermivorous_conv_Suspensorium.test,y=paca_Suspensorium$x[,1:4],variable=1,state=Vermivorous_vector)
Vermivorous_conv_Suspensorium.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#F98400"),lty=3),
                                                     points.args=list(pch=c(23,21),bg=c("#D9D0D3","#F98400")),
                                                     legend.args=list(pch=c(23,21),x="top"))
Vermivorous_conv_Suspensorium.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                                    polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                                    line.args=list(line.col="#F98400",lwd=4))

# Diet - Durophagous
# Skull
Durophagous_conv.plot <- plotConv(SC=Durophagous_conv.test,y=paca_test$x[,1:4],variable=1,state=Durophagous_vector)
Durophagous_conv.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#33A02C"),lty=3),
                             points.args=list(pch=c(23,21),bg=c("#D9D0D3","#33A02C")),
                             legend.args=list(pch=c(23,21),x="top"))
Durophagous_conv.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                            polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                            line.args=list(line.col="#33A02C",lwd=4))
# Snout
Durophagous_conv_Snout.plot <- plotConv(SC=Durophagous_conv_snout.test,y=paca_Snout$x[,1:4],variable=1,state=Durophagous_vector)
Durophagous_conv_Snout.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#33A02C"),lty=3),
                                       points.args=list(pch=c(23,21),bg=c("#D9D0D3","#33A02C")),
                                       legend.args=list(pch=c(23,21),x="top"))
Durophagous_conv_Snout.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                      polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                      line.args=list(line.col="#33A02C",lwd=4))
# Braincase
Durophagous_conv_Braincase.plot <- plotConv(SC=Durophagous_conv_Braincase.test,y=paca_Braincase$x[,1:4],variable=1,state=Durophagous_vector)
Durophagous_conv_Braincase.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#33A02C"),lty=3),
                                           points.args=list(pch=c(23,21),bg=c("#D9D0D3","#33A02C")),
                                           legend.args=list(pch=c(23,21),x="top"))
Durophagous_conv_Braincase.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                          polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                          line.args=list(line.col="#33A02C",lwd=4))
# Maxilla
Durophagous_conv_Maxilla.plot <- plotConv(SC=Durophagous_conv_Maxilla.test,y=paca_Maxilla$x[,1:4],variable=1,state=Durophagous_vector)
Durophagous_conv_Maxilla.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#33A02C"),lty=3),
                                         points.args=list(pch=c(23,21),bg=c("#D9D0D3","#33A02C")),
                                         legend.args=list(pch=c(23,21),x="top"))
Durophagous_conv_Maxilla.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                        polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                        line.args=list(line.col="#33A02C",lwd=4))
# Palatopterygoid Arch
Durophagous_conv_Palatoptery.plot <- plotConv(SC=Durophagous_conv_Palatoptery.test,y=paca_Palatoptery$x[,1:4],variable=1,state=Durophagous_vector)
Durophagous_conv_Palatoptery.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#33A02C"),lty=3),
                                             points.args=list(pch=c(23,21),bg=c("#D9D0D3","#33A02C")),
                                             legend.args=list(pch=c(23,21),x="top"))
Durophagous_conv_Palatoptery.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                            polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                            line.args=list(line.col="#33A02C",lwd=4))
# Mandible
Durophagous_conv_Mandible.plot <- plotConv(SC=Durophagous_conv_Mandible.test,y=paca_Mandible$x[,1:4],variable=1,state=Durophagous_vector)
Durophagous_conv_Mandible.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#33A02C"),lty=3),
                                          points.args=list(pch=c(23,21),bg=c("#D9D0D3","#33A02C")),
                                          legend.args=list(pch=c(23,21),x="top"))
Durophagous_conv_Mandible.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                         polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                         line.args=list(line.col="#33A02C",lwd=4))
# Suspensorium
Durophagous_conv_Suspensorium.plot <- plotConv(SC=Durophagous_conv_Suspensorium.test,y=paca_Suspensorium$x[,1:4],variable=1,state=Durophagous_vector)
Durophagous_conv_Suspensorium.plot$plotPChull(chull.args=list(border=c("#D9D0D3","#33A02C"),lty=3),
                                              points.args=list(pch=c(23,21),bg=c("#D9D0D3","#33A02C")),
                                              legend.args=list(pch=c(23,21),x="top"))
Durophagous_conv_Suspensorium.plot$plotPolar(polar.args=list(clockwise=T,start=0,rad.col="black",grid.col="black"),
                                             polygon.args=list(line.col="#D9D0D3",poly.col="#D9D0D3",lwd=3),
                                             line.args=list(line.col="#33A02C",lwd=4))

# Test for convergent evolution using convevol
# Habitat - Aquatic
aquatic_conv.sig <- convSig(All_Tree, paca_test$x[,1:4], Aquatic_Taxa, nsim = 1000)
aquatic_conv_snout.sig <- convSig(All_Tree, paca_Snout$x[,1:4], Aquatic_Taxa, nsim = 1000)
aquatic_conv_Braincase.sig <- convSig(All_Tree, paca_Braincase$x[,1:4], Aquatic_Taxa, nsim = 1000)
aquatic_conv_Maxilla.sig <- convSig(All_Tree, paca_Maxilla$x[,1:4], Aquatic_Taxa, nsim = 1000)
aquatic_conv_Palatoptery.sig <- convSig(All_Tree, paca_Palatoptery$x[,1:4], Aquatic_Taxa, nsim = 1000)
aquatic_conv_Mandible.sig <- convSig(All_Tree, paca_Mandible$x[,1:4], Aquatic_Taxa, nsim = 1000)
aquatic_conv_Suspensorium.sig <- convSig(All_Tree, paca_Suspensorium$x[,1:4], Aquatic_Taxa, nsim = 1000)
aquatic_convs.sig <- rbind(aquatic_conv.sig,aquatic_conv_snout.sig,aquatic_conv_Braincase.sig,aquatic_conv_Maxilla.sig,aquatic_conv_Palatoptery.sig,aquatic_conv_Mandible.sig,aquatic_conv_Suspensorium.sig)
write.csv(aquatic_convs.sig, file = "Convergence/Tables/Aquatic_Habitat_ConvSig.csv")

# Habitat - Terrestrial
Terrestrial_conv.sig <- convSig(All_Tree, paca_test$x[,1:4], Terrestrial_Taxa, nsim = 1000)
Terrestrial_conv_snout.sig <- convSig(All_Tree, paca_Snout$x[,1:4], Terrestrial_Taxa, nsim = 1000)
Terrestrial_conv_Braincase.sig <- convSig(All_Tree, paca_Braincase$x[,1:4], Terrestrial_Taxa, nsim = 1000)
Terrestrial_conv_Maxilla.sig <- convSig(All_Tree, paca_Maxilla$x[,1:4], Terrestrial_Taxa, nsim = 1000)
Terrestrial_conv_Palatoptery.sig <- convSig(All_Tree, paca_Palatoptery$x[,1:4], Terrestrial_Taxa, nsim = 1000)
Terrestrial_conv_Mandible.sig <- convSig(All_Tree, paca_Mandible$x[,1:4], Terrestrial_Taxa, nsim = 1000)
Terrestrial_conv_Suspensorium.sig <- convSig(All_Tree, paca_Suspensorium$x[,1:4], Terrestrial_Taxa, nsim = 1000)
Terrestrial_convs.sig <- rbind(Terrestrial_conv.sig,Terrestrial_conv_snout.sig,Terrestrial_conv_Braincase.sig,Terrestrial_conv_Maxilla.sig,Terrestrial_conv_Palatoptery.sig,Terrestrial_conv_Mandible.sig,Terrestrial_conv_Suspensorium.sig)
write.csv(Terrestrial_convs.sig, file = "Convergence/Tables/Terrestrial_Habitat_ConvSig.csv")

# Habitat - Fossorial
Fossorial_conv.sig <- convSig(All_Tree, paca_test$x[,1:4], Fossorial_Taxa, nsim = 1000)
Fossorial_conv_snout.sig <- convSig(All_Tree, paca_Snout$x[,1:4], Fossorial_Taxa, nsim = 1000)
Fossorial_conv_Braincase.sig <- convSig(All_Tree, paca_Braincase$x[,1:4], Fossorial_Taxa, nsim = 1000)
Fossorial_conv_Maxilla.sig <- convSig(All_Tree, paca_Maxilla$x[,1:4], Fossorial_Taxa, nsim = 1000)
Fossorial_conv_Palatoptery.sig <- convSig(All_Tree, paca_Palatoptery$x[,1:4], Fossorial_Taxa, nsim = 1000)
Fossorial_conv_Mandible.sig <- convSig(All_Tree, paca_Mandible$x[,1:4], Fossorial_Taxa, nsim = 1000)
Fossorial_conv_Suspensorium.sig <- convSig(All_Tree, paca_Suspensorium$x[,1:4], Fossorial_Taxa, nsim = 1000)
Fossorial_convs.sig <- rbind(Fossorial_conv.sig,Fossorial_conv_snout.sig,Fossorial_conv_Braincase.sig,Fossorial_conv_Maxilla.sig,Fossorial_conv_Palatoptery.sig,Fossorial_conv_Mandible.sig,Fossorial_conv_Suspensorium.sig)
write.csv(Fossorial_convs.sig, file = "Convergence/Tables/Fossorial_Habitat_ConvSig.csv")

# Diet - Generalist
Generalist_conv.sig <- convSig(All_Tree, paca_test$x[,1:4], Generalist_Taxa, nsim = 1000)
Generalist_conv_snout.sig <- convSig(All_Tree, paca_Snout$x[,1:4], Generalist_Taxa, nsim = 1000)
Generalist_conv_Braincase.sig <- convSig(All_Tree, paca_Braincase$x[,1:4], Generalist_Taxa, nsim = 1000)
Generalist_conv_Maxilla.sig <- convSig(All_Tree, paca_Maxilla$x[,1:4], Generalist_Taxa, nsim = 1000)
Generalist_conv_Palatoptery.sig <- convSig(All_Tree, paca_Palatoptery$x[,1:4], Generalist_Taxa, nsim = 1000)
Generalist_conv_Mandible.sig <- convSig(All_Tree, paca_Mandible$x[,1:4], Generalist_Taxa, nsim = 1000)
Generalist_conv_Suspensorium.sig <- convSig(All_Tree, paca_Suspensorium$x[,1:4], Generalist_Taxa, nsim = 1000)
Generalist_convs.sig <- rbind(Generalist_conv.sig,Generalist_conv_snout.sig,Generalist_conv_Braincase.sig,Generalist_conv_Maxilla.sig,Generalist_conv_Palatoptery.sig,Generalist_conv_Mandible.sig,Generalist_conv_Suspensorium.sig)
write.csv(Generalist_convs.sig, file = "Convergence/Tables/Generalist_Diet_ConvSig.csv")

# Diet - Aquatic Generalist
Aquatic_Generalist_conv.sig <- convSig(All_Tree, paca_test$x[,1:4], Aquatic_Generalist_Taxa, nsim = 1000)
Aquatic_Generalist_conv_snout.sig <- convSig(All_Tree, paca_Snout$x[,1:4], Aquatic_Generalist_Taxa, nsim = 1000)
Aquatic_Generalist_conv_Braincase.sig <- convSig(All_Tree, paca_Braincase$x[,1:4], Aquatic_Generalist_Taxa, nsim = 1000)
Aquatic_Generalist_conv_Maxilla.sig <- convSig(All_Tree, paca_Maxilla$x[,1:4], Aquatic_Generalist_Taxa, nsim = 1000)
Aquatic_Generalist_conv_Palatoptery.sig <- convSig(All_Tree, paca_Palatoptery$x[,1:4], Aquatic_Generalist_Taxa, nsim = 1000)
Aquatic_Generalist_conv_Mandible.sig <- convSig(All_Tree, paca_Mandible$x[,1:4], Aquatic_Generalist_Taxa, nsim = 1000)
Aquatic_Generalist_conv_Suspensorium.sig <- convSig(All_Tree, paca_Suspensorium$x[,1:4], Aquatic_Generalist_Taxa, nsim = 1000)
Aquatic_Generalist_convs.sig <- rbind(Aquatic_Generalist_conv.sig,Aquatic_Generalist_conv_snout.sig,Aquatic_Generalist_conv_Braincase.sig,Aquatic_Generalist_conv_Maxilla.sig,Aquatic_Generalist_conv_Palatoptery.sig,Aquatic_Generalist_conv_Mandible.sig,Aquatic_Generalist_conv_Suspensorium.sig)
write.csv(Aquatic_Generalist_convs.sig, file = "Convergence/Tables/Aquatic_Generalist_Diet_ConvSig.csv")

# Diet - Vermivorous
Vermivorous_conv.sig <- convSig(All_Tree, paca_test$x[,1:4], Vermivorous_Taxa, nsim = 1000)
Vermivorous_conv_snout.sig <- convSig(All_Tree, paca_Snout$x[,1:4], Vermivorous_Taxa, nsim = 1000)
Vermivorous_conv_Braincase.sig <- convSig(All_Tree, paca_Braincase$x[,1:4], Vermivorous_Taxa, nsim = 1000)
Vermivorous_conv_Maxilla.sig <- convSig(All_Tree, paca_Maxilla$x[,1:4], Vermivorous_Taxa, nsim = 1000)
Vermivorous_conv_Palatoptery.sig <- convSig(All_Tree, paca_Palatoptery$x[,1:4], Vermivorous_Taxa, nsim = 1000)
Vermivorous_conv_Mandible.sig <- convSig(All_Tree, paca_Mandible$x[,1:4], Vermivorous_Taxa, nsim = 1000)
Vermivorous_conv_Suspensorium.sig <- convSig(All_Tree, paca_Suspensorium$x[,1:4], Vermivorous_Taxa, nsim = 1000)
Vermivorous_convs.sig <- rbind(Vermivorous_conv.sig,Vermivorous_conv_snout.sig,Vermivorous_conv_Braincase.sig,Vermivorous_conv_Maxilla.sig,Vermivorous_conv_Palatoptery.sig,Vermivorous_conv_Mandible.sig,Vermivorous_conv_Suspensorium.sig)
write.csv(Vermivorous_convs.sig, file = "Convergence/Tables/Vermivorous_Diet_ConvSig.csv")

# Diet - Durophagous
Durophagous_conv.sig <- convSig(All_Tree, paca_test$x[,1:4], Durophagous_Taxa, nsim = 1000)
Durophagous_conv_snout.sig <- convSig(All_Tree, paca_Snout$x[,1:4], Durophagous_Taxa, nsim = 1000)
Durophagous_conv_Braincase.sig <- convSig(All_Tree, paca_Braincase$x[,1:4], Durophagous_Taxa, nsim = 1000)
Durophagous_conv_Maxilla.sig <- convSig(All_Tree, paca_Maxilla$x[,1:4], Durophagous_Taxa, nsim = 1000)
Durophagous_conv_Palatoptery.sig <- convSig(All_Tree, paca_Palatoptery$x[,1:4], Durophagous_Taxa, nsim = 1000)
Durophagous_conv_Mandible.sig <- convSig(All_Tree, paca_Mandible$x[,1:4], Durophagous_Taxa, nsim = 1000)
Durophagous_conv_Suspensorium.sig <- convSig(All_Tree, paca_Suspensorium$x[,1:4], Durophagous_Taxa, nsim = 1000)
Durophagous_convs.sig <- rbind(Durophagous_conv.sig,Durophagous_conv_snout.sig,Durophagous_conv_Braincase.sig,Durophagous_conv_Maxilla.sig,Durophagous_conv_Palatoptery.sig,Durophagous_conv_Mandible.sig,Durophagous_conv_Suspensorium.sig)
write.csv(Durophagous_convs.sig, file = "Convergence/Tables/Durophagous_Diet_ConvSig.csv")


clade <- extractVector(Thamno, "Clade")
genus <- extractVector(Thamno, "Genus")
gdf <- geomorph.data.frame(coords = Skull_Test, species = All_Tree$tip.label, diet = diet, habitat = habitat, genus = genus, clade = clade)
morphol.disparity(coords ~ 1, groups= ~ clade, data = gdf, iter = 999, print.progress = T)
morphol.disparity(coords ~ 1, groups= genus, data = gdf, iter = 999, print.progress = T)


clade_result <- dispRity.per.group(paca_test$x, list(semifossorial = c(1:3,16,50:51),
                                                 nerodia = c(4:12),
                                                 regina = c(13:14),
                                                 storeria = c(16:19),
                                                 ws_thamno = c(20:22,24,26:27,30,32:33,36:37,39,41,44:46),
                                                 mex_thamno = c(23,25,28:29,31,34:35,38,40,42:43,47:48)))
summary(clade_result) ; plot(clade_result)

genus_result <- dispRity.per.group(paca_test$x, list(Liodytes = c(2:3),
                                                     Nerodia = c(4:12),
                                                     Regina = c(13:14),
                                                     Storeria = c(16:19),
                                                     Thamnophis = c(20:48),
                                                     Virginia = c(50:51)))
summary(genus_result) ; plot(genus_result)

region_result <- dispRity.per.group(paca_test$x, list(east = c(1:17,19,21:22,37,41,44:45,49:51),
                                                     west = c(18,20,23:36,38:40,42:43,46:48)))
summary(region_result) ; plot(region_result)


color3<-colorRampPalette(c("#0c2c84","#225ea8","#31a354","#ffff00","#fe9929","#fc4e2a","red","darkred"))


color3 <- wes_palette("AsteroidCity1")
Thamno_Tree <- readNexus("BayesTraits/Data/Thamnophiini_WGS_51_Taxa.tree")
tree_BTraits_PACA<-BTRTools::rjpp(rjlog = "BayesTraits/Results/Skull/Full_Skull_PACA.txt.VarRates.txt",
                                  rjtrees = "BayesTraits/Results/Skull/Full_Skull_PACA.txt.Output.trees",
                                  tree = Thamno_Tree) #this is your time scaled tree that was used to input into bayestraits
eco <- read.csv("Thamno_Discrete_Eco_Traits.csv")
Thamno_Eco <- as.treedata.table(Thamno_Tree, eco, name_column = "Name_in_Tree")
eco <- as.data.frame(Thamno_Eco$dat)

# Read the reference tree
Thamno_Tree <- readNexus("BayesTraits/Data/Thamnophiini_WGS_51_Taxa.tree")

# Read the multiphylo object
inconsistent_trees <- read.nexus("BayesTraits/Results/Skull/Full_Skull_PACA.txt.Output.trees")

# Check if all trees have the same number of tips as the reference tree
if (!all(sapply(inconsistent_trees, function(tree) length(tree$tip.label) == length(Thamno_Tree$tip.label)))) {
  stop("Not all trees have the same number of tips as the reference tree.")
}

# Ensure consistency of tip labels for each tree
for(i in 1:length(inconsistent_trees)) {
  inconsistent_trees[[i]]$tip.label <- Thamno_Tree$tip.label
}


# Now use the modified tree for analysis
tree_BTraits_PACA <- BTRTools::rjpp(rjlog = "BayesTraits/Results/Skull/Full_Skull_PACA.txt.VarRates.txt",
                                    rjtrees = inconsistent_tree,
                                    tree = Thamno_Tree)


add_rjpp_to_tree <- function(rjpp_out){
  rjpp_data <- as_tibble(rjpp_out$data)
  timetree <- rjpp_out$meantree
  timetree$edge.length <- rjpp_data$orgBL[-1]
  timetree$root.time <- max(nodeHeights(timetree))
  rjpp_data_nodes <- rjpp_data %>% rename(., node=descNode) %>% mutate(., iters = rjpp_out$niter) %>% mutate(., ppRate = round(nOrgnNRate/iters,2))
  timetree <- treeio::as.treedata(timetree)
  treedata <- treeio::full_join(timetree, y = rjpp_data_nodes, by = "node")
  return(treedata)
}


habitat_cols <- c("#F21A00","#3B9AB2","#EBCC2A")
diet_cols <- c("#DD3C51","#46ACC8","#F98400","#33A02C")
geo_cols <- c("#10A494","#DE524A")
paca_tree_w_data <- add_rjpp_to_tree(tree_BTraits_PACA)
paca_tree_w_data <- full_join(paca_tree_w_data, mutate(eco, label = eco$tip.label), by = "label")
threshold <- .6 # the minimum posterior probability you want to plot a symbol for 
p<-ggtree(paca_tree_w_data, aes(color = log(meanRate)), size=2)+
  scale_colour_gradientn(colours = color3(100))+
  #geom_tippoint(aes(fill=Diet, x=x+5),color="black", shape=24)+
  #scale_fill_manual(breaks = c("Generalist","Aquatic_Generalist","Vermivorous","Durophagous"), values = diet_cols)+
  #geom_tippoint(aes(fill=Geography, x=x+3),color="black", shape=22)+
  #scale_fill_manual(breaks = c("East","West"), values = geo_cols)+
  theme(legend.position="top")+
  theme(legend.position=c(.10,.90),legend.direction = "horizontal",legend.box.background = element_rect(colour = "black",size =1))+
  scale_size(range = c(1,2))+ 
  labs(title="Full Skull Shape",
       color="log(Rate)")+
  geom_nodepoint(aes(subset=ppRate>threshold, size = ppRate),color='black',fill="grey", shape=24)+
  scale_size(range = c(1,2))+ 
  geom_tiplab(label= sub("_", " ",paca_tree_w_data@phylo$tip.label), size=3, color = "black", fontface="italic")+
  coord_cartesian(#xlim = c(0, 25), #you have to fiddle with these values to get your tip labels to show. the first value should be just before your root time, second value pads out space for tip labels
  ylim = c(-2,  52)) + #first value makes room for geo timescale, second value is vertical space and should be a few more than your number of tips
  #expand = FALSE) +
  scale_x_continuous(breaks=-epochs$max_age[c(1:5)], labels=epochs$max_age[c(1:5)]) +
  theme(panel.grid.major.x = element_line(colour="grey", size=0.5), legend.key.height =unit(.4,"cm"))#should also be modified based on your time scale limits
p <- revts(p);p
ptree1 <-  gggeo_scale(p, dat = "epochs", neg = TRUE, center_end_labels = TRUE, abbrv = FALSE, skip = c("Holocene"), height = unit(1, "line"), size=3)
pdffn = "BayesTraits/Results/Skull/Figures/Full_Skull_Shape_Tips_Var_Sites.pdf"
pdf(file=pdffn, width=16, height=14)
ptree1
dev.off()


tree_BTraits_PACA_braincase<-BTRTools::rjpp(rjlog = "BayesTraits/Results/Braincase/thamno_braincase_PACA.txt.VarRates.txt",
                                  rjtrees = "BayesTraits/Results/Braincase/thamno_braincase_PACA.txt.Output.trees",
                                  tree = Thamno_Tree) #this is your time scaled tree that was used to input into bayestraits
paca_tree_w_braincase <- add_rjpp_to_tree(tree_BTraits_PACA_braincase)
paca_tree_w_braincase <- full_join(paca_tree_w_braincase, mutate(eco, label = eco$tip.label), by = "label")
threshold <- .6 # the minimum posterior probability you want to plot a symbol for 
p<-ggtree(paca_tree_w_braincase, aes(color = log(meanRate)), size=2)+
  scale_colour_gradientn(colours = color3(100))+
  theme(legend.position="top")+
  theme(legend.position=c(.10,.90),legend.direction = "horizontal",legend.box.background = element_rect(colour = "black",size =1))+
  scale_size(range = c(1,2))+ 
  labs(title="Braincase Shape",
       color="log(Rate)")+
  geom_nodepoint(aes(subset=ppRate>threshold, size = ppRate),color='black',fill="grey", shape=24)+
  scale_size(range = c(1,2))+ 
  geom_tiplab(label= sub("_", " ",paca_tree_w_braincase@phylo$tip.label), size=3, color = "black", fontface="italic")+
  coord_cartesian(ylim = c(-2,  52)) +
  scale_x_continuous(breaks=-epochs$max_age[c(1:5)], labels=epochs$max_age[c(1:5)]) +
  theme(panel.grid.major.x = element_line(colour="grey", size=0.5), legend.key.height =unit(.4,"cm"))#should also be modified based on your time scale limits
p <- revts(p);p
ptree1 <-  gggeo_scale(p, dat = "epochs", neg = TRUE, center_end_labels = TRUE, abbrv = FALSE, skip = c("Holocene"), height = unit(1, "line"), size=3)
pdffn = "BayesTraits/Results/Braincase/Figures/Full_Skull_Shape_Tips_Var_Sites.pdf"
pdf(file=pdffn, width=16, height=14)
ptree1
dev.off()

tree_BTraits_PACA_mandible<-BTRTools::rjpp(rjlog = "BayesTraits/Results/Mandible/thamno_mandible_PACA.txt.VarRates.txt",
                                            rjtrees = "BayesTraits/Results/Mandible/thamno_mandible_PACA.txt.Output.trees",
                                            tree = Thamno_Tree) #this is your time scaled tree that was used to input into bayestraits
paca_tree_w_mandible <- add_rjpp_to_tree(tree_BTraits_PACA_mandible)
paca_tree_w_mandible <- full_join(paca_tree_w_mandible, mutate(eco, label = eco$tip.label), by = "label")
threshold <- .6 # the minimum posterior probability you want to plot a symbol for 
p<-ggtree(paca_tree_w_mandible, aes(color = log(meanRate)), size=2)+
  scale_colour_gradientn(colours = color3(100))+
  theme(legend.position="top")+
  theme(legend.position=c(.10,.90),legend.direction = "horizontal",legend.box.background = element_rect(colour = "black",size =1))+
  scale_size(range = c(1,2))+ 
  labs(title="Mandible Shape",
       color="log(Rate)")+
  geom_nodepoint(aes(subset=ppRate>threshold, size = ppRate),color='black',fill="grey", shape=24)+
  scale_size(range = c(1,2))+ 
  geom_tiplab(label= sub("_", " ",paca_tree_w_mandible@phylo$tip.label), size=3, color = "black", fontface="italic")+
  coord_cartesian(ylim = c(-2,  52)) +
  scale_x_continuous(breaks=-epochs$max_age[c(1:5)], labels=epochs$max_age[c(1:5)]) +
  theme(panel.grid.major.x = element_line(colour="grey", size=0.5), legend.key.height =unit(.4,"cm"))#should also be modified based on your time scale limits
p <- revts(p);p
ptree1 <-  gggeo_scale(p, dat = "epochs", neg = TRUE, center_end_labels = TRUE, abbrv = FALSE, skip = c("Holocene"), height = unit(1, "line"), size=3)
pdffn = "BayesTraits/Results/Mandible/Figures/Mandible_Shape_Tips_Var_Sites.pdf"
pdf(file=pdffn, width=16, height=14)
ptree1
dev.off()

tree_BTraits_PACA_maxilla<-BTRTools::rjpp(rjlog = "BayesTraits/Results/Maxilla/thamno_maxilla_PACA.txt.VarRates.txt",
                                           rjtrees = "BayesTraits/Results/Maxilla/thamno_maxilla_PACA.txt.Output.trees",
                                           tree = Thamno_Tree) #this is your time scaled tree that was used to input into bayestraits
paca_tree_w_maxilla <- add_rjpp_to_tree(tree_BTraits_PACA_maxilla)
paca_tree_w_maxilla <- full_join(paca_tree_w_maxilla, mutate(eco, label = eco$tip.label), by = "label")
threshold <- .6 # the minimum posterior probability you want to plot a symbol for 
p<-ggtree(paca_tree_w_maxilla, aes(color = log(meanRate)), size=2)+
  scale_colour_gradientn(colours = color3(100))+
  theme(legend.position="top")+
  theme(legend.position=c(.10,.90),legend.direction = "horizontal",legend.box.background = element_rect(colour = "black",size =1))+
  scale_size(range = c(1,2))+ 
  labs(title="Maxilla Shape",
       color="log(Rate)")+
  geom_nodepoint(aes(subset=ppRate>threshold, size = ppRate),color='black',fill="grey", shape=24)+
  scale_size(range = c(1,2))+ 
  geom_tiplab(label= sub("_", " ",paca_tree_w_maxilla@phylo$tip.label), size=3, color = "black", fontface="italic")+
  coord_cartesian(ylim = c(-2,  52)) +
  scale_x_continuous(breaks=-epochs$max_age[c(1:5)], labels=epochs$max_age[c(1:5)]) +
  theme(panel.grid.major.x = element_line(colour="grey", size=0.5), legend.key.height =unit(.4,"cm"))#should also be modified based on your time scale limits
p <- revts(p);p
ptree1 <-  gggeo_scale(p, dat = "epochs", neg = TRUE, center_end_labels = TRUE, abbrv = FALSE, skip = c("Holocene"), height = unit(1, "line"), size=3)
pdffn = "BayesTraits/Results/Maxilla/Figures/Maxilla_Shape_Tips_Var_Sites.pdf"
pdf(file=pdffn, width=16, height=14)
ptree1
dev.off()

tree_BTraits_PACA_palatoptery<-BTRTools::rjpp(rjlog = "BayesTraits/Results/Palatopterygoid/thamno_palatoptery_PACA.txt.VarRates.txt",
                                          rjtrees = "BayesTraits/Results/Palatopterygoid/thamno_palatoptery_PACA.txt.Output.trees",
                                          tree = Thamno_Tree) #this is your time scaled tree that was used to input into bayestraits
paca_tree_w_palatoptery <- add_rjpp_to_tree(tree_BTraits_PACA_palatoptery)
paca_tree_w_palatoptery <- full_join(paca_tree_w_palatoptery, mutate(eco, label = eco$tip.label), by = "label")
threshold <- .6 # the minimum posterior probability you want to plot a symbol for 
p<-ggtree(paca_tree_w_palatoptery, aes(color = log(meanRate)), size=2)+
  scale_colour_gradientn(colours = color3(100))+
  theme(legend.position="top")+
  theme(legend.position=c(.10,.90),legend.direction = "horizontal",legend.box.background = element_rect(colour = "black",size =1))+
  scale_size(range = c(1,2))+ 
  labs(title="Palatopterygoid Shape",
       color="log(Rate)")+
  geom_nodepoint(aes(subset=ppRate>threshold, size = ppRate),color='black',fill="grey", shape=24)+
  scale_size(range = c(1,2))+ 
  geom_tiplab(label= sub("_", " ",paca_tree_w_palatoptery@phylo$tip.label), size=3, color = "black", fontface="italic")+
  coord_cartesian(ylim = c(-2,  52)) +
  scale_x_continuous(breaks=-epochs$max_age[c(1:5)], labels=epochs$max_age[c(1:5)]) +
  theme(panel.grid.major.x = element_line(colour="grey", size=0.5), legend.key.height =unit(.4,"cm"))#should also be modified based on your time scale limits
p <- revts(p);p
ptree1 <-  gggeo_scale(p, dat = "epochs", neg = TRUE, center_end_labels = TRUE, abbrv = FALSE, skip = c("Holocene"), height = unit(1, "line"), size=3)
pdffn = "BayesTraits/Results/Palatopterygoid/Figures/Palatopterygoid_Shape_Tips_Var_Sites.pdf"
pdf(file=pdffn, width=16, height=14)
ptree1
dev.off()

tree_BTraits_PACA_snout<-BTRTools::rjpp(rjlog = "BayesTraits/Results/Snout/thamno_snout_PACA.txt.VarRates.txt",
                                              rjtrees = "BayesTraits/Results/Snout/thamno_snout_PACA.txt.Output.trees",
                                              tree = Thamno_Tree) #this is your time scaled tree that was used to input into bayestraits
paca_tree_w_snout <- add_rjpp_to_tree(tree_BTraits_PACA_snout)
paca_tree_w_snout <- full_join(paca_tree_w_snout, mutate(eco, label = eco$tip.label), by = "label")
threshold <- .6 # the minimum posterior probability you want to plot a symbol for 
p<-ggtree(paca_tree_w_snout, aes(color = log(meanRate)), size=2)+
  scale_colour_gradientn(colours = color3(100))+
  theme(legend.position="top")+
  theme(legend.position=c(.10,.90),legend.direction = "horizontal",legend.box.background = element_rect(colour = "black",size =1))+
  scale_size(range = c(1,2))+ 
  labs(title="Snout Shape",
       color="log(Rate)")+
  geom_nodepoint(aes(subset=ppRate>threshold, size = ppRate),color='black',fill="grey", shape=24)+
  scale_size(range = c(1,2))+ 
  geom_tiplab(label= sub("_", " ",paca_tree_w_snout@phylo$tip.label), size=3, color = "black", fontface="italic")+
  coord_cartesian(ylim = c(-2,  52)) +
  scale_x_continuous(breaks=-epochs$max_age[c(1:5)], labels=epochs$max_age[c(1:5)]) +
  theme(panel.grid.major.x = element_line(colour="grey", size=0.5), legend.key.height =unit(.4,"cm"))#should also be modified based on your time scale limits
p <- revts(p);p
ptree1 <-  gggeo_scale(p, dat = "epochs", neg = TRUE, center_end_labels = TRUE, abbrv = FALSE, skip = c("Holocene"), height = unit(1, "line"), size=3)
pdffn = "BayesTraits/Results/Snout/Figures/Snout_Shape_Tips_Var_Sites.pdf"
pdf(file=pdffn, width=16, height=14)
ptree1
dev.off()

tree_BTraits_PACA_suspensorium<-BTRTools::rjpp(rjlog = "BayesTraits/Results/Suspensorium/thamno_suspensorium_PACA.txt.VarRates.txt",
                                        rjtrees = "BayesTraits/Results/Suspensorium/thamno_suspensorium_PACA.txt.Output.trees",
                                        tree = Thamno_Tree) #this is your time scaled tree that was used to input into bayestraits
paca_tree_w_suspensorium <- add_rjpp_to_tree(tree_BTraits_PACA_suspensorium)
paca_tree_w_suspensorium <- full_join(paca_tree_w_suspensorium, mutate(eco, label = eco$tip.label), by = "label")
threshold <- .6 # the minimum posterior probability you want to plot a symbol for 
p<-ggtree(paca_tree_w_suspensorium, aes(color = log(meanRate)), size=2)+
  scale_colour_gradientn(colours = color3(100))+
  theme(legend.position="top")+
  theme(legend.position=c(.10,.90),legend.direction = "horizontal",legend.box.background = element_rect(colour = "black",size =1))+
  scale_size(range = c(1,2))+ 
  labs(title="Suspensorium Shape",
       color="log(Rate)")+
  geom_nodepoint(aes(subset=ppRate>threshold, size = ppRate),color='black',fill="grey", shape=24)+
  scale_size(range = c(1,2))+ 
  geom_tiplab(label= sub("_", " ",paca_tree_w_suspensorium@phylo$tip.label), size=3, color = "black", fontface="italic")+
  coord_cartesian(ylim = c(-2,  52)) +
  scale_x_continuous(breaks=-epochs$max_age[c(1:5)], labels=epochs$max_age[c(1:5)]) +
  theme(panel.grid.major.x = element_line(colour="grey", size=0.5), legend.key.height =unit(.4,"cm"))#should also be modified based on your time scale limits
p <- revts(p);p
ptree1 <-  gggeo_scale(p, dat = "epochs", neg = TRUE, center_end_labels = TRUE, abbrv = FALSE, skip = c("Holocene"), height = unit(1, "line"), size=3)
pdffn = "BayesTraits/Results/Suspensorium/Figures/suspensorium_Shape_Tips_Var_Sites.pdf"
pdf(file=pdffn, width=16, height=14)
ptree1
dev.off()
