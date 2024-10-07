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

source("Analyses_Scripts/Supplementary_Functions/read.markups.json.R")
meta_lms <- read.csv("Data/Landmarks_Bilat.csv")
lms_names <- meta_lms$Name

Female_coords <- readland.tps("Data/Female_Thamnophiini_Ana_Curves_UCE_CatNum.tps",specID = c("ID"))
Male_coords <- readland.tps("Data/Male_Thamnophiini_Ana_Curves_UCE_CatNum.tps",specID = c("ID"))
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
save(slid_coords, file = "Data/slid_coords.RData")
