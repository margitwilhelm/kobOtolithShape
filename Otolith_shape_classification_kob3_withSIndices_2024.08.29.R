# Script and data by Margit Wilhelm for the following manuscript:
# Otolith shape analysis as a tool for species identification and management 
# of cryptic congeners in the northern Benguela ocean warming hotspot
# Authors:
# M.R. Wilhelm 
# C.E. Jagger 
# N.M. Nghipangelwa
# B.A. Pringle 
# P.W. Shaw
# W.M. Potts
# R. Henriques 
# N.J. McKeown
#
# mwilhelm@unam.na
# 16 April 2024
# Updated 29 August 2024
#
# Originally Based on script by: Smoliński, S., Schade, F. M., & Berg, F. (2020). Assessing the performance of statistical classifiers to discriminate fish stocks using Fourier analysis of otolith shape. Canadian Journal of Fisheries and Aquatic Sciences, 77(4), 674-683.
#
# R version 4.3.0 (2023-04-21 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)

library(psych)
library(vegan)
library(dplyr)
library(tidyr) #Already used earlier

library(caret) #Used for LDA algorithm
library(MASS)  #Used for LDA

library(tidyverse) #already for extraction
library(Momocs) #already for extraction
library(doParallel) #already for extraction
library(foreach) #already for extraction
library(broom) #ALREADY FOR FOURIER

library(scales) #ALREADY FOR FOURIER -> Not yet
library(RColorBrewer) #ALREADY FOR FOURIER -> Not yet
library(ggspatial) #ALREADY FOR FOURIER -> NOt yet

select<-dplyr::select

theme_set(theme_bw())

# clean memory
gc()  
gc(reset=T)
rm(list = ls()) 

############################################################################
############################################################################
path <- dirname(rstudioapi::getActiveDocumentContext()$path)
path
setwd(path) # If not done in previous script and if not run at once

#IF NOT DONE TOGETHER
load("dataBIOAll5.RData")
load("dataSHAPEAll5_norm.RData")
load("eOTOAll5.RData")

############################################################################
############################################################################
## mean shapes
names(eOTO)
#This has the fish ID names

names(dataSHAPE_e)
#This has the Fourier descriptor names.
# [1] "A2"  "A6"  "A8"  "A12" "A13" "A15" "A16" "B6"  "B7"  "B9"  "B10" "B12" "B16" "C3"  "C5"  "C6"  "C7"  "C8" 
# [19] "C9"  "C10" "C11" "C12" "C13" "C14" "C15" "C17" "D1"  "D2"  "D4"  "D5"  "D6"  "D7"  "D10" "D11" "D12" "D13"
# [37] "D14" "D15" "D16" "D17"
#Removed 25, i.e. 40 left: 

names(data)

view(eOTO) #eOTO still has all descriptors included (including A1, B1 and C1)
view(data) #BIO DATA 
view(dataSHAPE_e)

ms_ <- MSHAPES(eOTO, fac = 'Species') 
ms_ <- ms_$shp
datams_<-rbind(data.frame(ms_$A.inodorus, Group="A.inodrous"),
               data.frame(ms_$A.coronus, Group="A.coronus"), 
               data.frame(ms_$Hybrid, Group="Hybrid"))
ms_
view(datams_) #360 observations for mean shapes

ggplot(datams_)+theme_bw()+geom_path( aes(x,y, color=Group, linetype=Group), linewidth=0.3)+theme_void()+theme(legend.position = c(0.6,0.6), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())+coord_equal()+
  theme(legend.key.size = unit(5, units = "mm"))
ggsave("meanshapesBySpp5.tiff", width=10, heigh=5, dpi = 600, compression = "lzw")


############################################################################
############################################################################
#PCA -> For PCA need library Momocs
############################################################################
pca_ <- PCA(eOTO) #Do PCA with all descriptors
names(eOTO)
view(eOTO)

names(data)
data$TLr <-round(data$TL)
data$Latr <-round(data$Lat)


tiff("PCASpp5.tiff", width=8, height = 5,units = "in", res = 600, compression = "lzw")
plot.new()
plot_PCA(pca_, ~Species, zoom =1.0, chull = F, eigen = F, morphospace = F, axesnames = F, axesvar = F, points = F)  %>% 
  layer_morphospace_PCA(., position = c("range"), nb = 12, nr = 6, nc = 5, rotate = 0,
                    size = 0.9, col = "#999999", flipx = FALSE, flipy = FALSE,
                    draw = TRUE)%>%
  layer_ellipses(conf = 0.9, lwd = 1, alpha = 0)%>%
  layer_axes(lwd = 1)%>%
  layer_axesvar(cex=1.5)%>%
  layer_ellipsesaxes (conf = 0.5,lwd=1.5)%>%
  layer_grid(col = "#999999", lty = 3, grid = 3)%>%
  layer_stars(alpha = 0.8)%>%
  layer_points(cex=1.3) %>% 
  layer_eigen(nb_max =5, cex = 1 )%>% 
  layer_legend(probs = seq(0, 1, 0.25), cex = 1)%>%
  layer_title(title = "", cex =1)

dev.off()

tiff("PCA1bSpp5.tiff", width=8, height = 5,units = "in", res = 600, compression = "lzw")
plot(pca_,~Species, pos="xy",xax=1,yax=2,points=T,lwd.shp=1,center=T,grid=T,nb.grids=5,box=T,axisnames=T)
dev.off()

tiff("PCA1cSpp5.tiff", width=8, height = 6,units = "in", res = 600, compression = "lzw")
plot(pca_,~Species,pos="xy",xax=1,yax=2,points=T,lwd.shp=1,center=T,grid=T,nb.grids=5,box=T,axisnames=T, chull = F,chull.filled = F, ellipses = T,ellipsesax = T , zoom = 1.1, stars = T )
dev.off()

tiff("PCA1dSpp5.tiff", width=8, height = 6,units = "in", res = 600, compression = "lzw")
plot(pca_,~Species,xax=1,yax=2,points=T,center=T,grid=T,box=T,axisnames=T, chull = F,chull.filled = F, zoom = 1.1, stars = T )
dev.off()

tiff("PCA1e-ID5.tiff", width=8, height=5,units = "in", res = 600, compression = "lzw")
plot(pca_,~ID,pos="xy",xax=1,yax=2,points=T,lwd.shp=1,center=T,grid=T,nb.grids=5,box=T,axisnames=T)
dev.off()

tiff("PCA1h-sex5.tiff", width=8, height = 6,units = "in", res = 600, compression = "lzw")
plot(pca_,~Sex,pos="xy",xax=1,yax=2,points=T,lwd.shp=1,center=T,grid=T,nb.grids=5,box=T,axisnames=T)
dev.off()

tiff("PCA1j-sppExt5.tiff", width=8, height =5,units = "in", res = 600, compression = "lzw")
plot(pca_,~SpeciesExt,pos="xy",xax=1,yax=2,points=T,lwd.shp=1,center=T,grid=T,nb.grids=5,box=T,axisnames=T)
dev.off()

glimpse(data)
table(data$Species)

############################################################################
############################################################################
load("dataBIOAll5.RData") #for data if not done before

names(data)
dataSHAPE = dataSHAPE_e
dataCLASSIFICATION <- cbind(pop=data$Species, dataSHAPE) 
view(dataCLASSIFICATION)
write.table(dataCLASSIFICATION, "dataShape5WithSpecies.csv")
write.table(dataSHAPE, "dataShape5.csv")

#####################
### CAP --> need library vegan and MASS and dplyr
#####################
library(vegan)
#library(MASS)
#library(dplyr)#uses %>%

names(dataCLASSIFICATION) #pop and 40 descriptors
glimpse(dataCLASSIFICATION) #needs libary(dplyr)
summary(dataCLASSIFICATION)
names(data)

capresults <- capscale(dataCLASSIFICATION%>%dplyr::select(-pop) ~ dataCLASSIFICATION$pop) #pop defined in line 159

capresults
anova(capresults, by = "terms", step = 1000) #Significant differences

eig = eigenvals(capresults, model = "constrained" ) 
eig.ratio = eig/sum(eig) 
eig
eig.ratio

capresults_sites <- scores(capresults)$sites[,1:2] %>% as.data.frame() %>% 
  mutate(ID2=rownames(.)) 

capresults_sites <-cbind(capresults_sites, data)
write.table(capresults_sites, "capresults_sites4.csv")
names(capresults_sites)

capresults_sites%>% 
  ggplot(aes(CAP1, CAP2, color=Species)) +geom_point()+
  stat_ellipse(aes(fill = Species), geom = "polygon", alpha = 0.2, level=0.95)+
  xlab(paste("CAP1:", round(eig.ratio[1]*100,2), "%"))+
  ylab(paste("CAP2:", round(eig.ratio[2]*100,2), "%"))
ggsave("capresultsBySpp_noLengthCorr4.tiff", width=10, heigh=5, dpi = 600, compression = "lzw")

####################
#### Are there any allometric effects?
####################
capresults_sites%>% 
  ggplot(aes(CAP1, CAP2, color=TL))+geom_point()+
  stat_ellipse(aes(fill = TL), geom = "path", alpha = 0.2, level=0.95)+
  xlab(paste("CAP1:", round(eig.ratio[1]*100,2), "%"))+
  ylab(paste("CAP2:", round(eig.ratio[2]*100,2), "%"))

capresults_sites %>% 
  ggplot(aes(TL,CAP1))+geom_point()+geom_smooth(method = "lm")

qplot(data$TL,dataSHAPE$D1)+geom_smooth(method = "lm")

data%>%ggplot(aes(TL))+geom_histogram()+facet_grid(Species~.) #Plots size distribution for the taxa

qplot(data$TL,dataSHAPE$D2, color=data$Spp)+geom_smooth(method = "lm")

###################
## Remove allometric effects (length-correction)
###################
view(dataSHAPE)
dataSHAPE_lc <- dataSHAPE%>%cbind(length=data$TL, pop = data$Species) 
View(dataSHAPE_lc)
write.table(dataSHAPE_lc, "kobshape5withSpp&Length.csv")

#test with ANCOVA Length significance in interaction with pop 
#function gather in tidyr package
library(tidyr)
library(broom)

ANCOVA <- dataSHAPE_lc%>%gather(Descriptor, VAL, -length, -pop)%>%group_by(Descriptor)%>%
  do(tidy(aov(lm(VAL ~ length*pop, data=.))))
view(ANCOVA)
write.csv(ANCOVA, "ancova5.csv")

#filter Descriptors with significant Length effect
sig_length_ <- ANCOVA%>%filter(term=="length"&p.value<0.05)%>%dplyr::select(Descriptor)%>%unlist() 
sig_length_
#20, i.e. also B10 included
# Descriptor1  Descriptor2  Descriptor3  Descriptor4  Descriptor5  Descriptor6  Descriptor7  Descriptor8 
# "A2"         "A6"         "A8"        "B10"         "B6"         "B7"         "B9"        "C15" 
# Descriptor9 Descriptor10 Descriptor11 Descriptor12 Descriptor13 Descriptor14 Descriptor15 Descriptor16 
# "C17"         "C3"         "C5"         "D1"        "D10"        "D12"        "D14"        "D16" 
# Descriptor17 Descriptor18 Descriptor19 Descriptor20 
# "D2"         "D5"         "D6"         "D7" 

view(dataSHAPE_lc) #This dataframe is not yet corrected, but needs to be corrected (B10 included)

#filter Descriptors with significant length and pop interaction and exclude them from the dataset
sig_pop_interaction_ <- ANCOVA%>%filter(Descriptor%in%sig_length_)%>%filter(term=="length:pop")%>%filter(p.value<0.05)%>%dplyr::select(Descriptor)%>%unlist() 
sig_pop_interaction_ #B10
length(sig_pop_interaction_) #1

#Take out those that need no conversion from length 
dataSHAPE_lcB <- dataSHAPE_lc[,-which(colnames(dataSHAPE_lc)%in%as.vector(sig_length_))]
view(dataSHAPE_lcB) #Those that don't need conversion (n=20, B10 excluded)

#Other variables will be corrected - and remove thos with significant length-population interaction- 
#leave this out if none
sig_length_ <- sig_length_[-which(sig_length_%in%sig_pop_interaction_)] 
sig_length_ 
length(sig_length_) 
#19 that have significant length  correction and no significant length-population interaction 
#i.e. Removed B10

LM<-dataSHAPE_lc%>%gather(Descriptor,VAL,-length,-pop)%>%filter(Descriptor%in%sig_length_)%>%group_by(Descriptor)%>%
  do(tidy((lm(VAL ~ length, data=.))))
LM<-LM%>%filter(Descriptor%in%sig_length_)%>%filter(term=="length")%>%dplyr::select(Descriptor,estimate)
LM
write.csv(LM, "LMdestriptors.csv")
#19

#Check relationships of randomly selected variable on plot
dataSHAPE_lc%>%ggplot(aes(length,D2))+geom_point()+geom_smooth(method="lm")

id<-which(colnames(dataSHAPE_lc)%in%sig_length_)
id #identifies the column names that show significant length

view(LM) #shows the coefficients / slopes for correction

dataSHAPE_lc1 <- dataSHAPE_lc[,id] #only the significant length ones -> need correction n=19 (B10 excluded)
dataSHAPE_lc1$id <- rownames(dataSHAPE_lc1) #Add fish ID - the row names
dataSHAPE_lc1$length <- dataSHAPE_lc$length #Add fish length
view(dataSHAPE_lc1) #n=19

dataSHAPE_lc1b <- dataSHAPE_lc1%>%gather(Descriptor, Value,-length,-id) #This makes it into a long sheet with four columns
dataSHAPE_lc1b <- dataSHAPE_lc1b%>%left_join(LM) #This adds the coefficients/slopes from LM for correction

view(dataSHAPE_lcB)#Those that don't need conversion (n=20, B10 excluded)

View(dataSHAPE_lc1b)
write.csv(dataSHAPE_lc1b, "dataShape1c1b.csv") #Correction matrix0

dataSHAPE_lc1b$Value <- dataSHAPE_lc1b$Value - dataSHAPE_lc1b$estimate*dataSHAPE_lc1b$length

dataSHAPE_lc1 <- dataSHAPE_lc1b%>%dplyr::select(Descriptor,Value, id)%>%spread(Descriptor, Value)%>%dplyr::select(-id)

dataSHAPE_lc <- cbind(dataSHAPE_lc1,dataSHAPE_lcB) #Put together corrected and uncorrected. 
#They are now in a different order but all corrected N = 39
view(dataSHAPE_lc)

#Check relationships of randomly selected variable on plot
dataSHAPE_lc%>%ggplot(aes(length,D2))+geom_point()+geom_smooth(method="lm") #length-corrected (flat line)

dataSHAPE_lc<-dataSHAPE_lc%>%dplyr::select(-length,-pop) #remove length and pop from data frame

view(dataSHAPE_lc) ### New matrix of Fourier descriptors, corrected for the length effects  
write.table(dataSHAPE_lc, "kobShape5length-corrected.csv")

######################################
######################################
##############  LDA and CAP on the length-corrected data -> No significant length interaction
######################################
load("dataBIOAll5.RData") #for data if not done before

dataSHAPE_lc = read.table("kobShape5length-corrected.csv")
view(dataSHAPE_lc)

dataCLASSIFICATION_lc <- cbind(pop=data$Species, dataSHAPE_lc)
view(dataCLASSIFICATION_lc)
write.table(dataCLASSIFICATION_lc, "dataCLASSIFICATION_lc.csv")

########################################################################################################################################
#############################
#Linear Discriminant analysis
#####  Discrimination using LDA ON LENGTH-CORRECTED DATA
#Smolinski, S., Schade, F. M., & Berg, F. (2020). Assessing the performance of statistical classifiers to discriminate fish stocks using Fourier analysis of otolith shape. Canadian Journal of Fisheries and Aquatic Sciences, 77(4), 674-683.
#############################
########################################################################################################################################
library(MASS)
dataSHAPE_lc = dataCLASSIFICATION_lc

#OR
dataSHAPE_lc = read.table("dataCLASSIFICATION_lc.csv")

view(dataSHAPE_lc)

tiff("visualisationAll_lengthcorr.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,2:40], 
             gap = 0,
             bg = c("red", "green", "black")[data$Species],
             pch = 21,smooth = F, )
dev.off()

tiff("visualisation1_lengthcorr2.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,2:6], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[data$Species],
             pch = 21,smooth = F, )
dev.off()

tiff("visualisation2_lengthcorr2.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,7:11], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[data$Species],
             pch = 21,smooth = F, )
dev.off()

tiff("visualisation3_lengthcorr2.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,12:16], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[data$Species],
             pch = 21,smooth = F, )
dev.off()

tiff("visualisation4_lengthcorr2.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,17:21], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[data$Species],
             pch = 21,smooth = F, )
dev.off()

tiff("visualisation5_lengthcorr2.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,22:26], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[data$Species],
             pch = 21,smooth = F, )
dev.off()

tiff("visualisation6_lengthcorr2.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,27:31], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[data$Species],
             pch = 21,smooth = F, )
dev.off()

tiff("visualisation7_lengthcorr2.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,32:36], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[data$Species],
             pch = 21,smooth = F, )
dev.off()

tiff("visualisation8_lengthcorr2.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,37:40], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[data$Species],
             pch = 21,smooth = F, )
dev.off()

########################################################################################################################################
#############################
#Linear Discriminant analysis
#############################
########################################################################################################################################
# Uses MASS package
library(MASS)

dataCLASSIFICATION_lc = read.table("dataCLASSIFICATION_lc.csv")

View(dataCLASSIFICATION_lc) #This one is produced in line 310 -> read in again in Line 556
dataCLASSIFICATION_lc$pop = as.factor(dataCLASSIFICATION_lc$pop)

linear <- lda(pop~., dataCLASSIFICATION_lc)
linear
#Shows the discriminatory power in otolith shape to distinguish between two groups = cross-validation
# 
# Call:
#   lda(pop ~ ., data = dataCLASSIFICATION_lc)
# 
# Prior probabilities of groups:
#   A.coronus A.inodorus     Hybrid 
# 0.13364055 0.82949309 0.03686636 
# 
# Group means:
#   A2            A6           A8           B6          B7          B9          C15          C17
# A.coronus  -0.002373862 -1.984144e-03 -0.001619536 -0.005401055 0.001922160 0.003031384 -0.001195380 0.0010613714
# A.inodorus -0.008947004 -9.507678e-05  0.001357534 -0.004868694 0.002078519 0.002373143 -0.001573055 0.0009505872
# Hybrid     -0.009794003  4.231251e-04  0.001070488 -0.003317501 0.002489585 0.001730408 -0.001704698 0.0008322896
# C3          C5         D1          D10         D12          D14         D16          D2
# A.coronus  0.018034243 -0.00861180 -0.6490027 0.0007116069 0.001354896 -0.002566099 0.001556684 -0.02556193
# A.inodorus 0.003889652 -0.00984747 -0.6182911 0.0002834177 0.001031590 -0.002361481 0.001431180 -0.04056551
# Hybrid     0.001108670 -0.01033993 -0.6236060 0.0010339669 0.001693191 -0.002723601 0.001311372 -0.03689989
# D5            D6          D7          A12          A13          A15           A16
# A.coronus  -0.01020897 -0.0003284411 0.003137023 0.0012594103 0.0011579355 0.0006438183  3.330405e-04
# A.inodorus -0.01016667  0.0057887461 0.006131121 0.0011201502 0.0007879199 0.0005358044 -7.790539e-05
# Hybrid     -0.01134217  0.0049452883 0.006433767 0.0009574997 0.0004133590 0.0007012752 -3.061124e-04
# B12           B16           C6           C7            C8           C9         C10
# A.coronus   0.0007714446 -1.719215e-05 -0.002850971 -0.006229655 -0.0008130481 0.0003551613 0.001735512
# A.inodorus -0.0008061839 -1.019612e-04 -0.004660560 -0.002160250  0.0006505090 0.0019494702 0.002139885
# Hybrid     -0.0015851687  1.924544e-05 -0.002268888 -0.003078627  0.0003843475 0.0025732062 0.003412662
# C11           C12           C13          C14           D4           D11           D13
# A.coronus  0.0014717377 -0.0004452489 -4.810031e-04 5.586864e-04 -0.001143296  0.0001186575 -5.482047e-04
# A.inodorus 0.0005574137 -0.0022416617  1.478451e-05 4.731837e-04  0.006331078 -0.0022519940 -5.629808e-05
# Hybrid     0.0016992730 -0.0026636242 -7.946901e-04 7.664099e-05  0.004235784 -0.0015931926 -1.276460e-04
# D15           D17
# A.coronus  -0.0007190158  3.162957e-04
# A.inodorus -0.0004120076  4.963305e-06
# Hybrid     -0.0001880225 -2.624064e-04
# 
# Coefficients of linear discriminants:
#   LD1        LD2
# A2    10.432404  -73.39666
# A6   -35.876234   88.03908
# A8   104.300414  531.18803
# B6   105.917518  251.04673
# B7  -238.798143  -59.24617
# B9   189.023004 -351.37789
# C15   13.649258  131.48726
# C17 -235.970136  212.82971
# C3   -65.867307  -16.46221
# C5   -34.002200  -94.63983
# D1    16.640950  -13.31141
# D10  -89.406661  344.59546
# D12  -94.347764  276.31296
# D14 -182.545221  -76.64182
# D16 -173.725090  -88.29249
# D2   -30.002506   19.74362
# D5     5.198687 -116.75314
# D6   143.244056 -155.01772
# D7   113.309632  -66.77389
# A12  -83.628076 -330.15461
# A13  -54.655993   16.79081
# A15   73.376251 -379.90556
# A16  -87.712508 -137.52604
# B12 -161.724792 -632.79978
# B16   49.994770   62.71384
# C6    65.907412   41.15658
# C7   226.521190 -113.96345
# C8    59.946289  -78.25881
# C9   191.197768  -40.83671
# C10  184.017630 -116.21253
# C11  260.467699  464.30025
# C12   67.782318   20.71257
# C13  334.256997   31.91792
# C14  -63.164505   56.99084
# D4    67.873559  -89.30710
# D11  -74.541769   81.81071
# D13 -167.341118  -14.16329
# D15 -358.985283 -167.15238
# D17  -43.932205 -372.89193
# 
# Proportion of trace:
#   LD1   LD2 
# 0.918 0.082 

plot(linear)

#To improve plot:
library(latticeExtra)

tiff("LDA1Spp5-lengthCorrN.tiff", width=5, heigh=5,units = "in", res = 600, compression = "lzw")
opar <- par(mar = c(4,5,1,1), oma = c(2,0,0,0))
plot(linear, xlab = "LD1: 91.8%", ylab = "LD2: 8.2%", col = dataCLASSIFICATION_lc$pop, 
     panel = points, pch = 4, lwd = 3, cex = 3, cex.lab = 1.5) 
dev.off()

#This is for the key in a separate plot. 
fkey = list(x = .95, y = 0.01, corner = c(1, 0),
            text = list(c("A. inodorus", 
                          "A. coronus",
                          "Hybrids")), 
            points = list(type = c("p", "p", "p"), 
                         col = c("2", "1", "3"),
                         pch = c(4, 4, 4), 
                         lwd = c(3, 3, 3)))
xyplot(A2 ~ A6, data=dataCLASSIFICATION_lc, key = fkey)

########################################################################################################################################
###### 4-fold Cross-Validation
########################################################################################################################################
cl <- makeCluster(detectCores())
registerDoParallel(cl)
cl

### info about libraries
#getModelInfo("lda")$lda$library

trcntr<-trainControl(method = "repeatedcv", number=4, repeats=100, allowParallel = T)
prepr=c('scale', 'center')

######################################
######################################
dataCLASSIFICATION_lc$pop <- as.factor(dataCLASSIFICATION_lc$pop)

LDA_results <- train(pop~., method='lda', preProcess=prepr, data=dataCLASSIFICATION_lc, trControl = trcntr)
#method = 'lda' (most commonly used) = fast

OriginalLDA_results <- LDA_results

#Original used LDA, but best according to Smolinski et al. 2019 is support Vector Machines (SVM)
#Least Squares Support Vector Machine with Radial Basis Function Kernel
#method = 'lssvmRadial'
#package required: kernlab 

LDA_resultsSVM <- train(pop~., method='lssvmRadial', preProcess=prepr, data=dataCLASSIFICATION_lc, trControl = trcntr)
#Slow method

######################
## Accuracy, Kappa
######################
#FOR ONLY FOURIER DESCRIPTORS
OriginalLDA_results

confusionMatrix(LDA_results)

LDA_resultsSVM

confusionMatrix(LDA_resultsSVM)

######################
## Variable importance
######################
varImp(LDA_results)
plot(varImp(LDA_results))

tiff("ImportanceSpp4_lengthcorr.tiff", width=6, heigh=12,units = "in", res = 600, compression = "lzw")
plot(varImp(LDA_results))
dev.off()

#####################
### CAP 
### Need library: vegan
#####################
library(vegan)

view(dataCLASSIFICATION_lc)

capresults <- capscale(dataCLASSIFICATION_lc%>%dplyr::select(-pop) ~ dataCLASSIFICATION_lc$pop) 
capresults
anova(capresults, by = "terms", step = 1000) #Sign. difference

eig = eigenvals(capresults, model = "constrained") 
eig.ratio = eig/sum(eig) 
eig.ratio

capresults_sites <- scores(capresults)$sites[,1:2] %>% as.data.frame() %>% 
  mutate(ID2=rownames(.)) 
view(capresults_sites) #Here only CAP1 and CAP2 for each otolith ID

####Now add data
capresults_sites <-cbind(capresults_sites, data)
write.csv(capresults_sites, "capresults_sitesCorrectedForLength4.csv")

names(capresults_sites)

tiff("Capresults_Spp5_LC.tiff", width=12, heigh=6,units = "in", res = 600, compression = "lzw")
capresults_sites%>% 
  ggplot(aes(CAP1, CAP2, color=Species)) +geom_point()+
  stat_ellipse(aes(fill = Species), geom = "polygon", alpha = 0.2, level=0.95)+
  xlab(paste("CAP1:", round(eig.ratio[1]*100,2), "%"))+
  ylab(paste("CAP2:", round(eig.ratio[2]*100,2), "%"))
dev.off()

tiff("Capresults_Lat5_LC.tiff", width=12, height=6,units = "in", res = 600, compression = "lzw")
capresults_sites%>% 
  ggplot(aes(CAP1, CAP2, color=Lat))+geom_point()+
  stat_ellipse(aes(fill = Lat), geom = "polygon", alpha = 0.2, level=0.95)+
  xlab(paste("CAP1:", round(eig.ratio[1]*100,2), "%"))+
  ylab(paste("CAP2:", round(eig.ratio[2]*100,2), "%"))
dev.off()

tiff("Capresults_Size4_LC.tiff", width=12, height=6,units = "in", res = 600, compression = "lzw")
capresults_sites%>% 
  ggplot(aes(CAP1, CAP2, color=TL))+geom_point()+
  stat_ellipse(aes(fill = TL), geom = "polygon", alpha = 0.2, level=0.95)+
  xlab(paste("CAP1:", round(eig.ratio[1]*100,2), "%"))+
  ylab(paste("CAP2:", round(eig.ratio[2]*100,2), "%"))
dev.off()

tiff("Capresults_mtDNA4_LC.tiff", width=12, height=6,units = "in", res = 600, compression = "lzw")
capresults_sites%>% 
  ggplot(aes(CAP1, CAP2, color=mtDNA))+geom_point()+
  stat_ellipse(aes(fill = mtDNA), geom = "polygon", alpha = 0.2, level=0.95)+
  xlab(paste("CAP1:", round(eig.ratio[1]*100,2), "%"))+
  ylab(paste("CAP2:", round(eig.ratio[2]*100,2), "%"))
dev.off()

tiff("Capresults_morph4_LC.tiff", width=12, height=6,units = "in", res = 600, compression = "lzw")
capresults_sites%>% 
  ggplot(aes(CAP1, CAP2, color=SpeciesExt))+geom_point()+
  stat_ellipse(aes(fill = SpeciesExt), geom = "polygon", alpha = 0.2, level=0.95)+
  xlab(paste("CAP1:", round(eig.ratio[1]*100,2), "%"))+
  ylab(paste("CAP2:", round(eig.ratio[2]*100,2), "%"))
dev.off()

########################################################################################################################################
##########################
# NOW ADD SHAPE INDICES TO FOURIER DESCRIPTORS AND RUN LDA AGAIN
# This can be done after extraction for any data using read.table to add data
##########################
########################################################################################################################################
#Shape index data are cleaned and obtained from SCRIPT: "anovaEtcOfShapeIndices_kob.R"
library(dplyr)

path <- dirname(rstudioapi::getActiveDocumentContext()$path)
path
setwd(path) # If not done in previous script and if not run at once

sapeInd = read.table("shapeIndices_forLDA.csv", header = TRUE) 
view(sapeInd)


dataCLASS = dataCLASSIFICATION_lc
#OR IF not done: 
dataCLASS = read.table("dataCLASSIFICATION_lc.csv", header = TRUE) # Do this if not done at once. 
view(dataCLASS) #Had ID as row names, and pop as column

dataCLASS$ID = sapeInd$ID #need ID as a column
view(dataCLASS)

sapeInd = sapeInd%>%dplyr::select(-pop) #Remove those that occur in both dataCLASS and sapeInd

dataCLASSIFICATION_lcFull = dataCLASS%>%left_join(sapeInd) #CROSS-JOIN REPEATS IT FOR POP!! Need ID in each!! If continuing from here - or "left_join" if ID in both

view(dataCLASSIFICATION_lcFull) #217 rows (fish); 46 columns

#Analysis including shape indices
dataCLASSIFICATION_lc = dataCLASSIFICATION_lcFull%>%dplyr::select(-length, -ID) #Also remove ID if necessary

#Analysis excluding shape indices
dataCLASSIFICATION_lc2 = dataCLASSIFICATION_lc%>%dplyr::select(-circ, -OCD.OH, -round, -ellip)


################################################################################
#############################
##### Discrimination using LDA ON LENGTH-CORRECTED DATA INCLUDING AND EXCLUDING SHAPE INDICES
#############################
################################################################################
library(MASS)

names(dataCLASSIFICATION_lc)
names(dataCLASSIFICATION_lc2)

dataCLASSIFICATION_lc$pop <- as.factor(dataCLASSIFICATION_lc$pop)
dataCLASSIFICATION_lc2$pop <- as.factor(dataCLASSIFICATION_lc2$pop)

linear <-  lda(pop~., dataCLASSIFICATION_lc)  #including shape indices
linear2 <- lda(pop~., dataCLASSIFICATION_lc2) #excluding shape indices

linear

# Call:
#   lda(pop ~ ., data = dataCLASSIFICATION_lc)
# 
# Prior probabilities of groups:
#   A.coronus A.inodorus     Hybrid 
# 0.13364055 0.82949309 0.03686636 
# 
# Group means:
#   A2            A6           A8           B6          B7          B9          C15          C17
# A.coronus  -0.002373862 -1.984144e-03 -0.001619536 -0.005401055 0.001922160 0.003031384 -0.001195380 0.0010613714
# A.inodorus -0.008947004 -9.507678e-05  0.001357534 -0.004868694 0.002078519 0.002373143 -0.001573055 0.0009505872
# Hybrid     -0.009794003  4.231251e-04  0.001070488 -0.003317501 0.002489585 0.001730408 -0.001704698 0.0008322896
# C3          C5         D1          D10         D12          D14         D16          D2
# A.coronus  0.018034243 -0.00861180 -0.6490027 0.0007116069 0.001354896 -0.002566099 0.001556684 -0.02556193
# A.inodorus 0.003889652 -0.00984747 -0.6182911 0.0002834177 0.001031590 -0.002361481 0.001431180 -0.04056551
# Hybrid     0.001108670 -0.01033993 -0.6236060 0.0010339669 0.001693191 -0.002723601 0.001311372 -0.03689989
# D5            D6          D7          A12          A13          A15           A16
# A.coronus  -0.01020897 -0.0003284411 0.003137023 0.0012594103 0.0011579355 0.0006438183  3.330405e-04
# A.inodorus -0.01016667  0.0057887461 0.006131121 0.0011201502 0.0007879199 0.0005358044 -7.790539e-05
# Hybrid     -0.01134217  0.0049452883 0.006433767 0.0009574997 0.0004133590 0.0007012752 -3.061124e-04
# B12           B16           C6           C7            C8           C9         C10
# A.coronus   0.0007714446 -1.719215e-05 -0.002850971 -0.006229655 -0.0008130481 0.0003551613 0.001735512
# A.inodorus -0.0008061839 -1.019612e-04 -0.004660560 -0.002160250  0.0006505090 0.0019494702 0.002139885
# Hybrid     -0.0015851687  1.924544e-05 -0.002268888 -0.003078627  0.0003843475 0.0025732062 0.003412662
# C11           C12           C13          C14           D4           D11           D13
# A.coronus  0.0014717377 -0.0004452489 -4.810031e-04 5.586864e-04 -0.001143296  0.0001186575 -5.482047e-04
# A.inodorus 0.0005574137 -0.0022416617  1.478451e-05 4.731837e-04  0.006331078 -0.0022519940 -5.629808e-05
# Hybrid     0.0016992730 -0.0026636242 -7.946901e-04 7.664099e-05  0.004235784 -0.0015931926 -1.276460e-04
# D15           D17     circ    OCD.OH     round     ellip
# A.coronus  -0.0007190158  3.162957e-04 13.81026 0.3565555 0.5803786 0.2822302
# A.inodorus -0.0004120076  4.963305e-06 14.71079 0.4849641 0.5342851 0.2860512
# Hybrid     -0.0001880225 -2.624064e-04 14.74290 0.4797834 0.5385447 0.2825055
# 
# Coefficients of linear discriminants:
#   LD1          LD2
# A2       19.9698945  -83.2239448
# A6      -27.1066957   59.2809859
# A8       -0.2791453  576.7315234
# B6       67.2465107  292.6403645
# B7     -198.8020833 -108.0722947
# B9      154.7462424 -387.1085435
# C15      58.2528175  192.9463045
# C17    -103.5026196  238.0998404
# C3      -53.3016037  -24.7372763
# C5       -6.2541000  -93.9409736
# D1       -2.6215548  -37.2687264
# D10     -99.3394926  355.3445123
# D12    -141.7820077  272.6429092
# D14    -113.7813382  -25.7199409
# D16     -36.5003856  -32.7881302
# D2      -17.5802570   22.5614678
# D5       11.7497169 -110.4872195
# D6       97.6872265 -179.7565413
# D7       77.1491800  -55.6748556
# A12      -0.6697420 -221.8614933
# A13     -53.7925343  -77.1601062
# A15     108.9420792 -412.0679467
# A16    -192.5532259 -111.6954213
# B12    -101.7373063 -715.9724027
# B16     -16.2937067   75.7532898
# C6       69.0686669   39.1726058
# C7      216.6759603 -113.8673929
# C8       29.5254401 -132.2448402
# C9      189.0252146   -8.7294662
# C10     171.4921581 -186.3828222
# C11     252.7018371  505.9459504
# C12      70.0001323   -0.8451842
# C13     313.2070264   85.3630524
# C14    -174.1063538   49.3793647
# D4       54.9577727 -101.1988338
# D11     -96.1894586  135.9405137
# D13    -177.7107230  -52.3148945
# D15    -267.5824401 -168.7364047
# D17      -0.8205623 -375.1650629
# circ      0.5518115   -0.9771813
# OCD.OH    7.9085201   -2.4075319
# round   -19.5922319  -42.0001578
# ellip   -22.0049338   -4.5767815
# 
# Proportion of trace:
#   LD1    LD2 
# 0.9234 0.0766 

linear2 #Same as before (exclusing shape indices)

plot.new()
op=par(mfrow = c(1,1), mar=c(4.2, 3.8, 0.2, 0.2))
plot(linear)
dev.off()

#Improve plot to points

tiff("LDA1Spp5-lengthCorr_SI4.tiff", width = 8, height = 6,units = "in", res = 600, compression = "lzw")
opar <- par(mar = c(4,5,1,1), oma = c(2,0,0,0))
plot(linear, xlab = "LD1: 92.3%", ylab = "LD2: 7.7%", col = dataCLASSIFICATION_lc$pop, 
     panel = points, pch = 4, lwd = 3, cex = 3, cex.lab = 1.5) 
#Done
#pop needs to be done 'as.factor'; otherwise Error "Invalid color name: 'A. inodorus'"
dev.off()

tiff("LDA1Spp5-lengthCorr_SI3ab2.tiff", width = 6, height = 8,units = "in", res = 600, compression = "lzw")
op=par(mfrow = c(2,1), mar=c(4.2, 3.8, 0.2, 0.2))
plot(linear)  #Compared with and without Shape indices
plot(linear2) #WITH, HYBRIDS stand out slightly better
dev.off()

########################################################################################################################################
###### Leave-one-out Cross-Validation
########################################################################################################################################
library(doParallel)

cl <- makeCluster(detectCores())
registerDoParallel(cl)

library(caret)

trcntr <- trainControl(method = "repeatedcv", number=4, repeats=100, allowParallel = T)
trcntr1 <- trainControl(method = "LOOCV", allowParallel = T)

prepr=c('scale', 'center')

#4-fold
LDA_results <- train(pop~., method='lda', preProcess=prepr, data=dataCLASSIFICATION_lc, trControl = trcntr)
LDA_results
confusionMatrix(LDA_results) #Can't compute confusion matrix - only use prediction

#LOOCV
LDA_results1 <- train(pop~., method='lda', preProcess=prepr, data=dataCLASSIFICATION_lc, trControl = trcntr1)
LDA_results1 #including SI - LOOCV
confusionMatrix(LDA_results1) #Can't compute confusion matrix - only use prediction

#4-fold
LDA_resultsSVM1 <- train(pop~., method='lssvmRadial', preProcess=prepr, data=dataCLASSIFICATION_lc, trControl = trcntr)
LDA_resultsSVM1
confusionMatrix(LDA_resultsSVM1)

dataCLASSIFICATION_lc_p = dataCLASSIFICATION_lc
dataCLASSIFICATION_lc_p$predLOOCV = predict(LDA_results1, newdata=dataCLASSIFICATION_lc)

view(dataCLASSIFICATION_lc_p)
write.csv(dataCLASSIFICATION_lc_p, "fourierDescrip_withpred3.csv")

getwd()

save(LDA_results1, file = "LDA_resultsLOOCV.RData")

########################################################################################################################################
###### K-fold Cross-Validation
########################################################################################################################################
trcntr5 <- trainControl(method = "repeatedcv", number=5, repeats=100, allowParallel = T) #library:caret 
trcntr10 <- trainControl(method = "repeatedcv", number=10, repeats=100, allowParallel = T) #library:caret 
trcntr100 <- trainControl(method = "repeatedcv", number=217, repeats = 100, allowParallel = T) #library:caret 
prepr=c('scale', 'center')

LDA_results5 <- train(pop~., method='lda', preProcess=prepr, data=dataCLASSIFICATION_lc, trControl = trcntr5)
LDA_results10 <- train(pop~., method='lda', preProcess=prepr, data=dataCLASSIFICATION_lc, trControl = trcntr10)
LDA_results100 <- train(pop~., method='lda', preProcess=prepr, data=dataCLASSIFICATION_lc, trControl = trcntr100)

confusionMatrix(LDA_results5)
confusionMatrix(LDA_results10)
confusionMatrix(LDA_results100)

########################################################################################################################################
###### STRATIFIED K-fold Cross-Validation
########################################################################################################################################
trcntr_stratified10 <- trainControl(method = "cv", number=10, classProbs = TRUE, allowParallel = T) #library:caret 
prepr=c('scale', 'center')

LDA_results_str10 <- train(pop~., method='lda', preProcess=prepr, data=dataCLASSIFICATION_lc2, trControl = trcntr_stratified10)

LDA_results_str10 
confusionMatrix(LDA_results_str10)

########################################################################################################################################
###### 4-fold Cross-Validation
########################################################################################################################################
library(doParallel)

cl <- makeCluster(detectCores())
registerDoParallel(cl)
cl

### info about libraries
getModelInfo("lda")$lda$library
library(caret)

trcntr <- trainControl(method = "repeatedcv", number=4, repeats=100, allowParallel = T) #library:caret 
prepr=c('scale', 'center')

######################################
######################################
LDA_results <- train(pop~., method='lda', preProcess=prepr, data=dataCLASSIFICATION_lc, trControl = trcntr)
LDA_results #including SI
confusionMatrix(LDA_results) #including SI

LDA_results2 <- train(pop~., method='lda', preProcess=prepr, data=dataCLASSIFICATION_lc2, trControl = trcntr)

LDA_results2 #excluding SI
confusionMatrix(LDA_results2)


#############################
##### Discrimination using SVM ON LENGTH-CORRECTED DATA
#############################
library(kernlab)

LDA_resultsSVM <- train(pop~., method='lssvmRadial', preProcess=prepr, data=dataCLASSIFICATION_lc, trControl = trcntr)

LDA_resultsSVM

confusionMatrix(LDA_resultsSVM)
# Cross-Validated (4 fold, repeated 100 times) Confusion Matrix 
# 
# (entries are percentual average cell counts across resamples)
# 
#               Reference
# Prediction   A.coronus A.inodorus Hybrid
# A.coronus       11.7        1.4    0.1
# A.inodorus       1.6       81.6    3.6
# Hybrid           0.0        0.0    0.0
# 
# Accuracy (average) : 0.9331

dataCLASSIFICATION_1p$predIncl = predict(LDA_results, newdata=dataCLASSIFICATION_lc)
dataCLASSIFICATION_1p$predExcl = predict(LDA_results2, newdata=dataCLASSIFICATION_lc2)
dataCLASSIFICATION_1p$predSVM = predict(LDA_resultsSVM, newdata=dataCLASSIFICATION_lc) 

view(dataCLASSIFICATION_1p)
write.csv(dataCLASSIFICATION_1p, "fourierDescrip_withpred3.csv")

save(LDA_resultsSVM, file = "LDA_resultsSVM.RData")

######################
## Variable importance
######################
varImp(LDA_results1)
# ROC curve variable importance
# 
# variables are sorted by maximum importance across the classes
# only 20 most important variables shown (out of 44)
# 
#         A.coronus A.inodorus Hybrid
# OCD.OH    100.00     100.00  89.81
# circ       98.32      95.25  98.32
# C3         97.37      85.92  97.37
# D2         90.76      89.03  90.76
# B12        86.03      60.92  86.03
# round      84.33      84.33  80.36
# C12        81.30      60.55  81.30
# D6         80.84      80.84  58.61
# A8         79.24      79.24  63.34
# D11        76.97      76.97  68.07
# C7         71.85      71.22  71.85
# D1         69.16      69.16  64.29
# D7         67.12      53.91  67.12
# D4         65.00      65.00  64.29
# A2         55.78      46.68  55.78
# B6         55.78      39.26  55.78
# B10        54.96      54.96  52.94
# C10        53.89      33.93  53.89
# A6         53.89      44.62  53.89
# B9         48.21      26.43  48.21

plot(varImp(LDA_results1))

tiff("ImportanceSpp4_lengthcorr_4.tiff", width=6, height = 12,units = "in", res = 600, compression = "lzw")
plot(varImp(LDA_results1))
dev.off()
