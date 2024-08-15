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
# Updated 05 August 2024
#
# Based on script by: Smoliński, S., Schade, F. M., & Berg, F. (2020). Assessing the performance of statistical classifiers to discriminate fish stocks using Fourier analysis of otolith shape. Canadian Journal of Fisheries and Aquatic Sciences, 77(4), 674-683.
#
# R version 4.3.0 (2023-04-21 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)

library(psych)
library(vegan)
library(dplyr)
library(tidyr) #already used earlier

library(caret) #used for lda alogorithm
library(MASS)  #USed for LDA

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

#iF NOT DONE TOGETHER
load("dataBIOAll5_norm.RData")
load("dataSHAPEAll5.RData")
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
#20
# Descriptor1  Descriptor2  Descriptor3  Descriptor4  Descriptor5  Descriptor6  Descriptor7  Descriptor8 
# "A2"         "A6"         "A8"        "B10"         "B6"         "B7"         "B9"        "C15" 
# Descriptor9 Descriptor10 Descriptor11 Descriptor12 Descriptor13 Descriptor14 Descriptor15 Descriptor16 
# "C17"         "C3"         "C5"         "D1"        "D10"        "D12"        "D14"        "D16" 
# Descriptor17 Descriptor18 Descriptor19 Descriptor20 
# "D2"         "D5"         "D6"         "D7" 


#filter Descriptors with significant length and pop interaction and exclude them from the dataset
sig_pop_interaction_ <- ANCOVA%>%filter(Descriptor%in%sig_length_)%>%filter(term=="length:pop")%>%filter(p.value<0.05)%>%dplyr::select(Descriptor)%>%unlist() 
sig_pop_interaction_ #B10
length(sig_pop_interaction_) #1

#Take out those that need no conversion from length 
dataSHAPE_lcB <- dataSHAPE_lc[,-which(colnames(dataSHAPE_lc)%in%as.vector(sig_length_))]
view(dataSHAPE_lcB) #Those that don't need conversion

#Other variables will be corrected - and remove thos with significant length-population interaction- 
#leave this out if none
sig_length_ <- sig_length_[-which(sig_length_%in%sig_pop_interaction_)] 
sig_length_ 
length(sig_length_) #19 that have significant length  correction and no significant length-population interaction 
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
view(dataSHAPE_lc) #This dataframe is not yet corrected, but need to be corrected
view(LM) #shows the coefficients / slopes for correction

dataSHAPE_lc1 <- dataSHAPE_lc[,id] #only the significant length ones -> need correction n=19
dataSHAPE_lc1$id <- rownames(dataSHAPE_lc1) #Add fish ID - the row names
dataSHAPE_lc1$length <- dataSHAPE_lc$length #Add fish length
view(dataSHAPE_lc1)

dataSHAPE_lc2 <- dataSHAPE_lc[,-id] #"not ID", no correction needed n = 20

dataSHAPE_lc1b <- dataSHAPE_lc1%>%gather(Descriptor, Value,-length,-id) #This makes it into a long sheet with four columns
dataSHAPE_lc1b <- dataSHAPE_lc1b%>%left_join(LM) #This adds the coefficients/slopes from LM for correction

View(dataSHAPE_lc1b)
write.csv(dataSHAPE_lc1b, "dataShape1c1b.csv") #Correction matrix0

dataSHAPE_lc1b$Value <- dataSHAPE_lc1b$Value - dataSHAPE_lc1b$estimate*dataSHAPE_lc1b$length

dataSHAPE_lc1 <- dataSHAPE_lc1b%>%dplyr::select(Descriptor,Value, id)%>%spread(Descriptor, Value)%>%dplyr::select(-id)

dataSHAPE_lc <- cbind(dataSHAPE_lc1,dataSHAPE_lc2) #Put together corrected and uncorrected. They are now in a different order but all corrected
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
view(dataSHAPE_lc)

tiff("visualisationAll_lengthcorr.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,1:39], 
             gap = 0,
             bg = c("red", "green", "black")[data$Species],
             pch = 21,smooth = F, )
dev.off()

tiff("visualisation1_lengthcorr2.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,1:5], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[data$Species],
             pch = 21,smooth = F, )
dev.off()

tiff("visualisation2_lengthcorr2.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,6:10], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[data$Species],
             pch = 21,smooth = F, )
dev.off()

tiff("visualisation3_lengthcorr2.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,11:15], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[data$Species],
             pch = 21,smooth = F, )
dev.off()

tiff("visualisation4_lengthcorr2.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,16:20], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[data$Species],
             pch = 21,smooth = F, )
dev.off()

tiff("visualisation5_lengthcorr2.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,21:25], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[data$Species],
             pch = 21,smooth = F, )
dev.off()

tiff("visualisation6_lengthcorr2.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,26:30], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[data$Species],
             pch = 21,smooth = F, )
dev.off()

tiff("visualisation7_lengthcorr2.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,31:35], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[data$Species],
             pch = 21,smooth = F, )
dev.off()

tiff("visualisation8_lengthcorr2.tiff", width=6, height=6 ,units = "in", res = 600, compression = "lzw")
pairs.panels(dataSHAPE_lc[,36:39], #package(psych) only use first four Fourier coefficients - visualise
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

View(dataCLASSIFICATION_lc)
linear <- lda(pop~., dataCLASSIFICATION_lc)
linear
#Shows the discriminatory power in otolith shape to distinguish between two groups = cross-validation

plot(linear)

#To improve plot:
library(latticeExtra)

tiff("LDA1Spp5-lengthCorrN.tiff", width=5, heigh=5,units = "in", res = 600, compression = "lzw")
opar <- par(mar = c(4,5,1,1), oma = c(2,0,0,0))
plot(linear, xlab = "LD1: 91.76%", ylab = "LD2: 8.24%", col = dataCLASSIFICATION_lc$pop, 
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
dataCLASSIFICATION_lc$pop<-as.factor(dataCLASSIFICATION_lc$pop)

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
#NOW ADD SHAPE INDICES TO FOURIER DESCRIPTORS AND RUN LDA AGAIN
# This can be done after extraction for any data using read.table to add data
##########################
########################################################################################################################################
#Shape index data are cleaned and obtained from SCRIPT: "anovaEtcOfShapeIndices_kob.R"
library(dplyr)

path <- dirname(rstudioapi::getActiveDocumentContext()$path)
path
setwd(path) # If not done in previous script and if not run at once

sapeInd = read.table("shapeIndices_forLDA.csv", header = TRUE) 
names(sapeInd)

dataCLASS = read.table("dataCLASSIFICATION_lc.csv", header = TRUE) # Do this if not done at once. 
names(dataCLASS) #Had ID as row names, and pop as column

sapeInd = sapeInd%>%dplyr::select(-pop) #Remove those that occur in both dataCLASS and sapeInd

dataCLASSIFICATION_lcFull = dataCLASS%>%cross_join(sapeInd) #If continuing from here - or "left_join" if ID in both

names(dataCLASSIFICATION_lcFull)

#Analysis including shape indices
dataCLASSIFICATION_lc = dataCLASSIFICATION_lcFull%>%dplyr::select(-length, -ID) #Also remove ID if necessary

#Analysis excluding shape indices
dataCLASSIFICATION_lc2 = dataCLASSIFICATION_lc%>%dplyr::select(-circ, -OCD.OH, -round, -ellip)

#Analysis excluding some variables (for SVM)
dataCLASSIFICATION_lc3 = dataCLASSIFICATION_lc2%>%dplyr::select(-C3, -C6, -C14, -D1, -D2, -A2, -B16, -A13)

################################################################################
#############################
##### Discrimination using LDA ON LENGTH-CORRECTED DATA INCLUDING AND EXCLUDING SHAPE INDICES
#############################
################################################################################
library(MASS)

names(dataCLASSIFICATION_lc)
names(dataCLASSIFICATION_lc2)
names(dataCLASSIFICATION_lc3)

dataCLASSIFICATION_lc$pop <- as.factor(dataCLASSIFICATION_lc$pop)
dataCLASSIFICATION_lc2$pop <- as.factor(dataCLASSIFICATION_lc2$pop)
dataCLASSIFICATION_lc3$pop <- as.factor(dataCLASSIFICATION_lc3$pop)

linear <- lda(pop~., dataCLASSIFICATION_lc) #including shape indices
linear2 <- lda(pop~., dataCLASSIFICATION_lc2) #excluding shape indices
linear3 <- lda(pop~., dataCLASSIFICATION_lc3) #excluding shape indices and some descriptors

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
# B10           B12           B16           C6           C7            C8           C9
# A.coronus   0.0011043734  0.0007714446 -1.719215e-05 -0.002850971 -0.006229655 -0.0008130481 0.0003551613
# A.inodorus -0.0003840593 -0.0008061839 -1.019612e-04 -0.004660560 -0.002160250  0.0006505090 0.0019494702
# Hybrid     -0.0005269451 -0.0015851687  1.924544e-05 -0.002268888 -0.003078627  0.0003843475 0.0025732062
# C10          C11           C12           C13          C14           D4           D11
# A.coronus  0.001735512 0.0014717377 -0.0004452489 -4.810031e-04 5.586864e-04 -0.001143296  0.0001186575
# A.inodorus 0.002139885 0.0005574137 -0.0022416617  1.478451e-05 4.731837e-04  0.006331078 -0.0022519940
# Hybrid     0.003412662 0.0016992730 -0.0026636242 -7.946901e-04 7.664099e-05  0.004235784 -0.0015931926
# D13           D15           D17     circ    OCD.OH     round     ellip
# A.coronus  -5.482047e-04 -0.0007190158  3.162957e-04 14.59162 0.4676125 0.5406021 0.2854098
# A.inodorus -5.629808e-05 -0.0004120076  4.963305e-06 14.59162 0.4676125 0.5406021 0.2854098
# Hybrid     -1.276460e-04 -0.0001880225 -2.624064e-04 14.59162 0.4676125 0.5406021 0.2854098
# 
# Coefficients of linear discriminants:
#   LD1           LD2
# A2      9.832593e+00 -7.192944e+01
# A6     -3.860824e+01  9.478158e+01
# A8      1.157654e+02  5.055456e+02
# B6      1.110955e+02  2.406384e+02
# B7     -2.397718e+02 -6.153758e+01
# B9      1.899036e+02 -3.514149e+02
# C15     1.079034e+01  1.396803e+02
# C17    -2.414114e+02  2.232481e+02
# C3     -6.590159e+01 -1.770230e+01
# C5     -3.357030e+01 -9.678460e+01
# D1      1.692160e+01 -1.377214e+01
# D10    -9.067934e+01  3.475345e+02
# D12    -9.540801e+01  2.783378e+02
# D14    -1.853268e+02 -7.323333e+01
# D16    -1.832162e+02 -6.736049e+01
# D2     -3.061207e+01  2.082616e+01
# D5      4.845990e+00 -1.161977e+02
# D6      1.413605e+02 -1.479929e+02
# D7      1.092528e+02 -5.434073e+01
# A12    -8.814481e+01 -3.213443e+02
# A13    -6.164569e+01  3.393436e+01
# A15     8.477460e+01 -4.095504e+02
# A16    -1.016434e+02 -1.036227e+02
# B10    -3.433754e+01  8.904876e+01
# B12    -1.466368e+02 -6.775199e+02
# B16     4.998138e+01  6.395473e+01
# C6      6.252155e+01  5.136411e+01
# C7      2.272362e+02 -1.119177e+02
# C8      5.626438e+01 -6.786742e+01
# C9      1.957558e+02 -4.914796e+01
# C10     1.901273e+02 -1.289818e+02
# C11     2.622264e+02  4.665642e+02
# C12     8.308582e+01 -1.759207e+01
# C13     3.378958e+02  2.902248e+01
# C14    -5.151068e+01  2.578023e+01
# D4      6.633411e+01 -8.436310e+01
# D11    -7.596400e+01  8.438986e+01
# D13    -1.734667e+02 -1.545171e+00
# D15    -3.657374e+02 -1.571894e+02
# D17    -4.152790e+01 -3.814363e+02
# circ    1.369232e-15  2.387063e-15
# OCD.OH  3.152354e-14  3.313689e-14
# round   1.118667e-13  1.048435e-13
# ellip   5.545731e-14  1.776415e-14
# 
# Proportion of trace:
#   LD1    LD2 
# 0.9176 0.0824 

linear2
# Call:
#   lda(pop ~ ., data = dataCLASSIFICATION_lc2)
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
# B10           B12           B16           C6           C7            C8           C9
# A.coronus   0.0011043734  0.0007714446 -1.719215e-05 -0.002850971 -0.006229655 -0.0008130481 0.0003551613
# A.inodorus -0.0003840593 -0.0008061839 -1.019612e-04 -0.004660560 -0.002160250  0.0006505090 0.0019494702
# Hybrid     -0.0005269451 -0.0015851687  1.924544e-05 -0.002268888 -0.003078627  0.0003843475 0.0025732062
# C10          C11           C12           C13          C14           D4           D11
# A.coronus  0.001735512 0.0014717377 -0.0004452489 -4.810031e-04 5.586864e-04 -0.001143296  0.0001186575
# A.inodorus 0.002139885 0.0005574137 -0.0022416617  1.478451e-05 4.731837e-04  0.006331078 -0.0022519940
# Hybrid     0.003412662 0.0016992730 -0.0026636242 -7.946901e-04 7.664099e-05  0.004235784 -0.0015931926
# D13           D15           D17
# A.coronus  -5.482047e-04 -0.0007190158  3.162957e-04
# A.inodorus -5.629808e-05 -0.0004120076  4.963305e-06
# Hybrid     -1.276460e-04 -0.0001880225 -2.624064e-04
# 
# Coefficients of linear discriminants:
#   LD1         LD2
# A2     9.832593  -71.929444
# A6   -38.608242   94.781579
# A8   115.765411  505.545641
# B6   111.095497  240.638424
# B7  -239.771816  -61.537579
# B9   189.903643 -351.414868
# C15   10.790345  139.680273
# C17 -241.411397  223.248130
# C3   -65.901586  -17.702304
# C5   -33.570296  -96.784601
# D1    16.921598  -13.772141
# D10  -90.679335  347.534517
# D12  -95.408005  278.337838
# D14 -185.326818  -73.233332
# D16 -183.216199  -67.360488
# D2   -30.612067   20.826156
# D5     4.845990 -116.197679
# D6   141.360472 -147.992897
# D7   109.252782  -54.340727
# A12  -88.144809 -321.344308
# A13  -61.645688   33.934362
# A15   84.774600 -409.550390
# A16 -101.643449 -103.622662
# B10  -34.337544   89.048762
# B12 -146.636774 -677.519939
# B16   49.981377   63.954732
# C6    62.521546   51.364112
# C7   227.236151 -111.917685
# C8    56.264381  -67.867421
# C9   195.755841  -49.147961
# C10  190.127335 -128.981811
# C11  262.226356  466.564230
# C12   83.085818  -17.592072
# C13  337.895839   29.022482
# C14  -51.510684   25.780227
# D4    66.334110  -84.363096
# D11  -75.963996   84.389859
# D13 -173.466738   -1.545171
# D15 -365.737410 -157.189432
# D17  -41.527904 -381.436263
# 
# Proportion of trace:
#   LD1    LD2 
# 0.9176 0.0824 
#These are identical in proportion of trace

linear3

plot.new()
op=par(mfrow = c(1,1), mar=c(4.2, 3.8, 0.2, 0.2))
plot(linear)
dev.off()

#Improve plot to points

tiff("LDA1Spp5-lengthCorr_SI3.tiff", width = 8, height = 6,units = "in", res = 600, compression = "lzw")
opar <- par(mar = c(4,5,1,1), oma = c(2,0,0,0))
plot(linear, xlab = "LD1: 91.76%", ylab = "LD2: 8.24%", col = dataCLASSIFICATION_lc$pop, 
     panel = points, pch = 4, lwd = 3, cex = 3, cex.lab = 1.5) 

#pop needs to be done 'as.factor'; otherwise Error "Invalid color name: 'A. inodorus'"

dev.off()

tiff("LDA1Spp5-lengthCorr_SI3ab.tiff", width = 6, height = 8,units = "in", res = 600, compression = "lzw")
op=par(mfrow = c(2,1), mar=c(4.2, 3.8, 0.2, 0.2))
plot(linear) #Compared with and without Shape indices
plot(linear2) #Identical plots
dev.off()

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
# Linear Discriminant Analysis 
# 
# 47089 samples
# 44 predictor
# 3 classes: 'A.coronus', 'A.inodorus', 'Hybrid' 
# 
# Pre-processing: scaled (44), centered (44) 
# Resampling: Cross-Validated (4 fold, repeated 100 times) 
# Summary of sample sizes: 35316, 35317, 35317, 35317, 35317, 35316, ... 
# Resampling results:
#   
#   Accuracy   Kappa    
# 0.9572363  0.8521588

confusionMatrix(LDA_results) #including SI
# Cross-Validated (4 fold, repeated 100 times) Confusion Matrix 
# 
# (entries are percentual average cell counts across resamples)
# 
# Reference
# Prediction   A.coronus A.inodorus Hybrid
# A.coronus       12.9        1.4    0.0
# A.inodorus       0.5       81.0    1.8
# Hybrid           0.0        0.6    1.8
# 
# Accuracy (average) : 0.9572

LDA_results2 <- train(pop~., method='lda', preProcess=prepr, data=dataCLASSIFICATION_lc2, trControl = trcntr)
LDA_results2 #excluding SI
confusionMatrix(LDA_results2)
# 
# Cross-Validated (4 fold, repeated 100 times) Confusion Matrix 
# 
# (entries are percentual average cell counts across resamples)
# 
# Reference
# Prediction   A.coronus A.inodorus Hybrid
# A.coronus       12.9        1.4    0.0
# A.inodorus       0.5       81.0    1.8
# Hybrid           0.0        0.6    1.8
# 
# Accuracy (average) : 0.957
#Identical accuary matrix

LDA_results3 <- train(pop~., method='lda', preProcess=prepr, data=dataCLASSIFICATION_lc3, trControl = trcntr)

LDA_results3

confusionMatrix(LDA_results3)
# Cross-Validated (4 fold, repeated 100 times) Confusion Matrix 
# 
# (entries are percentual average cell counts across resamples)
# 
# Reference
# Prediction   A.coronus A.inodorus Hybrid
# A.coronus       12.4        0.9    0.0
# A.inodorus       0.9       81.6    1.8
# Hybrid           0.0        0.5    1.8
# 
# Accuracy (average) : 0.9584

#############################
##### Discrimination using SVM ON LENGTH-CORRECTED DATA EXCLUDING SHAPE INDICES AND EXCLUDING ADDITIONAL 
###### (UNIMPORTANT) DESCRIPTORS
#############################
library(kernlab)

LDA_resultsSVM <- train(pop~., method='lssvmRadial', preProcess=prepr, data=dataCLASSIFICATION_lc3, trControl = trcntr)

LDA_resultsSVM
# Least Squares Support Vector Machine with Radial Basis Function Kernel 
# 
# 47089 samples
# 32 predictor
# 3 classes: 'A.coronus', 'A.inodorus', 'Hybrid' 
# 
# Pre-processing: scaled (32), centered (32) 
# Resampling: Cross-Validated (4 fold, repeated 100 times) 
# Summary of sample sizes: 35317, 35316, 35317, 35317, 35317, 35317, ... 
# Resampling results across tuning parameters:
#   
#   sigma        tau     Accuracy   Kappa    
# 0.009651078  0.0625  0.9677419  0.8853289
# 0.009651078  0.1250  0.9677419  0.8853289
# 0.009651078  0.2500  0.9539170  0.8298248
# 0.019832492  0.0625  0.9723502  0.9005743
# 0.019832492  0.1250  0.9723502  0.9005743
# 0.019832492  0.2500  0.9677419  0.8853289
# 0.030013907  0.0625  0.9769585  0.9182189
# 0.030013907  0.1250  0.9769585  0.9182189
# 0.030013907  0.2500  0.9769585  0.9182189
# 
# Accuracy was used to select the optimal model using the largest value.
# The final values used for the model were sigma = 0.03001391 and tau = 0.0625.

confusionMatrix(LDA_resultsSVM)
# cross-Validated (4 fold, repeated 100 times) Confusion Matrix 
# 
# (entries are percentual average cell counts across resamples)
# 
# Reference
# Prediction   A.coronus A.inodorus Hybrid
# A.coronus       12.9        0.5    0.0
# A.inodorus       0.5       82.5    1.4
# Hybrid           0.0        0.0    2.3
# 
# Accuracy (average) : 0.977

dataCLASSIFICATION_1p = dataCLASSIFICATION_lc2%>%dplyr::select(-pop)

dataCLASSIFICATION_1p$predIncl = predict(LDA_results, newdata=dataCLASSIFICATION_lc)
dataCLASSIFICATION_1p$predExcl = predict(LDA_results2, newdata=dataCLASSIFICATION_lc2)
dataCLASSIFICATION_1p$predSVM = predict(LDA_resultsSVM, newdata=dataCLASSIFICATION_lc) #pred2 = SVM

dataCLASSIFICATION_1p = cbind(dataCLASSIFICATION_1p, ID = dataCLASS$ID)
view(dataCLASSIFICATION_1p)
write.csv(dataCLASSIFICATION_1p, "fourierDescrip_withpred2.csv")

######################
## Variable importance
######################
varImp(LDA_results)
plot(varImp(LDA_results))

tiff("ImportanceSpp4_lengthcorr_3.tiff", width=6, height = 12,units = "in", res = 600, compression = "lzw")
plot(varImp(LDA_results))
dev.off()

#####################
### CAP Need library: vegan
#####################
library(vegan)

view(dataCLASSIFICATION_lc)

capresults <- capscale(dataCLASSIFICATION_lc%>%dplyr::select(-pop) ~ dataCLASSIFICATION_lc$pop) 
capresults
anova(capresults, by = "terms", step = 1000)

eig = eigenvals(capresults,model = "constrained") 
eig.ratio = eig/sum(eig) 
eig.ratio

capresults_sites <- scores(capresults)$sites[,1:2] %>% as.data.frame() %>% 
  mutate(ID2=rownames(.)) 
view(capresults_sites) #Here only CAP1 and CAP2 for each otolith ID

#Now add data
capresults_sites <-cbind(capresults_sites, data)
write.csv(capresults_sites, "capresults_sitesCorrectedForLength4.csv")

#Also add datashape_1c and data
testtosave <-cbind(capresults_sites, dataCLASSIFICATION_lc)
write.csv(testtosave, "allData_correctedForLength4.csv")

names(capresults_sites)

tiff("Capresults_Spp5_LC.tiff", width=12, heigh=6,units = "in", res = 600, compression = "lzw")
capresults_sites%>% 
  ggplot(aes(CAP1, CAP2, color=pop)) +geom_point()+
  stat_ellipse(aes(fill = pop), geom = "polygon", alpha = 0.2, level=0.95)+
  xlab(paste("CAP1:", round(eig.ratio[1]*100,2), "%"))+
  ylab(paste("CAP2:", round(eig.ratio[2]*100,2), "%"))
dev.off()

tiff("Capresults_Lat5_LC.tiff", width=12, height=6,units = "in", res = 600, compression = "lzw")
capresults_sites%>% 
  ggplot(aes(CAP1, CAP2, color=Latr))+geom_point()+
  stat_ellipse(aes(fill = Latr), geom = "polygon", alpha = 0.2, level=0.95)+
  xlab(paste("CAP1:", round(eig.ratio[1]*100,2), "%"))+
  ylab(paste("CAP2:", round(eig.ratio[2]*100,2), "%"))
dev.off()

tiff("Capresults_Area4_LC.tiff", width=12, height=6,units = "in", res = 600, compression = "lzw")
capresults_sites%>% 
  ggplot(aes(CAP1, CAP2, color=Area))+geom_point()+
  stat_ellipse(aes(fill = Area), geom = "polygon", alpha = 0.2, level=0.95)+
  xlab(paste("CAP1:", round(eig.ratio[1]*100,2), "%"))+
  ylab(paste("CAP2:", round(eig.ratio[2]*100,2), "%"))
dev.off()

tiff("Capresults_Size4_LC.tiff", width=12, height=6,units = "in", res = 600, compression = "lzw")
capresults_sites%>% 
  ggplot(aes(CAP1, CAP2, color=TLr))+geom_point()+
  stat_ellipse(aes(fill = TLr), geom = "polygon", alpha = 0.2, level=0.95)+
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

############################
##### Discrimination using LDA
#Smolinski, S., Schade, F. M., & Berg, F. (2020). Assessing the performance of statistical classifiers to discriminate fish stocks using Fourier analysis of otolith shape. Canadian Journal of Fisheries and Aquatic Sciences, 77(4), 674-683.
#############################
#############################
#need the following: library(psych); library(MASS)

library(psych)

pairs.panels(dataSHAPE[,21:25], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[data$Species],
             pch = 21,smooth = F, )


linear <- lda(pop~., dataCLASSIFICATION) #warns about variables that are collinear
linear #Gives group means wrt each F coefficient; COULD separate into different PCAs but here are only two groups
is_LDA(linear)


tiff("LDA1Spp3.tiff", width=8, heigh=5,units = "in", res = 600, compression = "lzw")
plot(linear) #Shows the discriminatory power in otolith shape to distinguish between two groups = cross-validation
dev.off()

########################################################################################################################################
###### 4-fold CROSS-VALIDATION -> uses library "caret" 
########################################################################################################################################
#Need the following: library(MASS); library(caret)
library(caret)
#install.packages("lava")

cl <- makeCluster(detectCores())
registerDoParallel(cl)
cl

### info about libraries -> MASS library
getModelInfo("lda")$lda$library

trcntr<-trainControl(method = "repeatedcv", number=4, repeats=100, allowParallel = T) #caret package
prepr=c('scale', 'center')

######################################
######################################
LDA_results<- train(pop~., method='lda', preProcess=prepr, data=dataCLASSIFICATION, trControl = trcntr)
#?train
# Original used LDA, but best according to Smolinski is support Vector Machines
#method = 'lda' (most commonly used)

#Least Squares Support Vector Machine with Radial Basis Function Kernel
#method = 'lssvmRadial'
#package required: kernlab => Smolinski

#Regularized Support Vector Machine (dual) with Linear Kernel 
#method = 'svmLinear3' 
#package(LiblineaR) needed

######################
## Accuracy, Kappa
######################
LDA_results
confusionMatrix(LDA_results)
#

######################
## Variable importance -> provides insight into which variables have most importance
######################
varImp(LDA_results)
tiff("ImportanceSpp3.tiff", width=6, heigh=12,units = "in", res = 600, compression = "lzw")
plot(varImp(LDA_results)) # How to plot them next to one another?
dev.off()
