# Script and data by Margit Wilhelm for the following manuscript:
# Otolith shape analysis as a tool for species identification and management 
# of cryptic congeners in the northern Benguela ocean warming hotspot
# Authors:
# M.R. Wilhelm 
#C.E. Jagger 
#N.M. Nghipangelwa
#B.A. Pringle 
#P.W. Shaw
#W.M. Potts
#R. Henriques 
#N.J. McKeown
#
# mwilhelm@unam.na
# 16 April 2024
# Updated 05 August 2024
#
#Based on script by: Smoliński, S., Schade, F. M., & Berg, F. (2020). Assessing the performance of statistical classifiers to discriminate fish stocks using Fourier analysis of otolith shape. Canadian Journal of Fisheries and Aquatic Sciences, 77(4), 674-683.
#
#R version 4.3.0 (2023-04-21 ucrt)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 19042)

library(tidyverse) #already for extraction
library(Momocs) #already for extraction
library(doParallel) #already for extraction
library(foreach) #already for extraction
library(broom) 
library(caret)
library(MASS)
library(scales)
library(RColorBrewer)
library(ggspatial)

# clean memory
gc()  
gc(reset=T)
rm(list = ls()) 

select <- dplyr::select #because of clashing functions "select"

theme_set(theme_bw()) #in ggplots theme of plots

############################################################################
############################################################################
path = dirname(rstudioapi::getActiveDocumentContext()$path)
path

setwd(path) #specify your working directory

load("otoKob_CJ4.RData") #From saved file, after "shape extraction", if not run at once
names(oto)

##Load the biological data associated with the otolith images
data <- read.csv(file = "kob_cj4.csv",header = TRUE )
view(data)
data$ID

#match biological data using names(oto) and data$ID
data <- data.frame(ID=names(oto)) %>% left_join(data, by = "ID")
data$ID
view(data)

gaps <- data$ID[is.na(data$X)]  #specify any variable where you expect some values missing
gaps
#present <- which(!is.na(data$X))
#oto <- oto[present] # subset of present otolith contours
#data <- data %>% filter(!is.na(X))

OTO <- Out(oto, fac = data) #Add oto with biological data
OTO
view(oto)
OTO$Species
names(OTO)

jpeg(filename = "OverviewAll4.jpg", width = 200, height = 200, res = 1000, units = "mm")
panel(OTO, col="grey", border=FALSE, names=TRUE, cex=0.3) #PLOTS the otoliths individual shapes in current folder
dev.off()

## Drop non existing levels of factor variables
OTO$fac<-OTO$fac%>%droplevels()

#Inspect the following object:
OTO
# 217 outlines, 7 classifiers

############################################################################
############################################################################
##Calibration and calculation of EFD
harmpow <-calibrate_harmonicpower_efourier(x=OTO,  nb.h = 30) #=manually specify; FIRST USE 30
#OTO is initial object & specify how many harmonics we want to extract initially

harmpow #plots sequence
harmpow$minh #prints sequence #Shows that 7 are required for 99% of variance and 17 for 99.9% variance. 

ggsave("calibrate_harmonicpower_efourier4.tiff", width=10, height=5, dpi = 600, compression = "lzw") 
#For printing figure, use only 9 (Line 90)

harm<-harmpow$minh[17]#if 99.9% of variance = 17; => specify number of harmonics as 17

####Visualisation of results = reconstructed outlines, based on different number of harmonics
rec<-calibrate_reconstructions_efourier(OTO, range = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)) 
rec$data$id<-as.factor(rec$data$id)
levels(rec$data$id)<-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)
rec

ggsave("calibrate_reconstructions_efourierAll4_99.9.tiff", width=10, heigh=5, dpi = 600, compression = "lzw")

##############################
#### Set landmark and calculate EFD
##############################
OTO <- coo_slidedirection(OTO, direction = c("right")) ##set "landmark" to the right from centroid (can be up or down)
OTO

eOTO <- efourier(OTO, nb.h = 17, smooth.it = 0, norm = T) #use 17 harmonics, don't smooth, normalization

#############################################
######## Shape descriptors
#############################################
dataSHAPE<-as.data.frame(eOTO$coe)%>%dplyr::select(-A1, -B1, -C1) #uses these three coefficients as starting points -> no additional information = will contain zero
view(dataSHAPE) #gives descriptors => can be saved as table => use this matrix in analysis in other statistical software

#############
### Histograms of the selected descriptors
#############

dataSHAPE%>%select(D1,A2,B2,C2,D2,A3,B3,C3,D3) %>% gather(Descriptor, Value)%>%ggplot(aes(Value))+geom_histogram()+facet_wrap(~Descriptor, scales="free") 
dataSHAPE%>%select(A3,B3,C3,D3,A4,B4,C4,D4) %>% gather(Descriptor, Value)%>%ggplot(aes(Value))+geom_histogram()+facet_wrap(~Descriptor, scales="free") 
dataSHAPE%>%select(A5,B5,C5,D5,A6,B6,C6,D6) %>% gather(Descriptor, Value)%>%ggplot(aes(Value))+geom_histogram()+facet_wrap(~Descriptor, scales="free") 
dataSHAPE%>%select(A7,B7,C7,D7,A8,B8,C8,D8) %>% gather(Descriptor, Value)%>%ggplot(aes(Value))+geom_histogram()+facet_wrap(~Descriptor, scales="free") 
dataSHAPE%>%select(A9,B9,C9,D9,A10,B10,C10,D10) %>% gather(Descriptor, Value)%>%ggplot(aes(Value))+geom_histogram()+facet_wrap(~Descriptor, scales="free") 
dataSHAPE%>%select(A11,B11,C11,D11,A12,B12,C12,D12) %>% gather(Descriptor, Value)%>%ggplot(aes(Value))+geom_histogram()+facet_wrap(~Descriptor, scales="free") 

#can be done for each descriptor = visual screening 
#too many bins; trouble with few datapoints (better with e.g. 50 observations)

#############
### Test for normality of descriptors
#############
dataSHAPE_Sp <- cbind(Sp=data$Species, dataSHAPE) 

dataSHAPE_SK <- dataSHAPE_Sp%>%filter(dataSHAPE_Sp$Sp == "A.inodorus")
dataSHAPE_DK <- dataSHAPE_Sp%>%filter(dataSHAPE_Sp$Sp == "A.coronus")
dataSHAPE_SK <- dataSHAPE_SK%>%select(-Sp)
dataSHAPE_DK <- dataSHAPE_DK%>%select(-Sp)

library(moments)

NORM <- dataSHAPE%>%
        gather(Descriptor, VAL)%>%
        group_by(Descriptor)%>%
        summarise(Mean = mean(VAL), 
                  Variance = var(VAL), 
                  ShapiroW = shapiro.test(VAL)$statistic,
                  ShapiroP = shapiro.test(VAL)$p.value,
                  Skew = moments::skewness(VAL), 
                  Kurt = moments::kurtosis(VAL))
view(NORM)
write.csv(NORM, "normDescriptors.csv")

view(dataSHAPE_SK)
NORM_SK <- dataSHAPE_SK%>%
  gather(Descriptor, VAL)%>%
  group_by(Descriptor)%>%
  summarise(Mean = mean(VAL), 
            Variance = var(VAL), 
            ShapiroW = shapiro.test(VAL)$statistic,
            ShapiroP = shapiro.test(VAL)$p.value,
            Skew = moments::skewness(VAL), 
            Kurt = moments::kurtosis(VAL))
view(NORM_SK)
write.csv(NORM_SK, "normDescriptorsSK.csv")

NORM_DK <- dataSHAPE_DK%>%
  gather(Descriptor, VAL)%>%
  group_by(Descriptor)%>%
  summarise(Mean = mean(VAL), 
            Variance = var(VAL), 
            ShapiroW = shapiro.test(VAL)$statistic,
            ShapiroP = shapiro.test(VAL)$p.value,
            Skew = moments::skewness(VAL), 
            Kurt = moments::kurtosis(VAL))
view(NORM_DK)
write.csv(NORM_DK, "normDescriptorsDK.csv")

#filter Descriptors with significant deviation from normality and exclude them from the dataset
sig_non_norm <- NORM%>%filter(ShapiroP<0.05)%>%dplyr::select(Descriptor)%>%unlist() 
sig_non_normSK <- NORM_SK%>%filter(ShapiroP<0.05)%>%dplyr::select(Descriptor)%>%unlist() 
sig_non_normDK <- NORM_DK%>%filter(ShapiroP<0.05)%>%dplyr::select(Descriptor)%>%unlist() 

sig_non_normSK 
# Descriptor1  Descriptor2  Descriptor3  Descriptor4  Descriptor5  Descriptor6  Descriptor7  Descriptor8 
# "A10"        "A11"        "A14"        "A17"         "A3"         "A4"         "A5"         "A7" 
# Descriptor9 Descriptor10 Descriptor11 Descriptor12 Descriptor13 Descriptor14 Descriptor15 Descriptor16 
# "A9"        "B11"         "B2"         "B3"         "B5"         "B8"         "C2"         "C4" 
# Descriptor17 
# "D8"
#ctrl-shift-c

sig_non_normDK 
# Descriptor1  Descriptor2  Descriptor3  Descriptor4  Descriptor5  Descriptor6  Descriptor7  Descriptor8 
# "A11"         "A3"        "B11"        "B13"        "B14"        "B15"        "B17"         "B2" 
# Descriptor9 Descriptor10 Descriptor11 Descriptor12 Descriptor13 Descriptor14 
# "B4"         "B5"        "C16"         "C2"         "D3"         "D9" 

length(sig_non_normSK) #For silver kob 17 (out of 65) = 48 left
length(sig_non_normDK) #For dusky kob 14 - in addition to SK = B13, B14, B15, B17, B4, C16, D3, D9 (8) = 40 left

#############
### Remove non-normal descriptors
#############
#List of those that are not normally distributed
excluded_vars <- which(colnames(dataSHAPE)%in%sig_non_normSK | colnames(dataSHAPE)%in%sig_non_normDK)
excluded_vars
#   2  3  4  6  8  9 10 13 16 17 18 19 20 23 26 28 29 30 32 33 35 47 51 56 57

dataSHAPE_e <- dataSHAPE[,-excluded_vars] #exclude those that have no normality for either pure-bred species
view(dataSHAPE_e)
#Remove 25, i.e. 40 left. 

view(dataSHAPE_e)

###########################
#Save analysis up to this point
save(dataSHAPE_e, file = "dataSHAPEAll5_norm.RData")
save(data, file = "dataBIOAll5.RData")
save(eOTO, file = "eOTOAll5.RData")

write.table(dataSHAPE_e, "kobshape5_norm.csv") #THIS IN DATA FILE, CAN BE READ IN AND MERGED WITH OTHER FILES TO PERFORM ANOVA ETC. 

##Move to Classification #Move to "anovaEtcOfShapeIndices_kob"