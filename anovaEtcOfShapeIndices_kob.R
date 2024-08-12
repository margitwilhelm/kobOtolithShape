# Script and data by Margit Wilhelm for the following manuscript:
#
# Otolith shape analysis as a tool for species identification and management 
# of cryptic congeners in the northern Benguela ocean warming hotspot
#
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
# Corresponding author: mwilhelm@unam.na
# 16 April 2024
# Updated 05 August 2024
#
#Based on script by: Smoliński, S., Schade, F. M., & Berg, F. (2020). Assessing the performance of statistical classifiers to discriminate fish stocks using Fourier analysis of otolith shape. Canadian Journal of Fisheries and Aquatic Sciences, 77(4), 674-683.
#
# R version 4.3.0 (2023-04-21 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)

library(latticeExtra)#for xyplots and coplots
library(vegan) 	#for descriptive community ecology
library(MASS) 	#for eigenvalues
library(shapeR)	#for cluster.plot

# clean memory
gc()  
gc(reset=T)
rm(list = ls()) 

path <- dirname(rstudioapi::getActiveDocumentContext()$path)
path
setwd(path)

#######KNOWN DATA SHAPE INDEX MEASUREMENTS
#file.show("kob_cj4_Measurements.csv")
alldata = read.csv("kob_cj4_Measurements.csv", header=T) 
  
alldata$FF 	= 4*pi*(alldata$OA) / (alldata$OP^2) #form factor perfectly negatively correlated with circularity
alldata$rect = alldata$OA/(alldata$OH * alldata$OL)
alldata$round = 4*alldata$OA/(pi*alldata$OL^2)
alldata$ellip = (alldata$OL - alldata$OH) / (alldata$OL + alldata$OH)
alldata$AR 	= (alldata$OL) / (alldata$OH) #Aspect ratio -> perfect correlation with ellipticity
alldata$circ = alldata$OP^2/ alldata$OA

alldata$logOA = log(alldata$OA)

names(alldata)

#########################
###test for normality
#########################

#############
library(tidyverse) 
#############

alldataSK <- alldata%>%filter(alldata$Species == "A.inodorus")
alldataDK <- alldata%>%filter(alldata$Species == "A.coronus")

NORM2sk <- alldataSK%>%
  gather(Descriptor, VAL, -ID, -Species, -TL)%>%
  group_by(Descriptor)%>%
  summarise(Mean = mean(VAL), 
            SD = sd(VAL), 
            ShapiroW = shapiro.test(VAL)$statistic,
            ShapiroP = shapiro.test(VAL)$p.value,
            Skew = moments::skewness(VAL), 
            Kurt = moments::kurtosis(VAL))

view(NORM2sk)
write.csv(NORM2sk, "normShapeIndicesSK.csv")

NORM2dk <- alldataDK%>%
  gather(Descriptor, VAL, -ID, -Species, -TL)%>%
  group_by(Descriptor)%>%
  summarise(Mean = mean(VAL), 
            SD = sd(VAL), 
            Var = var(VAL),
            ShapiroW = shapiro.test(VAL)$statistic,
            ShapiroP = shapiro.test(VAL)$p.value,
            Skew = moments::skewness(VAL), 
            Kurt = moments::kurtosis(VAL))

view(NORM2dk)
write.csv(NORM2dk, "normShapeIndicesDK.csv")

#filter Shape indices with significant deviation from normality and exclude them from the dataset
sig_non_norm2sk <- NORM2sk%>%filter(ShapiroP<0.05)%>%dplyr::select(Descriptor)%>%unlist() 
sig_non_norm2dk <- NORM2dk%>%filter(ShapiroP<0.05)%>%dplyr::select(Descriptor)%>%unlist() 

sig_non_norm2sk #REMOVE THESE FOLLOWING VARIABLES BECAUSE OF NON-NORMALITY WITH either species
#Descriptor1 Descriptor2 Descriptor3 Descriptor4 Descriptor5 
#"OA"       "OCD"        "OH"        "OL"        "OP" 

sig_non_norm2dk #REMOVE THESE FOLLOWING VARIABLES BECAUSE OF NON-NORMALITY WITH either species
#Descriptor 
#"rect" 

#############
### Remove non-normal descriptors
#############
#List of those that are not normally distributed
excluded_vars <- which(colnames(alldata)%in%sig_non_norm2sk | colnames(alldata)%in%sig_non_norm2dk)
excluded_vars
# 4 6 7 8 9 11

alldata_e <- alldata[,-excluded_vars] #exclude those that have no normality
view(alldata_e)
#Keep 7

write.table(alldata_e, "shapeIndiceskob2_norm.csv")

#############
library(tidyverse) 
#############

names(alldata_e)

##### DO COMPARISONS BEFORE LENGTH CORRECTION
m1=aov(OCD.OH ~ Species, data = alldata_e) #Can do this individually as before, or together 
summary(m1) # display Type I ANOVA table
drop1(m1,~.,test = "F") # type III SS and F Tests 
TukeyHSD(m1) 

####################################################
####Need library(vegan) for descriptive community ecology
####Need library(MASS) for eigenvals. And need library(shapeR) for cluster.plot
####Canonical analysis of principal coordinates
##NORMAL PCA - FIRST JUST ON UNSTANDARDIZED CO-ORDINATES - NOT SUCH GOOD SEPARATION
library(MASS)
library(vegan) #for capscale

cap.res3 = capscale(alldata_e[, 4:10] ~ Species, alldata_e)

plot(cap.res3)#Doesn't show much deviation. Need to remove correct for length-effects and combine with Fourier descriptors

eig3 = eigenvals(cap.res3, model = "constrained")
eig.ratio3 = eig3/sum(eig3)
eig.ratio3
eig3
cap.res3

#############################
#############################
#need the following: library(psych); library(MASS)
library(psych)

pairs.panels(alldata_e[,3:10], #package(psych) only use first four Fourier coefficients - visualise
             gap = 0,
             bg = c("red", "green")[shapeAge3$Species],
             pch = 21,smooth = F, ) #Some are co-linear, but first remove TL relationship

dataCLASSIFICATION2 <- alldata_e

linear <- lda(Species~., dataCLASSIFICATION2) #warns about variables that are collinear !! 
linear #Gives group means wrt each F coefficient; COULD separate into different PCAs but here only two groups

plot(linear) #Does bring out hybrids

#Shows the discriminatory power in otolith shape to distinguish between groups = cross-validation

#######
##STANDARDIZATION COEFFICIENTS USING LM
###################
## Remove allometric effects (length-correction)
###################

dataSHAPE_lc2 <- dataCLASSIFICATION2

View(dataSHAPE_lc2)

#test with ANCOVA Length significance in interaction with pop 
#function gather in tidyr package
library(tidyr)
library(broom)

ANCOVA2 <- dataSHAPE_lc2%>%gather(Descriptor, VAL,-ID, -TL, -Species)%>%group_by(Descriptor)%>%
  do(tidy(aov(lm(VAL ~ TL*Species, data=.))))
#view(ANCOVA2)
write.csv(ANCOVA2, "ancova_indices.csv")

#filter Descriptors with significant Length effect
sig_length_ <- ANCOVA2%>%filter(term=="TL"&p.value<0.05)%>%dplyr::select(Descriptor)%>%unlist() 
sig_length_ #4
#Descriptor1 Descriptor2 Descriptor3 Descriptor4 
#   "FF"    "OCD.OH"      "circ"     "logOA" 

#filter Descriptors with significant length and pop interaction and exclude them from the dataset
sig_pop_interaction_ <- ANCOVA2%>%filter(Descriptor%in%sig_length_)%>%filter(term=="TL:Species")%>%filter(p.value<0.05)%>%dplyr::select(Descriptor)%>%unlist() 
sig_pop_interaction_ 
#1 = logOA 
length(sig_pop_interaction_) #1

names(dataSHAPE_lc2)

dataSHAPE_lcB2 <- dataSHAPE_lc2[,-which(colnames(dataSHAPE_lc2)%in%as.vector(sig_length_))]
names(dataSHAPE_lcB2) #3 columns Those that don't need length conversion (round, ellip, AR) including ID, Species and TL

dataSHAPE_lcA2 <- dataSHAPE_lc2[,-which(colnames(dataSHAPE_lc2)%in%as.vector(sig_pop_interaction_))] #only those no sign pop interaction
names(dataSHAPE_lcA2) #Those that don't have sign population interaction (so logOA now dropped) 
#incl. ID, Species, TL

#Other variables will be corrected - 
sig_length_ <- sig_length_[-which(sig_length_%in%sig_pop_interaction_)] 
sig_length_ 
length(sig_length_) #3
# Descriptor1 Descriptor2 Descriptor3 
# "FF"    "OCD.OH"      "circ" 

LM <- dataSHAPE_lcA2%>%gather(Descriptor,VAL,-TL,-Species, -ID)%>%filter(Descriptor%in%sig_length_)%>%
  group_by(Descriptor)%>%
  do(tidy((lm(VAL ~ TL, data=.))))

LM <- LM%>%filter(Descriptor%in%sig_length_)%>%filter(term=="TL")%>%dplyr::select(Descriptor,estimate)

view(LM)
write.csv(LM, "LM_shapeInd.csv")

#Check relationships of randomly selected variable on plot
dataSHAPE_lcA2%>%ggplot(aes(TL, OCD.OH))+geom_point()+geom_smooth(method="lm")

id <- which(colnames(dataSHAPE_lcA2)%in%sig_length_)
id # Columns 4, 5, 9

dataSHAPE_lc12 <- dataSHAPE_lcA2[,id] #Choose only columns 4, 5, and 9 that need to be corrected (now has no ID or length)
names(dataSHAPE_lc12)

dataSHAPE_lc12$ID <- dataSHAPE_lc2$ID
dataSHAPE_lc12$length <- dataSHAPE_lc2$TL #Add ID and TL again

dataSHAPE_lcB2 <- dataSHAPE_lcB2%>%dplyr::select(-ID, -TL, -Species) #no correction needed - picked in line 202
names(dataSHAPE_lcB2) #Those that need no correction with ID, TL and Species removed
names(dataSHAPE_lc12) #Those that need correction (with ID and length included), but are not yet corrected

dataSHAPE_lc12 <- dataSHAPE_lc12%>%gather(Descriptor, Value, -length, -ID) #Made into long data frame
dataSHAPE_lc12 <- dataSHAPE_lc12%>%left_join(LM)
dataSHAPE_lc12$Value <- dataSHAPE_lc12$Value-dataSHAPE_lc12$estimate*dataSHAPE_lc12$length
view(dataSHAPE_lc12) #now they are corrected - but still in long data frame

dataSHAPE_lc12 <- dataSHAPE_lc12%>%dplyr::select(Descriptor, Value, ID, length)%>%
  spread(Descriptor, Value)#Now they become columns again

dataSHAPE_lc2 <- cbind(dataSHAPE_lc12, dataSHAPE_lcB2)
view(dataSHAPE_lc2) #All combined corrected with ID and length
write.csv(dataSHAPE_lc2, "shapeIndices_length-corr_kob.csv")

#Check relationships of randomly selected variable on plot
dataSHAPE_lc2%>%ggplot(aes(length, OCD.OH))+geom_point()+geom_smooth(method="lm") #length-corrected (flat line)
names(dataSHAPE_lc2)

dataSHAPE_lc2x <- dataSHAPE_lc2%>%dplyr::select(-length, -ID) #Remove length and ID for LDA

View(dataSHAPE_lc2x) ### New matrix of shape indices corrected for the length effects  
names(dataSHAPE_lc2x)

pairs.panels(dataSHAPE_lc2x[,1:6], #package(psych) - visualise
             gap = 0,
             bg = c("red", "green")[shapeAge3$Species],
             pch = 21,smooth = F, )

names(dataSHAPE_lc2)

###############AFTER STANDARDIZATION 
###############REDO ANALYSES
dataSHAPE_new2 <- dataSHAPE_lc2%>%cbind(pop=dataCLASSIFICATION2$Species)

shapeAge3 = dataSHAPE_new2 #length-corrected shape indices
names(shapeAge3)

#################################
#################################
### Boxplots & VIOPLOTS
library(vioplot)

#OCD, Stcomp, StFF
tiff("4boxplots3.tiff", width=6, height = 12, units = "in", res = 600, compression = "lzw")
op=par(mfrow = c(4,1), mar=c(4.2, 3.8, 3.8, 0.2)) #bottom, left, top, right
boxplot(round ~ pop, data=shapeAge3, xlab = "", ylab = "Otolith roundness", main = "A") 
boxplot(OCD.OH ~ pop, data=shapeAge3, xlab = "", ylab = "Standardised OCD-OH ratio", main ="B")
boxplot(circ ~ pop, data=shapeAge3, xlab = "", ylab = "Standardised otolith circularity", main ="C")
boxplot(ellip ~ pop, data = shapeAge3, xlab = "Species", ylab = "Otolith ellipticity", main = "D")
dev.off()

#VIOPLOTS OF LENGHT-CORRECTED DATA
tiff("4vioplots.tiff", width = 6, height = 12, units = "in", res = 600, compression = "lzw")
op=par(mfrow = c(4,1), mar=c(2, 4, 1, 0)) #bottom, left, top, right
vioplot(round ~ pop, data =shapeAge3, xlab= "", ylab = "Otolith roundness", 
        main = "A. Otolith roundness", 
        cex.axis = 1.5,
        col = c("lightgreen", "lightblue", "palevioletred"),
        lty = 1,
        lwd = 1, 
        pchMed = 19,
        colMed = "white", 
        cex = 1, #This is of the median dot
        cex.lab = 1.5, 
        cex.main = 1.5) 
    add_outliers(unlist(shapeAge3$round), shapeAge3$pop,
                     col = "black", fill = "red", bars = "grey85")
    legend("topright", legend=c("A. coronus", "A. inodorus", "Hybrids"),
               fill=c("lightgreen", "lightblue", "palevioletred"), cex = 1.5)
#    add_labels(unlist(shapeAge3$round), shapeAge3$pop, height = 0, cex = 1.5) #This adds sample sizes

vioplot(OCD.OH ~ pop, data=shapeAge3, xlab = "", ylab = "Standardised OCD-OH ratio", 
        main ="B. OCD-OH ratio", 
        cex.axis = 1.5, 
        col = c("lightgreen", "lightblue", "palevioletred"),
        lty = 1,
        lwd = 1, 
        pchMed = 19,
        colMed = "white", 
        cex = 1, #This is of the median dot
        cex.main = 1.5) 
  add_outliers(unlist(shapeAge3$OCD.OH), shapeAge3$pop,
             col = "black", fill = "red", bars = "grey85")
  
vioplot(circ ~ pop, data=shapeAge3, xlab = "", ylab = "Standardised otolith circularity", 
        main ="C. Otolith circularity", 
        cex.axis = 1.5, 
        col = c("lightgreen", "lightblue", "palevioletred"),
        lty = 1,
        lwd = 1, 
        pchMed = 19,
        colMed = "white", 
        cex = 1, #This is of the median dot
        cex.main = 1.5) 
add_outliers(unlist(shapeAge3$circ), shapeAge3$pop,
             col = "black", fill = "red", bars = "grey85")

vioplot(ellip ~ pop, data = shapeAge3, xlab = "Species", ylab = "Otolith ellipticity", 
        main = "D. Otolith ellipticity", 
        cex.axis = 1.5, 
        col = c("lightgreen", "lightblue", "palevioletred"),
        lty = 1,
        lwd = 1, 
        pchMed = 19,
        colMed = "white", 
        cex = 1, #This is of the median dot
        cex.main = 1.5) 
add_outliers(unlist(shapeAge3$ellip), shapeAge3$pop,
             col = "black", fill = "red", bars = "grey85")
dev.off()

names(alldata)
#Plot uncorrected
#VIOPLOTS OF RAW DATA
tiff("4vioplots_raw.tiff", width = 6, height = 12, units = "in", res = 600, compression = "lzw")

library(svglite)
svglite("4vioplots_raw.svg", width = 6, height = 12)

op=par(mfrow = c(4,1), mar=c(2, 4, 1, 0)) #bottom, left, top, right
vioplot(round ~ Species, data = dataCLASSIFICATION2, xlab= "", ylab = "Otolith roundness", 
        main = "A. Otolith roundness", 
        cex.axis = 1.5,
        col = c("lightgreen", "lightblue", "palevioletred"),
        lty = 1,
        lwd = 1, 
        pchMed = 19,
        colMed = "white", 
        cex = 1, #This is of the median dot
        cex.lab = 1.5, 
        cex.main = 1.5) 
add_outliers(unlist(dataCLASSIFICATION2$round), dataCLASSIFICATION2$Species,
             col = "black", fill = "red", bars = "grey85")

legend("topright", legend=c("A. coronus", "A. inodorus", "Hybrids"),
       fill=c("lightgreen", "lightblue", "palevioletred"), cex = 1.5)
#    add_labels(unlist(shapeAge3$round), shapeAge3$pop, height = 0, cex = 1.5) #This adds sample sizes

vioplot(OCD.OH ~ Species, data=dataCLASSIFICATION2, xlab = "", ylab = "Standardised OCD-OH ratio", 
        main ="B. OCD-OH ratio", 
        cex.axis = 1.5, 
        col = c("lightgreen", "lightblue", "palevioletred"),
        lty = 1,
        lwd = 1, 
        pchMed = 19,
        colMed = "white", 
        cex = 1, #This is of the median dot
        cex.main = 1.5)

add_outliers(unlist(dataCLASSIFICATION2$OCD.OH), dataCLASSIFICATION2$Species,
             col = "black", fill = "red", bars = "grey85")

vioplot(circ ~ Species, data=dataCLASSIFICATION2, xlab = "", ylab = "Standardised otolith circularity", 
        main ="C. Otolith circularity", 
        cex.axis = 1.5, 
        col = c("lightgreen", "lightblue", "palevioletred"),
        lty = 1,
        lwd = 1, 
        pchMed = 19,
        colMed = "white", 
        cex = 1, #This is of the median dot
        cex.main = 1.5) 
add_outliers(unlist(dataCLASSIFICATION2$circ), dataCLASSIFICATION2$Species,
             col = "black", fill = "red", bars = "grey85")

vioplot(ellip ~ Species, data = dataCLASSIFICATION2, xlab = "Species", ylab = "Otolith ellipticity", 
        main = "D. Otolith ellipticity", 
        cex.axis = 1.5, 
        col = c("lightgreen", "lightblue", "palevioletred"),
        lty = 1,
        lwd = 1, 
        pchMed = 19,
        colMed = "white", 
        cex = 1, #This is of the median dot
        cex.main = 1.5) 
add_outliers(unlist(dataCLASSIFICATION2$ellip), dataCLASSIFICATION2$Species,
             col = "black", fill = "red", bars = "grey85")
dev.off()



#OLD BOXPLOTS
tiff("6boxplots.tiff", width=6, height = 12, units = "in", res = 600, compression = "lzw")
op=par(mfrow = c(6,1), mar=c(4.2, 3.8, 0.2, 0.2))
boxplot(circ ~ Species, data=alldata, xlab = "", ylab = "Otolith circularity", main = "A")
boxplot(FF ~ Species, data=alldata, xlab = "", ylab = "Otolith form-factor", main = "B")
boxplot(OCD.OH ~ Species, data=alldata, xlab = "", ylab = "OCD-OH ratio", main ="C")
boxplot(round ~ pop, data=shapeAge3, xlab = "", ylab = "Otolith roundness", main = "D") 
boxplot(ellip ~ Species, data = alldata, xlab = "", ylab = "Otolith ellipticity", main = "E")
boxplot(AR ~ Species, data = alldata, xlab = "Species", ylab = "Otolith aspect ratio", main = "F")
dev.off()

tiff("4violots_uncorr.tiff", width=6, height = 12, units = "in", res = 600, compression = "lzw")
op=par(mfrow = c(6,1), mar=c(4.2, 3.8, 0.2, 0.2))
vioplot(circ ~ Species, data=alldata, xlab = "", ylab = "Otolith circularity", main = "A", cex.axis = 1.5, cex.lab = 1.5)
vioplot(OCD.OH ~ Species, data=alldata, xlab = "", ylab = "OCD-OH ratio", main ="B", cex.axis = 1.5)
vioplot(round ~ Species, data=alldata, xlab = "", ylab = "Otolith roundness", main = "C", cex.axis = 1.5) 
vioplot(ellip ~ Species, data = alldata, xlab = "Species", ylab = "Otolith ellipticity", main = "D", cex.axis = 1.5, cex.lab = 2)
dev.off()

####################################################
###### PLOT OCD AND OCD-OH against Fish length
######
##### XY plots for uncorrected OCD mm and OCD-OH against fish total length
library(latticeExtra)
names(alldata)

shapeAge3_uncorrected_inodorus = subset(alldata, alldata$Species=="A.inodorus")
shapeAge3_uncorrected_coronus = subset(alldata, alldata$Species=="A.coronus")
shapeAge3_uncorrected_hybrid = subset(alldata, alldata$Species=="Hybrid")

f_key = list(x = .95, y = 0.01, corner = c(1, 0),
            text = list(c("A. inodorus", 
                          "A. coronus",
                          "Hybrids")), 
            points = list(type = c("p", "p", "p"), 
                          col = c("red", "#009E73", "black"),
                          pch = c(1, 5, 3), 
#                          lwd = c(3, 3, 3), 
                          lty = c(5, 5, 5)))


p1a = xyplot(OCD ~ TL, data=shapeAge3_uncorrected_inodorus, xlab = "Fish total length (cm)", ylab = "OCD (mm)", main = "A",
             col.lab = "black", cex.lab = 20, 
             type = 'p', pch = 1, lwd = 3, col = "red", cex = 1, font.face = "bold",
             key = f_key, 
             scales=list(x=list(limits=c(30,100), at = seq(30, 100, by=10)), 
             y=list(limits=c(0,6), at = seq(0, 6, by = 1) ))
             ) #inodorus
p1a

p1b <- xyplot(OCD ~ TL, data=shapeAge3_uncorrected_coronus, type='p', pch = 5, col = "#009E73") #coronus
p1c <- xyplot(OCD ~ TL, data=shapeAge3_uncorrected_hybrid, type='p', pch = 3, col = "black") #Hybrid


p2a = xyplot(OCD.OH ~ TL, data=shapeAge3_uncorrected_inodorus, xlab = "Fish total length (cm)", ylab = "OCD / OH", main = "B",
             col.lab = "black", cex.lab = 20,
             type='p', pch = 1, lty = 5, lwd = 3, col = "red", key = f_key, 
             scales=list(x=list(limits=c(30,100), at = seq(30, 100, by=10)), 
                         y=list(limits=c(0.1,0.55), at = seq(0.1, 0.55, by = 0.05) ))
             ) #inodorus
p2b <- xyplot(OCD.OH ~ TL, data=shapeAge3_uncorrected_coronus, type='p', pch = 5, lty = 5, lwd = 3, col = "#009E73") #coronus
p2c <- xyplot(OCD.OH ~ TL, data=shapeAge3_uncorrected_hybrid, type='p', pch = 3, lty = 5, lwd = 3, col = "black") #Hybrid

A = p1a + p1b + p1c
B = p2a + p2b + p2c
A
B

tiff("OCD-TL.tiff", width=6, height = 8,units = "in", res = 600, compression = "lzw")
par(mar=c(4.2, 3.8, 0.2, 0.2))#bottom, left, top, right
print(A, split=c(1,1,1,2), more=TRUE)
print(B, split=c(1,2,1,2), more=TRUE)
dev.off()


library(svglite)
svglite("OCD-TL.svg", width = 6, height = 12)
par(mar=c(4.2, 3.8, 0.2, 0.2))#bottom, left, top, right
print(A, split=c(1,1,1,2), more=TRUE)
print(B, split=c(1,2,1,2), more=TRUE)
dev.off()


####################################################
##ANOVA ONLY DO THIS NOW FOR STANDARDIZED VALUES
#names(shapeAge3)
#view(shapeAge3)
ANOVA <- shapeAge3%>%gather(Descriptor, VAL,-length, -pop, -ID)%>%group_by(Descriptor)%>%
  do(tidy(aov(lm(VAL ~ pop, data=.))))
view(ANOVA)
write.csv(ANOVA, "anova_indices.csv")


m1 = aov(OCD.OH ~ pop, data = shapeAge3)
summary(m1) # display Type I ANOVA table
TukeyHSD(m1) #

m2 = aov(circ ~ pop, data = shapeAge3)
summary(m2) # display Type I ANOVA table
TukeyHSD(m2) #

m3 = aov(round ~ pop, data = shapeAge3)
summary(m3) # display Type I ANOVA table
TukeyHSD(m3) #
