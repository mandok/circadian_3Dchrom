#####################--------HiC--------#####################
##All analyses were carried out in R version 3.5.1 (2018-07-02)
#####################---Eigenvector scatterplots +ANOVA of compartments---#############
library(CHiCAGO)
library(gplots)
library(ggplot2)
library(reshape2)
library(hexbin)
library(gridExtra)C1569
library(CHiCAGO)
library(gplots)
library(ggplot2)
library(reshape2)
library(hexbin)
library(gridExtra)
setwd("HiC_analysis_Babraham/Juicebox/")
#readRDS("ZT0All_step2.Rds")
#setwd("/hdisk7/mandok/Temporal/HiC")
eig_06<-read.csv("ZT06sharedbins100kbeigenKR.bed", header = F, sep="\t", colClasses=c("NULL","NULL","NULL", NA, NA))
eig_012<-read.csv("ZT012sharedbins100kbeigenKR.bed", header = F, sep="\t",colClasses=c("NULL","NULL","NULL", NA, NA))
eig_018<-read.csv("ZT018sharedbins100kbeigenKR.bed", header = F, sep="\t",colClasses=c("NULL","NULL","NULL", NA, NA))
eig_612<-read.csv("ZT612sharedbins100kbeigenKR.bed", header = F, sep="\t",colClasses=c("NULL","NULL","NULL", NA, NA))
eig_618<-read.csv("ZT618sharedbins100kbeigenKR.bed", header = F, sep="\t",colClasses=c("NULL","NULL","NULL", NA, NA))
eig_1218<-read.csv("ZT1218sharedbins100kbeigenKR.bed", header = F, sep="\t",colClasses=c("NULL","NULL","NULL", NA, NA))



p1<-ggplot(eig_06) + geom_hex(aes(V4, V5), bins = 200) +
  theme_minimal()+xlab("ZT0")+ylab("ZT6")+ggtitle("Eigenvector")+
  scale_fill_gradientn("Counts", colours = rev(rainbow(10, end = 4/6)))

p2<-ggplot(eig_012) + geom_hex(aes(V4, V5), bins = 200) +
  theme_minimal()+xlab("ZT0")+ylab("ZT12")+ggtitle("Eigenvector")+
  scale_fill_gradientn("Counts", colours = rev(rainbow(10, end = 4/6)))

p3<-ggplot(eig_018) + geom_hex(aes(V4, V5), bins = 200) +
  theme_minimal()+xlab("ZT0")+ylab("ZT18")+ggtitle("Eigenvector")+
  scale_fill_gradientn("Counts", colours = rev(rainbow(10, end = 4/6)))

p4<-ggplot(eig_612) + geom_hex(aes(V4, V5), bins = 200) +
  theme_minimal()+xlab("ZT6")+ylab("ZT12")+ggtitle("Eigenvector")+
  scale_fill_gradientn("Counts", colours = rev(rainbow(10, end = 4/6)))

p5<-ggplot(eig_618) + geom_hex(aes(V4, V5), bins = 200) +
  theme_minimal()+xlab("ZT6")+ylab("ZT18")+ggtitle("Eigenvector")+
  scale_fill_gradientn("Counts", colours = rev(rainbow(10, end = 4/6)))

p6<-ggplot(eig_1218) + geom_hex(aes(V4, V5), bins = 200) +
  theme_minimal()+xlab("ZT12")+ylab("ZT18")+ggtitle("Eigenvector")+
  scale_fill_gradientn("Counts", colours = rev(rainbow(10, end = 4/6)))


#svglite::svglite("eigenvectors_alltimepoints_hexagon_200bins.svg")
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3)
#dev.off()


#svglite::svglite("eigenvectors_alltimepoints.svg")
#https://stackoverflow.com/questions/14271584/r-legend-for-color-density-scatterplot-produced-using-smoothscatter
fudgeit <- function(){
  xm <- get('xm', envir = parent.frame(1))
  ym <- get('ym', envir = parent.frame(1))
  z  <- get('dens', envir = parent.frame(1))
  colramp <- get('colramp', parent.frame(1))
  fields::image.plot(xm,ym,z, col = colramp(256), legend.only = T, add =F)
}

par(mar = c(5,4,4,5) + .1)
smoothScatter(x, nrpoints = 0, postPlotHook = fudgeit)

par(mfrow=(c(2,3)))
par(mar = c(5,4,4,5) + 1)
smoothScatter(eig_06[,1], eig_06[,2],
              colramp = colorRampPalette(colors = c("white", "blue", "cyan", "yellow","orange", "red")),
              xlab = "ZT0", ylab="ZT6", main="Eigenvector", nrpoints = 0, postPlotHook = fudgeit)

smoothScatter(eig_012[,1], eig_012[,2],
              colramp = colorRampPalette(colors = c("white", "blue", "cyan", "yellow","orange", "red")),
              xlab = "ZT0", ylab="ZT12", main="Eigenvector",
              nrpoints = 0)


smoothScatter(eig_018[,1], eig_018[,2],
              colramp = colorRampPalette(colors = c("white", "blue", "cyan", "yellow","orange", "red")),
              xlab = "ZT0", ylab="ZT18", main="Eigenvector",
              nrpoints = 0)

smoothScatter(eig_612[,1], eig_612[,2],
              colramp = colorRampPalette(colors = c("white", "blue", "cyan", "yellow","orange", "red")),
              xlab = "ZT6", ylab="ZT12", main="Eigenvector",
              nrpoints = 0)


smoothScatter(eig_618[,1], eig_618[,2],
              colramp = colorRampPalette(colors = c("white", "blue", "cyan", "yellow","orange", "red")),
              xlab = "ZT6", ylab="ZT18", main="Eigenvector",
              nrpoints = 0)


smoothScatter(eig_1218[,1], eig_1218[,2],
              colramp = colorRampPalette(colors = c("white", "blue", "cyan", "yellow","orange", "red")),
              xlab = "ZT12", ylab="ZT18", main="Eigenvector",
              nrpoints = 0)
#dev.off()
######################3
####Heatmaps
#Merged replicates
eig_061218<-read.csv("HiC_analysis_Babraham/Juicebox/Compartment_changes/Mergedreplicates_Comparmentschanges_ZT061218sharedbins100kbeigenKR_allcombinations", header = F, sep="\t", colClasses=c("NULL","NULL","NULL", NA, NA, NA, NA))
colnames(eig_061218)<-c("ZT0", "ZT6", "ZT12", "ZT18")
#svglite::svglite("/Users/andoku01/Dropbox (Cambridge University)/IFC/Figures2_Masami/HiC/Compartments/Mergedreplicates_changesincompartments.svg")
heatmap.2(as.matrix(eig_061218), key.xlab = "PCA1", key.title = "", Colv = F, trace="none", main="Compartment", density.info="none",col = colorRampPalette(c("blue","white","red")), labRow = NA)
#dev.off()

#Separated replicates
eig_061218<-read.csv("/Users/andoku01/Dropbox (Personal)/Masami/HiC_analysis_Babraham/Juicebox/Compartment_NOchanges/NochangesZT061218REPLICATESsharedbins100kbeigenKR.bed", header = F, sep="\t", colClasses=c("NULL","NULL","NULL", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
colnames(eig_061218)<-c("ZT0A","ZT0B", "ZT0C",  "ZT6A","ZT6B", "ZT6C", "ZT12A","ZT12B", "ZT12C","ZT18A","ZT18B", "ZT18C")
#svglite::svglite("/Users/andoku01/Dropbox (Cambridge University)/IFC/Figures2_Masami/HiC/Compartments/Replicatesseparated_nochangesincompartments.svg")
heatmap.2(as.matrix(eig_061218), key.xlab = "PCA1", key.title = "", Colv = F, trace="none", main="Compartment", density.info="none",col = colorRampPalette(c("blue","white","red")), labRow = NA)
#dev.off()

#Filter only different compartmnets

t<-eig_061218[1:25546,]<0 #B compartments 
tt<-eig_061218[1:25546,]>0 #A compartments 
t1<-apply(t, 1, sum) #Is the whole row B compartments
tt1<-apply(tt, 1, sum) #Is the whole row A compartments
t2<-eig_061218[(t1!=4) & (tt1!=4),] 
colnames(t2)<-c("ZT0", "ZT6", "ZT12", "ZT18")
heatmap.2(as.matrix(t2), key.xlab = "PCA1", key.title = "", Colv = F, trace="none", main="Compartment", density.info="none",col = colorRampPalette(c("blue","white","red")), labRow = NA)

####################################################################################
######################ANOVA
#Static Merged replicates---------------------------------------------------------------
eig_061218<-read.csv("/Users/andoku01/Dropbox (Personal)/Masami/HiC_analysis_Babraham/Juicebox/Compartment_NOchanges/NochangesZT061218sharedbins100kbeigenKR.bed", header = F, sep="\t", colClasses=c("NULL","NULL","NULL", NA, NA, NA, NA))
colnames(eig_061218)<-c("ZT0", "ZT6", "ZT12", "ZT18")

eig_061218anova<-as.data.frame(rbind(cbind(eig_061218[,1], colnames(eig_061218)[1]),
                                     cbind(eig_061218[,2], colnames(eig_061218)[2]),
                                     cbind(eig_061218[,3], colnames(eig_061218)[3]),
                                     cbind(eig_061218[,4], colnames(eig_061218)[4])))
eig_061218anova[,1]<-as.numeric(levels(eig_061218anova[,1]))[eig_061218anova[,1]]
colnames(eig_061218anova)<-c("pca1", "timepoint")
# Compute the analysis of variance
res.aov <- aov(pca1 ~ timepoint, data = eig_061218anova)
# Summary of the analysis
summary(res.aov)
#Df Sum Sq   Mean Sq F value Pr(>F)
#timepoint       3   0.00 0.0004232   0.585  0.625
#Residuals   84564  61.22 0.0007239 
#Normality assuption: Not normal
plot(res.aov, 2)
# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test
shapiro.test(x = sample(aov_residuals, 5000, replace = F) )
#Non-parametric alternative to one-way ANOVA test
kruskal.test(pca1 ~ timepoint, data = eig_061218anova)
#Kruskal-Wallis rank sum test
#
#data:  pca1 by timepoint
#Kruskal-Wallis chi-squared = 3.7245, df = 3, p-value = 0.2928

#CCC Merged replicates---------------------------------------------------------------
eig_061218<-read.csv("/Users/andoku01/Dropbox (Personal)/Masami/HiC_analysis_Babraham/Juicebox/Compartment_changes/Mergedreplicates_Comparmentschanges_ZT061218sharedbins100kbeigenKR_allcombinations", header = F, sep="\t", colClasses=c("NULL","NULL","NULL", NA, NA, NA, NA))
colnames(eig_061218)<-c("ZT0", "ZT6", "ZT12", "ZT18")

eig_061218anova<-as.data.frame(rbind(cbind(eig_061218[,1], colnames(eig_061218)[1]),
                                     cbind(eig_061218[,2], colnames(eig_061218)[2]),
                                     cbind(eig_061218[,3], colnames(eig_061218)[3]),
                                     cbind(eig_061218[,4], colnames(eig_061218)[4])))
eig_061218anova[,1]<-as.numeric(levels(eig_061218anova[,1]))[eig_061218anova[,1]]
colnames(eig_061218anova)<-c("pca1", "timepoint")
# Compute the analysis of variance
res.aov <- aov(pca1 ~ timepoint, data = eig_061218anova)
# Summary of the analysis
#summary(res.aov)
#Df Sum Sq Mean Sq F value Pr(>F)    
#timepoint       3  4.145  1.3816    3502 <2e-16 ***
#  Residuals   17612  6.947  0.0004                   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Normality assuption: Checked
plot(res.aov, 2)
#####
#Static Separate replicates---------------------------------------------------------------
eig_061218<-read.csv("HiC_analysis_Babraham/Juicebox/Compartment_NOchanges/NochangesZT061218REPLICATESsharedbins100kbeigenKR.bed", header = F, sep="\t", colClasses=c("NULL","NULL","NULL", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
colnames(eig_061218)<-c("ZT0A","ZT0B", "ZT0C",  "ZT6A","ZT6B", "ZT6C", "ZT12A","ZT12B", "ZT12C","ZT18A","ZT18B", "ZT18C")

eig_061218anova<-as.data.frame(rbind(cbind(eig_061218[,1], "ZT0", "A"),
                                     cbind(eig_061218[,2],  "ZT0", "B"),
                                     cbind(eig_061218[,3],  "ZT0", "C"),
                                     cbind(eig_061218[,4],  "ZT6", "A"),
                                     cbind(eig_061218[,5], "ZT6", "B"),
                                     cbind(eig_061218[,6], "ZT6", "C"),
                                     cbind(eig_061218[,7], "ZT12", "A"),
                                     cbind(eig_061218[,8], "ZT12", "B"),
                                     cbind(eig_061218[,9], "ZT12", "C"),
                                     cbind(eig_061218[,10], "ZT18", "A"),
                                     cbind(eig_061218[,11], "ZT18", "B"),
                                     cbind(eig_061218[,12], "ZT18", "C")
                                     ))
eig_061218anova[,1]<-as.numeric(levels(eig_061218anova[,1]))[eig_061218anova[,1]]
colnames(eig_061218anova)<-c("pca1", "timepoint", "replicate")
# Compute the analysis of variance
res.aov <- aov(pca1 ~ timepoint * replicate, data = eig_061218anova)
# Summary of the analysis
summary(res.aov)
#Df Sum Sq   Mean Sq F value Pr(>F)
#timepoint       3   0.00 0.0004232   0.585  0.625
#Residuals   84564  61.22 0.0007239 
#Normality assuption: Not normal
plot(res.aov, 2)
# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test
shapiro.test(x = sample(aov_residuals, 5000, replace = F) )
#Non-parametric alternative to one-way ANOVA test
kruskal.test(pca1 ~ timepoint, data = eig_061218anova)
#Kruskal-Wallis rank sum test
#
#data:  pca1 by timepoint
#Kruskal-Wallis chi-squared = 3.7245, df = 3, p-value = 0.2928

#CCC Merged replicates---------------------------------------------------------------
eig_061218<-read.csv("/Users/andoku01/Dropbox (Personal)/Masami/HiC_analysis_Babraham/Juicebox/Compartment_changes/ZT061218replicates_changeincompartments_allcombinations", header = F, sep="\t", colClasses=c("NULL","NULL","NULL", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
colnames(eig_061218)<-c("ZT0A","ZT0B", "ZT0C",  "ZT6A","ZT6B", "ZT6C", "ZT12A","ZT12B", "ZT12C","ZT18A","ZT18B", "ZT18C")

eig_061218anova<-as.data.frame(rbind(cbind(eig_061218[,1], "ZT0", "A"),
                                     cbind(eig_061218[,2],  "ZT0", "B"),
                                     cbind(eig_061218[,3],  "ZT0", "C"),
                                     cbind(eig_061218[,4],  "ZT6", "A"),
                                     cbind(eig_061218[,5], "ZT6", "B"),
                                     cbind(eig_061218[,6], "ZT6", "C"),
                                     cbind(eig_061218[,7], "ZT12", "A"),
                                     cbind(eig_061218[,8], "ZT12", "B"),
                                     cbind(eig_061218[,9], "ZT12", "C"),
                                     cbind(eig_061218[,10], "ZT18", "A"),
                                     cbind(eig_061218[,11], "ZT18", "B"),
                                     cbind(eig_061218[,12], "ZT18", "C")
))
eig_061218anova[,1]<-as.numeric(levels(eig_061218anova[,1]))[eig_061218anova[,1]]
colnames(eig_061218anova)<-c("pca1", "timepoint", "replicate")
# Compute the analysis of variance
res.aov <- aov(pca1 ~ timepoint * replicate, data = eig_061218anova)
# Summary of the analysis
summary(res.aov)
Df Sum Sq Mean Sq F value Pr(>F)    
timepoint       3  4.145  1.3816    3502 <2e-16 ***
  Residuals   17612  6.947  0.0004                   
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Normality assuption: Checked
plot(res.aov, 2)
