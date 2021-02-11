########## #Libraries
library(data.table)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(scales)
library(forcats)
library(ggrepel)
library(gridExtra)
library(grid)
library(extrafont)
library(dplyr)
library(gdata)
library(rtracklayer)
library(gplots)

#Load PCA1 values of compartments NO changes
nochanges_perreplicates<-fread("/inputs_for_scripts/HiC/Compartments/NochangesZT061218REPLICATESsharedbins100kbeigenKR.bed")
nochanges_mergedreplicates<-fread("/inputs_for_scripts/HiC/Compartments/NochangesZT061218sharedbins100kbeigenKR.bed")

colnames(nochanges_perreplicates)<-c("seqnames", "start", "end", "ZT0A", "ZT0B", "ZT0C" , "ZT6A", "ZT6B", "ZT6C", "ZT12A", "ZT12B", "ZT12C", "ZT18A", "ZT18B", "ZT18C")
colnames(nochanges_mergedreplicates)<-c("seqnames", "start", "end", "ZT0", "ZT6", "ZT12", "ZT18")

#Divide by A or B compartments; PCA1<0 -> B
Bcom_nochanges_mergedreplicates<-nochanges_mergedreplicates[nochanges_mergedreplicates$ZT0<0 &nochanges_mergedreplicates$ZT6<0 &nochanges_mergedreplicates$ZT12<0 &nochanges_mergedreplicates$ZT18<0]
Acom_nochanges_mergedreplicates<-nochanges_mergedreplicates[nochanges_mergedreplicates$ZT0>0 &nochanges_mergedreplicates$ZT6>0 &nochanges_mergedreplicates$ZT12>0 &nochanges_mergedreplicates$ZT18>0]

Bcom_nochanges_perreplicates<-nochanges_perreplicates[nochanges_perreplicates$ZT0A<0 &nochanges_perreplicates$ZT0B<0 &nochanges_perreplicates$ZT0C<0 &
                                                        nochanges_perreplicates$ZT6A<0 &nochanges_perreplicates$ZT6B<0 &nochanges_perreplicates$ZT6C<0 &
                                                        nochanges_perreplicates$ZT12A<0 &nochanges_perreplicates$ZT12B<0 &nochanges_perreplicates$ZT12C<0 &
                                                        nochanges_perreplicates$ZT18A<0 &nochanges_perreplicates$ZT18B<0 & nochanges_perreplicates$ZT18C<0]
Acom_nochanges_perreplicates<-nochanges_perreplicates[nochanges_perreplicates$ZT0A>0 &nochanges_perreplicates$ZT0B>0 &nochanges_perreplicates$ZT0C>0 &
                                                        nochanges_perreplicates$ZT6A>0 &nochanges_perreplicates$ZT6B>0 &nochanges_perreplicates$ZT6C>0 &
                                                        nochanges_perreplicates$ZT12A>0 &nochanges_perreplicates$ZT12B>0 &nochanges_perreplicates$ZT12C>0 &
                                                        nochanges_perreplicates$ZT18A>0 &nochanges_perreplicates$ZT18B>0 & nochanges_perreplicates$ZT18C>0]



#summary(Acom_nochanges_mergedreplicates[,-c(1:3)])
#ZT0                 ZT6                 ZT12                ZT18          
#Min.   :4.184e-05   Min.   :1.641e-05   Min.   :0.0000067   Min.   :5.490e-06  
#1st Qu.:2.556e-02   1st Qu.:2.570e-02   1st Qu.:0.0251910   1st Qu.:2.499e-02  
#Median :3.214e-02   Median :3.290e-02   Median :0.0326327   Median :3.261e-02  
#Mean   :2.985e-02   Mean   :3.048e-02   Mean   :0.0301650   Mean   :3.018e-02  
#3rd Qu.:3.689e-02   3rd Qu.:3.775e-02   3rd Qu.:0.0375342   3rd Qu.:3.768e-02  
#Max.   :5.348e-02   Max.   :5.425e-02   Max.   :0.0550385   Max.   :5.453e-02  

#summary(Bcom_nochanges_mergedreplicates[,-c(1:3)])
#ZT0                  ZT6                  ZT12                 ZT18           
#Min.   :-6.231e-02   Min.   :-6.150e-02   Min.   :-6.305e-02   Min.   :-6.242e-02  
#1st Qu.:-3.390e-02   1st Qu.:-3.359e-02   1st Qu.:-3.373e-02   1st Qu.:-3.364e-02  
#Median :-2.607e-02   Median :-2.574e-02   Median :-2.563e-02   Median :-2.561e-02  
#Mean   :-2.421e-02   Mean   :-2.398e-02   Mean   :-2.405e-02   Mean   :-2.402e-02  
#3rd Qu.:-1.459e-02   3rd Qu.:-1.442e-02   3rd Qu.:-1.448e-02   3rd Qu.:-1.451e-02  
#Max.   :-3.200e-07   Max.   :-1.487e-05   Max.   :-7.330e-06   Max.   :-2.970e-06 

#########################-----------------------
########################## Divide compartments according to sd quantiles

Acom_nochanges_mergedreplicatesSD<-(apply(Acom_nochanges_mergedreplicates[,-c(1:3)], 1, sd))
summary(Acom_nochanges_mergedreplicatesSD)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#2.913e-05 5.286e-04 7.591e-04 8.959e-04 1.092e-03 6.001e-03 
Bcom_nochanges_mergedreplicatesSD<-(apply(Bcom_nochanges_mergedreplicates[,-c(1:3)], 1, sd))
summary(Bcom_nochanges_mergedreplicatesSD)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0000237 0.0005400 0.0009122 0.0011402 0.0014939 0.0077577 
#First quantile means low variance between timepoints; 4th quantile means higher changes between timepoint. This does not take inconsideration strong A
sd_merged_replicates<-as.data.frame(cbindX(as.matrix(Acom_nochanges_mergedreplicatesSD),as.matrix(Bcom_nochanges_mergedreplicatesSD)))
colnames(sd_merged_replicates)<-c("Acompartments", "Bcompartments")
svglite::svglite("/Users/andoku01/Dropbox (Cambridge University)/IFC/Figures2_Masami/HiC/Compartments/Boxplot_StandardDeviation_nochangescompartments_mergedreplicates.svg", height = 5, width = 4)
ggplot(melt(sd_merged_replicates),aes(x=variable, y=value, fill=variable) )+geom_boxplot(notch=T)+ylab("Standard Deviation of PCA1 in all timepoints")+theme_bw()+xlab("")+
  theme(axis.text.x = element_text(vjust = 1, size = 12, color="black"),
        axis.text.y = element_text(size =12, color="black"),
        legend.position = "none",
         panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )
dev.off()


Acom_nochanges_mergedreplicates$SDbetweentimepoints<-NA
Acom_nochanges_mergedreplicates$SDbetweentimepoints[Acom_nochanges_mergedreplicatesSD<=summary(Acom_nochanges_mergedreplicatesSD)[2]]<-"1stquantile"
Acom_nochanges_mergedreplicates$SDbetweentimepoints[Acom_nochanges_mergedreplicatesSD>summary(Acom_nochanges_mergedreplicatesSD)[2] & Acom_nochanges_mergedreplicatesSD<=summary(Acom_nochanges_mergedreplicatesSD)[3]]<-"2ndquantile"
Acom_nochanges_mergedreplicates$SDbetweentimepoints[Acom_nochanges_mergedreplicatesSD>summary(Acom_nochanges_mergedreplicatesSD)[3]& Acom_nochanges_mergedreplicatesSD<=summary(Acom_nochanges_mergedreplicatesSD)[5]]<-"3rdquantile"
Acom_nochanges_mergedreplicates$SDbetweentimepoints[Acom_nochanges_mergedreplicatesSD>summary(Acom_nochanges_mergedreplicatesSD)[5]]<-"4thquantile"
Acom_nochanges_mergedreplicates$SDbetweentimepoints<-factor(Acom_nochanges_mergedreplicates$SDbetweentimepoints, levels=c("1stquantile","2ndquantile","3rdquantile","4thquantile" ))

Bcom_nochanges_mergedreplicates$SDbetweentimepoints<-NA
Bcom_nochanges_mergedreplicates$SDbetweentimepoints[Bcom_nochanges_mergedreplicatesSD<=summary(Bcom_nochanges_mergedreplicatesSD)[2]]<-"1stquantile"
Bcom_nochanges_mergedreplicates$SDbetweentimepoints[Bcom_nochanges_mergedreplicatesSD>summary(Bcom_nochanges_mergedreplicatesSD)[2] & Bcom_nochanges_mergedreplicatesSD<=summary(Bcom_nochanges_mergedreplicatesSD)[3]]<-"2ndquantile"
Bcom_nochanges_mergedreplicates$SDbetweentimepoints[Bcom_nochanges_mergedreplicatesSD>summary(Bcom_nochanges_mergedreplicatesSD)[3]& Bcom_nochanges_mergedreplicatesSD<=summary(Bcom_nochanges_mergedreplicatesSD)[5]]<-"3rdquantile"
Bcom_nochanges_mergedreplicates$SDbetweentimepoints[Bcom_nochanges_mergedreplicatesSD>summary(Bcom_nochanges_mergedreplicatesSD)[5]]<-"4thquantile"

#quant1<-sapply(1:4, function(timepoint){
#  cutoff<-as.numeric(gsub(pattern = "1st Qu.:", replacement = "", x = summary(Acom_nochanges_mergedreplicates[,-c(1:3)])[2,timepoint]));
#  
#  Acom_nochanges_mergedreplicates[Acom_nochanges_mergedreplicates[,c("ZT0", "ZT6", "ZT12", "ZT18")[timepoint]]<=cutoff[timepoint], ]
#       })

#########################-----------------------
#########################Boxplots
#svglite::svglite("/Figures/HiC/Compartments/Boxplot_PCA1_dividedbySD_Acomnochangescompartments_mergedreplicates_2.svg", height = 5, width = 5)
#ggplot(melt(Acom_nochanges_mergedreplicates[,-c(1:3)]), aes(x=SDbetweentimepoints,fill=variable, y=value))+geom_boxplot(notch = T)+
#  #facet_wrap(~SDbetweentimepoints)+
#  theme_bw()+
#  theme( panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(size =10, color="black"),
#         axis.text.y = element_text(size =10, color="black"))+ylab("PCA1")+ggtitle("A-Acompartments")+xlab("")
#dev.off()

#svglite::svglite("/Figures/HiC/Compartments/Boxplot_PCA1_dividedbySD_Bcomnochangescompartments_mergedreplicates_2.svg", height = 5, width = 5#)
#ggplot(melt(Bcom_nochanges_mergedreplicates[,-c(1:3)]), aes(x=SDbetweentimepoints,fill=variable, y=value))+geom_boxplot(notch = T)+
#  #facet_wrap(~SDbetweentimepoints)+
#  theme_bw()+
#  theme( panel.grid.major = element_blank(),
#         axis.text.x = element_text(size =10, color="black"),
#         axis.text.y = element_text(size =10, color="black"))+ylab("PCA1")+ggtitle("B-Bcompartments")+xlab("")
#dev.off()


#########################Geom_lines
#Acom_nochanges_mergedreplicates$bin<-paste(Acom_nochanges_mergedreplicates$seqnames, Acom_nochanges_mergedreplicates$start, Acom_nochanges_mergedreplicates$end, sep = "_")
#ggplot(melt(Acom_nochanges_mergedreplicates[,-c(1:3)],id.vars  =  c("bin","SDbetweentimepoints")), aes(x=variable,y=value))+geom_line(aes(group=bin))+ stat_smooth(aes(group = 1)) +facet_wrap(~SDbetweentimepoints)
#ggplot(melt(Acom_nochanges_mergedreplicates[,-c(1:3)],id.vars  =  c("bin","SDbetweentimepoints")), aes(x=variable,y=value))+geom_jitter()+facet_wrap(~SDbetweentimepoints)
#+ stat_smooth(aes(group = 1)) +facet_wrap(~SDbetweentimepoints)

#Acom_nochanges_mergedreplicates$bin<-NULL
#########################-----------------------
#####################Export files
write.table(x = Acom_nochanges_mergedreplicates[Acom_nochanges_mergedreplicates$SDbetweentimepoints=="1stquantile"], file = "/inputs_for_scripts/HiC/Compartments/Acomp_nochanges_mergedreplicates_1stquant.bed", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = Acom_nochanges_mergedreplicates[Acom_nochanges_mergedreplicates$SDbetweentimepoints=="2ndquantile"], file = "/inputs_for_scripts/HiC/Compartments/Acomp_nochanges_mergedreplicates_2ndquant.bed", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = Acom_nochanges_mergedreplicates[Acom_nochanges_mergedreplicates$SDbetweentimepoints=="3rdquantile"], file = "/inputs_for_scripts/HiC/Compartments/Acomp_nochanges_mergedreplicates_3rdquant.bed", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = Acom_nochanges_mergedreplicates[Acom_nochanges_mergedreplicates$SDbetweentimepoints=="4thquantile"], file = "/inputs_for_scripts/HiC/Compartments/Acomp_nochanges_mergedreplicates_4thquant.bed", quote = F, sep = "\t", row.names = F, col.names = T)


write.table(x = Bcom_nochanges_mergedreplicates[Bcom_nochanges_mergedreplicates$SDbetweentimepoints=="1stquantile"], file = "/inputs_for_scripts/HiC/Compartments/Bcomp_nochanges_mergedreplicates_1stquant.bed", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = Bcom_nochanges_mergedreplicates[Bcom_nochanges_mergedreplicates$SDbetweentimepoints=="2ndquantile"], file = "/inputs_for_scripts/HiC/Compartments/Bcomp_nochanges_mergedreplicates_2ndquant.bed", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = Bcom_nochanges_mergedreplicates[Bcom_nochanges_mergedreplicates$SDbetweentimepoints=="3rdquantile"], file = "/inputs_for_scripts/HiC/Compartments/Bcomp_nochanges_mergedreplicates_3rdquant.bed", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = Bcom_nochanges_mergedreplicates[Bcom_nochanges_mergedreplicates$SDbetweentimepoints=="4thquantile"], file = "/inputs_for_scripts/HiC/Compartments/Bcomp_nochanges_mergedreplicates_4thquant.bed", quote = F, sep = "\t", row.names = F, col.names = T)

#########################Heatmap A compartments
#All quantiles together
#tempo<-Acom_nochanges_mergedreplicates[order(Acom_nochanges_mergedreplicates$SDbetweentimepoints),]
#df<-tempo[,-c(1:3,8)]
### cluster genes
#clusters <- kmeans(dist(df),2)
#order_genes = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")#[sort(clusters$cluster,index.return=TRUE)$ix]

#df$Bins = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")
#meltdf <- melt(df)
#meltdf$Bins = factor(meltdf$Bins,levels=order_genes)
#df$Bins = NULL

#fixedvalue = max(abs(df))
#ggplot(data = meltdf, aes(variable, Bins, fill = value))+
#  geom_tile(size=0)+ 
 # scale_fill_gradient(low = "white", high = "red", na.value = NA,
                      limits=c(0,fixedvalue),name="PCA1")+
  #scale_fill_gradientn(colors=c("#1D06AF","#0614F0","#0658F0","#0689F0","#06BBF0","white","#EDAC09","#F29208","#F26408","#FC2503","#AB0707"),limits=c(-fixedvalue,fixedvalue),values=rescale(c(-fixedvalue,qnorm(0.10,0,fixedvalue/3),qnorm(0.20,0,fixedvalue/3),qnorm(0.30,0,fixedvalue/3),qnorm(0.40,0,fixedvalue/3),qnorm(0.60,0,fixedvalue/3),qnorm(0.70,0,fixedvalue/3),qnorm(0.80,0,fixedvalue/3),qnorm(0.90,0,fixedvalue/3),fixedvalue),to=c(0,1)),name="logFC") +
#  theme_minimal()+ 
#  scale_y_discrete(expand=c(0,0))+
#  scale_x_discrete(expand=c(0,0))+
#  theme(axis.text.x = element_text(vjust = 1, size = 8, color="black"),
#        #axis.text.y = element_text(size =8, color="black")
#        axis.text.y=element_blank(),axis.ticks.y=element_blank()
#  )+ggtitle("A-Acompartment 4th quantile")

#4thQUANT
tempo<-Acom_nochanges_mergedreplicates[Acom_nochanges_mergedreplicates$SDbetweentimepoints=="4thquantile",]
df<-tempo[,-c(1:3,8)]
### cluster genes
clusters <- kmeans(dist(df),2)
order_genes = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")[sort(clusters$cluster,index.return=TRUE)$ix]

df$Bins = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")
meltdf <- melt(df)
meltdf$Bins = factor(meltdf$Bins,levels=order_genes)
df$Bins = NULL

fixedvalue = max(abs(df))
p4<-ggplot(data = meltdf, aes(variable, Bins, fill = value))+
  geom_tile(size=0)+ 
  scale_fill_gradient(low = "white", high = "red", na.value = NA,
                      limits=c(0,fixedvalue),name="PCA1")+
  #scale_fill_gradientn(colors=c("#1D06AF","#0614F0","#0658F0","#0689F0","#06BBF0","white","#EDAC09","#F29208","#F26408","#FC2503","#AB0707"),limits=c(-fixedvalue,fixedvalue),values=rescale(c(-fixedvalue,qnorm(0.10,0,fixedvalue/3),qnorm(0.20,0,fixedvalue/3),qnorm(0.30,0,fixedvalue/3),qnorm(0.40,0,fixedvalue/3),qnorm(0.60,0,fixedvalue/3),qnorm(0.70,0,fixedvalue/3),qnorm(0.80,0,fixedvalue/3),qnorm(0.90,0,fixedvalue/3),fixedvalue),to=c(0,1)),name="logFC") +
  theme_minimal()+ 
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  theme(axis.text.x = element_text(vjust = 1, size = 8, color="black"),
        #axis.text.y = element_text(size =8, color="black")
        axis.text.y=element_blank(),axis.ticks.y=element_blank()
  )+ggtitle("A-Acompartment 4th quantile")


#3rd quant
tempo<-Acom_nochanges_mergedreplicates[Acom_nochanges_mergedreplicates$SDbetweentimepoints=="3rdquantile",]
df<-tempo[,-c(1:3,8)]

### cluster genes
clusters <- kmeans(dist(df),2)
order_genes = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")[sort(clusters$cluster,index.return=TRUE)$ix]

df$Bins = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")
meltdf <- melt(df)
meltdf$Bins = factor(meltdf$Bins,levels=order_genes)
df$Bins = NULL


p3<-ggplot(data = meltdf, aes(variable, Bins, fill = value))+
  geom_tile(size=0)+ 
  scale_fill_gradient(low = "white", high = "red", na.value = NA,
                      limits=c(0,fixedvalue),name="PCA1")+
  #scale_fill_gradientn(colors=c("#1D06AF","#0614F0","#0658F0","#0689F0","#06BBF0","white","#EDAC09","#F29208","#F26408","#FC2503","#AB0707"),limits=c(-fixedvalue,fixedvalue),values=rescale(c(-fixedvalue,qnorm(0.10,0,fixedvalue/3),qnorm(0.20,0,fixedvalue/3),qnorm(0.30,0,fixedvalue/3),qnorm(0.40,0,fixedvalue/3),qnorm(0.60,0,fixedvalue/3),qnorm(0.70,0,fixedvalue/3),qnorm(0.80,0,fixedvalue/3),qnorm(0.90,0,fixedvalue/3),fixedvalue),to=c(0,1)),name="logFC") +
  theme_minimal()+ 
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  theme(axis.text.x = element_text(vjust = 1, size = 8, color="black"),
        #axis.text.y = element_text(size =8, color="black")
        axis.text.y=element_blank(),axis.ticks.y=element_blank()
  )+ggtitle("A-Acompartment 3rd quantile")


#2nd quant
tempo<-Acom_nochanges_mergedreplicates[Acom_nochanges_mergedreplicates$SDbetweentimepoints=="2ndquantile",]
df<-tempo[,-c(1:3,8)]

### cluster genes
clusters <- kmeans(dist(df),2)
order_genes = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")[sort(clusters$cluster,index.return=TRUE)$ix]

df$Bins = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")
meltdf <- melt(df)
meltdf$Bins = factor(meltdf$Bins,levels=order_genes)
df$Bins = NULL


p2<-ggplot(data = meltdf, aes(variable, Bins, fill = value))+
  geom_tile(size=0)+ 
  scale_fill_gradient(low = "white", high = "red", na.value = NA,
                      limits=c(0,fixedvalue),name="PCA1")+
  #scale_fill_gradientn(colors=c("#1D06AF","#0614F0","#0658F0","#0689F0","#06BBF0","white","#EDAC09","#F29208","#F26408","#FC2503","#AB0707"),limits=c(-fixedvalue,fixedvalue),values=rescale(c(-fixedvalue,qnorm(0.10,0,fixedvalue/3),qnorm(0.20,0,fixedvalue/3),qnorm(0.30,0,fixedvalue/3),qnorm(0.40,0,fixedvalue/3),qnorm(0.60,0,fixedvalue/3),qnorm(0.70,0,fixedvalue/3),qnorm(0.80,0,fixedvalue/3),qnorm(0.90,0,fixedvalue/3),fixedvalue),to=c(0,1)),name="logFC") +
  theme_minimal()+ 
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  theme(axis.text.x = element_text(vjust = 1, size = 8, color="black"),
        #axis.text.y = element_text(size =8, color="black")
        axis.text.y=element_blank(),axis.ticks.y=element_blank()
  )+ggtitle("A-Acompartment 2nd quantile")

#1st quant
tempo<-Acom_nochanges_mergedreplicates[Acom_nochanges_mergedreplicates$SDbetweentimepoints=="1stquantile",]
df<-tempo[,-c(1:3,8)]

### cluster genes
clusters <- kmeans(dist(df),2)
order_genes = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")[sort(clusters$cluster,index.return=TRUE)$ix]

df$Bins = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")
meltdf <- melt(df)
meltdf$Bins = factor(meltdf$Bins,levels=order_genes)
df$Bins = NULL


p1<-ggplot(data = meltdf, aes(variable, Bins, fill = value))+
  geom_tile(size=0)+ 
  scale_fill_gradient(low = "white", high = "red", na.value = NA,
                      limits=c(0,fixedvalue),name="PCA1")+
  #scale_fill_gradientn(colors=c("#1D06AF","#0614F0","#0658F0","#0689F0","#06BBF0","white","#EDAC09","#F29208","#F26408","#FC2503","#AB0707"),limits=c(-fixedvalue,fixedvalue),values=rescale(c(-fixedvalue,qnorm(0.10,0,fixedvalue/3),qnorm(0.20,0,fixedvalue/3),qnorm(0.30,0,fixedvalue/3),qnorm(0.40,0,fixedvalue/3),qnorm(0.60,0,fixedvalue/3),qnorm(0.70,0,fixedvalue/3),qnorm(0.80,0,fixedvalue/3),qnorm(0.90,0,fixedvalue/3),fixedvalue),to=c(0,1)),name="logFC") +
  theme_minimal()+ 
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  theme(axis.text.x = element_text(vjust = 1, size = 8, color="black"),
        #axis.text.y = element_text(size =8, color="black")
        axis.text.y=element_blank(),axis.ticks.y=element_blank()
  )+ggtitle("A-Acompartment 1st quantile")

pdf("/Figures/HiC/Compartments/Heatmap_quantiles_Acom.pdf",height=4,width=12)
grid.arrange(p1,p2,p3,p4 , ncol=4)
dev.off()
#########################-----------------------
#########################B compartments
#4thQUANT
tempo<-Bcom_nochanges_mergedreplicates[Bcom_nochanges_mergedreplicates$SDbetweentimepoints=="4thquantile",]
df<-tempo[,-c(1:3,8)]
### cluster genes
clusters <- kmeans(dist(df),2)
order_genes = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")[sort(clusters$cluster,index.return=TRUE)$ix]

df$Bins = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")
meltdf <- melt(df)
meltdf$Bins = factor(meltdf$Bins,levels=order_genes)
df$Bins = NULL

fixedvalue = unlist(df[order(df)][1,1])

p4<-ggplot(data = meltdf, aes(variable, Bins, fill = value))+
  geom_tile(size=0)+ 
  scale_fill_gradient(low = "blue", high = "white", na.value = NA,
                      limits=c(fixedvalue,0),
                      name="PCA1")+
  #scale_fill_gradientn(colors=c("#1D06AF","#0614F0","#0658F0","#0689F0","#06BBF0","white","#EDAC09","#F29208","#F26408","#FC2503","#AB0707"),limits=c(-fixedvalue,fixedvalue),values=rescale(c(-fixedvalue,qnorm(0.10,0,fixedvalue/3),qnorm(0.20,0,fixedvalue/3),qnorm(0.30,0,fixedvalue/3),qnorm(0.40,0,fixedvalue/3),qnorm(0.60,0,fixedvalue/3),qnorm(0.70,0,fixedvalue/3),qnorm(0.80,0,fixedvalue/3),qnorm(0.90,0,fixedvalue/3),fixedvalue),to=c(0,1)),name="logFC") +
  theme_minimal()+ 
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  theme(axis.text.x = element_text(vjust = 1, size = 8, color="black"),
        #axis.text.y = element_text(size =8, color="black")
        axis.text.y=element_blank(),axis.ticks.y=element_blank()
  )+ggtitle("B-Bcompartment 4th quantile")


#3rd quant
tempo<-Bcom_nochanges_mergedreplicates[Bcom_nochanges_mergedreplicates$SDbetweentimepoints=="3rdquantile",]
df<-tempo[,-c(1:3,8)]

### cluster genes
clusters <- kmeans(dist(df),2)
order_genes = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")[sort(clusters$cluster,index.return=TRUE)$ix]

df$Bins = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")
meltdf <- melt(df)
meltdf$Bins = factor(meltdf$Bins,levels=order_genes)
df$Bins = NULL


p3<-ggplot(data = meltdf, aes(variable, Bins, fill = value))+
  geom_tile(size=0)+ 
  scale_fill_gradient(low = "blue", high = "white", na.value = NA,
                      limits=c(fixedvalue,0),name="PCA1")+
  #scale_fill_gradientn(colors=c("#1D06AF","#0614F0","#0658F0","#0689F0","#06BBF0","white","#EDAC09","#F29208","#F26408","#FC2503","#AB0707"),limits=c(-fixedvalue,fixedvalue),values=rescale(c(-fixedvalue,qnorm(0.10,0,fixedvalue/3),qnorm(0.20,0,fixedvalue/3),qnorm(0.30,0,fixedvalue/3),qnorm(0.40,0,fixedvalue/3),qnorm(0.60,0,fixedvalue/3),qnorm(0.70,0,fixedvalue/3),qnorm(0.80,0,fixedvalue/3),qnorm(0.90,0,fixedvalue/3),fixedvalue),to=c(0,1)),name="logFC") +
  theme_minimal()+ 
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  theme(axis.text.x = element_text(vjust = 1, size = 8, color="black"),
        #axis.text.y = element_text(size =8, color="black")
        axis.text.y=element_blank(),axis.ticks.y=element_blank()
  )+ggtitle("B-Bcompartment 3rd quantile")


#2nd quant
tempo<-Bcom_nochanges_mergedreplicates[Bcom_nochanges_mergedreplicates$SDbetweentimepoints=="2ndquantile",]
df<-tempo[,-c(1:3,8)]

### cluster genes
clusters <- kmeans(dist(df),2)
order_genes = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")[sort(clusters$cluster,index.return=TRUE)$ix]

df$Bins = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")
meltdf <- melt(df)
meltdf$Bins = factor(meltdf$Bins,levels=order_genes)
df$Bins = NULL


p2<-ggplot(data = meltdf, aes(variable, Bins, fill = value))+
  geom_tile(size=0)+ 
  scale_fill_gradient(low = "blue", high = "white", na.value = NA,
                      limits=c(fixedvalue,0),name="PCA1")+
  #scale_fill_gradientn(colors=c("#1D06AF","#0614F0","#0658F0","#0689F0","#06BBF0","white","#EDAC09","#F29208","#F26408","#FC2503","#AB0707"),limits=c(-fixedvalue,fixedvalue),values=rescale(c(-fixedvalue,qnorm(0.10,0,fixedvalue/3),qnorm(0.20,0,fixedvalue/3),qnorm(0.30,0,fixedvalue/3),qnorm(0.40,0,fixedvalue/3),qnorm(0.60,0,fixedvalue/3),qnorm(0.70,0,fixedvalue/3),qnorm(0.80,0,fixedvalue/3),qnorm(0.90,0,fixedvalue/3),fixedvalue),to=c(0,1)),name="logFC") +
  theme_minimal()+ 
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  theme(axis.text.x = element_text(vjust = 1, size = 8, color="black"),
        #axis.text.y = element_text(size =8, color="black")
        axis.text.y=element_blank(),axis.ticks.y=element_blank()
  )+ggtitle("B-Bcompartment 2nd quantile")

#1st quant
tempo<-Bcom_nochanges_mergedreplicates[Bcom_nochanges_mergedreplicates$SDbetweentimepoints=="1stquantile",]
df<-tempo[,-c(1:3,8)]

### cluster genes
clusters <- kmeans(dist(df),2)
order_genes = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")[sort(clusters$cluster,index.return=TRUE)$ix]

df$Bins = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")
meltdf <- melt(df)
meltdf$Bins = factor(meltdf$Bins,levels=order_genes)
df$Bins = NULL


p1<-ggplot(data = meltdf, aes(variable, Bins, fill = value))+
  geom_tile(size=0)+ 
  scale_fill_gradient(low = "blue", high = "white", na.value = NA,
                      limits=c(fixedvalue,0),name="PCA1")+
  #scale_fill_gradientn(colors=c("#1D06AF","#0614F0","#0658F0","#0689F0","#06BBF0","white","#EDAC09","#F29208","#F26408","#FC2503","#AB0707"),limits=c(-fixedvalue,fixedvalue),values=rescale(c(-fixedvalue,qnorm(0.10,0,fixedvalue/3),qnorm(0.20,0,fixedvalue/3),qnorm(0.30,0,fixedvalue/3),qnorm(0.40,0,fixedvalue/3),qnorm(0.60,0,fixedvalue/3),qnorm(0.70,0,fixedvalue/3),qnorm(0.80,0,fixedvalue/3),qnorm(0.90,0,fixedvalue/3),fixedvalue),to=c(0,1)),name="logFC") +
  theme_minimal()+ 
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  theme(axis.text.x = element_text(vjust = 1, size = 8, color="black"),
        #axis.text.y = element_text(size =8, color="black")
        axis.text.y=element_blank(),axis.ticks.y=element_blank()
  )+ggtitle("B-Bcompartment 1st quantile")

pdf("/Figures/HiC/Compartments/Heatmap_quantiles_Bcom.pdf",height=4,width=12)
grid.arrange(p1,p2,p3,p4 , ncol=4)
dev.off()
#########################-----------------------
##Use the absolute change
#Acomp
#ZT0 - ZTX
Acom_nochanges_mergedreplicates$ZT0minusZT6<-abs(Acom_nochanges_mergedreplicates$ZT0- Acom_nochanges_mergedreplicates$ZT6)
Acom_nochanges_mergedreplicates$ZT0minusZT12<-abs(Acom_nochanges_mergedreplicates$ZT0- Acom_nochanges_mergedreplicates$ZT12)
Acom_nochanges_mergedreplicates$ZT0minusZT18<-abs(Acom_nochanges_mergedreplicates$ZT0- Acom_nochanges_mergedreplicates$ZT18)

#ZT6 - ZTX
Acom_nochanges_mergedreplicates$ZT6minusZT12<-abs(Acom_nochanges_mergedreplicates$ZT6- Acom_nochanges_mergedreplicates$ZT12)
Acom_nochanges_mergedreplicates$ZT6minusZT18<-abs(Acom_nochanges_mergedreplicates$ZT6- Acom_nochanges_mergedreplicates$ZT18)

#ZT12- ZTX
Acom_nochanges_mergedreplicates$ZT12minusZT18<-abs(Acom_nochanges_mergedreplicates$ZT12- Acom_nochanges_mergedreplicates$ZT18)

svglite::svglite("/Users/andoku01/Dropbox (Cambridge University)/IFC/Figures2_Masami/HiC/Compartments/Boxplot_absolutedifferenceofPCA1_dividedbySD_Acomnochangescompartments_mergedreplicates.svg", height = 5, width = 10)
ggplot(melt(Acom_nochanges_mergedreplicates[,-c(1:7)]), aes(x=variable ,fill=SDbetweentimepoints, y=value))+geom_boxplot(notch = T)+
  #facet_wrap(~SDbetweentimepoints)+
  theme_bw()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.x = element_text(size =10, color="black"),
         axis.text.y = element_text(size =10, color="black"))+ylab("Absolute difference between timepoints PCA1")+ggtitle("A-Acompartments")+xlab("")
dev.off()
summary(Acom_nochanges_mergedreplicates[Acom_nochanges_mergedreplicates$SDbetweentimepoints=="4thquantile", c("ZT0minusZT6", "ZT0minusZT12", "ZT0minusZT18", "ZT6minusZT12", "ZT6minusZT18", "ZT12minusZT18")])
#ZT0minusZT6         ZT0minusZT12        ZT0minusZT18        ZT6minusZT12        ZT6minusZT18       ZT12minusZT18      
#Min.   :7.400e-08   Min.   :7.050e-06   Min.   :3.800e-08   Min.   :0.0000064   Min.   :4.643e-06   Min.   :7.800e-08  
#1st Qu.:1.148e-03   1st Qu.:1.166e-03   1st Qu.:7.336e-04   1st Qu.:0.0008364   1st Qu.:9.389e-04   1st Qu.:7.251e-04  
#Median :2.201e-03   Median :2.104e-03   Median :1.603e-03   Median :0.0017466   Median :1.976e-03   Median :1.540e-03  
#Mean   :2.278e-03   Mean   :2.263e-03   Mean   :1.748e-03   Mean   :0.0019052   Mean   :2.133e-03   Mean   :1.812e-03  
#3rd Qu.:3.030e-03   3rd Qu.:2.902e-03   3rd Qu.:2.506e-03   3rd Qu.:0.0026897   3rd Qu.:2.973e-03   3rd Qu.:2.575e-03  
#Max.   :1.073e-02   Max.   :1.278e-02   Max.   :9.325e-03   Max.   :0.0140363   Max.   :1.080e-02   Max.   :1.062e-02 

#Use above the 3rd quantile of the absolute difference of the 4th quantile of the SD for each timepoint comparison
#1078
#Acom_nochanges_mergedreplicates<-Acom_nochanges_mergedreplicates[Acom_nochanges_mergedreplicates$ZT0minusZT6>=summary(Acom_nochanges_mergedreplicates$ZT0minusZT6[Acom_nochanges_mergedreplicates$SDbetweentimepoints=="4thquantile"])[5]| Acom_nochanges_mergedreplicates$ZT0minusZT12>=summary(Acom_nochanges_mergedreplicates$ZT0minusZT12[Acom_nochanges_mergedreplicates$SDbetweentimepoints=="4thquantile"])[5]| Acom_nochanges_mergedreplicates$ZT0minusZT18>=summary(Acom_nochanges_mergedreplicates$ZT0minusZT18[Acom_nochanges_mergedreplicates$SDbetweentimepoints=="4thquantile"])[5]| Acom_nochanges_mergedreplicates$ZT6minusZT12>=summary(Acom_nochanges_mergedreplicates$ZT6minusZT12[Acom_nochanges_mergedreplicates$SDbetweentimepoints=="4thquantile"])[5]| Acom_nochanges_mergedreplicates$ZT6minusZT18>=summary(Acom_nochanges_mergedreplicates$ZT6minusZT18[Acom_nochanges_mergedreplicates$SDbetweentimepoints=="4thquantile"])[5]| Acom_nochanges_mergedreplicates$ZT12minusZT18>=summary(Acom_nochanges_mergedreplicates$ZT12minusZT18[Acom_nochanges_mergedreplicates$SDbetweentimepoints=="4thquantile"])[5]]
Acom_nochanges_mergedreplicates<-Acom_nochanges_mergedreplicates[Acom_nochanges_mergedreplicates$SDbetweentimepoints=="4thquantile"]
#Filter by max PCA by timepoint
maxZT0_Acom_nochanges_mergedreplicates<-Acom_nochanges_mergedreplicates[apply(Acom_nochanges_mergedreplicates[,c(4:7)], 1, which.max)==1,]
maxZT6_Acom_nochanges_mergedreplicates<-Acom_nochanges_mergedreplicates[apply(Acom_nochanges_mergedreplicates[,c(4:7)], 1, which.max)==2,]
maxZT12_Acom_nochanges_mergedreplicates<-Acom_nochanges_mergedreplicates[apply(Acom_nochanges_mergedreplicates[,c(4:7)], 1, which.max)==3,]
maxZT18_Acom_nochanges_mergedreplicates<-Acom_nochanges_mergedreplicates[apply(Acom_nochanges_mergedreplicates[,c(4:7)], 1, which.max)==4,]

#svglite::svglite("/Users/andoku01/Dropbox (Cambridge University)/IFC/Figures2_Masami/HiC/Compartments/Boxplot_differenceofPCA1_dividedbySD_Acomnochangescompartments_mergedreplicates.svg", height = 5, width = 10)
#ggplot(melt(Acom_nochanges_mergedreplicates[,-c(1:7)]), aes(x=variable ,fill=SDbetweentimepoints, y=value))+geom_boxplot(notch = T)+
  #facet_wrap(~SDbetweentimepoints)+
#  theme_bw()+
#  geom_hline(yintercept = 0, linetype="dashed")+
#  theme( panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(size =10, color="black"),
#         axis.text.y = element_text(size =10, color="black"))+ylab("Absolute difference between timepoints PCA1")+ggtitle("A-Acompartments")+xlab("")
#dev.off()

#########################-----------------------
#Heatmaps of bins with highest ZT0
tempo<-maxZT18_Acom_nochanges_mergedreplicates
df<-tempo[,-c(1:3,8:14)]

### cluster genes
clusters <- kmeans(dist(df),2)
order_genes = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")[sort(clusters$cluster,index.return=TRUE)$ix]

df$Bins = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")
meltdf <- melt(df)
meltdf$Bins = factor(meltdf$Bins,levels=order_genes)
df$Bins = NULL


fixedvalue = max(abs(df))
minval = min(abs(df))
p4<-ggplot(data = meltdf, aes(variable, Bins, fill = value))+
  geom_tile(size=0)+ theme_minimal()+
   scale_fill_steps2(low = "#E2E6BD", mid = "#E99A2C", high = "#D33F6A", midpoint = 0.02,limits=c(0.001,0.04),name="PC1"     )+
  #scale_fill_stepsn(colours = terrain.colors(10))+ 
  #scale_fill_gradient(   low = "grey", high = "red", na.value = NA,
  #                    limits=c(minval,fixedvalue),name="PC1")+
  #scale_fill_gradientn(colors=c("#1D06AF","#0614F0","#0658F0","#0689F0","#06BBF0","white","#EDAC09","#F29208","#F26408","#FC2503","#AB0707"),limits=c(-fixedvalue,fixedvalue),values=rescale(c(-fixedvalue,qnorm(0.10,0,fixedvalue/3),qnorm(0.20,0,fixedvalue/3),qnorm(0.30,0,fixedvalue/3),qnorm(0.40,0,fixedvalue/3),qnorm(0.60,0,fixedvalue/3),qnorm(0.70,0,fixedvalue/3),qnorm(0.80,0,fixedvalue/3),qnorm(0.90,0,fixedvalue/3),fixedvalue),to=c(0,1)),name="logFC") +
   
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  theme(axis.text.x = element_text(vjust = 1, size = 8, color="black"),
        #axis.text.y = element_text(size =8, color="black")
        axis.text.y=element_blank(),axis.ticks.y=element_blank()
  )+ggtitle("A-Acompartment 4th quantile:\n Max PC1 ZT18")+xlab("")

svglite::svglite("/Figures/HiC/Compartments/Boxplot_PCA1maxZT18_Acomnochangescompartments_mergedreplicates.svg", height = 5, width = 4)
ggplot(meltdf, aes(x=variable, y=value))+geom_boxplot(notch = T)+
  #facet_wrap(~SDbetweentimepoints)+
  theme_bw()+
  theme( panel.grid.major = element_blank(),
         axis.text.x = element_text(size =10, color="black"),
         axis.text.y = element_text(size =10, color="black"))+ylab("PC1")+ggtitle("A-Acompartments\n Max PC1 ZT18")+xlab("")
dev.off()

pdf("/Figures/HiC/Compartments/Heatmap_4thquantile_Acom_maxPC1pertimepoint.pdf",height=4,width=12)
grid.arrange(p1,p2,p3,p4, ncol=4)
dev.off()
#########################-----------------------
#####################Export files
write.table(x = maxZT0_Acom_nochanges_mergedreplicates[,1:8], file = "/inputs_for_scripts/HiC/Compartments/ZT0maxPC1_Acomp_nochanges_mergedreplicates_4thquantSD.bed", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = maxZT6_Acom_nochanges_mergedreplicates[,1:8], file = "/inputs_for_scripts/HiC/Compartments/ZT6maxPC1_Acomp_nochanges_mergedreplicates_4thquantSD.bed", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x =maxZT12_Acom_nochanges_mergedreplicates[,1:8] , file = "/inputs_for_scripts/HiC/Compartments/ZT12maxPC1_Acomp_nochanges_mergedreplicates_4thquantSD.bed", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x =maxZT18_Acom_nochanges_mergedreplicates[,1:8], file = "/inputs_for_scripts/HiC/Compartments/ZT18maxPC1_Acomp_nochanges_mergedreplicates_4thquantSD.bed", quote = F, sep = "\t", row.names = F, col.names = T)

#########################-----------------------

#Bcomp
#ZT0 - ZTX
Bcom_nochanges_mergedreplicates$ZT0minusZT6<-abs(Bcom_nochanges_mergedreplicates$ZT0- Bcom_nochanges_mergedreplicates$ZT6)
Bcom_nochanges_mergedreplicates$ZT0minusZT12<-abs(Bcom_nochanges_mergedreplicates$ZT0- Bcom_nochanges_mergedreplicates$ZT12)
Bcom_nochanges_mergedreplicates$ZT0minusZT18<-abs(Bcom_nochanges_mergedreplicates$ZT0- Bcom_nochanges_mergedreplicates$ZT18)

#ZT6 - ZTX
Bcom_nochanges_mergedreplicates$ZT6minusZT12<-abs(Bcom_nochanges_mergedreplicates$ZT6- Bcom_nochanges_mergedreplicates$ZT12)
Bcom_nochanges_mergedreplicates$ZT6minusZT18<-abs(Bcom_nochanges_mergedreplicates$ZT6- Bcom_nochanges_mergedreplicates$ZT18)

#ZT12- ZTX
Bcom_nochanges_mergedreplicates$ZT12minusZT18<-abs(Bcom_nochanges_mergedreplicates$ZT12- Bcom_nochanges_mergedreplicates$ZT18)

svglite::svglite("/Users/andoku01/Dropbox (Cambridge University)/IFC/Figures2_Masami/HiC/Compartments/Boxplot_absolutedifferenceofPCA1_dividedbySD_Bcomnochangescompartments_mergedreplicates.svg", height = 5, width = 10)
ggplot(melt(Bcom_nochanges_mergedreplicates[,-c(1:7)]), aes(x=variable ,fill=SDbetweentimepoints, y=value))+geom_boxplot(notch = T)+
  #facet_wrap(~SDbetweentimepoints)+
  theme_bw()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.x = element_text(size =10, color="black"),
         axis.text.y = element_text(size =10, color="black"))+ylab("Absolute difference between timepoints PC1")+ggtitle("B-Bcompartments")+xlab("")
dev.off()
summary(Bcom_nochanges_mergedreplicates[Bcom_nochanges_mergedreplicates$SDbetweentimepoints=="4thquantile", c("ZT0minusZT6", "ZT0minusZT12", "ZT0minusZT18", "ZT6minusZT12", "ZT6minusZT18", "ZT12minusZT18")])
#ZT0minusZT6         ZT0minusZT12        ZT0minusZT18        ZT6minusZT12        ZT6minusZT18       ZT12minusZT18      
#Min.   :4.700e-07   Min.   :2.060e-07   Min.   :1.586e-06   Min.   :1.219e-06   Min.   :1.229e-06   Min.   :1.469e-06  
#1st Qu.:1.651e-03   1st Qu.:1.501e-03   1st Qu.:1.138e-03   1st Qu.:1.048e-03   1st Qu.:1.170e-03   1st Qu.:1.402e-03  
#Median :2.990e-03   Median :2.945e-03   Median :2.353e-03   Median :2.162e-03   Median :2.437e-03   Median :2.629e-03  
#Mean   :3.253e-03   Mean   :3.132e-03   Mean   :2.845e-03   Mean   :2.347e-03   Mean   :2.593e-03   Mean   :2.793e-03  
#3rd Qu.:4.263e-03   3rd Qu.:4.242e-03   3rd Qu.:3.715e-03   3rd Qu.:3.383e-03   3rd Qu.:3.672e-03   3rd Qu.:3.801e-03  
#Max.   :1.667e-02   Max.   :1.836e-02   Max.   :1.797e-02   Max.   :1.396e-02   Max.   :1.122e-02   Max.   :1.279e-02  
#Keep only 4th quant
Bcom_nochanges_mergedreplicates<-Bcom_nochanges_mergedreplicates[Bcom_nochanges_mergedreplicates$SDbetweentimepoints=="4thquantile"]
#Filter by max PCA by timepoint
maxZT0_Bcom_nochanges_mergedreplicates<-Bcom_nochanges_mergedreplicates[apply(abs(Bcom_nochanges_mergedreplicates[,c(4:7)]), 1, which.max)==1,]
maxZT6_Bcom_nochanges_mergedreplicates<-Bcom_nochanges_mergedreplicates[apply(abs(Bcom_nochanges_mergedreplicates[,c(4:7)]), 1, which.max)==2,]
maxZT12_Bcom_nochanges_mergedreplicates<-Bcom_nochanges_mergedreplicates[apply(abs(Bcom_nochanges_mergedreplicates[,c(4:7)]), 1, which.max)==3,]
maxZT18_Bcom_nochanges_mergedreplicates<-Bcom_nochanges_mergedreplicates[apply(abs(Bcom_nochanges_mergedreplicates[,c(4:7)]), 1, which.max)==4,]

#####################Export files

write.table(x = maxZT0_Bcom_nochanges_mergedreplicates[,1:8], file = "/inputs_for_scripts/HiC/Compartments/ZT0maxPC1_Bcomp_nochanges_mergedreplicates_4thquantSD.bed", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = maxZT6_Bcom_nochanges_mergedreplicates[,1:8], file = "/inputs_for_scripts/HiC/Compartments/ZT6maxPC1_Bcomp_nochanges_mergedreplicates_4thquantSD.bed", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x =maxZT12_Bcom_nochanges_mergedreplicates[,1:8] , file = "/inputs_for_scripts/HiC/Compartments/ZT12maxPC1_Bcomp_nochanges_mergedreplicates_4thquantSD.bed", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x =maxZT18_Bcom_nochanges_mergedreplicates[,1:8], file = "/inputs_for_scripts/HiC/Compartments/ZT18maxPC1_Bcomp_nochanges_mergedreplicates_4thquantSD.bed", quote = F, sep = "\t", row.names = F, col.names = T)
#########################-----------------------
#Heatmaps of bins with highest ZT0
tempo<-maxZT12_Bcom_nochanges_mergedreplicates
df<-tempo[,-c(1:3,8:14)]

### cluster genes
clusters <- kmeans(dist(df),2)
order_genes = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")[sort(clusters$cluster,index.return=TRUE)$ix]

df$Bins = paste(tempo$seqnames, tempo$start,tempo$end,sep = "_")
meltdf <- melt(df)
meltdf$Bins = factor(meltdf$Bins,levels=order_genes)
df$Bins = NULL


fixedvalue = max(abs(df))
minval = min(abs(df))
p4<-ggplot(data = meltdf, aes(variable, Bins, fill = value))+
  geom_tile(size=0)+ theme_minimal()+
  #scale_fill_steps2(low = "#E2E6BD", mid = "#E99A2C", high = "#D33F6A", midpoint = 0.02,limits=c(0.001,0.04),name="PC1"     )+
  scale_fill_steps()+
      #midpoint = 0.02,limits=c(0.001,0.04),name="PC1" )+ 
  #scale_fill_gradient(   low = "grey", high = "red", na.value = NA,
  #                    limits=c(minval,fixedvalue),name="PC1")+
  #scale_fill_gradientn(colors=c("#1D06AF","#0614F0","#0658F0","#0689F0","#06BBF0","white","#EDAC09","#F29208","#F26408","#FC2503","#AB0707"),limits=c(-fixedvalue,fixedvalue),values=rescale(c(-fixedvalue,qnorm(0.10,0,fixedvalue/3),qnorm(0.20,0,fixedvalue/3),qnorm(0.30,0,fixedvalue/3),qnorm(0.40,0,fixedvalue/3),qnorm(0.60,0,fixedvalue/3),qnorm(0.70,0,fixedvalue/3),qnorm(0.80,0,fixedvalue/3),qnorm(0.90,0,fixedvalue/3),fixedvalue),to=c(0,1)),name="logFC") +
  
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  theme(axis.text.x = element_text(vjust = 1, size = 8, color="black"),
        #axis.text.y = element_text(size =8, color="black")
        axis.text.y=element_blank(),axis.ticks.y=element_blank()
  )+ggtitle("B-Bcompartment 4th quantile:\n Max PC1 ZT18")+xlab("")

svglite::svglite("/Figures/HiC/Compartments/Boxplot_PCA1maxZT12_Bcomnochangescompartments_mergedreplicates.svg", height = 5, width = 4)
ggplot(meltdf, aes(x=variable, y=value))+geom_boxplot(notch = T)+
  #facet_wrap(~SDbetweentimepoints)+
  theme_bw()+
  theme( panel.grid.major = element_blank(),
         axis.text.x = element_text(size =10, color="black"),
         axis.text.y = element_text(size =10, color="black"))+ylab("PCA1")+ggtitle("B-Bcompartments\n Max PC1 ZT12")+xlab("")
dev.off()
pdf("/Figures/HiC/Compartments/Heatmap_4thquantile_Bcom_maxPC1pertimepoint.pdf",height=4,width=12)
grid.arrange(p1,p2,p3,p4, ncol=4)
dev.off()

