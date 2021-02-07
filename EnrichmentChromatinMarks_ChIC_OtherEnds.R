#####################--------P-CHiC--------#####################
##All analyses were carried out in R version 3.5.1 (2018-07-02)
#####################---Libraries---####################

library(hexbin)
library(GenomicRanges)
library(parallel)
library(gtools)
library(gdata)
library(doParallel)
library(grid)
#library(biomaRt)
library(gtools)
library(ggplot2)
library(gridExtra)
library(reshape2)
#library(cowplot)
library(gplots)
library(RColorBrewer)
library(rtracklayer)
library(gtools)
library(preprocessCore)
library(edgeR)
require(VennDiagram)
library(viridis)
library(plyr)

#####################---Functions---####################
#Function to convert to a GRanges object
columnnames<-c('chr','start','end')
bedfile<-function(feature, columnnames){
  colnames(feature) <- c('chr','start','end')
  bed<-with(feature, GRanges(chr, IRanges(start+1, end)))
  return(bed)
}

#Function that gets the index of the bait, whose OE interacts with feature of interest
retrieveindex_fromchic<-function(bedfile_feature, chic_otherends_bed, chic_bait_bed){
  t_1<-subjectHits(findOverlaps(bedfile_feature, chic_otherends_bed))
  #Only the indexes of the baits that overlap with a feature
  #Retrieve the baits that correspond to the unique indexes of the otherends
  return(chic_bait_bed[t_1])
}

#Creates a list of several files contaning "name" pattern 
createvars<-function(list, name, phases){
  names(list) <- paste(name, phases, sep = "")
list2env(list , envir = .GlobalEnv)
remove(list)}

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


#####################---Enrichments:chromatin marks and coreclock tfs---##########

###############Open datasets

#CHiC all genes
chic<-read.csv("/inputs_for_scripts/All_step2_washU_text_1.txt", header=F, sep=" ")

#eRNAs, histone modifications and methylation
ernas<-read.csv("/inputs_for_scripts/eRNAs_de_novo_oscillating.txt", header=T, sep="\t") 

H3K27ac_BR<-read.csv("/inputs_for_scripts/circmarks/RenLab-H3K27ac-Liver.peaks", header=F, sep="\t")
H3K27ac_BR<-H3K27ac_BR[,1:3]
H3K4me1_BR<-read.csv("/inputs_for_scripts/circmarks/RenLab-H3K4me1-liver-peaks.txt", header=F, sep="\t")
H3K4me1_BR<-H3K4me1_BR[,1:3]
H3K4me3_BR<-read.csv("/inputs_for_scripts/circmarks/RenLab-H3K4me3-liver-peaks.txt", header=F, sep="\t")
H3K4me3_BR<-H3K4me3_BR[,1:3]
CTCF_BR<-read.csv("/inputs_for_scripts/circmarks/RenLab-CTCF-liver-peaks.txt", header=F, sep="\t")
CTCF_BR<-CTCF_BR[,1:3]
H3K27me3_BR<-read.csv("/inputs_for_scripts/circmarks/RenLab-H3K27me3-liver.peaks", header=F, sep="\t")
supen_liver<-read.csv("/inputs_for_scripts/circmarks/superenhancers_Liver.bed", sep="\t", header = F)
supen_liver_bed<-bedfile(supen_liver)

#Add circproms from MFM RNAseq
circpromphases<-read.csv("/inputs_for_scripts/circproms/HindIIIfragments_circadiangenes_MFMRNAseq_phases.bed", header = F, sep="\t")
circpromphasesINTRONS<-read.csv("/inputs_for_scripts/circproms/HindIIIfragments_circadiangenesonlyintrons_MFMRNAseq_phases.bed", header = F, sep="\t")
allINTRONSchic<-read.csv("/inputs_for_scripts/circproms/CHiC_ov_circproms_onlyintrons.bedpe", header = F, sep="\t")
fangcircproms<- read.csv("/inputs_for_scripts/circproms/Updated/HindIIIfragments_circadiangenes_fang_phases.bed", header = F, sep="\t")

#Chipseq circiadian tf
setwd("/ChIP-seq/")
temp<- list.files(pattern="*.txt")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i], header=T, sep="\t" )[,2:4])

#Divide hic in baits and otherends
chic1<-(chic[,1:3])
chic2<-(chic[,4:6])
colnames(chic1) <- c('chr','start','end')
colnames(chic2) <- c('chr','start','end')
#chic with only introns as baits
chicINTRONS1<-(allINTRONSchic[,1:3])
chicINTRONS2<-(allINTRONSchic[,4:6])
colnames(chicINTRONS1) <- c('chr','start','end')
colnames(chicINTRONS2) <- c('chr','start','end')

#-------------------------------------------------------------------
##############Create separate Granges objects for baits and otherends
chic_bed_baits<-with(chic1, GRanges(chr, IRanges(start+1, end)))
chic_bed_otherends<-with(chic2, GRanges(chr, IRanges(start+1, end)))
chic_bedINTRONS_baits<-with(chicINTRONS1, GRanges(chr, IRanges(start+1, end)))
chic_bedINTRONS_otherends<-with(chicINTRONS2, GRanges(chr, IRanges(start+1, end)))

#ernas
ernas_bed<- bedfile(ernas, columnnames)

#neuro2d
#neuro2d_bed<-bedfile(neuro2d_neg, columnnames)
#H3K4me3 int
#h3k4me3_int_bed<-bedfile(H3K4me3_int, columnnames)

#H3K4me3 prom
#h3k4me3_prom_bed<- bedfile(H3K4me3_prom, columnnames)

#enhancers
#enhancersgen_bed<- bedfile(enhancersgen, columnnames)

#H3K27ac_BR
H3K27ac_BR_bed<- bedfile(H3K27ac_BR, columnnames)

#H3K4me1_BR
H3K4me1_BR_bed<- bedfile(H3K4me1_BR, columnnames)

#H3K4me3_BR
H3K4me3_BR_bed<- bedfile(H3K4me3_BR, columnnames)

#CTCF_BR
CTCF_BR_bed<- bedfile(CTCF_BR, columnnames)

#H3K27me3_BR
H3K27me3_BR_bed<- bedfile(H3K27me3_BR, columnnames)

#circproms
circproms_MFM<- bedfile(circpromphases[1:3], columnnames)
circproms_MFM_introns<- bedfile(circpromphasesINTRONS[1:3], columnnames)

#Chipseq cirdian tf
peaks<-ls(pattern = ".txt")
peaks_bedfiles<-sapply(peaks, function(x){
  t<-get(x)
  bedfile(t, columnnames)
})

list2env(peaks_bedfiles, envir = .GlobalEnv)
#-------------------------------------------------------------------
#Count observed overlaps with only other ends
OLernas_otherends<-subjectHits(findOverlaps(ernas_bed, chic_bed_otherends))
OLH3K27ac_BR_bed_otherends<-subjectHits(findOverlaps(H3K27ac_BR_bed, chic_bed_otherends))
OLH3K4me1_BR_bed_otherends<-subjectHits(findOverlaps(H3K4me1_BR_bed, chic_bed_otherends))
OLH3K4me3_BR_bed_otherends<-subjectHits(findOverlaps(H3K4me3_BR_bed, chic_bed_otherends))
OLCTCF_BR_bed_otherends<-subjectHits(findOverlaps(CTCF_BR_bed, chic_bed_otherends))
OLH3K27me3_BR_bed_otherends<-subjectHits(findOverlaps(H3K27me3_BR_bed, chic_bed_otherends))
OLSUPERENH_bed_otherends<-subjectHits(findOverlaps(supen_liver_bed, chic_bed_otherends))
OLcircpromsmfm_bed_otherends<-subjectHits(findOverlaps(circproms_MFM, chic_bed_otherends))
OLcircpromsmfmonlyintrons_bed_otherends<-subjectHits(findOverlaps(circproms_MFM_introns, chic_bed_otherends))
OLdhs_bed_otherends<-subjectHits(findOverlaps(dhs, chic_bed_otherends))


#subset with otherends, only first appearance
length_peaks<-sapply(peaks, function(x){
  t<-get(x)
  l<-length(subsetByOverlaps(t, chic_bed_otherends))
  return(l)
})

l1_otherends<-length(OLernas_otherends)
l5_otherends<-length(OLH3K27ac_BR_bed_otherends)
l6_otherends<-length(OLH3K4me1_BR_bed_otherends)
l7_otherends<-length(OLH3K4me3_BR_bed_otherends)
l8_otherends<-length(OLCTCF_BR_bed_otherends)
l9_otherends<-length(OLSUPERENH_bed_otherends)
l10_otherends<-length(OLH3K27me3_BR_bed_otherends)
l11_otherends<-length(OLcircpromsmfm_bed_otherends)
l12_otherends<-length(OLcircpromsmfmonlyintrons_bed_otherends)
l13_otherends<-length(OLdhs_bed_otherends)

t_otherends<-c("ernas"=l1_otherends, 
               "H3K27ac"=l5_otherends, "H3K4me1"=l6_otherends, "H3K4me3"=l7_otherends, "CTCF"= l8_otherends, "Superenhancers"= l9_otherends, "H3K27me3"=l10_otherends,"DHS"=l13_otherends,"Circproms_MFMRNAseq"=l10_otherends, "CircpromsINTRONS_MFMRNAseq"=l11_otherends)


#-------------------------------------------------------------------
############Overlap with only other ends

#1) Divide circroms into static and circadian #No overlap between the indexes of the baits corresponding to circ or static, but there's 52960 overlaps in the otherends, meaning the circ and static proms share OE
chicOE_mfmcircpcroms<-chic_bed_otherends[queryHits(findOverlaps(query = chic_bed_baits, subject = circproms_MFM, type = "equal"))]
chicOE_mfmstaticpcroms<-chic_bed_otherends[setdiff(1:248967,queryHits(findOverlaps(query = chic_bed_baits, subject = circproms_MFM, type = "equal")))]

#Filter whole chic into one with only introns># onlt 15 OE overlap between circ and static introns
chicOE_mfmcircpcromsINTRONS<-chic_bedINTRONS_otherends[queryHits(findOverlaps(query = chic_bedINTRONS_baits, subject = circproms_MFM_introns))]
chicOE_mfmstaticpcromsINTRONS<-chic_bedINTRONS_otherends[setdiff(1:5634,queryHits(findOverlaps(query = chic_bedINTRONS_baits, subject = circproms_MFM_introns, type = "equal")))]

#2) Get observed number of overlaps:
countobsOV<-function(featureslist, OEsetcircadian,OEsetstatic ){
  t<-lapply(featureslist, function(feat){
    circ<-length(subjectHits(findOverlaps(feat, OEsetcircadian)))
    stat<-length(subjectHits(findOverlaps(feat, OEsetstatic)))
    return(c("circOv_OE"=circ, "statOv_OE"=stat))})
  #names(t)<-names(featureslist)
  final<-data.frame(matrix(unlist(t), ncol=length(featureslist), byrow=F))
  colnames(final)<-names(featureslist)
  rownames(final)<-c("circOv_OE", "statOv_OE")
  return(final)
  
#  return(t)
  
}

listoffeatures<-list(ernas_bed, H3K27ac_BR_bed, H3K4me1_BR_bed,H3K4me3_BR_bed, CTCF_BR_bed, H3K27me3_BR_bed,supen_liver_bed, dhs )
names(listoffeatures)<-c("ernas_bed", "H3K27ac_BR_bed", "H3K4me1_BR_bed","H3K4me3_BR_bed", "CTCF_BR_bed", "H3K27me3_BR_bed","supen_liver_bed" , "DHS" )
#3) Apply function for chromatin marks
allpromsMFM_obscountsOE<-countobsOV(featureslist = listoffeatures,OEsetcircadian = chicOE_mfmcircpcroms, OEsetstatic = chicOE_mfmstaticpcroms )
allINTRONSMFM_obscountsOE<-countobsOV(featureslist = listoffeatures,OEsetcircadian = chicOE_mfmcircpcromsINTRONS, OEsetstatic = chicOE_mfmstaticpcromsINTRONS )

#OLcircpromsmfmonlyintrons_bed_otherends<-subjectHits(findOverlaps(circproms_MFM_introns, chic_bed_otherends))

#4) apply function for TFs
listoffeaturesTFs<-lapply(peaks , get)
names(listoffeaturesTFs)<-peaks
allpromsMFM_obscountsOE_TFs<-countobsOV(featureslist = listoffeaturesTFs,OEsetcircadian = chicOE_mfmcircpcroms, OEsetstatic = chicOE_mfmstaticpcroms )
allINTRONSMFM_obscountsOE_TFs<-countobsOV(featureslist = listoffeaturesTFs,OEsetcircadian = chicOE_mfmcircpcromsINTRONS, OEsetstatic = chicOE_mfmstaticpcromsINTRONS )

#t_otherends<-c("ernas"=l1_otherends, "H3K27ac"=l5_otherends, "H3K4me1"=l6_otherends, "H3K4me3"=l7_otherends, "CTCF"= l8_otherends, "Superenhancers"= l9_otherends,"H3K27me3"=l10_otherends, "Circproms_MFMRNAseq"=l11_otherends, "CircpromsINTRONS_MFMRNAseq"=l12_otherends)
#barplot(t_otherends,cex.names = 0.65, las=1, col = "red", main = "Overlap with Capture-hiC fragments")

#5) Function of the Enrichment analysis
#1) It identifies the distance between fragments, then takes the length of the corresponding chr and calculated the max available coordinate for that pair of interactions
#2) Starts iterating per chromosome
#3) #The expected number is calculated by take 
chrsizes<-read.csv("/Users/andoku01/Dropbox (Cambridge University)/IFC/inputs_for_scripts/chrsizes_mm9.txt", header = F, sep = "\t")
chrsizes$V1[20]<-"chrMT"
#--------------
#---------------
#This is the latest version of enrichment tool like the one used in chicago! Use this

enrichoverlap<-function(feature_bed, rep, chic1, chic2)
{
  Et2<-vector()
  sapply(1:rep, function(x){
    #First part is to build a random ChiC dataset with the same number of fragments 
    #Retrieve random chromosomes
    t1_baits<-factor(sample(chic1[,1], length(chic1[,1]), replace=T))
   # (sapply(strsplit(x = t1_baits, "chr"), function(x){x[2]}))
    #t1_baits<-factor(t1_baits, levels(t1_baits)[c(1, 12:19, 2:11,"X", "Y")])
    #Retrieve random start positions
    #Max pos of chr, chrsizes
    tempo<-unlist(sapply(1:length(levels(t1_baits)), function(x){rep(chrsizes[x,2],table(t1_baits)[x])}))- abs(chic2[,3]-chic1[,2])
    t2_baits_startpos<-sapply(tempo, function(limit){sample(limit, 1, replace = F)})
    #Distance between baits and oe to determine start position of otherend fragment
    #t3_baits_endpos<-unlist(t2_baits_startpos)+(chic1[,3]-chic1[,2])
    t4_otherends_startpos<-unlist(t2_baits_startpos)+ abs(chic1[,3]-chic1[,2])+abs(chic2[,2]-chic1[,3])
    t5_otherends_endpos<-t4_otherends_startpos+abs(chic2[,3]-chic2[,2])
    
    t8_otherends<-as.data.frame(cbind((as.character(t1_baits)), t4_otherends_startpos,t5_otherends_endpos))
    
    colnames(t8_otherends) <- c('chr','start','end')
    enrich_bed_otherends<- sort(GRanges(seqnames = paste("chr",t8_otherends[,1], sep = ""), ranges = IRanges(as.numeric(as.character(t8_otherends[,2])),as.numeric(as.character(t8_otherends[,3])))))
       #Second part is to make the overlap using the genomic ranges function subsetByOverlap
    #USE #Subsetbyoverlaps
    EOL<-subjectHits(findOverlaps(feature_bed, enrich_bed_otherends))
    El1_1<-length(EOL)
    Et2<-c(Et2, El1_1)
  }
  
  )
}
}
#-------------------------------------------
# 2019-02-09 Enrichment Obs in circadian vs expected in static
#Control set of static interactions
#Non circ interactions:chicOE_mfmstaticpcroms
#To compare with circproms- take 19829 interactions from the static set, count overlaps, repeat 100 times
#To compare with introns- take 5193 interactions from the static set, count overlaps, repeat 100 times

#Function to get expected counts
countexpOV<-function(featureslist, OEsetstatic , rep, numbercircproms){
  #final<-data.frame(matrix(vector(), rep, 7,dimnames=list(c(), names(listoffeatures))),
 #                 stringsAsFactors=F)
  counter<-1
  t<-lapply(featureslist, function(feat){
    counts<-vector()
    t1<-sapply(1:rep, function(x){
        randomindexes<-sample(229647, numbercircproms, replace = F)
        counts_tempo<-length(subjectHits(findOverlaps(feat, OEsetstatic[randomindexes])))
        counts<-c(counts, counts_tempo)
    })
    return(t1)
  })
  final<-data.frame(matrix(unlist(t), nrow=rep, byrow=F))
  colnames(final)<-names(featureslist)
  return(final)
  }

#Apply function for histone marks

no_cores <- 3 ## change how many cores to use
cl<-makeCluster(no_cores)
registerDoParallel(cl)
expallproms_staticinterascontrol<-countexpOV(featureslist = listoffeatures, OEsetstatic = chicOE_mfmstaticpcroms, rep = 100, numbercircproms = 19829)
expintrons_staticinterascontrol<-countexpOV(featureslist = listoffeatures, OEsetstatic = chicOE_mfmstaticpcroms, rep = 100, numbercircproms = 5193)


#Apply function for TFs
expallproms_staticinterascontrol_TFs<-countexpOV(featureslist = listoffeaturesTFs, OEsetstatic = chicOE_mfmstaticpcroms, rep = 100, numbercircproms = 19829)
expintrons_staticinterascontrol_TFs<-countexpOV(featureslist = listoffeaturesTFs, OEsetstatic = chicOE_mfmstaticpcroms, rep = 100, numbercircproms = 5193)

stopCluster(cl)

#---------------------
#Run enrichment on random set of interactions
chic1$chr<-gsub(pattern = "chr", replacement = "", chic1$chr)
OELernas1<-enrichoverlap(ernas_bed, 100,  chic1, chic2)
#EOLh3k4me3_int1<-enrichoverlap(chic_bed_otherends, h3k4me3_int_bed, 100, chic2)
#EOLh3k4me3_prom1<-enrichoverlap(chic_bed_otherends, h3k4me3_prom_bed, 100, chic2)
#EOLenhancer1<-enrichoverlap(chic_bed_otherends, enhancersgen_bed, 100, chic2)

EOLH3K27ac_BR_bed<-enrichoverlap(H3K27ac_BR_bed, 100,  chic1, chic2)
EOLH3K4me1_BR_bed<-enrichoverlap( H3K4me1_BR_bed, 100,  chic1, chic2)
EOLH3K4me3_BR_bed<-enrichoverlap(H3K4me3_BR_bed, 100,  chic1, chic2)
EOLCTCF_BR_bed<-enrichoverlap(CTCF_BR_bed, 100,  chic1, chic2)
EOLSUPERENH_bed<-enrichoverlap(supen_liver_bed, 100,  chic1, chic2)
EOLH3K27me3_bed<-enrichoverlap(H3K27me3_BR_bed, 100,  chic1, chic2)
EOLDHS_bed<-enrichoverlap(dhs, 100,  chic1, chic2)
EOLcircpromsmfmrnaseq_bed<-enrichoverlap(circproms_MFM, 100,  chic1, chic2)
EOLcircpromsmfmrnaseqintrons_bed<-enrichoverlap(circproms_MFM_introns, 100,  chic1, chic2)
#enrichment for peaks
enrichment_peaks<-sapply(peaks, function(x){
  t<-get(x)
  l<-(enrichoverlap(t, 100,  chic1, chic2))
  return(l)
})

#---------------------


#library(gmodels)
#library(ggpubr)
ObsExp_histmarksallproms<-allpromsMFM_obscountsOE[1,]/apply(expallproms_staticinterascontrol, 2, mean)
ObsExp_histmarksintrons<-allINTRONSMFM_obscountsOE[1,]/apply(expintrons_staticinterascontrol, 2, mean)

ObsExp_TFsallproms<-allpromsMFM_obscountsOE_TFs[1,]/apply(expallproms_staticinterascontrol_TFs, 2, mean)
ObsExp_TFsintrons<-allINTRONSMFM_obscountsOE_TFs[1,]/apply(expintrons_staticinterascontrol_TFs, 2, mean)


t<-melt(as.matrix(rbind("Allproms"=ObsExp_histmarksallproms, "Introns"=ObsExp_histmarksintrons)))
tt<-melt(as.matrix(rbind("Allproms"=allpromsMFM_obscountsOE[1,]/(apply(expallproms_staticinterascontrol, 2, mean)+apply(expallproms_staticinterascontrol, 2, sd)), 
                          "Introns"=allINTRONSMFM_obscountsOE[1,]/(apply(expintrons_staticinterascontrol, 2, mean)+apply(expintrons_staticinterascontrol, 2, sd)))))
ttt<-melt(as.matrix(rbind("Allproms"=allpromsMFM_obscountsOE[1,]/(apply(expallproms_staticinterascontrol, 2, mean)-apply(expallproms_staticinterascontrol, 2, sd)), 
                         "Introns"=allINTRONSMFM_obscountsOE[1,]/(apply(expintrons_staticinterascontrol, 2, mean)-apply(expintrons_staticinterascontrol, 2, sd)))))


p1<-ggplot(t, aes(x=Var2, y=value, fill=Var1))+geom_bar(stat="identity", position = "dodge")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("Histone marks")+ylab("Obs/Exp")+
  geom_errorbar(aes(ymin=tt$value, ymax=ttt$value), width=.2,position=position_dodge(.9))+
  scale_fill_grey(start=0.8, end=0.2, name="")+ geom_hline(yintercept=1, linetype="dashed", color = "red")
  
  #stat_compare_means( aes(label = ..p.signif..),label.x = 1.5, label.y = 5)
  #stat_compare_means(aes(group = Var1), method = "t.test", label = "p.signif", label.y = 5)

t1<-melt(as.matrix(rbind("Allproms"=ObsExp_TFsallproms, "Introns"=ObsExp_TFsintrons)))

tt1<-melt(as.matrix(rbind("Allproms"=allpromsMFM_obscountsOE_TFs[1,]/(apply(expallproms_staticinterascontrol_TFs, 2, mean)+apply(expallproms_staticinterascontrol_TFs, 2, sd)), 
                         "Introns"=allINTRONSMFM_obscountsOE_TFs[1,]/(apply(expintrons_staticinterascontrol_TFs, 2, mean)+apply(expintrons_staticinterascontrol_TFs, 2, sd)))))
ttt1<-melt(as.matrix(rbind("Allproms"=allpromsMFM_obscountsOE_TFs[1,]/(apply(expallproms_staticinterascontrol_TFs, 2, mean)-apply(expallproms_staticinterascontrol_TFs, 2, sd)), 
                          "Introns"=allINTRONSMFM_obscountsOE_TFs[1,]/(apply(expintrons_staticinterascontrol_TFs, 2, mean)-apply(expintrons_staticinterascontrol_TFs, 2, sd)))))



p2<-ggplot(t1, aes(x=Var2, y=value, fill=Var1))+geom_bar(stat="identity", position = "dodge")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("TFs")+ylab("Obs/Exp")+
  geom_errorbar(aes(ymin=tt1$value, ymax=ttt1$value), width=.2,position=position_dodge(.9))+
  scale_fill_grey(start=0.8, end=0.2,  name="")+ geom_hline(yintercept=1, linetype="dashed", color = "red")

#svglite::svglite("/Enrichment/Enrichment_epigeneticmarks_ObsinCircinteractions_ExpStaticinteractions.svg")
gridExtra::grid.arrange(p1, p2, ncol=1)
sapply(1:8, function(x){t.test(expallproms_staticinterascontrol[,x], mu = allpromsMFM_obscountsOE[1,x])$p.value})

sapply(1:6, function(x){t.test(expallproms_staticinterascontrol_TFs[,x], mu = allpromsMFM_obscountsOE_TFs[1,x])$p.value})
#dev.off()

#Pvals:
#ernas_bed H3K27ac_BR_bed H3K4me1_BR_bed H3K4me3_BR_bed CTCF_BR_bed H3K27me3_BR_bed  supen_liver_bed DHS              
#[1] 5.755666e-129 7.475513e-143 1.927844e-142  2.478375e-96  2.157781e-47 4.045414e-103
#[7] 2.259606e-130 1.835416e-137

#Plot epigenetic marks
Enrich_features<-cbind("ernas"=OELernas1, 
                       "H3K27ac"=EOLH3K27ac_BR_bed, "H3K4me1"=EOLH3K4me1_BR_bed, "H3K4me3"=EOLH3K4me3_BR_bed, "CTCF"=EOLCTCF_BR_bed, "Superenhancers"=EOLSUPERENH_bed, "H3K27me3"= EOLH3K27me3_bed, "DHS"=EOLDHS_bed,"circpromsmfmrnaseq"=EOLcircpromsmfmrnaseq_bed, "circpromsmfmrnaseqINTRONS"=EOLcircpromsmfmrnaseqintrons_bed)
Enrich_features_means<-apply(Enrich_features, 2, mean) #means of the replicates of each feature
Enrich_features_sd<-apply(Enrich_features, 2, sd) #sd of the replicates of each feature
Enrich_features_LCI<-Enrich_features_means-Enrich_features_sd# confidence intervals  of the replicates of each feature
Enrich_features_HCI<-Enrich_features_means+Enrich_features_sd# confidence intervals  of the replicates of each feature

#barplot(Enrich_features_means, cex.names = 0.65, las=1, col = "blue", main = "Overlap with Random-generated Capture-hiC fragments (only other ends)")

#####################PLOT epigenetic marks, error intervals
#barplot(Enrich_features_means, cex.names = 0.65,  las=1, col = "blue", main = "Overlap with Random-generated Capture-hiC fragments (only other ends)")
#plot(Enrich_features_means, type="n", ylim=c(0,max(t_otherends)),  xaxt = "n", xlab="",ylab="Overlaps with OE",  main = "Overlap with Random-generated Capture-hiC fragments (only other ends)")
#points(Enrich_features_means, col = "blue",pch=19,lwd = 2.5,  xaxt = "n", xlab="")
#points(t_otherends,pch=19,lwd = 2.5,xaxt = "n", xlab="")
#axis(1, at=1:5, labels=names(Enrich_features_means))
#arrows(1:5, Enrich_features_means-Enrich_features_sd, 1:5, Enrich_features_means+Enrich_features_sd, length=0.2, angle=90, code=3)

#legend("topright", c("Observed","Expected"), pch=16,col=c("black","blue"), bty = "n",
#       cex = 1.2, y.intersp=0.2)

#Barplots
#Plot real counts+enrichment of epigenetic features
#svglite::svglite("/Enrichment/Enrichment_epigeneticmarks.svg")
tempo<-barplot((t_otherends/Enrich_features_means)[-(9:10)], ylim=c(0,6),ylab="Obs/Exp Overlaps with C-HiC",
              names.arg = c("eRNAs", "H3K27ac", "H3K4me1", "H3K4me3", "CTCF", "Superenh", "H3K27me3","DHS"),
                            #"Circproms\nMFMRNAseq", "CircpromsINTRONS\nMFMRNAseq"), 
              las=2, col="gray30", border = F)
abline(h = 1, col = "black", lty = 3)
arrows(tempo, (t_otherends[-(9:10)]/(Enrich_features_means[-(9:10)]-Enrich_features_sd[-(9:10)])) ,tempo, (t_otherends[-(9:10)]/(Enrich_features_means[-(9:10)]+Enrich_features_sd[-(9:10)])), length=0.2, angle=90, code=3)
dev.off()
sapply(1:8, function(x){t<-t.test(Enrich_features[,x],mu=t_otherends[x]); t$p.value})
#ernas        H3K27ac        H3K4me1        H3K4me3           CTCF Superenhancers       H3K27me3 
#3.441496e-230  8.946165e-249  1.629007e-261  4.319663e-254  2.108128e-252  1.522489e-213  1.383905e-172 
#DHS 
#4.286069e-250 

####Plot chipseq circadian tfs
###############Plot circadian proms
##Plot real counts+enrichment of circadian tf

#Calculate means, sd and intervals
stats_enrichment<-function(Enrich_features){
  Enrich_features_means<-apply(Enrich_features, 2, mean) #means of the replicates of each feature
  Enrich_features_sd<-apply(Enrich_features, 2, sd) #sd of the replicates of each feature
  Enrich_features_LCI<-Enrich_features_means-Enrich_features_sd# confidence intervals  of the replicates of each feature
  Enrich_features_HCI<-Enrich_features_means+Enrich_features_sd# confidence intervals  of the replicates of each feature
  return(c("mean"=Enrich_features_means,"sd"=Enrich_features_sd, "LCI"=Enrich_features_LCI,"HCI"=Enrich_features_HCI ))
}

#Labels
t<-gsub("mean.", "", names(peak_enrichment_means))
t<-gsub("_peaklist.txt", "", t)

#Plot
peak_enrichment_means<-stats_enrichment(enrichment_peaks)[1:7]
peak_enrichment_sd<-stats_enrichment(enrichment_peaks)[8:14]

#svglite::svglite("/Enrichment/Enrichment_circadianTFs.svg")
tempo<-barplot((length_peaks/peak_enrichment_means),ylim=c(0,6), ylab="Obs/Exp Overlaps with C-HiC",
              names.arg = t,
              #"Circproms\nMFMRNAseq", "CircpromsINTRONS\nMFMRNAseq"), 
              las=2, col="gray30", border = F)
abline(h = 1, col = "black", lty = 3)
arrows(tempo, (length_peaks/(peak_enrichment_means-peak_enrichment_sd)) ,tempo, (length_peaks/(peak_enrichment_means+peak_enrichment_sd)), length=0.2, angle=90, code=3)
dev.off()
sapply(1:7, function(x){t<-t.test(enrichment_peaks[,x],mu=length_peaks[x]); t$p.value})
#BMAL1_peaklist.txt CLOCK_peaklist.txt  CRY1_peaklist.txt  CRY2_peaklist.txt NPAS2_peaklist.txt 
#6.001212e-216      7.439541e-213      6.641199e-234      3.622151e-231      3.870992e-205 
#PER1_peaklist.txt  PER2_peaklist.txt 
#3.597387e-216      5.348654e-234
#BMAL1_peaklist.txt: 6.001212e-216 CLOCK_peaklist.txt: 7.439541e-213  CRY1_peaklist.txt:6.641199e-234
#CRY2_peaklist.txt:3.622151e-231 NPAS2_peaklist.txt:3.870992e-205
#PER1_peaklist.txt:3.597387e-216  PER2_peaklist.txt:5.348654e-234
      

##############Plot circadian proms


tempo_peaks<-rbind(length_peaks,peak_enrichment_means) #plot the means of the overlaps of the real and random dataset
barplot(tempo_peaks,cex.names = 1, beside = TRUE, las=1, main = "Circadian TFs:Overlap with Capture-hiC fragments (only other ends)", xlab = "", names.arg = t)
legend("topright", c("Observed","Expected"),
       pch=16 ,col=c("black","gray"), bty = "n",
       cex = 1.2, y.intersp=0.2)


plot(peak_enrichment_means, type="n", ylim=c(0,max(length_peaks)),  xaxt = "n", xlab="",ylab="Overlaps with OE",  main = "Overlap with Random-generated Capture-hiC fragments (only other ends)")
points(peak_enrichment_means, col = "blue",pch=19,lwd = 2.5,  xaxt = "n", xlab="")
points(length_peaks,pch=19,lwd = 2.5,xaxt = "n", xlab="")
axis(1, at=1:7, labels=t)
arrows(1:7, peak_enrichment_means-stats_enrichment(enrichment_peaks)[8:14], 1:7, peak_enrichment_means+stats_enrichment(enrichment_peaks)[8:14], length=0.2, angle=90, code=3)

legend("topright", c("Observed","Expected"), pch=16,col=c("black","blue"), bty = "n",
       cex = 1.2, y.intersp=0.2)

#####################################################################33
###     P-values
##################################################################

#For epigenetic marks
#First calculate zcores, mean and sd from enrichment data. This will be the null distribution 


#pvals
zscores_epigenmarks_exp<-sapply(1:4, function(x){
  sapply(Enrich_features[,x], function(y){
    zscore<-(y-Enrich_features_means[x])/Enrich_features_sd[x]
    #p<-1-pnorm(abs(zscore))
    return(zscore)
  })
})
hist(zscores_epigenmarks_exp[,1], xlab = "zscores", main="Null distribution")# normalized data

#2) Calculate zscore for observed overlaps
zscore_epigenmarks_obs<-sapply(1:4, function(x){(t_otherends[x]-Enrich_features_means[x])/(Enrich_features_sd[x])})

#3) Pvals


pvals_expected<-sapply(1:4, function(x){
  sapply(zscores_epigenmarks_exp[,x], function(y){
    pnorm(abs(y))
  })
})
hist(pvals_expected)
#La distribución muestra que los pvalores se encuentran cercanos al 1, indicando que suficiente evidencia para no rechazar a la hipotesis nulas
length(pvals_expected[pvals_expected[,1]>0.05,1])


pvals_observed<-sapply(zscore_epigenmarks_obs, function(y){
    pnorm(abs(y))
}) 

#TFs
zscore_tfs_obs<-sapply(1:7, function(x){(length_peaks[x]-peak_enrichment_means[x])/(peak_enrichment_sd[x])})
