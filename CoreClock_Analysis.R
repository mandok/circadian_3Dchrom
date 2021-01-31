#####################---Core-clock interactions---##########
#Open coreclock gene positions
coreclock<-read.csv("Core_clock_short_fragments.txt", sep = "\t", header = T)
coreclock_bed<-GRanges(seqnames = paste("chr", coreclock[,1], sep = ""), ranges = IRanges(coreclock[,2], coreclock[,3]), name=coreclock[,4])
#### Open chic merge file

################################  
#### Open chic merge file s#####
################################  
#I use the symmetrical set
chic<-read.csv("inputs_for_scripts/All_step2_washU_text_symm.txt", sep = "\t", header = F)

#Symmetric set is the complete one
#chic_baits<-read.csv("chic_baits_mm9.txt", sep = "\t", header = F)
chic_baits<-chic[,1:3]
#chic_otherends<-read.csv("chic_otherends_mm9.txt", sep = "\t", header = F)
chic_otherends<-chic[,4:6]

#Bedfiles, GRanges Objects
chic_otherends_bed<-GRanges(seqnames= Rle(chic_otherends[,1]), ranges = IRanges(chic_otherends[,2], chic_otherends[,3]))
chic_bait_bed<-GRanges(seqnames= Rle(chic_baits[,1]), ranges = IRanges(chic_baits[,2], chic_baits[,3]))

###Open eRNAs
#Open TOTAL eRNAs 
total_ernas<-read.csv("eRNAs_de_novo_oscillating.txt", sep = "\t", header = T)
total_ernas<-bedfile(total_ernas, columnnames)

#Osc eRNA
osc_ernas<-read.csv("eRNAs_de_novo_oscillating_phases.txt", header = T, sep = "\t")
osc_ernas<-bedfile(osc_ernas, columnnames)

#Retrieve bait index
BI_oscernas<-retrieveindex_fromchic(total_ernas, chic_otherends_bed, chic_bait_bed) #9730 interactions; sym 109244

coreclock_bed<-circpromphases[circpromphases$name%in% coreclock$Feature]
#Count interactions associated with core clock genes
int_coreclock_OE<-sapply(1:17, function(coreclockgene){
  length(findOverlaps(query = coreclock_bed[coreclockgene], subject = chic_bait_bed))
})
#Count interactions associated with core clock genes associated with eRNAs
int_coreclock_OE_witheRNAs<-sapply(1:17, function(coreclockgene){
  length(findOverlaps(query = coreclock_bed[coreclockgene], subject = BI_oscernas))
})


#other circadian genes, not core clock

#set working directory
circpromphases<-read.csv("inputs_for_scripts/circproms/MFM_RNAseq/HindIIIfragments_circadiangenes_MFMRNAseq_phases.bed", header = F, sep="\t")

circpromphasesINTRONS<-read.csv("inputs_for_scripts/circproms/MFM_RNAseq/HindIIIfragments_circadiangenesonlyintrons_MFMRNAseq_phases.bed", header = F, sep="\t")
#Obtain positions from biomart
#circpromphases_pos<-getBM(attributes=c('chromosome_name','start_position','end_position', 'refseq_mrna'), filters ='refseq_mrna', value=circpromphases[,1],  mart = ensembl)
#circpromphases_pos[,1]<-paste("chr", circpromphases_pos[,1], sep = "")
circpromall<-bedfile(circpromphases, columnnames)
circpromallintrons<-bedfile(circpromphasesINTRONS, columnnames)
circprom_notcoreclock<-circpromphases[!(1:1195 %in% queryHits(findOverlaps( query = circpromphases, coreclock_bed)))]
circpromINTRONS_notcoreclock<-circpromallintrons[!(1:271 %in% queryHits(findOverlaps( circpromallintrons, coreclock_bed)))]

#Core clock Genes that are not in the list of mfm
mcols(coreclock_bed[!(1:25 %in% subjectHits(findOverlaps( circpromall, coreclock_bed)))])
####1   Ccrn4l
#2   Arntl2
#3     Rora
#4   Csnk1d
#5    Fbxl3
#6   Csnk1e
#7  Csnk1a1
#8     Rorb
 
#Functions
#all circproms
#Calculates the distribution of expected interactions of 25 noncoreclock genes, iterating at least 100 times
#The function below consider either all the circgenes detected by the MFMRNAseq or only the introns and also whether the circadian proms interact with an osc erna
notcoreclockgenes_interactions<-function(iterations){
  finalcount<-sapply(1:iterations, function(y){
  rand<-sample(1:1178, 17, replace = F)
  int_notcoreclock_OE<-sapply(rand, function(notcoreclockgene){
    l<-length(findOverlaps(query = circprom_notcoreclock[notcoreclockgene], subject = chic_bait_bed))
    })
  return(int_notcoreclock_OE)
  })
  return(finalcount)
}
notcoreclockgenes_interactions_ernas<-function(iterations){
  finalcount<-sapply(1:iterations, function(y){
    rand<-sample(1:1178, 17, replace = F)
    int_notcoreclock_OE<-sapply(rand, function(notcoreclockgene){
      l<-length(findOverlaps(query = circprom_notcoreclock[notcoreclockgene], subject = BI_oscernas))
    })
    return(int_notcoreclock_OE)
  })
  return(finalcount)
}

circpromall

#introns
notcoreclockgenesINTRONS_interactions<-function(iterations){
  finalcount<-sapply(1:iterations, function(y){
    rand<-sample(1:257, 25, replace = F)
    int_notcoreclock_OE<-sapply(rand, function(notcoreclockgene){
      l<-length(findOverlaps(query = circpromINTRONS_notcoreclock[notcoreclockgene], subject = chic_bait_bed))
    })
    return(int_notcoreclock_OE)
  })
  return(finalcount)
}

notcoreclockgenesINTRONS_interactions_ernas<-function(iterations){
  finalcount<-sapply(1:iterations, function(y){
    rand<-sample(1:257, 25, replace = F)
    int_notcoreclock_OE<-sapply(rand, function(notcoreclockgene){
      l<-length(findOverlaps(query = circpromINTRONS_notcoreclock[notcoreclockgene], subject = BI_oscernas))
    })
    return(int_notcoreclock_OE)
  })
  return(finalcount)
}

#Apply functions
#aLL
tempo_all<-notcoreclockgenes_interactions(100)
tempo_ernas<-notcoreclockgenes_interactions_ernas(100)
#All NO SUBSETTING
int_notcoreclock_obs<-as.vector(table(queryHits(findOverlaps(query = circprom_notcoreclock, subject = chic_bait_bed))))
int_notcoreclock_eRNAs_obs<-as.vector(table(queryHits(findOverlaps(query = circprom_notcoreclock, subject = BI_oscernas))))

mean_iter<-apply(tempo_all, 1, mean)
mean_iter_ernas<-apply(tempo_ernas, 1, mean)
#iNTRONS
tempo_alliNTRONS<-notcoreclockgenesINTRONS_interactions(100)
tempo_ernasINTRONS<-notcoreclockgenesINTRONS_interactions_ernas(100)

mean_iterINTRONS<-apply(tempo_alliNTRONS, 1, mean)
mean_iter_ernasINTRONS<-apply(tempo_ernasINTRONS, 1, mean)

#Plot interactions: Only shows the number of interactions made by the set of core clock genes, including those that have interactions with an osceRNA.
par(las=2)
plot(int_coreclock_OE, ylim=c(0, 51), type="l", col="red", ylab="Number of interactions", xlab="Core clock genes", xaxt="n", main="Interactions of Core clock genes according to CHiC")+axis(1, at=1:25, labels=unlist(as.list(mcols(coreclock_bed))))
lines(int_coreclock_OE_witheRNAs)
lines(mean_iter, col="blue")
lines(mean_iter_ernas, col="green")
legend(x=1, y=50, xpd = TRUE, y.intersp=0.7, x.intersp=0.1, bty = "n", c("Core clock genes-All interactions", "Core clock genes-Interactions with an eRNA", "Equivalent non-core clock genes-All interactions", "Equivalent non-core clock genes-Interactions with an eRNA"), fill= c("red", "black", "blue", "green"), cex=0.8)

###Plot interactions of cc vs all chic and only vs ernas
#All circproms
#svglite::svglite("Coreclockinteractions_expusingallcircproms.svg")
boxplot(cbindX(as.data.frame(int_coreclock_OE), as.data.frame(int_notcoreclock_obs), as.data.frame(int_coreclock_OE_witheRNAs), as.data.frame(int_notcoreclock_eRNAs_obs)),outline=F  ,names=c("Coreclock", "No \nCoreclock", "Coreclock-\neRNAs", "No \nCoreclock-\neRNAs" ), las=2, ylab="Interactions", main="All circproms", notch=T)
boxplot(cbindX(as.data.frame(int_coreclock_OE), as.data.frame(mean_iter), as.data.frame(int_coreclock_OE_witheRNAs), as.data.frame(mean_iter_ernas)),outline=F  ,names=c("Coreclock", "No \nCoreclock", "Coreclock-\neRNAs", "No \nCoreclock-\neRNAs" ), las=2, ylab="Interactions", main="All circproms")
#dev.off()
#Intrnos
#svglite::svglite("Coreclockinteractions_expusingintrons.svg")
boxplot(cbind(int_coreclock_OE, mean_iterINTRONS, int_coreclock_OE_witheRNAs, mean_iter_ernasINTRONS),outline=F  ,names=c("Coreclock", "No \nCoreclock", "Coreclock-\neRNAs", "No \nCoreclock-\neRNAs" ), las=2, ylab="Interactions", main="Introns")
#dev.off()
#See distrbutions, normality test shapiro 
#Circproms
shapiro.test(mean_iter) #Normal; p-value = 0.6261
shapiro.test(int_coreclock_OE)#not normal p-value =0.004195 
#Circproms-ernas
shapiro.test(int_coreclock_OE_witheRNAs) #Normal; p-value =  0.0001595
shapiro.test(mean_iter_ernas)#not normal p-value = 0.3345
#Introns
shapiro.test(mean_iterINTRONS) #Normal; p-value = 0.2721

#Circproms-ernas
shapiro.test(mean_iter_ernasINTRONS)#not normal p-value = 0.2626



#Statistical test
wilcox.test(mean_iter, int_coreclock_OE)[[3]]
#Noncoreclock genes vs coreclock genes: all chic-OE interactions
#pval 0.05463609

wilcox.test(mean_iter_ernas, int_coreclock_OE_witheRNAs)[[3]]
##Noncoreclock genes vs coreclock genes: ernas overlap with chic-OE interactions
#pval 0.330835

wilcox.test(int_coreclock_OE, mean_iterINTRONS)[[3]]
#pval  0.02035716


wilcox.test(int_coreclock_OE_witheRNAs, mean_iter_ernasINTRONS)[[3]]
#pval 0.2634396

#Plot results as barplots
#All circproms
interactions<-cbind(mean(mean_iter), mean(int_coreclock_OE))
sd_interactions<-cbind(sd(mean_iter), sd(int_coreclock_OE))
minus<-interactions[1]-(sd_interactions)
minus[2]<-0
plus<-interactions[1]+sd_interactions[1]
#svglite::svglite("coreclock_vs_allcircproms.svg")
b<-barplot(t(interactions),  las=1, main = "All cirproms", ylab="Interactions", xlab = "CHiC interactions", cex.axis = .9, col = c("dodgerblue4", "darkorange"), beside=T, border="white", cex.lab=.9, space=c(0,0), cex.names = .9, ylim=c(0,22))
error.bar(b, c(interactions[1],0), c(sd_interactions[1], 0))
legend("topright", legend = c("Non core-clock genes", "Core clock genes"),  xpd=T, pch=16 ,col=c("dodgerblue4","darkorange"), bty = "n", cex = .9)
#dev.off()
#Introns  
interactions<-cbind(mean(mean_iterINTRONS), mean(int_coreclock_OE))
sd_interactions<-cbind(sd(mean_iterINTRONS), sd(int_coreclock_OE))
minus<-interactions[1]-(sd_interactions)
minus[2]<-0
plus<-interactions[1]+sd_interactions[1]
#svglite::svglite("coreclockinter_vs_introns.svg")
b<-barplot(t(interactions),  las=1, main = "Introns", ylab="Interactions", xlab = "CHiC interactions", cex.axis = .9, col = c("dodgerblue4", "darkorange"), beside=T, border="white", cex.lab=.9, space=c(0,0), cex.names = .9, ylim=c(0,25))
error.bar(b, c(interactions[1],0), c(sd_interactions[1], 0))
legend("topright", legend = c("Non core-clock genes", "Core clock genes"),  xpd=T, pch=16 ,col=c("dodgerblue4","darkorange"), bty = "n", cex = .9)
#dev.off()

############################################################################################
###########################################
# Dynamic interactions

load("inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_diffs_readcount_above150kb2018")
load("inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_diffs_readcount_below150kb2018")
load("inputs_for_scripts/circtablesinter")
merge_uniquedifinter<-function(tablecircadianinteractions, outputadamscript){
  sapply(1:2, function(updown){
    t1<-tablecircadianinteractions[as.numeric(rownames(outputadamscript[[1]][[updown]])),]
    t2<-tablecircadianinteractions[as.numeric(rownames(outputadamscript[[2]][[updown]])),]
    t3<-tablecircadianinteractions[as.numeric(rownames(outputadamscript[[3]][[updown]])),]
    t4<-tablecircadianinteractions[as.numeric(rownames(outputadamscript[[4]][[updown]])),]
    t5<-tablecircadianinteractions[as.numeric(rownames(outputadamscript[[5]][[updown]])),]
    t6<-tablecircadianinteractions[as.numeric(rownames(outputadamscript[[6]][[updown]])),]
    alldiffinter<-rbind(t1, t2, t3, t4, t5, t6)
    return(alldiffinter)
  })
  #alldiffinter<-unique(alldiffinter[order(alldiffinter$V1, alldiffinter$V2, alldiffinter$V3, alldiffinter$V4, alldiffinter$V5, alldiffinter$V6),])
}


#>150kb
allreadsabove_difinter_timepoint<-merge_uniquedifinter(circadian_tableinteractions, diffs_readcount_above150kb)
allreadsabove_difinter_timepoint<-rbind(as.data.frame(allreadsabove_difinter_timepoint[,1]), as.data.frame(allreadsabove_difinter_timepoint[,2]))#transform into one large data frame
allreadsabove_difinter_timepoint<-unique(allreadsabove_difinter_timepoint[order(allreadsabove_difinter_timepoint$V1, allreadsabove_difinter_timepoint$V2, allreadsabove_difinter_timepoint$V3, allreadsabove_difinter_timepoint$V4, allreadsabove_difinter_timepoint$V5, allreadsabove_difinter_timepoint$V6),]) #unique 

#<150kb
allreadsbelow_difinter_timepoint<-merge_uniquedifinter(circadian_tableinteractions, diffs_readcount_below150kb)
allreadsbelow_difinter_timepoint<-rbind(as.data.frame(allreadsbelow_difinter_timepoint[,1]), as.data.frame(allreadsbelow_difinter_timepoint[,2]))#transform into one large data frame
allreadsbelow_difinter_timepoint<-unique(allreadsbelow_difinter_timepoint[order(allreadsbelow_difinter_timepoint$V1, allreadsbelow_difinter_timepoint$V2, allreadsbelow_difinter_timepoint$V3, allreadsbelow_difinter_timepoint$V4, allreadsbelow_difinter_timepoint$V5, allreadsbelow_difinter_timepoint$V6),]) #unique 

#Overlap of 982 rows
#Read nondifferential inters: static; This set includes all the interactions in circadian-table
allNONdiffinter_abovebelow150Kb<-read.csv("inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_allNONdiffinter_abovebelow150kb2020includingnotMFMcps.bedpe", sep="\t", header=F)
#Bedfiles, GRanges Objects

STATICchic_bait_bed<-GRanges(seqnames= Rle(allNONdiffinter_abovebelow150Kb[,1]), ranges = IRanges(allNONdiffinter_abovebelow150Kb[,2], allNONdiffinter_abovebelow150Kb[,3]))
STATICchic_otherends_bed<-GRanges(seqnames= Rle(allNONdiffinter_abovebelow150Kb[,4]), ranges = IRanges(allNONdiffinter_abovebelow150Kb[,5], allNONdiffinter_abovebelow150Kb[,6]))

#t1<-queryHits(findOverlaps(query = STATICchic_bait_bed, c(BELOWgrangesconservedinalltimepoints_B, ABOVEgrangesconservedinalltimepoints_B), type="equal"))
#t2<-queryHits(findOverlaps(query = STATICchic_otherends_bed, c(BELOWgrangesconservedinalltimepoints_OE, ABOVEgrangesconservedinalltimepoints_OE), type="equal"))
#STATICchic_bait_bed<-STATICchic_bait_bed[1:length(STATICchic_bait_bed) %in% intersect(t1, t2)]
######Observed: Diffinters overlapping coreclock function / static interactions

#coreclock_difinters<-function(difinters,coreclock_bed ){
  grangesconservedinalltimepoints_B<-GRanges(seqnames = difinters[,1], ranges = IRanges(difinters[,2], difinters[,3]), ZT0_difint= difinters[,7], ZT6_difint= difinters[,8], ZT12_difint= difinters[,9], ZT18_difint= difinters[,10])
  grangesconservedinalltimepoints_OE<-GRanges(seqnames = difinters[,4], ranges = IRanges(difinters[,5], difinters[,6]),ZT0_difint= difinters[,7], ZT6_difint= difinters[,8], ZT12_difint= difinters[,9], ZT18_difint= difinters[,10] )
  #overlap cc with dynamic inters
  ov_cc_dyna<-length(findOverlaps(grangesconservedinalltimepoints_B, coreclock_bed))
  ov_cc_dyna_index_baits<-grangesconservedinalltimepoints_B[queryHits(findOverlaps(grangesconservedinalltimepoints_B, coreclock_bed))]
  ov_cc_dyna_index_oe<-grangesconservedinalltimepoints_OE[queryHits(findOverlaps(grangesconservedinalltimepoints_B, coreclock_bed))]
  
  #overlap with static inter
  ov_cc_static<-length(findOverlaps(STATICchic_bait_bed, coreclock_bed))
  ov_cc_static_bait<-STATICchic_bait_bed[queryHits(findOverlaps(STATICchic_bait_bed, coreclock_bed))]
  ov_cc_static_oe<-STATICchic_otherends_bed[queryHits(findOverlaps(STATICchic_bait_bed, coreclock_bed))]
  
  #how many cc-dynamicinters are in cc-staticinters
  dynaoverstatic_baitindex<-(queryHits(findOverlaps(ov_cc_static_bait, ov_cc_dyna_index_baits)))
  dynaoverstatic_oeindex<-(queryHits(findOverlaps(ov_cc_static_oe, ov_cc_dyna_index_oe)))
  
    t<-c(length(dynaoverstatic_oeindex), ov_cc_static, ov_cc_dyna)
  names(t)<-c("Dynamicfromstatic", "Coreclock_with_staticinteractions", "Coreclock_with_dynamicinteractions")
  return(t)
}
#obs_cc_above<-coreclock_difinters(allreadsabove_difinter_timepoint, coreclock_bed )#
#obs_cc_above
#Dynamicfromstatic  Coreclock_with_staticinteractions Coreclock_with_dynamicinteractions 
#6                                297                                  6                               372                                  6 
#obs_cc_below<-coreclock_difinters(allreadsbelow_difinter_timepoint, coreclock_bed )
#Dynamicfromstatic  Coreclock_with_staticinteractions Coreclock_with_dynamicinteractions 
#70                                297                                 72 


#####################Separately CC's

mergebaitandoe<-function(tableinters_bait,tableinters_OE,  coreclock){
  t<-tableinters_bait$Freq
  names(t)<-coreclock$name[tableinters_bait$Var1]
  t2<-tableinters_OE$Freq
  names(t2)<-coreclock$name[tableinters_OE$Var1]
  #Merge bait and OE numbers
  t3<-sapply(coreclock$name, function(x){sum(c(t,t2)[names(c(t,t2))%in% x])})
  names(t3)<-coreclock$name
  return(t3)
}

coreclock_separately_difinters<-function(difintersbelow,difintersabove, coreclock_bed){
  BELOWgrangesconservedinalltimepoints_B<-GRanges(seqnames = difintersbelow[,1], ranges = IRanges(difintersbelow[,2], difintersbelow[,3]), ZT0_difint= difintersbelow[,7], ZT6_difint= difintersbelow[,8], ZT12_difint= difintersbelow[,9], ZT18_difint= difintersbelow[,10])
  BELOWgrangesconservedinalltimepoints_OE<-GRanges(seqnames = difintersbelow[,4], ranges = IRanges(difintersbelow[,5], difintersbelow[,6]),ZT0_difint= difintersbelow[,7], ZT6_difint= difintersbelow[,8], ZT12_difint= difintersbelow[,9], ZT18_difint= difintersbelow[,10] )

  ABOVEgrangesconservedinalltimepoints_B<-GRanges(seqnames = difintersabove[,1], ranges = IRanges(difintersabove[,2], difintersabove[,3]), ZT0_difint= difintersabove[,7], ZT6_difint= difintersabove[,8], ZT12_difint= difintersabove[,9], ZT18_difint= difintersabove[,10])
  ABOVEgrangesconservedinalltimepoints_OE<-GRanges(seqnames = difintersabove[,4], ranges = IRanges(difintersabove[,5], difintersabove[,6]),ZT0_difint= difintersabove[,7], ZT6_difint= difintersabove[,8], ZT12_difint= difintersabove[,9], ZT18_difint= difintersabove[,10] )
  
  #Overlaps cc with total inters  
  ov_cc_total<-as.data.frame(table(subjectHits(findOverlaps(query = chic_bait_bed,subject =  coreclock_bed))))
  ov_cc_total$Var1<-as.integer(levels(ov_cc_total$Var1))[ov_cc_total$Var1]
  #ov_cc_total$genename<-centralcc$genename[ov_cc_total$Var1]
  ov_cc_total_1<-ov_cc_total$Freq
  names(ov_cc_total_1)<-coreclock_bed$name[ov_cc_total$Var1]
  
  #overlap
  #overlap cc with dynamic inters
  #ov_cc_dyna<-length(findOverlaps(grangesconservedinalltimepoints_B, coreclock_bed))
  #Overlaps of baits and other ends with core clock, to make it symmetrical
  ov_cc_dynaB<-as.data.frame(table(subjectHits(findOverlaps(query = (c(BELOWgrangesconservedinalltimepoints_B, ABOVEgrangesconservedinalltimepoints_B)),subject =  coreclock_bed))))
  ov_cc_dynaOE<-as.data.frame(table(subjectHits(findOverlaps(query = (c(BELOWgrangesconservedinalltimepoints_OE, ABOVEgrangesconservedinalltimepoints_OE)),subject =  coreclock_bed))))
  ov_cc_dynaB$Var1<-as.integer(levels(ov_cc_dynaB$Var1))[ov_cc_dynaB$Var1]
  #ov_cc_dynaB$genename<-coreclock_bed$genename[ov_cc_dynaB$Var1]
  ov_cc_dynaOE$Var1<-as.integer(levels(ov_cc_dynaOE$Var1))[ov_cc_dynaOE$Var1]
  ov_cc_dyna<-mergebaitandoe(ov_cc_dynaB,ov_cc_dynaOE,  coreclock_bed)
  
  #ov_cc_dynaOE$genename<-coreclock_bed$genename[ov_cc_dynaOE$Var1]
  #ov_cc_dyna_index_baits<-grangesconservedinalltimepoints_B[queryHits(findOverlaps(grangesconservedinalltimepoints_B, coreclock_bed))]
  #ov_cc_dyna_index_oe<-grangesconservedinalltimepoints_OE[queryHits(findOverlaps(grangesconservedinalltimepoints_B, coreclock_bed))]
  #findOverlaps(c(BELOWgrangesconservedinalltimepoints_B, ABOVEgrangesconservedinalltimepoints_B), coreclock_bed)[which(subjectHits(findOverlaps(c(BELOWgrangesconservedinalltimepoints_B, ABOVEgrangesconservedinalltimepoints_B), centralcc))==2)]
  #grangesconservedinalltimepoints_OE[4450:4452]
  #c(BELOWgrangesconservedinalltimepoints_B, ABOVEgrangesconservedinalltimepoints_B)[141:146]
  #c(BELOWgrangesconservedinalltimepoints_OE, ABOVEgrangesconservedinalltimepoints_OE)[141:146]
  
  #overlap with static inter
  #Overlaps of baits and other ends with core clock, to make it symmetrical
  ov_cc_staticB<-as.data.frame(table(subjectHits(findOverlaps(STATICchic_bait_bed, coreclock_bed))))
  ov_cc_staticB$Var1<-as.integer(levels(ov_cc_staticB$Var1))[ov_cc_staticB$Var1]
  #ov_cc_staticB$genename<-centralcc$genename[ov_cc_staticB$Var1]
 
  ov_cc_staticOE<-as.data.frame(table(subjectHits(findOverlaps(STATICchic_otherends_bed, coreclock_bed))))
  ov_cc_staticOE$Var1<-as.integer(levels(ov_cc_staticOE$Var1))[ov_cc_staticOE$Var1]
  #ov_cc_staticOE$genename<-centralcc$genename[ov_cc_staticOE$Var1]
  
  ov_cc_static<-mergebaitandoe(ov_cc_staticB,ov_cc_staticOE,  coreclock_bed)
  
  #ov_cc_static_bait<-chic_bait_bed[queryHits(findOverlaps(chic_bait_bed, coreclock_bed))]
 # ov_cc_static_oe<-chic_otherends_bed[queryHits(findOverlaps(chic_bait_bed, coreclock_bed))]
  
  #how many cc-dynamicinters are in cc-staticinters
  #dynaoverstatic_baitindex<-(queryHits(findOverlaps(ov_cc_static_bait, ov_cc_dyna_index_baits)))
  #dynaoverstatic_baitindex<-as.data.frame(table(subjectHits(findOverlaps(ov_cc_static_bait, ov_cc_dyna_index_baits))))$Freq
  
  #dynaoverstatic_oeindex<-(queryHits(findOverlaps(ov_cc_static_oe, ov_cc_dyna_index_oe)))
  
  #t<-merge(ov_cc_staticB[, c("Freq", "genename")], ov_cc_dynaB[, c("Freq", "genename")], by= "genename")
  t<-list(ov_cc_total_1,ov_cc_dyna,ov_cc_static)
  names(t)<-c("All_inters", "Dynamic", "Static")
  return(t)}
#obs_cc_sep_above<-coreclock_separately_difinters(allreadsabove_difinter_timepoint, coreclock_bed )#
#Run function
#coreclock_bed$genename<-as.character(levels(coreclock_bed$genename))[coreclock_bed$genename]

obs_cc_sep_belowabove<-coreclock_separately_difinters(allreadsbelow_difinter_timepoint,allreadsabove_difinter_timepoint, circpromphases[circpromphases$name%in% coreclock_bed$name] )#
library(plyr)
#Put in a df
listtodf<-function(listovresults){
  t<-as.data.frame(cbind("Names"=names(listovresults$All_inters), "All_inters"=listovresults$All_inters))
  t2<-as.data.frame(cbind("Names"=names(listovresults$Dynamic), "Dynamic"=listovresults$Dynamic))
  t3<-as.data.frame(cbind("Names"=names(listovresults$Static), "Static"=listovresults$Static))
  t$All_inters<-as.integer(levels(t$All_inters))[t$All_inters]
  t2$Dynamic<-as.integer(levels(t2$Dynamic))[t2$Dynamic]
  t3$Static<-as.integer(levels(t3$Static))[t3$Static]
  t$Names<-as.character(levels(t$Names))[t$Names]
  t2$Names<-as.character(levels(t2$Names))[t2$Names]
  t3$Names<-as.character(levels(t3$Names))[t3$Names]
  
  #t4<-as.data.frame(cbind(t[,-2], t2[,-2], t3[,-2] ))
  #row.names(t4)<-t3$Name
  #colnames(t4)<-c("All_inters", "Dynamic", "Static")
  t4<-join_all(dfs = list( t2,t, t3),by = "Names",match = "all" )
  t4<-t4[,c("Names", "All_inters", "Dynamic", "Static")]
  return(t4)
}
listtodf(obs_cc_sep_belowabove)
#t2<-coreclock_separately_difinters(allreadsbelow_difinter_timepoint,allreadsabove_difinter_timepoint,  coreclock_bed )#
#apply(t2[t2$genename %in% c("Arntl", "Clock", "Cry2", "Npas2", "Nr1d1", "Nr1d2", "Per1", "Per2", "Rorc"), -1], 2, median)
#apply(t2[!(t2$genename %in% c("Arntl", "Clock", "Cry2", "Npas2", "Nr1d1", "Nr1d2", "Per1", "Per2", "Rorc")), -1], 2, mean)
#apply(t2[, -1], 2, mean)
#apply(t2[, -1], 2, median)

##### Expected with updated static set
apply(exp_cc_above, 1, median)
apply(exp_cc_below, 1, median)

#Export results as a table:

t1<-rbind(exp_cc_below[-3,],exp_cc_above[-3,])
t1<-t(t1)

colnames(t1)<-c("Dynamicinters_below150kb", "Staticinters_below150kb", "Dynamicinters_above150kb", "Staticinters_above150kb")
t1<-t1[,c("Dynamicinters_below150kb","Dynamicinters_above150kb", "Staticinters_below150kb", "Staticinters_above150kb") ]
apply(t1, 2, median)
write.table(x = t1, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/Expected_interactions_coreclock_100iter.txt", quote = F, sep = "\t", row.names = F, col.names = T)


#Extra: Dina vs static of genes Arntl", "Per1", "Per2", "Clock", "Npas2", "Cry1","Cry2", "Nr1d1", "Nr1d2", "Rorc
t<-(circpromphases[circpromphases$V4%in% c("Arntl", "Per1", "Per2", "Clock", "Npas2", "Cry1","Cry2", "Nr1d1", "Nr1d2", "Rorc"),])
colnames(t)<-c("seqnames", "start", "end", "genename", "transcriptional_phase")
centralcc<-makeGRangesFromDataFrame(t, keep.extra.columns = T)
#2020-07-06 all coreclocks that appeared in MFMRNAseq
t<-circpromphases[circpromphases$V4%in% coreclock_bed$genename,]

#########2020-02-10 Expected: Diffinters overlapping circproms not coreclock  / static interactions + overlap with eRNAs

#expected_coreclock_difinters_oveRNAs<-function(difinters,circprom_notcoreclock, iterations ){
  t<-sapply(1:iterations, function(x){
    grangesconservedinalltimepoints_B<-GRanges(seqnames = difinters[,1], ranges = IRanges(difinters[,2], difinters[,3]), ZT0_difint= difinters[,7], ZT6_difint= difinters[,8], ZT12_difint= difinters[,9], ZT18_difint= difinters[,10])
    grangesconservedinalltimepoints_OE<-GRanges(seqnames = difinters[,4], ranges = IRanges(difinters[,5], difinters[,6]),ZT0_difint= difinters[,7], ZT6_difint= difinters[,8], ZT12_difint= difinters[,9], ZT18_difint= difinters[,10] )
    #overlap cc with dynamic inters
    rd<-sample(1:length(circprom_notcoreclock), 17, replace = F)
    ov_cc_dyna<-length(findOverlaps(grangesconservedinalltimepoints_B, circprom_notcoreclock[rd]))
    ov_cc_dyna_index_baits<-grangesconservedinalltimepoints_B[queryHits(findOverlaps(grangesconservedinalltimepoints_B, circprom_notcoreclock[rd]))]
    ov_cc_dyna_index_oe<-grangesconservedinalltimepoints_OE[queryHits(findOverlaps(grangesconservedinalltimepoints_B, circprom_notcoreclock[rd]))]
    
    #overlap with static inter
    ov_cc_static<-length(findOverlaps(STATICchic_bait_bed, circprom_notcoreclock[rd]))
    ov_cc_static_bait<-STATICchic_bait_bed[queryHits(findOverlaps(STATICchic_bait_bed, circprom_notcoreclock[rd]))]
    ov_cc_static_oe<-STATICchic_otherends_bed[queryHits(findOverlaps(STATICchic_bait_bed, circprom_notcoreclock[rd]))]
    
    #overlap with eRNAs
    ov_cc_static_oe_eRNAs<- ov_cc_static_oe[queryHits(findOverlaps(query =ov_cc_static_oe , subject = total_ernas))]
    
    #how many cc-dynamicinters are in cc-staticinters
    dynaoverstatic_baitindex<-(queryHits(findOverlaps(query = ov_cc_static_bait, subject = ov_cc_dyna_index_baits)))
    dynaoverstatic_oeindex<-(queryHits(findOverlaps(query = ov_cc_static_oe, subject = ov_cc_dyna_index_oe)))
    dynaoverstatic_oeindex_ov_ernas<-(queryHits(findOverlaps(query = ov_cc_static_oe_eRNAs, subject = ov_cc_dyna_index_oe)))
    
    t<-c(length(dynaoverstatic_oeindex_ov_ernas),length(dynaoverstatic_oeindex), ov_cc_static, ov_cc_dyna)
    names(t)<-c("Dynamicfromstatic_oveRNAs", "Dynamicfromstatic", "Coreclock_with_staticinteractions", "Coreclock_with_dynamicinteractions")
    return(t)
    
  })
  return(t)
}
#exp_cc_above_oveRNAs<-expected_coreclock_difinters_oveRNAs(allreadsabove_difinter_timepoint, circprom_notcoreclock, 100)
#exp_cc_below_oveRNAs<-expected_coreclock_difinters_oveRNAs(allreadsbelow_difinter_timepoint, circprom_notcoreclock, 100)


#Export results as a table:

#t11<-rbind(exp_cc_below_oveRNAs[-4,],exp_cc_above_oveRNAs[-4,])
#t11<-t(t11)
#colnames(t11)<-c("Dynamicfromstatic_oveRNAs_below150kb", "Dynamicinters_below150kb", "Staticinters_below150kb","Dynamicfromstatic_oveRNAs_above150kb" ,"Dynamicinters_above150kb", "Staticinters_above150kb")
#t11<-t11[,c("Dynamicfromstatic_oveRNAs_below150kb", "Dynamicinters_below150kb", "Staticinters_below150kb","Dynamicfromstatic_oveRNAs_above150kb" ,"Dynamicinters_above150kb", "Staticinters_above150kb")]
#apply(t11, 2, median)
#write.table(x = t11, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/Expected_interactions_coreclock_100iter_overlapwitheRNAs.txt", quote = F, sep = "\t", row.names = F, col.names = T)

##############2020-07-06 Overlap with eRNAs

osc_ernas<-bedfile(read.csv("eRNAs_de_novo_oscillating_phases.txt", header = T, sep = "\t"),columnnames)
total_ernas<-bedfile(-read.csv("eRNAs_de_novo_oscillating.txt", sep = "\t", header = T), columnnames)

coreclock_separately_difinters_overnas<-function(difintersbelow,difintersabove, coreclock_bed, total_ernas){
  BELOWgrangesconservedinalltimepoints_B<-GRanges(seqnames = difintersbelow[,1], ranges = IRanges(difintersbelow[,2], difintersbelow[,3]), ZT0_difint= difintersbelow[,7], ZT6_difint= difintersbelow[,8], ZT12_difint= difintersbelow[,9], ZT18_difint= difintersbelow[,10])
  BELOWgrangesconservedinalltimepoints_OE<-GRanges(seqnames = difintersbelow[,4], ranges = IRanges(difintersbelow[,5], difintersbelow[,6]),ZT0_difint= difintersbelow[,7], ZT6_difint= difintersbelow[,8], ZT12_difint= difintersbelow[,9], ZT18_difint= difintersbelow[,10] )
  
  ABOVEgrangesconservedinalltimepoints_B<-GRanges(seqnames = difintersabove[,1], ranges = IRanges(difintersabove[,2], difintersabove[,3]), ZT0_difint= difintersabove[,7], ZT6_difint= difintersabove[,8], ZT12_difint= difintersabove[,9], ZT18_difint= difintersabove[,10])
  ABOVEgrangesconservedinalltimepoints_OE<-GRanges(seqnames = difintersabove[,4], ranges = IRanges(difintersabove[,5], difintersabove[,6]),ZT0_difint= difintersabove[,7], ZT6_difint= difintersabove[,8], ZT12_difint= difintersabove[,9], ZT18_difint= difintersabove[,10] )
  
  #Overlaps cc with total inters  
  #Take chic that overlaps eRNAs
  ov_cc_total<-as.data.frame(table(subjectHits(findOverlaps(query = chic_bait_bed[queryHits(findOverlaps(query = chic_otherends_bed, subject = total_ernas))],subject =  coreclock_bed))))
  ov_cc_total$Var1<-as.integer(levels(ov_cc_total$Var1))[ov_cc_total$Var1]
  #ov_cc_total$genename<-centralcc$genename[ov_cc_total$Var1]
  ov_cc_total_1<-ov_cc_total$Freq
  names(ov_cc_total_1)<-coreclock_bed$genename[ov_cc_total$Var1]
  
  #overlap
  #overlap cc with dynamic inters
  #Overlaps of baits and other ends with core clock, to make it symmetrical
  ov_cc_dynaB<-as.data.frame(table(subjectHits(findOverlaps(query = c(BELOWgrangesconservedinalltimepoints_B, ABOVEgrangesconservedinalltimepoints_B)[queryHits(findOverlaps(query = c(BELOWgrangesconservedinalltimepoints_OE, ABOVEgrangesconservedinalltimepoints_OE), subject = total_ernas))],subject =  coreclock_bed))))
  ov_cc_dynaOE<-as.data.frame(table(subjectHits(findOverlaps(query = c(BELOWgrangesconservedinalltimepoints_OE, ABOVEgrangesconservedinalltimepoints_OE)[queryHits(findOverlaps(query = c(BELOWgrangesconservedinalltimepoints_B, ABOVEgrangesconservedinalltimepoints_B), subject = total_ernas))],subject =  coreclock_bed))))
  ov_cc_dynaB$Var1<-as.integer(levels(ov_cc_dynaB$Var1))[ov_cc_dynaB$Var1]
  #ov_cc_dynaB$genename<-coreclock_bed$genename[ov_cc_dynaB$Var1]
  ov_cc_dynaOE$Var1<-as.integer(levels(ov_cc_dynaOE$Var1))[ov_cc_dynaOE$Var1]
  ov_cc_dyna<-mergebaitandoe(ov_cc_dynaB,ov_cc_dynaOE,  coreclock_bed)
  
  #overlap with static inter
  #Overlaps of baits and other ends with core clock, to make it symmetrical
  ov_cc_staticB<-as.data.frame(table(subjectHits(findOverlaps(STATICchic_bait_bed[queryHits(findOverlaps(query = STATICchic_otherends_bed, subject = total_ernas))], coreclock_bed))))
  ov_cc_staticB$Var1<-as.integer(levels(ov_cc_staticB$Var1))[ov_cc_staticB$Var1]
  #ov_cc_staticB$genename<-centralcc$genename[ov_cc_staticB$Var1]
  
  ov_cc_staticOE<-as.data.frame(table(subjectHits(findOverlaps(STATICchic_otherends_bed[queryHits(findOverlaps(query = STATICchic_bait_bed, subject = total_ernas))], coreclock_bed))))
  ov_cc_staticOE$Var1<-as.integer(levels(ov_cc_staticOE$Var1))[ov_cc_staticOE$Var1]
  #ov_cc_staticOE$genename<-centralcc$genename[ov_cc_staticOE$Var1]
  
  ov_cc_static<-mergebaitandoe(ov_cc_staticB,ov_cc_staticOE,  coreclock_bed)
  
 
  #how many cc-dynamicinters are in cc-staticinters
  t<-list(ov_cc_total_1,ov_cc_dyna,ov_cc_static)
  names(t)<-c("All_inters", "Dynamic", "Static")
  return(t)}


obs_cc_sep_totalernas_belowabove<-coreclock_separately_difinters_overnas(allreadsbelow_difinter_timepoint,allreadsabove_difinter_timepoint, coreclock_bed, total_ernas)

obs_cc_sep_oscernas_belowabove<-coreclock_separately_difinters_overnas(allreadsbelow_difinter_timepoint,allreadsabove_difinter_timepoint, coreclock_bed, osc_ernas)


#########Expected number of interactions of non coreclock circadian genes:
#circpromall<-bedfile(circpromphases, columnnames)
circpromphases<-import("inputs_for_scripts/circproms/MFM_RNAseq/HindIIIfragments_circadiangenes_MFMRNAseq.bed", format = "BED")
circpromallintrons<-bedfile(circpromphasesINTRONS, columnnames)
circprom_notcoreclock<-circpromphases[!(circpromphases$name%in%coreclock_bed$name)]
circpromINTRONS_notcoreclock<-circpromallintrons[!(1:271 %in% queryHits(findOverlaps( circpromallintrons, coreclock_bed)))]


mergebaitandoe_forexp<-function(tableinters_bait,tableinters_OE,  coreclock){
  t<-tableinters_bait$Freq
  names(t)<-coreclock$name[tableinters_bait$Var1]
  t2<-tableinters_OE$Freq
  names(t2)<-coreclock$name[tableinters_OE$Var1]
  #Merge bait and OE numbers
  t3<-sapply(coreclock$name, function(x){sum(c(t,t2)[names(c(t,t2))%in% x])})
  return(t3)
}
coreclock_separately_difinters_expected<-function(difintersbelow,difintersabove, circprom_notcoreclock, iterations, numberofcc){
  t<-sapply(1:iterations, function(x){
    rd<-sample(1:length(circprom_notcoreclock), numberofcc, replace = F)
    BELOWgrangesconservedinalltimepoints_B<-GRanges(seqnames = difintersbelow[,1], ranges = IRanges(difintersbelow[,2], difintersbelow[,3]), ZT0_difint= difintersbelow[,7], ZT6_difint= difintersbelow[,8], ZT12_difint= difintersbelow[,9], ZT18_difint= difintersbelow[,10])
    BELOWgrangesconservedinalltimepoints_OE<-GRanges(seqnames = difintersbelow[,4], ranges = IRanges(difintersbelow[,5], difintersbelow[,6]),ZT0_difint= difintersbelow[,7], ZT6_difint= difintersbelow[,8], ZT12_difint= difintersbelow[,9], ZT18_difint= difintersbelow[,10] )
    
    ABOVEgrangesconservedinalltimepoints_B<-GRanges(seqnames = difintersabove[,1], ranges = IRanges(difintersabove[,2], difintersabove[,3]), ZT0_difint= difintersabove[,7], ZT6_difint= difintersabove[,8], ZT12_difint= difintersabove[,9], ZT18_difint= difintersabove[,10])
    ABOVEgrangesconservedinalltimepoints_OE<-GRanges(seqnames = difintersabove[,4], ranges = IRanges(difintersabove[,5], difintersabove[,6]),ZT0_difint= difintersabove[,7], ZT6_difint= difintersabove[,8], ZT12_difint= difintersabove[,9], ZT18_difint= difintersabove[,10] )
    
    #Overlaps cc with total inters  
    ov_cc_total<-as.data.frame(table(subjectHits(findOverlaps(query = chic_bait_bed,subject =  circprom_notcoreclock[rd]))))
    ov_cc_total$Var1<-as.integer(levels(ov_cc_total$Var1))[ov_cc_total$Var1]
    ov_cc_total_1<-ov_cc_total$Freq
    #names(ov_cc_total_1)<-coreclock_bed$genename[ov_cc_total$Var1]
    
    #overlap cc with dynamic inters
    #Overlaps of baits and other ends with core clock, to make it symmetrical
    ov_cc_dynaB<-as.data.frame(table(subjectHits(findOverlaps(query = (c(BELOWgrangesconservedinalltimepoints_B, ABOVEgrangesconservedinalltimepoints_B)),subject =  circprom_notcoreclock[rd]))))
    ov_cc_dynaOE<-as.data.frame(table(subjectHits(findOverlaps(query = (c(BELOWgrangesconservedinalltimepoints_OE, ABOVEgrangesconservedinalltimepoints_OE)),subject =  circprom_notcoreclock[rd]))))
    ov_cc_dynaB$Var1<-as.integer(levels(ov_cc_dynaB$Var1))[ov_cc_dynaB$Var1]
    ov_cc_dynaOE$Var1<-as.integer(levels(ov_cc_dynaOE$Var1))[ov_cc_dynaOE$Var1]
    ov_cc_dyna<-mergebaitandoe_forexp(ov_cc_dynaB,ov_cc_dynaOE,  circprom_notcoreclock[rd])
    #t1<-ov_cc_dynaB$Var1
    #names(t1)<-circprom_notcoreclock[rd, "genename"][ov_cc_dynaB$Var1]
    #t2<-tableinters_OE$Var1
    #names(t2)<-circprom_notcoreclock[rd, "genename"][tableinters_OE$Var1]
    #ov_cc_dyna<-sapply(circprom_notcoreclock[rd, "genename"], function(x){sum(c(t,t2)[names(c(t,t2))%in% x])})
    #overlap with static inter
    #Overlaps of baits and other ends with core clock, to make it symmetrical
    ov_cc_staticB<-as.data.frame(table(subjectHits(findOverlaps(STATICchic_bait_bed, circprom_notcoreclock[rd]))))
    ov_cc_staticB$Var1<-as.integer(levels(ov_cc_staticB$Var1))[ov_cc_staticB$Var1]
    ov_cc_staticOE<-as.data.frame(table(subjectHits(findOverlaps(STATICchic_otherends_bed, circprom_notcoreclock[rd]))))
    ov_cc_staticOE$Var1<-as.integer(levels(ov_cc_staticOE$Var1))[ov_cc_staticOE$Var1]
    ov_cc_static<-mergebaitandoe_forexp(ov_cc_staticB,ov_cc_staticOE,  circprom_notcoreclock[rd])
    #t1<-ov_cc_staticB$Var1
    #names(t1)<-circprom_notcoreclock[rd, "genename"][ov_cc_staticB$Var1]
    #t2<-ov_cc_staticOE$Var1
    #names(t2)<-circprom_notcoreclock[rd, "genename"][ov_cc_staticOE$Var1]
    #ov_cc_static<-sapply(circprom_notcoreclock[rd, "genename"], function(x){sum(c(t,t2)[names(c(t,t2))%in% x])})
    
    #how many cc-dynamicinters are in cc-staticinters
    t<-list(ov_cc_total_1,ov_cc_dyna,ov_cc_static)
    names(t)<-c("All_inters", "Dynamic", "Static")
    return(t)})
  return(t)}

exp_cc_sep_belowabove<-coreclock_separately_difinters_expected(allreadsbelow_difinter_timepoint,allreadsabove_difinter_timepoint, circprom_notcoreclock, iterations = 100, numberofcc = 17)

#Report median of 100 iters for 17 genes
#All inters
#NOTE: if we use the sum of interactions of 17 genes and the divide by 17 we get a median of ~18 genes
##But if I get the median of 17 genes of all iterations then it becomes ~13
median(sapply(1:100, function(iter){median(unlist(exp_cc_sep_belowabove[1,iter]))}))
median(sapply(1:100, function(x){sum(unlist(exp_cc_sep_belowabove[1,x]))}))
median(sapply(1:100, function(x){sum(unlist(exp_cc_sep_belowabove[1,x]))}))/17
#Dynamic
median(sapply(1:100, function(iter){median(unlist(exp_cc_sep_belowabove[2,iter]))}))
median(sapply(1:100, function(x){sum(unlist(exp_cc_sep_belowabove[2,x]))}))
median(sapply(1:100, function(x){sum(unlist(exp_cc_sep_belowabove[2,x]))}))/17

#Static
median(sapply(1:100, function(iter){median(unlist(exp_cc_sep_belowabove[3,iter]))}))
median(sapply(1:100, function(x){sum(unlist(exp_cc_sep_belowabove[3,x]))}))
median(sapply(1:100, function(x){sum(unlist(exp_cc_sep_belowabove[3,x]))}))/17


#Export:
allinters<-do.call(rbind.data.frame, exp_cc_sep_belowabove[1,])
colnames(allinters)<-NULL
write.table(x = allinters, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/2020_07_08_Expected_interactions_coreclock_100iter_Totalinters.txt", quote = F, sep = "\t", row.names = F, col.names = F)

Dynamic<-do.call(rbind.data.frame, exp_cc_sep_belowabove[2,])
colnames(Dynamic)<-NULL
write.table(x = Dynamic, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/2020_07_08_Expected_interactions_coreclock_100iter_Dynamic.txt", quote = F, sep = "\t", row.names = F, col.names = F)

Static<-do.call(rbind.data.frame, exp_cc_sep_belowabove[3,])
colnames(Static)<-NULL
write.table(x = Static, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/2020_07_06_Expected_interactions_coreclock_100iter_Static.txt", quote = F, sep = "\t", row.names = F, col.names = F)

#Expected per random 17 genes
#exp_cc_sep_belowabove_allinters<-read.csv("inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/2020_07_06_Expected_interactions_coreclock_100iter_Totalinters.txt",sep = "\t",header = T)
#colnames(exp_cc_sep_belowabove_allinters)<-NULL
#median(as.matrix(as.numeric(apply(exp_cc_sep_belowabove_allinters, 2, function(x){median(x, na.rm=T)}))))
#exp_cc_sep_belowabove_dynamic<-read.csv("inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/2020_07_06_Expected_interactions_coreclock_100iter_Dynamic.txt",sep = "\t",header = T)
#median(as.matrix(as.numeric(apply(exp_cc_sep_belowabove_dynamic, 2, function(x){median(x, na.rm=T)}))))
#exp_cc_sep_belowabove_static<-read.csv("inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/2020_07_06_Expected_interactions_coreclock_100iter_Static.txt",sep = "\t",header = T)
#median(as.matrix(as.numeric(apply(exp_cc_sep_belowabove_static, 2, function(x){median(x, na.rm=T)}))))



#Expected for eRNAs

coreclock_separately_difinters_overnas_expected<-function(difintersbelow,difintersabove, circprom_notcoreclock, iterations, total_ernas){
  t<-sapply(1:iterations, function(x){
    rd<-sample(1:length(circprom_notcoreclock), 17)
    BELOWgrangesconservedinalltimepoints_B<-GRanges(seqnames = difintersbelow[,1], ranges = IRanges(difintersbelow[,2], difintersbelow[,3]), ZT0_difint= difintersbelow[,7], ZT6_difint= difintersbelow[,8], ZT12_difint= difintersbelow[,9], ZT18_difint= difintersbelow[,10])
    BELOWgrangesconservedinalltimepoints_OE<-GRanges(seqnames = difintersbelow[,4], ranges = IRanges(difintersbelow[,5], difintersbelow[,6]),ZT0_difint= difintersbelow[,7], ZT6_difint= difintersbelow[,8], ZT12_difint= difintersbelow[,9], ZT18_difint= difintersbelow[,10] )
    
    ABOVEgrangesconservedinalltimepoints_B<-GRanges(seqnames = difintersabove[,1], ranges = IRanges(difintersabove[,2], difintersabove[,3]), ZT0_difint= difintersabove[,7], ZT6_difint= difintersabove[,8], ZT12_difint= difintersabove[,9], ZT18_difint= difintersabove[,10])
    ABOVEgrangesconservedinalltimepoints_OE<-GRanges(seqnames = difintersabove[,4], ranges = IRanges(difintersabove[,5], difintersabove[,6]),ZT0_difint= difintersabove[,7], ZT6_difint= difintersabove[,8], ZT12_difint= difintersabove[,9], ZT18_difint= difintersabove[,10] )
    
    #Overlaps cc with total inters  
    ov_cc_total<-as.data.frame(table(subjectHits(findOverlaps(query = chic_bait_bed[queryHits(findOverlaps(query = chic_otherends_bed, subject = total_ernas))],subject =  circprom_notcoreclock[rd]))))
    ov_cc_total$Var1<-as.integer(levels(ov_cc_total$Var1))[ov_cc_total$Var1]
    ov_cc_total_1<-ov_cc_total$Freq
    #names(ov_cc_total_1)<-coreclock_bed$genename[ov_cc_total$Var1]
    
    #overlap cc with dynamic inters
    #Overlaps of baits and other ends with core clock, to make it symmetrical
    ov_cc_dynaB<-as.data.frame(table(subjectHits(findOverlaps(query = c(BELOWgrangesconservedinalltimepoints_B, ABOVEgrangesconservedinalltimepoints_B)[queryHits(findOverlaps(query = c(BELOWgrangesconservedinalltimepoints_OE, ABOVEgrangesconservedinalltimepoints_OE), subject = total_ernas))],subject =  circprom_notcoreclock[rd]))))
    ov_cc_dynaOE<-as.data.frame(table(subjectHits(findOverlaps(query = c(BELOWgrangesconservedinalltimepoints_OE, ABOVEgrangesconservedinalltimepoints_OE)[queryHits(findOverlaps(query = c(BELOWgrangesconservedinalltimepoints_B, ABOVEgrangesconservedinalltimepoints_B), subject = total_ernas))],subject =  circprom_notcoreclock[rd]))))
    ov_cc_dynaB$Var1<-as.integer(levels(ov_cc_dynaB$Var1))[ov_cc_dynaB$Var1]
    ov_cc_dynaOE$Var1<-as.integer(levels(ov_cc_dynaOE$Var1))[ov_cc_dynaOE$Var1]
    ov_cc_dyna<-mergebaitandoe_forexp(ov_cc_dynaB,ov_cc_dynaOE,  circprom_notcoreclock[rd])
    
    #overlap with static inter
    #Overlaps of baits and other ends with core clock, to make it symmetrical
    ov_cc_staticB<-as.data.frame(table(subjectHits(findOverlaps(STATICchic_bait_bed[queryHits(findOverlaps(query = STATICchic_otherends_bed, subject = total_ernas))], circprom_notcoreclock[rd]))))
    ov_cc_staticB$Var1<-as.integer(levels(ov_cc_staticB$Var1))[ov_cc_staticB$Var1]
    ov_cc_staticOE<-as.data.frame(table(subjectHits(findOverlaps(STATICchic_otherends_bed[queryHits(findOverlaps(query = STATICchic_bait_bed, subject = total_ernas))], circprom_notcoreclock[rd]))))
    ov_cc_staticOE$Var1<-as.integer(levels(ov_cc_staticOE$Var1))[ov_cc_staticOE$Var1]
    ov_cc_static<-mergebaitandoe_forexp(ov_cc_staticB,ov_cc_staticOE,  circprom_notcoreclock[rd])
    
    #how many cc-dynamicinters are in cc-staticinters
    t<-list(ov_cc_total_1,ov_cc_dyna,ov_cc_static)
    names(t)<-c("All_inters", "Dynamic", "Static")
    return(t)})
  return(t)}


exp_cc_totalernas_belowabove<-coreclock_separately_difinters_overnas_expected(allreadsbelow_difinter_timepoint,allreadsabove_difinter_timepoint, circprom_notcoreclock, 100, total_ernas)

exp_cc_oscernas_belowabove<-coreclock_separately_difinters_overnas_expected(allreadsbelow_difinter_timepoint,allreadsabove_difinter_timepoint, circprom_notcoreclock, 100, osc_ernas)

#Report median of 100 iters for 17 genes
#All inters
median(sapply(1:100, function(iter){median(unlist(exp_cc_totalernas_belowabove[1,iter]))}))
median(sapply(1:100, function(iter){median(unlist(exp_cc_oscernas_belowabove[1,iter]))}))

#Dynamic
median(sapply(1:100, function(iter){median(unlist(exp_cc_totalernas_belowabove[2,iter]))}))
median(sapply(1:100, function(iter){median(unlist(exp_cc_oscernas_belowabove[2,iter]))}))

#Static
median(sapply(1:100, function(iter){median(unlist(exp_cc_totalernas_belowabove[3,iter]))}))
median(sapply(1:100, function(iter){median(unlist(exp_cc_oscernas_belowabove[3,iter]))}))

#Export:
#Total eRNAs
allinters<-do.call(rbind.data.frame, exp_cc_totalernas_belowabove[1,])
write.table(x = allinters, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/2020_07_06_Expected_interactions_coreclock_100iter_ovtotalernas_Totalinters.txt", quote = F, sep = "\t", row.names = F, col.names = T)

Dynamic<-do.call(rbind.data.frame, exp_cc_totalernas_belowabove[2,])
write.table(x = Dynamic, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/2020_07_06_Expected_interactions_coreclock_100iter_ovtotalernas_Dynamic.txt", quote = F, sep = "\t", row.names = F, col.names = T)

Static<-do.call(rbind.data.frame, exp_cc_totalernas_belowabove[3,])
write.table(x = Static, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/2020_07_06_Expected_interactions_coreclock_100iter_ovtotalernas_Static.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#osc eRNAs
allinters<-do.call(rbind.data.frame, exp_cc_oscernas_belowabove[1,])
write.table(x = allinters, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/2020_07_06_Expected_interactions_coreclock_100iter_ovoscernas_Totalinters.txt", quote = F, sep = "\t", row.names = F, col.names = T)

Dynamic<-do.call(rbind.data.frame, exp_cc_oscernas_belowabove[2,])
write.table(x = Dynamic, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/2020_07_06_Expected_interactions_coreclock_100iter_ovoscernas_Dynamic.txt", quote = F, sep = "\t", row.names = F, col.names = T)

Static<-do.call(rbind.data.frame, exp_cc_oscernas_belowabove[3,])
write.table(x = Static, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/2020_07_06_Expected_interactions_coreclock_100iter_ovoscernas_Static.txt", quote = F, sep = "\t", row.names = F, col.names = T)



#all the previous blots show that the non core genes have slightly higher interactions, but the expected interactions from the table show even one less interactoins
library(tidyverse)
library(gdata)
t<-read.csv("inputs_for_scripts/countsinteractionspercircproms_table.txt", sep="\t", header=F)
ccinters<-t[t$V4 %in%coreclock_bed$name,c(4:5)] %>% group_by(V4) %>% summarize(interactions = sum(V5))

nonccinters<-t[!(t$V4 %in%coreclock_bed$name),c(4:5)] %>% group_by(V4) %>% summarize(interactions = sum(V5))

t2<-cbindX(as.matrix(ccinters[,-1]), as.matrix(nonccinters[,-1]))
colnames(t2)<-c("Coreclock\ngenes", "Circadian\nnoncoreclock")
boxplot(t2, las=2,ylab="Number of total interactions")


t3<-cbindX(as.matrix(ccinters), as.matrix(nonccinters))
colnames(t3)<-c("Coreclockgenes", "NumberofTotalInteractions","Circadiannoncoreclock","NumberofTotalInteractions")
write.table(x = t3, file = "inputs_for_scripts/countsinteractionspercircproms_table_summarized.txt", quote = F, sep = "\t" ,row.names = F, col.names = T)

#########################
####2020-07-8 Repeat with fang GROseq to confirm results: 1098 fang; 13 cc
fangcircproms<- read.csv("/inputs_for_scripts/circproms/Updated/HindIIIfragments_circadiangenes_fang_phases.bed", header = F, sep="\t")
fangcircproms<- bedfile(fangcircproms[,-4], columnnames)
t<- read.csv("/inputs_for_scripts/circproms/Updated/HindIIIfragments_circadiangenes_fang_phases.bed", header = F, sep="\t")
mcols(fangcircproms)$name<-t[,5]
fangcircproms$name<-as.character(levels(fangcircproms$name))[fangcircproms$name]
fangcircproms_coreclock<-fangcircproms[unique(match(coreclock_bed$name,fangcircproms$name ))[!is.na(unique(match(coreclock_bed$name,fangcircproms$name )))]]
fangcircproms_notcoreclock<-fangcircproms[fangcircproms$name%in% setdiff(fangcircproms$name,coreclock_bed$name),]

obs_cc_sep_belowabove_FANG<-coreclock_separately_difinters(allreadsbelow_difinter_timepoint,allreadsabove_difinter_timepoint, fangcircproms_coreclock )#

#Put in a df
n<-names(obs_cc_sep_belowabove_FANG$All_inters)
obs_cc_sep_belowabove_FANG<-do.call(cbind.data.frame, obs_cc_sep_belowabove_FANG)
rownames(obs_cc_sep_belowabove_FANG)<-n
obs_cc_sep_belowabove_FANG #Same because they re the same fragments


exp_cc_sep_belowaboveFANG<-coreclock_separately_difinters_expected(allreadsbelow_difinter_timepoint,allreadsabove_difinter_timepoint, fangcircproms_notcoreclock,iterations = 100, numberofcc = 17)

#Report median of 100 iters for 17 genes
#All inters
#NOTE: if we use the sum of interactions of 17 genes and the divide by 17 we get a median of ~18 genes
##But if I get the median of 17 genes of all iterations then it becomes ~13
median(sapply(1:100, function(iter){median(unlist(exp_cc_sep_belowaboveFANG[1,iter]))}))
median(sapply(1:100, function(x){sum(unlist(exp_cc_sep_belowaboveFANG[1,x]))}))
median(sapply(1:100, function(x){sum(unlist(exp_cc_sep_belowaboveFANG[1,x]))}))/17
#Dynamic
median(sapply(1:100, function(iter){median(unlist(exp_cc_sep_belowaboveFANG[2,iter]))}))
median(sapply(1:100, function(x){sum(unlist(exp_cc_sep_belowaboveFANG[2,x]))}))
median(sapply(1:100, function(x){sum(unlist(exp_cc_sep_belowaboveFANG[2,x]))}))/17

#Static
median(sapply(1:100, function(iter){median(unlist(exp_cc_sep_belowaboveFANG[3,iter]))}))
median(sapply(1:100, function(x){sum(unlist(exp_cc_sep_belowaboveFANG[3,x]))}))
median(sapply(1:100, function(x){sum(unlist(exp_cc_sep_belowaboveFANG[3,x]))}))/17


#Export:
allinters<-do.call(rbind.data.frame, exp_cc_sep_belowaboveFANG[1,])
colnames(allinters)<-NULL
write.table(x = allinters, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/2020_07_08_Expected_interactions_coreclock_100iter_TotalintersFANG.txt", quote = F, sep = "\t", row.names = F, col.names = F)

Dynamic<-do.call(rbind.data.frame, exp_cc_sep_belowaboveFANG[2,])
colnames(Dynamic)<-NULL
write.table(x = Dynamic, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/2020_07_08_Expected_interactions_coreclock_100iter_DynamicFANG.txt", quote = F, sep = "\t", row.names = F, col.names = F)

Static<-do.call(rbind.data.frame, exp_cc_sep_belowaboveFANG[3,])
colnames(Static)<-NULL
write.table(x = Static, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/2020_07_06_Expected_interactions_coreclock_100iter_StaticFANG.txt", quote = F, sep = "\t", row.names = F, col.names = F)

####Count interactions per circproms for fang


#all the previous blots show that the non core genes have slightly higher interactions, but the expected interactions from the table show even one less interactoins
library(tidyverse)
library(gdata)
t<-read.csv("inputs_for_scripts/countsinteractionspercircproms_tableFANG.txt", sep="\t", header=F)
ccinters<-t[t$V5 %in%coreclock_bed$name,c(5,7)] %>% group_by(V5) %>% summarize(interactions = sum(V7))

nonccinters<-t[!(t$V5 %in%coreclock_bed$name),c(5,7)] %>% group_by(V5) %>% summarize(interactions = sum(V7))

t2<-cbindX(as.matrix(ccinters[,-1]), as.matrix(nonccinters[,-1]))
colnames(t2)<-c("Coreclock\ngenes", "Circadian\nnoncoreclock")
boxplot(t2, las=2,ylab="Number of total interactions")


t3<-cbindX(as.matrix(ccinters), as.matrix(nonccinters))
colnames(t3)<-c("Coreclockgenes", "NumberofTotalInteractions","Circadiannoncoreclock","NumberofTotalInteractions")
write.table(x = t3, file = "inputs_for_scripts/countsinteractionspercircproms_table_summarizedFANG.txt", quote = F, sep = "\t" ,row.names = F, col.names = T)


#########2020-07-11 Try different set of core clocks
library(data.table)
classiccc<-c("Clock", "Npas2", "Arntl", "Arntl2", "Per1", "Per2", "Per3", "Cry1", "Cry2","Nr1d1", "Nr1d2", "Rora","Rorb", "Rorc", "Fbxl3", "Csnk1d", "Csnk1e")
ml_cc<-c("Clock", "Npas2", "Arntl", "Arntl2", "Per1", "Per2", "Per3", "Cry1", "Cry2","Nr1d1", "Nr1d2", "Rora","Rorb", "Rorc", "Fbxl3", "Csnk1d", "Csnk1e","Tars","Nfil3", "Wee1", "Psen2","Clpx","Bnip3", "Hdac11", "Gm129", "Tef", "Caprin1")
HindBaits<-bedfile(fread("inputs_for_scripts/circproms/Updated/HindIIIfragments_ov_baits.bed")[,1:3], columnnames)
Baits<-fread("inputs_for_scripts/Baits.txt")

Baits<-GRanges(seqnames = Baits$V1, ranges = IRanges(Baits$V2, Baits$V3),genename=Baits$V5)

BCC<-(sort(Baits[unlist(sapply(classiccc, function(x){grep(x, Baits$genename)}))]) )
bhcc<-HindBaits[subjectHits(findOverlaps(query = BCC, subject = HindBaits))]
bhcc$name<-sub(pattern = "-\\d0\\d",replacement = "",lapply(strsplit(unlist(mcols(BCC[queryHits(findOverlaps(query = BCC, subject = HindBaits))])), ","), function(x){x[1]}))
bhcc<-unique(bhcc)

BMLCC<-(sort(Baits[unlist(sapply(ml_cc, function(x){grep(fixed = T,x, Baits$genename)}))]) )
bhmlcc<-HindBaits[subjectHits(findOverlaps(query = BMLCC, subject = HindBaits))]
bhmlcc$name<-sub(pattern = "-\\d0\\d",replacement = "",lapply(strsplit(unlist(mcols(BMLCC[queryHits(findOverlaps(query = BMLCC, subject = HindBaits))])), ","), function(x){x[1]}))
bhmlcc<-unique(bhmlcc)
bhmlcc<-bhmlcc[bhmlcc$name%in%ml_cc]
#Do the analysis by filtering: All genes in chic and genes in MFMrnaseq+chic
#Total inters (total, dynamic, static; dynamic based from diffinters that have MFMRNAseq circproms)
#---------------------#Obs#---------------------#
#classicsallchic_obs_cc_sep_belowabove<-coreclock_separately_difinters(difintersbelow = allreadsbelow_difinter_timepoint,difintersabove = allreadsabove_difinter_timepoint, coreclock_bed = )
#11/17 ClassicCC on MFM:Done
classicsonMFMRNAseq_obs_cc_sep_belowabove<-coreclock_separately_difinters(difintersbelow = allreadsbelow_difinter_timepoint,difintersabove = allreadsabove_difinter_timepoint, coreclock_bed =circpromphases[circpromphases$name %in% classiccc] )
listtodf(classicsonMFMRNAseq_obs_cc_sep_belowabove)
#16/17 ClassicCC ALL: Done
classicsall_obs_cc_sep_belowabove<-coreclock_separately_difinters(difintersbelow = allreadsbelow_difinter_timepoint,difintersabove = allreadsabove_difinter_timepoint, coreclock_bed =bhcc )
listtodf(classicsall_obs_cc_sep_belowabove)

#19/27 MachineLearning on MFM: Done
machlearnonMFMRNAseq_obs_cc_sep_belowabove<-coreclock_separately_difinters(difintersbelow = allreadsbelow_difinter_timepoint,difintersabove = allreadsabove_difinter_timepoint, coreclock_bed =circpromphases[circpromphases$name %in% ml_cc] )
listtodf(machlearnonMFMRNAseq_obs_cc_sep_belowabove)

#26/27 MachineLearning ALL: Done
machlearnall_obs_cc_sep_belowabove<-coreclock_separately_difinters(difintersbelow = allreadsbelow_difinter_timepoint,difintersabove = allreadsabove_difinter_timepoint, coreclock_bed =bhmlcc )
listtodf(machlearnall_obs_cc_sep_belowabove)


#---------------------#Expected#---------------------#
#---11/17 ClassicCC on MFM: Done
circprom_notcoreclock<-circpromphases[!(circpromphases$name%in%circpromphases[circpromphases$name %in% classiccc]$name)]
classicsonMFMRNAse_exp_cc_sep_belowabove<-coreclock_separately_difinters_expected(allreadsbelow_difinter_timepoint,allreadsabove_difinter_timepoint, circprom_notcoreclock, iterations = 100, numberofcc = length(circpromphases[circpromphases$name %in% classiccc]))
  #Export:
allinters<-do.call(rbind.data.frame, classicsonMFMRNAse_exp_cc_sep_belowabove[1,])
colnames(allinters)<-NULL
write.table(x = allinters, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/Anafi_2014/2020_07_11_Expected_interactions_coreclock_100iter_Totalinters_classicsonMFMRNAse.txt", quote = F, sep = "\t", row.names = F, col.names = F)

Dynamic<-do.call(rbind.data.frame, classicsonMFMRNAse_exp_cc_sep_belowabove[2,])
colnames(Dynamic)<-NULL
write.table(x = Dynamic, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/Anafi_2014/2020_07_11_Expected_interactions_coreclock_100iter_Dynamic_classicsonMFMRNAse.txt", quote = F, sep = "\t", row.names = F, col.names = F)

Static<-do.call(rbind.data.frame, classicsonMFMRNAse_exp_cc_sep_belowabove[3,])
colnames(Static)<-NULL
write.table(x = Static, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/Anafi_2014/2020_07_11_Expected_interactions_coreclock_100iter_Static_classicsonMFMRNAse.txt", quote = F, sep = "\t", row.names = F, col.names = F)
summarize_Exp(allinters,Dynamic,Static)

#---16/17 ClassicCC all: Done
circprom_notcoreclock<-circpromphases[!(circpromphases$name%in%bhcc$name)]
classicsall_exp_cc_sep_belowabove<-coreclock_separately_difinters_expected(allreadsbelow_difinter_timepoint,allreadsabove_difinter_timepoint, circprom_notcoreclock, iterations = 100, numberofcc = length(bhcc))
#Export:
allinters<-do.call(rbind.data.frame, classicsall_exp_cc_sep_belowabove[1,])
colnames(allinters)<-NULL
write.table(x = allinters, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/Anafi_2014/2020_07_11_Expected_interactions_coreclock_100iter_Totalinters_classicsall.txt", quote = F, sep = "\t", row.names = F, col.names = F)

Dynamic<-do.call(rbind.data.frame, classicsall_exp_cc_sep_belowabove[2,])
colnames(Dynamic)<-NULL
write.table(x = Dynamic, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/Anafi_2014/2020_07_11_Expected_interactions_coreclock_100iter_Dynamic_classicsall.txt", quote = F, sep = "\t", row.names = F, col.names = F)

Static<-do.call(rbind.data.frame, classicsall_exp_cc_sep_belowabove[3,])
colnames(Static)<-NULL
write.table(x = Static, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/Anafi_2014/2020_07_11_Expected_interactions_coreclock_100iter_Static_classicsall.txt", quote = F, sep = "\t", row.names = F, col.names = F)
summarize_Exp(allinters,Dynamic,Static)

#--------19/27 Machine Learning on MFM: Done
circprom_notcoreclock<-circpromphases[!(circpromphases$name%in%circpromphases[circpromphases$name %in% ml_cc]$name)]
machlearnonMFMRNAseq_exp_cc_sep_belowabove<-coreclock_separately_difinters_expected(allreadsbelow_difinter_timepoint,allreadsabove_difinter_timepoint, circprom_notcoreclock, iterations = 100, numberofcc = length(circpromphases[circpromphases$name %in% ml_cc]))
#Export:
allinters<-do.call(rbind.data.frame, machlearnonMFMRNAseq_exp_cc_sep_belowabove[1,])
colnames(allinters)<-NULL
write.table(x = allinters, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/Anafi_2014/2020_07_11_Expected_interactions_coreclock_100iter_Totalinters_machlearnonMFMRNAseq.txt", quote = F, sep = "\t", row.names = F, col.names = F)

Dynamic<-do.call(rbind.data.frame, machlearnonMFMRNAseq_exp_cc_sep_belowabove[2,])
colnames(Dynamic)<-NULL
write.table(x = Dynamic, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/Anafi_2014/2020_07_11_Expected_interactions_coreclock_100iter_Dynamic_machlearnonMFMRNAseq.txt", quote = F, sep = "\t", row.names = F, col.names = F)

Static<-do.call(rbind.data.frame, machlearnonMFMRNAseq_exp_cc_sep_belowabove[3,])
colnames(Static)<-NULL
write.table(x = Static, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/Anafi_2014/2020_07_11_Expected_interactions_coreclock_100iter_Static_machlearnonMFMRNAseq.txt", quote = F, sep = "\t", row.names = F, col.names = F)
summarize_Exp(allinters,Dynamic,Static)


#--------26/27 Machine Learning all: Done
circprom_notcoreclock<-circpromphases[!(circpromphases$name%in%bhmlcc$name)]
machlearnall_exp_cc_sep_belowabove<-coreclock_separately_difinters_expected(allreadsbelow_difinter_timepoint,allreadsabove_difinter_timepoint, circprom_notcoreclock, iterations = 100, numberofcc = length(bhmlcc))
#Export:
allinters<-do.call(rbind.data.frame, machlearnall_exp_cc_sep_belowabove[1,])
colnames(allinters)<-NULL
write.table(x = allinters, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/Anafi_2014/2020_07_11_Expected_interactions_coreclock_100iter_Totalinters_machlearnall.txt", quote = F, sep = "\t", row.names = F, col.names = F)

Dynamic<-do.call(rbind.data.frame, machlearnall_exp_cc_sep_belowabove[2,])
colnames(Dynamic)<-NULL
write.table(x = Dynamic, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/Anafi_2014/2020_07_11_Expected_interactions_coreclock_100iter_Dynamic_machlearnall.txt", quote = F, sep = "\t", row.names = F, col.names = F)

Static<-do.call(rbind.data.frame, machlearnall_exp_cc_sep_belowabove[3,])
colnames(Static)<-NULL
write.table(x = Static, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/Core_clock/Anafi_2014/2020_07_11_Expected_interactions_coreclock_100iter_Static_machlearnall.txt", quote = F, sep = "\t", row.names = F, col.names = F)

summarize_Exp(allinters,Dynamic,Static)

#inters with eRNAs 
#obs_cc_sep_totalernas_belowabove<-coreclock_separately_difinters_overnas(allreadsbelow_difinter_timepoint,allreadsabove_difinter_timepoint, coreclock_bed, total_ernas)
#obs_cc_sep_oscernas_belowabove<-coreclock_separately_difinters_overnas(allreadsbelow_difinter_timepoint,allreadsabove_difinter_timepoint, coreclock_bed, osc_ernas)
summarize_Exp<-function(allinters, Dynamic,Static){
  c("Allinters"=median(apply(allinters, 1, sum)/dim(allinters)[2]),
    "Dynamic"=median(apply(Dynamic, 1, sum)/dim(Dynamic)[2]),
    "Static"=median(apply(Static, 1, sum)/dim(Static)[2]))
}

