#####################---Differential interactions---##########
#######Open data
#Open eRNAs
osc_ernas<-read.csv("eRNAs_de_novo_oscillating_phases.txt", header = T, sep = "\t")

#Divide per phases
phases<-unique(osc_ernas[,19])
#Binning per 3 hours
phases_bin<-split(phases, ceiling(seq_along(phases)/3))

#Divide the ernas per phases, create separe variables

ernas_divided_list<-lapply(phases_bin, function(x){
  tempo1<-osc_ernas[osc_ernas[,19]==x[1],c(1,2,3,8:15, 19)]
  tempo2<-osc_ernas[osc_ernas[,19]==x[2],c(1,2,3,8:15,19)]
  tempo3<-osc_ernas[osc_ernas[,19]==x[3],c(1,2,3,8:15,19)]
  tempo4<-rbind(tempo1, tempo2, tempo3)
  return(tempo4)
})
#8 items in list with the ernas per phase (each 3 hrs)
#names for the new phases
phases_bin_names<-sapply(1:8, function(x){paste(phases_bin[[x]][1], phases_bin[[x]][2], phases_bin[[x]][3], sep = "_")})

names(ernas_divided_list) <- paste("df_ernas", phases_bin_names, sep = "")
list2env(ernas_divided_list , envir = .GlobalEnv)
remove(ernas_divided_list)

#Transform to Granges objects
list_ernasperphase<-sapply(mixedsort(ls(pattern = "df_ernas", sorted = F), decreasing = F), function(x){list(x)})

list_ernasperphase_bed<-lapply(list_ernasperphase, function(x){
  t<-get(x)
  t1<-GRanges(seqnames= Rle(t[,1]), ranges = IRanges(t[,2], t[,3]), phase=t[,12])
  return(t1)
})



#Open circproms


circpromphases<-read.csv("inputs_for_scripts/circproms/Updated/HindIIIfragments_circadiangenes_MFMRNAseq_phases.bed", header = F, sep="\t")

circpromphasesINTRONS<-read.csv("inputs_for_scripts/circproms/Updated/HindIIIfragments_circadiangenesonlyintrons_MFMRNAseq_phases.bed", header = F, sep="\t")

fangcircproms<- read.csv("inputs_for_scripts/circproms/Updated/HindIIIfragments_circadiangenes_fang_phases.bed", header = F, sep="\t")
fangcircproms<-fangcircproms[,-4]

#Divide per phases
circpromphases<-circpromphases[order(circpromphases[,5]),]
phases<-unique(circpromphases[,5])
circpromphasesfang<-fangcircproms[order(fangcircproms[,5]),]
phasesfang<-unique(circpromphasesfang[,5])
phases_bin<-split(phasesfang, ceiling(seq_along(phasesfang)/3))
#Add separaterly last phase, otherwise only last bin would have 1 phase
phases_bin[[8]]<-c(phases_bin[[8]], 23.5)
phases_bin[[9]]<-NULL


#all proms
circpromphases<-circpromphases[order(circpromphases[,5]),]
circprom_divided_list<-lapply(phases, function(x){circpromphases[circpromphases[,5]==x,]})


#introns
circpromphasesINTRONS<-circpromphasesINTRONS[order(circpromphasesINTRONS[,5]),]

#Introns
circpromINTRONS_divided_list<-lapply(phases, function(x){circpromphasesINTRONS[circpromphasesINTRONS[,5]==x,]})

#Fang circproms

circpromFANG_divided_list<-lapply(phases_bin, function(x){
  tempo1<-fangcircproms[fangcircproms[,5]==x[1],]
  tempo2<-fangcircproms[fangcircproms[,5]==x[2],]
  tempo3<-fangcircproms[fangcircproms[,5]==x[3],]
  tempo4<-fangcircproms[fangcircproms[,5]==x[4],]
  tempo5<-rbind(tempo1, tempo2, tempo3, tempo4)
  return(tempo5)
})

#Remove NAs
circprom_divided_list<-lapply(circprom_divided_list, function(x){
  t<-x
  t<-t[complete.cases(t),]
})
  #Introns
circpromINTRONS_divided_list<-lapply(circpromINTRONS_divided_list, function(x){
  t<-x
  t<-t[complete.cases(t),]
})
  #Fang
circpromFANG_divided_list<-lapply(circpromFANG_divided_list, function(x){
  t<-x
  t<-t[complete.cases(t),]
})


#8 items in list with the ernas per phase (each 3 hrs)
#names for the new phases
#phases_bin_names<-sapply(1:8, function(x){paste(phases_bin[[x]][1], phases_bin[[x]][2], phases_bin[[x]][3], phases_bin[[x]][4], sep = "_")})

#phases_bin_names<-gsub("_NA", "", phases_bin_names)
names(circprom_divided_list)<-gsub(pattern = "ZT", replacement = "" , phases)
createvars(circprom_divided_list, "df_circprom", gsub(pattern = "ZT", replacement = "" , phases))


names(circpromINTRONS_divided_list)<-gsub(pattern = "ZT", replacement = "" , phases)
createvars(circpromINTRONS_divided_list, "df_INTRONScircprom", gsub(pattern = "ZT", replacement = "" , phases))


phases_bin_names<-sapply(1:8, function(x){paste(phases_bin[[x]][1], phases_bin[[x]][2], phases_bin[[x]][3], sep = "_")})
names(circpromFANG_divided_list)<-phases_bin_names
createvars(circpromFANG_divided_list, "df_FANGcircprom", gsub(pattern = "ZT", replacement = "" , phases_bin_names))

#Transform into granges object
#biomartpos_expr replace for 
list_circpromsperphase_bed<-lapply(sapply(mixedsort(ls(pattern = "df_circprom")), list), function(x){
  t<-get(x) 
  t1<-GRanges(seqnames= Rle(t[,1]), ranges = IRanges(t[,2], t[,3]), phase=t[,5])
  return(t1)
})
  #introns
list_circpromsperphaseINTRONS_bed<-lapply(sapply(mixedsort(ls(pattern = "df_INTRONScircprom")), list), function(x){
  t<-get(x) 
  t1<-GRanges(seqnames= Rle(t[,1]), ranges = IRanges(t[,2], t[,3]), phase=t[,5])
  return(t1)
})

  #FANG
list_circpromsperphaseFANG_bed<-lapply(sapply(mixedsort(ls(pattern = "df_FANGcircprom")), list), function(x){
  t<-get(x) 
  t1<-GRanges(seqnames= Rle(t[,1]), ranges = IRanges(t[,2], t[,3]), phase=t[,5])
  return(t1)
})




#To compare with chic, change bin phases for ernas timepoints 8 to 4 timepoints

list_ernasperphase_bed_mod<-list()
list_ernasperphase_bed_mod[[1]]<-c(list_ernasperphase_bed[[1]], list_ernasperphase_bed[[2]])
list_ernasperphase_bed_mod[[3]]<-c(list_ernasperphase_bed[[5]], list_ernasperphase_bed[[6]])
list_ernasperphase_bed_mod[[2]]<-c(list_ernasperphase_bed[[3]], list_ernasperphase_bed[[4]])
list_ernasperphase_bed_mod[[4]]<-c(list_ernasperphase_bed[[7]], list_ernasperphase_bed[[8]])


list_circpromsperphaseFANG_bed_mod<-list()
list_circpromsperphaseFANG_bed_mod[[1]]<-c(list_circpromsperphaseFANG_bed[[1]],list_circpromsperphaseFANG_bed[[2]])
list_circpromsperphaseFANG_bed_mod[[3]]<-c(list_circpromsperphaseFANG_bed[[5]], list_circpromsperphaseFANG_bed[[6]])
list_circpromsperphaseFANG_bed_mod[[2]]<-c(list_circpromsperphaseFANG_bed[[3]], list_circpromsperphaseFANG_bed[[4]])
list_circpromsperphaseFANG_bed_mod[[4]]<-c(list_circpromsperphaseFANG_bed[[7]], list_circpromsperphaseFANG_bed[[8]])



#main function


#6) Differential interactions: separate UP and DOWN
#Main funtion modification
#Modify the script to retrieve up and down interactions in each interaction

difinteractions_updown<-function(inter_table){
  
  
  #### b2g edgeR differential calls
  counts <- cbind(inter_table[,7:10])
  
  # Counts filter - b2g
  circ_dist <- abs(rowMeans(inter_table[,2:3])-rowMeans(inter_table[,5:6]))/1e3
  
  #Join replicates
  #NOt for circadian
  
  #repMean = cbind(rowMeans(inter_table[,c(10,13)]),rowMeans(inter_table[,c(11,14)]),rowMeans(inter_table[,c(12,15)]))#mean of replicates
  circ_maxN <- apply(inter_table[,7:10],1,max) #max val of each 3 conditions
  circ_counts <- inter_table[circ_maxN>15 ,7:10] #modify max to 5 if using scores
  #Modify this only if using below 150kb & circ_dist>=150 ,7:10]  # filter for low N contacts
  
  ###################################################
  ###################################################
  
  # qNorm the filtered counts (take care to preserve the row names)
  countsq = normalize.quantiles(as.matrix(circ_counts))
  rownames(countsq) = rownames(circ_counts)
  circ_counts<-countsq
  
  # Run edgeR on two reps for each day
  y_circ <- DGEList(circ_counts, group=c(1,2,3,4), remove.zeros = T)
  #y_circ <- estimateCommonDisp(y_circ) #dont work, i need replicates
  #y_circ <- estimateTagwiseDisp(y_circ)
  
  # Compare two groups
  et_circ<-apply(combn(1:4, 2), 2,function(x){
    pairs<-as.character(x)
    et_circ <- exactTest(y_circ, pair=pairs, dispersion = 1e-100)
    return(et_circ)
  })
  #et_circ <- exactTest(y_circ, pair=c("1", "3"), dispersion = 0) #why only 1 and 3, make all the combinations
  edge<-lapply(et_circ, function(timepoints){
    edge_circ <- as.data.frame(topTags(timepoints, n= dim(circ_counts)[1])) #Extracts the top DE tags in a data frame for a given pair of groups, ranked by p-value or absolute log-fold change.
    edge_padj_circ <- edge_circ[edge_circ$FDR < 0.1, ]
    
    edge_padj_up_circ <- edge_padj_circ[edge_padj_circ$logFC > 1,]
    edge_padj_down_circ <- edge_padj_circ[edge_padj_circ $logFC < -1,]
    #dim(edge_padj_up_circ)
    #dim(edge_padj_down_circ)
    return(list(edge_padj_up_circ,edge_padj_down_circ))
  })
  return(edge)
  ###################################################
  ###################################################

}


#Apply function

#6.1) Apply function
load("inputs_for_scripts/circtablesinter")
#1) Filter only interactions with a circadian promoter, only work with reads
mfm_mrnas<-bedfile(read.csv("inputs_for_scripts/circproms/MFM_RNAseq/HindIIIfragments_circadiangenes_MFMRNAseq_phases.bed", header = F, sep = "\t"))

mfm_intronsmrnas<-bedfile(read.csv("inputs_for_scripts/circproms/MFM_RNAseq/HindIIIfragments_circadiangenesonlyintrons_MFMRNAseq_phases.bed", header = F, sep = "\t"))

fang_mrnas<-bedfile(read.csv("inputs_for_scripts/circproms/Updated/HindIIIfragments_circadiangenes_fang_phases.bed", header = F, sep = "\t"))

#Make symmertic? NO
intmerged_baits_bed<-GRanges(seqnames= Rle(circadian_tableinteractions[,1]), ranges = IRanges(circadian_tableinteractions[,2], circadian_tableinteractions[,3]))
intmerged_otherends_bed<-GRanges(seqnames= Rle(circadian_tableinteractions[,4]), ranges = IRanges(circadian_tableinteractions[,5], circadian_tableinteractions[,6]))

#filter only circadian interactions
intMerge_circbaits<-queryHits(findOverlaps(intmerged_baits_bed, mfm_mrnas))
circadian_tableinteractions_filteredr<-circadian_tableinteractions[intMerge_circbaits, ]

intMerge_circbaitsINTRONS<-queryHits(findOverlaps(intmerged_baits_bed, mfm_intronsmrnas))
circadian_tableinteractions_INTRONS_filteredr<-circadian_tableinteractions[intMerge_circbaitsINTRONS, ]

intMerge_circbaitsFANG<-queryHits(findOverlaps(intmerged_baits_bed, fang_mrnas))
circadian_tableinteractions_FANG_filteredr<-circadian_tableinteractions[intMerge_circbaitsFANG, ]


#Divide interactions into two categories depending on the distance between the fragments
below150kb<-abs(circadian_tableinteractions_filteredr[,5]-circadian_tableinteractions_filteredr[,2])<=150000
above150kb<-abs(circadian_tableinteractions_filteredr[,5]-circadian_tableinteractions_filteredr[,2])>150000

below150kbINTRONS<-abs(circadian_tableinteractions_INTRONS_filteredr[,5]-circadian_tableinteractions_INTRONS_filteredr[,2])<=150000
above150kbINTRONS<-abs(circadian_tableinteractions_INTRONS_filteredr[,5]-circadian_tableinteractions_INTRONS_filteredr[,2])>150000

below150kbFANG<-abs(circadian_tableinteractions_FANG_filteredr[,5]-circadian_tableinteractions_FANG_filteredr[,2])<=150000
above150kbFANG<-abs(circadian_tableinteractions_FANG_filteredr[,5]-circadian_tableinteractions_FANG_filteredr[,2])>150000


diffs_readcount_below150kb<-difinteractions_updown(circadian_tableinteractions_filteredr[below150kb,])

diffs_readcount_below150kbINTRONS<-difinteractions_updown(circadian_tableinteractions_INTRONS_filteredr[below150kbINTRONS,])

diffs_readcount_below150kbFANG<-difinteractions_updown(circadian_tableinteractions_FANG_filteredr[below150kbFANG,])

###If I change dispersion 0 to 0.01 works ? 
diffs_readcount_above150kb<-difinteractions_updown(circadian_tableinteractions_filteredr[above150kb,])
#diff_score_updown<-difinteractions_updown(circadian_tableinteractions_filtereds, 5)
diffs_readcount_above150kbINTRONS<-difinteractions_updown(circadian_tableinteractions_INTRONS_filteredr[above150kbINTRONS,])

diffs_readcount_above150kbFANG<-difinteractions_updown(circadian_tableinteractions_FANG_filteredr[above150kbFANG,])



#Save as objects. 
save(diffs_readcount_above150kb, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_diffs_readcount_above150kb2018")
save(diffs_readcount_below150kb, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_diffs_readcount_below150kb2018")

save(diffs_readcount_above150kbINTRONS , file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_diffs_readcount_above150kb2020INTRONS")
save(diffs_readcount_below150kbINTRONS, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_diffs_readcount_below150kb2020INTRONS")

save(diffs_readcount_above150kbFANG , file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_diffs_readcount_above150kb2020FANG")
save(diffs_readcount_below150kbFANG, file = "inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_diffs_readcount_below150kb2020FANG")

#load("inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_diffs_readcount_above150kb2018")
#load("inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_diffs_readcount_below150kb2018")

#Run code in 1388 to get the merge dataset of difinter for each distance
## 9.3 Correlation indicating the timepoint where the interactions is stronger

allreadsabove_difinter_timepoint
allreadsbelow_difinter_timepoint



#Boxplots

#6.2) Boxplots
#Create merge dataset of difinter
merge_uniquedifinter_upordown<-function(tablecircadianinteractions, outputadamscript){
  sapply(1:2, function(upordown){
    alldiffinter<-rbind(tablecircadianinteractions[as.numeric(rownames(outputadamscript[[1]][[upordown]])),],
                        tablecircadianinteractions[as.numeric(rownames(outputadamscript[[2]][[upordown]])),], 
                        tablecircadianinteractions[as.numeric(rownames(outputadamscript[[3]][[upordown]])),],
                        tablecircadianinteractions[as.numeric(rownames(outputadamscript[[4]][[upordown]])),],
                        tablecircadianinteractions[as.numeric(rownames(outputadamscript[[5]][[upordown]])),],
                        tablecircadianinteractions[as.numeric(rownames(outputadamscript[[6]][[upordown]])),])
    alldiffinter<-unique(alldiffinter[order(alldiffinter$V1, alldiffinter$V2, alldiffinter$V3, alldiffinter$V4, alldiffinter$V5, alldiffinter$V6),])
    #return(alldiffinter)
    heatmap.2(as.matrix(alldiffinter[,7:10]), trace = "none", density.info = "none", dendrogram="row",labRow =  "", labCol  = c("ZT0", "ZT6", "ZT12", "ZT18"), xlab="" , Colv="NA", main= as.character(dim(alldiffinter)[1]), breaks = seq(0,50,1))
    
  })

}


boxplotupordowninteractions<-function(tablecircadianinteractions, outputadamscript, labely){
  sapply(1:2, function(upordown){
  par(mfrow=c(2,3))
  sub<-c("Up", "Down")
  sapply(1:6, function(combination){
    title<-c("ZT0-6", "ZT0-12", "ZT0-18", "ZT6-12", "ZT6-18", "ZT12-18")
    boxplot(tablecircadianinteractions[as.numeric(rownames(outputadamscript[[combination]][[upordown]])),7:10], outline=F, notch=T, 
            names= c("ZT0", "ZT6", "ZT12", "ZT18"), 
            main=title[combination], col=brewer.pal(4, "BrBG"),
            ylab=labely)
    mtext(sub[upordown], 3, cex=.8)
    legend("topleft", paste("n=", dim(tablecircadianinteractions[as.numeric(rownames(outputadamscript[[combination]][[upordown]])),7:10])[1]), bty="n")
    
  })
})
}


pdf("MFMINTRONSreadsbelow150kbUPDOWNboxplots.pdf")
boxplotupordowninteractions(circadian_tableinteractions, diffs_readcount_below150kb, "")
dev.off()

pdf("MFMINTRONSreadsabove150kbUPDOWNboxplots.pdf")
boxplotupordowninteractions(circadian_tableinteractions, diffs_readcount_above150kb, "")
dev.off()

#For scores, as one combination is zero
pdf("scoresUPDOWNboxplots.pdf")

sapply(1:2, function(upordown){
  par(mfrow=c(2,3))
  sub<-c("Up", "Down")
  sapply(1:5, function(combination){
    title<-c("ZT0-6",  "ZT0-18", "ZT6-12", "ZT6-18", "ZT12-18")
    boxplot(circadian_tableinteractions_SCORES[as.numeric(rownames(diff_score_updown[-2][[combination]][[upordown]])),7:10],
            names= c("ZT0", "ZT6", "ZT12", "ZT18"), 
            main=title[combination], col=brewer.pal(6, "RdYlBu"),
            ylab="CHiCAGO score")
    mtext(sub[upordown], 3, cex=.6)
    mtext(paste("n=", dim(circadian_tableinteractions_SCORES[as.numeric(rownames(diff_score_updown[-2][[combination]][[upordown]])),7:10])[1]), 3, las=1, cex=0.6)
    
  })
})
dev.off()


#9.1  #Heatmaps


#Up BELOW 150kbs
#zero<-matrix(c(1,2,2,2,3,2), byrow = T, ncol = 2)
#six<-matrix(c(1,1,4,2,5,2), byrow = T, ncol = 2)
#twel<-matrix(c(2,1,4,1,6,2), byrow = T, ncol = 2)
#eight<-matrix(c(3,1,5,1,6,1), byrow = T, ncol = 2)

#Merge interactions in all timepoints, dont include information about the timiepoints where the inter was detected


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

  #INTRONS
allreadsabove_difinter_timepointINTRONS<-merge_uniquedifinter(circadian_tableinteractions, diffs_readcount_above150kbINTRONS)
allreadsabove_difinter_timepointINTRONS<-rbind(as.data.frame(allreadsabove_difinter_timepointINTRONS[,1]), as.data.frame(allreadsabove_difinter_timepointINTRONS[,2]))#transform into one large data frame
allreadsabove_difinter_timepointINTRONS<-unique(allreadsabove_difinter_timepointINTRONS[order(allreadsabove_difinter_timepointINTRONS$V1, allreadsabove_difinter_timepointINTRONS$V2, allreadsabove_difinter_timepointINTRONS$V3, allreadsabove_difinter_timepointINTRONS$V4, allreadsabove_difinter_timepointINTRONS$V5, allreadsabove_difinter_timepointINTRONS$V6),]) #unique 
  #FANG
allreadsabove_difinter_timepointFANG<-merge_uniquedifinter(circadian_tableinteractions, diffs_readcount_above150kbFANG)
allreadsabove_difinter_timepointFANG<-rbind(as.data.frame(allreadsabove_difinter_timepointFANG[,1]), as.data.frame(allreadsabove_difinter_timepointFANG[,2]))#transform into one large data frame
allreadsabove_difinter_timepointFANG<-unique(allreadsabove_difinter_timepointFANG[order(allreadsabove_difinter_timepointFANG$V1, allreadsabove_difinter_timepointFANG$V2, allreadsabove_difinter_timepointFANG$V3, allreadsabove_difinter_timepointFANG$V4, allreadsabove_difinter_timepointFANG$V5, allreadsabove_difinter_timepointFANG$V6),]) #unique 

#<150kb
allreadsbelow_difinter_timepoint<-merge_uniquedifinter(circadian_tableinteractions, diffs_readcount_below150kb)
allreadsbelow_difinter_timepoint<-rbind(as.data.frame(allreadsbelow_difinter_timepoint[,1]), as.data.frame(allreadsbelow_difinter_timepoint[,2]))#transform into one large data frame
allreadsbelow_difinter_timepoint<-unique(allreadsbelow_difinter_timepoint[order(allreadsbelow_difinter_timepoint$V1, allreadsbelow_difinter_timepoint$V2, allreadsbelow_difinter_timepoint$V3, allreadsbelow_difinter_timepoint$V4, allreadsbelow_difinter_timepoint$V5, allreadsbelow_difinter_timepoint$V6),]) #unique 

  #INTRONS
allreadsbelow_difinter_timepointINTRONS<-merge_uniquedifinter(circadian_tableinteractions, diffs_readcount_below150kbINTRONS)
allreadsbelow_difinter_timepointINTRONS<-rbind(as.data.frame(allreadsbelow_difinter_timepointINTRONS[,1]), as.data.frame(allreadsbelow_difinter_timepointINTRONS[,2]))#transform into one large data frame
allreadsbelow_difinter_timepointINTRONS<-unique(allreadsbelow_difinter_timepointINTRONS[order(allreadsbelow_difinter_timepointINTRONS$V1, allreadsbelow_difinter_timepointINTRONS$V2, allreadsbelow_difinter_timepointINTRONS$V3, allreadsbelow_difinter_timepointINTRONS$V4, allreadsbelow_difinter_timepointINTRONS$V5, allreadsbelow_difinter_timepointINTRONS$V6),]) #unique 
  #FANG
allreadsbelow_difinter_timepointFANG<-merge_uniquedifinter(circadian_tableinteractions, diffs_readcount_below150kbFANG)
allreadsbelow_difinter_timepointFANG<-rbind(as.data.frame(allreadsbelow_difinter_timepointFANG[,1]), as.data.frame(allreadsbelow_difinter_timepointFANG[,2]))#transform into one large data frame
allreadsbelow_difinter_timepointFANG<-unique(allreadsbelow_difinter_timepointFANG[order(allreadsbelow_difinter_timepointFANG$V1, allreadsbelow_difinter_timepointFANG$V2, allreadsbelow_difinter_timepointFANG$V3, allreadsbelow_difinter_timepointFANG$V4, allreadsbelow_difinter_timepointFANG$V5, allreadsbelow_difinter_timepointFANG$V6),]) #unique 

##Phase cor

phasecor_diffinteralltimepoints_pe_vectorphases<-function(intertableall, list_circproms){
  grangesconservedinalltimepoints_B<-GRanges(seqnames = intertableall[,1], ranges = IRanges(intertableall[,2], intertableall[,3]))
  
  grangesconservedinalltimepoints_OE<-GRanges(seqnames = intertableall[,4], ranges = IRanges(intertableall[,5], intertableall[,6]))
      
    c3<-sapply(1:4, function(x){
      index<-findOverlaps(grangesconservedinalltimepoints_OE, list_ernasperphase_bed_mod[[x]])
      mcols(list_ernasperphase_bed_mod[[x]][subjectHits(index)])
      t3<-sapply(1:4, function(y){
        index2<-findOverlaps(grangesconservedinalltimepoints_B[queryHits(index)], list_circproms[[y]])
      t1<-mcols(list_ernasperphase_bed_mod[[x]][subjectHits(index)[queryHits(index2)]])
      t2<-mcols(list_circproms[[y]][subjectHits(index2)])
      #Only for MFM datastset run the next line:
      t2<-gsub(pattern = "ZT", replacement = "", as.data.frame(t2)[,1])
      #return(cbind(as.data.frame(t1)[,1], as.data.frame(t2)[,1]))
      return(cbind(as.data.frame(t1)[,1],t2))

         })
      t4<-rbind(t3[[1]], t3[[2]], t3[[3]], t3[[4]])
      return(t4)
        })
     c4<-rbind(c3[[1]], c3[[2]], c3[[3]], c3[[4]])
    return(c4)
}

v1<-phasecor_diffinteralltimepoints_pe_vectorphases(allreadsabove_difinter_timepoint, list_circpromsperphase_bed)

v2<-phasecor_diffinteralltimepoints_pe_vectorphases(allreadsbelow_difinter_timepoint,  list_circpromsperphase_bed)

#Heatmap
my_palette <- colorRampPalette(c("black","blue"))(n = 51)

#heatmap.2(as.matrix(circadian_tableinteractions[as.numeric(rownames(diff_score[[1]])),7:10]), 
#          trace = "none", density.info = "none", dendrogram="row",labRow =  "", 
#          labCol  = c("ZT0", "ZT6", "ZT12", "ZT18"), xlab="" ,
#          col=my_palette, Colv="NA", main="ZT0-ZT6", na.rm = T)

  ##Merge set of all differential interactions function

#there are only 61 circadian promoters shared in all the timepoints
heatmap.2(as.matrix(allreadsabove_difinter_timepoint[,7:10]), trace = "none", density.info = "none", dendrogram="row",labRow =  "", key.title = "", key.xlab = "Reads", labCol  = c("ZT0", "ZT6", "ZT12", "ZT18"), xlab="" ,col = viridis(50), breaks = seq(0,50, 1), Colv="NA", main= paste(">150kb n=", as.character(dim(allreadsabove_difinter_timepoint)[1]), sep = ""))

heatmap.2(as.matrix(allreadsbelow_difinter_timepoint[,7:10]), trace = "none", density.info = "none", dendrogram="row",labRow =  "", labCol  = c("ZT0", "ZT6", "ZT12", "ZT18"), xlab="" , Colv="NA", main=  paste("<150kb n=", as.character(dim(allreadsbelow_difinter_timepoint)[1]), sep = ""),  key.title = "", key.xlab = "Reads",col = viridis(50), breaks = seq(0,50, 1))
Heatmap_above150KbDI




#Count the number of differential interactions- using dual comparisons

##Positions in the list 2 corresponds to ZT0-ZT12 and index 5 to ZT6-ZT18
#---------Above 150Kb-------
ZT0_12_presentinZT12above150Kb<-circadian_tableinteractions[as.numeric(rownames(diffs_readcount_above150kb[[2]][[1]])),]
ZT0_12_presentinZT0above150Kb<-circadian_tableinteractions[as.numeric(rownames(diffs_readcount_above150kb[[2]][[2]])),]

ZT6_18_presentinZT18above150Kb<-circadian_tableinteractions[as.numeric(rownames(diffs_readcount_above150kb[[5]][[1]])),]
ZT6_18_presentinZT6above150Kb<-circadian_tableinteractions[as.numeric(rownames(diffs_readcount_above150kb[[5]][[2]])),]

#---------Below 150Kb-------
ZT0_12_presentinZT12below150Kb<-circadian_tableinteractions[as.numeric(rownames(diffs_readcount_below150kb[[2]][[1]])),]
ZT0_12_presentinZT0below150Kb<-circadian_tableinteractions[as.numeric(rownames(diffs_readcount_below150kb[[2]][[2]])),]

ZT6_18_presentinZT18below150Kb<-circadian_tableinteractions[as.numeric(rownames(diffs_readcount_below150kb[[5]][[1]])),]
ZT6_18_presentinZT6below150Kb<-circadian_tableinteractions[as.numeric(rownames(diffs_readcount_below150kb[[5]][[2]])),]



##Function to count the circproms that have differential interactions
countdifinter_circpromsperphase<-function(difinters, indexcircproms1, indexcircproms2){
  grangesconservedinalltimepoints_B<-GRanges(seqnames = difinters[,1], ranges = IRanges(difinters[,2], difinters[,3]), ZT0= difinters[,7], ZT6= difinters[,8], ZT12= difinters[,9], ZT18= difinters[,10])
  grangesconservedinalltimepoints_OE<-GRanges(seqnames = difinters[,4], ranges = IRanges(difinters[,5], difinters[,6]),ZT0= difinters[,7], ZT6= difinters[,8], ZT12= difinters[,9], ZT18= difinters[,10] )
  c1<-length(findOverlaps(grangesconservedinalltimepoints_B, list_circpromsperphase_bed[[indexcircproms1]]))
  c2<-length(findOverlaps(grangesconservedinalltimepoints_B, list_circpromsperphase_bed[[indexcircproms2]]))
  return(c(c1, c2))
}

#---------Above 150Kb-------
countdifinter_circpromsperphase(ZT0_12_presentinZT12above150Kb, 1, 3) #zt0 proms: 48 zt12 proms:131
countdifinter_circpromsperphase(ZT0_12_presentinZT0above150Kb, 1, 3) #zt0 proms: 56 zt12 proms: 96

countdifinter_circpromsperphase(ZT6_18_presentinZT18above150Kb, 2, 4) #zt6 proms:78 zt18 proms:37
countdifinter_circpromsperphase(ZT6_18_presentinZT6above150Kb, 2, 4)#zt6 proms: 104 zt18 proms: 55

#---------Below 150Kb-------
countdifinter_circpromsperphase(ZT0_12_presentinZT12below150Kb, 1, 3) #zt0 proms: 315 zt12 proms:811
countdifinter_circpromsperphase(ZT0_12_presentinZT0below150Kb, 1, 3) #zt0 proms: 284 zt12 proms: 434

countdifinter_circpromsperphase(ZT6_18_presentinZT18below150Kb, 2, 4) #zt6 proms:435 zt18 proms:315
countdifinter_circpromsperphase(ZT6_18_presentinZT6below150Kb, 2, 4)#zt6 proms: 355 zt18 proms: 198



##Function to count the eRNAs that have differential interactions: Doesnt work, not many eRNAs overlap with chic

countdifinter_ernasperphase<-function(difinters, indexcircproms1, indexcircproms2){
  grangesconservedinalltimepoints_B<-GRanges(seqnames = difinters[,1], ranges = IRanges(difinters[,2], difinters[,3]), ZT0= difinters[,7], ZT6= difinters[,8], ZT12= difinters[,9], ZT18= difinters[,10])
  grangesconservedinalltimepoints_OE<-GRanges(seqnames = difinters[,4], ranges = IRanges(difinters[,5], difinters[,6]),ZT0= difinters[,7], ZT6= difinters[,8], ZT12= difinters[,9], ZT18= difinters[,10] )
  c1<-length(findOverlaps(grangesconservedinalltimepoints_OE, list_ernasperphase_bed_mod[[indexcircproms1]]))
  c2<-length(findOverlaps(grangesconservedinalltimepoints_OE, list_ernasperphase_bed_mod[[indexcircproms2]]))
  return(c(c1, c2))
}

#---------Above 150Kb-------
countdifinter_ernasperphase(ZT0_12_presentinZT12above150Kb, 1, 3) #zt0 ernas: 4 zt12 ernas:1
countdifinter_ernasperphase(ZT0_12_presentinZT0above150Kb, 1, 3) #zt0 ernas: 3 zt12 ernas: 0

countdifinter_ernasperphase(ZT6_18_presentinZT18above150Kb, 2, 4) #zt6 ernas:0 zt18 ernas:13
countdifinter_ernasperphase(ZT6_18_presentinZT6above150Kb, 2, 4)#zt6 ernas: 1 zt18 ernas: 1

#---------Below 150Kb-------
countdifinter_ernasperphase(ZT0_12_presentinZT12below150Kb, 1, 3) #zt0 ernas: 57 zt12 ernas:10
countdifinter_ernasperphase(ZT0_12_presentinZT0below150Kb, 1, 3) #zt0 ernas: 21 zt12 ernas: 9

countdifinter_ernasperphase(ZT6_18_presentinZT18below150Kb, 2, 4) #zt6 ernas:9 zt18 ernas:43
countdifinter_ernasperphase(ZT6_18_presentinZT6below150Kb, 2, 4)#zt6 ernas: 13 zt18 ernas: 22


##Prom-eRNAs phase
difinter_circpromsernasperphase<-function(difinters, indexcircproms1, indexcircproms2){
  grangesconservedinalltimepoints_B<-GRanges(seqnames = difinters[,1], ranges = IRanges(difinters[,2], difinters[,3]), ZT0= difinters[,7], ZT6= difinters[,8], ZT12= difinters[,9], ZT18= difinters[,10])
  grangesconservedinalltimepoints_OE<-GRanges(seqnames = difinters[,4], ranges = IRanges(difinters[,5], difinters[,6]),ZT0= difinters[,7], ZT6= difinters[,8], ZT12= difinters[,9], ZT18= difinters[,10] )
      c3<-sapply(1:4, function(x){
      index<-findOverlaps(grangesconservedinalltimepoints_OE, list_ernasperphase_bed_mod[[x]])
      mcols(list_ernasperphase_bed_mod[[x]][subjectHits(index)])
      t3<-sapply(1:4, function(y){
        index2<-findOverlaps(grangesconservedinalltimepoints_B[queryHits(index)], list_circpromsperphase_bed[[y]])
      t1<-grangesconservedinalltimepoints_OE[[subjectHits(index)[queryHits(index2)]]]
      t2<-grangesconservedinalltimepoints_B[subjectHits(index2),]
      return(cbind(as.data.frame(t2), as.data.frame(t1)))

         })
      t4<-rbind(t3[[1]], t3[[2]], t3[[3]], t3[[4]])
      return(t4)
        })
  
})



#Count the number of interactions (complete ChIC per phase)- using circproms per phase that make dif inters


# 1) Starting from circproms per phase that make differential interactions, I count the number of interactions they make by viewing the normal set of CHiC by timepoints: Tables (x: PROM ZTX y: ChIC ZT0; differentials ChIC ZT0; enhancers; differentials with enhancers; enhancers phase?)


#Function to extract the circproms that make dif inters per phase
circpromsperphase_fromdifinters<-function(difinters,list_circpromsperphase_bed ){
  grangesconservedinalltimepoints_B<-GRanges(seqnames = difinters[,1], ranges = IRanges(difinters[,2], difinters[,3]), ZT0_difint= difinters[,7], ZT6_difint= difinters[,8], ZT12_difint= difinters[,9], ZT18_difint= difinters[,10])
  grangesconservedinalltimepoints_OE<-GRanges(seqnames = difinters[,4], ranges = IRanges(difinters[,5], difinters[,6]),ZT0_difint= difinters[,7], ZT6_difint= difinters[,8], ZT12_difint= difinters[,9], ZT18_difint= difinters[,10] )
  c3<-sapply(1:4, function(x){
    index<-subjectHits(findOverlaps(grangesconservedinalltimepoints_B, list_circpromsperphase_bed[[x]]))
    t<-list_circpromsperphase_bed[[x]][index]
    return(t)
    })
  return(c3)
}

circpromsabove<-circpromsperphase_fromdifinters(allreadsabove_difinter_timepoint, list_circpromsperphase_bed)
circpromsbelow<-circpromsperphase_fromdifinters(allreadsbelow_difinter_timepoint, list_circpromsperphase_bed)
  #Introns
circpromsINTRONSabove<-circpromsperphase_fromdifinters(allreadsabove_difinter_timepointINTRONS, list_circpromsperphaseINTRONS_bed)
circpromsINTRONSbelow<-circpromsperphase_fromdifinters(allreadsbelow_difinter_timepointINTRONS, list_circpromsperphaseINTRONS_bed)
  #Fang
circpromsFANGabove<-circpromsperphase_fromdifinters(allreadsabove_difinter_timepointFANG, list_circpromsperphaseFANG_bed_mod)
circpromsFANGbelow<-circpromsperphase_fromdifinters(allreadsbelow_difinter_timepointFANG, list_circpromsperphaseFANG_bed_mod)

#Merge circproms of each distance per pahse
circproms_makingdifinter<-sapply(1:4, function(x){
  unique(c(circpromsabove[[x]], circpromsbelow[[x]]))
})
  #Introns
circpromsINTRONS_makingdifinter<-sapply(1:4, function(x){
  unique(c(circpromsINTRONSabove[[x]], circpromsINTRONSbelow[[x]]))
})
  #FANG
circpromsFANG_makingdifinter<-sapply(1:4, function(x){
  unique(c(circpromsFANGabove[[x]], circpromsFANGbelow[[x]]))
})

#Open Chicago files 
ZT0<-read.csv("inputs_for_scripts/ZT0All_step2_washU_text_2.txt", sep="\t", header = F)
ZT6<-read.csv("inputs_for_scripts/ZT6All_step2_washU_text_2.txt", sep="\t", header = F)
ZT12<-read.csv("inputs_for_scripts/ZT12All_step2_washU_text_2.txt", sep="\t", header = F)
ZT18<-read.csv("inputs_for_scripts/ZT18All_step2_washU_text_2.txt", sep="\t", header = F)

#Overlap of circproms with chic and count numer of interactions per chicphase
countintersfromchics<-function(ZT0, ZT6, ZT12, ZT18, circproms_makingdifinter){
  countsallphases<-sapply(1:4, function(circpromphase){
       grangesconservedinalltimepoints_B<-GRanges(seqnames = ZT0[,1], ranges = IRanges(ZT0[,2], ZT0[,3]))
  grangesconservedinalltimepoints_OE<-GRanges(seqnames = ZT0[,4], ranges = IRanges(ZT0[,5], ZT0[,6]) )
  zt0chic<-countOverlaps(circproms_makingdifinter[[circpromphase]], grangesconservedinalltimepoints_B)
      grangesconservedinalltimepoints_B<-GRanges(seqnames = ZT6[,1], ranges = IRanges(ZT6[,2], ZT6[,3]))
  zt6chic<-countOverlaps(circproms_makingdifinter[[circpromphase]], grangesconservedinalltimepoints_B)  
        grangesconservedinalltimepoints_B<-GRanges(seqnames = ZT12[,1], ranges = IRanges(ZT12[,2], ZT12[,3]))
  zt12chic<-countOverlaps(circproms_makingdifinter[[circpromphase]], grangesconservedinalltimepoints_B) 
          grangesconservedinalltimepoints_B<-GRanges(seqnames = ZT18[,1], ranges = IRanges(ZT18[,2], ZT18[,3]))
  zt18chic<-countOverlaps(circproms_makingdifinter[[circpromphase]], grangesconservedinalltimepoints_B) 
  
    counts<-cbind(zt0chic, zt6chic, zt12chic, zt18chic)
    rownames(counts)<-paste(as.data.frame(circproms_makingdifinter[[circpromphase]])[,1], as.data.frame(circproms_makingdifinter[[circpromphase]])[,2], as.data.frame(circproms_makingdifinter[[circpromphase]])[,3], sep = "_")
    return(counts)
  })
  return(countsallphases)
}


#Apply function
countinter_circpromsperphasedifinters<-countintersfromchics(ZT0, ZT6, ZT12, ZT18, circproms_makingdifinter)
names(countinter_circpromsperphasedifinters)<-c("ZT0", "ZT6", "ZT12", "ZT18") 
  #Introns
countinter_circpromsINTRONSperphasedifinters<-countintersfromchics(ZT0, ZT6, ZT12, ZT18, circpromsINTRONS_makingdifinter)
names(countinter_circpromsINTRONSperphasedifinters)<-c("ZT0", "ZT6", "ZT12", "ZT18") 
  #FANG
countinter_circpromsFANGperphasedifinters<-countintersfromchics(ZT0, ZT6, ZT12, ZT18, circpromsFANG_makingdifinter)
names(countinter_circpromsFANGperphasedifinters)<-c("ZT0", "ZT6", "ZT12", "ZT18") 

#Export lists
write.table(countinter_circpromsperphasedifinters[[1]], "inputs_for_scripts/differential_interactions/MFM_RNAseq/tablecount_interactionscircpromsZT0dectectedindifinters.txt", quote = F, sep = "\t")
write.table(countinter_circpromsperphasedifinters[[2]], "inputs_for_scripts/differential_interactions/MFM_RNAseq/tablecount_interactionscircpromsZT6dectectedindifinters.txt", quote = F, sep = "\t")
write.table(countinter_circpromsperphasedifinters[[3]], "inputs_for_scripts/differential_interactions/MFM_RNAseq/tablecount_interactionscircpromsZT12dectectedindifinters.txt", quote = F, sep = "\t")
write.table(countinter_circpromsperphasedifinters[[4]], "inputs_for_scripts/differential_interactions/MFM_RNAseq/tablecount_interactionscircpromsZT18dectectedindifinters.txt", quote = F, sep = "\t")

#In which timepoint I can find the max and min per circprom;
##ZT0
#2) What proportion of ZT0 promoters have the maximum interactions to ZT0?  53/183
table(apply(countinter_circpromsperphasedifinters[[1]], 1, which.max)==1)
# 3) What proportion of ZT0 promoters have the least interactions with ZT12? 33/183
table(apply(countinter_circpromsperphasedifinters[[1]], 1, which.min)==3)
# 4) What proportion of ZT0 promoters have the maximum at 0 and the minimum at 12? 16/183
table((apply(countinter_circpromsperphasedifinters[[1]], 1, which.max)==1)+(apply(countinter_circpromsperphasedifinters[[1]], 1, which.min)==3))
##ZT6
# 2) What proportion of ZT6 promoters have the maximum interactions to ZT6? 81/283= 0.2862191
table(apply(countinter_circpromsperphasedifinters[[2]], 1, which.max)==2)
# 3) What proportion of ZT6 promoters have the least interactions to ZT18? 39/283=0.1378092
table(apply(countinter_circpromsperphasedifinters[[2]], 1, which.min)==4)
# 4) What proportion of ZT6 promoters have the maximum at 6 and the minimum at 18? 18/283=0.06360424
table((apply(countinter_circpromsperphasedifinters[[2]], 1, which.max)==2)+(apply(countinter_circpromsperphasedifinters[[2]], 1, which.min)==4))
##ZT12
# 2) What proportion of ZT12 promoters have the maximum interactions to ZT12? 110/364= 0.3021978
table(apply(countinter_circpromsperphasedifinters[[3]], 1, which.max)==3)
# 3) What proportion of ZT12 promoters have the least interactions to ZT0?  190/364=0.521978
table(apply(countinter_circpromsperphasedifinters[[3]], 1, which.min)==1)
# 4) What proportion of ZT12 promoters have the maximum at 12 and the minimum at 0? 66/364=0.1813187
table((apply(countinter_circpromsperphasedifinters[[3]], 1, which.max)==3)+(apply(countinter_circpromsperphasedifinters[[3]], 1, which.min)==1))
##ZT18
# 2) What proportion of ZT18 promoters have the maximum interactions to ZT18? 42/160= 0.2625
table(apply(countinter_circpromsperphasedifinters[[4]], 1, which.max)==4)
# 3) What proportion of ZT18 promoters have the least interactions to ZT6? 54/160= 0.3375
table(apply(countinter_circpromsperphasedifinters[[4]], 1, which.min)==2)
# 4) What proportion of ZT18 promoters have the maximum at 18 and the minimum at 6? 14/160=0.0875
table((apply(countinter_circpromsperphasedifinters[[4]], 1, which.max)==4)+(apply(countinter_circpromsperphasedifinters[[4]], 1, which.min)==2))


#General max and min table 
#count max and mins
countmaxandmins<-function(countinter_difinters){
  t<-sapply(1:4, function(phase){
  max<-sum(unlist(apply(countinter_difinters, 1, function(x){if(sum(x)>0){return(which.max(x)==phase)}})))
  min<-sum(unlist(apply(countinter_difinters, 1, function(x){if(sum(x)>0){return(which.min(x)==phase)}})))
  return(c(max, min))
  })

  rownames(t)<-c("Number of Max interactions","Number of Min interactions" )
colnames(t)<-c("ZT0CHiC", "ZT6CHiC", "ZT12CHiC", "ZT18CHiC") 
  return(t)
}

#count max and mins
countmaxandmins(countinter_circpromsperphasedifinters[[1]])
countmaxandmins(countinter_circpromsperphasedifinters[[2]])
countmaxandmins(countinter_circpromsperphasedifinters[[3]])
countmaxandmins(countinter_circpromsperphasedifinters[[4]])


########2020-01-26-INTRONS & FANG
    ##INTRONS
#count max and mins
countmaxandmins(countinter_circpromsINTRONSperphasedifinters[[1]])
countmaxandmins(countinter_circpromsINTRONSperphasedifinters[[2]])
countmaxandmins(countinter_circpromsINTRONSperphasedifinters[[3]])
countmaxandmins(countinter_circpromsINTRONSperphasedifinters[[4]])

    ##FANG
#count max and mins
countmaxandmins(countinter_circpromsFANGperphasedifinters[[1]])
countmaxandmins(countinter_circpromsFANGperphasedifinters[[2]])
countmaxandmins(countinter_circpromsFANGperphasedifinters[[3]])
countmaxandmins(countinter_circpromsFANGperphasedifinters[[4]])

#-----------------------------------------
#Count the enhancers, this functions counts the osceRNAs
countintersfromchicALLENH<-function(ZT0, ZT6, ZT12, ZT18, circproms_makingdifinter){
  countsallphases<-sapply(1:4, function(transcriptionalphase){
    #ZT0
    grangesconservedinalltimepoints_B<-GRanges(seqnames = ZT0[,1], ranges = IRanges(ZT0[,2], ZT0[,3]))
    grangesconservedinalltimepoints_OE<-GRanges(seqnames = ZT0[,4], ranges = IRanges(ZT0[,5], ZT0[,6]) )
    oetempo<-subjectHits(findOverlaps(circproms_makingdifinter[[transcriptionalphase]], grangesconservedinalltimepoints_B)) # retrieve the index OE of the Baits that overlap with circproms
    zt0chicenh<-sapply(1:4, function(ernasphase){queryHits(findOverlaps(grangesconservedinalltimepoints_OE[oetempo], list_ernasperphase_bed_mod[[ernasphase]]))}) #overlap of the OE with the eRNAs, retrieve the OE indexes
    zt0chicenh_retrievecircrpomsindex<-countOverlaps(circproms_makingdifinter[[transcriptionalphase]], grangesconservedinalltimepoints_B[oetempo[unlist(zt0chicenh)]]) #Again, overlap of the baits [index of OE that overlaped with erNAS] and the circrproms 
    
    #ZT6
    grangesconservedinalltimepoints_B<-GRanges(seqnames = ZT6[,1], ranges = IRanges(ZT6[,2], ZT6[,3]))
    grangesconservedinalltimepoints_OE<-GRanges(seqnames = ZT6[,4], ranges = IRanges(ZT6[,5], ZT6[,6]) )
    oetempo<-subjectHits(findOverlaps(circproms_makingdifinter[[transcriptionalphase]], grangesconservedinalltimepoints_B)) # retrieve the index OE of the Baits that overlap with circproms
    zt6chicenh<-sapply(1:4, function(ernasphase){queryHits(findOverlaps(grangesconservedinalltimepoints_OE[oetempo], list_ernasperphase_bed_mod[[ernasphase]]))}) #overlap of the OE with the eRNAs, retrieve the OE indexes
    zt6chicenh_retrievecircrpomsindex<-countOverlaps(circproms_makingdifinter[[transcriptionalphase]], grangesconservedinalltimepoints_B[oetempo[unlist(zt6chicenh)]]) #Again, overlap of the baits [index of OE that overlaped with erNAS] and the circrproms 
    
    #ZT12
    grangesconservedinalltimepoints_B<-GRanges(seqnames = ZT12[,1], ranges = IRanges(ZT12[,2], ZT12[,3]))
    grangesconservedinalltimepoints_OE<-GRanges(seqnames = ZT12[,4], ranges = IRanges(ZT12[,5], ZT12[,6]) )
    oetempo<-subjectHits(findOverlaps(circproms_makingdifinter[[transcriptionalphase]], grangesconservedinalltimepoints_B)) # retrieve the index OE of the Baits that overlap with circproms
    zt12chicenh<-sapply(1:4, function(ernasphase){queryHits(findOverlaps(grangesconservedinalltimepoints_OE[oetempo], list_ernasperphase_bed_mod[[ernasphase]]))}) #overlap of the OE with the eRNAs, retrieve the OE indexes
    zt12chicenh_retrievecircrpomsindex<-countOverlaps(circproms_makingdifinter[[transcriptionalphase]], grangesconservedinalltimepoints_B[oetempo[unlist(zt12chicenh)]]) #Again, overlap of the baits [index of OE that overlaped with erNAS] and the circrproms 
    
    
    #ZT18
    grangesconservedinalltimepoints_B<-GRanges(seqnames = ZT18[,1], ranges = IRanges(ZT18[,2], ZT18[,3]))
    grangesconservedinalltimepoints_OE<-GRanges(seqnames = ZT18[,4], ranges = IRanges(ZT18[,5], ZT18[,6]) )
    oetempo<-subjectHits(findOverlaps(circproms_makingdifinter[[transcriptionalphase]], grangesconservedinalltimepoints_B)) # retrieve the index OE of the Baits that overlap with circproms
    zt18chicenh<-sapply(1:4, function(ernasphase){queryHits(findOverlaps(grangesconservedinalltimepoints_OE[oetempo], list_ernasperphase_bed_mod[[ernasphase]]))}) #overlap of the OE with the eRNAs, retrieve the OE indexes
    zt18chicenh_retrievecircrpomsindex<-countOverlaps(circproms_makingdifinter[[transcriptionalphase]], grangesconservedinalltimepoints_B[oetempo[unlist(zt18chicenh)]]) #Again, overlap of the baits [index of OE that overlaped with erNAS] and the circrproms 
    
    counts<-cbind(zt0chicenh_retrievecircrpomsindex, zt6chicenh_retrievecircrpomsindex, zt12chicenh_retrievecircrpomsindex, zt18chicenh_retrievecircrpomsindex)
    colnames(counts)<-c("ZT0CHiC", "ZT6CHiC", "ZT12CHiC", "ZT18CHiC") 
    rownames(counts)<-paste(as.data.frame(circproms_makingdifinter[[transcriptionalphase]])[,1], as.data.frame(circproms_makingdifinter[[transcriptionalphase]])[,2], as.data.frame(circproms_makingdifinter[[transcriptionalphase]])[,3], sep = "_")
    return(counts)
  })
  return(countsallphases)
}

#Apply function; the result is a list (4) where each position is the transcriptional phase of the circproms, inside each list there is a data.frame with 4 cols that are the number of interactions with eRNAs accoridng to the chics and the rownames is the circproms genomic coordinates. 
countinter_eRNAsassociatedwithcircpromsperphasedifinters<-countintersfromchicALLENH(ZT0, ZT6, ZT12, ZT18, circproms_makingdifinter)
names(countinter_eRNAsassociatedwithcircpromsperphasedifinters)<-c("ZT0circproms", "ZT6circproms", "ZT12circproms", "ZT18circproms") 

#count max and mins
countmaxandmins(countinter_eRNAsassociatedwithcircpromsperphasedifinters[[1]])
countmaxandmins(countinter_eRNAsassociatedwithcircpromsperphasedifinters[[2]])
countmaxandmins(countinter_eRNAsassociatedwithcircpromsperphasedifinters[[3]])
countmaxandmins(countinter_eRNAsassociatedwithcircpromsperphasedifinters[[4]])


#How many interactions are made in total by the circproms
ZT0circproms that make differential interactions 183; from these how many interactions are detected in the chicZT0 2492, chicZT6 2531, chicZT12 2607, chicZT18 2699. 
ZT6circproms that make differential interactions 238
ZT12circproms that make differential interactions 364
ZT18circproms that make differential interactions 160; 

#-----------------------------------------
#2018-08-20 Count again the number of enhancers but using the total list
#1)Open total eRNAs
total_ernas<-read.csv("/Users/andoku01/Dropbox/Lab_IFC_Mayra/eRNAs/eRNAs_de_novo_oscillating.txt", header = T, sep = "\t")
total_ernas<-bedfile(total_ernas, columnnames)

#2) Count the enhancers, this functions counts the total eRNAs

countintersfromchicALLENH_TOTAL<-function(ZT0, ZT6, ZT12, ZT18, circproms_makingdifinter){
  countsallphases<-sapply(1:4, function(transcriptionalphase){
    #ZT0
    grangesconservedinalltimepoints_B<-GRanges(seqnames = ZT0[,1], ranges = IRanges(ZT0[,2], ZT0[,3]))
    grangesconservedinalltimepoints_OE<-GRanges(seqnames = ZT0[,4], ranges = IRanges(ZT0[,5], ZT0[,6]) )
    oetempo<-subjectHits(findOverlaps(circproms_makingdifinter[[transcriptionalphase]], grangesconservedinalltimepoints_B)) # retrieve the index OE of the Baits that overlap with circproms
    zt0chicenh<-queryHits(findOverlaps(grangesconservedinalltimepoints_OE[oetempo], total_ernas)) #overlap of the OE with the eRNAs, retrieve the OE indexes
    zt0chicenh_retrievecircrpomsindex<-countOverlaps(circproms_makingdifinter[[transcriptionalphase]], grangesconservedinalltimepoints_B[oetempo[zt0chicenh]]) #Again, overlap of the baits [index of OE that overlaped with erNAS] and the circrproms 
    
    #ZT6
    grangesconservedinalltimepoints_B<-GRanges(seqnames = ZT6[,1], ranges = IRanges(ZT6[,2], ZT6[,3]))
    grangesconservedinalltimepoints_OE<-GRanges(seqnames = ZT6[,4], ranges = IRanges(ZT6[,5], ZT6[,6]) )
    oetempo<-subjectHits(findOverlaps(circproms_makingdifinter[[transcriptionalphase]], grangesconservedinalltimepoints_B)) # retrieve the index OE of the Baits that overlap with circproms
    zt6chicenh<-queryHits(findOverlaps(grangesconservedinalltimepoints_OE[oetempo], total_ernas)) #overlap of the OE with the eRNAs, retrieve the OE indexes
    zt6chicenh_retrievecircrpomsindex<-countOverlaps(circproms_makingdifinter[[transcriptionalphase]], grangesconservedinalltimepoints_B[oetempo[zt6chicenh]]) #Again, overlap of the baits [index of OE that overlaped with erNAS] and the circrproms 
    
    #ZT12
    grangesconservedinalltimepoints_B<-GRanges(seqnames = ZT12[,1], ranges = IRanges(ZT12[,2], ZT12[,3]))
    grangesconservedinalltimepoints_OE<-GRanges(seqnames = ZT12[,4], ranges = IRanges(ZT12[,5], ZT12[,6]) )
    oetempo<-subjectHits(findOverlaps(circproms_makingdifinter[[transcriptionalphase]], grangesconservedinalltimepoints_B)) # retrieve the index OE of the Baits that overlap with circproms
    zt12chicenh<-queryHits(findOverlaps(grangesconservedinalltimepoints_OE[oetempo], total_ernas)) #overlap of the OE with the eRNAs, retrieve the OE indexes
    zt12chicenh_retrievecircrpomsindex<-countOverlaps(circproms_makingdifinter[[transcriptionalphase]], grangesconservedinalltimepoints_B[oetempo[zt12chicenh]]) #Again, overlap of the baits [index of OE that overlaped with erNAS] and the circrproms 
    
    
    #ZT18
    grangesconservedinalltimepoints_B<-GRanges(seqnames = ZT18[,1], ranges = IRanges(ZT18[,2], ZT18[,3]))
    grangesconservedinalltimepoints_OE<-GRanges(seqnames = ZT18[,4], ranges = IRanges(ZT18[,5], ZT18[,6]) )
    oetempo<-subjectHits(findOverlaps(circproms_makingdifinter[[transcriptionalphase]], grangesconservedinalltimepoints_B)) # retrieve the index OE of the Baits that overlap with circproms
    zt18chicenh<-queryHits(findOverlaps(grangesconservedinalltimepoints_OE[oetempo], total_ernas)) #overlap of the OE with the eRNAs, retrieve the OE indexes
    zt18chicenh_retrievecircrpomsindex<-countOverlaps(circproms_makingdifinter[[transcriptionalphase]], grangesconservedinalltimepoints_B[oetempo[zt18chicenh]]) #Again, overlap of the baits [index of OE that overlaped with erNAS] and the circrproms 
    
    counts<-cbind(zt0chicenh_retrievecircrpomsindex, zt6chicenh_retrievecircrpomsindex, zt12chicenh_retrievecircrpomsindex, zt18chicenh_retrievecircrpomsindex)
    colnames(counts)<-c("ZT0CHiC", "ZT6CHiC", "ZT12CHiC", "ZT18CHiC") 
    rownames(counts)<-paste(as.data.frame(circproms_makingdifinter[[transcriptionalphase]])[,1], as.data.frame(circproms_makingdifinter[[transcriptionalphase]])[,2], as.data.frame(circproms_makingdifinter[[transcriptionalphase]])[,3], sep = "_")
    return(counts)
  })
  return(countsallphases)
}

#Apply function; the result is a list (4) where each position is the transcriptional phase of the circproms, inside each list there is a data.frame with 4 cols that are the number of interactions with eRNAs accoridng to the chics and the rownames is the circproms genomic coordinates. 
countinter_TOTALeRNAsassociatedwithcircpromsperphasedifinters<-countintersfromchicALLENH_TOTAL(ZT0, ZT6, ZT12, ZT18, circproms_makingdifinter)
names(countinter_TOTALeRNAsassociatedwithcircpromsperphasedifinters)<-c("ZT0circproms", "ZT6circproms", "ZT12circproms", "ZT18circproms") 

#count max and mins
countmaxandmins(countinter_TOTALeRNAsassociatedwithcircpromsperphasedifinters[[1]])
countmaxandmins(countinter_TOTALeRNAsassociatedwithcircpromsperphasedifinters[[2]])
countmaxandmins(countinter_TOTALeRNAsassociatedwithcircpromsperphasedifinters[[3]])
countmaxandmins(countinter_TOTALeRNAsassociatedwithcircpromsperphasedifinters[[4]])


#How many interactions are made in total by the circproms
ZT0circproms that make differential interactions 183; from these how many interactions are detected in the chicZT0 2492, chicZT6 2531, chicZT12 2607, chicZT18 2699. 
ZT6circproms that make differential interactions 238
ZT12circproms that make differential interactions 364
ZT18circproms that make differential interactions 160; 

#-----------------------------------------
#2018-08-20 Filter only the core clock genes and report all the min/max counts :)
#1) Open core clock genes; #25 gene list and divide per phases
coreclock<-read.csv("core_clock_genes/Core_clock_short_fragments.txt", header = T, sep = "\t")
coreclock[,1]<-paste("chr", coreclock[,1], sep = "")
coreclock<-makeGRangesFromDataFrame(coreclock, keep.extra.columns = T)

#Divide per phases; only 17 are in our circproms ? 
coreclock<-sapply(1:4, function(circpromsphase){
  t<-coreclock[queryHits(findOverlaps(query = coreclock,subject = list_circpromsperphase_bed[[circpromsphase]] , type = "within"))]
t$timepoint<-as.data.frame(mcols(list_circpromsperphase_bed[[circpromsphase]][subjectHits(findOverlaps(query = coreclock,subject = list_circpromsperphase_bed[[circpromsphase]] , type = "within"))]))
return(unique(t))
})


#2)#Function to extract the CORE-CLOCK circproms that make dif inters per phase
CCcircpromsabove<-circpromsperphase_fromdifinters(allreadsabove_difinter_timepoint, coreclock)
CCcircpromsbelow<-circpromsperphase_fromdifinters(allreadsbelow_difinter_timepoint, coreclock)

#3)Merge circproms of each distance per pahse
CCcircproms_makingdifinter<-sapply(1:4, function(x){
  unique(c(CCcircpromsabove[[x]], CCcircpromsbelow[[x]]))
})

#4) Count diff inters according to chic
#Apply function
countinter_CCcircpromsperphasedifinters<-countintersfromchics(ZT0, ZT6, ZT12, ZT18, CCcircproms_makingdifinter)
names(countinter_CCcircpromsperphasedifinters)<-c("ZT0", "ZT6", "ZT12", "ZT18") 
#Export lists
write.table(countinter_CCcircpromsperphasedifinters[[1]], "inputs_for_scripts/differential_interactions/MFM_RNAseq/tablecount_interactionsCORECLOCKcircpromsZT0dectectedindifinters.txt", quote = F, sep = "\t")
write.table(countinter_CCcircpromsperphasedifinters[[2]], "inputs_for_scripts/differential_interactions/MFM_RNAseq/tablecount_interactionsCORECLOCKcircpromsZT6dectectedindifinters.txt", quote = F, sep = "\t")
write.table(countinter_CCcircpromsperphasedifinters[[3]], "inputs_for_scripts/differential_interactions/MFM_RNAseq/tablecount_interactionsCORECLOCKcircpromsZT12dectectedindifinters.txt", quote = F, sep = "\t")
write.table(countinter_CCcircpromsperphasedifinters[[4]], "inputs_for_scripts/differential_interactions/MFM_RNAseq/tablecount_interactionsCORECLOCKcircpromsZT18dectectedindifinters.txt", quote = F, sep = "\t")

#count max and mins
countmaxandmins(countinter_CCcircpromsperphasedifinters[[1]])
countmaxandmins(countinter_CCcircpromsperphasedifinters[[2]])
countmaxandmins(countinter_CCcircpromsperphasedifinters[[3]])
countmaxandmins(countinter_CCcircpromsperphasedifinters[[4]])

#5) Count diff inters according to chic with total and osc eRNAs
#----total
#Apply function; the result is a list (4) where each position is the transcriptional phase of the circproms, inside each list there is a data.frame with 4 cols that are the number of interactions with eRNAs accoridng to the chics and the rownames is the circproms genomic coordinates. 
countinter_TOTALeRNAsassociatedwithCCcircpromsperphasedifinters<-countintersfromchicALLENH_TOTAL(ZT0, ZT6, ZT12, ZT18, CCcircproms_makingdifinter)
names(countinter_TOTALeRNAsassociatedwithCCcircpromsperphasedifinters)<-c("ZT0circproms", "ZT6circproms", "ZT12circproms", "ZT18circproms") 

#count max and mins
countmaxandmins(countinter_TOTALeRNAsassociatedwithCCcircpromsperphasedifinters[[1]])
countmaxandmins(countinter_TOTALeRNAsassociatedwithCCcircpromsperphasedifinters[[2]])
countmaxandmins(countinter_TOTALeRNAsassociatedwithCCcircpromsperphasedifinters[[3]])
countmaxandmins(countinter_TOTALeRNAsassociatedwithCCcircpromsperphasedifinters[[4]])

#----osceRNAs
#Apply function; the result is a list (4) where each position is the transcriptional phase of the circproms, inside each list there is a data.frame with 4 cols that are the number of interactions with eRNAs accoridng to the chics and the rownames is the circproms genomic coordinates. 
countinter_eRNAsassociatedwithCCcircpromsperphasedifinters<-countintersfromchicALLENH(ZT0, ZT6, ZT12, ZT18, CCcircproms_makingdifinter)
names(countinter_eRNAsassociatedwithCCcircpromsperphasedifinters)<-c("ZT0circproms", "ZT6circproms", "ZT12circproms", "ZT18circproms") 

#count max and mins
countmaxandmins(countinter_eRNAsassociatedwithCCcircpromsperphasedifinters[[1]])
countmaxandmins(countinter_eRNAsassociatedwithCCcircpromsperphasedifinters[[2]])
countmaxandmins(countinter_eRNAsassociatedwithCCcircpromsperphasedifinters[[3]])
countmaxandmins(countinter_eRNAsassociatedwithCCcircpromsperphasedifinters[[4]])

#2018-08-20 Identificar en la table de counts el circproms en donde el máximo transcripcional y el máximo de interacciones coinciden así como identificar el set en dónde el máximo transcripcional coincide con el máximo de interacciones y el mínimo transcripcional coincide con el mínimo de interacciones.
#Modification of countmaxandmins
countmaxandmins_addcolumstotablecounts<-function(countinter_difinters, maxtimepoint, mintimepoint){
  max<-apply(countinter_difinters, 1, function(x){if(sum(x)>0){return(which.max(x)==maxtimepoint)}})
  t<-cbind(countinter_difinters,"maxinter_in_maxtranscriptonalphase"=max)
  min<-apply(countinter_difinters, 1, function(x){if(sum(x)>0){return(which.min(x)==mintimepoint)}})
  w<-cbind(t, "mininter_in_mintranscriptonalphase"=min)
  return(w)

#  rownames(t)<-c("Number of Max interactions","Number of Min interactions" )
#colnames(t)<-c("ZT0CHiC", "ZT6CHiC", "ZT12CHiC", "ZT18CHiC") 
 # return(t)
}


countinter_circpromsperphasedifinters[[1]]<-countmaxandmins_addcolumstotablecounts(countinter_circpromsperphasedifinters[[1]], maxtimepoint = 1, mintimepoint = 3)
countinter_circpromsperphasedifinters[[2]]<-countmaxandmins_addcolumstotablecounts(countinter_circpromsperphasedifinters[[2]], maxtimepoint = 2, mintimepoint = 4)
countinter_circpromsperphasedifinters[[3]]<-countmaxandmins_addcolumstotablecounts(countinter_circpromsperphasedifinters[[3]], maxtimepoint = 3, mintimepoint = 1)
countinter_circpromsperphasedifinters[[4]]<-countmaxandmins_addcolumstotablecounts(countinter_circpromsperphasedifinters[[4]], maxtimepoint = 4, mintimepoint = 2)

#Export lists again and add gene name

#------------------------------------------------------
#2019-01-16 Report number using introns

circpromsINTRONSabove<-circpromsperphase_fromdifinters(allreadsabove_difinter_timepoint, list_circpromsperphaseINTRONS_bed)
circpromsINTRONSbelow<-circpromsperphase_fromdifinters(allreadsbelow_difinter_timepoint, list_circpromsperphaseINTRONS_bed)

#Merge circproms of each distance per pahse
circpromsINTRONS_makingdifinter<-sapply(1:4, function(x){
  unique(c(circpromsINTRONSabove[[x]], circpromsINTRONSbelow[[x]]))
})


#
#Apply function
countinter_circpromsINTRONSperphasedifinters<-countintersfromchics(ZT0, ZT6, ZT12, ZT18, circpromsINTRONS_makingdifinter)
names(countinter_circpromsINTRONSperphasedifinters)<-c("ZT0", "ZT6", "ZT12", "ZT18") 
#Export lists
write.table(countinter_circpromsINTRONSperphasedifinters[[1]], "inputs_for_scripts/differential_interactions/MFM_RNAseq/tablecount_interactionscircpromsINTRONSZT0dectectedindifinters.txt", quote = F, sep = "\t")
write.table(countinter_circpromsINTRONSperphasedifinters[[2]], "inputs_for_scripts/differential_interactions/MFM_RNAseq/tablecount_interactionscircpromsINTRONSZT6dectectedindifinters.txt", quote = F, sep = "\t")
write.table(countinter_circpromsINTRONSperphasedifinters[[3]], "inputs_for_scripts/differential_interactions/MFM_RNAseq/tablecount_interactionscircpromsZT12INTRONSdectectedindifinters.txt", quote = F, sep = "\t")
write.table(countinter_circpromsINTRONSperphasedifinters[[4]], "inputs_for_scripts/differential_interactions/MFM_RNAseq/tablecount_interactionscircpromsZT18INTRONSdectectedindifinters.txt", quote = F, sep = "\t")

#In which timepoint I can find the max and min per circprom;
##ZT0
#2) What proportion of ZT0 promoters have the maximum interactions to ZT0?  16/52
table(apply(countinter_circpromsINTRONSperphasedifinters[[1]], 1, which.max)==1)
# 3) What proportion of ZT0 promoters have the least interactions with ZT12? 7/52
table(apply(countinter_circpromsINTRONSperphasedifinters[[1]], 1, which.min)==3)
# 4) What proportion of ZT0 promoters have the maximum at 0 and the minimum at 12? 4/52
table((apply(countinter_circpromsINTRONSperphasedifinters[[1]], 1, which.max)==1)+(apply(countinter_circpromsINTRONSperphasedifinters[[1]], 1, which.min)==3))
##ZT6
# 2) What proportion of ZT6 promoters have the maximum interactions to ZT6? 7/22
table(apply(countinter_circpromsINTRONSperphasedifinters[[2]], 1, which.max)==2)
# 3) What proportion of ZT6 promoters have the least interactions to ZT18? 5/22
table(apply(countinter_circpromsINTRONSperphasedifinters[[2]], 1, which.min)==4)
# 4) What proportion of ZT6 promoters have the maximum at 6 and the minimum at 18? 3/22
table((apply(countinter_circpromsINTRONSperphasedifinters[[2]], 1, which.max)==2)+(apply(countinter_circpromsINTRONSperphasedifinters[[2]], 1, which.min)==4))
##ZT12
# 2) What proportion of ZT12 promoters have the maximum interactions to ZT12? 40/101
table(apply(countinter_circpromsINTRONSperphasedifinters[[3]], 1, which.max)==3)
# 3) What proportion of ZT12 promoters have the least interactions to ZT0? 50/101
table(apply(countinter_circpromsINTRONSperphasedifinters[[3]], 1, which.min)==1)
# 4) What proportion of ZT12 promoters have the maximum at 12 and the minimum at 0? 24/101
table((apply(countinter_circpromsINTRONSperphasedifinters[[3]], 1, which.max)==3)+(apply(countinter_circpromsINTRONSperphasedifinters[[3]], 1, which.min)==1))
##ZT18
# 2) What proportion of ZT18 promoters have the maximum interactions to ZT18? 6/43
table(apply(countinter_circpromsINTRONSperphasedifinters[[4]], 1, which.max)==4)
# 3) What proportion of ZT18 promoters have the least interactions to ZT6? 18/43
table(apply(countinter_circpromsINTRONSperphasedifinters[[4]], 1, which.min)==2)
# 4) What proportion of ZT18 promoters have the maximum at 18 and the minimum at 6? 3/43
table((apply(countinter_circpromsINTRONSperphasedifinters[[4]], 1, which.max)==4)+(apply(countinter_circpromsINTRONSperphasedifinters[[4]], 1, which.min)==2))

#General max and min table 

#count max and mins
countmaxandmins(countinter_circpromsINTRONSperphasedifinters[[1]])
countmaxandmins(countinter_circpromsINTRONSperphasedifinters[[2]])
countmaxandmins(countinter_circpromsINTRONSperphasedifinters[[3]])
countmaxandmins(countinter_circpromsINTRONSperphasedifinters[[4]])

#Count interactions with total eRNAs

#Apply function; the result is a list (4) where each position is the transcriptional phase of the circproms, inside each list there is a data.frame with 4 cols that are the number of interactions with eRNAs accoridng to the chics and the rownames is the circproms genomic coordinates. 
countinter_TOTALeRNAsassociatedwithcircpromsINTRONSperphasedifinters<-countintersfromchicALLENH_TOTAL(ZT0, ZT6, ZT12, ZT18, circpromsINTRONS_makingdifinter)
names(countinter_TOTALeRNAsassociatedwithcircpromsINTRONSperphasedifinters)<-c("ZT0circproms", "ZT6circproms", "ZT12circproms", "ZT18circproms") 

#count max and mins
countmaxandmins(countinter_TOTALeRNAsassociatedwithcircpromsINTRONSperphasedifinters[[1]])
countmaxandmins(countinter_TOTALeRNAsassociatedwithcircpromsINTRONSperphasedifinters[[2]])
countmaxandmins(countinter_TOTALeRNAsassociatedwithcircpromsINTRONSperphasedifinters[[3]])
countmaxandmins(countinter_TOTALeRNAsassociatedwithcircpromsINTRONSperphasedifinters[[4]])

#Count the enhancers, this functions counts the osceRNAs

#Apply function; the result is a list (4) where each position is the transcriptional phase of the circproms, inside each list there is a data.frame with 4 cols that are the number of interactions with eRNAs accoridng to the chics and the rownames is the circproms genomic coordinates. 
countinter_eRNAsassociatedwithcircpromsINTRONSperphasedifinters<-countintersfromchicALLENH(ZT0, ZT6, ZT12, ZT18, circpromsINTRONS_makingdifinter)
names(countinter_eRNAsassociatedwithcircpromsINTRONSperphasedifinters)<-c("ZT0circproms", "ZT6circproms", "ZT12circproms", "ZT18circproms") 

#count max and mins
countmaxandmins(countinter_eRNAsassociatedwithcircpromsINTRONSperphasedifinters[[1]])
countmaxandmins(countinter_eRNAsassociatedwithcircpromsINTRONSperphasedifinters[[2]])
countmaxandmins(countinter_eRNAsassociatedwithcircpromsINTRONSperphasedifinters[[3]])
countmaxandmins(countinter_eRNAsassociatedwithcircpromsINTRONSperphasedifinters[[4]])






#Skip
#9.2 Export list of interactions, circprom-oscenh


#9.2 Export list of interactions, circprom-oscenh


phasecor_diffinteralltimepoints_pe_positionspairs<-function(intertableall, list_circproms){
  grangesconservedinalltimepoints_B<-GRanges(seqnames = intertableall[,1], ranges = IRanges(intertableall[,2], intertableall[,3]))
  
  grangesconservedinalltimepoints_OE<-GRanges(seqnames = intertableall[,4], ranges = IRanges(intertableall[,5], intertableall[,6]))
      
    c3<-sapply(1:4, function(x){
      index<-findOverlaps(grangesconservedinalltimepoints_OE, list_ernasperphase_bed_mod[[x]])
      mcols(list_ernasperphase_bed_mod[[x]][subjectHits(index)])
      t3<-sapply(1:4, function(y){
        index2<-findOverlaps(grangesconservedinalltimepoints_B[queryHits(index)], list_circproms[[y]])
      t1<-grangesconservedinalltimepoints_OE[[subjectHits(index)[queryHits(index2)]],]
      t2<-grangesconservedinalltimepoints_B[subjectHits(index2),]
      return(cbind(as.data.frame(t2), as.data.frame(t1)))

         })
      t4<-rbind(t3[[1]], t3[[2]], t3[[3]], t3[[4]])
      return(t4)
        })
     c4<-rbind(c3[[1]], c3[[2]], c3[[3]], c3[[4]])
    return(c4)
}


phasecor_diffinteralltimepoints_pe_positionspairs<-function(intertableall){ grangesconservedinalltimepoints_B<-GRanges(seqnames = intertableall[,1], ranges = IRanges(intertableall[,2], intertableall[,3]))
  
  grangesconservedinalltimepoints_OE<-GRanges(seqnames = intertableall[,4], ranges = IRanges(intertableall[,5], intertableall[,6]))
      
    c3<-sapply(1:4, function(x){
      index<-findOverlaps(grangesconservedinalltimepoints_OE, list_ernasperphase_bed_mod[[x]])
      mcols(list_ernasperphase_bed_mod[[x]][subjectHits(index)])
      index2<-findOverlaps(grangesconservedinalltimepoints_B[queryHits(index)], list_circpromsperphase_bed[[x]])
      t1<-list_ernasperphase_bed_mod[[x]][subjectHits(index)[queryHits(index2)]]
      t2<-list_circpromsperphase_bed[[x]][subjectHits(index2)]
      return(cbind(as.character(levels(as.data.frame(t2)[,1]))[as.data.frame(t2)[,1]], as.data.frame(t2)[,2:3],as.character(levels(as.data.frame(t1)[,1]))[as.data.frame(t1)[,1]], as.data.frame(t1)[,2:3]))
      })
    return(c3)
}
v1<-phasecor_diffinteralltimepoints_pe_positionspairs(allreadsabove_difinter_timepoint, list_circpromsperphaseFANG_bed_mod)
v2<-phasecor_diffinteralltimepoints_pe_positionspairs(allreadsbelow_difinter_timepoint)
#1-3 only circproms, if 4-6 cols they would be th ernas
v1<-rbind(cbind(as.character(levels(v1[,1][[1]]))[v1[,1][[1]]], v1[,1][[2]], v1[,1][[3]],as.character(levels(v1[,1][[4]]))[v1[,1][[4]]], v1[,1][[5]], v1[,1][[6]]),cbind(as.character(levels(v1[,2][[1]]))[v1[,2][[1]]], v1[,2][[2]], v1[,2][[3]],as.character(levels(v1[,2][[4]]))[v1[,2][[4]]], v1[,2][[5]], v1[,2][[6]]), cbind(as.character(levels(v1[,3][[1]]))[v1[,3][[1]]], v1[,3][[2]], v1[,3][[3]], as.character(levels(v1[,3][[4]]))[v1[,3][[4]]], v1[,3][[5]], v1[,3][[6]]),cbind(as.character(levels(v1[,4][[1]]))[v1[,4][[1]]], v1[,4][[2]],v1[,4][[3]],as.character(levels(v1[,4][[4]]))[v1[,4][[4]]], v1[,4][[5]], v1[,4][[6]]))

v2<-rbind(cbind(as.character(levels(v2[,1][[1]]))[v2[,1][[1]]], v2[,1][[2]], v2[,1][[3]],as.character(levels(v2[,1][[4]]))[v2[,1][[4]]], v2[,1][[5]], v2[,1][[6]]),cbind(as.character(levels(v2[,2][[1]]))[v2[,2][[1]]], v2[,2][[2]], v2[,2][[3]],as.character(levels(v2[,2][[4]]))[v2[,2][[4]]], v2[,2][[5]], v2[,2][[6]]), cbind(as.character(levels(v2[,3][[1]]))[v2[,3][[1]]], v2[,3][[2]], v2[,3][[3]], as.character(levels(v2[,3][[4]]))[v2[,3][[4]]], v2[,3][[5]], v2[,3][[6]]),cbind(as.character(levels(v2[,4][[1]]))[v2[,4][[1]]], v2[,4][[2]],v2[,4][[3]],as.character(levels(v2[,4][[4]]))[v2[,4][[4]]], v2[,4][[5]], v2[,4][[6]]))

setwd("/hdisk7/mandok/")
write.table(unique(v1),"inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_promenhdiffinter_above150kb.bedpe",quote=F,sep="\t",row.names=F,col.names=F) 
write.table(unique(v2),"inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_promenhdiffinter_below150kb.bedpe",quote=F,sep="\t",row.names=F,col.names=F) 


#9.2.1 Export list of interactions, circprom-circprom 2019-01-14


#9.2.1 Export list of interactions, circprom-circprom 2019-01-14


phasecor_diffinteralltimepoints_pp_positionspairs<-function(intertableall, list_circproms){
  grangesconservedinalltimepoints_B<-GRanges(seqnames = intertableall[,1], ranges = IRanges(intertableall[,2], intertableall[,3]))
  
  grangesconservedinalltimepoints_OE<-GRanges(seqnames = intertableall[,4], ranges = IRanges(intertableall[,5], intertableall[,6]))
      
    c3<-sapply(1:4, function(x){
      index<-findOverlaps(grangesconservedinalltimepoints_OE, list_circpromsperphase_bed[[x]])
      mcols(list_circpromsperphase_bed[[x]][subjectHits(index)])
      t3<-sapply(1:4, function(y){
        index2<-findOverlaps(grangesconservedinalltimepoints_B[queryHits(index)], list_circpromsperphase_bed[[y]])
      t1<-grangesconservedinalltimepoints_OE[[subjectHits(index)[queryHits(index2)]],]
      t2<-grangesconservedinalltimepoints_B[subjectHits(index2),]
      return(cbind(as.data.frame(t2), as.data.frame(t1)))

         })
      t4<-rbind(t3[[1]], t3[[2]], t3[[3]], t3[[4]])
      return(t4)
        })
     c4<-rbind(c3[[1]], c3[[2]], c3[[3]], c3[[4]])
    return(c4)
}


phasecor_diffinteralltimepoints_pp_positionspairs<-function(intertableall){ grangesconservedinalltimepoints_B<-GRanges(seqnames = intertableall[,1], ranges = IRanges(intertableall[,2], intertableall[,3]))
  
  grangesconservedinalltimepoints_OE<-GRanges(seqnames = intertableall[,4], ranges = IRanges(intertableall[,5], intertableall[,6]))
      
    c3<-sapply(1:4, function(x){
      index<-findOverlaps(grangesconservedinalltimepoints_OE, list_ernasperphase_bed_mod[[x]])
      mcols(list_ernasperphase_bed_mod[[x]][subjectHits(index)])
      index2<-findOverlaps(grangesconservedinalltimepoints_B[queryHits(index)], list_circpromsperphase_bed[[x]])
      t1<-list_ernasperphase_bed_mod[[x]][subjectHits(index)[queryHits(index2)]]
      t2<-list_circpromsperphase_bed[[x]][subjectHits(index2)]
      return(cbind(as.character(levels(as.data.frame(t2)[,1]))[as.data.frame(t2)[,1]], as.data.frame(t2)[,2:3],as.character(levels(as.data.frame(t1)[,1]))[as.data.frame(t1)[,1]], as.data.frame(t1)[,2:3]))
      })
    return(c3)
}
v1<-phasecor_diffinteralltimepoints_pe_positionspairs(allreadsabove_difinter_timepoint, list_circpromsperphaseFANG_bed_mod)
v2<-phasecor_diffinteralltimepoints_pe_positionspairs(allreadsbelow_difinter_timepoint)
#1-3 only circproms, if 4-6 cols they would be th ernas
v1<-rbind(cbind(as.character(levels(v1[,1][[1]]))[v1[,1][[1]]], v1[,1][[2]], v1[,1][[3]],as.character(levels(v1[,1][[4]]))[v1[,1][[4]]], v1[,1][[5]], v1[,1][[6]]),cbind(as.character(levels(v1[,2][[1]]))[v1[,2][[1]]], v1[,2][[2]], v1[,2][[3]],as.character(levels(v1[,2][[4]]))[v1[,2][[4]]], v1[,2][[5]], v1[,2][[6]]), cbind(as.character(levels(v1[,3][[1]]))[v1[,3][[1]]], v1[,3][[2]], v1[,3][[3]], as.character(levels(v1[,3][[4]]))[v1[,3][[4]]], v1[,3][[5]], v1[,3][[6]]),cbind(as.character(levels(v1[,4][[1]]))[v1[,4][[1]]], v1[,4][[2]],v1[,4][[3]],as.character(levels(v1[,4][[4]]))[v1[,4][[4]]], v1[,4][[5]], v1[,4][[6]]))

v2<-rbind(cbind(as.character(levels(v2[,1][[1]]))[v2[,1][[1]]], v2[,1][[2]], v2[,1][[3]],as.character(levels(v2[,1][[4]]))[v2[,1][[4]]], v2[,1][[5]], v2[,1][[6]]),cbind(as.character(levels(v2[,2][[1]]))[v2[,2][[1]]], v2[,2][[2]], v2[,2][[3]],as.character(levels(v2[,2][[4]]))[v2[,2][[4]]], v2[,2][[5]], v2[,2][[6]]), cbind(as.character(levels(v2[,3][[1]]))[v2[,3][[1]]], v2[,3][[2]], v2[,3][[3]], as.character(levels(v2[,3][[4]]))[v2[,3][[4]]], v2[,3][[5]], v2[,3][[6]]),cbind(as.character(levels(v2[,4][[1]]))[v2[,4][[1]]], v2[,4][[2]],v2[,4][[3]],as.character(levels(v2[,4][[4]]))[v2[,4][[4]]], v2[,4][[5]], v2[,4][[6]]))


write.table(unique(v1),"inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_promenhdiffinter_above150kb.bedpe",quote=F,sep="\t",row.names=F,col.names=F) 
write.table(unique(v2),"inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_promenhdiffinter_below150kb.bedpe",quote=F,sep="\t",row.names=F,col.names=F) 


# 9.3 Correlation indicating the timepoint where the interactions is stronger


#9.3 Correlation indicating the timepoint where the interactions is stronger

merge_uniquedifinter_updowntimepoint<-function(tablecircadianinteractions, outputadamscript){
 sapply(1:2, function(updown){
   t1<-cbind(tablecircadianinteractions[as.numeric(rownames(outputadamscript[[1]][[updown]])),], "Timepoint"=c(0, 6)[updown])
   t2<-cbind(tablecircadianinteractions[as.numeric(rownames(outputadamscript[[2]][[updown]])),], "Timepoint"=c(0, 12)[updown])
   t3<-cbind(tablecircadianinteractions[as.numeric(rownames(outputadamscript[[3]][[updown]])),], "Timepoint"=c(0, 18)[updown])
   t4<-cbind(tablecircadianinteractions[as.numeric(rownames(outputadamscript[[4]][[updown]])),], "Timepoint"=c(6, 12)[updown])
   t5<-cbind(tablecircadianinteractions[as.numeric(rownames(outputadamscript[[5]][[updown]])),], "Timepoint"=c(6, 18)[updown])
   t6<-cbind(tablecircadianinteractions[as.numeric(rownames(outputadamscript[[6]][[updown]])),], "Timepoint"=c(12, 18)[updown])
   alldiffinter<-rbind(t1, t2, t3, t4, t5, t6)
   return(alldiffinter)
 })
 #alldiffinter<-unique(alldiffinter[order(alldiffinter$V1, alldiffinter$V2, alldiffinter$V3, alldiffinter$V4, alldiffinter$V5, alldiffinter$V6),])
}
 
#>150kb
allreadsabove_difinter_timepoint<-merge_uniquedifinter_updowntimepoint(circadian_tableinteractions, diffs_readcount_above150kb)
allreadsabove_difinter_timepoint<-rbind(as.data.frame(allreadsabove_difinter_timepoint[,1]), as.data.frame(allreadsabove_difinter_timepoint[,2]))#transform into one large data frame
allreadsabove_difinter_timepoint<-unique(allreadsabove_difinter_timepoint[order(allreadsabove_difinter_timepoint$V1, allreadsabove_difinter_timepoint$V2, allreadsabove_difinter_timepoint$V3, allreadsabove_difinter_timepoint$V4, allreadsabove_difinter_timepoint$V5, allreadsabove_difinter_timepoint$V6),]) #unique interactions
#<150kb
allreadsbelow_difinter_timepoint<-merge_uniquedifinter_updowntimepoint(circadian_tableinteractions, diffs_readcount_below150kb)
allreadsbelow_difinter_timepoint<-rbind(as.data.frame(allreadsbelow_difinter_timepoint[,1]), as.data.frame(allreadsbelow_difinter_timepoint[,2]))#transform into one large data frame
allreadsbelow_difinter_timepoint<-unique(allreadsbelow_difinter_timepoint[order(allreadsbelow_difinter_timepoint$V1, allreadsbelow_difinter_timepoint$V2, allreadsbelow_difinter_timepoint$V3, allreadsbelow_difinter_timepoint$V4, allreadsbelow_difinter_timepoint$V5, allreadsbelow_difinter_timepoint$V6),]) #unique interactions

#Export data
write.table(allreadsabove_difinter_timepoint,"inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_alldiffinter_above150kb.bedpe",quote=F,sep="\t",row.names=F,col.names=F) 
write.table(allreadsbelow_difinter_timepoint,"inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_alldiffinter_below150kb.bedpe",quote=F,sep="\t",row.names=F,col.names=F) 

phasecor_diffinteralltimepoints_pe_vectorphases_timepoint<-function(intertableall){
  grangesconservedinalltimepoints_B<-GRanges(seqnames = intertableall[,1], ranges = IRanges(intertableall[,2], intertableall[,3]), "Timepoint"=intertableall[,11])
  
  grangesconservedinalltimepoints_OE<-GRanges(seqnames = intertableall[,4], ranges = IRanges(intertableall[,5], intertableall[,6]),  "Timepoint"=intertableall[,11])
      
    c3<-sapply(1:4, function(x){
      index<-findOverlaps(grangesconservedinalltimepoints_OE, list_ernasperphase_bed_mod[[x]])
      #mcols(list_ernasperphase_bed_mod[[x]][subjectHits(index)])
      index2<-findOverlaps(grangesconservedinalltimepoints_B[queryHits(index)], list_circpromsperphase_bed_mod[[x]])
      t1<-mcols(list_ernasperphase_bed_mod[[x]][subjectHits(index)[queryHits(index2)]])
      t2<-mcols(list_circpromsperphase_bed_mod[[x]][subjectHits(index2)])
      t3<-mcols(grangesconservedinalltimepoints_B[queryHits(index)][queryHits(index2)] )
      return(cbind(as.data.frame(t1)[,1], as.data.frame(t2)[,1], as.data.frame(t3)[,1]))
      })
    return(c3)
}

vv1<-phasecor_diffinteralltimepoints_pe_vectorphases_timepoint(allreadsabove_difinter_timepoint)
vv1<-rbind(vv1[[1]], vv1[[2]], vv1[[3]], vv1[[4]])

vv2<-phasecor_diffinteralltimepoints_pe_vectorphases_timepoint(allreadsbelow_difinter_timepoint)
vv2<-rbind(vv2[[1]], vv2[[2]], vv2[[3]], vv2[[4]])

#Heatmaps
my_palette <- colorRampPalette(c("black","blue"))(n = 51)

#there are only 61 circadian promoters shared in all the timepoints
heatmap.2(as.matrix(allreadsabove_difinter_timepoint[,7:10]), trace = "none", density.info = "none", dendrogram="row",labRow =  "", key.title = "", key.xlab = "Reads", labCol  = c("ZT0", "ZT6", "ZT12", "ZT18"), xlab="" ,col = colorRampPalette(c("white", "dodgerblue2"), 50), breaks = seq(0,50, 1), Colv="NA", main= paste(">150kb n=", as.character(dim(allreadsabove_difinter_timepoint)[1]), sep = ""))

heatmap.2(as.matrix(allreadsbelow_difinter_timepoint[,7:10]), trace = "none", density.info = "none", dendrogram="row",labRow =  "", labCol  = c("ZT0", "ZT6", "ZT12", "ZT18"), xlab="" , Colv="NA", main=  paste("<150kb n=", as.character(dim(allreadsbelow_difinter_timepoint)[1]), sep = ""),  key.title = "", key.xlab = "Reads",col = colorRampPalette(c("white", "dodgerblue2"), 50), breaks = seq(0,50, 1))






#Plot 
#Plots
#cor(v1[,1], v1[,2], method = "pearson")
#cor(v2[,1], v2[,2], method = "pearson")
colnames(vv1)<-c("Circadian promoters", "Oscillatory eRNAs", "HiC Timepoint")
colnames(vv2)<-c("Circadian promoters", "Oscillatory eRNAs", "HiC Timepoint")
p1<-ggplot(melt(vv1), aes(x=Var1, y=value, group=Var2, col=Var2))+geom_line(aes(linetype=Var2))+xlab("Pair prom-enh")+ylab("Phase")+theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12), axis.text.y = element_text(size = 12),axis.title.x = element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_blank(), legend.text =  element_text(size=12),  legend.position="none")+ scale_colour_brewer(palette="Dark2")+  annotate("text", x =3, y = 22, label ="")+ggtitle("> 150 Kb Alldifinter")+scale_linetype_manual(values = c(1,1,3))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2<-ggplot(melt(vv2), aes(x=Var1, y=value, group=Var2, col=Var2))+geom_line(aes(linetype=Var2))+xlab("Pair prom-enh")+ylab("Phase")+theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12), axis.text.y = element_text(size = 12),axis.title.x = element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_blank(), legend.text =  element_text(size=12), legend.justification=c(1,0), legend.position=c(1,0))+ scale_colour_brewer(palette="Dark2")+  annotate("text", x = 15, y = 22, label ="")+ggtitle("< 150 Kb Alldifinter")+scale_linetype_manual(values = c(1,1,3))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))

tiff("Figures/Paper/Differential interactions/Correction_of_circproms/Corphases_allinteractions_HiCtimepoint.tiff",height = 10, width = 25, units = 'cm', compression = "lzw", res = 300)
grid.arrange(p1,p2, ncol=2)
dev.off()




###############3
# 9.5 Repeat the script using the whole list of genes, not only circadian filtered table: Negative control 

#Divide interactions into two categories depending on the distance between the fragments
below150kb<-abs(circadian_tableinteractions[,5]-circadian_tableinteractions[,2])<=150000
above150kb<-abs(circadian_tableinteractions[,5]-circadian_tableinteractions[,2])>150000

#Apply function
ALLgenesdiffs_readcount_below150kb<-difinteractions_updown(circadian_tableinteractions[below150kb,])
###If I change dispersion 0 to 0.01 works ? why?
ALLgenesdiffs_readcount_above150kb<-difinteractions_updown(circadian_tableinteractions[above150kb,])
#diff_score_updown<-difinteractions_updown(circadian_tableinteractions_filtereds, 5)


#Merged set of dif inter; to compare with the difinter using only circproms
rownnumbermerge_uniquedifinter<-function(outputadamscript){
 sapply(1:2, function(updown){
   t1<-as.numeric(rownames(outputadamscript[[1]][[updown]]))
   t2<-as.numeric(rownames(outputadamscript[[2]][[updown]]))
   t3<-as.numeric(rownames(outputadamscript[[3]][[updown]]))
   t4<-as.numeric(rownames(outputadamscript[[4]][[updown]]))
   t5<-as.numeric(rownames(outputadamscript[[5]][[updown]]))
   t6<-as.numeric(rownames(outputadamscript[[6]][[updown]]))
   alldiffinter<-unique(c(t1, t2, t3, t4, t5, t6))
   return(alldiffinter)
 })
 #alldiffinter<-unique(alldiffinter[order(alldiffinter$V1, alldiffinter$V2, alldiffinter$V3, alldiffinter$V4, alldiffinter$V5, alldiffinter$V6),])
}

#>150kb
ALLgenesrownumberabove_difinter_timepoint<-rownnumbermerge_uniquedifinter(ALLgenesdiffs_readcount_above150kb)
circpromsrownumberabove_difinter_timepoint<-rownnumbermerge_uniquedifinter(diffs_readcount_above150kb)

#<150kb
ALLgenesrownumberbelow_difinter_timepoint<-rownnumbermerge_uniquedifinter(ALLgenesdiffs_readcount_below150kb)
circpromsrownumberbelow_difinter_timepoint<-rownnumbermerge_uniquedifinter(diffs_readcount_below150kb)

###Compare if the circprom-difinter are contained in the allgenes-difinters
#>150kb
table(unique(c(ALLgenesrownumberabove_difinter_timepoint[[1]], ALLgenesrownumberabove_difinter_timepoint[[2]])) %in%unique(c(circpromsrownumberabove_difinter_timepoint[[1]], circpromsrownumberabove_difinter_timepoint[[2]])))

#<150kb
table(unique(c(ALLgenesrownumberbelow_difinter_timepoint[[1]], ALLgenesrownumberbelow_difinter_timepoint[[2]])) %in%unique(c(circpromsrownumberabove_difinter_timepoint[[1]], circpromsrownumberbelow_difinter_timepoint[[2]])))

####Merge difinters: result will be genomic coordinates
#>150kb
ALLgeneallreadsabove_difinter_timepoint<-merge_uniquedifinter(circadian_tableinteractions, ALLgenesdiffs_readcount_above150kb)
ALLgeneallreadsabove_difinter_timepoint<-rbind(as.data.frame(ALLgeneallreadsabove_difinter_timepoint[,1]), as.data.frame(ALLgeneallreadsabove_difinter_timepoint[,2]))#transform into one large data frame
ALLgeneallreadsabove_difinter_timepoint<-unique(ALLgeneallreadsabove_difinter_timepoint[order(ALLgeneallreadsabove_difinter_timepoint$V1, ALLgeneallreadsabove_difinter_timepoint$V2, ALLgeneallreadsabove_difinter_timepoint$V3, ALLgeneallreadsabove_difinter_timepoint$V4, ALLgeneallreadsabove_difinter_timepoint$V5, ALLgeneallreadsabove_difinter_timepoint$V6),]) #unique 

#<150kb
ALLgeneallreadsbelow_difinter_timepoint<-merge_uniquedifinter(circadian_tableinteractions, ALLgenesdiffs_readcount_below150kb)
ALLgeneallreadsbelow_difinter_timepoint<-rbind(as.data.frame(ALLgeneallreadsbelow_difinter_timepoint[,1]), as.data.frame(ALLgeneallreadsbelow_difinter_timepoint[,2]))#transform into one large data frame
ALLgeneallreadsbelow_difinter_timepoint<-unique(ALLgeneallreadsbelow_difinter_timepoint[order(ALLgeneallreadsbelow_difinter_timepoint$V1, ALLgeneallreadsbelow_difinter_timepoint$V2, ALLgeneallreadsbelow_difinter_timepoint$V3, ALLgeneallreadsbelow_difinter_timepoint$V4, ALLgeneallreadsbelow_difinter_timepoint$V5, ALLgeneallreadsbelow_difinter_timepoint$V6),]) #unique 


#####Count the number of interactions according to the chic
#Merge the difinter of both distance criteria 
dataframedifinter_togranges<-function(difinters){
  grangesconservedinalltimepoints_B<-GRanges(seqnames = difinters[,1], ranges = IRanges(difinters[,2], difinters[,3]), ZT0_difint= difinters[,7], ZT6_difint= difinters[,8], ZT12_difint= difinters[,9], ZT18_difint= difinters[,10])
  grangesconservedinalltimepoints_OE<-GRanges(seqnames = difinters[,4], ranges = IRanges(difinters[,5], difinters[,6]),ZT0_difint= difinters[,7], ZT6_difint= difinters[,8], ZT12_difint= difinters[,9], ZT18_difint= difinters[,10] )
  return(grangesconservedinalltimepoints_B)
}

ALLgeneallreadsabovebelow_difinter<-unique(c(dataframedifinter_togranges(ALLgeneallreadsabove_difinter_timepoint), dataframedifinter_togranges(ALLgeneallreadsbelow_difinter_timepoint)))

#Function to count the overlaps the difinters to the chic
#Overlap of circproms with chic and count numer of interactions per chicphase
countintersfromchics_mod<-function(ZT0, ZT6, ZT12, ZT18, circproms_makingdifinter){
       grangesconservedinalltimepoints_B<-GRanges(seqnames = ZT0[,1], ranges = IRanges(ZT0[,2], ZT0[,3]))
  grangesconservedinalltimepoints_OE<-GRanges(seqnames = ZT0[,4], ranges = IRanges(ZT0[,5], ZT0[,6]) )
  zt0chic<-countOverlaps(circproms_makingdifinter, grangesconservedinalltimepoints_B)
      grangesconservedinalltimepoints_B<-GRanges(seqnames = ZT6[,1], ranges = IRanges(ZT6[,2], ZT6[,3]))
  zt6chic<-countOverlaps(circproms_makingdifinter, grangesconservedinalltimepoints_B)  
        grangesconservedinalltimepoints_B<-GRanges(seqnames = ZT12[,1], ranges = IRanges(ZT12[,2], ZT12[,3]))
  zt12chic<-countOverlaps(circproms_makingdifinter, grangesconservedinalltimepoints_B) 
          grangesconservedinalltimepoints_B<-GRanges(seqnames = ZT18[,1], ranges = IRanges(ZT18[,2], ZT18[,3]))
  zt18chic<-countOverlaps(circproms_makingdifinter, grangesconservedinalltimepoints_B) 
  
    counts<-cbind(zt0chic, zt6chic, zt12chic, zt18chic)
    rownames(counts)<-paste(as.data.frame(circproms_makingdifinter)[,1], as.data.frame(circproms_makingdifinter)[,2], as.data.frame(circproms_makingdifinter)[,3], sep = "_")
    return(counts)
}



#Apply function
countinter_difinters<-countintersfromchics_mod(ZT0, ZT6, ZT12, ZT18, ALLgeneallreadsabovebelow_difinter)
#names(countinter_difinters)<-c("ZT0", "ZT6", "ZT12", "ZT18") 

#Export lists
write.table(countinter_difinters, "inputs_for_scripts/differential_interactions/MFM_RNAseq/tablecount_interactionsALLBAITSdifinters.txt", quote = F, sep = "\t")

#Count proportions
#To look at all the max and mins counts in all the transcriptional phases :P
t<-sapply(1:4, function(phase){
  max<-sum(unlist(apply(countinter_difinters, 1, function(x){if(sum(x)>0){return(which.max(x)==phase)}})))
  min<-sum(unlist(apply(countinter_difinters, 1, function(x){if(sum(x)>0){return(which.min(x)==phase)}})))
  return(c(max, min))
  })
rownames(t)<-c("Number of Max interactions","Number of Min interactions" )
colnames(t)<-c("circpromZT0","circpromZT6", "circpromZT12", "circpromZT18" )
t

###2018-08-12  Specific questions:

#-How many promoters make differential interactions of the approximately 22,000 total genes? From the proportions they report, there does not seem to be any preference for the time point, which is different from the circproms.
Allbaits<-read.csv("inputs_for_scripts/Baits.txt", header = F, sep = "\t")
grangesAllbaits<-GRanges(seqnames = Allbaits[,1], ranges = IRanges(Allbaits[,2], Allbaits[,3]), strand = Allbaits[,4], Genename=Allbaits[,5])
remove(Allbaits)
length(unique(subjectHits(findOverlaps(grangesAllbaits, ALLgeneallreadsabovebelow_difinter))))#13172 
dim(unique(mcols(grangesAllbaits[queryHits(findOverlaps(grangesAllbaits, ALLgeneallreadsabovebelow_difinter))]))) #13171/22226

#-Of the 1,256 circadian proms, how many make differential interactions? 949/1195
length(unique(sort(subjectHits(findOverlaps(ALLgeneallreadsabovebelow_difinter, c(list_circpromsperphase_bed[[1]], list_circpromsperphase_bed[[2]], list_circpromsperphase_bed[[3]], list_circpromsperphase_bed[[4]])))))) #949/1195


#what is the proportion of the set of genes that make differentials that fall into the category of transcriptionally active.
activegenesliver<-read.csv("RNAseq_Liver/active_genes_liver.txt", header = F, sep = "\t")

activegenesliverGR<-GRanges(seqnames = paste("chr", activegenesliver[,2], sep=""), ranges = IRanges(activegenesliver[,3], activegenesliver[,4]), strand = activegenesliver[,5], Genename=activegenesliver[,1])
remove(activegenesliver)
length(unique(subjectHits(findOverlaps(query = activegenesliverGR, subject = ALLgeneallreadsabovebelow_difinter))))#7468/ 13235

dim(unique(mcols(activegenesliverGR[queryHits(findOverlaps(query = activegenesliverGR, subject =  ALLgeneallreadsabovebelow_difinter))]))) #7025/9418


#Diffinters export files:Updated 2020-07-06

#Merge above and below 150kb dif inters
 library(dplyr)
t<-unique(rbind(allreadsabove_difinter_timepoint,allreadsbelow_difinter_timepoint))
#Ctl
#nrow(match_df(circadian_tableinteractions, t))#6047
#Only interactions mediated by circproms# 242920 = 6047
nondifinterall<-(anti_join(circadian_tableinteractions, t))

write.table(nondifinterall, "inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_allNONdiffinter_abovebelow150kb2020includingnotMFMcps.bedpe", quote = F, row.names = F, col.names = F, sep = "\t")
