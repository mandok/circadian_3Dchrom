#####################---Enrichments:eRNAs static vs osceRNAs---##############

#Create a random set of no osc ernas, keeping the size of osc ernas and see overlap with circprom
chic<-read.csv("/inputs_for_scripts/All_step2_washU_text_symm.txt", sep = "\t", header = F)

#Symmetric set is the complete one
#chic_baits<-read.csv("chic_baits_mm9.txt", sep = "\t", header = F)
chic_baits<-chic[,1:3]
#chic_otherends<-read.csv("chic_otherends_mm9.txt", sep = "\t", header = F)
chic_otherends<-chic[,4:6]
#Same number of lines! Allows mapping bait<->otherend

#Bedfiles, GRanges Objects
chic_otherends_bed<-GRanges(seqnames= Rle(chic_otherends[,1]), ranges = IRanges(chic_otherends[,2], chic_otherends[,3]))
chic_bait_bed<-GRanges(seqnames= Rle(chic_baits[,1]), ranges = IRanges(chic_baits[,2], chic_baits[,3]))

#Bed eRNAs 
osc_ernas<-read.csv("eRNAs_de_novo_oscillating_phases.txt", header = T, sep = "\t")
osc_ernas<-osc_ernas[,1:3]
osc_ernas<-bedfile(osc_ernas, columnnames)

###Control randompicked, use no osc eRNAs
poolofenas<-bedfile(read.csv("eRNAs_de_novo_oscillating.txt", header = T, sep = "\t"), columnnames)
noosc_ernas<-poolofenas[!(1:19086 %in% queryHits(findOverlaps(poolofenas, (osc_ernas), type="equal")))]


enrichoverlap_noosc_ernas_to_chicEO<-function(chic_bed_otherends, noosc_ernas, osc_ernas, rep)
{
  Et2<-vector()
  sapply(1:rep, function(x){
    #First part is to build a random noosc dataset with the same number of ernas as osc_ernas 
    random_noosernas<-sample(noosc_ernas, length(osc_ernas), replace = F)
    
    #Second part is to make the overlap using the genomic ranges function subsetByOverlaps
    #USE #Subsetbyoverlaps
    El1_1<-retrieveindex_fromchic(random_noosernas, chic_otherends_bed, chic_bait_bed)
    Et2<-c(Et2, length(El1_1))
  }
  
  )
}


EO_oscernas_tochic_randomset<-enrichoverlap_noosc_ernas_to_chicEO(chic_bed_otherends, noosc_ernas, osc_ernas, 500)
EO_oscernas_tochic_randomset_mean<-mean(EO_oscernas_tochic_randomset)
EO_oscernas_tochic_randomset_sd<-sd(EO_oscernas_tochic_randomset)
EO_oscernas_tochic_randomset_HI<-EO_oscernas_tochic_randomset_mean+EO_oscernas_tochic_randomset_sd
EO_oscernas_tochic_randomset_LI<-EO_oscernas_tochic_randomset_mean-EO_oscernas_tochic_randomset_sd

osc_ernas_BI<-subjectLength(findOverlaps(query = chic_otherends_bed, subject = osc_ernas))

EO_oscernas_tochic_randomset_plot<-rbind("oscernas_observed"=(osc_ernas_BI),"oscenas_expected"=EO_oscernas_tochic_randomset_mean)
barplot(EO_oscernas_tochic_randomset_plot,cex.names = 1.2, beside = TRUE, las=1, main = "Overlap with Capture-hiC \n fragments (only other ends)", xlab = "", names.arg = "Oscillatory eRNAs", ylab="Number of overlaps", width = 0.1, border = "white", space = 0.5)
legend("topright", c("Observed","Expresed"),pch=16 ,col=c("black","gray"), bty = "n",cex = 1.2, y.intersp=0.2)
arrows(1:2, c(0,EO_oscernas_tochic_randomset_LI), 1:2, c(0,EO_oscernas_tochic_randomset_HI), length=0.5, angle=90, code=3,col = "red")

#No need to make bins, ernas are always 500bp; only uses fang data set; enrichment of association osc_erna with a circadian promoter
enrichoverlap_noosc_ernas<-function(chic_bed_otherends, noosc_ernas, osc_ernas, rep, mfm_mrnas)
{
  Et2<-vector()
  sapply(1:rep, function(x){
    #First part is to build a random noosc dataset with the same number of ernas as osc_ernas 
    random_noosernas<-sample(noosc_ernas, length(osc_ernas), replace = F)
    
    #Second part is to make the overlap using the genomic ranges function subsetByOverlaps
    #USE #Subsetbyoverlaps
    indexes<-retrieveindex_fromchic(random_noosernas, chic_otherends_bed, chic_bait_bed)
    El1_1<-count_matches_circprom_ernasindexes(indexes,mfm_mrnas)
    Et2<-c(Et2, El1_1)
  }
  
  )
}
#Apply function to calculate the expectd number of interactions between eRNAs and cirpcproms MFM; 100 iterations
EO_nooscernas_randomset<-enrichoverlap_noosc_ernas(chic_bed_otherends, noosc_ernas, osc_ernas, 100, mfm_mrnasINTRONS)
#Calculate the mean, sd, higher interval and lower interval for all the iterations
EO_nooscernas_randomset_mean<-mean(EO_nooscernas_randomset)
EO_nooscernas_randomset_sd<-sd(EO_nooscernas_randomset)
EO_nooscernas_randomset_HI<-EO_nooscernas_randomset_mean+EO_nooscernas_randomset_sd
EO_nooscernas_randomset_LI<-EO_nooscernas_randomset_mean-EO_nooscernas_randomset_sd

##################Otherends from differential interactions
#Run differentialinteractions_adamscript, use all the unique pairs
load("/inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_diffs_readcount_above150kb2018")
load("/inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_diffs_readcount_below150kb2018")
load("/inputs_for_scripts/circtablesinter")

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


difOE_above_symm<-unique(sort(GRanges(seqnames = paste("chr", rbind(allreadsabove_difinter_timepoint[,1], allreadsabove_difinter_timepoint[,4]), sep = ""), ranges = IRanges(rbind(allreadsabove_difinter_timepoint[,2], allreadsabove_difinter_timepoint[,5]),rbind(allreadsabove_difinter_timepoint[,3], allreadsabove_difinter_timepoint[,6])))))
difOE_below_symm<-unique(sort(GRanges(seqnames = paste("chr", rbind(allreadsbelow_difinter_timepoint[,1], allreadsbelow_difinter_timepoint[,4]), sep = ""), ranges = IRanges(rbind(allreadsbelow_difinter_timepoint[,2], allreadsbelow_difinter_timepoint[,5]),rbind(allreadsbelow_difinter_timepoint[,3], allreadsbelow_difinter_timepoint[,6])))))

difOE_above<-unique(sort(GRanges(seqnames =   allreadsabove_difinter_timepoint[,4], 
                                 ranges = IRanges( allreadsabove_difinter_timepoint[,5],allreadsabove_difinter_timepoint[,6]))))
difOE_below<-unique(sort(GRanges(seqnames =  allreadsbelow_difinter_timepoint[,4], 
                                 ranges = IRanges( allreadsbelow_difinter_timepoint[,5], allreadsbelow_difinter_timepoint[,6]))))


#List of OE that are nt  differential interactions
t<-!(1:dim(circadian_tableinteractions)[1] %in% rownames(allreadsabove_difinter_timepoint))
nondifOE_above_symm<-unique(sort(GRanges(seqnames = paste("chr", rbind(circadian_tableinteractions[t,1], circadian_tableinteractions[t,4]), sep = ""), ranges = IRanges(rbind(circadian_tableinteractions[t,2], circadian_tableinteractions[t,5]),rbind(circadian_tableinteractions[t,3], circadian_tableinteractions[t,6])))))
nondifOE_above<-unique(sort(GRanges(seqnames =  circadian_tableinteractions[t,4], 
                                         ranges = IRanges( circadian_tableinteractions[t,5],circadian_tableinteractions[t,6]))))

t<-!(1:dim(circadian_tableinteractions)[1] %in% rownames(allreadsbelow_difinter_timepoint))
nondifOE_below_symm<-unique(sort(GRanges(seqnames = paste("chr", rbind(circadian_tableinteractions[t,1], circadian_tableinteractions[t,4]), sep = ""), ranges = IRanges(rbind(circadian_tableinteractions[t,2], circadian_tableinteractions[t,5]),rbind(circadian_tableinteractions[t,3], circadian_tableinteractions[t,6])))))
nondifOE_below<-unique(sort(GRanges(seqnames =  circadian_tableinteractions[t,4], 
                                    ranges = IRanges( circadian_tableinteractions[t,5],circadian_tableinteractions[t,6]))))


#Obs
obs_difoe_ab<-length(queryHits(findOverlaps(difOE_above, osc_ernas)))
obs_difoe_ab_noosc<-length(queryHits(findOverlaps(difOE_above, noosc_ernas)))
obs_difoe_be<-length(queryHits(findOverlaps(difOE_below, osc_ernas)))
obs_difoe_be_noosc<-length(queryHits(findOverlaps(difOE_below, noosc_ernas)))

###########Exp randompick nondif
enrichoverlap_oscernas_nondifOE<-function(nondifOE_above_symm, difOE_above_symm, osc_ernas, rep)
{
  Et2<-vector()
  sapply(1:rep, function(x){
    #First part is to build a random noosc dataset with the same number of ernas as osc_ernas 
    random_noosernas<-sample(nondifOE_above_symm, length(difOE_above_symm), replace = F)
    #Second part is to make the overlap using the genomic ranges function subsetByOverlaps
    #USE #Subsetbyoverlaps
    El1_1<-length(unique(subjectHits(findOverlaps(osc_ernas, random_noosernas))))
    Et2<-c(Et2, El1_1)
  }
  
  )
}
exp_difoe_ab<-enrichoverlap_oscernas_nondifOE(nondifOE_above, difOE_above, osc_ernas, 100)
exp_difoe_ab_noosc<-enrichoverlap_oscernas_nondifOE(nondifOE_above, difOE_above, noosc_ernas, 100)

exp_difoe_be<-enrichoverlap_oscernas_nondifOE(nondifOE_below, difOE_below, osc_ernas, 100)
exp_difoe_be_noosc<-enrichoverlap_oscernas_nondifOE(nondifOE_below, difOE_below, noosc_ernas, 100)

############Stat test
t.test(exp_difoe_ab, mu = obs_difoe_ab)
t.test(exp_difoe_be, mu = obs_difoe_be)
t.test(exp_difoe_ab_noosc, mu = obs_difoe_ab_noosc)
t.test(exp_difoe_be_noosc, mu = obs_difoe_be_noosc)
#########Plot
t<-(rbind(cbind("> 150 Kb"= obs_difoe_ab/ mean(exp_difoe_ab), "> 150 Kb "=obs_difoe_ab_noosc /mean(exp_difoe_ab_noosc)),
         cbind(  "< 150 Kb "=obs_difoe_be/mean(exp_difoe_be), "< 150 Kb "=obs_difoe_be_noosc/ mean(exp_difoe_be_noosc))))
t1<-c("> 150 Kb"= obs_difoe_ab, "< 150 Kb "=obs_difoe_be,"> 150 Kb "=obs_difoe_ab_noosc ,
       "< 150 Kb "=obs_difoe_be_noosc)
t2<-c("> 150 Kb"= mean(exp_difoe_ab),"< 150 Kb "=mean(exp_difoe_be), "> 150 Kb "=mean(exp_difoe_ab_noosc),
       "< 150 Kb "=mean(exp_difoe_be_noosc))
t3<-c("> 150 Kb"= sd(exp_difoe_ab),  "< 150 Kb "=sd(exp_difoe_be),"> 150 Kb "=sd(exp_difoe_ab_noosc),
      "< 150 Kb "=sd(exp_difoe_be_noosc))
#svglite::svglite("enrichmenternas_oe.svg")

b<-barplot(t,beside = T,  ylim=c(0,4),names.arg=c("Osc eRNAs", "No Osc eRNAs"),col= c("gray30", " lightsteelblue3"), ylab = "Differential PIRs/ non differential PIRs", main="Enrichment of eRNAs in PIRs",border = "white")
arrows( b, (t1/(t2+t3)) ,b, (t1/(t2-t3)), length=0.2, angle=90, code=3)
#error.bar(b, (t1/(t2+t3)),(t1/(t2-t3)))
legend("topright",legend =  c("> 150 Kb", "< 150 Kb"),fill= c("gray30", " lightsteelblue3"), bty = "n", border="white" )
abline(h=1, lty=3, col = "black")

#dev.off()


