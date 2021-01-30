#######################################
#2019-01-15 using MFM RNAseq proms; and differentila interactions
#Random pick, p-e; NO INTRONS


#------------
#RANDOM ERNAS, now try random proms
#Extra: EXPECTED
#Repeat 100 times
#Take a random list everytime 
#----------------------------
#Obtain new chic files in granges format from the differential interactions

library(rtracklayer)

load("/inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_diffs_readcount_above150kb2018")
load("/inputs_for_scripts/differential_interactions/MFM_RNAseq/MFMcircproms_diffs_readcount_below150kb2018")

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

Baits_diffinter_GR<-(c(makeGRangesFromDataFrame(allreadsbelow_difinter_timepoint[,1:3], seqnames.field = "V1", start.field = "V2", end.field = "V3"), makeGRangesFromDataFrame(allreadsabove_difinter_timepoint[,1:3], seqnames.field = "V1", start.field = "V2", end.field = "V3") ))

OE_diffinter_GR<-(c(makeGRangesFromDataFrame(allreadsbelow_difinter_timepoint[,4:6], seqnames.field = "V4", start.field = "V5", end.field = "V6"), makeGRangesFromDataFrame(allreadsabove_difinter_timepoint[,4:6], seqnames.field = "V4", start.field = "V5", end.field = "V6") ))
----------------------------------
 
#Take a random list everytime 
#Random no osc ernas
expbothphases_randompickernas<-sapply(1:100, function(rep){
  RANlist_ernasperphase_bed<-sapply(1:8, function(phase){
    pos<-sample(list_ernasperphase_bed[[phase]] ,160, replace = F)
    })
                #lapply(list_ernasperphase_bed, length)[[phase]], replace = F)
#Random no osc proms
  #RANlist_circpromperphase_bed<-sapply(1:8, function(phase){
  #pos<-sample(list_circpromsperphase_bed[[phase]],109, replace = F)
  #})
#Retrieve indexes of eRNAs
  RANindexes_baitsofernas_perphases<-lapply(RANlist_ernasperphase_bed, function(x){
  tempo<-retrieveindex_fromchic(x, OE_diffinter_GR, Baits_diffinter_GR)
  })
  
  #Bait(index ernas) overlap with proms

  RANbothphases_circadianprom_ernas_counts<-
  lapply(list_circpromsperphase_bed, function(tempo_circ){# each circprom
    t<-lapply(RANindexes_baitsofernas_perphases, function(tempo_ernas){
    tempo_ernas<-length(subjectHits(findOverlaps(tempo_ernas, tempo_circ)))
    })
})
  names(RANbothphases_circadianprom_ernas_counts)<-circpromlist_phases
  
  
  #names(RANoverlaps_and_correlations_bothphases_circadianprom_ernas)<-circpromlist_phases
  return(RANbothphases_circadianprom_ernas_counts)

})
expbothphases_randompickproms<-sapply(1:100, function(rep){
  #RANlist_ernasperphase_bed<-sapply(1:8, function(phase){
  #  pos<-sample(list_ernasperphase_bed[[phase]] ,160, replace = F)
  #  })
                #lapply(list_ernasperphase_bed, length)[[phase]], replace = F)
#Random no osc proms
  RANlist_circpromperphase_bed<-sapply(1:4, function(phase){ # 4 groups for mfm;8 for fang
  pos<-sample(list_circpromsperphase_bed[[phase]],193, replace = F) #103 for fang ; 193 for MFM
  })
#Retrieve indexes of eRNAs
  indexes_baitsofernas_perphases<-lapply(list_ernasperphase_bed, function(x){
  tempo<-retrieveindex_fromchic(x, OE_diffinter_GR, Baits_diffinter_GR)
  })
  
  #Bait(index ernas) overlap with proms

  RANbothphases_circadianprom_ernas_counts<-
  lapply(RANlist_circpromperphase_bed, function(tempo_circ){# each circprom
    t<-lapply(indexes_baitsofernas_perphases  , function(tempo_ernas){
    tempo_ernas<-length(subjectHits(findOverlaps(tempo_ernas, tempo_circ)))
    })
})
  names(RANbothphases_circadianprom_ernas_counts)<-circpromlist_phases
  
  
  #names(RANoverlaps_and_correlations_bothphases_circadianprom_ernas)<-circpromlist_phases
  return(RANbothphases_circadianprom_ernas_counts)

})

#Obtain the mean, 
mean_exp_bothphases_randompickproms<-sapply(1:4, function(cirpromphase){sapply(1:8, function(ernasphase){mean(unlist(sapply(1:100, function(rep){expbothphases_randompickproms[cirpromphase,][[rep]][ernasphase]})))})})

mean_exp_bothphases_randompickernas<-sapply(1:8, function(ernasphase){sapply(1:4, function( cirpromphase){mean(unlist(sapply(1:100, function(rep){expbothphases_randompickernas[cirpromphase,][[rep]][ernasphase]})))})})

 


mean_exp_bothphases_randompickproms
mean_exp_bothphases_randompickernas

#------------
#Plot
###DAY
tempoclock<-sapply(1:4, function(x){
  unlist(bothphases_circadianprom_ernas_counts[[x]])
})
tempoclock_dayRP<-mean_exp_bothphases_randompickernas[,1:4]
rownames(tempoclock_dayRP)<-c("0", "6", "12", "18")
colnames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12")
tempoclock_nightRP<-mean_exp_bothphases_randompickernas[,5:8]
rownames(tempoclock_nightRP)<-c("0", "6", "12", "18")
colnames(tempoclock_nightRP)<-c("13-15", "16-18", "19-21", "22-24")
tempoclock_dayRP<-mean_exp_bothphases_randompickproms[,1:2]
rownames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_dayRP)<-c("0", "6")
tempoclock_nightRP<-mean_exp_bothphases_randompickproms[,3:4]
rownames(tempoclock_nightRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_nightRP)<-c("12", "18")


###Control randompicked, use no osc eRNAs
poolofenas<-bedfile(read.csv("/eRNAs_de_novo_oscillating.txt", header = T, sep = "\t"), columnnames)
noosc_ernas<-poolofenas[!(1:19086 %in% queryHits(findOverlaps(poolofenas, bedfile(osc_ernas), type="equal")))]

#Retrieve non circadian promoters, 
#First retrieve only baits
onlybaits<-bedfile(read.csv("/Baits_Gen_Res/Baits.txt", header = T, sep="\t"), columnnames)
t<-chic_bait_bed[unique(queryHits(findOverlaps(chic_bait_bed, onlybaits)))]
noncircproms<-t[!(1:247943 %in% queryHits(findOverlaps(t,fang_mrnas)))]
#247943 289336

expbothphases_randompicknooscproms_control<-sapply(1:100, function(rep){
  #RANlist_ernasperphase_bed<-sapply(1:8, function(phase){
   # pos<-sample(noosc_ernas ,160, replace = F)
    #})
                #lapply(list_ernasperphase_bed, length)[[phase]], replace = F)
#Random no osc proms
  RANlist_circpromperphase_bed<-sapply(1:4, function(phase){
  pos<-sample(unique(noncircproms), 193, replace=F)})#103 for fang ; 193 for MFM
    #sample(list_circpromsperphase_bed[[phase]],109, replace = F)})
#Retrieve indexes of eRNAs
  RANindexes_baitsofernas_perphases<-lapply(list_ernasperphase_bed, function(x){
  tempo<-retrieveindex_fromchic(x,OE_diffinter_GR, Baits_diffinter_GR)
  })
  
  #Bait(index ernas) overlap with proms

  RANbothphases_circadianprom_ernas_counts<-
  lapply(RANlist_circpromperphase_bed, function(tempo_circ){# each circprom
    t<-lapply(RANindexes_baitsofernas_perphases, function(tempo_ernas){
    tempo_ernas<-length(subjectHits(findOverlaps(tempo_ernas, tempo_circ)))
    })
})
  names(RANbothphases_circadianprom_ernas_counts)<-circpromlist_phases
  
  
  #names(RANoverlaps_and_correlations_bothphases_circadianprom_ernas)<-circpromlist_phases
  return(RANbothphases_circadianprom_ernas_counts)

})


expbothphases_randompicknooscernas_control<-sapply(1:100, function(rep){
  RANlist_ernasperphase_bed<-sapply(1:8, function(phase){
    pos<-sample(noosc_ernas ,160, replace = F)
    })
                #lapply(list_ernasperphase_bed, length)[[phase]], replace = F)
#Random no osc proms
  #RANlist_circpromperphase_bed<-sapply(1:8, function(phase){
  #pos<-sample(unique(noncircproms), 109, replace=F)})
    #sample(list_circpromsperphase_bed[[phase]],109, replace = F)})
#Retrieve indexes of eRNAs
  RANindexes_baitsofernas_perphases<-lapply(RANlist_ernasperphase_bed, function(x){
  tempo<-retrieveindex_fromchic(x,OE_diffinter_GR, Baits_diffinter_GR)
  })
  
  #Bait(index ernas) overlap with proms

  RANbothphases_circadianprom_ernas_counts<-
  lapply(list_circpromsperphase_bed, function(tempo_circ){# each circprom
    t<-lapply(RANindexes_baitsofernas_perphases, function(tempo_ernas){
    tempo_ernas<-length(subjectHits(findOverlaps(tempo_ernas, tempo_circ)))
    })
})
  names(RANbothphases_circadianprom_ernas_counts)<-circpromlist_phases
  
  
  #names(RANoverlaps_and_correlations_bothphases_circadianprom_ernas)<-circpromlist_phases
  return(RANbothphases_circadianprom_ernas_counts)

})



#Obtain the mean, 
mean_exp_bothphases_randompicknooscproms_control<-sapply(1:4, function(cirpromphase){sapply(1:8, function(ernasphase){mean(unlist(sapply(1:100, function(rep){expbothphases_randompicknooscproms_control[cirpromphase,][[rep]][ernasphase]})))})})

mean_exp_bothphases_randompicknooscernas_control<-sapply(1:8, function(ernasphase){sapply(1:4, function(cirpromphase){mean(unlist(sapply(1:100, function(rep){expbothphases_randompicknooscernas_control[cirpromphase,][[rep]][ernasphase]})))})})

 
mean_exp_bothphases_randompicknooscproms_control
mean_exp_bothphases_randompicknooscernas_control

#Plot day night and control, general not per hour
#mean_exp_bothphases_randompickproms
#mean_exp_bothphases_randompickernas
tempoclock_dayRP<-mean_exp_bothphases_randompickernas[,1:4]
#rownames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
#colnames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12")
#tempoclock_nightRP<-mean_exp_bothphases_randompickernas[,5:8]
#rownames(tempoclock_nightRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
#colnames(tempoclock_nightRP)<-c("13-15", "16-18", "19-21", "22-24")


#Proms anchors

tempoclock_dayRP<-mean_exp_bothphases_randompickproms[,1:2] #1:2 MFM and 1:4 fang #mean_exp_bothphases_randompickproms[,1:2] for MFM and mean_exp_bothphases_randompickproms[,1:4] for fang
rownames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_dayRP)<-c("0", "6") #c("0", "6") for MFM and c("1-3", "4-6", "7-9", "10-12") for fang
tempoclock_nightRP<-mean_exp_bothphases_randompickproms[,3:4] #3:4 MFM and 5:8 for fang
rownames(tempoclock_nightRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_nightRP)<-c("12", "18") # c("12", "18") for MFM and c("13-15", "16-18", "19-21", "22-24")

tempoclock_nightRP_control<-mean_exp_bothphases_randompicknooscproms_control[,3:4]
rownames(tempoclock_nightRP_control)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_nightRP_control)<-c("12", "18") #c("12", "18") for MFM and c("13-15", "16-18", "19-21", "22-24") for fang

tempoclock_dayRP_control<-mean_exp_bothphases_randompicknooscproms_control[,1:2]
rownames(tempoclock_dayRP_control)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_dayRP_control)<-c("0", "6")#c("0", "6") for MFM and c("1-3", "4-6", "7-9", "10-12") for fang


t<-cbind(apply(tempoclock_dayRP, 1, mean)/sum(apply(tempoclock_dayRP, 1, mean)), apply(tempoclock_nightRP, 1, mean)/sum(apply(tempoclock_nightRP, 1, mean)), apply(mean_exp_bothphases_randompicknooscproms_control, 1, mean)/sum(apply(mean_exp_bothphases_randompicknooscproms_control, 1, mean)))

colnames(t)<-c("Diurnal", "Nocturnal","Control")


p1<-ggplot(melt(t), aes(x=Var2, y=value, fill=Var1))+geom_bar(stat = "identity")+
  theme( panel.background = element_blank(), axis.text=element_text(size=12, vjust=0.8), axis.line.x = element_line(color="white", size = 1),axis.line.y = element_line(color="white", size = 1), axis.text.y = element_blank(),axis.title=element_text(size=12), axis.ticks.y = element_blank())+
    scale_fill_brewer(palette="RdBu", name="Phases")+
  #scale_x_discrete(labels= c("1", "4", "7", "10", "13", "16", "19", "22"))+
  ylab("")+xlab("")+ggtitle("Circadian Promoters")
p1<-p1+theme(axis.text.x = element_text(angle = 90, hjust = 1))


#ernas anchors
tempoclock_dayRP<-mean_exp_bothphases_randompickernas[,1:4]
rownames(tempoclock_dayRP)<-c("0", "6", "12", "18")#c("0", "6", "12", "18") for MFM and c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24") for fang
colnames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12")
tempoclock_nightRP<-mean_exp_bothphases_randompickernas[,5:8]
rownames(tempoclock_nightRP)<-c("0", "6", "12", "18")#c("0", "6", "12", "18") for MFM and c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24") for fang
colnames(tempoclock_nightRP)<-c("13-15", "16-18", "19-21", "22-24")

tempoclock_nightRP_control<-mean_exp_bothphases_randompicknooscernas_control[,5:8]
rownames(tempoclock_nightRP)<-c("ZT0", "ZT6", "ZT12", "ZT18")#c("ZT0", "ZT6", "ZT12", "ZT18") for MFM and c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24") for fang
colnames(tempoclock_nightRP)<-c("13-15", "16-18", "19-21", "22-24")

tempoclock_dayRP_control<-mean_exp_bothphases_randompicknooscernas_control[,1:4]
rownames(tempoclock_dayRP)<-c("ZT0", "ZT6", "ZT12", "ZT18")#c("ZT0", "ZT6", "ZT12", "ZT18") for MFM and c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24") for fang
colnames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12")


t<-cbind(apply(tempoclock_dayRP, 1, mean)/sum(apply(tempoclock_dayRP, 1, mean)), apply(tempoclock_nightRP, 1, mean)/sum(apply(tempoclock_nightRP, 1, mean)), apply(mean_exp_bothphases_randompicknooscernas_control, 1, mean)/sum(apply(mean_exp_bothphases_randompicknooscernas_control, 1, mean)))
colnames(t)<-c("Diurnal", "Nocturnal","Control")

p2<-ggplot(melt(t), aes(x=Var2, y=value, fill=Var1))+geom_bar(stat = "identity")+
  theme( panel.background = element_blank(), axis.text=element_text(size=12, vjust=0.8), axis.line.x = element_line(color="white", size = 1),axis.line.y = element_line(color="white", size = 1), axis.text.y = element_blank(),axis.title=element_text(size=12), axis.ticks.y = element_blank())+
    scale_fill_brewer(palette="RdBu", name="Phases")+
  #scale_x_discrete(labels= c("1", "4", "7", "10", "13", "16", "19", "22"))+
  ylab("")+xlab("")+ggtitle("eRNAs")
p2<-p2+theme(axis.text.x = element_text(angle = 90, hjust = 1))



#svglite::svglite("Analysis_per_phases/MFM_RNAseqrandompicked_daynightcontrol_anchorprom_anchorpromintron_anchorernas_proportions.svg")
grid.arrange(p1, p2, p3, ncol=3)
#dev.off()

#fang control only p1 and p2
#svglite::svglite("Differential_interactions/Analysis_per_phases/MFMRNAseq_GROseqrandompicked_daynightcontrol_anchorprom_anchorernas_proportions.svg")
grid.arrange(p1, p2,ncol=2)
#dev.off()

#Stats

#Statistical test to see if there is ss difference between D-C N-C and D-N

#mean_exp_bothphases_randompickproms
#mean_exp_bothphases_randompickernas
tempoclock_dayRP<-mean_exp_bothphases_randompickproms[,1:2]
rownames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_dayRP)<-c("0", "6")
tempoclock_nightRP<-mean_exp_bothphases_randompickproms[,3:4]
rownames(tempoclock_nightRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_nightRP)<-c("12", "18")

#For fishers test, prom anchors no proportions
t<-cbind(apply(tempoclock_dayRP, 1, mean), apply(tempoclock_nightRP, 1, mean), apply(mean_exp_bothphases_randompicknooscproms_control, 1, mean))
colnames(t)<-c("Diurnal", "Nocturnal", "Control")
#fisher.test(matrix(c(sum(t[1:4,1]), sum(t[5:8, 1]), sum(t[1:4, 3]), sum(t[5:8,3])), byrow = T, ncol=2))#p-value = 0.5471
#fisher.test(matrix(c(sum(t[1:4,2]), sum(t[5:8, 2]), sum(t[1:4, 3]), sum(t[5:8,3])), byrow = T, ncol=2))#p-value =  0.1304
#fisher.test(matrix(c(sum(t[1:4,1]), sum(t[5:8, 1]), sum(t[1:4, 2]), sum(t[5:8,2])),  byrow = T, ncol=2)) #the only one significant, p-value = 0.005327


wilcox.test(t[,1], t[,3], paired = T) #Diurnal vs control p-value = 0.007813
wilcox.test(t[,2], t[,3], paired = T) #Nocturnal vs control p-value =  0.007813
wilcox.test(t[,1], t[,2], paired = T) #Diurnal vs Nocturnal p-value = 0.1484

#svglite::svglite("Differential_interactions/Analysis_per_phases/MFM_RNAseqrandompicked_daynightcontrol_anchorprom_anchorproms_rawDIFFERENTIALinteractionscounts.svg")
ggplot(melt(t), aes(x=Var1, y=value, group=Var2, col=Var2))+geom_line()+xlab("eRNAs phase")+ylab("Number of interaction with eRNAs")+ggtitle("Circproms")+scale_color_brewer(palette="Set1")+scale_color_brewer(palette="Set1")+guides(color=guide_legend(title="Circproms phase"))+theme_classic()
#dev.off()

#For fishers test, ernas anchor no proportions
tempoclock_dayRP<-mean_exp_bothphases_randompickernas[,1:4]
rownames(tempoclock_dayRP)<-c("0", "6", "12", "18")
colnames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12")
tempoclock_nightRP<-mean_exp_bothphases_randompickernas[,5:8]
rownames(tempoclock_nightRP)<-c("0", "6", "12", "18")


t<-cbind(apply(tempoclock_dayRP, 1, mean), apply(tempoclock_nightRP, 1, mean), apply(mean_exp_bothphases_randompicknooscernas_control, 1, mean))
colnames(t)<-c("Diurnal", "Nocturnal", "Control")
#fisher.test(matrix(c(sum(t[1:4,1]), sum(t[5:8, 1]), sum(t[1:4, 3]), sum(t[5:8,3])), byrow = T, ncol=2))#p-value = 0.4959
#fisher.test(matrix(c(sum(t[1:4,2]), sum(t[5:8, 2]), sum(t[1:4, 3]), sum(t[5:8,3])), byrow = T, ncol=2))#p-value = 0.4935
#fisher.test(matrix(c(sum(t[1:4,1]), sum(t[5:8, 1]), sum(t[1:4, 2]), sum(t[5:8,2])),  byrow = T, ncol=2))#closest to significant  0.1297


wilcox.test(t[,1], t[,3], paired = T) #Diurnal vs control p-value = 0.625
wilcox.test(t[,2], t[,3], paired = T) #Nocturnal vs control p-value =  0.125
wilcox.test(t[,1], t[,2], paired = T) #Diurnal vs Nocturnal p-value = 0.875

#svglite::svglite("Differential_interactions/Analysis_per_phases/MFM_RNAseqrandompicked_daynightcontrol_anchorprom_anchorernas_rawDIFFERENTIALinteractionscounts.svg")
ggplot(melt(t), aes(x=Var1, y=value, group=Var2, col=Var2))+geom_line()+xlab("Circproms phase")+ylab("Number of interaction with circproms")+ggtitle("eRNAs")+scale_color_brewer(palette="Set1")+guides(color=guide_legend(title="eRNAs phase"))+theme_classic()
#dev.off()

    ###----all circproms pvals
#1) Mean of diurnals ZT0 Y 6 for each permutation, I will have 100 values for each of the ZT ernas for the diurnal and noc Circproms
#2) With two distributions for each timepoint i can compare. In the end I will have a vector of 8 pvals 

zt0circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickproms[1,perm])[ernasphase]})})

zt6circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickproms[2,perm])[ernasphase]})})

zt12circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickproms[3,perm])[ernasphase]})})

zt18circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickproms[4,perm])[ernasphase]})})


#mean of diurnal and noc
diurcircproms<-(zt0circ+zt6circ)/2
noccircproms<-(zt12circ+zt18circ)/2

#Check if distributions are normal or not
nd<-sapply(c(1:8), function(x){shapiro.test(diurcircproms[,x])[[2]]}) #remove col 5 because its an identical distribution
nn<-sapply(1:8, function(x){shapiro.test(noccircproms[,x])[[2]]})
par(mfrow=c(1,2))
plot(nn, ylim=c(0,1))+abline(h = .001, col="red")
plot(nd)+abline(h = .001, col="red")
#dev.off()
#normal: 1.2,3 7.8
 #not normal: 4,5,6
pvalsallcircproms<-c(wilcox.test(diurcircproms[,1],noccircproms[,1] )[[3]],
                wilcox.test(diurcircproms[,2],noccircproms[,2] )[[3]], 
                wilcox.test(diurcircproms[,3],noccircproms[,3] )[[3]], 
                wilcox.test(diurcircproms[,4],noccircproms[,4] )[[3]], 
                wilcox.test(diurcircproms[,5],noccircproms[,5] )[[3]], 
                wilcox.test(diurcircproms[,6],noccircproms[,6] )[[3]], 
                wilcox.test(diurcircproms[,7],noccircproms[,7] )[[3]], 
                wilcox.test(diurcircproms[,7],noccircproms[,7] )[[3]])
 
 names(pvalsallcircproms)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
 
#svglite::svglite("Differential_interactions/Analysis_per_phases/MFM_RNAseqrandompicked_daynightcontrol_pvals_diurvsnocALLPROMS.svg")
ggplot(melt(log10(pvalsallcircproms)), aes(x=c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24"), y=value))+geom_point(size=4,aes(col=rownames(melt(log10(pvalsallcircproms)))))+
    geom_point(shape = 1,size=4,colour = "gray49")+ scale_y_continuous(name="-log10 pval")+xlab("Phases")+
   scale_x_discrete( limits=c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24"))+ geom_hline(yintercept = 0.001, color="black", linetype="dashed")+
  ggtitle("Diurnal vs Nocturnal: All circproms")+theme_bw()+
  scale_color_manual(values=brewer.pal(9,"RdBu"), name="Phases")
#dev.off()



  ####################3#Introns
tempoclock_dayRP<-mean_exp_bothphases_randompickpromsINTRONS[,1:2]
rownames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_dayRP)<-c("0", "6")
tempoclock_nightRP<-mean_exp_bothphases_randompickpromsINTRONS[,3:4]
rownames(tempoclock_nightRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_nightRP)<-c("12", "18")

#For fishers test, prom anchors no proportions
t<-cbind(apply(tempoclock_dayRP, 1, mean), apply(tempoclock_nightRP, 1, mean), apply(mean_expbothphases_randompicknooscpromsINTRONS_control, 1, mean))
colnames(t)<-c("Diurnal", "Nocturnal", "Control")
#fisher.test(matrix(c(sum(t[1:4,1]), sum(t[5:8, 1]), sum(t[1:4, 3]), sum(t[5:8,3])), byrow = T, ncol=2))#p-value = 0.5471
#fisher.test(matrix(c(sum(t[1:4,2]), sum(t[5:8, 2]), sum(t[1:4, 3]), sum(t[5:8,3])), byrow = T, ncol=2))#p-value =  0.1304
#fisher.test(matrix(c(sum(t[1:4,1]), sum(t[5:8, 1]), sum(t[1:4, 2]), sum(t[5:8,2])),  byrow = T, ncol=2)) #the only one significant, p-value = 0.005327


wilcox.test(t[,1], t[,3], paired = T) #Diurnal vs control p-value = 0.007813
wilcox.test(t[,2], t[,3], paired = T) #Nocturnal vs control p-value =  0.01563
wilcox.test(t[,1], t[,2], paired = T) #Diurnal vs Nocturnal p-value = 0.3125

    ###----introns pvals
#1) Mean of diurnals ZT0 Y 6 for each permutation, I will have 100 values for each of the ZT ernas for the diurnal and noc Circproms
#2) With two distributions for each timepoint i can compare. In the end I will have a vector of 8 pvals 

zt0circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickpromsINTRONS[1,perm])[ernasphase]})})

zt6circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickpromsINTRONS[2,perm])[ernasphase]})})

zt12circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickpromsINTRONS[3,perm])[ernasphase]})})

zt18circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickpromsINTRONS[4,perm])[ernasphase]})})


#mean of diurnal and noc
diurcircproms<-(zt0circ+zt6circ)/2
noccircproms<-(zt12circ+zt18circ)/2

#Check if distributions are normal or not
nd<-sapply(c(1:4, 6:8), function(x){shapiro.test(diurcircproms[,x])[[2]]}) #remove col 5 because its an identical distribution
nn<-sapply(1:8, function(x){shapiro.test(noccircproms[,x])[[2]]})
par(mfrow=c(1,2))
plot(nn)+abline(h = .001, col="red")
plot(nd)+abline(h = .001, col="red")
#dev.off()
#normal: 1.7.8
#not normal: 2,3,4,5,6
pvalsintrons<-c(t.test(diurcircproms[,1],noccircproms[,1] )[[3]],
                wilcox.test(diurcircproms[,2],noccircproms[,2] )[[3]], 
                wilcox.test(diurcircproms[,3],noccircproms[,3] )[[3]], 
                wilcox.test(diurcircproms[,4],noccircproms[,4] )[[3]], 
                wilcox.test(diurcircproms[,5],noccircproms[,5] )[[3]], 
                wilcox.test(diurcircproms[,6],noccircproms[,6] )[[3]], 
                t.test(diurcircproms[,7],noccircproms[,7] )[[3]], 
                t.test(diurcircproms[,7],noccircproms[,7] )[[3]])
 
 names(pvalsintrons)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")

 #svglite::svglite("Analysis_per_phases/MFM_RNAseqrandompicked_daynightcontrol_pvals_diurvsnocINTRONS.svg")
ggplot(melt(log10(pvalsintrons)), aes(x=c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24"), y=value))+geom_point(size=4, aes(col=rownames(melt(log10(pvalsintrons)))))+  scale_y_continuous(name="-log10 pval")+xlab("Phases")+
    geom_point(shape = 1,size=4,colour = "gray49")+
   scale_x_discrete( limits=c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24"))+ geom_hline(yintercept = 0.001, color="black", linetype="dashed")+
  ggtitle("Diurnal vs Nocturnal: Introns")+theme_bw()+scale_color_manual(values=brewer.pal(9,"RdBu"), name="Phases")
#dev.off()

# ###############results were not st  sig, try only comparing the distrbution
# #proms anchor
 shapiro.test(t[,1]) #p-value = 0.1826
shapiro.test(t[,2]) #p-value = 0.008647
shapiro.test(t[,3]) #p-value = 0.08065
t.test(t[,1], t[,3]) #p-value = 0.1344
t.test(t[,2], t[,3]) #p-value = 0.1183
t.test(t[,1], t[,2])#p-value = 0.977
#ernas anchor
shapiro.test(t[,1]) #p-value = 0.655
shapiro.test(t[,2]) #p-value = 0.2955
shapiro.test(t[,3]) #p-value = 0.242
t.test(t[,1], t[,3]) #p-value = 0.1445
t.test(t[,2], t[,3]) #p-value = 0.202
t.test(t[,1], t[,2], )#p-value = 0.9469

##############Try using the proportions
#proms
zt0circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickproms[1,perm])[ernasphase]})})
zt6circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickproms[2,perm])[ernasphase]})})
zt12circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickproms[3,perm])[ernasphase]})})
zt18circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickproms[4,perm])[ernasphase]})})
diurcircproms<-(zt0circ+zt6circ)/2
noccircproms<-(zt12circ+zt18circ)/2
diurcircproms<-t(apply(diurcircproms, 1, function(x){x/sum(x)}))
noccircproms<-t(apply(noccircproms, 1, function(x){x/sum(x)}))

rownames(diurcircproms)<-NULL
colnames(diurcircproms)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
rownames(diurcircproms)<-rep("Diurnal", 100)
colnames(noccircproms)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
rownames(noccircproms)<-NULL
rownames(noccircproms)<-rep("Nocturnal", 100)
#Check if distributions are normal or not
nd<-sapply(c(1:8), function(x){shapiro.test(diurcircproms[,x])[[2]]}) 
nn<-sapply(1:8, function(x){shapiro.test(noccircproms[,x])[[2]]})
par(mfrow=c(1,2))
plot(nn)+abline(h = .001, col="red")
plot(nd)+abline(h = .001, col="red")
#dev.off()
 
#plot distributions
t<-apply(diurcircproms, 2, density)
par(mfrow=c(2,4))
lapply(t, plot)

t<-apply(noccircproms, 2, density)
par(mfrow=c(2,4))
lapply(t, plot)
#They all look normal

#Since they are normal, check their variances
#Check variances of distributions, 
vnd<-sapply(c(1:8), function(x){var.test(diurcircproms[,x],noccircproms[,x], alternative = "two.sided")$p.value}) 

plot(vnd)+abline(h = .001, col="red")
plot(nd)+abline(h = .001, col="red")
#dev.off()


#vnd>.001 if pval is >.001 I dont have enough evidence to reject h0 which asumes the variances are equivalent, so h0==h1

[1]  TRUE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE
#t.test: 1.3.7
#welch: 2,4,5,6,8
pvalsintrons<-c(t.test(diurcircproms[,1],noccircproms[,1], var.equal = T )[[3]],
                t.test(diurcircproms[,2],noccircproms[,2], var.equal = F  )[[3]], 
                t.test(diurcircproms[,3],noccircproms[,3], var.equal = T  )[[3]], 
                t.test(diurcircproms[,4],noccircproms[,4], var.equal = F )[[3]], 
                t.test(diurcircproms[,5],noccircproms[,5] , var.equal = F)[[3]], 
                t.test(diurcircproms[,6],noccircproms[,6], var.equal = F )[[3]], 
                t.test(diurcircproms[,7],noccircproms[,7], var.equal = T  )[[3]], 
                t.test(diurcircproms[,7],noccircproms[,7], var.equal = F )[[3]])
 
 names(pvalsintrons)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")

#svglite::svglite("Analysis_per_phases/MFM_RNAseqrandompicked_daynightcontrol_pvals_diurvsnocALLPROMS_proportions.svg")
  ggplot(melt(log10(pvalsintrons)), aes(x=c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24"), y=value))+geom_point(size=4, aes(col=rownames(melt(log10(pvalsintrons)))))+  scale_y_continuous(name="-log10 pval")+xlab("Phases")+
    geom_point(shape = 1,size=4,colour = "gray49")+
   scale_x_discrete( limits=c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24"))+ geom_hline(yintercept = 0.001, color="black", linetype="dashed")+
  ggtitle("Diurnal vs Nocturnal: All proms ")+theme_bw()+scale_color_manual(values=brewer.pal(9,"RdBu"), name="Phases")
#dev.off()
#////////
#introns
zt0circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickpromsINTRONS[1,perm])[ernasphase]})})
zt6circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickpromsINTRONS[2,perm])[ernasphase]})})
zt12circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickpromsINTRONS[3,perm])[ernasphase]})})
zt18circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickpromsINTRONS[4,perm])[ernasphase]})})
diurcircproms<-(zt0circ+zt6circ)/2
noccircproms<-(zt12circ+zt18circ)/2
diurcircproms<-t(apply(diurcircproms, 1, function(x){x/sum(x)}))
noccircproms<-t(apply(noccircproms, 1, function(x){x/sum(x)}))

rownames(diurcircproms)<-NULL
colnames(diurcircproms)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
rownames(diurcircproms)<-rep("Diurnal", 100)
colnames(noccircproms)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
rownames(noccircproms)<-NULL
rownames(noccircproms)<-rep("Nocturnal", 100)
#Check if distributions are normal or not
nd<-sapply(c(1:8), function(x){shapiro.test(diurcircproms[,x])[[2]]}) 
nn<-sapply(1:8, function(x){shapiro.test(noccircproms[,x])[[2]]})
par(mfrow=c(1,2))
plot(nn)+abline(h = .001, col="red")
plot(nd)+abline(h = .001, col="red")
nn>.001
[1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
nd>.001
[1]  TRUE FALSE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE

 
#plot distributions
t<-apply(diurcircproms, 2, density)
par(mfrow=c(2,4))
lapply(t, plot)

#Not all of them look normal
t<-apply(noccircproms, 2, density)
par(mfrow=c(2,4))
lapply(t, plot)
#They all look normal

#Since they are normal, check their variances
#Check variances of distributions, 
vnd<-sapply(c(1,3,5:8), function(x){var.test(diurcircproms[,x],noccircproms[,x], alternative = "two.sided")$p.value}) 
plot(vnd)+abline(h = .001, col="red")#All variances all different, so use welc

#vnd>.001 if pval is >.001 I dont have enough evidence to reject h0 which asumes the variances are equivalent, so h0==h1

[1] FALSE FALSE FALSE FALSE FALSE FALSE
#mann-whitney: 2,4
#welch: 1,3,5:8
pvalsintrons<-c(t.test(diurcircproms[,1],noccircproms[,1], var.equal = F )[[3]],
                wilcox.test(diurcircproms[,2],noccircproms[,2], paired = F  )[[3]], 
                t.test(diurcircproms[,3],noccircproms[,3],var.equal = F  )[[3]], 
                wilcox.test(diurcircproms[,4],noccircproms[,4],paired = F )[[3]], 
                t.test(diurcircproms[,5],noccircproms[,5] ,var.equal = F )[[3]], 
                t.test(diurcircproms[,6],noccircproms[,6], var.equal = F  )[[3]], 
                t.test(diurcircproms[,7],noccircproms[,7], var.equal = F   )[[3]], 
                t.test(diurcircproms[,7],noccircproms[,7], var.equal = F )[[3]])
 
 names(pvalsintrons)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")

#svglite::svglite("Analysis_per_phases/MFM_RNAseqrandompicked_daynightcontrol_pvals_diurvsnocINTRONS_proportions.svg")
ggplot(melt(log10(pvalsintrons)), aes(x=c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24"), y=value))+geom_point(size=4, aes(col=rownames(melt(log10(pvalsintrons)))))+  scale_y_continuous(name="-log10 pval")+xlab("Phases")+
    geom_point(shape = 1,size=4,colour = "gray49")+
   scale_x_discrete( limits=c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24"))+ geom_hline(yintercept = -3, color="black", linetype="dashed")+
  ggtitle("Diurnal vs Nocturnal: Introns ")+theme_bw()+scale_color_manual(values=brewer.pal(9,"RdBu"), name="Phases")
#dev.off()




#boxplots

#All proms
zt0circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickproms[1,perm])[ernasphase]})})
zt6circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickproms[2,perm])[ernasphase]})})
zt12circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickproms[3,perm])[ernasphase]})})
zt18circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickproms[4,perm])[ernasphase]})})
diurcircproms<-(zt0circ+zt6circ)/2
noccircproms<-(zt12circ+zt18circ)/2

rownames(diurcircproms)<-NULL
colnames(diurcircproms)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
rownames(diurcircproms)<-rep("Diurnal", 100)
colnames(noccircproms)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
rownames(noccircproms)<-NULL
rownames(noccircproms)<-rep("Nocturnal", 100)

#svglite::svglite("Differential_interactions/Analysis_per_phases/MFM_RNAseqrandompicked_daynigh_promehn_allproms_boxplotofDIFFERENTIALinteractions.svg")
#colnames(matrix.combined)<-rep(c("diurnal", "nocturnal"), 8)
ggplot(melt(rbind(diurcircproms, noccircproms)), aes(x=as.factor(Var2), y=value, fill=as.factor(Var1)))+geom_boxplot(notch = T, notchwidth = 0.5, coef=10)+
  ylab("Interactions")+
  scale_x_discrete(name="eRNAs phases")+theme_bw()+
  scale_fill_manual(values = c("brown3" , "steelblue", "brown3" , "steelblue", "brown3" , "steelblue", "brown3",   "steelblue", "brown3", "steelblue", "brown3", "steelblue", "brown3", "steelblue", "brown3",  "steelblue"),name="Promoter phases")+ggtitle("Prom-enh: All proms")
#dev.off()


#Introns
zt0circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickpromsINTRONS[1,perm])[ernasphase]})})
zt6circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickpromsINTRONS[2,perm])[ernasphase]})})
zt12circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickpromsINTRONS[3,perm])[ernasphase]})})
zt18circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickpromsINTRONS[4,perm])[ernasphase]})})
diurcircproms<-(zt0circ+zt6circ)/2
noccircproms<-(zt12circ+zt18circ)/2
rownames(diurcircproms)<-NULL
colnames(diurcircproms)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
rownames(diurcircproms)<-rep("Diurnal", 100)
colnames(noccircproms)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
rownames(noccircproms)<-NULL
rownames(noccircproms)<-rep("Nocturnal", 100)

#svglite::svglite("Analysis_per_phases/MFM_RNAseqrandompicked_daynigh_promehn_introns_boxplotofinteractions.svg")
#colnames(matrix.combined)<-rep(c("diurnal", "nocturnal"), 8)
ggplot(melt(rbind(diurcircproms, noccircproms)), aes(x=as.factor(Var2), y=value, fill=as.factor(Var1)))+geom_boxplot(notch = T, notchwidth = 0.5, coef=10)+
  ylab("Interactions")+
  scale_x_discrete(name="eRNAs phases")+theme_bw()+
  scale_fill_manual(values = c("brown3" , "steelblue", "brown3" , "steelblue", "brown3" , "steelblue", "brown3",   "steelblue", "brown3", "steelblue", "brown3", "steelblue", "brown3", "steelblue", "brown3",  "steelblue"),name="Promoter phases")+ggtitle("Prom-enh: Introns")
#dev.off()

#####################----------------------------------------------------------
#Proportions

#All proms
zt0circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickproms[1,perm])[ernasphase]})})
zt6circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickproms[2,perm])[ernasphase]})})
zt12circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickproms[3,perm])[ernasphase]})})
zt18circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickproms[4,perm])[ernasphase]})})
diurcircproms<-(zt0circ+zt6circ)/2
noccircproms<-(zt12circ+zt18circ)/2
diurcircproms<-t(apply(diurcircproms, 1, function(x){x/sum(x)}))
noccircproms<-t(apply(noccircproms, 1, function(x){x/sum(x)}))

rownames(diurcircproms)<-NULL
colnames(diurcircproms)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
rownames(diurcircproms)<-rep("Diurnal", 100)
colnames(noccircproms)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
rownames(noccircproms)<-NULL
rownames(noccircproms)<-rep("Nocturnal", 100)

#svglite::svglite("Differential_interactions/Analysis_per_phases/MFM_RNAseqrandompicked_daynigh_promehn_allproms_boxplotofinteractions_proportions.svg")
#colnames(matrix.combined)<-rep(c("diurnal", "nocturnal"), 8)
ggplot(melt(rbind(diurcircproms, noccircproms)), aes(x=as.factor(Var2), y=value, fill=as.factor(Var1)))+geom_boxplot(notch = T, notchwidth = 0.5, coef=10)+
  ylab("Interactions")+
  scale_x_discrete(name="eRNAs phases")+theme_bw()+
  scale_fill_manual(values = c("brown3" , "steelblue", "brown3" , "steelblue", "brown3" , "steelblue", "brown3",   "steelblue", "brown3", "steelblue", "brown3", "steelblue", "brown3", "steelblue", "brown3",  "steelblue"),name="Promoter phases")+ggtitle("Proportions Prom-enh: All proms")
#dev.off()

#introns
zt0circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickpromsINTRONS[1,perm])[ernasphase]})})
zt6circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickpromsINTRONS[2,perm])[ernasphase]})})
zt12circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickpromsINTRONS[3,perm])[ernasphase]})})
zt18circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickpromsINTRONS[4,perm])[ernasphase]})})
diurcircproms<-(zt0circ+zt6circ)/2
noccircproms<-(zt12circ+zt18circ)/2
diurcircproms<-t(apply(diurcircproms, 1, function(x){x/sum(x)}))
noccircproms<-t(apply(noccircproms, 1, function(x){x/sum(x)}))

rownames(diurcircproms)<-NULL
colnames(diurcircproms)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
rownames(diurcircproms)<-rep("Diurnal", 100)
colnames(noccircproms)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
rownames(noccircproms)<-NULL
rownames(noccircproms)<-rep("Nocturnal", 100)

#svglite::svglite("Analysis_per_phases/MFM_RNAseqrandompicked_daynigh_promehn_introns_boxplotofinteractions_proportions.svg")
#colnames(matrix.combined)<-rep(c("diurnal", "nocturnal"), 8)
ggplot(melt(rbind(diurcircproms, noccircproms)), aes(x=as.factor(Var2), y=value, fill=as.factor(Var1)))+geom_boxplot(notch = T, notchwidth = 0.5, coef=10)+
  ylab("% Interactions")+
  scale_x_discrete(name="eRNAs phases")+theme_bw()+
  scale_fill_manual(values = c("brown3" , "steelblue", "brown3" , "steelblue", "brown3" , "steelblue", "brown3",   "steelblue", "brown3", "steelblue", "brown3", "steelblue", "brown3", "steelblue", "brown3",  "steelblue"),name="Promoter phases")+ggtitle("Proportions Prom-enh: Introns")
#dev.off()
