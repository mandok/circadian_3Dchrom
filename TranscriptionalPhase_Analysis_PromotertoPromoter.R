#######################################################

#Repeat analysis only for p-p interactions


################################  
#### Open chic merge file s#####
################################  
#I use the symmetrical set
chic<-read.csv("/inputs_for_scripts/All_step2_washU_text_1.txt", sep = " ", header = F)
#chic_baits<-read.csv("chic_baits_mm9.txt", sep = "\t", header = F)
chic_baits<-chic[,1:3]
#chic_otherends<-read.csv("chic_otherends_mm9.txt", sep = "\t", header = F)
chic_otherends<-chic[,4:6]
#Same number of lines! Allows mapping bait<->otherend

#Bedfiles, GRanges Objects
chic_otherends_bed<-GRanges(seqnames= Rle(chic_otherends[,1]), ranges = IRanges(chic_otherends[,2], chic_otherends[,3]))
chic_bait_bed<-GRanges(seqnames= Rle(chic_baits[,1]), ranges = IRanges(chic_baits[,2], chic_baits[,3]))

#--------------------------------------------------------------------
##Overlap of circproms to baits
#skip

#Open circadian promoters

fang_mrnas<-read.csv("/Users/andoku01/Dropbox/Masami/inputs_for_scripts/circproms/MFM_RNAseq/HindIIIfragments_circadiangenes_MFMRNAseq_phases.bed", sep = "\t", header = F)
fang_mrnas<-bedfile(fang_mrnas, columnnames)
circpromlist<-list("fang_mrnas")




#cirproms


circpromphases<-read.csv("/Users/andoku01/Dropbox/Masami/inputs_for_scripts/circproms/MFM_RNAseq/HindIIIfragments_circadiangenes_MFMRNAseq_phases.bed", sep = "\t", header = F)
circpromphasesINTRONS<-read.csv("/Users/andoku01/Dropbox/Masami/inputs_for_scripts/circproms/MFM_RNAseq/HindIIIfragments_circadiangenesonlyintrons_MFMRNAseq_phases.bed", sep = "\t", header = F)
#Divide per phases
#List of the phases
#Divide per phases
#List of the phases
#phases<-unique(circpromphases[,17])
#SORT FILE ACCORDING TO PHASES
circpromphases<-circpromphases[order(circpromphases[,5]),]
phases<-unique(circpromphases[,5])
#introns
circpromphasesINTRONS<-circpromphasesINTRONS[order(circpromphasesINTRONS[,5]),]
#phases_bin<-split(phases, ceiling(seq_along(phases)/3))
#Add separaterly last phase, otherwise only last bin would have 1 phase
#phases_bin[[8]]<-c(phases_bin[[8]], 23.5)
#phases_bin[[9]]<-NULL
#Divide the circprom per phases, create separe variables
#circprom_divided_list<-lapply(phases, function(x){
#  tempo<-circpromphases[circpromphases[,17]==x,c(1,17)]
#  return(tempo)
#  })


circprom_divided_list<-lapply(phases, function(x){
  tempo1<-circpromphases[circpromphases[,5]==x[1],]
  tempo2<-circpromphases[circpromphases[,5]==x[2],]
  tempo3<-circpromphases[circpromphases[,5]==x[3],]
  tempo4<-circpromphases[circpromphases[,5]==x[5],]
  tempo5<-rbind(tempo1, tempo2, tempo3, tempo4)
  return(tempo5)
})

#Introns
circpromINTRONS_divided_list<-lapply(phases, function(x){
  tempo1<-circpromphasesINTRONS[circpromphasesINTRONS[,5]==x[1],]
  tempo2<-circpromphasesINTRONS[circpromphasesINTRONS[,5]==x[2],]
  tempo3<-circpromphasesINTRONS[circpromphasesINTRONS[,5]==x[3],]
  tempo4<-circpromphasesINTRONS[circpromphasesINTRONS[,5]==x[5],]
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


#8 items in list with the ernas per phase (each 3 hrs)
#names for the new phases
#phases_bin_names<-sapply(1:8, function(x){paste(phases_bin[[x]][1], phases_bin[[x]][2], phases_bin[[x]][3], phases_bin[[x]][4], sep = "_")})

#phases_bin_names<-gsub("_NA", "", phases_bin_names)
names(circprom_divided_list)<-gsub(pattern = "ZT", replacement = "" , phases)
createvars(circprom_divided_list, "df_circprom", gsub(pattern = "ZT", replacement = "" , phases))


names(circpromINTRONS_divided_list)<-gsub(pattern = "ZT", replacement = "" , phases)
createvars(circpromINTRONS_divided_list, "df_INTRONScircprom", gsub(pattern = "ZT", replacement = "" , phases))

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

createvars(list_circpromsperphase_bed, "bed_circprom",  gsub(pattern = "ZT", replacement = "" , phases))

createvars(list_circpromsperphaseINTRONS_bed, "bed_intronscircprom",  gsub(pattern = "ZT", replacement = "" , phases))




#Overlap with CHiC

####################Overlap with ChiC
circpromlist_phases<-sapply(mixedsort(ls(pattern = "bed_circprom", sorted = F), decreasing = F), function(x){list(x)})
    #intron
circpromlistintrons_phases<-sapply(mixedsort(ls(pattern = "bed_intronscircprom", sorted = F), decreasing = F), function(x){list(x)})


#####################Retrieve index from bait, proms per phases
phases_circprom_index_baits<-lapply(circpromlist_phases, function(x){
  tempo<-retrieveindex_fromchic(get(x), chic_otherends_bed, chic_bait_bed)
})

names(phases_circprom_index_baits) <- paste("BI_proms", phases, sep = "")#phases MFM and phases_bin_names for fang
list2env(phases_circprom_index_baits , envir = .GlobalEnv)
#remove(phases_eRNAs_index_baits)

    #introns
phases_circpromintrons_index_baits<-lapply(circpromlistintrons_phases, function(x){
  tempo<-retrieveindex_fromchic(get(x), chic_otherends_bed, chic_bait_bed)
})

names(phases_circpromintrons_index_baits) <- paste("BI_intronproms", phases, sep = "")
list2env(phases_circpromintrons_index_baits , envir = .GlobalEnv)
#remove(phases_eRNAs_index_baits)




###### Both circproms and circprom per phase analysis



#Compare each circadian promoter vs all the ernas_phases,
#Rerun chunk in line 104
bothphases_circadianprom_counts<-
  lapply(circpromlist_phases, function(y){# each circprom
    tempo_circ<-get(y);
    t<-lapply(phases_circprom_index_baits, function(tempo_circ2){
    tempo_circ2<-length(subjectHits(findOverlaps(tempo_circ2, tempo_circ)))
    })
    
})


names(bothphases_circadianprom_counts)<-circpromlist_phases


    #introns
bothphases_circadianpromintron_counts<-
  lapply(circpromlistintrons_phases, function(y){# each circprom
    tempo_circ<-get(y);
    t<-lapply(phases_circpromintrons_index_baits, function(tempo_circ2){
    tempo_circ2<-length(subjectHits(findOverlaps(tempo_circ2, tempo_circ)))
    })
    
})


names(bothphases_circadianpromintron_counts)<-circpromlistintrons_phases



#Plot as clocks

#Plot as clocks, proms as anchors


###########DAY
tempoclock<-sapply(1:2, function(x){
  unlist(bothphases_circadianprom_counts[[x]])
})
#tempoclock_day<-apply(tempoclock, 1, function(y){ mean(y)})

#temporal<-cbind(1:8, "Value"=tempoclock_day)
 # temporal<-as.data.frame(temporal)
#  temporal[,2]<-(temporal[,2]/sum(temporal[,2]))#as.numeric(levels(temporal[,#2]))[temporal[,2]]
  
  
d<-ggplot(melt(sapply(1:2, function(x){tempoclock[,x]/apply(tempoclock, 2, sum)[x]})), aes(x=Var1, y=value, fill=(Var2)))+geom_bar(width = 1, stat = "identity")+ 
  scale_y_discrete(limits="")+
  theme( panel.background = element_blank(), axis.text=element_text(size=20, vjust=0.8, face = "bold"), axis.line.x = element_line(color="white", size = 1), axis.line.y = element_line(color="white", size = 1), axis.title=element_text(size=20))+
  scale_fill_distiller(palette = "OrRd", breaks=c(1,2), labels=c(1,2), name="Phases")+
  scale_x_discrete(labels= c("0", "6", "12", "18"))+
  ylab("")+xlab("Diurnal promoters")+
  coord_polar() 
d

 
  ###########NIGHT
tempoclock<-sapply(3:4, function(x){
  unlist(bothphases_circadianprom_counts[[x]])
})
#tempoclock_day<-apply(tempoclock, 1, function(y){  mean(y)})

#temporal<-cbind(1:8, "Value"=tempoclock_day)temporal<-as.data.frame(temporal)temporal[,2]<-(temporal[,2]/sum(temporal[,2]))#as.numeric(levels(temporal[,2]))[temporal[,2]]
n<-ggplot(melt(sapply(1:2, function(x){tempoclock[,x]/apply(tempoclock, 2, sum)[x]})), aes(x=Var1, y=value, fill=Var2))+geom_bar(width = 1, stat = "identity")+
  scale_y_discrete(limits="")+
  theme( panel.background = element_blank(), axis.text=element_text(size=20, vjust=0.8, face = "bold"), axis.line.x = element_line(color="white", size = 1),axis.line.y = element_line(color="white", size = 1), axis.title=element_text(size=20))+
    scale_fill_distiller(palette="PuBu", breaks=c(1,2), labels=c(1,2),name="Phases")+ scale_x_discrete(labels= c("0", "6", "12", "18"))+
  ylab("")+xlab("Nocturnal promoters")+
  coord_polar() 
n
tiff("clocks_day_circpromsB.tiff", height = 12, width = 12, units = 'cm', 
     compression = "lzw", res = 300)
d
#dev.off()

tiff("clocks_night_circpromsB.tiff", height = 12, width = 12, units = 'cm', 
     compression = "lzw", res = 300)
n
#dev.off()

    #introns
###########DAY
tempoclock<-sapply(1:2, function(x){
  unlist(bothphases_circadianpromintron_counts[[x]])
})
#tempoclock_day<-apply(tempoclock, 1, function(y){ mean(y)})

#temporal<-cbind(1:8, "Value"=tempoclock_day)
 # temporal<-as.data.frame(temporal)
#  temporal[,2]<-(temporal[,2]/sum(temporal[,2]))#as.numeric(levels(temporal[,#2]))[temporal[,2]]
  
  
d<-ggplot(melt(sapply(1:2, function(x){tempoclock[,x]/apply(tempoclock, 2, sum)[x]})), aes(x=Var1, y=value, fill=(Var2)))+geom_bar(width = 1, stat = "identity")+ 
  scale_y_discrete(limits="")+
  theme( panel.background = element_blank(), axis.text=element_text(size=20, vjust=0.8, face = "bold"), axis.line.x = element_line(color="white", size = 1), axis.line.y = element_line(color="white", size = 1), axis.title=element_text(size=20))+
  scale_fill_distiller(palette = "OrRd", breaks=c(1,2), labels=c(1,2), name="Phases")+
  scale_x_discrete(labels= c("0", "6", "12", "18"))+
  ylab("")+xlab("Diurnal promoters")+
  coord_polar() 
d

 
  ###########NIGHT
tempoclock<-sapply(3:4, function(x){
  unlist(bothphases_circadianpromintron_counts[[x]])
})
#tempoclock_day<-apply(tempoclock, 1, function(y){  mean(y)})

#temporal<-cbind(1:8, "Value"=tempoclock_day)temporal<-as.data.frame(temporal)temporal[,2]<-(temporal[,2]/sum(temporal[,2]))#as.numeric(levels(temporal[,2]))[temporal[,2]]
n<-ggplot(melt(sapply(1:2, function(x){tempoclock[,x]/apply(tempoclock, 2, sum)[x]})), aes(x=Var1, y=value, fill=Var2))+geom_bar(width = 1, stat = "identity")+
  scale_y_discrete(limits="")+
  theme( panel.background = element_blank(), axis.text=element_text(size=20, vjust=0.8, face = "bold"), axis.line.x = element_line(color="white", size = 1),axis.line.y = element_line(color="white", size = 1), axis.title=element_text(size=20))+
    scale_fill_distiller(palette="PuBu", breaks=c(1,2), labels=c(1,2),name="Phases")+ scale_x_discrete(labels= c("0", "6", "12", "18"))+
  ylab("")+xlab("Nocturnal promoters")+
  coord_polar() 
n
#tiff("clocks_day_circpromsintronsB.tiff", height = 12, width = 12, units = 'cm', 
     compression = "lzw", res = 300)
d
#dev.off()

#tiff("clocks_night_circpromsintronsB.tiff", height = 12, width = 12, units = 'cm', 
     compression = "lzw", res = 300)
n
#dev.off()


#######################################
#2018-09-09 using MFM RNAseq proms
#Random pick, p-p; 


#------------
#RANDOM ERNAS, now try random proms
#Extra: EXPECTED
#Repeat 100 times
#Take a random list everytime; use general list of genes and introns 
ppgenelist<-read.csv("/inputs_for_scripts/circproms/MFM_RNAseq/HindIIIfragments_circadiangenes_MFMRNAseq_phases.bed", header = F, sep="\t")
ppgenelist<-read.csv("/inputs_for_scripts/circproms/Updated/HindIIIfragments_circadiangenes_fang_phases.bed", header = F, sep="\t")
ppgenelist<-ppgenelist[,-4]
ppintronslist<-read.csv("/inputs_for_scripts/circproms/MFM_RNAseq/HindIIIfragments_circadiangenesonlyintrons_MFMRNAseq_phases.bed", header = F, sep="\t")

#Divide per phases
#List of the phases
#phases<-unique(circpromphases[,17])
#SORT FILE ACCORDING TO PHASES
ppgenelist<-ppgenelist[order(ppgenelist[,5]),]
ppintronslist<-ppintronslist[order(ppintronslist[,5]),]

phases<-unique(ppintronslist[,5])

phases_bin<-split(unique(ppgenelist[,5]),ceiling(seq_along(unique(ppgenelist[,5]))/3))#For Fang
phases_bin<-phases#For MFM-RNAseq
#Add separaterly last phase, otherwise only last bin would have 1 phase
phases_bin[[8]]<-c(phases_bin[[8]], 23.5)
phases_bin[[9]]<-NULL
phases_bin_names<-sapply(1:8, function(x){paste(phases_bin[[x]][1], phases_bin[[x]][2], phases_bin[[x]][3], phases_bin[[x]][4], sep = "_")})

phases_bin_names<-gsub("_NA", "", phases_bin_names)

ppcircprom_divided_list<-lapply(phases_bin, function(x){ #for fang phases_bin; phases for MFM
  tempo1<-ppgenelist[ppgenelist[,5]==x[1],]
  tempo2<-ppgenelist[ppgenelist[,5]==x[2],]
  tempo3<-ppgenelist[ppgenelist[,5]==x[3],]
  tempo4<-ppgenelist[ppgenelist[,5]==x[4],] # x[5]
  tempo5<-rbind(tempo1, tempo2, tempo3, tempo4) #Add tempo4 for fang
  return(tempo5)
})

ppintroncircprom_divided_list<-lapply(phases, function(x){
  tempo1<-ppintronslist[ppintronslist[,5]==x[1],]
  tempo2<-ppintronslist[ppintronslist[,5]==x[2],]
  tempo3<-ppintronslist[ppintronslist[,5]==x[3],]
  tempo4<-ppintronslist[ppintronslist[,5]==x[4],]
  tempo5<-rbind(tempo1, tempo2, tempo3, tempo4)
  return(tempo5)
})


#Remove NAs
ppcircprom_divided_list<-lapply(ppcircprom_divided_list, function(x){
  t<-x
  t<-t[complete.cases(t),]
})

ppintroncircprom_divided_list<-lapply(ppintroncircprom_divided_list, function(x){
  t<-x
  t<-t[complete.cases(t),]
})


#8 items in list with the ernas per phase (each 3 hrs)
#names for the new phases
#phases_bin_names<-sapply(1:8, function(x){paste(phases_bin[[x]][1], phases_bin[[x]][2], phases_bin[[x]][3], phases_bin[[x]][4], sep = "_")})

#phases_bin_names<-gsub("_NA", "", phases_bin_names)
names(ppcircprom_divided_list)<-gsub(pattern = "ZT", replacement = "" , phases_bin_names) #phases_bin_names Fang; phases MFM
createvars(ppcircprom_divided_list, "df_ppcircprom", gsub(pattern = "ZT", replacement = "" , phases_bin_names))#phases_bin_names for fang; phases MFM
#For fang; run ver 3.0 to separate per phases
names(phases)<- phases_bin_names #phases
createvars(circprom_divided_list, "df_ppcircprom", phases_bin_names)


names(ppintroncircprom_divided_list)<-gsub(pattern = "ZT", replacement = "" , phases)
createvars(ppintroncircprom_divided_list, "df_ppintronscircprom", gsub(pattern = "ZT", replacement = "" , phases))

#Transform into granges object
#biomartpos_expr replace for 
list_ppcircpromsperphase_bed<-lapply(sapply(mixedsort(ls(pattern = "df_ppcircprom")), list), function(x){
  t<-get(x) 
  t1<-GRanges(seqnames= Rle(t[,1]), ranges = IRanges(t[,2], t[,3]), phase=t[,5])
  return(t1)
})
# MFM: lapply(list_ppcircpromsperphase_bed, length) $df_ppcircprom0: 208; $df_ppcircprom6: 322; $df_ppcircprom12: 472; $df_ppcircprom18: 193
#FANG: lapply(list_ppcircpromsperphase_bed, length) $df_ppcircprom1_2.5_3.25: 165; $df_ppcircprom4_5.5_6.25: 111; $df_ppcircprom7_8.5_9.25: 114; $df_ppcircprom10_10.75_11.5: 103; $df_ppcircprom12.25_13_13.75: 107; $df_ppcircprom14.5_15.25_16: 175; $df_ppcircprom17.5_18.25_19: 130; $df_ppcircprom20.5_21.25_22_23.5: 193

list_ppintroncircpromsperphase_bed<-lapply(sapply(mixedsort(ls(pattern = "df_ppintronscircprom")), list), function(x){
  t<-get(x) 
  t1<-GRanges(seqnames= Rle(t[,1]), ranges = IRanges(t[,2], t[,3]), phase=t[,5])
  return(t1)
})
#lapply(list_ppintroncircpromsperphase_bed, length)  $df_ppcircprom0: 58; $df_ppcircprom6: 32; $df_ppcircprom12: 128; $df_ppcircprom18: 53


#createvars(list_ppcircpromsperphase_bed, "bed_ppcircprom",  gsub(pattern = "ZT", replacement = "" , phases))
#createvars(list_ppcircpromsperphase_bed, "bed_ppcircprom", phases_bin_names) for fang

#@@@@@@@@@@@@@@@
 expbothphases_randompickcircproms1<-sapply(1:100, function(rep){
  RANlist_circpromsperphase_bed<-sapply(1:8, function(phase){
    pos<-sample(list_ppcircpromsperphase_bed[[phase]] ,103, replace = F)#193 for MFM and 103 for fang; 
    })
                #lapply(list_ernasperphase_bed, length)[[phase]], replace = F)
#Random no osc proms
  #RANlist_circpromperphase_bed<-sapply(1:8, function(phase){
  #pos<-sample(list_circpromsperphase_bed[[phase]],109, replace = F)
  #})
#Retrieve indexes of circproms
  RANindexes_baitsofcircproms_perphases<-lapply(RANlist_circpromsperphase_bed, function(x){
  tempo<-retrieveindex_fromchic(x, chic_otherends_bed, chic_bait_bed)
  })
  
  #Bait(index ernas) overlap with proms

  RANbothphases_circadianprom_circadianproms_counts<-
  lapply(list_ppcircpromsperphase_bed, function(tempo_circ){# each circprom
    t<-lapply(RANindexes_baitsofcircproms_perphases, function(tempo_ernas){
    tempo_ernas<-length(subjectHits(findOverlaps(tempo_ernas, tempo_circ)))
    })
})
  names(RANbothphases_circadianprom_circadianproms_counts)<-phases_bin_names #phases_bin_names for fang, phases for MFM
  
  
  #names(RANoverlaps_and_correlations_bothphases_circadianprom_ernas)<-circpromlist_phases
  return(RANbothphases_circadianprom_circadianproms_counts)

})

    #introns
expbothphases_randompickcircpromsintrons1<-sapply(1:100, function(rep){
  RANlist_circpromsperphase_bed<-sapply(1:4, function(phase){
    pos<-sample(list_ppintroncircpromsperphase_bed[[phase]] ,32, replace = F)
    })
                #lapply(list_ernasperphase_bed, length)[[phase]], replace = F)
#Random no osc proms
  #RANlist_circpromperphase_bed<-sapply(1:8, function(phase){
  #pos<-sample(list_circpromsperphase_bed[[phase]],109, replace = F)
  #})
#Retrieve indexes of circproms
  RANindexes_baitsofcircproms_perphases<-lapply(RANlist_circpromsperphase_bed, function(x){
  tempo<-retrieveindex_fromchic(x, chic_otherends_bed, chic_bait_bed)
  })
  
  #Bait(index ernas) overlap with proms

  RANbothphases_circadianprom_circadianproms_counts<-
  lapply(list_ppintroncircpromsperphase_bed, function(tempo_circ){# each circprom
    t<-lapply(RANindexes_baitsofcircproms_perphases, function(tempo_ernas){
    tempo_ernas<-length(subjectHits(findOverlaps(tempo_ernas, tempo_circ)))
    })
})
  names(RANbothphases_circadianprom_circadianproms_counts)<-phases
  
  
  #names(RANoverlaps_and_correlations_bothphases_circadianprom_ernas)<-circpromlist_phases
  return(RANbothphases_circadianprom_circadianproms_counts)

})


#Retrieve indexes of eRNAs
  RANindexes_baitsofcircproms_perphases<-lapply(list_circpromsperphase_bed, function(x){
  tempo<-retrieveindex_fromchic(x, chic_otherends_bed, chic_bait_bed)
  })
  
  #Bait(index ernas) overlap with proms

  RANbothphases_circadianprom_circadianproms_counts<-
  lapply(RANlist_circpromperphase_bed, function(tempo_circ){# each circprom
    t<-lapply(RANindexes_baitsofcircproms_perphases, function(tempo_ernas){
    tempo_ernas<-length(subjectHits(findOverlaps(tempo_ernas, tempo_circ)))
    })
})
  names(RANbothphases_circadianprom_circadianproms_counts)<-circpromlist_phases
  
  
  #names(RANoverlaps_and_correlations_bothphases_circadianprom_ernas)<-circpromlist_phases
  return(RANbothphases_circadianprom_circadianproms_counts)

})

#Obtain the mean, 
mean_expbothphases_randompickcircproms1<-sapply(1:8, function(cirpromphase){sapply(1:8, function(circpromphase1){mean(unlist(sapply(1:100, function(rep){expbothphases_randompickcircproms1[cirpromphase,][[rep]][circpromphase1]})))})})

mean_expbothphases_randompickcircpromsintrons1<-sapply(1:4, function(cirpromphase){sapply(1:4, function(circpromphase1){mean(unlist(sapply(1:100, function(rep){expbothphases_randompickcircpromsintrons1[cirpromphase,][[rep]][circpromphase1]})))})})

#Skip
#mean_expbothphases_randompickcircproms2<-sapply(1:8, function(circpromphase1){sapply(1:8, function( cirpromphase){mean(unlist(sapply(1:100, function(rep){expbothphases_randompickcircproms2[cirpromphase,][[rep]][circpromphase1]})))})})

mean_expbothphases_randompickcircproms1
mean_expbothphases_randompickcircpromsintrons1
#mean_expbothphases_randompickcircproms2
#------------
#Plot
###DAY
tempoclock<-sapply(1:2, function(x){
  unlist(expbothphases_randompickcircproms1[[x]])
})
tempoclock_dayRP<-mean_expbothphases_randompickcircproms1[,1:2]
rownames(tempoclock_dayRP)<-c("0", "6", "12", "18")
colnames(tempoclock_dayRP)<-c("0", "6")


###NIGHT
tempoclock<-sapply(3:4, function(x){
  unlist(expbothphases_randompickcircproms1[[x]])
})
tempoclock_nightRP<-mean_expbothphases_randompickcircproms1[,3:4]
rownames(tempoclock_nightRP)<-c("0", "6", "12", "18")
colnames(tempoclock_nightRP)<-c("12", "18")

#Skip and go to 2274
####Plot as proportions

#ernas anchor
t<-sapply(1:4, function(x){tempoclock_dayRP[,x]/apply(tempoclock_dayRP, 2, sum)[x]})
colnames(t)<-c("1-3", "4-6", "7-9", "10-12")
p1<-ggplot(melt(t), aes(x=Var1, y=value, fill=Var2))+geom_bar(width = 1, stat = "identity")+
  theme( panel.background = element_blank(), axis.text=element_text(size=12, vjust=0.8), axis.line.x = element_line(color="white", size = 1), axis.text.x = element_text(angle = 90, hjust = 1) ,axis.line.y = element_line(color="white", size = 1), axis.title=element_text(size=12))+
  scale_fill_brewer(palette="OrRd", name="Phases")+
  ylab("")+xlab("Promoter Phases")+ggtitle("Diurnal eRNAs")


t<-sapply(1:4, function(x){tempoclock_nightRP[,x]/apply(tempoclock_nightRP, 2, sum)[x]})
colnames(t)<-c("13-15", "16-18", "19-21", "22-24")

p2<-ggplot(melt(t), aes(x=Var1, y=value, fill=Var2))+geom_bar(width = 1, stat = "identity")+
  theme( panel.background = element_blank(), axis.text=element_text(size=12, vjust=0.8), axis.line.x = element_line(color="white", size = 1),axis.text.x = element_text(angle = 90, hjust = 1),axis.line.y = element_line(color="white", size = 1), axis.title=element_text(size=12))+
  scale_fill_brewer(palette="PuBu", name="Phases")+
  ylab("")+xlab("Promoter Phases")+ggtitle("Nocturnal eRNAs")

#Proms anchors
t<-sapply(1:4, function(x){tempoclock_dayRP[,x]/apply(tempoclock_dayRP, 2, sum)[x]})
colnames(t)<-c("1-3", "4-6", "7-9", "10-12")
p3<-ggplot(melt(t), aes(x=Var1, y=value, fill=Var2))+geom_bar(width = 1, stat = "identity")+
  theme( panel.background = element_blank(), axis.text=element_text(size=12, vjust=0.8), axis.line.x = element_line(color="white", size = 1),axis.text.x = element_text(angle = 90, hjust = 1),axis.line.y = element_line(color="white", size = 1), axis.title=element_text(size=12))+
    scale_fill_brewer(palette="OrRd", name="Phases")+
  #scale_x_discrete(labels= c("1", "4", "7", "10", "13", "16", "19", "22"))+
  ylab("")+xlab("eRNAs Phases")+ggtitle("Diurnal Promoters")


t<-sapply(1:4, function(x){tempoclock_nightRP[,x]/apply(tempoclock_nightRP, 2, sum)[x]})
colnames(t)<-c("13-15", "16-18", "19-21", "22-24")
p4<-ggplot(melt(t), aes(x=Var1, y=value, fill=Var2))+geom_bar(width = 1, stat = "identity")+
  theme( panel.background = element_blank(), axis.text=element_text(size=12, vjust=0.8), axis.line.x = element_line(color="white", size = 1),axis.text.x = element_text(angle = 90, hjust = 1),axis.line.y = element_line(color="white", size = 1), axis.title=element_text(size=12))+
    scale_fill_brewer(palette="PuBu", name="Phases")+
  #scale_x_discrete(labels= c("1", "4", "7", "10", "13", "16", "19", "22"))+
  ylab("")+xlab("eRNAs Phases")+ggtitle("Nocturnal Promoters")

grid.arrange(p1, p2, p3, p4)

###Control randompicked, use no circ proms
onlybaits<-bedfile(read.csv("/Users/andoku01/Dropbox/Masami/Baits_Gen_Res/Baits.txt", header = T, sep="\t"), columnnames)
t<-chic_bait_bed[unique(queryHits(findOverlaps(chic_bait_bed, onlybaits)))]
noncircproms<-t[!(1:247943 %in% queryHits(findOverlaps(t,bedfile(circpromphases))))]


expbothphases_randompicknooscproms_controlPP<-sapply(1:100, function(rep){
  #RANlist_ernasperphase_bed<-sapply(1:8, function(phase){
   # pos<-sample(noosc_ernas ,160, replace = F)
    #})
                #lapply(list_ernasperphase_bed, length)[[phase]], replace = F)
#Random no osc proms
  RANlist_circpromperphase_bed<-sapply(1:8, function(phase){#1:4 mfm and 1:8
  pos<-sample(unique(noncircproms), 103, replace=F)})#68 for mfm and 103 for fang
    #sample(list_circpromsperphase_bed[[phase]],109, replace = F)})
#Retrieve indexes of eRNAs
  RANindexes_baitsofernas_perphases<-lapply(list_ppcircpromsperphase_bed, function(x){
  tempo<-retrieveindex_fromchic(x, chic_otherends_bed, chic_bait_bed)
  })
  
  #Bait(index ernas) overlap with proms

  RANbothphases_circadianprom_ernas_counts<-
  lapply(RANlist_circpromperphase_bed, function(tempo_circ){# each circprom
    t<-lapply(RANindexes_baitsofernas_perphases, function(tempo_ernas){
    tempo_ernas<-length(subjectHits(findOverlaps(tempo_ernas, tempo_circ)))
    })
})
  names(RANbothphases_circadianprom_ernas_counts)<-phases_bin_names #phases_bin_names for fang and phases for MFM
  
  
  #names(RANoverlaps_and_correlations_bothphases_circadianprom_ernas)<-circpromlist_phases
  return(RANbothphases_circadianprom_ernas_counts)

})

    #introns
expbothphases_randompicknooscpromsintrons_controlPP<-sapply(1:100, function(rep){
  #RANlist_ernasperphase_bed<-sapply(1:8, function(phase){
   # pos<-sample(noosc_ernas ,160, replace = F)
    #})
                #lapply(list_ernasperphase_bed, length)[[phase]], replace = F)
#Random no osc proms
  RANlist_circpromperphase_bed<-sapply(1:4, function(phase){
  pos<-sample(unique(noncircproms), 32, replace=F)})
    #sample(list_circpromsperphase_bed[[phase]],109, replace = F)})
#Retrieve indexes of eRNAs
  RANindexes_baitsofernas_perphases<-lapply(list_ppintroncircpromsperphase_bed, function(x){
  tempo<-retrieveindex_fromchic(x, chic_otherends_bed, chic_bait_bed)
  })
  
  #Bait(index ernas) overlap with proms

  RANbothphases_circadianprom_ernas_counts<-
  lapply(RANlist_circpromperphase_bed, function(tempo_circ){# each circprom
    t<-lapply(RANindexes_baitsofernas_perphases, function(tempo_ernas){
    tempo_ernas<-length(subjectHits(findOverlaps(tempo_ernas, tempo_circ)))
    })
})
  names(RANbothphases_circadianprom_ernas_counts)<-phases
  
  
  #names(RANoverlaps_and_correlations_bothphases_circadianprom_ernas)<-circpromlist_phases
  return(RANbothphases_circadianprom_ernas_counts)

})

#Obtain the mean, 
mean_exp_bothphases_randompicknooscproms_controlPP<-sapply(1:8, function(cirpromphase){sapply(1:8, function(cirpromphase2){mean(unlist(sapply(1:100, function(rep){expbothphases_randompicknooscproms_controlPP[cirpromphase,][[rep]][cirpromphase2]})))})})

mean_exp_bothphases_randompicknooscpromsintrons_controlPP<-sapply(1:4, function(cirpromphase){sapply(1:4, function(cirpromphase2){mean(unlist(sapply(1:100, function(rep){expbothphases_randompicknooscpromsintrons_controlPP[cirpromphase,][[rep]][cirpromphase2]})))})})

mean_exp_bothphases_randompicknooscproms_controlPP
mean_exp_bothphases_randompicknooscpromsintrons_controlPP

tempoclock_nightRP_control<-mean_exp_bothphases_randompicknooscproms_controlPP[,5:8]
rownames(tempoclock_nightRP_control)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
#c("0", "6", "12", "18") for MFM and c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24") for fang
colnames(tempoclock_nightRP_control)<-c( "13-15", "16-18", "19-21", "22-24")

tempoclock_dayRP_control<-mean_exp_bothphases_randompicknooscproms_controlPP[,1:4]
rownames(tempoclock_dayRP_control)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")#c("0", "6", "12", "18")for MFM and c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24") for fang
colnames(tempoclock_dayRP_control)<-c("1-3","4-6","7-9", "10-12")

#Plot day night and control, general not per hour
#mean_exp_bothphases_randompickproms
#mean_exp_bothphases_randompickernas

#Proms anchors
tempoclock_dayRP<-mean_expbothphases_randompickcircproms1[,1:4]
rownames(tempoclock_dayRP)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24") #c("0", "6", "12", "18")for MFM and c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24") for fang
colnames(tempoclock_dayRP)<-c("1-3","4-6","7-9", "10-12")


tempoclock_nightRP<-mean_expbothphases_randompickcircproms1[,5:8]
rownames(tempoclock_nightRP)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24") #c("0", "6", "12", "18") for MFM and c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24") for fang
colnames(tempoclock_nightRP)<-c("13-15", "16-18", "19-21", "22-24")

t<-cbind(apply(tempoclock_dayRP, 1, mean)/sum(apply(tempoclock_dayRP, 1, mean)), apply(tempoclock_nightRP, 1, mean)/sum(apply(tempoclock_nightRP, 1, mean)), apply(mean_exp_bothphases_randompicknooscproms_controlPP, 1, mean)/sum(apply(mean_exp_bothphases_randompicknooscproms_controlPP, 1, mean)))

colnames(t)<-c("Diurnal", "Nocturnal","Control")


p1<-ggplot(melt(t), aes(x=Var2, y=value, fill=factor(Var1)))+geom_bar(stat = "identity")+
  theme( panel.background = element_blank(), axis.text=element_text(size=12, vjust=0.8), axis.line.x = element_line(color="white", size = 1),axis.line.y = element_line(color="white", size = 1), axis.text.y = element_blank(),axis.title=element_text(size=12), axis.ticks.y = element_blank())+
    scale_fill_brewer(palette="RdBu", name="Phases")+
  #scale_x_discrete(labels= c("1", "4", "7", "10", "13", "16", "19", "22"))+
  ylab("")+xlab("")+ggtitle("Circadian Promoters")
p1<-p1+theme(axis.text.x = element_text(angle = 90, hjust = 1))

    #introns
tempoclock_dayRP<-mean_expbothphases_randompickcircpromsintrons1[,1:2]
rownames(tempoclock_dayRP)<-c("0", "6", "12", "18")
colnames(tempoclock_dayRP)<-c("0", "6")

tempoclock_nightRP<-mean_expbothphases_randompickcircpromsintrons1[,3:4]
rownames(tempoclock_nightRP)<-c("0", "6", "12", "18")
colnames(tempoclock_nightRP)<-c("12", "18")

t<-cbind(apply(tempoclock_dayRP, 1, mean)/sum(apply(tempoclock_dayRP, 1, mean)), apply(tempoclock_nightRP, 1, mean)/sum(apply(tempoclock_nightRP, 1, mean)), apply(mean_exp_bothphases_randompicknooscpromsintrons_controlPP, 1, mean)/sum(apply(mean_exp_bothphases_randompicknooscpromsintrons_controlPP, 1, mean)))

colnames(t)<-c("Diurnal", "Nocturnal","Control")


p2<-ggplot(melt(t), aes(x=Var2, y=value, fill=factor(Var1)))+geom_bar(stat = "identity")+
  theme( panel.background = element_blank(), axis.text=element_text(size=12, vjust=0.8), axis.line.x = element_line(color="white", size = 1),axis.line.y = element_line(color="white", size = 1), axis.text.y = element_blank(),axis.title=element_text(size=12), axis.ticks.y = element_blank())+
    scale_fill_brewer(palette="RdBu", name="Phases")+
  #scale_x_discrete(labels= c("1", "4", "7", "10", "13", "16", "19", "22"))+
  ylab("")+xlab("")+ggtitle("Circadian Promoters\n Introns")
p2<-p2+theme(axis.text.x = element_text(angle = 90, hjust = 1))

#svglite::svglite("Prom-prom_interactions/MFM_RNAseqrandompicked_daynightcontrol_anchorpromTOPROM_introns_proportions.svg")
grid.arrange(p1,p2, ncol=2)
#dev.off()


#svglite::svglite("Prom-prom_interactions/Fang_randompicked_daynightcontrol_anchorpromTOPROM_proportions.svg")
p1
#dev.off()

#####Boxplots
#FanG
zt0circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[1,perm])[ernasphase]})})
zt4circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[2,perm])[ernasphase]})})
zt6circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[3,perm])[ernasphase]})})
zt8circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[4,perm])[ernasphase]})})
zt12circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[5,perm])[ernasphase]})})
zt14circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[6,perm])[ernasphase]})})
zt16circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[7,perm])[ernasphase]})})
zt18circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[8,perm])[ernasphase]})})

diurcircproms<-(zt0circ+zt4circ+zt6circ+zt8circ)/4
noccircproms<-(zt12circ+zt14circ+zt16circ+zt18circ)/4
diurcircproms<-t(apply(diurcircproms, 1, function(x){x/sum(x)}))
noccircproms<-t(apply(noccircproms, 1, function(x){x/sum(x)}))

rownames(diurcircproms)<-NULL
colnames(diurcircproms)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
rownames(diurcircproms)<-rep("Diurnal", 100)
colnames(noccircproms)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
rownames(noccircproms)<-NULL
rownames(noccircproms)<-rep("Nocturnal", 100)



#svglite::svglite("Prom-prom_interactions/Fang_randompicked_daynightcontrol_anchorpromTOPROM_boxplotofinteractions_proportions.svg")
#colnames(matrix.combined)<-rep(c("diurnal", "nocturnal"), 8)
ggplot(melt(rbind(diurcircproms, noccircproms)), aes(x=as.factor(Var2), y=value, fill=as.factor(Var1)))+geom_boxplot(notch = T, notchwidth = 0.5, coef=10)+
  ylab("Proportion of Interactions")+
  scale_x_discrete(name="Circproms phases")+theme_bw()+
  scale_fill_manual(values = c("brown3" , "steelblue", "brown3" , "steelblue", "brown3" , "steelblue", "brown3",   "steelblue", "brown3", "steelblue", "brown3", "steelblue", "brown3", "steelblue", "brown3",  "steelblue"),name="Promoter phases")+ggtitle("Prom-prom: Fang")
#dev.off()



#Stats

#STATISTICAL TESTS
#For fishers test, prom anchors no proportions
##DAY
tempoclock_dayRP<-mean_expbothphases_randompickcircproms1[,1:2]
rownames(tempoclock_dayRP)<-c("0", "6", "12", "18")
colnames(tempoclock_dayRP)<-c("0", "6")
###NIGHT

tempoclock_nightRP<-mean_expbothphases_randompickcircproms1[,3:4]
rownames(tempoclock_nightRP)<-c("0", "6", "12", "18")
colnames(tempoclock_nightRP)<-c("12", "18")

t<-cbind(apply(tempoclock_dayRP, 1, mean), apply(tempoclock_nightRP, 1, mean), apply(mean_exp_bothphases_randompicknooscproms_controlPP, 1, mean))
#fisher.test(matrix(c(sum(t[1:4,1]), sum(t[5:8, 1]), sum(t[1:4, 3]), sum(t[5:8,3])), byrow = T, ncol=2))#p-value = 0.6262 D vs control
#fisher.test(matrix(c(sum(t[1:4,2]), sum(t[5:8, 2]), sum(t[1:4, 3]), sum(t[5:8,3])), byrow = T, ncol=2))#p-value = 0.4686 N vs control
#fisher.test(matrix(c(sum(t[1:4,1]), sum(t[5:8, 1]), sum(t[1:4, 2]), sum(t[5:8,2])),  byrow = T, ncol=2)) #the only one significant, p-value = 0.1297 D vs N

wilcox.test(t[,1], t[,3], paired = T) #Diurnal vs control p-value = 0.125
wilcox.test(t[,2], t[,3], paired = T) #Nocturnal vs control p-value =  0.125
wilcox.test(t[,1], t[,2], paired = T) #Diurnal vs Nocturnal p-value = 0.875

    #introns
##DAY
tempoclock_dayRP<-mean_expbothphases_randompickcircpromsintrons1[,1:2]
rownames(tempoclock_dayRP)<-c("0", "6", "12", "18")
colnames(tempoclock_dayRP)<-c("0", "6")
###NIGHT

tempoclock_nightRP<-mean_expbothphases_randompickcircpromsintrons1[,3:4]
rownames(tempoclock_nightRP)<-c("0", "6", "12", "18")
colnames(tempoclock_nightRP)<-c("12", "18")

t<-cbind(apply(tempoclock_dayRP, 1, mean), apply(tempoclock_nightRP, 1, mean), apply(mean_exp_bothphases_randompicknooscpromsintrons_controlPP, 1, mean))
#fisher.test(matrix(c(sum(t[1:4,1]), sum(t[5:8, 1]), sum(t[1:4, 3]), sum(t[5:8,3])), byrow = T, ncol=2))#p-value = 0.6262 D vs control
#fisher.test(matrix(c(sum(t[1:4,2]), sum(t[5:8, 2]), sum(t[1:4, 3]), sum(t[5:8,3])), byrow = T, ncol=2))#p-value = 0.4686 N vs control
#fisher.test(matrix(c(sum(t[1:4,1]), sum(t[5:8, 1]), sum(t[1:4, 2]), sum(t[5:8,2])),  byrow = T, ncol=2)) #the only one significant, p-value = 0.1297 D vs N

wilcox.test(t[,1], t[,3], paired = T) #Diurnal vs control p-value = 0.625
wilcox.test(t[,2], t[,3], paired = T) #Nocturnal vs control p-value =  0.25
wilcox.test(t[,1], t[,2], paired = T) #Diurnal vs Nocturnal p-value = 0.625



###############results were not st  sig, try only comparing the distrbution
#proms anchor
shapiro.test(t[,1]) #p-value = 0.5905
shapiro.test(t[,2]) #p-value = 0.07562
shapiro.test(t[,3]) #p-value = 0.1351
t.test(t[,1], t[,3]) #p-value =  0.06135
t.test(t[,2], t[,3]) #p-value =0.1947
t.test(t[,1], t[,2])#p-value = 0.9109


    ###----all circproms pvals
#1) Mean of diurnals ZT0 Y 6 for each permutation, I will have 100 values for each of the ZT ernas for the diurnal and noc Circproms
#2) With two distributions for each timepoint i can compare. In the end I will have a vector of 8 pvals 

zt0circ<-sapply(1:4, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[1,perm])[ernasphase]})})

zt6circ<-sapply(1:4, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[2,perm])[ernasphase]})})

zt12circ<-sapply(1:4, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[3,perm])[ernasphase]})})

zt18circ<-sapply(1:4, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[4,perm])[ernasphase]})})


#mean of diurnal and noc
diurcircproms<-(zt0circ+zt6circ)/2
noccircproms<-(zt12circ+zt18circ)/2

#Check if distributions are normal or not
nd<-sapply(c(1:3), function(x){shapiro.test(diurcircproms[,x])[[2]]}) #remove col 5 because its an identical distribution
nn<-sapply(1:3, function(x){shapiro.test(noccircproms[,x])[[2]]})
par(mfrow=c(1,2))
plot(nn, ylim=c(0,1))+abline(h = .001, col="red")
plot(nd, ylim=c(0,1))+abline(h = .001, col="red")
#dev.off()
#normal: 1.2,3 7.8
#not normal: 4,5,6
pvalsallcircproms<-c(t.test(diurcircproms[,1],noccircproms[,1] )[[3]],
                t.test(diurcircproms[,2],noccircproms[,2] )[[3]], 
                t.test(diurcircproms[,3],noccircproms[,3] )[[3]], 
                wilcox.test(diurcircproms[,4],noccircproms[,4] )[[3]])
 
 names(pvalsallcircproms)<-c("0", "6", "12", "18")

 #svglite::svglite("Prom-prom_interactions/MFM_RNAseqrandompicked_daynightcontrol_pvals_diurvsnocALLPROMS.svg")
ggplot(melt(log10(pvalsallcircproms)), aes(x=c("0", "6", "12", "18"), y=value))+geom_point(size=4,aes(col=rownames(melt(log10(pvalsallcircproms)))))+
    geom_point(shape = 1,size=4,colour = "gray49")+ scale_y_continuous(name="-log10 pval")+xlab("Phases")+
   scale_x_discrete( limits=c("0", "6", "12", "18"))+ geom_hline(yintercept = 0.001, color="black", linetype="dashed")+
  ggtitle("Diurnal vs Nocturnal: All circproms")+theme_bw()+
  scale_color_manual(values=brewer.pal(9,"RdBu"), name="Phases")
#dev.off()


    ###----introns circproms pvals
#1) Mean of diurnals ZT0 Y 6 for each permutation, I will have 100 values for each of the ZT ernas for the diurnal and noc Circproms
#2) With two distributions for each timepoint i can compare. In the end I will have a vector of 8 pvals 

zt0circ<-sapply(1:4, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircpromsintrons1[1,perm])[ernasphase]})})

zt6circ<-sapply(1:4, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircpromsintrons1[2,perm])[ernasphase]})})

zt12circ<-sapply(1:4, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircpromsintrons1[3,perm])[ernasphase]})})

zt18circ<-sapply(1:4, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircpromsintrons1[4,perm])[ernasphase]})})


#mean of diurnal and noc
diurcircproms<-(zt0circ+zt6circ)/2
noccircproms<-(zt12circ+zt18circ)/2

#Check if distributions are normal or not
nd<-sapply(c(1,4), function(x){shapiro.test(diurcircproms[,x])[[2]]}) #remove col 5 because its an identical distribution
nn<-sapply(c(1, 3:4), function(x){shapiro.test(noccircproms[,x])[[2]]})
par(mfrow=c(1,2))
plot(nn, ylim=c(0,1))+abline(h = .001, col="red")
plot(nd, ylim=c(0,1))+abline(h = .001, col="red")
#dev.off()
#normal: 1.2,3 7.8
#not normal: 4,5,6
pvalsallcircproms<-c(wilcox.test(diurcircproms[,1],noccircproms[,1] )[[3]],
               wilcox.test(diurcircproms[,2],noccircproms[,2] )[[3]], 
               wilcox.test(diurcircproms[,3],noccircproms[,3] )[[3]], 
                wilcox.test(diurcircproms[,4],noccircproms[,4] )[[3]])
 
 names(pvalsallcircproms)<-c("0", "6", "12", "18")
 
 
#svglite::svglite("Prom-prom_interactions/MFM_RNAseqrandompicked_daynightcontrol_pvals_diurvsnocINTRONS.svg")
ggplot(melt(log10(pvalsallcircproms)), aes(x=c("0", "6", "12", "18"), y=value))+geom_point(size=4,aes(col=rownames(melt(log10(pvalsallcircproms)))))+
    geom_point(shape = 1,size=4,colour = "gray49")+ scale_y_continuous(name="-log10 pval")+xlab("Phases")+
   scale_x_discrete( limits=c("0", "6", "12", "18"))+ geom_hline(yintercept = 0.001, color="black", linetype="dashed")+
  ggtitle("Diurnal vs Nocturnal: All circproms")+theme_bw()+
  scale_color_manual(values=brewer.pal(9,"RdBu"), name="Phases")
#dev.off()


#-------------
#Try proportions
#proms MFM
zt0circ<-sapply(1:4, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[1,perm])[ernasphase]})})
zt6circ<-sapply(1:4, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[2,perm])[ernasphase]})})
zt12circ<-sapply(1:4, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[3,perm])[ernasphase]})})
zt18circ<-sapply(1:4, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[4,perm])[ernasphase]})})
diurcircproms<-(zt0circ+zt6circ)/2
noccircproms<-(zt12circ+zt18circ)/2
diurcircproms<-t(apply(diurcircproms, 1, function(x){x/sum(x)}))
noccircproms<-t(apply(noccircproms, 1, function(x){x/sum(x)}))

rownames(diurcircproms)<-NULL
colnames(diurcircproms)<-c("0","6","12","18")
rownames(diurcircproms)<-rep("Diurnal", 100)
colnames(noccircproms)<-c("0","6","12","18")
rownames(noccircproms)<-NULL
rownames(noccircproms)<-rep("Nocturnal", 100)
#Check if distributions are normal or not
nd<-sapply(c(1:4), function(x){shapiro.test(diurcircproms[,x])[[2]]}) #4th col is ident
nn<-sapply(1:4, function(x){shapiro.test(noccircproms[,x])[[2]]})#4th col is ident
par(mfrow=c(1,2))
plot(nn)+abline(h = .001, col="red")
plot(nd)+abline(h = .001, col="red")
#dev.off()
 
#plot distributions
t<-apply(diurcircproms, 2, density)
par(mfrow=c(1,4))
lapply(t, plot)

t<-apply(noccircproms, 2, density)
par(mfrow=c(1,4))
lapply(t, plot)
#They all look normal

#Since they are normal, check their variances
#Check variances of distributions, 
vnd<-sapply(c(1:4), function(x){var.test(diurcircproms[,x],noccircproms[,x], alternative = "two.sided")$p.value}) 

plot(vnd)+abline(h = .001, col="red")


#vnd>.001 if pval is >.001 I dont have enough evidence to reject h0 which asumes the variances are equivalent, so h0==h1

#[1] TRUE TRUE TRUE   NA
#t.test: 1.2,3
#4th?
pvalsintrons<-c(t.test(diurcircproms[,1],noccircproms[,1], var.equal = T )[[3]],
                t.test(diurcircproms[,2],noccircproms[,2],  var.equal = T  )[[3]], 
                t.test(diurcircproms[,3],noccircproms[,3],  var.equal = T )[[3]] ,
                t.test(diurcircproms[,4],noccircproms[,4], var.equal = T )[[3]])
 
 names(pvalsintrons)<-c("0","6","12", "18")
 
#svglite::svglite("Prom-prom_interactions/MFM_RNAseqrandompicked_daynightcontrol_pvals_diurvsnocALLPROMS_proportions.svg")
ggplot(melt(log10(pvalsintrons)), aes(x=c("0","6","12", "18"), y=value))+geom_point(size=4, aes(col=rownames(melt(log10(pvalsintrons)))))+  scale_y_continuous(name="-log10 pval")+xlab("Phases")+
    geom_point(shape = 1,size=4,colour = "gray49")+
   scale_x_discrete( limits=c("0","6","12", "18"))+ geom_hline(yintercept = -3, color="black", linetype="dashed")+
  ggtitle("Diurnal vs Nocturnal: All proms ")+theme_bw()+scale_color_manual(values=brewer.pal(4,"RdBu"), name="Phases")
#dev.off()
#2020-03-16 Proms FANG----------
zt0circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[1,perm])[ernasphase]})})
zt4circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[2,perm])[ernasphase]})})
zt6circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[3,perm])[ernasphase]})})
zt8circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[4,perm])[ernasphase]})})
zt12circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[5,perm])[ernasphase]})})
zt14circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[6,perm])[ernasphase]})})
zt16circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[7,perm])[ernasphase]})})
zt18circ<-sapply(1:8, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircproms1[8,perm])[ernasphase]})})

diurcircproms<-(zt0circ+zt4circ+zt6circ+zt8circ)/4
noccircproms<-(zt12circ+zt14circ+zt16circ+zt18circ)/4
diurcircproms<-t(apply(diurcircproms, 1, function(x){x/sum(x)}))
noccircproms<-t(apply(noccircproms, 1, function(x){x/sum(x)}))

rownames(diurcircproms)<-NULL
colnames(diurcircproms)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
rownames(diurcircproms)<-rep("Diurnal", 100)
colnames(noccircproms)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
rownames(noccircproms)<-NULL
rownames(noccircproms)<-rep("Nocturnal", 100)
#Check if distributions are normal or not
nd<-sapply(c(1:8), function(x){shapiro.test(diurcircproms[,x])[[2]]}) #4th col is ident
nn<-sapply(1:8, function(x){shapiro.test(noccircproms[,x])[[2]]})#4th col is ident
par(mfrow=c(1,2))
plot(nn)+abline(h = .001, col="red")
plot(nd)+abline(h = .001, col="red")
#dev.off() #All look normal
 
#plot distributions
t<-apply(diurcircproms, 2, density)
par(mfrow=c(1,4))
lapply(t, plot)

t<-apply(noccircproms, 2, density)
par(mfrow=c(1,4))
lapply(t, plot)
#They all look normal

#Since they are normal, check their variances
#Check variances of distributions, 
vnd<-sapply(c(1:8), function(x){var.test(diurcircproms[,x],noccircproms[,x], alternative = "two.sided")$p.value}) 

plot(vnd)+abline(h = .001, col="red")


#vnd>.001 if pval is >.001 I dont have enough evidence to reject h0 which asumes the variances are equivalent, so h0==h1

#Var equal: TRUE FALSE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE
#t.test: 1.2,3 V
#4th?
pvalsallcircproms_fang<-c(t.test(diurcircproms[,1],noccircproms[,1], var.equal = T )[[3]],
                t.test(diurcircproms[,2],noccircproms[,2],  var.equal = F  )[[3]], 
                t.test(diurcircproms[,3],noccircproms[,3],  var.equal = T )[[3]] ,
                t.test(diurcircproms[,4],noccircproms[,4], var.equal = F )[[3]],
                t.test(diurcircproms[,5],noccircproms[,5], var.equal = T )[[3]],
                t.test(diurcircproms[,6],noccircproms[,6], var.equal = T )[[3]],
                t.test(diurcircproms[,7],noccircproms[,7], var.equal = T )[[3]],
                t.test(diurcircproms[,8],noccircproms[,8], var.equal = T )[[3]])
 
 names(pvalsallcircproms_fang)<-c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24")
 
#svglite::svglite("Prom-prom_interactions/FANGrandompicked_daynightcontrol_pvals_diurvsnoc_proportions.svg")
ggplot(melt(log10(pvalsallcircproms_fang)), aes(x=c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24"), y=value))+geom_point(size=4, aes(col=rownames(melt(log10(pvalsallcircproms_fang)))))+  scale_y_continuous(name="-log10 pval")+xlab("Phases")+
    geom_point(shape = 1,size=4,colour = "gray49")+
   scale_x_discrete( limits=c("1-3","4-6","7-9", "10-12","13-15", "16-18", "19-21", "22-24"))+ geom_hline(yintercept = -3, color="black", linetype="dashed")+
  ggtitle("Diurnal vs Nocturnal: All proms FANG ")+theme_bw()+scale_color_manual(values=brewer.pal(8,"RdBu"), name="Phases")
#dev.off()

#--------------------

#introns
zt0circ<-sapply(1:4, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircpromsintrons1[1,perm])[ernasphase]})})
zt6circ<-sapply(1:4, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircpromsintrons1[2,perm])[ernasphase]})})
zt12circ<-sapply(1:4, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircpromsintrons1[3,perm])[ernasphase]})})
zt18circ<-sapply(1:4, function(ernasphase){sapply(1:100, function(perm){unlist(expbothphases_randompickcircpromsintrons1[4,perm])[ernasphase]})})
diurcircproms<-(zt0circ+zt6circ)/2
noccircproms<-(zt12circ+zt18circ)/2
diurcircproms<-t(apply(diurcircproms, 1, function(x){x/sum(x)}))
noccircproms<-t(apply(noccircproms, 1, function(x){x/sum(x)}))

rownames(diurcircproms)<-NULL
colnames(diurcircproms)<-c("0","6","12","18")
rownames(diurcircproms)<-rep("Diurnal", 100)
colnames(noccircproms)<-c("0","6","12","18")
rownames(noccircproms)<-NULL
rownames(noccircproms)<-rep("Nocturnal", 100)
#Check if distributions are normal or not
nd<-sapply(c(1:2, 4), function(x){shapiro.test(diurcircproms[,x])[[2]]}) #4th col is ident
nn<-sapply(c(1,3:4), function(x){shapiro.test(noccircproms[,x])[[2]]})#4th col is ident
par(mfrow=c(1,2))
plot(nn)+abline(h = .001, col="red")
plot(nd)+abline(h = .001, col="red")
#dev.off()

nn>.001
#[1] TRUE TRUE NA TRUE
nd>.001
#[1] FALSE NA FALSE  TRUE
 
#plot distributions
t<-apply(diurcircproms, 2, density)
par(mfrow=c(1,4))
lapply(t, plot)

t<-apply(noccircproms, 2, density)
par(mfrow=c(1,4))
lapply(t, plot)


#Since 3rd, check their variances
#Check variances of distributions, 
vnd<-sapply(c(4), function(x){var.test(diurcircproms[,x],noccircproms[,x], alternative = "two.sided")$p.value}) 

plot(vnd)+abline(h = .001, col="red")


#vnd>.001 if pval is >.001 I dont have enough evidence to reject h0 which asumes the variances are equivalent, so h0==h1

#[1] TRUE TRUE TRUE   NA
#t.test: 1.2,3
#4th?
pvalsintrons<-c(wilcox.test(diurcircproms[,1],noccircproms[,1],paired = F )[[3]],
                wilcox.test(diurcircproms[,2],noccircproms[,2],  paired = F  )[[3]], 
                wilcox.test(diurcircproms[,3],noccircproms[,3], paired = F )[[3]], 
                t.test(diurcircproms[,4],noccircproms[,4], var.equal = F )[[3]])
 
 names(pvalsintrons)<-c("0","6","12", "18")
 
#svglite::svglite("Prom-prom_interactions/MFM_RNAseqrandompicked_daynightcontrol_pvals_diurvsnocINTRONS_proportions.svg")
ggplot(melt(log10(pvalsintrons)), aes(x=c("0","6","12", "18"), y=value))+geom_point(size=4, aes(col=rownames(melt(log10(pvalsintrons)))))+  scale_y_continuous(name="-log10 pval")+xlab("Phases")+
    geom_point(shape = 1,size=4,colour = "gray49")+
   scale_x_discrete( limits=c("0","6","12", "18"))+ geom_hline(yintercept = -3, color="black", linetype="dashed")+
  ggtitle("Diurnal vs Nocturnal: Introns ")+theme_bw()+scale_color_manual(values=brewer.pal(4,"RdBu"), name="Phases")
#dev.off()
