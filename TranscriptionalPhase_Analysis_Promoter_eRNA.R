#####################---osceRNAs-cirproms interactions per phase---##############
chic<-read.csv("/inputs_for_scripts/All_step2_washU_text_symm.txt", sep = "\t", header = F)
chic_baits<-chic[,1:3]
chic_otherends<-chic[,4:6]

#Bedfiles, GRanges Objects
chic_otherends_bed<-GRanges(seqnames= Rle(chic_otherends[,1]), ranges = IRanges(chic_otherends[,2], chic_otherends[,3]))
chic_bait_bed<-GRanges(seqnames= Rle(chic_baits[,1]), ranges = IRanges(chic_baits[,2], chic_baits[,3]))



#--------------------------------------------------------------------
##Overlap of eRNAs per phases

#open ernas data. Fang 2014
osc_ernas<-read.csv("/inputs_for_scripts/eRNAs_de_novo_oscillating_phases.txt", header = T, sep = "\t")

#Group per phases
#List of the phases
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

names(list_ernasperphase_bed) <- paste(phases_bin_names, "bed_ernas",  sep = "")
list2env(list_ernasperphase_bed , envir = .GlobalEnv)

####################Overlap with ChiC
ernaslist_phases<-sapply(mixedsort(ls(pattern = "bed_ernas", sorted = F), decreasing = F), function(x){list(x)})

#####################Retrieve index from bait, ernas per phases
phases_eRNAs_index_baits<-lapply(ernaslist_phases, function(x){
  tempo<-retrieveindex_fromchic(get(x), chic_otherends_bed, chic_bait_bed)
})

names(phases_eRNAs_index_baits) <- paste("BI_ernas", phases_bin_names, sep = "")
list2env(phases_eRNAs_index_baits , envir = .GlobalEnv)
#phases_bin_names for fang and 
#remove(phases_eRNAs_index_baits)


#Open circadian promoters
fang_mrnas<-read.csv("/inputs_for_scripts/circproms/Updated/HindIIIfragments_circadiangenes_fang_phases.bed", sep = "\t", header = F)
fang_mrnas<-bedfile(fang_mrnas, columnnames)
circpromlist<-list("fang_mrnas")

circadianprom_phases_ernasbaits_counts<-
  lapply(circpromlist, function(y){
    tempo_circ<-get(y);
    t<-lapply(phases_eRNAs_index_baits, function(tempo_ernas){
    tempo_ernas<-length(findOverlaps(tempo_ernas,tempo_circ))
    return(tempo_ernas)
    })
  })



names(circadianprom_phases_ernasbaits_counts)<-circpromlist

##---------------------------------------------------------------------------
##Overlap of circprom per phases
circpromphases<-read.csv("/inputs_for_scripts/circproms/MFM_RNAseq/HindIIIfragments_circadiangenes_MFMRNAseq_phases.bed", header = F, sep="\t")
#Run only for control
circpromphases<-read.csv("/inputs_for_scripts/circproms/Updated/HindIIIfragments_circadiangenes_fang_phases.bed", header = F, sep="\t")
circpromphases<-circpromphases[,-(4)]
circpromphasesINTRONS<-read.csv("/inputs_for_scripts/circproms/MFM_RNAseq/HindIIIfragments_circadiangenesonlyintrons_MFMRNAseq_phases.bed", header = F, sep="\t")

#Divide per phases
#List of the phases
#phases<-unique(circpromphases[,17])
#SORT FILE ACCORDING TO PHASES
circpromphases<-circpromphases[order(circpromphases[,5]),]
phases<-unique(circpromphases[,5])
#introns
circpromphasesINTRONS<-circpromphasesINTRONS[order(circpromphasesINTRONS[,5]),]
#For FANG
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
  tempo4<-circpromphases[circpromphases[,5]==x[4],]
  tempo5<-rbind(tempo1, tempo2, tempo3, tempo4)
  return(tempo5)
})

#Introns
circpromINTRONS_divided_list<-lapply(phases, function(x){
  tempo1<-circpromphasesINTRONS[circpromphasesINTRONS[,5]==x[1],]
  tempo2<-circpromphasesINTRONS[circpromphasesINTRONS[,5]==x[2],]
  tempo3<-circpromphasesINTRONS[circpromphasesINTRONS[,5]==x[3],]
  tempo4<-circpromphasesINTRONS[circpromphasesINTRONS[,5]==x[4],]
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




circpromlist_phases<-sapply(mixedsort(ls(pattern = "bed_circprom", sorted = F), decreasing = F), function(x){list(x)})
  #introns
circpromlistintrons_phases<-sapply(mixedsort(ls(pattern = "bed_intronscircprom", sorted = F), decreasing = F), function(x){list(x)})

#Bed eRNAs 
osc_ernas<-osc_ernas[,1:3]
osc_ernas<-bedfile(osc_ernas, columnnames)


#Overlap all osc ernas to chic and retrieve their interacting bait
ernaslist_phases<-list("osc_ernas")
circprom_phases_eRNAs_index_baits<-lapply(ernaslist_phases, function(x){
  tempo<-get(x);
  tempo<-retrieveindex_fromchic(tempo, chic_otherends_bed, chic_bait_bed)
}) #only a granges object with the baits is recovered, 9730 baits


#Count the interactions of the baits that interact with circproms per phase
circadianprom_phases_ernasbaits_counts<-lapply(circpromlist_phases, function(y){
    tempo_circ<-get(y);
    t<-lapply(circprom_phases_eRNAs_index_baits, function(tempo_ernas){
      tempo_ernas<-length(findOverlaps(tempo_ernas,tempo_circ))
      })
  })

names(circadianprom_phases_ernasbaits_counts)<-circpromlist_phases

    ###Introns------------------------
#Count the interactions of the baits that interact with circproms per phase
circadianpromintrons_phases_ernasbaits_counts<-lapply(circpromlistintrons_phases, function(y){
    tempo_circ<-get(y);
    t<-lapply(circprom_phases_eRNAs_index_baits, function(tempo_ernas){
      tempo_ernas<-length(findOverlaps(tempo_ernas,tempo_circ))
      })
  })

names(circadianpromintrons_phases_ernasbaits_counts)<-circpromlistintrons_phases


#Plot results
#num_ernas_percircpromphases<-as.data.frame(cbind(c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24"), unlist(circadianprom_phases_ernasbaits_counts)))
#num_ernas_percircpromphases[,2]<-as.numeric(levels(num_ernas_percircpromphases[,2]))[num_ernas_percircpromphases[,2]]

#pp4<-ggplot(num_ernas_percircpromphases, aes(x=V1, y=V2, group=1 ))+theme_bw()+geom_point()+geom_line()+xlab("Promoter Phases")+ylab("Associated oscillatory eRNAs")+scale_x_discrete(limits=c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24"))+theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 20), axis.text.y = element_text(size = 20), axis.title.x = element_text(size=20), axis.title.y=element_text(size=20))

#tiff("promsperphases.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
#pp4
##dev.off()

#--------------------------------------------------------------------------

###### Both ernas and circprom per phase analysis
#Compare each circadian promoter vs all the ernas_phases,

#Compare each circadian promoter vs all the ernas_phases,
#Rerun chunk in line 104
bothphases_circadianprom_ernas_counts<-
  lapply(circpromlist_phases, function(y){# each circprom
    tempo_circ<-get(y);
    t<-lapply(phases_eRNAs_index_baits, function(tempo_ernas){
    tempo_ernas<-length(subjectHits(findOverlaps(tempo_ernas, tempo_circ)))
    })
    
})


names(bothphases_circadianprom_ernas_counts)<-circpromlist_phases

# intrnos
bothphases_circadianpromintrons_ernas_counts<-
  lapply(circpromlistintrons_phases, function(y){# each circprom
    tempo_circ<-get(y);
    t<-lapply(phases_eRNAs_index_baits, function(tempo_ernas){
    tempo_ernas<-length(subjectHits(findOverlaps(tempo_ernas, tempo_circ)))
    })
    
})


names(bothphases_circadianpromintrons_ernas_counts)<-circpromlistintrons_phases



#Create files of pairs circproms-osceRNAs diurnal and nocturnal

#Compare each circadian promoter vs all the ernas_phases,
bothphases_circadianprom_ernas_counts<-
  lapply(circpromlist_phases, function(y){# each circprom
    tempo_circ<-get(y);
    t<-lapply(phases_eRNAs_index_baits, function(tempo_ernas){
    tempo_ernas<-length(subjectHits(findOverlaps(tempo_ernas, tempo_circ)))
    })
    
})


names(bothphases_circadianprom_ernas_counts)<-circpromlist_phases

# introns
bothphases_circadianpromintrons_ernas_counts<-
  lapply(circpromlistintrons_phases, function(y){# each circprom
    tempo_circ<-get(y);
    t<-lapply(phases_eRNAs_index_baits, function(tempo_ernas){
    tempo_ernas<-length(subjectHits(findOverlaps(tempo_ernas, tempo_circ)))
    })
    
})


names(bothphases_circadianpromintrons_ernas_counts)<-circpromlistintrons_phases

#######################################
#Random pick, but anchor only one either the proms or the eRNas
#------------
#RANDOM ERNAS, now try random proms
#Extra: EXPECTED
#Repeat 100 times
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
  tempo<-retrieveindex_fromchic(x, chic_otherends_bed, chic_bait_bed)
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
  tempo<-retrieveindex_fromchic(x, chic_otherends_bed, chic_bait_bed)
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
  #Introns
expbothphases_randompickpromsINTRONS<-sapply(1:100, function(rep){
  #RANlist_ernasperphase_bed<-sapply(1:8, function(phase){
  #  pos<-sample(list_ernasperphase_bed[[phase]] ,160, replace = F)
  #  })
                #lapply(list_ernasperphase_bed, length)[[phase]], replace = F)
#Random no osc proms
  RANlist_circpromperphase_bed<-sapply(1:4, function(phase){
  pos<-sample(list_circpromsperphaseINTRONS_bed[[phase]],32, replace = F)
  })
#Retrieve indexes of eRNAs
  indexes_baitsofernas_perphases<-lapply(list_ernasperphase_bed, function(x){
  tempo<-retrieveindex_fromchic(x, chic_otherends_bed, chic_bait_bed)
  })
  
  #Bait(index ernas) overlap with proms

  RANbothphases_circadianprom_ernas_counts<-
  lapply(RANlist_circpromperphase_bed, function(tempo_circ){# each circprom
    t<-lapply(indexes_baitsofernas_perphases  , function(tempo_ernas){
    tempo_ernas<-length(subjectHits(findOverlaps(tempo_ernas, tempo_circ)))
    })
})
  names(RANbothphases_circadianprom_ernas_counts)<-circpromlistintrons_phases
  
  
  #names(RANoverlaps_and_correlations_bothphases_circadianprom_ernas)<-circpromlist_phases
  return(RANbothphases_circadianprom_ernas_counts)

})

#Obtain the mean, 
mean_exp_bothphases_randompickproms<-sapply(1:4, function(cirpromphase){sapply(1:8, function(ernasphase){mean(unlist(sapply(1:100, function(rep){expbothphases_randompickproms[cirpromphase,][[rep]][ernasphase]})))})})

mean_exp_bothphases_randompickernas<-sapply(1:8, function(ernasphase){sapply(1:4, function( cirpromphase){mean(unlist(sapply(1:100, function(rep){expbothphases_randompickernas[cirpromphase,][[rep]][ernasphase]})))})})

  #Introns
mean_exp_bothphases_randompickpromsINTRONS<-sapply(1:4, function(cirpromphase){sapply(1:8, function(ernasphase){mean(unlist(sapply(1:100, function(rep){expbothphases_randompickpromsINTRONS[cirpromphase,][[rep]][ernasphase]})))})})



mean_exp_bothphases_randompickproms
mean_exp_bothphases_randompickernas
mean_exp_bothphases_randompickpromsINTRONS

#------------
#Plot
###DAY
tempoclock<-sapply(1:4, function(x){
  unlist(bothphases_circadianprom_ernas_counts[[x]])
})
tempoclock_dayRP<-mean_exp_bothphases_randompickernas[,1:4]
rownames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12")
tempoclock_nightRP<-mean_exp_bothphases_randompickernas[,5:8]
rownames(tempoclock_nightRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_nightRP)<-c("13-15", "16-18", "19-21", "22-24")
tempoclock_dayRP<-mean_exp_bothphases_randompickproms[,1:4]
rownames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12")
tempoclock_nightRP<-mean_exp_bothphases_randompickproms[,5:8]
rownames(tempoclock_nightRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_nightRP)<-c("13-15", "16-18", "19-21", "22-24")
#SKIP AND GO TO 984
#Plot Obs and Exp
DOE<-ggplot(melt(tempoclock_dayRP), aes(x=Var1, y=value, group=Var2, col=Var2))+theme_bw()+
  geom_point()+geom_line()+xlab("Promoter Phases")+ylab("Number of Interactions")+theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12), axis.text.y = element_text(size = 12),axis.title.x = element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_blank(), legend.text =  element_text(size=12))+scale_colour_brewer(palette = "Spectral")+
  ggtitle("Diurnal eRNAs")
  
  #geom_hline(yintercept = 1, linetype =2)

#  annotate("text", x = 1:8, y = 505, label = paste("n=", unlist(lapply(list_ernasperphase_bed, length)), sep = ""))
 #tiff("randompicked_diurnalernas.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
#DOE
##dev.off()


###NIGHT
tempoclock<-sapply(5:8, function(x){
  unlist(bothphases_circadianprom_ernas_counts[[x]])
})
tempoclock_nightRP<-mean_exp_bothphases_randompickernas[,5:8]
rownames(tempoclock_nightRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_nightRP)<-c("13-15", "16-18", "19-21", "22-24")

#Plot Obs and Exp
NOE<-ggplot(melt(tempoclock_nightRP), aes(x=Var1, y=value, group=Var2, col=Var2))+theme_bw()+
  geom_point()+geom_line()+xlab("Promoter Phases")+ylab("Number of interactions")+theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12), axis.text.y = element_text(size = 12),axis.title.x = element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_blank(), legend.text =  element_text(size=12))+scale_colour_brewer(palette = "Dark2")+ggtitle("Nocturnal Promoters")
#  annotate("text", x = 1:8, y = 505, label = paste("n=", unlist(lapply(list_ernasperphase_bed, length)), sep = ""))
tiff("randompicked_nocturnalernas.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
NOE
#dev.off()

####Plot as proportions

#ernas anchor
tempoclock_dayRP<-mean_exp_bothphases_randompickernas[,1:4]
rownames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12")
tempoclock_nightRP<-mean_exp_bothphases_randompickernas[,5:8]
rownames(tempoclock_nightRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_nightRP)<-c("13-15", "16-18", "19-21", "22-24")


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
tempoclock_dayRP<-mean_exp_bothphases_randompickproms[,1:4]
rownames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12")
tempoclock_nightRP<-mean_exp_bothphases_randompickproms[,5:8]
rownames(tempoclock_nightRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_nightRP)<-c("13-15", "16-18", "19-21", "22-24")


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

###Control randompicked, use no osc eRNAs
poolofenas<-bedfile(read.csv("/inputs_for_scripts/eRNAs_de_novo_oscillating.txt", header = T, sep = "\t"), columnnames)
noosc_ernas<-poolofenas[!(1:19086 %in% queryHits(findOverlaps(poolofenas, (osc_ernas))))]

#Retrieve non circadian promoters, 
#First retrieve only baits
onlybaits<-bedfile(read.csv("/inputs_for_scripts/Baits.txt", header = T, sep="\t"), columnnames)
t<-chic_bait_bed[unique(queryHits(findOverlaps(chic_bait_bed, onlybaits)))]
noncircproms<-t[!(1:247943 %in% queryHits(findOverlaps(t,bedfile(circpromphases))))]
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
  tempo<-retrieveindex_fromchic(x, chic_otherends_bed, chic_bait_bed)
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
  tempo<-retrieveindex_fromchic(x, chic_otherends_bed, chic_bait_bed)
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

####Introns
expbothphases_randompicknooscpromsINTRONS_control<-sapply(1:100, function(rep){
  #RANlist_ernasperphase_bed<-sapply(1:8, function(phase){
   # pos<-sample(noosc_ernas ,160, replace = F)
    #})
                #lapply(list_ernasperphase_bed, length)[[phase]], replace = F)
#Random no osc proms
  RANlist_circpromINTRONSperphase_bed<-sapply(1:8, function(phase){
  pos<-sample(unique(noncircproms), 193, replace=F)})
    #sample(list_circpromsperphase_bed[[phase]],109, replace = F)})
#Retrieve indexes of eRNAs
  RANindexes_baitsofernas_perphases<-lapply(list_ernasperphase_bed, function(x){
  tempo<-retrieveindex_fromchic(x, chic_otherends_bed, chic_bait_bed)
  })
  
  #Bait(index ernas) overlap with proms

  RANbothphases_circadianpromINTRONS_ernas_counts<-
  lapply(RANlist_circpromINTRONSperphase_bed, function(tempo_circ){# each circprom
    t<-lapply(RANindexes_baitsofernas_perphases, function(tempo_ernas){
    tempo_ernas<-length(subjectHits(findOverlaps(tempo_ernas, tempo_circ)))
    })
})
  names(RANbothphases_circadianpromINTRONS_ernas_counts)<-circpromlistintrons_phases
  
  
  #names(RANoverlaps_and_correlations_bothphases_circadianprom_ernas)<-circpromlist_phases
  return(RANbothphases_circadianpromINTRONS_ernas_counts)

})


#Obtain the mean, 
mean_exp_bothphases_randompicknooscproms_control<-sapply(1:4, function(cirpromphase){sapply(1:8, function(ernasphase){mean(unlist(sapply(1:100, function(rep){expbothphases_randompicknooscproms_control[cirpromphase,][[rep]][ernasphase]})))})})

mean_exp_bothphases_randompicknooscernas_control<-sapply(1:8, function(ernasphase){sapply(1:4, function(cirpromphase){mean(unlist(sapply(1:100, function(rep){expbothphases_randompicknooscernas_control[cirpromphase,][[rep]][ernasphase]})))})})

    #Introns
mean_expbothphases_randompicknooscpromsINTRONS_control<-sapply(1:4, function(cirpromphase){sapply(1:8, function(ernasphase){mean(unlist(sapply(1:100, function(rep){expbothphases_randompicknooscpromsINTRONS_control[cirpromphase,][[rep]][ernasphase]})))})})

mean_exp_bothphases_randompicknooscproms_control
mean_exp_bothphases_randompicknooscernas_control
mean_expbothphases_randompicknooscpromsINTRONS_control

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
colnames(tempoclock_dayRP)<-c("0", "6")#c("0", "6") for MFM and c("1-3", "4-6", "7-9", "10-12") for fang
tempoclock_nightRP<-mean_exp_bothphases_randompickproms[,3:4] #3:4 MFM and 5:8 for fang
rownames(tempoclock_nightRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_nightRP)<-c("12", "18") # c("12", "18") for MFM and c("13-15", "16-18", "19-21", "22-24")

tempoclock_nightRP_control<-mean_exp_bothphases_randompicknooscproms_control[,3:4]
rownames(tempoclock_nightRP_control)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_nightRP_control)<-c("12", "18") #c("12", "18") for MFM and c("13-15", "16-18", "19-21", "22-24") for fang

tempoclock_dayRP_control<-mean_exp_bothphases_randompicknooscproms_control[,1:2]
rownames(tempoclock_dayRP_control)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_dayRP_control)<-c("0", "6") #c("0", "6") for MFM and c("1-3", "4-6", "7-9", "10-12") for fang


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
rownames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")#c("0", "6", "12", "18") for MFM and c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24") for fang
colnames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12")
tempoclock_nightRP<-mean_exp_bothphases_randompickernas[,5:8]
rownames(tempoclock_nightRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")#c("0", "6", "12", "18") for MFM and c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24") for fang
colnames(tempoclock_nightRP)<-c("13-15", "16-18", "19-21", "22-24")

tempoclock_nightRP_control<-mean_exp_bothphases_randompicknooscernas_control[,5:8]
rownames(tempoclock_nightRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")#c("ZT0", "ZT6", "ZT12", "ZT18") for MFM and c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24") for fang
colnames(tempoclock_nightRP)<-c("13-15", "16-18", "19-21", "22-24")

tempoclock_dayRP_control<-mean_exp_bothphases_randompicknooscernas_control[,1:4]
rownames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")#c("ZT0", "ZT6", "ZT12", "ZT18") for MFM and c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24") for fang
colnames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12")


t<-cbind(apply(tempoclock_dayRP, 1, mean)/sum(apply(tempoclock_dayRP, 1, mean)), apply(tempoclock_nightRP, 1, mean)/sum(apply(tempoclock_nightRP, 1, mean)), apply(mean_exp_bothphases_randompicknooscernas_control, 1, mean)/sum(apply(mean_exp_bothphases_randompicknooscernas_control, 1, mean)))
colnames(t)<-c("Diurnal", "Nocturnal","Control")

p2<-ggplot(melt(t), aes(x=Var2, y=value, fill=Var1))+geom_bar(stat = "identity")+
  theme( panel.background = element_blank(), axis.text=element_text(size=12, vjust=0.8), axis.line.x = element_line(color="white", size = 1),axis.line.y = element_line(color="white", size = 1), axis.text.y = element_blank(),axis.title=element_text(size=12), axis.ticks.y = element_blank())+
    scale_fill_brewer(palette="RdBu", name="Phases")+
  #scale_x_discrete(labels= c("1", "4", "7", "10", "13", "16", "19", "22"))+
  ylab("")+xlab("")+ggtitle("eRNAs")
p2<-p2+theme(axis.text.x = element_text(angle = 90, hjust = 1))


######introns
#Proms-introns anchors

tempoclock_dayRP<-mean_exp_bothphases_randompickpromsINTRONS[,1:2]
rownames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_dayRP)<-c("0", "6")
tempoclock_nightRP<-mean_exp_bothphases_randompickpromsINTRONS[,3:4]
rownames(tempoclock_nightRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_nightRP)<-c("12", "18")

tempoclock_nightRP_control<-mean_expbothphases_randompicknooscpromsINTRONS_control[,3:4]
rownames(tempoclock_nightRP_control)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_nightRP_control)<-c("12", "18")

tempoclock_dayRP_control<-mean_expbothphases_randompicknooscpromsINTRONS_control[,1:2]
rownames(tempoclock_dayRP_control)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_dayRP_control)<-c("0", "6")


t<-cbind(apply(tempoclock_dayRP, 1, mean)/sum(apply(tempoclock_dayRP, 1, mean)), apply(tempoclock_nightRP, 1, mean)/sum(apply(tempoclock_nightRP, 1, mean)), apply(mean_expbothphases_randompicknooscpromsINTRONS_control, 1, mean)/sum(apply(mean_expbothphases_randompicknooscpromsINTRONS_control, 1, mean)))

colnames(t)<-c("Diurnal", "Nocturnal","Control")


p3<-ggplot(melt(t), aes(x=Var2, y=value, fill=Var1))+geom_bar(stat = "identity")+
  theme( panel.background = element_blank(), axis.text=element_text(size=12, vjust=0.8), axis.line.x = element_line(color="white", size = 1),axis.line.y = element_line(color="white", size = 1), axis.text.y = element_blank(),axis.title=element_text(size=12), axis.ticks.y = element_blank())+
    scale_fill_brewer(palette="RdBu", name="Phases")+
  #scale_x_discrete(labels= c("1", "4", "7", "10", "13", "16", "19", "22"))+
  ylab("")+xlab("")+ggtitle("Circadian \n Promoters Introns")
p3<-p3+theme(axis.text.x = element_text(angle = 90, hjust = 1))

##svglite::svglite("/Analysis_per_phases/MFM_RNAseqrandompicked_daynightcontrol_anchorprom_anchorpromintron_anchorernas_proportions.svg")
grid.arrange(p1, p2, p3, ncol=3)
##dev.off()

#fang control only p1 and p2
##svglite::svglite("/Analysis_per_phases/Fang_GROseqrandompicked_daynightcontrol_anchorprom_anchorernas_proportions.svg")
grid.arrange(p1, p2,ncol=2)
##dev.off()

#Stats

#Statistical test to see if there is ss difference between D-C N-C and D-N

#mean_exp_bothphases_randompickproms
#mean_exp_bothphases_randompickernas
#mean_exp_bothphases_randompickpromsINTRONS
tempoclock_dayRP<-mean_exp_bothphases_randompickproms[,1:2]
rownames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_dayRP)<-c("1-3", "4-6", "7-9", "10-12")
tempoclock_nightRP<-mean_exp_bothphases_randompickproms[,3:4]
rownames(tempoclock_nightRP)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
colnames(tempoclock_nightRP)<-c("13-15", "16-18", "19-21", "22-24")

#For fishers test, prom anchors no proportions
t<-cbind(apply(tempoclock_dayRP, 1, mean), apply(tempoclock_nightRP, 1, mean), apply(mean_exp_bothphases_randompicknooscproms_control, 1, mean))
colnames(t)<-c("Diurnal", "Nocturnal", "Control")
#fisher.test(matrix(c(sum(t[1:4,1]), sum(t[5:8, 1]), sum(t[1:4, 3]), sum(t[5:8,3])), byrow = T, ncol=2))#p-value = 0.5471
#fisher.test(matrix(c(sum(t[1:4,2]), sum(t[5:8, 2]), sum(t[1:4, 3]), sum(t[5:8,3])), byrow = T, ncol=2))#p-value =  0.1304
#fisher.test(matrix(c(sum(t[1:4,1]), sum(t[5:8, 1]), sum(t[1:4, 2]), sum(t[5:8,2])),  byrow = T, ncol=2)) #the only one significant, p-value = 0.005327


wilcox.test(t[,1], t[,3], paired = T) #Diurnal vs control p-value = 0.007813
wilcox.test(t[,2], t[,3], paired = T) #Nocturnal vs control p-value =  0.007813
wilcox.test(t[,1], t[,2], paired = T) #Diurnal vs Nocturnal p-value = 0.1484
#Fang: Diurnal vs control p-value = 0.007813, Nocturnal vs control p-value =  0.007813, Diurnal vs Nocturnal p-value = 0.1484 

##svglite::svglite("/Analysis_per_phases/MFM_RNAseqrandompicked_daynightcontrol_anchorprom_anchorproms_rawinteractionscounts.svg")
ggplot(melt(t), aes(x=Var1, y=value, group=Var2, col=Var2))+geom_line()+xlab("eRNAs phase")+ylab("Number of interaction with eRNAs")+ggtitle("Circproms")+scale_color_brewer(palette="Set1")+scale_color_brewer(palette="Set1")+guides(color=guide_legend(title="Circproms phase"))+theme_classic()
##dev.off()

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

##svglite::svglite("Analysis_per_phases/MFM_RNAseqrandompicked_daynightcontrol_anchorprom_anchorernas_rawinteractionscounts.svg")
ggplot(melt(t), aes(x=Var1, y=value, group=Var2, col=Var2))+geom_line()+xlab("Circproms phase")+ylab("Number of interaction with circproms")+ggtitle("eRNAs")+scale_color_brewer(palette="Set1")+guides(color=guide_legend(title="eRNAs phase"))+theme_classic()
##dev.off()

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
pvalsallcircproms<-c(t.test(diurcircproms[,1],noccircproms[,1] )[[3]],
                t.test(diurcircproms[,2],noccircproms[,2] )[[3]], 
                t.test(diurcircproms[,3],noccircproms[,3] )[[3]], 
                wilcox.test(diurcircproms[,4],noccircproms[,4] )[[3]], 
                wilcox.test(diurcircproms[,5],noccircproms[,5] )[[3]], 
                wilcox.test(diurcircproms[,6],noccircproms[,6] )[[3]], 
                t.test(diurcircproms[,7],noccircproms[,7] )[[3]], 
                t.test(diurcircproms[,7],noccircproms[,7] )[[3]])
 
 names(pvalsallcircproms)<-c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24")
 
##svglite::svglite("Analysis_per_phases/MFM_RNAseqrandompicked_daynightcontrol_pvals_diurvsnocALLPROMS.svg")
ggplot(melt(log10(pvalsallcircproms)), aes(x=c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24"), y=value))+geom_point(size=4,aes(col=rownames(melt(log10(pvalsallcircproms)))))+
    geom_point(shape = 1,size=4,colour = "gray49")+ scale_y_continuous(name="-log10 pval")+xlab("Phases")+
   scale_x_discrete( limits=c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24"))+ geom_hline(yintercept = 0.001, color="black", linetype="dashed")+
  ggtitle("Diurnal vs Nocturnal: All circproms")+theme_bw()+
  scale_color_manual(values=brewer.pal(9,"RdBu"), name="Phases")
##dev.off()



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

##svglite::svglite("Analysis_per_phases/MFM_RNAseqrandompicked_daynightcontrol_pvals_diurvsnocINTRONS.svg")

ggplot(melt(log10(pvalsintrons)), aes(x=c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24"), y=value))+geom_point(size=4, aes(col=rownames(melt(log10(pvalsintrons)))))+  scale_y_continuous(name="-log10 pval")+xlab("Phases")+
    geom_point(shape = 1,size=4,colour = "gray49")+
   scale_x_discrete( limits=c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24"))+ geom_hline(yintercept = 0.001, color="black", linetype="dashed")+
  ggtitle("Diurnal vs Nocturnal: Introns")+theme_bw()+scale_color_manual(values=brewer.pal(9,"RdBu"), name="Phases")
##dev.off()

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

#[1]  TRUE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE
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

##svglite::svglite("Analysis_per_phases/MFM_RNAseqrandompicked_daynightcontrol_pvals_diurvsnocALLPROMS_proportions.svg")
  
ggplot(melt(log10(pvalsintrons)), aes(x=c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24"), y=value))+geom_point(size=4, aes(col=rownames(melt(log10(pvalsintrons)))))+  scale_y_continuous(name="-log10 pval")+xlab("Phases")+
    geom_point(shape = 1,size=4,colour = "gray49")+
   scale_x_discrete( limits=c("1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24"))+ geom_hline(yintercept = 0.001, color="black", linetype="dashed")+
  ggtitle("Diurnal vs Nocturnal: All proms ")+theme_bw()+scale_color_manual(values=brewer.pal(9,"RdBu"), name="Phases")
##dev.off()
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
#nn>.001
#[1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
#nd>.001
#[1]  TRUE FALSE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE

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

#[1] FALSE FALSE FALSE FALSE FALSE FALSE
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

##svglite::svglite("Analysis_per_phases/MFM_RNAseqrandompicked_daynigh_promehn_allproms_boxplotofinteractions.svg")
#colnames(matrix.combined)<-rep(c("diurnal", "nocturnal"), 8)
ggplot(melt(rbind(diurcircproms, noccircproms)), aes(x=as.factor(Var2), y=value, fill=as.factor(Var1)))+geom_boxplot(notch = T, notchwidth = 0.5, coef=10)+
  ylab("Interactions")+
  scale_x_discrete(name="eRNAs phases")+theme_bw()+
  scale_fill_manual(values = c("brown3" , "steelblue", "brown3" , "steelblue", "brown3" , "steelblue", "brown3",   "steelblue", "brown3", "steelblue", "brown3", "steelblue", "brown3", "steelblue", "brown3",  "steelblue"),name="Promoter phases")+ggtitle("Prom-enh: All proms")
##dev.off()


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

##svglite::svglite("Analysis_per_phases/MFM_RNAseqrandompicked_daynigh_promehn_introns_boxplotofinteractions.svg")
#colnames(matrix.combined)<-rep(c("diurnal", "nocturnal"), 8)
ggplot(melt(rbind(diurcircproms, noccircproms)), aes(x=as.factor(Var2), y=value, fill=as.factor(Var1)))+geom_boxplot(notch = T, notchwidth = 0.5, coef=10)+
  ylab("Interactions")+
  scale_x_discrete(name="eRNAs phases")+theme_bw()+
  scale_fill_manual(values = c("brown3" , "steelblue", "brown3" , "steelblue", "brown3" , "steelblue", "brown3",   "steelblue", "brown3", "steelblue", "brown3", "steelblue", "brown3", "steelblue", "brown3",  "steelblue"),name="Promoter phases")+ggtitle("Prom-enh: Introns")
##dev.off()

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

#svglite::svglite("Analysis_per_phases/MFM_RNAseqrandompicked_daynigh_promehn_allproms_boxplotofinteractions_proportions.svg")
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

