#####################---Plot metaplots in R---####################
#Code: Metaplots_codeforgraph_submission


library(gplots)
library(RColorBrewer)
library(gtools)

##########CTCF MESC

#CTCF mesc zt0 chrx chr19
ctcfM0<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/ctcf_mesc_zt0chrX.txt", sep=",", header = F)
ctcfM0[lower.tri(ctcfM0)]<-NA
##tiff("/metaplot_ctcfmesc_zt0_chrX.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
heatmap.2(as.matrix(ctcfM0), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="mESCs CTCF \n ZT0 chrX", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA)
#dev.off()

ctcfM0<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/1117CTCFmesc_zt0chr19", sep=",", header = F)
ctcfM0[lower.tri(ctcfM0)]<-NA
#tiff("metaplot_ctcfmesc_zt0_chr19.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
heatmap.2(as.matrix(ctcfM0), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="mESCs CTCF \n ZT0 chr19", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA)
#dev.off()

#CTCF mesc zt12 chrx 19
ctcfM12<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/ctcf_mesc_zt12chrX.txt", sep=",", header = F)
ctcfM12[lower.tri(ctcfM12)]<-NA
#tiff("metaplot_ctcfmesc_zt12_chrX.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
heatmap.2(as.matrix(ctcfM12), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="mESCs CTCF \n ZT12 chrX", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA)
#dev.off()

ctcfM12<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/1117CTCFmesc_zt12chr19.txt", sep=",", header = F)
ctcfM12[lower.tri(ctcfM12)]<-NA
#tiff("metaplot_ctcfmesc_zt12_chr19.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
heatmap.2(as.matrix(ctcfM12), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="mESCs CTCF \n ZT12 chr19", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA)
#dev.off()


################
#CTCF-liver  ENCODE and ours

#CTCF ENCONDE
#HiC ZT0 - chr19
ctcfliveENCODEzt0chr19<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/0318ctcfliverENCODE_zt0allchr19_median.txt", sep=",", header = F)
ctcfliveENCODEzt0chr19[lower.tri(ctcfliveENCODEzt0chr19)]<-NA
#tiff("metaplot_ctcfliverENCODE_zt0_chr19.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
heatmap.2(as.matrix(ctcfliveENCODEzt0chr19), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver ENCODE CTCF \n ZT0 chr19", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA)
#dev.off()

#HiC ZT12 - all chrs
ctcfliveENCODEzt12allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/0318CTCFliverENCODE1000rp_zt12allchrs_median.txt", sep=",", header = F)

#Remove cells of "genebody"
ctcfliveENCODEzt12allchrs_1000RP<-ctcfliveENCODEzt12allchrs_1000RP[c(1:11,16:25),c(1:11,16:25)]

ctcfliveENCODEzt12allchrs_1000RP[lower.tri(ctcfliveENCODEzt12allchrs_1000RP)]<-NA


#tiff("metaplot_ctcfliverENCODE_zt12_allchrs.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
heatmap.2(as.matrix(ctcfliveENCODEzt12allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver ENCODE CTCF \n ZT12 allchrs", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA)
#dev.off()

#CTCF LIVER MFM
#HiC ZT12 - chr19
ctcfliverMFM20FE_zt12chr19_500RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/0318CTCFliverMFM20FE_zt12chr19_median.txt", sep=",", header = F)
#Remove cells of "genebody"
ctcfliverMFM20FE_zt12chr19_500RP<-ctcfliverMFM20FE_zt12chr19_500RP[c(1:11,16:25),c(1:11,16:25)]

ctcfliverMFM20FE_zt12chr19_500RP[lower.tri(ctcfliverMFM20FE_zt12chr19_500RP)]<-NA

#tiff("/Figures/HiC/Metaplots/metaplot_ctcfliverMFM20EF_zt12_chr19_differentscale.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

##tiff("metaplot_ctcfliverMFM20EF_zt12_chr19.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
heatmap.2(as.matrix(ctcfliverMFM20FE_zt12chr19_500RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver MFM CTCF 500RP \n ZT12 chr19", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks=seq(.7, 1.6, .1))
#dev.off()

#HiC ZT0 - all chrs
ctcfliverMFM20FE_zt0allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/0318CTCFliverMFM20FE1000rp_zt0allchrs_median.txt", sep=",", header = F)
#Remove cells of "genebody"
ctcfliverMFM20FE_zt0allchrs_1000RP<-ctcfliverMFM20FE_zt0allchrs_1000RP[c(1:11,16:25),c(1:11,16:25)]

ctcfliverMFM20FE_zt0allchrs_1000RP[lower.tri(ctcfliverMFM20FE_zt0allchrs_1000RP)]<-NA

#tiff("/Figures/HiC/Metaplots/metaplot_ctcfliverMFM20EF_zt0_allchrs.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
##tiff("metaplot_ctcfliverMFM20EF_zt0_allchrs.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
heatmap.2(as.matrix(ctcfliverMFM20FE_zt0allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver MFM CTCF 1000RP \n ZT0 all chrs", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks = seq(.7, 1.035, .035))
#dev.off()


#HiC ZT12 - all chrs
ctcfliverMFM20FE_zt12allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/0318CTCFliverMFM20FE1000rp_zt12allchrs_median.txt", sep=",", header = F)
#Remove cells of "genebody"
ctcfliverMFM20FE_zt12allchrs_1000RP<-ctcfliverMFM20FE_zt12allchrs_1000RP[c(1:11,16:25),c(1:11,16:25)]

ctcfliverMFM20FE_zt12allchrs_1000RP[lower.tri(ctcfliverMFM20FE_zt12allchrs_1000RP)]<-NA

#tiff("/Figures/HiC/Metaplots/metaplot_ctcfliverMFM20EF_zt12_allchrs.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
##tiff("metaplot_ctcfliverMFM20EF_zt12_allchrs.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
heatmap.2(as.matrix(ctcfliverMFM20FE_zt12allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver MFM CTCF 1000RP \n ZT12 all chrs", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA,  breaks = seq(.7, 1.035, .035))
#dev.off()


#ZT12-ZT0 ALL CHRS
#tiff("/Figures/HiC/Metaplots/metaplot_ctcfliverMFM20EF_zt12-zt0_allchrs.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
heatmap.2(as.matrix(ctcfliverMFM20FE_zt12allchrs_1000RP-ctcfliverMFM20FE_zt0allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver MFM CTCF 1000RP \n ZT12-ZT0 all chrs", density.info="none",col = rev(colorRampPalette(brewer.pal(9,"RdBu"))(9)), 
          labRow = NA, labCol = NA, seq(-.5, .5, .1))
#dev.off()

##




#--------------------
#CTCF divided by timepoint

#--------CTCF ZT0--------#

#HiC ZT0 - all chrs
ctcfliverMFMZT0_zt0allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/CTCF/0518CTCFZT0MFM1000RP_HiCZT0allchrs_median.txt", sep=",", header = F)
#Remove cells of "genebody"
ctcfliverMFMZT0_zt0allchrs_1000RP<-ctcfliverMFMZT0_zt0allchrs_1000RP[c(1:11,16:25),c(1:11,16:25)]

ctcfliverMFMZT0_zt0allchrs_1000RP[lower.tri(ctcfliverMFMZT0_zt0allchrs_1000RP)]<-NA

#svg("metaplot_ctcfliverMFMZT0rp1000_hiczt0allchrs.svg")
par(cex.main=0.8)
heatmap.2(as.matrix(ctcfliverMFMZT0_zt0allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver MFM CTCF ZT0 1000RP \n  HiC ZT0 all chrs", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks = seq(.7, 1.035, .035))
#dev.off()


#HiC ZT12 - all chrs
ctcfliverMFMZT0_zt12allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/CTCF/0518CTCFZT0MFM1000RP_HiCZT12allchrs_median.txt", sep=",", header = F)
#Remove cells of "genebody"
ctcfliverMFMZT0_zt12allchrs_1000RP<-ctcfliverMFMZT0_zt12allchrs_1000RP[c(1:11,16:25),c(1:11,16:25)]

ctcfliverMFMZT0_zt12allchrs_1000RP[lower.tri(ctcfliverMFMZT0_zt12allchrs_1000RP)]<-NA

#svg("metaplot_ctcfliverMFMZT0rp1000_hiczt12allchrs.svg")
#par(cex.main=0.8)
heatmap.2(as.matrix(ctcfliverMFMZT0_zt12allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver MFM CTCF ZT0 1000RP \n  HiC ZT12 all chrs", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks = seq(.7, 1.035, .035))
#dev.off()

#ZT12-ZT0 ALL CHRS
#tiff("metaplot_ctcfliverMFMZT0rp1000_hiczt12overzt0allchrs.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
par(cex.main=0.8)
heatmap.2(as.matrix(ctcfliverMFMZT0_zt12allchrs_1000RP/ctcfliverMFMZT0_zt0allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver MFM CTCF ZT0 1000RP \n  HiC ZT12/ZT0 all chrs", density.info="none",col = rev(colorRampPalette(brewer.pal(9,"RdBu"))(9)), 
          labRow = NA, labCol = NA, breaks = seq(.95,1.05,length.out=10))

#seq(.95,1.05,length.out=10)
#seq(-.05, .05, length.out = 10)
#dev.off()


#--------CTCF ZT6--------#

#HiC ZT6 - all chrs
ctcfliverMFMZT6_zt6allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/CTCF/0518CTCFZT6MFM1000RP_HiCZT6allchrs_median.txt", sep=",", header = F)
#Remove cells of "genebody"
ctcfliverMFMZT6_zt6allchrs_1000RP<-ctcfliverMFMZT6_zt6allchrs_1000RP[c(1:11,16:25),c(1:11,16:25)]

ctcfliverMFMZT6_zt6allchrs_1000RP[lower.tri(ctcfliverMFMZT6_zt6allchrs_1000RP)]<-NA

#tiff("CTCF_ZT6_rp1000_usingHiCresol10Kb/  metaplot_ctcfliverMFMZT6rp1000_hiczt6allchrs.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
par(cex.main=0.8)
heatmap.2(as.matrix(ctcfliverMFMZT6_zt6allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver MFM CTCF ZT6 1000RP \n  HiC ZT6 all chrs", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks = seq(.7, 1.035, .035))
#dev.off()


#HiC ZT18 - all chrs
ctcfliverMFMZT6_zt18allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/CTCF/0518CTCFZT6MFM1000RP_HiCZT18allchrs_median.txt", sep=",", header = F)
#Remove cells of "genebody"
ctcfliverMFMZT6_zt18allchrs_1000RP<-ctcfliverMFMZT6_zt18allchrs_1000RP[c(1:11,16:25),c(1:11,16:25)]

ctcfliverMFMZT6_zt18allchrs_1000RP[lower.tri(ctcfliverMFMZT6_zt18allchrs_1000RP)]<-NA

#tiff("CTCF_ZT6_rp1000_usingHiCresol10Kb/  metaplot_ctcfliverMFMZT6rp1000_hiczt18allchrs.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
par(cex.main=0.8)
heatmap.2(as.matrix(ctcfliverMFMZT6_zt18allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver MFM CTCF ZT6 1000RP \n  HiC ZT18 all chrs", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks = seq(.7, 1.035, .035))
#dev.off()

#ZT18-ZT6 ALL CHRS
#tiff("CTCF_ZT6_rp1000_usingHiCresol10Kb/metaplot_ctcfliverMFMZT6rp1000_hiczt18-zt6allchrs.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
par(cex.main=0.8)
heatmap.2(as.matrix(ctcfliverMFMZT6_zt18allchrs_1000RP-ctcfliverMFMZT6_zt6allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver MFM CTCF ZT6 1000RP \n  HiC ZT18-ZT6 all chrs", density.info="none",col = rev(colorRampPalette(brewer.pal(9,"RdBu"))(9)), 
          labRow = NA, labCol = NA, breaks =seq(-.05, .05, length.out = 10) )

#seq(.95,1.05,length.out=10)
#seq(-.05, .05, length.out = 10)
#dev.off()







#--------CTCF ZT12--------#

#HiC ZT0 - all chrs
ctcfliverMFMZT12_zt0allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/CTCF/0518CTCFZT12MFM1000RP_HiCZT0allchrs_median.txt", sep=",", header = F)
#Remove cells of "genebody"
ctcfliverMFMZT12_zt0allchrs_1000RP<-ctcfliverMFMZT12_zt0allchrs_1000RP[c(1:11,16:25),c(1:11,16:25)]

ctcfliverMFMZT12_zt0allchrs_1000RP[lower.tri(ctcfliverMFMZT12_zt0allchrs_1000RP)]<-NA

#svg("CTCF_ZT12_rp1000_usingHiCresol10Kb/  metaplot_ctcfliverMFMZT12rp1000_hiczt0allchrs.svg")
par(cex.main=0.8)
heatmap.2(as.matrix(ctcfliverMFMZT12_zt0allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver MFM CTCF ZT12 1000RP \n  HiC ZT0 all chrs", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks = seq(.7, 1.035, .035))
#dev.off()


#HiC ZT12 - all chrs
ctcfliverMFMZT12_zt12allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/CTCF/0518CTCFZT12MFM1000RP_HiCZT12allchrs_median.txt", sep=",", header = F)
#Remove cells of "genebody"
ctcfliverMFMZT12_zt12allchrs_1000RP<-ctcfliverMFMZT12_zt12allchrs_1000RP[c(1:11,16:25),c(1:11,16:25)]

ctcfliverMFMZT12_zt12allchrs_1000RP[lower.tri(ctcfliverMFMZT12_zt12allchrs_1000RP)]<-NA

#svg("CTCF_ZT12_rp1000_usingHiCresol10Kb/  metaplot_ctcfliverMFMZT12rp1000_hiczt12allchrs.svg")
par(cex.main=0.8)
heatmap.2(as.matrix(ctcfliverMFMZT12_zt12allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver MFM CTCF ZT12 1000RP \n  HiC ZT12 all chrs", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks = seq(.7, 1.035, .035))
#dev.off()

#ZT12-ZT0 ALL CHRS
#tiff("CTCF_ZT12_rp1000_usingHiCresol10Kb/metaplot_ctcfliverMFMZT12rp1000_hiczt12-zt0allchrs.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
par(cex.main=0.8)
heatmap.2(as.matrix(ctcfliverMFMZT12_zt12allchrs_1000RP-ctcfliverMFMZT12_zt0allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver MFM CTCF ZT12 1000RP \n  HiC ZT12-ZT0 all chrs", density.info="none",col = rev(colorRampPalette(brewer.pal(9,"RdBu"))(9)), 
          labRow = NA, labCol = NA, breaks = seq(-.05, .05, length.out = 10))

#seq(.95,1.05,length.out=10)
#seq(-.05, .05, length.out = 10)
#dev.off()






#--------CTCF ZT18--------#

#HiC ZT6 - all chrs
ctcfliverMFMZT18_zt6allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/CTCF/0518CTCFZT18MFM1000RP_HiCZT6allchrs_median.txt", sep=",", header = F)
#Remove cells of "genebody"
ctcfliverMFMZT18_zt6allchrs_1000RP<-ctcfliverMFMZT18_zt6allchrs_1000RP[c(1:11,16:25),c(1:11,16:25)]

ctcfliverMFMZT18_zt6allchrs_1000RP[lower.tri(ctcfliverMFMZT18_zt6allchrs_1000RP)]<-NA

#tiff("CTCF_ZT18_rp1000_usingHiCresol10Kb/  metaplot_ctcfliverMFMZT18rp1000_hiczt6allchrs.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
par(cex.main=0.8)
heatmap.2(as.matrix(ctcfliverMFMZT18_zt6allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver MFM CTCF ZT18 1000RP \n  HiC ZT6 all chrs", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks = seq(.7, 1.035, .035))
#dev.off()


#HiC ZT18 - all chrs
ctcfliverMFMZT18_zt18allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/CTCF/0518CTCFZT18MFM1000RP_HiCZT18allchrs_median.txt", sep=",", header = F)
#Remove cells of "genebody"
ctcfliverMFMZT18_zt18allchrs_1000RP<-ctcfliverMFMZT18_zt18allchrs_1000RP[c(1:11,16:25),c(1:11,16:25)]

ctcfliverMFMZT18_zt18allchrs_1000RP[lower.tri(ctcfliverMFMZT18_zt18allchrs_1000RP)]<-NA

#tiff("CTCF_ZT18_rp1000_usingHiCresol10Kb/  metaplot_ctcfliverMFMZT18rp1000_hiczt18allchrs.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
par(cex.main=0.8)
heatmap.2(as.matrix(ctcfliverMFMZT18_zt18allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver MFM CTCF ZT18 1000RP \n  HiC ZT18 all chrs", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks = seq(.7, 1.035, .035))
#dev.off()

#ZT18-ZT6 ALL CHRS
#tiff("CTCF_ZT18_rp1000_usingHiCresol10Kb//metaplot_ctcfliverMFMZT18rp1000_hiczt18overzt6allchrs.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
par(cex.main=0.8)
heatmap.2(as.matrix(ctcfliverMFMZT18_zt18allchrs_1000RP/ctcfliverMFMZT18_zt6allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="Liver MFM CTCF ZT6 1000RP \n  HiC ZT18/ZT6 all chrs", density.info="none",col = rev(colorRampPalette(brewer.pal(9,"RdBu"))(9)), 
          labRow = NA, labCol = NA, breaks = seq(.95,1.05,length.out=10))

#seq(.95,1.05,length.out=10)
#seq(-.05, .05, length.out = 10)
#dev.off()




#------------------

# TADs
########## = ---------------------------------
#ZT0-TADs 50Kb Resolution 200Kb window

#HiC ZT0 - chr19
#TADsZT0MFM_zt0chr19_1000RP<-read.csv("/inputs_for_scripts/Metaplots/inputs_for_R_finalfigs/TADs/50Kbres_200Kbwindow/0418TADsZT0MFM1000rp_zt0chr19_median.txt", sep=",", header = F)

#TADsZT0MFM_zt0chr19_1000RP[lower.tri(TADsZT0MFM_zt0chr19_1000RP)]<-NA

##tiff("/Figures/HiC/Metaplots/TADs/50Kbres_200Kbwindow/metaplot_TADsZT0MFM_zt0chr19_1000RP.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

##tiff("metaplot_ctcfliverMFM20EF_zt12_chr19.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

#par(cex.main=0.8)
#heatmap.2(as.matrix(TADsZT0MFM_zt0chr19_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", main="ZT0 TADs  50Kb res\n 200Kb window MFM 100RP \n ZT0 chr19", density.info="none",col = brewer.pal(9, "BuPu"), labRow = NA, labCol = NA, breaks=seq(.7, 1.6, .1))

##dev.off()

#HiC ZT0 - all chrs
TADsZT0MFM_zt0allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/TADs/50Kbres_200Kbwindow/0418TADsZT10_50KbresolZT0HiCallchrs_median.txt", sep=",", header = F)

TADsZT0MFM_zt0allchrs_1000RP[lower.tri(TADsZT0MFM_zt0allchrs_1000RP)]<-NA

#svglite::svglite(filename = "metaplot_TADsZT0_zt050kballchrs_1000RP.svg")

par(cex.main=0.8)
heatmap.2(as.matrix(TADsZT0MFM_zt0allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="ZT0 TADs  50Kb res\n 200Kb window MFM 1000RP \n ZT0 allchrs 50 Kb resol", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks=seq(.65, 1.1, .05))
#dev.off()


##HiC ZT12 - chr19
#TADsZT0MFM_zt12chr19_1000RP<-read.csv("/inputs_for_scripts/Metaplots/inputs_for_R_finalfigs/TADs/50Kbres_200Kbwindow/0418TADsZT0MFM1000rp_zt12chr19_median.txt", sep=",", header = F)

#TADsZT0MFM_zt12chr19_1000RP[lower.tri(TADsZT0MFM_zt12chr19_1000RP)]<-NA

##tiff("/Figures/HiC/Metaplots/TADs/50Kbres_200Kbwindow/metaplot_TADsZT0MFM_zt12chr19_1000RP.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

##tiff("metaplot_ctcfliverMFM20EF_zt12_chr19.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

#par(cex.main=0.8)
#heatmap.2(as.matrix(TADsZT0MFM_zt12chr19_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", main="ZT0 TADs  50Kb res\n 200Kb window MFM 100RP \n ZT12 chr19", density.info="none",col = brewer.pal(9, "BuPu"), labRow = NA, labCol = NA, breaks=seq(.7, 1.6, .1))

##dev.off()

##HiC ZT12 - all chrs
TADsZT0MFM_zt12allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/TADs/50Kbres_200Kbwindow/0418TADsZT10_50KbresolZT12HiCallchrs_median.txt", sep=",", header = F)

TADsZT0MFM_zt12allchrs_1000RP[lower.tri(TADsZT0MFM_zt12allchrs_1000RP)]<-NA

#svglite::svglite("metaplot_TADsZT0_zt1250kballchrs_1000RP.svg")

##tiff("metaplot_ctcfliverMFM20EF_zt12_chr19.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

par(cex.main=0.8)
heatmap.2(as.matrix(TADsZT0MFM_zt12allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="ZT0 TADs  50Kb res\n 200Kb window MFM 1000RP \n ZT12 allchrs 50 Kb resol", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks=seq(.65, 1.1, .05))

#dev.off()

#Difference  hic ZT12-ZT0 ALL CHRS
#tiff("metaplot_TADsZT0_zt12-zt050kballchrs_1000RP.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

##tiff("metaplot_ctcfliverMFM20EF_zt12_chr19.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

par(cex.main=0.8)
heatmap.2(as.matrix(TADsZT0MFM_zt12allchrs_1000RP-TADsZT0MFM_zt0allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="ZT0 TADs  50Kb res\n 200Kb window MFM 1000RP \n  HiC ZT12-ZT0 allchrs 50 Kb resol", density.info="none",col = rev(colorRampPalette(brewer.pal(9,"RdBu"))(9)), 
          labRow = NA, labCol = NA, breaks =  seq(-.05, .05, length.out = 10) )
            #seq(.95,1.05,length.out=10) seq(-.05, .05, length.out = 10)
          

#dev.off()


#ZT16-TADs 50Kb Resolution 200Kb window


#HiC ZT6 - allchrs

TADsZT6MFM_zt6allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/TADs/50Kbres_200Kbwindow/0418TADsZT6_50KbresolZT6HiCallchrs_median.txt", sep=",", header = F)

TADsZT6MFM_zt6allchrs_1000RP[lower.tri(TADsZT6MFM_zt6allchrs_1000RP)]<-NA

#tiff("metaplot_TADsZT6_zt650kballchrs_1000RP.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

#svglite::svglite("HiC_50Kbresolution/metaplot_TADsZT6_zt650kballchrs_1000RP.svg")
par(cex.main=0.8)
heatmap.2(as.matrix(TADsZT6MFM_zt6allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="ZT6 TADs  50Kb res\n 200Kb window MFM 1000RP \n ZT6 allchrs 50 Kb resol", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks=seq(.65, 1.1, .05))

#dev.off()


#HiC ZT18 - allchrs

TADsZT6MFM_zt18allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/TADs/50Kbres_200Kbwindow/0418TADsZT6_50KbresolZT18HiCallchrs_median.txt", sep=",", header = F)

TADsZT6MFM_zt18allchrs_1000RP[lower.tri(TADsZT6MFM_zt18allchrs_1000RP)]<-NA

#tiff("metaplot_TADsZT6_zt1850kballchrs_1000RP.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)
#svglite::svglite("HiC_50Kbresolution/metaplot_TADsZT6_zt1850kballchrs_1000RP.svg")

par(cex.main=0.8)
heatmap.2(as.matrix(TADsZT6MFM_zt18allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="ZT6 TADs  50Kb res\n 200Kb window MFM 1000RP \n ZT18 allchrs 50 Kb resol", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks=seq(.65, 1.1, .05))

#dev.off()


# Difference HiC ZT18 - ZT6 allchrs
#tiff("metaplot_TADsZT6_ZT18overZT650kballchrs_1000RP.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

par(cex.main=0.8)
heatmap.2(as.matrix(TADsZT6MFM_zt18allchrs_1000RP/TADsZT6MFM_zt6allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="ZT6 TADs  50Kb res\n 200Kb window MFM 1000RP \n ZT18/ZT6 allchrs 50 Kb resol", density.info="none",col = rev(colorRampPalette(brewer.pal(9,"RdBu"))(9)),
          labRow = NA, labCol = NA, breaks=seq(.95,1.05,length.out=10))

#seq(.95,1.05,length.out=10)
#seq(-.05, .05, length.out = 10)

#dev.off()



#ZT12-TADs 50Kb Resolution 200Kb window

#HiC ZT0 - chr19
#TADsZT12MFM_zt0chr19_1000RP<-read.csv("/inputs_for_scripts/Metaplots/inputs_for_R_finalfigs/TADs/50Kbres_200Kbwindow/0418TADsZT12MFM1000rp_zt0chr19_median.txt", sep=",", header = F)

#TADsZT12MFM_zt0chr19_1000RP[lower.tri(TADsZT12MFM_zt0chr19_1000RP)]<-NA

##tiff("/Figures/HiC/Metaplots/TADs/50Kbres_200Kbwindow/metaplot_TADsZT12MFM_zt0chr19_1000RP.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

##tiff("metaplot_ctcfliverMFM20EF_zt12_chr19.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

#par(cex.main=0.8)
#heatmap.2(as.matrix(TADsZT12MFM_zt0chr19_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", main="ZT12 TADs  50Kb res\n 200Kb window MFM 100RP \n ZT0 chr19", density.info="none",col = brewer.pal(9, "BuPu"), labRow = NA, labCol = NA, breaks=seq(.7, 1.6, .1))
##dev.off()

#HiC ZT0 - allchrs
TADsZT12MFM_zt0allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/TADs/50Kbres_200Kbwindow/0418TADsZT12_50KbresolZT0HiCallchrs_median.txt", sep=",", header = F)

TADsZT12MFM_zt0allchrs_1000RP[lower.tri(TADsZT12MFM_zt0allchrs_1000RP)]<-NA

#svg("HiC_50Kbresolution/metaplot_TADsZT12_zt050kballchrs_1000RP.svg")

##tiff("metaplot_ctcfliverMFM20EF_zt12_chr19.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

par(cex.main=0.8)
heatmap.2(as.matrix(TADsZT12MFM_zt0allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="ZT12 TADs  50Kb res\n 200Kb window MFM 1000RP \n ZT0 allchrs 50 Kb resol", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks=seq(.65, 1.1, .05))

#dev.off()



##HiC ZT12 - chr19
#TADsZT12MFM_zt12chr19_1000RP<-read.csv("/inputs_for_scripts/Metaplots/inputs_for_R_finalfigs/TADs/50Kbres_200Kbwindow/0418TADsZT12MFM1000rp_zt12chr19_median.txt", sep=",", header = F)

#TADsZT12MFM_zt12chr19_1000RP[lower.tri(TADsZT12MFM_zt12chr19_1000RP)]<-NA

##tiff("/Figures/HiC/Metaplots/TADs/50Kbres_200Kbwindow/metaplot_TADsZT12MFM_zt12chr19_1000RP.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

##tiff("metaplot_ctcfliverMFM20EF_zt12_chr19.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

#par(cex.main=0.8)
#heatmap.2(as.matrix(TADsZT12MFM_zt12chr19_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", main="ZT0 TADs  50Kb res\n 200Kb window MFM 100RP \n ZT12 chr19", density.info="none",col = brewer.pal(9, "BuPu"),labRow = NA, labCol = NA, breaks=seq(.7, 1.6, .1))

##dev.off()

#HiC ZT12 - allchrs
TADsZT12MFM_zt12allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/TADs/50Kbres_200Kbwindow/0418TADsZT12_50KbresolZT12HiCallchrs_median.txt", sep=",", header = F)

TADsZT12MFM_zt12allchrs_1000RP[lower.tri(TADsZT12MFM_zt12allchrs_1000RP)]<-NA

#svg("HiC_50Kbresolution/metaplot_TADsZT12_zt1250kballchrs_1000RP.svg")

##tiff("metaplot_ctcfliverMFM20EF_zt12_chr19.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

par(cex.main=0.8)
heatmap.2(as.matrix(TADsZT12MFM_zt12allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="ZT12 TADs  50Kb res\n 200Kb window MFM 1000RP \n ZT12 allchrs 50 Kb resol", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks=seq(.65, 1.1, .05))

#dev.off()

#Difference  hic ZT12-ZT0 ALL CHRS
#tiff("metaplot_TADsZT12_zt12-zt050kballchrs_1000RP.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

##tiff("metaplot_ctcfliverMFM20EF_zt12_chr19.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

par(cex.main=0.8)
heatmap.2(as.matrix(TADsZT12MFM_zt12allchrs_1000RP-TADsZT12MFM_zt0allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="ZT12 TADs  50Kb res\n 200Kb window MFM 1000RP \n  HiC ZT12-ZT0 allchrs 50 Kb resol", density.info="none",col = rev(colorRampPalette(brewer.pal(9,"RdBu"))(9)), 
          labRow = NA, labCol = NA, breaks = seq(-.05, .05, length.out = 10))

#seq(.95,1.05,length.out=10)
#dev.off()





#ZT18-TADs 50Kb Resolution 200Kb window


#HiC ZT6 - allchrs
TADsZT18MFM_zt6allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/TADs/50Kbres_200Kbwindow/0418TADsZT18_50KbresolZT6HiCallchrs_median.txt", sep=",", header = F)

TADsZT18MFM_zt6allchrs_1000RP[lower.tri(TADsZT18MFM_zt6allchrs_1000RP)]<-NA

#tiff("metaplot_TADsZT18_zt650kballchrs_1000RP.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

#svglite::svglite("HiC_50Kbresolution/metaplot_TADsZT18_zt650kballchrs_1000RP.svg")

##tiff("metaplot_ctcfliverMFM20EF_zt12_chr19.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

par(cex.main=0.8)
heatmap.2(as.matrix(TADsZT18MFM_zt6allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="ZT18 TADs  50Kb res\n 200Kb window MFM 1000RP \n ZT6 allchrs 50 Kb resol", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks=seq(.65, 1.1, .05))

#dev.off()


#HiC ZT18 - allchrs
TADsZT18MFM_zt18allchrs_1000RP<-read.csv("/inputs_for_scripts/HiC/Metaplots/inputs_for_R_finalfigs/TADs/50Kbres_200Kbwindow/0418TADsZT18_50KbresolZT18HiCallchrs_median.txt", sep=",", header = F)

TADsZT18MFM_zt18allchrs_1000RP[lower.tri(TADsZT18MFM_zt18allchrs_1000RP)]<-NA

#tiff("metaplot_TADsZT18_zt1850kballchrs_1000RP.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)


#svglite::svglite("HiC_50Kbresolution/metaplot_TADsZT18_zt1850kballchrs_1000RP.svg")


##tiff("metaplot_ctcfliverMFM20EF_zt12_chr19.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

par(cex.main=0.8)
heatmap.2(as.matrix(TADsZT18MFM_zt18allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="ZT18 TADs  50Kb res\n 200Kb window MFM 1000RP \n ZT18 allchrs 50 Kb resol", density.info="none",col = brewer.pal(9, "BuPu"), 
          labRow = NA, labCol = NA, breaks=seq(.65, 1.1, .05))

#dev.off()

#Difference  hic ZT18-ZT6 ALL CHRS
#tiff("metaplot_TADsZT18_zt18overzt650kballchrs_1000RP.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

##tiff("metaplot_ctcfliverMFM20EF_zt12_chr19.tiff", height = 12, width = 12, units = 'cm', compression = "lzw", res = 300)

par(cex.main=0.8)
heatmap.2(as.matrix(TADsZT18MFM_zt18allchrs_1000RP/TADsZT18MFM_zt6allchrs_1000RP), key.xlab = "Median OE", key.title = "", Colv = F, Rowv = F,  trace="none", 
          main="ZT18 TADs  50Kb res\n 200Kb window MFM 1000RP \n  HiC ZT18/ZT6 allchrs 50 Kb resol", density.info="none",col = rev(colorRampPalette(brewer.pal(9,"RdBu"))(9)), 
          labRow = NA, labCol = NA, breaks = seq(.95,1.05,length.out=10))

#  seq(-.05, .05, length.out = 10)
#dev.off()

