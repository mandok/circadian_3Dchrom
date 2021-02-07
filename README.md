# circadian_3Dchrom
This repository contains all the in-house scripts and inputs for this publication: https://doi.org/10.1101/2020.07.23.217992


The following lines are needed to run before executing any script
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
