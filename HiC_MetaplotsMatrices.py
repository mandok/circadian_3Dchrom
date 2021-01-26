#####################---Metaplots---#############
#Scripts ran in Jupyter notebook, Python 2.7
##Libraries

%matplotlib inline
import sys, getopt
import math
import numpy as np
import pandas as pd
from scipy import stats, integrate,sparse
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import coo_matrix,csr_matrix
#import straw
import time
from itertools import combinations
from collections import Counter
import itertools
from scipy import stats
import operator
import csv
from itertools import groupby
from math import floor

#### Definition of the functions used in the script

#To calculate the metaplot we segmented the triangular matrix in 4 areas: 
#1) area of the left wing (x:1-10,y:1-10) with the  retrievehiclwing function
#2) area of the right wing (x:16-20,y:16-20) with the  retrievehicrwing function
#3) area of top (x:1-11,y:16-20) with retrievehicfeattopright
#4) are of regions of interest (x:11-15,y:11-15) with retrievehicfeatmiddle
#5) area of region of interest vs the wing (x:1-15,y:11-20) with retrievehicfeatmiwings


def retrievehiclwing(minposlwing , maxposlwing):
    tlwing = ((hiczt0[(hiczt0.bin1 >= minposlwing ) & (hiczt0.bin1 <= maxposlwing ) & (hiczt0.bin2 <= maxposlwing) &  (hiczt0.bin2 >= minposlwing) ])) # filter the hic bins that are within the range of the positions of the defined leftwing
    bin1 = (tlwing.bin1).values - minposlwing #change the genomic position into the index inside the matrix; e.g.  genomic positions of 10Kb bins  becomes the x in the matrix
    bin2 = (tlwing.bin2).values - minposlwing #change the genomic position into the index inside the matrix; e.g.  genomic positions of 10Kb bins  becomes the y in the matrix
    oe = (tlwing.OE).values #values retrived from the hic correspoding to the x,y positions (bin1 bin2)
    
    for cell in range(len(oe)): # bin1 and bin2 now are in the range of 1-10, hence the matrix can be filled
        matrixmeta [bin1[cell]][bin2[cell]] =  oe[cell]
    
    return(matrixmeta)


def retrievehicrwing(minposrwing , maxposrwing):
        trwing = ((hiczt0[(hiczt0.bin1 >= minposrwing ) & (hiczt0.bin1 <= maxposrwing ) & (hiczt0.bin2 <= maxposrwing) &  (hiczt0.bin2 >= minposrwing) ])) # filter the hic bins that are within the range of the positions of the defined rightwing
        #print(tlwing)
        bin1 = (trwing.bin1).values - minposrwing + 15 #change the genomic position into the index inside the matrix; e.g.  genomic positions of 10Kb bins  becomes the x in the matrix
        bin2 = (trwing.bin2).values - minposrwing + 15 #change the genomic position into the index inside the matrix; e.g.  genomic positions of 10Kb bins  becomes the y in the matrix
        #print([bin1, bin2])   
        oe = (trwing.OE).values #values retrived from the hic correspoding to the x,y positions (bin1 bin2)

        for cell in range(len(oe)): # bin1 and bin2 now are in the range of 15-20, hence the matrix can be filled
            matrixmeta [bin1[cell]][bin2[cell]] =  oe[cell]
        return(matrixmeta)

def retrievehicfeattopright(pos, minposlwing, maxposrwing, minposrwing):
        #print((hiczt0.bin1))
        t = hiczt0[((hiczt0.bin1 < pos[0] ) & (hiczt0.bin1 >= minposlwing )& (hiczt0.bin2 >= minposrwing ) & (hiczt0.bin2 <= maxposrwing ) )] #!!!!Modify

        bin1 = (t.bin1).values #- minposlwing
        bin2 = (t.bin2).values #- minposlwing

        bin1 = bin1- minposlwing
        bin2 = bin2 - minposlwing
        bin2 = bin2 - (minposrwing - minposlwing -15)
        
        oe = (t.OE).values
        #print([bin1, bin2, oe]) 
        
        for cell in range(len(oe)):
            matrixmeta [bin1[cell]][bin2[cell]] =  oe[cell]
        return(matrixmeta)

def retrievehicfeatmiddle(pos, minposlwing, maxposrwing, minposrwing):
        i=0
        while(i <=4):
            #Only diag
            t = hiczt0[((hiczt0.bin1 >= pos[i] ) & (hiczt0.bin1 <= pos[i+1] )& (hiczt0.bin2 >= pos[i] ) & (hiczt0.bin2 <= pos[i+1]) )] #!!!!Modify
            oe = (t.OE.median())
            matrixmeta [10+i][10+i] =  oe
            i= i+1

        #other combinations 
        for j in range(1,5, 1):
            i=0
            while(i < (len(pos)-j -1)):
                t = hiczt0[((hiczt0.bin1 >= pos[i] ) & (hiczt0.bin1 <= pos[i+1] )& (hiczt0.bin2 > pos[i+j] ) & (hiczt0.bin2 <= pos[i+j+1]) )]
                oe = (t.OE.median())
                matrixmeta [10+i][10+i+j] =  oe
                i= i+1


        return(matrixmeta) 

def retrievehicfeatmiwings(pos, minposlwing, maxposlwing, minposrwing, maxposrwing):
        intervalslw=(np.linspace(minposlwing, maxposlwing, 10)).astype(int)
        intervalsrw=(np.linspace(minposrwing, maxposrwing, 10)).astype(int)

        j=0
        for index in range(10):
            i=0
            while(i< 5):
                #midtop
                t = hiczt0[((hiczt0.bin1 == intervalslw[index] ) )& (hiczt0.bin2 >= pos[i] ) & (hiczt0.bin2 <= pos[i+1] )  ]
                oe = (t.OE.median())
                matrixmeta [j][10+i] =  oe
                #midleft
                t = hiczt0[((hiczt0.bin2 == intervalsrw[index] ) )& (hiczt0.bin1 >= pos[i] ) & (hiczt0.bin1 <= pos[i+1] )  ]
                oe = (t.OE.median())
                matrixmeta [10+i][15+j] =  oe
                i= i+1
            j=j+1
            
        return(matrixmeta)
#Open region of interest/anchors. Format: chr start end sign
tads="TADs/0100_ZT12_ABC_TADs_TADtool_ws200kb_co140_mod.bed"
#ctcf="CTCF/Peaks/Shared_BUTNOT_ZT18_rep1_rep2_modmetaplotinput.txt"

#numebr of anchors (minimum 2)
anchornum = 2

#number of bins between each anchor
binnum = 5

#chromosome size file
chrsizefile = "inputs_for_scripts/chrsizes_mm9"


#resolution of the interaction file in kb
resolution = 10
strres = str(resolution)
res = resolution * 1000

#option for toggling wings
wings = True

chrlist = []
anchorlist = []

#totalbinnum = binnum*(anchornum-1)

newresults = []

pos = []
posrwing = []
poslwing = []
chrsizes = []
sign = []

#Create list of chrsnames
chrlist = open(chrsizefile, 'r')
chrsizes = [ chrs.split()[0] for chrs in chrlist]

tic = time.clock()
#Iterate over chr name
for chrs in chrsizes[:19]:
    print(chrs)
    anchors = open(tads,'r')
    for line in anchors:
        li = line.split()

        
        if(li[0]==chrs):       #Create the bins for each feature 
            ini = int(li[1])//res
            fin = int(li[2])//res 
            #t = (stats.binned_statistic([f for f in range(ini, fin+1, 1)], [f for f in range(ini, fin+1, 1)], bins=binnum-1)[1])
            #t = [int(f) for f in t]
            #Bins for the 
            intervals=(np.linspace(ini, fin, 6)).astype(int)
            pos.append(intervals)
            
            
            lwing= [ini-w for w in range(10, 0, -1)] #10 bins as left wing  10*res = number of bp for each wing 
            poslwing.append(lwing)
            rwing= [fin+w for w in range(1,11, 1)] #10 bins as right wing  10*res = number of bp for each wing   
            posrwing.append(rwing)
            #pos.append(rwing+feat+lwing)
        
            #counter+=1
            #Store the sign of the feature
            sign.append(li[3])
    #Open hic files; output from juicer bin1 bin2 value(obs/exp); bins are in 10Kb bins; header! bin1 bin2 OE  
    ##IMPORTANT: to change of timepoint of the HiC change manually this part
    # ../Metaplots/ZT0/ZT0ABC_chr_arm_%.... to the corresponding timepoint e.g. ../Metaplots/ZT6/ZT6ABC_chr_arm_%....
    hiczt0name = 'inputs_for_scripts/Metaplots/ZT12/ZT12ABC_chr_arm_%s_10kbKR_oe_binnedmtplotinput.txt' % chrs

    hiczt0 = pd.read_csv(hiczt0name, header=0, sep="\t") #Open hic OE ZT0 for the current chr
   

    for feat in range(sum([1 for p in pos])): #Iterate over feature positions
    #for feat in range(2):
        matrixmeta = [[0 for x in range(25)] for y in range(25)] #Create matrix
        col = 0
        row = 0
        minposlwing = min(poslwing[feat])
        maxposlwing = max(poslwing[feat])
        minposrwing = min(posrwing[feat])
        maxposrwing = max(posrwing[feat])

        #Fill lefwing
        #print([minposlwing , maxposlwing ,maxposlwing-minposlwing ])
        retrievehiclwing(minposlwing, maxposlwing)


        #Fill rightwing
        #print([minposrwing , maxposrwing ,maxposrwing-minposrwing ])
        retrievehicrwing(minposrwing , maxposrwing)
        #Fill feature
        #print([pos[feat] ,pos[feat][4]-pos[feat][0] ])
        retrievehicfeattopright(pos[feat], minposlwing, maxposrwing, minposrwing)
        retrievehicfeatmiddle(pos[feat], minposlwing, maxposrwing, minposrwing)
        retrievehicfeatmiwings(pos[feat], minposlwing, maxposlwing, minposrwing, maxposrwing)
        
        #Flip matrix if it is 
        if(sign[feat]=="-"):
            t = np.rot90(np.flip(matrixmeta, 1),3)
            newresults.append(t) #Save results
        else:  
            newresults.append(matrixmeta) #Save results in a list, the lenght of the list represents the number of regions of interest. Inside list[0] there will the positions of the matrix and the hic value
        
toc = time.clock()
print(toc - tic )

##Calculate the mean, sum and median all of the matrices to build the metamatrix

#np.median([newresults[feat][0][0] for feat in range(57)])
finalresults_mean = [[0 for x in range(25)] for y in range(25)]
finalresults_median = [[0 for x in range(25)] for y in range(25)]
finalresults_sum = [[0 for x in range(25)] for y in range(25)]
for feat in range(sum([1 for p in pos])):
#for feat in range(2):
    for col in range(25):
        for row in range(25):   
            #finalresults_median[col][row]=np.median([newresults[feat][col][row] for feat in range(2)])

            finalresults_mean[col][row]=np.nanmean([newresults[feat][col][row] for feat in range(sum([1 for p in pos]))], )
            finalresults_median[col][row]=np.nanmedian([newresults[feat][col][row] for feat in range(sum([1 for p in pos]))])
            finalresults_sum[col][row]=np.nansum([newresults[feat][col][row] for feat in range(sum([1 for p in pos]))])

#Quick plot to look at the results of the median metamatrix
sns.heatmap(finalresults_median, cmap="YlGnBu")
plt.axvline(x=10, color="w")
plt.axvline(x=15, color="w")
plt.axhline(y=10, color="w")
plt.axhline(y=15, color="w")


#Save results; change manually the name of the file and the variable (either finalresults_median, finalresults_mean or finalresults_sum)

with open("TADsZT12MFMcicpromsintronsZT12MFM_ZT12HiCallchrs_median.txt", "w") as f:
    writer = csv.writer(f)
    writer.writerows(finalresults_median)

