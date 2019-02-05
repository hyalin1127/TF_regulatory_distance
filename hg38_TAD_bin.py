import os
import pandas as pd
import math
import pickle
import itertools
import numpy as np
import h5py
import scipy
from collections import Counter, defaultdict
import operator
import itertools
from scipy import stats

target_path = os.getcwd()

def hg38_TAD_100bpbin_conversion():
    TAD_bin=defaultdict(list)
    bin_TAD=dict()

    hg38_100bpbin_dict = pickle.load(open("hg38_window100bp_dict.p",'rb'))
    file=open("hg38_H1_domain.txt")
    count=0
    for line in file:
        count+=1
        elements=line.strip().split("\t")
        chromosome=elements[0]
        start=int(((int(elements[1]))/100)+1)
        end=int(((int(elements[2]))/100))
        for i in range(start,end):
            bins=hg38_100bpbin_dict[chromosome][3]+i
            TAD_bin["hg38_H1_%s" %(str(count))].append(bins)
            bin_TAD[bins]="hg38_H1_%s" %(str(count))

    pickle.dump(bin_TAD,file=open("hg38_H1_bin100bp_TAD.p",'wb'))
    pickle.dump(TAD_bin,file=open("hg38_H1_TAD_bin100bp.p",'wb'))

def main():
    hg38_TAD_100bpbin_conversion()

main()
