#!/usr/bin/env python
import os
import sys
import pandas as pd
import pickle
import itertools
import numpy as np
import h5py
from collections import Counter, defaultdict
import operator
import itertools
from scipy import stats
from optparse import OptionParser

def Info(infoStr):
    print("[%s] %s" %(time.strftime('%H:%M:%S'), infoStr))

def within_TAD_RP_model(RP_distance_constant,bindings_bins,hg38_TSS_bin_list,H1_bin100bp_TAD):
    gene_RP_record =defaultdict(list)

    binding_bins_in_TAD = defaultdict(list)
    for bins in bindings_bins:
        try:
            TAD = H1_bin100bp_TAD[bins]
            binding_bins_in_TAD[TAD].append(float(bins)/RP_distance_constant)
        except:
            pass

    for gene,bins,TAD in hg38_TSS_bin_list:
        bins = float(bins)/RP_distance_constant
        within_TAD_bindings_bins = binding_bins_in_TAD[TAD]
        RP = np.sum([np.exp(-abs(bindings_bin - bins)) for bindings_bin in within_TAD_bindings_bins])
        gene_RP_record[gene].append(RP)
    return_gene_RP_record = {i:np.mean(j) for i,j in gene_RP_record.items()}
    return(return_gene_RP_record)

def prepare_optparser():
    usage = "usage: python %prog ......"
    description = "Input a TF ChIPseq file and TF-gene expression correlation, output correlation between RP and TF-gene correlation."
    description = "Demo command: python RP_TAD_model.py -T TEAD1_YY1_ranked_occupancy_100bp_bin.hdf5 -C pearson_correlation_record_df.p -S hg38_TSS_bin.csv -A hg38_H1_bin100bp_TAD.p"
    optparser = OptionParser(version="%prog v1.00", description=description, usage=usage, add_help_option=False)

    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-T","--TFChipSeq",dest="TF_occupancy_file",type="string",
                         help="Name of TF_occupancy_file")
    optparser.add_option("-C","--Correlation",dest="correlation_file",type="string",
                         help="Name of TF-gene icorrelation_file")
    optparser.add_option("-S","--TSS",dest="hg38_TSS_info",type="string",
                         help="Name of hg38_TSS_info file")
    optparser.add_option("-A","--TAD",dest="hg38_bin100bp_TAD_info",type="string",
                         help="hg38_bin100bp_TAD_info")

    (options,args) = optparser.parse_args()
    return(options)

class RP():
    def __init__(self, options):
        self.TF_occupancy_100bp_bin = h5py.File(options.TF_occupancy_file, "r")
        self.pearson_correlation_record_df = pd.read_pickle(options.correlation_file)
        self.hg38_TSS_bin_df = pd.read_csv(options.hg38_TSS_info,sep="\t",header=0,index_col=None)
        self.hg38_TSS_bin_list = (self.hg38_TSS_bin_df[["gene","bin","TAD"]]).values.tolist()
        self.H1_bin100bp_TAD = pickle.load(open(options.hg38_bin100bp_TAD_info,'rb'))
        self.RP_distance_constants = [1,5,10,30,50,80,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,
        1800,1900,2000,2200,2300,2500,2700,3000,3200,3600,3800,4000,4200,4400,4600,5000,5500,6000,7000,8000,10000,15000,20000]
        self.target_path = os.getcwd()

    def TF_RP_vs_gene_complete(self):
        samples = list(self.TF_occupancy_100bp_bin.keys())
        for RP_distance_constant in self.RP_distance_constants:
            total_gene_RP_record = dict()

            for sample in samples:
                bindings_bins = list(self.TF_occupancy_100bp_bin[sample][()])
                total_gene_RP_record[sample] = within_TAD_RP_model(RP_distance_constant,bindings_bins,self.hg38_TSS_bin_list,self.H1_bin100bp_TAD)

            total_gene_RP_record_df = pd.DataFrame(total_gene_RP_record)
            pickle.dump(total_gene_RP_record_df,file=open("/%s/total_gene_RP_record_df_constant_%s_complete.p" %(self.target_path,str(RP_distance_constant)),'wb'))

    def TF_RP_vs_gene_complete_processing_all_TADs(self):
        genes_in_pearson_correlation_record_df = set(self.pearson_correlation_record_df.index.tolist())

        stat_record = defaultdict(dict)
        for RP_distance_constant in self.RP_distance_constants:
            total_gene_RP_record_df = pickle.load(open("/%s/total_gene_RP_record_df_constant_%s_complete.p" %(self.target_path,str(RP_distance_constant)),'rb'))
            samples = list(total_gene_RP_record_df.columns.values)
            common_genes = [i for i in total_gene_RP_record_df.index.tolist() if i in genes_in_pearson_correlation_record_df]
            selected_total_gene_RP_record_df = total_gene_RP_record_df.loc[common_genes]
            selected_pearson_correlation_record_df = self.pearson_correlation_record_df.loc[common_genes]

            for sample in samples:
                TF = (sample.split("_"))[0]
                if TF in selected_pearson_correlation_record_df.columns.values:
                    expression_correlation = list(selected_pearson_correlation_record_df[TF])
                    RP = list(selected_total_gene_RP_record_df[sample])
                    temp_df = pd.DataFrame()
                    temp_df["exp"] = expression_correlation
                    temp_df["RP"] = RP
                    temp_df.replace([np.inf, -np.inf], np.nan,inplace=True)
                    temp_df.dropna(axis=0,how="any",inplace=True)

                    stat_record[sample][RP_distance_constant*100*np.log(2)] = stats.pearsonr(list(temp_df["RP"]), list(temp_df["exp"]))[0]

        stat_record_df = pd.DataFrame(stat_record)
        stat_record_df.to_csv("/%s/TF_RP_all_TADs_correlation_CCLE.csv" %(self.target_path))

def main():
    opts = prepare_optparser()
    g = RP(opts)
    g.TF_RP_vs_gene_complete()
    g.TF_RP_vs_gene_complete_processing_all_TADs()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) Bye!\n")
        sys.exit(0)
