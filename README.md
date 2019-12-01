This script contains demosntrative codes to calculate:
1) regulatory potentials 
2) decay distances-depndent correlation between RP and TF-gene associations.

# Requirements:
Python >= 3.6  
Pandas >= 0.24.2  
h5py >= 2.0.0  
Scipy >= 1.1.0  
Numpy >= 1.14.2  

# Installation:
git clone https://hyalin1127@github.org/tf_regulatory_distance.git

# Required data preprocessing:

1. All coordinates in required files are transformed into 100bp bin with hg38 assembly. 
2. hg38_H1_bin100bp_TAD.p can be derived using provided scripts with provided files([hg38_window100bp_dict.p](https://bitbucket.org/liulab/tf_regulatory_distance/downloads/hg38_window100bp_dict.p), [hg38_H1_domain.txt](https://bitbucket.org/liulab/tf_regulatory_distance/downloads/hg38_H1_domain.txt))

#### Usage:
python hg38_TAD_bin.py 

# TF regulatory distance derivation #

#### Explanation:
This script will perform the following tasks:  
1) calculates regulatory potential (RP) of transcription faactors.  
2) derive the regulatory distance-dependent correlation between RP and gene-TF correlation.  

#### Required files (Demo):

1. [hg38_TSS_bin.csv](https://bitbucket.org/liulab/tf_regulatory_distance/downloads/hg38_TSS_bin.csv): TSS coordinates and corresponding bin / TAD
2. [ChIPseq_occupancy_100bp_bin.hdf5](https://bitbucket.org/liulab/tf_regulatory_distance/downloads/TEAD1_YY1_ranked_occupancy_100bp_bin.hdf5): each key/value contains sample_name/peak_bins. What provided is a demonstrative file containing TEAD1 and YY1 ChIPseqs.
3. [pearson_correlation_record_df](https://bitbucket.org/liulab/tf_regulatory_distance/downloads/pearson_correlation_record_df.p): TF-gene expression correlations
4. hg38_H1_bin100bp_TAD.p: Dictionary with bin as keys and TAD as values. Can be derived in 'Data Preprocessing' section.

#### Arguments: 

```
-T TF_OCCUPANCY_FILE, --TFChipSeq=TF_OCCUPANCY_FILE
                        Name of TF_occupancy_file
-C CORRELATION_FILE, --Correlation=CORRELATION_FILE
                        Name of TF-gene icorrelation_file
-S HG38_TSS_INFO, --TSS=HG38_TSS_INFO
                        Name of hg38_TSS_info file
-A HG38_BIN100BP_TAD_INFO, --TAD=HG38_BIN100BP_TAD_INFO
                        hg38_bin100bp_TAD_info
```

#### Usage (Demo):
python RP_TAD_model.py -T TEAD1_YY1_ranked_occupancy_100bp_bin.hdf5 -C pearson_correlation_record_df.p -S hg38_TSS_bin.csv -A hg38_H1_bin100bp_TAD.p

#### Output:

CSV file: index as recay distances, columns as TF samples, and cells as pearson correlation between RP and gene-TF expression correlation.

### Citations:


