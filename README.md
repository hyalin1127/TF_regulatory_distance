This script contains demosntrative codes to calculate regulatory potential conditioning on different decay distances.

# Data preprocessing: 
1. All coordinates in required files are transformed into 100bp bin with hg38 assembly. 
2. hg38_H1_bin100bp_TAD.p can be derived using provided scripts with provided files([hg38_window100bp_dict.p](https://bitbucket.org/liulab/tf_regulatory_distance/downloads/hg38_window100bp_dict.p), [hg38_H1_domain.txt](https://bitbucket.org/liulab/tf_regulatory_distance/downloads/hg38_H1_domain.txt))

#### Usage:
python Bitbucket_hg38_TAD_bin_conversion.py

# TF regulatory distance derivation

#### Required files (Demo)
1. [hg38_TSS_bin.csv](https://bitbucket.org/liulab/tf_regulatory_distance/downloads/hg38_TSS_bin.csv): TSS coordinates and corresponding bin / TAD
2. [ChIPseq_occupancy_100bp_bin.hdf5](https://bitbucket.org/liulab/tf_regulatory_distance/downloads/TEAD1_YY1_ranked_occupancy_100bp_bin.hdf5): each key/value contains sample_name/peak_bins. What provided is a demonstrative file containing TEAD1 and YY1 ChIPseqs.
3. [pearson_correlation_record_df](https://bitbucket.org/liulab/tf_regulatory_distance/downloads/pearson_correlation_record_df.p): TF-gene expression correlations
4. hg38_H1_bin100bp_TAD.p

#### Usage:
python RP_TAD_model/RP_TAD_model_bitbucket.py 

#### Output:
CSV file: index as recay distances, columns as TF samples, and cells as pearson correlation between RP and gene-TF correlation.