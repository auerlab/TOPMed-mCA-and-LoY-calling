#!/bin/bash

############################################################################################################
#
#   Objective: Integrated TopMed WGS Loss of Y mCA calling pipeline from raw vcf data to mCA instances using MoChA.
#
#              This script was designed to filter the mCA for Loss of Y.
#              
#   NOTE: Please refer to https://github.com/freeseek/mocha for GRCh38 resources.
#         In this pipeline, we assume they are ready to be applied.
#
#         This pipeline is contingent to be tweaked according to your research requirements.
#
#         All required/updated softwares should be installed before running this pipeline; otherwise,
#         it will cause errors. 
#
#         bcftools(1.11 or later), shapeit4(4.2.2 or later, if VCF needs phasing), bcftools-mocha(1.15 or later), 
#         bedtools(2.30 or later)
#          
#
############################################################################################################
#   
#   Revision History: 
#
#   Date              Author               Title                       Modification
#   09-24-2023        Mr.Xiaolong Ma       Bioinformatics Analyst      v1
#   
#   Contact                       Institute
#   xima@mcw.edu      Medical College of Wisconsin
#
############################################################################################################
#
#   Directory Hierarchy
#
#   your_proj_working_directory
#   .
#   ├── GRCh38              ----> to save the prepared GRCh38 resources 
#   ├── scripts             ----> to save the bash scripts  
#   ├── mCA                 ----> to save the mChA calling result files
#   |── raw_data            ----> to store the input vcf/bcf raw data 
#   |── logs                ----> to save logs for checkup
#   ├── mis                 ----> to save sample list, metadata and miscellaneous           
#
#
#   Input: Each individual mCA phenotype table .tsv. Make sure the gender of the samples is male only. 
#          It was generated from step1_LoY_mCA_calling: mCA_calls.$pfx.tsv
#
############################################################################################################

# Loss of Y mCA Filter: 
# keep1: lod_lrr_baf >= 5 && rev_col <= 2.9 || (bdev <= 0.16 && rev_col <=2.5)
# keep2: length >= 2000

# Keep: ChrX && Type = Loss for LoY mCA check

############ mCA table HEADER: ##############

# 1:sample_id   2:computed_gender       3:chrom 4:beg_GRCh38    5:end_GRCh38    6:length        7:p_arm 8:q_arm 9:n_sites       10:n_hets
# 11:n50_hets   12:bdev 13:bdev_se      14:rel_cov      15:rel_cov_se   16:lod_lrr_baf  17:lod_baf_phase        18:n_flips      19:baf_conc     20:lod_baf_conc 21:type 22:cf

echo -e "script:step2_LoY filter started at $(date)\n"

# specify your wd
wd="your_proj_working_directory"
echo -e "Working Directory: $wd\n"

mCA_tsv="$wd/mCA/mCA_results/mCA_calls"

# combine all individual mCA called table into a single .tsv for further steps
# save the result in the upper directory
awk '(NR == 1) || (FNR > 1)' ${mCA_tsv}/mCA_calls.*.tsv > $wd/mCA/mCA_results/all.raw.LoY.mCA.tsv 

# if the old filtered LoY mCA list exists, delete it first.
if [ -a $wd/mCA/mCA_results/all.filtered.LoY.mCA.tsv ]
then

echo -e "Old all.filtered.LoY.mCA.tsv exists, delete it first.\n"
rm $wd/mCA/mCA_results/all.filtered.LoY.mCA.tsv

fi

# filtering #
echo "Only Keep lod_baf_phase >= 5 && rev_cov <= 2.9 || (bdev <= 0.16 && rev_cov <=2.5)"
echo -e "Only keep chrX and type = Loss\n"

# lod_lrr_baf >= 5 && rev_col <= 2.9 || (bdev <= 0.16 && rev_col <=2.5) 
# length >= 2000
awk '(($17 >= 5) || (NR==1)) {print $0}'  $wd/mCA/mCA_results/all.raw.LoY.mCA.tsv  \
|  awk '((NR==1) || ($14 <= 2.9) || ($12 <= 0.16 && $14 <= 2.5)) {print $0}'| \
awk '(($6 >= 2000) || (NR==1)) {print $0}'| \
awk '(NR==1)) || (($21== "Loss") && ($3 =="chrX")) {print $0}'
 > $wd/mCA/mCA_results/all.filtered.LoY.mCA.tsv


echo "Done! step2_LoY_filter."
echo "Ending at $(date)"
