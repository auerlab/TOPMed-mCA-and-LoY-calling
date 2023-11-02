#!/bin/bash 

############################################################################################################
#
#   Objective: Integrated TopMed WGS mCA calling pipeline from raw vcf data to mCA instances using MoChA.
#
#              This script is optional. It was designed to compile all mCA phenotype table into a single .tsv
#              The output .tsv will be easiler to be filtered by means of step3 R script.
#              
#   NOTE: Please refer to https://github.com/freeseek/mocha for GRCh38 resources.
#         In this pipeline, we assume they are ready to be applied.
#
#   	  This pipeline is contingent to be tweaked according to your research requirements.
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
#   Contact			  Institute
#   xima@mcw.edu      Medical College of Wisconsin
#
############################################################################################################
#
#   Directory Hierarchy
#
#   your_proj_working_directory
#   .
#   ├── GRCh38              ----> to save the prepared GRCh38 resources 
#   ├── scripts 			----> to save the bash scripts  
#   ├── mCA             	----> to save the mChA calling result files
#   |── raw_data            ----> to store the input vcf/bcf raw data 
#	|── logs				----> to save logs for checkup
#   ├── mis    				----> to save sample list, metadata and miscellaneous		
#
#
#   Input: The absolute path of the directory where the mCA phenotype tables saved.
#
############################################################################################################

echo -e "script:optional_step2_mCA_compile started at $(date)\n"

wd="your_proj_working_directory"
mCA_dir="$wd/mCA/mCA_results/mCA_calls"
output="$wd/mCA/mCA_results"

# Project Name 
proj_pfx="proj1"

# combine all the mCA phenotypes from each individual mCA table to a single .tsv file
# NR == 1 --> pick up the header at a time
# FNR > 1 --> pick up every mCA phenotype from each mCA_calls .tsv file 

awk '(NR == 1) || (FNR > 1)' ${mCA_dir}/mCA_calls.*.tsv > ${output}/${proj_pfx}.all.mCA.tsv
echo -e "All mCA phenotypes are compiled into a single .tsv and saved under: ${output}/${proj_pfx}.all.mCA.tsv\n"


echo "Done! optional_step2_mCA_compile.sh"
echo "Ending at $(date)"