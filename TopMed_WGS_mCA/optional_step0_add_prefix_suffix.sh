#!/bin/bash

############################################################################################################
#
#   Objective: Integrated TopMed WGS mCA calling pipeline from raw vcf data to mCA instances using MoChA.
#
#              This script is optional. It was designed to add the suffix and prefix for each TopMed sample, and
#              eventually make a VCF file with standard name to proceed. 
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
#   09-22-2023        Mr.Xiaolong Ma       Bioinformatics Analyst      v1
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
#   Input: A VCF sample list .txt file, each row should be a single sample.
#   
#   Like below:
#   			NWD123654
#   			NWD456987
#   			NWD159753
#   			......
#   			......
#
############################################################################################################
echo -e "script:optional_step0_add_prefix_suffix started at $(date)/n"

# specify your wd
wd="your_proj_working_directory"
echo -e "Working Directory: $wd\n"

# Make sure the directories are ready
mkdir -p $wd/{logs,raw_data,mis,scripts,GRCh38,mCA}

# Make sure the TopMed sample list exists, it should be one row one sample ID
sample_dir="$wd/mis"

# Add the prefix and suffix for each TopMed sample 
# Such as file: combined.NWD100014-ad.vcf.xz
# pfx="combined."
# sfx="-ad.vcf.xz" 
# ID="NWD100014"

pfx="combined."
sfx="-ad.vcf.xz"

while read ID
do
echo "${pfx}${ID}${sfx}"

done < ${sample_dir}/TopMed_sample_list >> ${sample_dir}/raw_VCF_list
echo -e "Raw VCF files were listed out in: ${sample_dir}/raw_VCF_list\n"

echo -e "script:optional_step0_add_prefix_suffix ended at $(date)/n"