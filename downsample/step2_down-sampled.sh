#!/bin/bash


############################################################################################################
#
#   Objective: Down-sample method aimed to balance the variant markers density across difference race or ancestries 
#
#              This script serves an optional and general down-sampling method for any study. 
#              
#              
#   NOTE: Please refer to https://github.com/auerlab/biolibc for biolibc package installation.
#         Please refer to https://github.com/auerlab/biolibc-tools for biolibc-tools package installation.
#         Down-sample is one of the tools offered by biolibc-tools.
#
#
#         bcftools(1.11 or later), biolibc(0.2.4 or later), biolibc-tools
#
############################################################################################################
#   
#   Revision History: 
#
#   Date              Author               Title                       Modification
#   10-02-2023        Mr.Xiaolong Ma       Bioinformatics Analyst      v1
#   
#   Contact                       Institute
#   xima@mcw.edu      Medical College of Wisconsin
#
#
#   Input: A list generated from step1_select_sample, which has VCF sample ID (col1) and number of 
#          adjusted markers (col2) for down-sampling.
# 
#   vcf.list:
#             NWD112233 765432
#			  NWD102030 123456
#			  NWD998877 908070
#             ......
#             ......               
#
############################################################################################################


module load bcftools/1.11
module load biolibc/0.2.4
module load biolibc

# specify you working directory
wd="Your_path_here"

# set up the downsample list and output directory
downsample_list="$wd/male.het.count.black.list" # Take black male as an example 

echo -e "Down-sampling each VCF from ${downsample_list}...\n"

while read id hets # id is NWD ID; hets is adjusted number of markers across race or ancestries
do

bcf="${wd}/raw_vcf/${id}.bcf"
bcfout="${wd}/down-sample/down-sampled.${id}.bcf"

echo -e "$id is processing...\n"
bcftools view ${bcf} -Ov | vcf-downsample ${hets} | bcftools view -Ob --output ${bcfout}
bcftools index -f ${bcfout}
echo -e "$id is down-sampled.\n"

done < ${downsample_list}

echo -e "step2_down-sampled.sh is done.\n"





