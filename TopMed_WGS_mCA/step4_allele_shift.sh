#!/bin/bash

#SBATCH --job-name=AS
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=24:00:00
#SBATCH --account=your_acc
#SBATCH --array=1-5
#SBATCH --partition=normal
#SBATCH --output=../logs/%u.%x.%j.out
#SBATCH --error=../logs/%u.%x.%j.err

############################################################################################################
#
#   Objective: Integrated TopMed WGS mCA calling pipeline from raw vcf data to mCA instances using MoChA.
#
#              This script was designed to make allele shift analysis for MoChA output mCA VCF.
#              This script focused on gene ATM & MPL shift. 
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
#   09-25-2023        Mr.Xiaolong Ma       Bioinformatics Analyst      v1
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
#   Input: 1: The compiled mCA phenotype .tsv (the all.mCA.tsv), which was generated from optional_step2_mCA_compile.sh. 
#          2: Import allelic shift information from the MoChA output VCF into a single VCF file, which means all the output
#             MoChA VCFs have to be combined together.
#
############################################################################################################
module load bcftools/1.16
module load bcftools-mocha/1.15

# set up the threads, it is equal to #SBATCH --ntasks=20
threads=20

echo -e "script:step4_allele_shift started at $(date)\n"
echo -e "Job name: ${SLURM_JOB_NAME}, Job ID: ${SLURM_JOB_ID}\n"

# specify your wd
wd="your_proj_working_directory"
mCA_vcf="$wd/mCA/mCA_results/mCA_vcf"
echo -e "Working Directory: $wd\n"

mkdir -p $wd/mCA/AS

# Keep the chr1 & chr11 (where the MPL & ATM located) CN-LOH phenotypes 

# Header Information for MoChA mCA call .tsv
#-------------------------------------------------------------------------
# 1:sample_id	2:computed_gender	3:chrom	4:beg_GRCh38	5:end_GRCh38	6:length	7:p_arm	8:q_arm	9:n_sites	10:n_hets	
# 11:n50_hets	12:bdev	13:bdev_se	14:rel_cov	15:rel_cov_se	16:lod_lrr_baf	17:lod_baf_phase	18:n_flips	
# 19:baf_conc	20:lod_baf_conc	21:type	22:cf
#-------------------------------------------------------------------------

# keep chr1 & p_arm=T & CN-LOH --> MPL: chromosome 1 (1p)
# keep chr11 & q_arm=T & CN-LOH  ---> ATM:chromosome 11 (11q)

as="$wd/mCA/AS"
all_mCA="$wd//mCA/mCA_results/all.mCA.tsv"

# For MPL sample ID 
awk '{if($3 == "chr1" && $7 == "T" && $21 == "CN-LOH") print $1}' ${all_mCA} > $as/MPL.sample.id

# For ATM sample ID 
awk '{if($3 == "chr11" && $8 == "T" && $21 == "CN-LOH") print $1}' ${all_mCA} > $as/ATM.sample.id

# make the index file for each mCA VCF for both MPL & ATM 
while read id
do 
	bcftools index --force ${mCA_vcf}/mCA_vcf.$id.bcf
	echo "${mCA_vcf}/mCA_vcf.$id.bcf" > ${mCA_vcf}/MPL.vcf.list 
done < $as/MPL.sample.id

while read id
do 
	bcftools index --force ${mCA_vcf}/mCA_vcf.$id.bcf
	echo "${mCA_vcf}/mCA_vcf.$id.bcf" > ${mCA_vcf}/ATM.vcf.list 
done < $as/ATM.sample.id

# Make the file list for VCF combination 
MPL="${mCA_vcf}/MPL.vcf.list"
ATM="${mCA_vcf}/ATM.vcf.list"

##### merge all MPL vcf together ####
echo -e "Combining MPL(chr1.p_arm.CN-LOH) bcf files into a single one\n\n"
bcftools merge -l ${MPL} -Ob -o $wd/mCA/AS/AS.MPL_merged.bcf
bcftools index --force $wd/mCA/AS/AS.MPL_merged.bcf 
echo -e "All MPL mCA VCF were combined.\n"

##### merge all ATM vcf together ####
echo -e "Combining ATM(chr11.q_arm.CN-LOH) bcf files into a single one\n\n"
bcftools merge -l ${ATM} -Ob -o $wd/mCA/AS/AS.ATM_merged.bcf
bcftools index --force $wd/mCA/AS/AS.ATM_merged.bcf
echo -e "All ATM mCA VCF were combined.\n"

echo -e "Allele Shift Analysis Begins...\n"

# MPL AS
bcftools +extendFMT \
  --no-version -Ou \
  --format AS \
  --phase \
  --dist 500000 \
  $wd/mCA/AS/AS.MPL_merged.bcf | \
bcftools +mochatools \
  --no-version \
  --output-type b \
  -- --summary AS \
  --test AS \
  --drop-genotypes | \
  tee $wd/mCA/AS/Result.AS.MPL.bcf | \
  bcftools index --force --output $wd/mCA/AS/Result.AS.MPL.bcf.csi

echo -e "77 ID AS analysis is done.\n"

# ATM AS
bcftools +extendFMT \
  --no-version -Ou \
  --format AS \
  --phase \
  --dist 500000 \
  $wd/mCA/AS/AS.ATM_merged.bcf | \
bcftools +mochatools \
  --no-version \
  --output-type b \
  -- --summary AS \
  --test AS \
  --drop-genotypes | \
  tee  $wd/mCA/AS/Result.AS.ATM.bcf  | \
  bcftools index --force --output $wd/mCA/AS/Result.AS.ATM.bcf.csi

echo -e "MPL & ATM AS is done.\n"




