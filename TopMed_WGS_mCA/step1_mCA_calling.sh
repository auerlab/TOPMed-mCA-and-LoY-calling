#!/bin/bash

#SBATCH --job-name=mCA
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
#              This script is required. It was designed to call mCA phenotypes for input VCF using MoChA.
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
#   Input: A VCF list .txt file, each row should be a single VCF. Alternatively, you can list out the sample ID in a .txt file, 
#          and manually add up the suffix and prefix as their filenames respectively.
#
#          This script assume you have multiple samples for calling; therefore, it has been applied to SLURM array setting for parellel proocessing.
#          
#          If your VCF does not have header, this script will help you concatenate the header and main body together as an intact VCF.
#          Please make sure the general header.txt file is saved under your $wd/mis
#
#          First of all, you have to specify your account and array number in the top SLURM setting accordingly. 
#   
#   VCF list is like below:
#   						proj1_NWD123654.vcf
#   						proj1_NWD456987.vcf
#   						proj1_NWD159753.vcf
#   						......
#   						......
#
############################################################################################################
module load bcftools/1.16
module load bcftools-mocha/1.15

# set up the threads, it is equal to #SBATCH --ntasks=20
threads=20

echo -e "script:step1_mCA_calling started at $(date)\n"
echo -e "Job name: ${SLURM_JOB_NAME}, Job ID: ${SLURM_JOB_ID}\n"

# specify your wd
wd="your_proj_working_directory"
echo -e "Working Directory: $wd\n"

# Make sure the directories are ready
mkdir -p $wd/{logs,raw_data,mis,scripts,GRCh38,mCA}
mkdir -p $wd/mCA/{gc,mCA_results} # to save mCA results
mkdir -p $wd/mCA/mCA_results/{mCA_calls,mCA_stats,mCA_vcf} # to save mCA results

# Make sure the VCF file list exists, it should be one row one VCF file
sample_dir="$wd/mis"

# ## #SBATCH --array=1-5, supposing you have 5 VCFs in such a -----> vcf.list.txt
vcf=$(head -n ${SLURM_ARRAY_TASK_ID} ${sample_dir}/vcf.list.txt | tail -n +${SLURM_ARRAY_TASK_ID})   

# sample ID extract
id=$(basename $vcf ".vcf") # vcf=proj1_NWD123654.vcf --> id=proj1_NWD123654
echo -e "$vcf is processing, sample_id=$id\n"

# Make sure your reference fasta exists. Please refer to https://github.com/freeseek/mocha for GRCh38 resources
# Supposing you are using hg38 for such as an analysis
ref="$wd/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
echo -e "GRCh38 reference: ${ref} is applied.\n"

# If your VCF does not have header, the following steps will help you concatenate the header before calling.

############ The following section is optional ######################

# make a specific header for specific sample 
cp ${wd}/mis/head.txt ${wd}/raw_data/head.${id}.txt
echo -e "Adding up header back to VCF and cut off unnecessary format info.\n"

# label the sample ID
sed -i "s/SAMPLE/$id/g" ${wd}/raw_data/head.${id}.txt  

# concatenate the header with vcf main body 
cat ${wd}/raw_data/head.${id}.txt ${wd}/raw_data/$vcf > ${wd}/raw_data/head.$vcf

# chop off the unnecessary format and only keep the ones MoChA required
sed -Ei "s/,[[:digit:]]*:[[:digit:]]*$//g" ${wd}/raw_data/head.$vcf
sed -Ei "s/:DP//g" ${wd}/raw_data/head.$vcf

#####################################################################

# Keep the variants with AD >=5

bcftools view -i 'MIN(FMT/AD)>=5' --threads ${threads} -Ou ${wd}/raw_data/head.$vcf | bcftools norm --no-version -Ob -o ${wd}/raw_data/${id}.bcf -d none -f $ref
bcftools index --threads ${threads} -f ${wd}/raw_data/${id}.bcf
echo -e "Keep the variants with AD>=5\n"

# remove intermediate files 
rm ${wd}/raw_data/head.${id}.txt
rm ${wd}/raw_data/head.$vcf

echo -e "$id sample: GC annotation is processing\n"
gcbcf="${wd}/mCA/gc/gc.${id}.bcf"

bcftools +mochatools ${wd}/raw_data/${id}.bcf --threads ${threads} -Ob -o $gcbcf -- -t GC --fasta-ref $ref
bcftools index --threads ${threads} -f $gcbcf
echo -e "GC annotation is done. Output saved under: ${wd}/mCA/gc/gc.${id}.bcf\n"

echo -e "$id sample: mCA Call is processing\n"

#output setup
outdir1="${wd}/mCA/mCA_results/mCA_calls"
outdir2="${wd}/mCA/mCA_results/mCA_stats"
outdir3="${wd}/mCA/mCA_results/mCA_vcf"

otsv1=${outdir1}/mCA_calls.$id.tsv
otsv2=${outdir2}/mCA_stats.$id.tsv
ovcf=${outdir3}/mCA_vcf.$id.bcf

bcftools +mocha --threads ${threads} -g GRCh38 $gcbcf -z $otsv2 -c $otsv1 --min-dist 1000 --LRR-weight 0.0 --bdev-LRR-BAF 6.0 -o $ovcf -Ob


echo "Done! step1_mCA_calling for sample:${id}"
echo "Ending at $(date)"











