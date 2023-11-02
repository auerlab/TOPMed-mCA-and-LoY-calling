

############################################################################################################
#
#   Objective: Down-sample method aimed to balance the variant markers density across difference race or ancestries 
#
#              This script serves an optional and general down-sampling method for any study. 
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
#   Input: GenomeWideHetCounts table for downsample information of TopMed samples  
#
############################################################################################################

setwd("/Your_working_dir/")

# set up the gender and race variables
genders<-c("female", "male")
pops<-c("asians", "black", "hispanic", "white")

#sex="female" # you can specify the gender and race of your interest
#anc="asians"

# make adjusted markers across races 
# GenomeWideHetCounts should at least have two column: Sample ID & adjusted het counts
for(anc in pops){
for(sex in genders){
file=paste0("./GenomeWideHetCounts/new.",sex,".het.count.",anc,".txt")
ofile=paste0(sex,".het.count.",anc,".list")

da<-read.csv(file, header=T,sep="\t")
write.table(cbind(da$NWDID,da$adjusted_tot_hets), file=ofile, col.names=F, row.names=F, quote=F)
}}
