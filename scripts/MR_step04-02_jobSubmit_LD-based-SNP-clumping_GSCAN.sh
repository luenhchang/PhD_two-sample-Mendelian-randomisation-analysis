#!/bin/bash
## file name: MR_step04-02_jobSubmit_LD-based-SNP-clumping_GSCAN.sh
## old file name: 
## modified from: 
## date created: 20180719
## purpose: select independent SNPs
## Run dependency: 
## How to run this script: . ${locScripts}/MR_step04-02_jobSubmit_LD-based-SNP-clumping_GSCAN.sh > ${locHistory}/jobSubmitted_20190408_LD-based-clumping_GSCAN-GWAS_sample-excluded-ICC-23andMe-QIMR-UKB

## Time 	Change
##--------------------------------------------------------------------------------------------------------------
## 20190408	qsub jobs 7296750,7296751,7296754,7296755,7296758,7296759,7296763-7296766 (10 jobs)
## 20190403 qsub jobs 7275558-7275567
## 20190305 qsub jobs 7069251-7069260 (set p1 to 5e-8 and p2 to 1e-6)
## 20181119	qsub jobs 6535209-6535219
## 20180729 qsub jobs 6025935-6025944 (10 jobs)
## reclump smoking-related GWAS with less stringent --clump-p1
### --clump-p1 1e-6 
### --clump-p2 1e-6 
## 20180719	qsub jobs 5954723-5954732 (10 jobs)
## --clump-p1 5e-8 [index variant p-value threshold] This is the p value threshold for choosing the best SNP per haplotype
## --clump-p2 1e-6 [clumped variant p-value threshold] Once the best SNP is chosen, other SNPs within the same haplotype with p value < p2 and r2 0.01 are chosen (i.e. clumped)
##--------------------------------------------------------------------------------------------------------------

# Type 	File
#--------------------------------------------------------------------------------------------------------------
# Input	${locInput}/LD-SNP-clumping-criteria.txt
# Input ${locQC3}/ai_noICC.ambiguSNPRemoved
# Input ${locQC3}/cpd_noICC.ambiguSNPRemoved
# Input ${locQC3}/dpw_noICC.ambiguSNPRemoved
# Input ${locQC3}/sc_noICC.ambiguSNPRemoved
# Input ${locQC3}/si_noICC.ambiguSNPRemoved
# Outpu ${locOutput}/*_noICC_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped (10 files)
#---------------------------------------------------------------------------------------------------------------
## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/MR_ICC_GSCAN_201806";
jobScriptFilePath="${locScripts}/MR_step04-01-01_jobScript_LD-based-SNP-clumping.sh";
locHistory="${homeDir}/history";
locGSCAN="${homeDir}/LabData/Lab_NickM/lunC/GSCAN";
locICC="${homeDir}/LabData/Lab_NickM/lunC/international_cannabis_consortium";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locMR="${workingDir}/MR_ICC_GSCAN_201806/data";
locGSCAN="$locMR/noICC_results";
locQC3=$locGSCAN/QC3_remove_ambiguousSNPs_indel;

locLDBasedClumping=$locGSCAN/LD-based_SNP_clumping;
locInput=$locLDBasedClumping/input;
locOutput=$locLDBasedClumping/output;
locArchive=$locOutput/archive;
pbs_output_dir=$locLDBasedClumping/pbs_output;
pbs_archive=$pbs_output_dir/archive;
mkdir -p $locLDBasedClumping $locInput $locOutput $pbs_output_dir $locArchive $pbs_archive;

bfile="/mnt/lustre/reference/data/UKBB_500k/versions/lab_stuartma/LD_reference/LD_ref_201803_rsOnly";
email_address="luenhchang@gmail.com";

#------------------------------------------------------------------------------------------------------------#
#------------------------Archive existing files before repeating same analyses-------------------------------#
#------------------------------------------------------------------------------------------------------------#
# Archive every files in $locOutput except for the archive folder for backup purposes
# cd $locOutput
# mv !(archive) $locArchive/;

# Archive every files in $pbs_output_dir except for the archive folder for backup purposes
# cd $pbs_output_dir
# mv !(archive) $pbs_archive/;

#------------------------------------------------------------------------------------------------------------#
#------------------------Perform SNP clumping for discovery sample GSCAN ------------------------------------#
#------------------------------------------------------------------------------------------------------------#
# Store file paths of input files as a file
realpath $locQC3/*_noICC.ambiguousSNPRemoved > $locQC3/filePath # 5 file paths

# Loop thru each line of the comma-separated GWAS information file, skipping 1st row (header)
## Number of iteration: 5 (files) * 2 (sets of clumping criteria)
IFS=$'\n';
count=0;
for filePath in `cat $locQC3/filePath`;do
shortFileName=`basename $filePath | cut -d"." -f1`; 
clumpFileName="SNP_PValue_${shortFileName}";
# Get SNP and P columns from each input file
awk 'BEGIN {print "SNP P"} (NR !=1) {print $3, $7}' $filePath > $locInput/${clumpFileName};
for line in `tail -n+2 ${locInput}/LD-SNP-clumping-criteria.txt`; do
count=$((${count}+1)); 
echo "========================================== iteration ${count} ============================================"; 
outputDir=$locOutput;
clumping_r2_threshold=$(echo $line | cut -d" " -f1);
LDWindow=$(echo $line | cut -d" " -f2);
p_threshold_1=$(echo $line | cut -d" " -f3);
p_threshold_2=$(echo $line | cut -d" " -f4);
outputFileName=${shortFileName}_LDWindow-kb-${LDWindow}_R2-${clumping_r2_threshold}_p1-${p_threshold_1}_p2-${p_threshold_2};
jobName="LDBasedSNPclumping_${shortFileName}";
echo "shortFileName=$shortFileName"; 
echo "clumpFileName=$clumpFileName"; 
echo "clumping_r2_threshold=$clumping_r2_threshold";
echo "LDWindow=$LDWindow";
echo "p_threshold_1=$p_threshold_1";
echo "p_threshold_2=$p_threshold_2";
echo "qsub -N ${jobName} -m bea -M $email_address -v v_bfilePath=${bfile},v_clumpFilePath=${locInput}/${clumpFileName},v_clumping_r2_threshold=${clumping_r2_threshold},v_clumping_LD_window=${LDWindow},v_p_threshold_1=${p_threshold_1},v_p_threshold_2=${p_threshold_2},v_outputDir=${outputDir},v_outputResultFileName=${outputFileName} -e ${pbs_output_dir}/${outputFileName}.pbs.err -o ${pbs_output_dir}/${outputFileName}.pbs.out -l ncpus=1,walltime=01:00:00,mem=2500mb ${jobScriptFilePath}";
qsub -N ${jobName} -m bea -M $email_address -v v_bfilePath=${bfile},v_clumpFilePath=${locInput}/${clumpFileName},v_clumping_r2_threshold=${clumping_r2_threshold},v_clumping_LD_window=${LDWindow},v_p_threshold_1=${p_threshold_1},v_p_threshold_2=${p_threshold_2},v_outputDir=${outputDir},v_outputResultFileName=${outputFileName} -e ${pbs_output_dir}/${outputFileName}.pbs.err -o ${pbs_output_dir}/${outputFileName}.pbs.out -l ncpus=1,walltime=01:00:00,mem=2500mb ${jobScriptFilePath};
done	
done

#cp -n ${locScripts}/MR_step04-02_jobSubmit_LD-based-SNP-clumping_GSCAN.sh ${locScripts}/MR_step04-03_jobSubmit_LD-based-SNP-clumping_UKB.sh

##---------------------------------This is the end of this file-------------------------------------------------##