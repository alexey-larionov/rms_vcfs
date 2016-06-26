#!/bin/bash

# a01_read_config.sh
# Parse config file for merging rms dbGAP vcfs
# Alexey Larionov, 10Jun2016

# Function for reading parameters
function get_parameter()
{
	local parameter="${1}"
  local line
	line=$(awk -v p="${parameter}" 'BEGIN { FS=":" } $1 == p {print $2}' "${job_file}") 
	echo ${line} # return value
}

# === Data location and analysis settings === # 

data_server=$(get_parameter "data_server") # e.g. admin@mgqnap.medschl.cam.ac.uk
project_location=$(get_parameter "project_location") # e.g. /share/mae

project=$(get_parameter "project") # e.g. rms_dbGAP 
input_folder=$(get_parameter "input_vcfs_folder") # e.g. s04_filtered_vcfs
suffix=$(get_parameter "input_vcfs_suffix") # e.g. filt
output_folder=$(get_parameter "output_folder") # e.g. s05_merged_vcf_v1

remove_project_folder=$(get_parameter "Remove_project_folder_from_HPC_scratch_after_run") # e.g. no

# ============= mgqnap settings ============= #

mgqnap_user=$(get_parameter "mgqnap_user") # e.g. mae
mgqnap_group=$(get_parameter "mgqnap_group") # e.g. mtgroup

# =============== HPC settings ============== #

working_folder=$(get_parameter "working_folder") # e.g. /scratch/medgen/users/mae 
 
account_to_use=$(get_parameter "Account_to_use_on_HPC") # e.g. TISCHKOWITZ-SL2 
time_to_request=$(get_parameter "Max_time_to_request_(hrs.min.sec)") # e.g. 12.00.00 
time_to_request=${time_to_request//./:} # substitute dots to colons  

# ============ Standard settings ============ #

# ----------- Tools ---------- #

tools_folder=$(get_parameter "tools_folder") # e.g. /scratch/medgen/tools

java7=$(get_parameter "java7") # e.g. java/jre1.7.0_76/bin/java
java7="${tools_folder}/${java7}"

gatk=$(get_parameter "gatk") # e.g. gatk/gatk-3.4-46/GenomeAnalysisTK.jar 
gatk="${tools_folder}/${gatk}" 
 
bcftools=$(get_parameter "bcftools") # e.g. bcftools/bcftools-1.2/bin/bcftools 
bcftools="${tools_folder}/${bcftools}" 

plot_vcfstats=$(get_parameter "plot_vcfstats") # e.g. bcftools/bcftools-1.2/bin/plot-vcfstats 
plot_vcfstats="${tools_folder}/${plot_vcfstats}" 

python_bin=$(get_parameter "python_bin") # e.g. python/python_2.7.10/bin/ 
python_bin="${tools_folder}/${python_bin}" 
PATH="${python_bin}:${PATH}" # for updated version of matplotlib library for plot-vcfstats

r_folder=$(get_parameter "r_folder") # e.g. r/R-3.2.0/bin
r_folder="${tools_folder}/${r_folder}"
PATH="${r_folder}:${PATH}" # for GATK plots

r_bin_folder=$(get_parameter "r_bin_folder") # e.g. r/R-3.2.2/bin/
r_bin_folder="${tools_folder}/${r_bin_folder}" # for Rmarkdown reports

r_lib_folder=$(get_parameter "r_lib_folder") # e.g. r/R-3.2.2/lib64/R/library
r_lib_folder="${tools_folder}/${r_lib_folder}" # for Rmarkdown reports

# ----------- Resources ---------- #

resources_folder=$(get_parameter "resources_folder") # e.g. /scratch/medgen/resources

decompressed_bundle_folder=$(get_parameter "decompressed_bundle_folder") # e.g. gatk_bundle/hg19/decompressed
decompressed_bundle_folder="${resources_folder}/${decompressed_bundle_folder}"

ref_genome=$(get_parameter "ref_genome") # e.g. ucsc.hg19.fasta
ref_genome="${decompressed_bundle_folder}/${ref_genome}"

hapmap=$(get_parameter "hapmap") # e.g. hapmap_3.3.hg19.sites.vcf
hapmap="${decompressed_bundle_folder}/${hapmap}"

omni=$(get_parameter "omni") # e.g. 1000G_omni2.5.hg19.sites.vcf
omni="${decompressed_bundle_folder}/${omni}"

phase1_1k_hc=$(get_parameter "phase1_1k_hc") # e.g. 1000G_phase1.snps.high_confidence.hg19.sites.vcf
phase1_1k_hc="${decompressed_bundle_folder}/${phase1_1k_hc}"

dbsnp_138=$(get_parameter "dbsnp_138") # e.g. dbsnp_138.hg19.vcf
dbsnp_138="${decompressed_bundle_folder}/${dbsnp_138}"

dbsnp_138_sites129=$(get_parameter "dbsnp_138_sites129") # e.g. dbsnp_138.hg19.excluding_sites_after_129.vcf
dbsnp_138_sites129="${decompressed_bundle_folder}/${dbsnp_138_sites129}"

mills=$(get_parameter "mills") # e.g. Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
mills="${decompressed_bundle_folder}/${mills}"

broad_exomes_intervals=$(get_parameter "broad_exomes_intervals") # e.g. broad_human_exome_hg19.interval_list
broad_exomes_intervals="${decompressed_bundle_folder}/${broad_exomes_intervals}"

ph3_1k_folder=$(get_parameter "ph3_1k_folder") # e.g. phase3_1k_release20130502/vcfs
ph3_1k_folder="${resources_folder}/${ph3_1k_folder}"

ph3_1k_vcf=$(get_parameter "ph3_1k_vcf") # e.g. ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.hg19.sites.vcf
ph3_1k_vcf="${ph3_1k_folder}/${ph3_1k_vcf}"

fa_mask_folder=$(get_parameter "fa_mask_folder") # e.g. phase3_1k_release20130502/variants_with_frequent_ALT_allele_b37
fa_mask_folder="${resources_folder}/${fa_mask_folder}"

fa_mask_90=$(get_parameter "fa_mask_90") # e.g. FAA_mask_90.vcf
fa_mask_90="${fa_mask_folder}/${fa_mask_90}"

fa_mask_95=$(get_parameter "fa_mask_95") # e.g. FAA_mask_95.vcf
fa_mask_95="${fa_mask_folder}/${fa_mask_95}"

fa_mask_99=$(get_parameter "fa_mask_99") # e.g. FAA_mask_99.vcf
fa_mask_99="${fa_mask_folder}/${fa_mask_99}"

fa_mask_100=$(get_parameter "fa_mask_100") # e.g. FAA_mask_100.vcf
fa_mask_100="${fa_mask_folder}/${fa_mask_100}"

# ----------- Working folders ---------- #
 
project_folder="${working_folder}/${project}" # e.g. rms_dbGAP 

output_folder="${project_folder}/${output_folder}" # e.g. s05_merged_vcf_v1

logs_folder=$(get_parameter "logs_folder") # e.g. logs
logs_folder="${output_folder}/${logs_folder}" 

vqsr_folder=$(get_parameter "vqsr_folder") # e.g. vqsr
vqsr_folder="${output_folder}/${vqsr_folder}" 

vcfstats_folder=$(get_parameter "vcfstats_folder") # e.g. vcfstats
vcfstats_folder="${output_folder}/${vcfstats_folder}" 
all_vcfstats_folder="${vcfstats_folder}/all" 
cln_vcfstats_folder="${vcfstats_folder}/cln" 

histograms_folder=$(get_parameter "histograms_folder") # e.g. histograms
histograms_folder="${output_folder}/${histograms_folder}" 
