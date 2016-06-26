#!/bin/bash

# s01_rms_merge_vcfs.sh
# Merging rms dbGAP vcfs
# Started: Alexey Larionov, 10Jun2016
# Last updated: AL, 23Jun2016

# --- sbatch section --- #

#SBATCH -J merge_rms_vcfs
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH -p sandybridge
##SBATCH --qos=INTR

# Set initial working folder
cd "${SLURM_SUBMIT_DIR}"

# --- Modules section --- #

# Standard cluster's modules
. /etc/profile.d/modules.sh
module purge
module load default-impi

# Additional modules for knitr-rmarkdown (used for histograms)
module load gcc/5.2.0
module load boost/1.50.0
module load texlive/2015
module load pandoc/1.15.2.1

# --- Set environment --- #

# Stop at any error
set -e

# Command line arguments
job_file="${1}"
scripts_folder="${2}"

# Read and report setings
source "${scripts_folder}/a01_read_config.sh"

echo "Started s01_rms_merge_vcfs: $(date +%d%b%Y_%H:%M:%S)"
echo ""
echo "Job name: ${SLURM_JOB_NAME}"
echo "Allocated node: $(hostname)"
echo ""
echo "Initial working folder:"
echo "${SLURM_SUBMIT_DIR}"
echo ""
echo "====================== Settings ======================"
echo ""

source "${scripts_folder}/a02_report_settings.sh"

echo "=================== Pipeline steps ==================="
echo ""

# Make folders
tmp_folder="${output_folder}/tmp"
mkdir -p "${tmp_folder}"
mkdir -p "${vqsr_folder}"
mkdir -p "${all_vcfstats_folder}"
mkdir -p "${cln_vcfstats_folder}"
mkdir -p "${histograms_folder}"

dataset="${project}"

# Go to working folder (= the output folder)
init_dir="$(pwd)"
cd "${output_folder}"

# --- Copy source vcfs to cluster --- #

# Progress report
echo "Started copying source data: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# Copy and read list of samples
rsync -ae "ssh -x" "${data_server}:${project_location}/${project}/${input_folder}/samples.txt" "${tmp_folder}/"
samples=$(<"${tmp_folder}/samples.txt")

# For each sample
for sample in ${samples}
do

  # Copy vcf
  if ! rsync -ae "ssh -x" "${data_server}:${project_location}/${project}/${input_folder}/${sample}.${suffix}.vcf" "${tmp_folder}/"
  then
    echo ""
    echo "Failed getting source vcf from NAS"
    echo "Sample: ${sample}"
    echo "Script terminated"
    echo ""
    exit
  fi

  # Copy idx
  if ! rsync -ae "ssh -x" "${data_server}:${project_location}/${project}/${input_folder}/${sample}.${suffix}.vcf.idx" "${tmp_folder}/"
  then
    echo ""
    echo "Failed getting source idx from NAS"
    echo "Sample: ${sample}"
    echo "Script terminated"
    echo ""
    exit
  fi

  # Progress report
  echo "${sample}"

done # next sample
echo ""

# Progress report
echo "Completed copying source data"
echo ""

# --- Merge vcfs --- #

# Progress report
echo "Started merging vcfs: $(date +%d%b%Y_%H:%M:%S)"

# File names
merged_raw_vcf="${tmp_folder}/${dataset}_merged_raw.vcf"
merging_log="${logs_folder}/${dataset}_merging_vcfs.log"

# make list of source vcfs
vcf_files=""
for sample in ${samples}
do
  if [ -z "${vcf_files}" ]; then
    vcf_files="-V ${tmp_folder}/${sample}.${suffix}.vcf"
  else 
    vcf_files="${vcf_files} -V ${tmp_folder}/${sample}.${suffix}.vcf"
  fi
done

# Merge vcf files
"${java7}" -Xmx60g -jar "${gatk}" \
  -T CombineVariants \
  -R "${ref_genome}" \
  -L "${broad_exomes_intervals}" -ip 100 \
  ${vcf_files} \
  -o "${merged_raw_vcf}" \
  -nt 14 &>  "${merging_log}"

# Progress report
echo "Completed merging vcfs"
echo ""

# --- Trim the variants --- #
# Removes variants and alleles that have not been detected in any genotype

# Progress report
echo "Started trimming variants: $(date +%d%b%Y_%H:%M:%S)"

# File names
trim_vcf="${tmp_folder}/${dataset}_merged_trim.vcf"
trim_log="${logs_folder}/${dataset}_trimming_vcfs.log"

"${java7}" -Xmx60g -jar "${gatk}" \
  -T SelectVariants \
  -R "${ref_genome}" \
  -L "${broad_exomes_intervals}" -ip 100 \
  -V "${merged_raw_vcf}" \
  -o "${trim_vcf}" \
  --excludeNonVariants \
  --removeUnusedAlternates \
  -nt 14 &>  "${trim_log}"

# Note: 
# This trimming looks excessive in our pipeline

# Progress report
echo "Completed trimming variants: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Add variants IDs to INFO field --- #
# To trace variants during the later steps

# Progress report
echo "Started adding variants IDs to INFO field"

# File name
trim_id_vcf="${tmp_folder}/${dataset}_trim_id.vcf"

# Compile names for temporary files
tmp1=$(mktemp --tmpdir="${tmp_folder}" "${dataset}_tmp1".XXXXXX)
tmp2=$(mktemp --tmpdir="${tmp_folder}" "${dataset}_tmp2".XXXXXX)
tmp3=$(mktemp --tmpdir="${tmp_folder}" "${dataset}_tmp3".XXXXXX)
tmp4=$(mktemp --tmpdir="${tmp_folder}" "${dataset}_tmp4".XXXXXX)

# Prepare data witout header
grep -v "^#" "${trim_vcf}" > "${tmp1}"
awk '{printf("RawVarID=var%09d\t%s\n", NR, $0)}' "${tmp1}" > "${tmp2}"
awk 'BEGIN {OFS="\t"} ; { $9 = $9";"$1 ; print}' "${tmp2}" > "${tmp3}"
cut -f2- "${tmp3}" > "${tmp4}"

# Prepare header
grep "^##" "${trim_vcf}" > "${trim_id_vcf}"
echo '##INFO=<ID=RawVarID,Number=1,Type=String,Description="Raw Variant ID">' >> "${trim_id_vcf}"
grep "^#CHROM" "${trim_vcf}" >> "${trim_id_vcf}"

# Append data to header in the output file
cat "${tmp4}" >> "${trim_id_vcf}"

# Clean-up
rm "${tmp1}" "${tmp2}" "${tmp3}" "${tmp4}"

# Progress report
echo "Completed adding variants IDs to INFO field: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Make mask for multiallelic variants --- #

# Progress report
echo "Started making mask for multiallelic variants"

# File names
trim_id_ma_mask_vcf="${tmp_folder}/${dataset}_trim_id_ma_mask.vcf"
trim_id_ma_mask_log="${logs_folder}/${dataset}_trim_id_ma_mask.log"

# Make mask
"${java7}" -Xmx60g -jar "${gatk}" \
  -T SelectVariants \
  -R "${ref_genome}" \
  -L "${broad_exomes_intervals}" -ip 100 \
  -V "${trim_id_vcf}" \
  -o "${trim_id_ma_mask_vcf}" \
  -restrictAllelesTo MULTIALLELIC \
  -nt 14 &>  "${trim_id_ma_mask_log}"

# Progress report
echo "Completed making mask: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Add flag for multiallelic variants --- #

# Progress report
echo "Started adding flag for multiallelic variants"

# File names
trim_id_ma_vcf="${tmp_folder}/${dataset}_trim_id_ma.vcf"
trim_id_ma_log="${logs_folder}/${dataset}_trim_id_ma.log"

# Add info
"${java7}" -Xmx60g -jar "${gatk}" \
  -T VariantAnnotator \
  -R "${ref_genome}" \
  -L "${broad_exomes_intervals}" -ip 100 \
  -V "${trim_id_vcf}" \
  -comp:MultiAllelic "${trim_id_ma_mask_vcf}" \
  -o "${trim_id_ma_vcf}" \
  -nt 14 &>  "${trim_id_ma_log}"

# Progress report
echo "Completed adding flag for multiallelic variants: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Add flags for variants with frequent alt allele in 1k ph3 (hg19) --- #

# Progress report
echo "Started adding flags for variants with frequent alt allele in 1k ph3 (hg19)"

# File names
trim_id_ma_fa_vcf="${tmp_folder}/${dataset}_trim_id_ma_fa.vcf"
trim_id_ma_fa_log="${logs_folder}/${dataset}_trim_id_ma_fa.log"

# Add info
"${java7}" -Xmx60g -jar "${gatk}" \
  -T VariantAnnotator \
  -R "${ref_genome}" \
  -L "${broad_exomes_intervals}" -ip 100 \
  -V "${trim_id_ma_vcf}" \
  -comp:ALT_frequency_in_1k_90 "${fa_mask_90}" \
  -comp:ALT_frequency_in_1k_95 "${fa_mask_95}" \
  -comp:ALT_frequency_in_1k_99 "${fa_mask_99}" \
  -comp:ALT_frequency_in_1k_100 "${fa_mask_100}" \
  -o "${trim_id_ma_fa_vcf}" \
  -nt 14 &>  "${trim_id_ma_fa_log}"

# Progress report
echo "Completed adding flag: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Add actual AF values from 1k ph3 (hg19) --- #

# Progress report
echo "Add AF values from 1k ph3 (hg19)"

# File names
trim_id_ma_fa_1k_vcf="${tmp_folder}/${dataset}_trim_id_ma_fa_1k.vcf"
trim_id_ma_fa_1k_log="${logs_folder}/${dataset}_trim_id_ma_fa_1k.log"

"${java7}" -Xmx60g -jar "${gatk}" \
  -T VariantAnnotator \
  -R "${ref_genome}" \
  -L "${broad_exomes_intervals}" -ip 100 \
  -V "${trim_id_ma_fa_vcf}" \
  --resource:ph3_1k "${ph3_1k_vcf}" \
  --expression ph3_1k.AF \
  -o "${trim_id_ma_fa_1k_vcf}" \
  -nt 14 &>  "${trim_id_ma_fa_1k_log}"

#  --resourceAlleleConcordance true \

# Progress report
echo "Completed adding 1k AF values: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Train vqsr snp model --- #

# Progress report
echo "Started training vqsr snp model"

# File names
recal_snp="${vqsr_folder}/${dataset}_snp.recal"
plots_snp="${vqsr_folder}/${dataset}_snp_plots.R"
tranches_snp="${vqsr_folder}/${dataset}_snp.tranches"
log_train_snp="${logs_folder}/${dataset}_snp_train.log"

# Train vqsr snp model
"${java7}" -Xmx60g -jar "${gatk}" \
  -T VariantRecalibrator \
  -R "${ref_genome}" \
  -L "${broad_exomes_intervals}" -ip 100 \
  -input "${trim_id_ma_fa_1k_vcf}" \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 "${hapmap}" \
  -resource:omni,known=false,training=true,truth=true,prior=12.0 "${omni}" \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 "${phase1_1k_hc}" \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "${dbsnp_138}" \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS \
  -recalFile "${recal_snp}" \
  -tranchesFile "${tranches_snp}" \
  -rscriptFile "${plots_snp}" \
  --target_titv 3.2 \
  -mode SNP \
  -tranche 100.0 -tranche 99.0 -tranche 97.0 -tranche 95.0 -tranche 90.0 \
  -nt 14 &>  "${log_train_snp}"

# Excluded SOR and InbreedingCoeff because these annotations where missed 
# in all training variants in the dataset (-an SOR -an InbreedingCoeff)

# Progress report
echo "Completed training vqsr snp model: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Apply vqsr snp model --- #

# Progress report
echo "Started applying vqsr snp model"

# File names
vqsr_snp_vcf="${tmp_folder}/${dataset}_snp_vqsr.vcf"
log_apply_snp="${logs_folder}/${dataset}_snp_apply.log"

# Apply vqsr snp model
"${java7}" -Xmx60g -jar "${gatk}" \
  -T ApplyRecalibration \
  -R "${ref_genome}" \
  -L "${broad_exomes_intervals}" -ip 100 \
  -input "${trim_id_ma_fa_1k_vcf}" \
  -recalFile "${recal_snp}" \
  -tranchesFile "${tranches_snp}" \
  -o "${vqsr_snp_vcf}" \
  -mode SNP \
  -nt 14 &>  "${log_apply_snp}"  

# Progress report
echo "Completed applying vqsr snp model: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Train vqsr indel model --- #

# Progress report
echo "Started training vqsr indel model"

# File names
recal_indel="${vqsr_folder}/${dataset}_indel.recal"
plots_indel="${vqsr_folder}/${dataset}_indel_plots.R"
tranches_indel="${vqsr_folder}/${dataset}_indel.tranches"
log_train_indel="${logs_folder}/${dataset}_indel_train.log"

# Train model
"${java7}" -Xmx60g -jar "${gatk}" \
  -T VariantRecalibrator \
  -R "${ref_genome}" \
  -L "${broad_exomes_intervals}" -ip 100 \
  -input "${vqsr_snp_vcf}" \
  -resource:mills,known=false,training=true,truth=true,prior=12.0 "${mills}" \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "${dbsnp_138}" \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS \
  -recalFile "${recal_indel}" \
  -tranchesFile "${tranches_indel}" \
  -rscriptFile "${plots_indel}" \
  -tranche 100.0 -tranche 99.0 -tranche 97.0 -tranche 95.0 -tranche 90.0 \
  --maxGaussians 4 \
  -mode INDEL \
  -nt 14 &>  "${log_train_indel}"

# Excluded SOR and InbreedingCoeff because these annotations where missed 
# in all training variants in the dataset (-an SOR -an InbreedingCoeff)

# Progress report
echo "Completed training vqsr indel model: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Apply vqsr indel model --- #

# Progress report
echo "Started applying vqsr indel model"

# File names
out_vcf="${output_folder}/${dataset}_raw.vcf"
out_vcf_md5="${output_folder}/${dataset}_raw.md5"
log_apply_indel="${logs_folder}/${dataset}_indel_apply.log"

# Apply vqsr indel model
"${java7}" -Xmx60g -jar "${gatk}" \
  -T ApplyRecalibration \
  -R "${ref_genome}" \
  -L "${broad_exomes_intervals}" -ip 100 \
  -input "${vqsr_snp_vcf}" \
  -recalFile "${recal_indel}" \
  -tranchesFile "${tranches_indel}" \
  -o "${out_vcf}" \
  -mode INDEL \
  -nt 14 &>  "${log_apply_indel}"  

# Make md5 file
md5sum $(basename "${out_vcf}") $(basename "${out_vcf}.idx") > "${out_vcf_md5}"

# Progress report
echo "Completed applying vqsr indel model: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Prepare data for histograms --- #

# Progress report
echo "Started preparing data for histograms"

# File names
histograms_data_txt="${histograms_folder}/${dataset}_histograms_data.txt"
histograms_data_log="${logs_folder}/${dataset}_histograms_data.log"

# Prepare data
"${java7}" -Xmx60g -jar "${gatk}" \
  -T VariantsToTable \
  -R "${ref_genome}" \
  -L "${broad_exomes_intervals}" -ip 100 \
  -V "${out_vcf}" \
  -F RawVarID -F FILTER -F TYPE -F MultiAllelic \
  -F ALT_frequency_in_1k_90 -F ALT_frequency_in_1k_95 -F ALT_frequency_in_1k_99 -F ALT_frequency_in_1k_100 \
  -F CHROM -F POS -F REF -F ALT -F DP -F QUAL -F VQSLOD \
  -o "${histograms_data_txt}" \
  -AMD -raw &>  "${histograms_data_log}"  

# -AMD allow missed data
# -raw keep filtered

# Progress report
echo "Completed preparing data for histograms: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Generate histograms using R markdown script --- #

# Progress report
echo "Started making histograms"

# File names
histograms_report_pdf="${histograms_folder}/${dataset}_histograms_report.pdf"
histograms_report_html="${histograms_folder}/${dataset}_histograms_report.html"
histograms_plot_log="${logs_folder}/${dataset}_histograms_plot.log"

# Prepare R scripts
latex_dataset="${dataset//_/-}" # Underscores have special meaning in LaTex, so they should be avoided in PDF output

pdf_script="library('rmarkdown', lib='"${r_lib_folder}"'); render('"${scripts_folder}"/r01_make_pdf.Rmd', params=list(dataset='"${latex_dataset}-raw"' , data_file='"${histograms_data_txt}"'), output_file='"${histograms_report_pdf}"')"

html_script="library('rmarkdown', lib='"${r_lib_folder}"'); render('"${scripts_folder}"/r02_make_html.Rmd', params=list(dataset='"${dataset}-raw"' , working_folder='"${histograms_folder}"/' , data_file='"${histograms_data_txt}"'), output_file='"${histograms_report_html}"')"

# Execute R scripts
# Notes:
# Path to R was added to environment and modules required for 
# R with knitr were loaded in s01_genotype_gvcfs.sb.sh:
# module load gcc/5.2.0
# module load boost/1.50.0
# module load texlive/2015
# module load pandoc/1.15.2.1

# Underscore within ${dataset} may cause problem during rendering of pdf report

echo "-------------- Preparing pdf report -------------- " > "${histograms_plot_log}"
echo "" >> "${histograms_plot_log}"
"${r_bin_folder}/R" -e "${pdf_script}" &>> "${histograms_plot_log}"

echo "-------------- Preparing html report -------------- " >> "${histograms_plot_log}"
echo "" >> "${histograms_plot_log}"
"${r_bin_folder}/R" -e "${html_script}" &>> "${histograms_plot_log}"

echo "" >> "${histograms_plot_log}"

# Progress report
echo "Completed making histograms: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Calculating vcfstats for full data emitted by HC --- #

# Progress report
echo "Started vcfstats"
echo ""

# File name
vcf_stats="${all_vcfstats_folder}/${dataset}_raw.vchk"

# Calculate vcf stats
"${bcftools}" stats -F "${ref_genome}" "${out_vcf}" > "${vcf_stats}" 
#To be done: explore -R option to focus stats on targets:
# -R "${nextera_targets_bed}" ?? 

# Plot the stats
"${plot_vcfstats}" "${vcf_stats}" -p "${all_vcfstats_folder}/"
echo ""

# Progress report
echo "Completed vcfstats: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Make a "clean" copy of vcf without filtered variants --- #

# Progress report
echo "Started making clean vcf vithout HC filtered variants"

# File names
cln_vcf="${output_folder}/${dataset}_raw_cln.vcf"
cln_vcf_md5="${output_folder}/${dataset}_raw_cln.md5"
cln_vcf_log="${logs_folder}/${dataset}_raw_cln.log"

# Exclude filtered variants
"${java7}" -Xmx60g -jar "${gatk}" \
  -T SelectVariants \
  -R "${ref_genome}" \
  -L "${broad_exomes_intervals}" -ip 100 \
  -V "${out_vcf}" \
  -o "${cln_vcf}" \
  --excludeFiltered \
  -nt 14 &>  "${cln_vcf_log}"

# Make md5 file
md5sum $(basename "${cln_vcf}") $(basename "${cln_vcf}.idx") > "${cln_vcf_md5}"

# Completion message to log
echo "Completed making clean vcf: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Calculating vcfstats after minimal HC and VQSR filters--- #

# Progress report
echo "Started vcfstats on clean data"
echo ""

# File name
vcf_stats="${cln_vcfstats_folder}/${dataset}_cln.vchk"

# Calculate vcf stats
"${bcftools}" stats -F "${ref_genome}" "${cln_vcf}" > "${vcf_stats}" 
#To be done: explore -R option to focus stats on targets:
# -R "${nextera_targets_bed}" ?? 

# Plot the stats
"${plot_vcfstats}" "${vcf_stats}" -p "${cln_vcfstats_folder}/"
echo ""

# Progress report
echo "Completed vcfstats on clean data: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Copy output back to NAS --- #

# Progress report
echo "Started copying results to NAS"

# Remove temporary files from cluster
rm -fr "${tmp_folder}"

# Copy files to NAS
rsync -thrqe "ssh -x" "${output_folder}" "${data_server}:${project_location}/${project}/" 
exit_code="${?}"

# Stop if copying failed
if [ "${exit_code}" != "0" ]  
then
  echo ""
  echo "Failed copying results to NAS"
  echo "Script terminated"
  echo ""
  exit
fi

# Change ownership on nas (to allow user manipulating files later w/o administrative privileges)
ssh -x "${data_server}" "chown -R ${mgqnap_user}:${mgqnap_group} ${project_location}/${project}/$(basename ${output_folder})"
ssh -x "${data_server}" "chown -R ${mgqnap_user}:${mgqnap_group} ${project_location}/${project}" # just in case...
ssh -x "${data_server}" "chown -R ${mgqnap_user}:${mgqnap_group} ${project_location}" # just in case...

# Progress report
log_on_nas="${project_location}/${project}/$(basename ${output_folder})/logs/s01_rms_merge_vcfs.log"

ssh -x "${data_server}" "echo \"Completed copying results to NAS: $(date +%d%b%Y_%H:%M:%S)\" >> ${log_on_nas}"
ssh -x "${data_server}" "echo \"\" >> ${log_on_nas}"

# Remove results from cluster
#rm -fr "${logs_folder}"
rm -fr "${vqsr_folder}"
rm -fr "${histograms_folder}"
rm -fr "${vcfstats_folder}"

rm -f "${out_vcf}"
rm -f "${out_vcf}.idx"
rm -f "${out_vcf_md5}"

rm -f "${cln_vcf}"
rm -f "${cln_vcf}.idx"
rm -f "${cln_vcf_md5}"

ssh -x "${data_server}" "echo \"Removed results from cluster\" >> ${log_on_nas}"
ssh -x "${data_server}" "echo \"\" >> ${log_on_nas}"

# Return to the initial folder
cd "${init_dir}"

# Remove project folder (if requested)
if [ "${remove_project_folder}" == "yes" ] || [ "${remove_project_folder}" == "Yes" ] 
then 
  rm -fr "${project_folder}"
  ssh -x "${data_server}" "echo \"Removed working folder from cluster\" >> ${log_on_nas}"
else
  ssh -x "${data_server}" "echo \"Working folder is emptied and left on cluster\" >> ${log_on_nas}"
fi 
