#!/bin/bash
# s01_filter_rms_vcf.sh
# Alexey Larionov, 03Jun2016
# Last updated: 09Jun2016

# Envirnoment
set -e

java8="/media/ajme/ajme/tools/java/jre1.8.0_91/bin/java"
picard241="/media/ajme/ajme/tools/picard/picard-tools-2.4.1/picard.jar"
gatk36="/media/ajme/ajme/tools/gatk/gatk3.6/GenomeAnalysisTK.jar"

ref_genome="/media/ajme/ajme/resources/ref_genomes/hg19/ucsc.hg19.fasta"
ref_genome_dict="/media/ajme/ajme/resources/ref_genomes/hg19/ucsc.hg19.dict"

source_vcfs_folder="/media/ajme/ajme/external/rms_db_vcf/rms120_vcf"
reordered_vcfs_folder="/media/ajme/ajme/external/rms_db_vcf/s01_reordered_vcfs"
sorted_vcfs_folder="/media/ajme/ajme/external/rms_db_vcf/s02_sorted_vcfs"
trimmed_vcfs_folder="/media/ajme/ajme/external/rms_db_vcf/s03_trimmed_vcfs"
filtered_vcfs_folder="/media/ajme/ajme/external/rms_db_vcf/s04_filtered_vcfs"
scripts_folder="/media/ajme/ajme/external/rms_db_vcf/scripts"

samples_file="rms_files_samples_AL.MGcorr.07Jun2016.txt"
samples_file="${data_folder}/${samples_file}"

samples_file="rms_files_samples_AL.MGcorr.07Jun2016.txt"
samples_file="${source_vcfs_folder}/${samples_file}"

mkdir -p "${reordered_vcfs_folder}"
mkdir -p "${sorted_vcfs_folder}"
mkdir -p "${trimmed_vcfs_folder}"
mkdir -p "${filtered_vcfs_folder}"

# Progress report
echo "Filter rms vcfs"
echo "Started: $(date +%d%b%Y_%H:%M:%S)"
echo ""
echo "Retains variants present in blood and tumour"
echo "Additional filters: DP>10 & QUAL >20"
echo ""
echo "samples_file: ${samples_file}"
echo ""
echo "source_vcfs_folder: ${source_vcfs_folder}"
echo "reordered_vcfs_folder: ${reordered_vcfs_folder}"
echo "sorted_vcfs_folder: ${sorted_vcfs_folder}"
echo "trimmed_vcfs_folder: ${trimmed_vcfs_folder}"
echo "filtered_vcfs_folder: ${filtered_vcfs_folder}"
echo "scripts_folder: ${scripts_folder}"
echo ""
echo "java8: ${java8}"
echo "picard241: ${picard241}"
echo "gatk36: ${gatk36}"
echo ""
echo "ref_genome: ${ref_genome}"
echo "ref_genome_dict: ${ref_genome_dict}"
echo ""
echo "--------------- Started filtering ---------------"
echo ""

# For each line in the samples file
while read file blood tumour patient
do
  
  # Skip 1st line
  if [ "${file}" == "file" ]
  then
    continue
  fi

  # Report settings
  echo "file: ${file}"
  echo "blood: ${blood}"
  echo "tumour: ${tumour}"
  echo "patient: ${patient}"
  echo ""

  # Compile file names
  vcf_in="${source_vcfs_folder}/${file}"
  vcf_reod="${reordered_vcfs_folder}/${patient}.reod.vcf"
  vcf_sort="${sorted_vcfs_folder}/${patient}.sort.vcf"
  vcf_trim="${trimmed_vcfs_folder}/${patient}.trim.vcf"
  vcf_filt="${filtered_vcfs_folder}/${patient}.filt.vcf"
  
  vcf_tmp1=$(mktemp --suffix=".tmp1.vcf" ${patient}.XXXXX)
  vcf_tmp2=$(mktemp --suffix=".tmp2.vcf" ${patient}.XXXXX)
  
  # Report file names
  echo "vcf_in: ${vcf_in}"
  echo "vcf_reod: ${vcf_reod}"
  echo "vcf_sort: ${vcf_sort}"
  echo "vcf_trim: ${vcf_trim}"  
  echo "vcf_filt: ${vcf_filt}"
  echo ""

  # Reorder contigs (for compartibility with GATK h19 ref genome)
  "${scripts_folder}/a01_reorder_rms_vcf.sh" \
    "${vcf_in}" \
    "${vcf_reod}"
  
  # Sort variants and validate against reference genome dictionary using picard
  echo ""
  echo "### Sort variants and validate against reference genome dictionary using picard ###"
  echo ""
  "${java8}" -Xms20g -Xmx60g -jar "${picard241}" \
    SortVcf \
    I="${vcf_reod}" \
    O="${vcf_sort}" \
    SEQUENCE_DICTIONARY="${ref_genome_dict}" \
    CREATE_INDEX=false # Index for sorted vcfs will be created by GATK at the next step

  # Remove alt alleles non-supported by genotypes
  echo ""
  echo "### Remove alt alleles non-supported by genotypes ###"
  echo ""
  "${java8}" -Xms20g -Xmx60g -jar "${gatk36}" \
    -T SelectVariants \
    -R "${ref_genome}" \
    -V "${vcf_sort}" \
    -o "${vcf_trim}" \
    --excludeNonVariants \
    --removeUnusedAlternates

  # Validate source vcf
  echo ""
  echo "### Validate preprocessed vcf file by GATK ###"
  echo ""
  "${java8}" -Xms20g -Xmx60g -jar "${gatk36}" \
     -T ValidateVariants \
     -R "${ref_genome}" \
     -V "${vcf_trim}"
  
  # Select variants
  
  echo ""
  echo "### Select variants present in blood ###"
  echo ""
  "${java8}" -Xms20g -Xmx60g -jar "${gatk36}" \
     -T SelectVariants \
     -R "${ref_genome}" \
     -V "${vcf_trim}" \
     -o "${vcf_tmp1}" \
     -select 'vc.getGenotype("'"${blood}"'").isHomVar() || vc.getGenotype("'"${blood}"'").isHet()'

  echo ""
  echo "### Select variants present in tumour ###"
  echo ""
  "${java8}" -Xms20g -Xmx60g -jar "${gatk36}" \
     -T SelectVariants \
     -R "${ref_genome}" \
     -V "${vcf_tmp1}" \
     -o "${vcf_tmp2}" \
     -select 'vc.getGenotype("'"${tumour}"'").isHomVar() || vc.getGenotype("'"${tumour}"'").isHet()'
  
  echo ""
  echo "### Select variants with DP >= 10 && QUAL >= 20.0 ###"
  echo ""
  "${java8}" -Xms20g -Xmx60g -jar "${gatk36}" \
     -T SelectVariants \
     -R "${ref_genome}" \
     -V "${vcf_tmp2}" \
     -o "${vcf_filt}" \
     -select 'DP >= 10 && QUAL >= 20.0'
  
  # Remove temporary files
  rm "${vcf_tmp1}" "${vcf_tmp1}.idx" "${vcf_tmp2}" "${vcf_tmp2}.idx"
  
  # Progress reoport
  echo ""
  echo "--------------- done ${patient} ---------------"
  echo ""
  
done < "${samples_file}" # next line in the samples file

# Completipon message
echo "Completed all samples: $(date +%d%b%Y_%H:%M:%S)"
echo ""
