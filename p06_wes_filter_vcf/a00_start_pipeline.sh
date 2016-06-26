#!/bin/bash

# a00_start_pipeline.sh
# Start filtering rms vcf by DP, QUAL and VQSLOD
# Alexey Larionov, 06Feb2016
# Last updated: 24Jun2016

## Read parameter
job_file="${1}"
scripts_folder="${2}"

# Main log name and location
working_folder=$(awk '$1=="working_folder:" {print $2}' "${job_file}") # e.g. /scratch/medgen/users/mae 
project=$(awk '$1=="project:" {print $2}' "${job_file}") # e.g. rms_dbGAP
output_folder=$(grep ^"output_folder" "${job_file}" | awk '{print $2}') # e.g. s06_merged_vcf_v1_flt1
logs_folder=$(awk '$1=="logs_folder:" {print $2}' "${job_file}") # e.g. logs
fa_threshold=$(awk '$1=="Apply_1k_ALT_frequency_threshold:" {print $2}' "${job_file}") # e.g. no

logs_folder="${working_folder}/${project}/${output_folder}/${logs_folder}"
mkdir -p "${logs_folder}"
log="${logs_folder}/s01_filter_vcf.log"

# Requested time and account
time_to_request=$(grep ^"Max_time_to_request_(hrs.min.sec):" "${job_file}" | awk '{print $2}')
time_to_request=${time_to_request//./:} # substitute dots to colons
account_to_use=$(grep ^"Account_to_use_on_HPC" "${job_file}" | awk '{print $2}')

# Check the value of AF threshold
if [ "${fa_threshold}" != "no" ] && \
   [ "${fa_threshold}" != "90" ] && \
   [ "${fa_threshold}" != "95" ] && \
   [ "${fa_threshold}" != "99" ] && \
   [ "${fa_threshold}" != "100" ]
then
  echo "" 
  echo "Unexpected value for 1k ALT frequency threshold: ${fa_threshold}"
  echo ""
  echo "Allowed values: no, 90, 95, 99, 100"
  echo ""
  echo "Script terminated"
  echo "" 
  exit 1
fi

# Submit job
sbatch \
  --output="${log}" \
  --time="${time_to_request}" \
  --account="${account_to_use}" \
  "${scripts_folder}/s01_filter_vcf.sh" \
    "${job_file}" \
    "${scripts_folder}"
