#!/bin/bash

# a00_start_pipeline.sh
# Start merging rms dbGAP vcfs
# Alexey Larionov, 10Jun2016
# Last updated: 24Jun2016

# Command line arguments
job_file="${1}"
scripts_folder="${2}"

# Main log name and location
working_folder=$(awk '$1=="working_folder:" {print $2}' "${job_file}") # e.g. /scratch/medgen/users/mae 
project=$(awk '$1=="project:" {print $2}' "${job_file}") # e.g. rms_dbGAP 
output_folder=$(grep ^"output_folder" "${job_file}" | awk '{print $2}') # e.g. s05_merged_vcf_v1
logs_folder=$(awk '$1=="logs_folder:" {print $2}' "${job_file}") # e.g. logs

logs_folder="${working_folder}/${project}/${output_folder}/${logs_folder}"
mkdir -p "${logs_folder}"
log="${logs_folder}/s01_rms_merge_vcfs.log"

# Requested time and account
time_to_request=$(grep ^"Max_time_to_request_(hrs.min.sec):" "${job_file}" | awk '{print $2}')
time_to_request=${time_to_request//./:} # substitute dots to colons
account_to_use=$(grep ^"Account_to_use_on_HPC" "${job_file}" | awk '{print $2}')

# Submit job
sbatch \
  --output="${log}" \
  --time="${time_to_request}" \
  --account="${account_to_use}" \
  "${scripts_folder}/s01_rms_merge_vcfs.sh" \
  "${job_file}" \
  "${scripts_folder}"
