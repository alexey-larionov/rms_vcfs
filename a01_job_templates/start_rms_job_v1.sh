#!/bin/bash

# start_rms_job.sh
# Start rms job described in the job file
# Alexey Larionov, 10Jun2016
# Version 1

# Use: 
# start_rms_job.sh rms_job_file

# Get job file
job="${1}"

# ------------ Check input ------------ #

# Check that an argument has been provided
if [ -z "${job}" ]
then
  echo "" 
  echo "No job file given"
  echo "" 
  echo "Use:"
  echo "start_rms_job.sh job_file"
  echo ""  
  echo "Script terminated"
  echo ""
  exit 1
fi

# Help message
if [ "${job}" == "-h" ] || [ "${job}" == "--help" ]
then
  echo ""
  echo "Start data analysis described in the job file"
  echo ""
  echo "Use:"
  echo "start_rms_job.sh rms_job_file"
  echo ""  
  exit
fi

# Make full file name for the job description file
job_file="$(pwd)/${job}"

# Check that job description file exists
if [ ! -e "${job_file}" ]
then
  echo "" 
  echo "Job file ${job_file} does not exist"
  echo ""  
  echo "Use:"
  echo "start_rms_job.sh rms_job_file"
  echo ""
  echo "Script terminated"
  echo ""
  exit 1
fi

# Check the job description file format (just check the first line only)
read line1 < "${job_file}"
if [ "${line1}" != "Job description file for merging rms dbGAP vcfs" ] && \
   [ "${line1}" != "Job description file for filtering rms vcf" ] && \
   [ "${line1}" != "Job description file for rms vep annotation" ] && \
   [ "${line1}" != "Job description file for rms data export" ]
then
  echo "" 
  echo "Unexpected format of the job file ${job_file}"
  echo ""
  echo "Use:"
  echo "start_job.sh job_file"
  echo ""
  echo "Script terminated"
  echo "" 
  exit 1
fi

# ------------ Start pipeline ------------ #

# Get start script name and folder from the job file
scripts_folder=$(awk '$1=="scripts_folder:" {print $2}' "${job_file}")
start_script=$(awk '$1=="start_script:" {print $2}' "${job_file}")

# Ask user to confirm the job before launching

echo ""
echo "Requested job:"
echo ""
echo "Pipeline launching script: ${scripts_folder}/${start_script}"
echo "Job description file: ${job_file}"
echo ""
echo "Start this job? y/n"
read user_choice

if [ "${user_choice}" != "y" ]
then
  echo ""
  echo "Script terminated"
  echo ""
  exit
fi

# Start the job
echo ""
"${scripts_folder}/${start_script}" "${job_file}" "${scripts_folder}" 
echo ""
