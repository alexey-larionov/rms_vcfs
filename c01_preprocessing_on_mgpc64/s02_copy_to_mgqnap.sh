#!/bin/bash
# s02_copy_to_mgqnap.sh
# Started: Alexey Larionov, 10Jun2016
# Last updated: AL, 10Jun2016

# Copy pre-processed rms vcfs to mgqnap

# Stop at any error
set -e 

# Start message
echo "Copy pre-processed rms vcfs to mgqnap"
echo "Started: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# Folders
source1="/media/ajme/ajme/external/rms_db_vcf/rms120_vcf"
source2="/media/ajme/ajme/external/rms_db_vcf/s03_trimmed_vcfs"
source3="/media/ajme/ajme/external/rms_db_vcf/s04_filtered_vcfs"
target="admin@mgqnap.medschl.cam.ac.uk:/share/mae/rms_dbGAP"

# Copying

echo "Raw data"
echo ""
rsync -avhe "ssh -x" "${source1}" "${target}/"

echo "Reshaped data"
echo ""
rsync -avhe "ssh -x" "${source2}" "${target}/"

echo "Filtered data"
echo ""
rsync -avhe "ssh -x" "${source3}" "${target}/"

# Completion message
echo ""
echo "Completed: $(date +%d%b%Y_%H:%M:%S)"
