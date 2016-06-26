#!/bin/bash
# a01_reorder_rms_vcf.sh
# Alexey Larionov, 03Jun2016
# Last updated: 08Jun2016

# Stop at any error
set -e 

# Read file names
vcf_in="${1}"
vcf_out="${2}"

# Check input
if [ -z "${1}" ] | [ -z "${2}" ]
then
  echo "In-house script for vcf re-ordering and header editing"
  echo "Wrong input"
  echo "Use: a01_reorder_rms_vcf.sh vcf_in vcf_out"
  exit
fi

if [ ! -e "${1}" ]
then
  echo "In-house script for vcf re-ordering and header editing"
  echo "Input vcf file does not exist"
  echo "Use: a01_reorder_rms_vcf.sh vcf_in vcf_out"
  exit
fi

# Print info message
echo "#------------------------------------------------------------------#"
echo "#      In-house script for vcf re-ordering and header editing      #"
echo "#------------------------------------------------------------------#"
echo ""
echo "vcf_in: ${vcf_in}"
echo "vcf_out: ${vcf_out}"
echo ""

# Tools and resourses
java8="/media/ajme/ajme/tools/java/jre1.8.0_91/bin/java"
picard241="/media/ajme/ajme/tools/picard/picard-tools-2.4.1/picard.jar"
gatk36="/media/ajme/ajme/tools/gatk/gatk3.6/GenomeAnalysisTK.jar"
ref_genome_dict="/media/ajme/ajme/resources/ref_genomes/hg19/ucsc.hg19.dict"

# Chromosomes contigs in GATK-bundle compartible order 
contigs_list="chrM chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

# --- Reshape the header --- #

# Progress report
echo "Processing header"
echo ""

# Remove mis-ordered reference contigs from the header
grep "^##" "${vcf_in}" | \
  grep -v "^##contig=" | \
  grep -v "^##reference=" > "${vcf_out}"

# Check and add chromosome contigs in the given order
for chr in ${contigs_list}
do
  if grep "^##contig=<ID=${chr}," "${vcf_in}" >> "${vcf_out}"
  then
    echo "${chr} - OK"
  else
    echo "${chr} - Contig missed"
  fi
done
echo ""

# Add patches and fixes in required order
# This is done for files formate compartibility ionly
# These contigs are present in RMS dataset anyway 
echo "##contig=<ID=chr1_gl000191_random,length=106433,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr1_gl000192_random,length=547496,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr4_ctg9_hap1,length=590426,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr4_gl000193_random,length=189789,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr4_gl000194_random,length=191469,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr6_apd_hap1,length=4622290,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr6_cox_hap2,length=4795371,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr6_dbb_hap3,length=4610396,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr6_mann_hap4,length=4683263,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr6_mcf_hap5,length=4833398,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr6_qbl_hap6,length=4611984,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr6_ssto_hap7,length=4928567,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr7_gl000195_random,length=182896,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr8_gl000196_random,length=38914,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr8_gl000197_random,length=37175,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr9_gl000198_random,length=90085,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr9_gl000199_random,length=169874,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr9_gl000200_random,length=187035,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr9_gl000201_random,length=36148,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr11_gl000202_random,length=40103,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr17_ctg5_hap1,length=1680828,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr17_gl000203_random,length=37498,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr17_gl000204_random,length=81310,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr17_gl000205_random,length=174588,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr17_gl000206_random,length=41001,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr18_gl000207_random,length=4262,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr19_gl000208_random,length=92689,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr19_gl000209_random,length=159169,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chr21_gl000210_random,length=27682,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000211,length=166566,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000212,length=186858,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000213,length=164239,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000214,length=137718,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000215,length=172545,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000216,length=172294,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000217,length=172149,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000218,length=161147,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000219,length=179198,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000220,length=161802,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000221,length=155397,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000222,length=186861,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000223,length=180455,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000224,length=179693,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000225,length=211173,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000226,length=15008,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000227,length=128374,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000228,length=129120,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000229,length=19913,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000230,length=43691,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000231,length=27386,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000232,length=40652,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000233,length=45941,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000234,length=40531,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000235,length=34474,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000236,length=41934,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000237,length=45867,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000238,length=39939,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000239,length=33824,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000240,length=41933,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000241,length=42152,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000242,length=43523,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000243,length=43341,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000244,length=39929,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000245,length=36651,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000246,length=38154,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000247,length=36422,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000248,length=39786,assembly=hg19>" >> "${vcf_out}"
echo "##contig=<ID=chrUn_gl000249,length=38502,assembly=hg19>" >> "${vcf_out}"
echo "##ReferenceContigsSource=InHouseScript" >> "${vcf_out}"

# The #chrom line
grep "^#CHROM[[:blank:]]" "${vcf_in}" >> "${vcf_out}"

# Reorder the contigs in the vcf file body
echo "Processing variants"
echo ""

for chr in ${contigs_list}
do
  if grep "^${chr}[[:blank:]]" "${vcf_in}" >> "${vcf_out}"
  then
    echo "${chr} - OK"
  else
    echo "${chr} - No variants found"
  fi
done
echo ""

# Completion message
echo "#------------------------------------------------------------------#"
echo "# Completed in-house script for vcf re-ordering and header editing #"
echo "#------------------------------------------------------------------#"
