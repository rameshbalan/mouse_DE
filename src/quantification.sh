#!/bin/bash

# Indexing the reference transcriptome
salmon index -t reference_rna.fa -i reference_rna_index

declare -a arr=(
"raw_data/C1M3L2"
"raw_data/C1M3L2"
"raw_data/C1M3R2"
"raw_data/C1M3R2"
"raw_data/C1M5L2"
"raw_data/C1M5L2"
"raw_data/C1M5R2"
"raw_data/C1M5R2"
"raw_data/C2M1L2"
"raw_data/C2M1L2"
"raw_data/C2M1R2"
"raw_data/C2M1R2"
"raw_data/C2M3L2"
"raw_data/C2M3L2"
"raw_data/C2M3R2"
"raw_data/C2M3R2"
"raw_data/C2M4L2"
"raw_data/C2M4L2"
"raw_data/C2M4R2"
"raw_data/C2M4R2"
"raw_data/C2M5L2"
"raw_data/C2M5L2"
"raw_data/C2M5R2"
"raw_data/C2M5R2"
)

# Looping through each file to quantify 
for fn in "${arr[@]}";
do
samp=`basename ${fn}`
echo "Processing sample ${fn}"
salmon quant -i reference_rna_index -l A \
         -1 ${fn}_1.fq.gz \
         -2 ${fn}_2.fq.gz \
         -p 16 -o quants/${samp}_quant
done 
