#!/bin/sh
#Calculates the linkage disequilibrium(LD) between a top snp and the variants within a specific window using PLINK
# Usage: bash extractLD_bfile.sh plink_path bfile_prefix keep_file chr pos_start pos_end top_snp out_file


plink_path=$1
bfile_prefix=$2
chr=$3
pos_start=$4
pos_end=$5
top_snp=$6
window_kb=$7
out_file=$8

#Create the output directory if it does not exist
mkdir -p ld_files

#Run PLINK to calculate LD

${plink_path}/plink \
    --bfile ${bfile_prefix} \
    --allow-extra-chr \
    --chr ${chr} \
    --from-bp ${pos_start} \
    --to-bp ${pos_end} \
    --r2 \
    --ld-window 100000000 \
    --ld-window-kb {window_kb}\
    --ld-window-r2 0 \
    --ld-snp ${top_snp} \
    --keep-allele-order \
    --out ld_files/${out_file}
