#!/bin/bash
cd data/peritoneal_cavity/

ref_path="/home/test/software/refdata-gex-GRCm39-2024-A/"

# 直接指定样本列表
parallel -j 2 --eta --joblog cellranger_jobs.log \
    cellranger count \
    --id={} \
    --transcriptome="$ref_path" \
    --fastqs="$HOME/下载/mouse/peritoneal_cavity/{}" \
    --sample={} \
    --localcores=64 \
    --localmem=64 \
    --create-bam false \
    --nosecondary \
    ::: "HC_1" "HC_2"


