#!/bin/bash
cd ~/Dropbox/Research/taguchi_sensei/yoshimoto_sensei/250612_variant_call

# “.hard-filtered.vcf.gz” をコピー
mkdir vcfs
find /data2/F7317_vx0041/DragenEnrichment -type f -name '*.hard-filtered.vcf.gz' | while read -r gz; do
  base=$(basename "${gz}")
  echo "[INFO] Copy and index: ${base}"
  cp "${gz}" "./vcfs/${base}"
done
gunzip -d ./vcfs/*.gz

# vcf to maf
mkdir mafs
ls ./vcfs/*.hard-filtered.vcf | parallel -j 16 '
  sample=$(basename {} .hard-filtered.vcf)
  ~/sw/miniconda3/bin/vcf2maf.pl \
    --input-vcf {} \
    --output-maf ./mafs/${sample}.maf \
    --ref-fasta ~/Dropbox/Research/taguchi_sensei/yoshimoto_sensei/250610_variant_call/files/GRCh38.primary_assembly.genome.fa \
    --vep-path ~/sw/miniconda3/bin \
    --vep-data ~/Dropbox/Research/taguchi_sensei/yoshimoto_sensei/250610_variant_call/files/vep_cache \
    --ncbi-build GRCh38 \
    --species homo_sapiens \
    --tumor-id "${sample}"
'

