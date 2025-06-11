#!/bin/bash
set -euo pipefail

BASE_DIR=~/Dropbox/Research/taguchi_sensei/yoshimoto_sensei/250611_variant_call
INPUT_DIR=/data2/F7317_vx0041/DragenEnrichment
OUT_DIR="${BASE_DIR}/vcfs/"

cd "${BASE_DIR}"
mkdir -p "${OUT_DIR}" "${IDX_DIR}"

bcftools norm \
  --check-ref e \
  -f ./files/GRCh38.primary_assembly.genome.fa \
  ./vcfs/Cx1.hard-filtered.vcf.gz \
  -Ou \
  -o /dev/null

# “.hard-filtered.vcf.gz” をまとめてコピーし，直接 tabix インデックス
find "${INPUT_DIR}" -type f -name '*.hard-filtered.vcf.gz' | while read -r gz; do
  base=$(basename "${gz}")
  echo "[INFO] Copy and index: ${base}"
  cp "${gz}" "${OUT_DIR}/${base}"
  tabix -p vcf "${OUT_DIR}/${base}"
done

genes="^(TAP1|TAP2|PSMB8|PSMB9|LMP7|ERAP1|TAPBP|PDIA3|CALR)$"

# BED 出力
awk '$3=="gene" \
    && match($0, /gene_name "([^"]+)"/, m) \
    && m[1] ~ "'"${genes}"'" \
  { print $1 "\t" $4-1 "\t" $5 "\t" m[1] }' \
  ./files/gencode.v44.annotation.gtf \
> ./files/target.bed

#annotation by VEP
# プロジェクト直下の files ディレクトリへ移動
cd ~/Dropbox/Research/taguchi_sensei/yoshimoto_sensei/250611_variant_call/files

# キャッシュ格納用ディレクトリを作成
mkdir -p vep_cache && cd vep_cache

# indexed_vep_cache（事前にtabixインデックス済み）を取得
wget ftp://ftp.ensembl.org/pub/release-114/variation/indexed_vep_cache/homo_sapiens_merged_vep_114_GRCh38.tar.gz

# 展開
tar zxvf homo_sapiens_merged_vep_114_GRCh38.tar.gz

# （不要ならアーカイブを削除）
rm homo_sapiens_merged_vep_114_GRCh38.tar.gz

cd ../files/vep_cache/ExAC
wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/subsets/ExAC_nonTCGA.r0.3.1.sites.vep.vcf
bgzip -f ExAC_nonTCGA.r0.3.1.sites.vep.vcf
tabix -p vcf ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz

rm -r target
mkdir -p target
for gz in ./vcfs/*.vcf.gz; do
  sample=$(basename "$gz" .vcf.gz)
  bcftools view \
    -R ./files/target.bed \
    -Ov \
    "$gz" \
  > ./target/${sample}.vcf
done

mkdir -p mafs
for vcf in ~/Dropbox/Research/taguchi_sensei/yoshimoto_sensei/250611_variant_call/target/*.hard-filtered.vcf; do
    sample=$(basename "$vcf" .hard-filtered.vcf)
    ~/sw/miniconda3/bin/vcf2maf.pl \
    --input-vcf "$vcf" \
    -output-maf ~/Dropbox/Research/taguchi_sensei/yoshimoto_sensei/250610_variant_call/mafs/${sample}.maf \
    --ref-fasta ./files/GRCh38.primary_assembly.genome.fa \
    --vep-path ~/sw/miniconda3/bin \
    --vep-data ./files/vep_cache \
    --ncbi-build GRCh38 \
    --species homo_sapiens \
    --tumor-id "${sample}"
done
