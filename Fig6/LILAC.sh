#!/bin/bash

samples=$(find /data2/WES/data/ -maxdepth 1 -type d -name "Cx*" -printf "%f\n")

# 処理を関数化
run_lilac() {
    local sample=$1
    java -jar ~/sw/lilac_v1.5.2.jar \
        -sample "${sample}" \
        -ref_genome /data1/genomes/GRCh38_ENCODE/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
        -ref_genome_version V38 \
        -resource_dir /data1/HLA/241027_Lilac/resource/v6_0/ref/38/immune \
        -tumor_bam "/data2/WES/data/${sample}/enrichment_seq/${sample}_tumor.bam" \
        -output_dir "./res/${sample}/"
}

export -f run_lilac

# 並列実行
parallel run_lilac ::: "${samples[@]}"
