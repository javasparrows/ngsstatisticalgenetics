#!/bin/bash
# GVCF→VCF変換（染色体ごとに分割実行）

set -e
echo "=== Running joint genotyping (GVCF to VCF) ==="

# 染色体リスト
chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)

for chr in "${chromosomes[@]}"; do
  echo "Processing chromosome: $chr"
  
  # GenomicsDBImport
  gatk --java-options "-Xmx4g" GenomicsDBImport \
    --variant /workspace/intermediate/son.g.vcf \
    --variant /workspace/intermediate/father.g.vcf \
    --variant /workspace/intermediate/mother.g.vcf \
    --reference /workspace/materials/JG.fa \
    --genomicsdb-workspace-path /workspace/intermediate/genomics_database.${chr} \
    --intervals ${chr}
  
  # GenotypeGVCFs
  gatk --java-options "-Xmx4g" GenotypeGVCFs \
    --reference /workspace/materials/JG.fa \
    --variant gendb:///workspace/intermediate/genomics_database.${chr} \
    --output /workspace/intermediate/joint_genotyped.${chr}.vcf \
    --intervals ${chr}
done

echo "Merging chromosome-wise VCFs..."
# 染色体ごとのVCFを統合する
input_files=""
for chr in "${chromosomes[@]}"; do
  input_files+=" --INPUT /workspace/intermediate/joint_genotyped.${chr}.vcf"
done

gatk --java-options "-Xmx4G" MergeVcfs \
  --OUTPUT /workspace/results/merged.vcf.gz \
  ${input_files}

echo "=== Joint genotyping completed ==="