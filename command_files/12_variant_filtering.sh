#!/bin/bash
# コード17: ハードフィルタリング

set -e
echo "=== Running variant filtering ==="

echo "Splitting variants into SNVs and INDELs..."
# SNVの抽出
gatk --java-options "-Xmx4g" SelectVariants \
  --reference /workspace/materials/JG.fa \
  --variant /workspace/results/merged.vcf.gz \
  --select-type-to-include SNP \
  --output /workspace/intermediate/merged.snv.vcf.gz

# Indelの抽出
gatk --java-options "-Xmx4g" SelectVariants \
  --reference /workspace/materials/JG.fa \
  --variant /workspace/results/merged.vcf.gz \
  --select-type-to-include INDEL \
  --output /workspace/intermediate/merged.indel.vcf.gz

echo "Applying hard filtering for SNVs..."
gatk --java-options "-Xmx4g" VariantFiltration \
  --reference /workspace/materials/JG.fa \
  --variant /workspace/intermediate/merged.snv.vcf.gz \
  --output /workspace/intermediate/merged.snv.hardfiltering.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "SOR > 3.0" --filter-name "SOR3" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5"

echo "Applying hard filtering for INDELs..."
gatk --java-options "-Xmx4g" VariantFiltration \
  --reference /workspace/materials/JG.fa \
  --variant /workspace/intermediate/merged.indel.vcf.gz \
  --output /workspace/intermediate/merged.indel.hardfiltering.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 200.0" --filter-name "FS200" \
  --filter-expression "SOR > 10.0" --filter-name "SOR10" \
  --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

echo "Merging filtered VCFs..."
gatk --java-options "-Xmx4G" MergeVcfs \
  --OUTPUT /workspace/results/merged.hardfiltering.vcf.gz \
  --INPUT /workspace/intermediate/merged.snv.hardfiltering.vcf.gz \
  --INPUT /workspace/intermediate/merged.indel.hardfiltering.vcf.gz

echo "=== Variant filtering completed ==="