#!/bin/bash
# コード15: バリアントコーリング（BAM→GVCF）

set -e
echo "=== Running variant calling (BAM to GVCF) ==="

echo "Running HaplotypeCaller for Son..."
gatk --java-options "-Xmx4g" HaplotypeCaller \
  --reference /workspace/materials/JG.fa \
  --emit-ref-confidence GVCF \
  --input /workspace/results/son.sort.markdup.bam \
  --output /workspace/intermediate/son.g.vcf

echo "Running HaplotypeCaller for Father..."
gatk --java-options "-Xmx4g" HaplotypeCaller \
  --reference /workspace/materials/JG.fa \
  --emit-ref-confidence GVCF \
  --input /workspace/results/father.sort.markdup.bam \
  --output /workspace/intermediate/father.g.vcf

echo "Running HaplotypeCaller for Mother..."
gatk --java-options "-Xmx4g" HaplotypeCaller \
  --reference /workspace/materials/JG.fa \
  --emit-ref-confidence GVCF \
  --input /workspace/results/mother.sort.markdup.bam \
  --output /workspace/intermediate/mother.g.vcf

echo "=== Variant calling completed ==="