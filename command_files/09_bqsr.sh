#!/bin/bash
# コード14: BQSR（Base Quality Score Recalibration）

set -e
echo "=== Running BQSR ==="

echo "Running BaseRecalibrator for Son..."
gatk BaseRecalibrator \
  --reference /workspace/materials/JG.fa \
  --input /workspace/results/son.sort.markdup.bam \
  --known-sites /workspace/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.dbsnp_138.b37.success.sorted.vcf.gz \
  --known-sites /workspace/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.Mills_and_1000G_gold_standard.indels.b37.success.sorted.vcf.gz \
  --output /workspace/intermediate/son_recal_data.table

echo "Running BaseRecalibrator for Father..."
gatk BaseRecalibrator \
  --reference /workspace/materials/JG.fa \
  --input /workspace/results/father.sort.markdup.bam \
  --known-sites /workspace/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.dbsnp_138.b37.success.sorted.vcf.gz \
  --known-sites /workspace/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.Mills_and_1000G_gold_standard.indels.b37.success.sorted.vcf.gz \
  --output /workspace/intermediate/father_recal_data.table

echo "Running BaseRecalibrator for Mother..."
gatk BaseRecalibrator \
  --reference /workspace/materials/JG.fa \
  --input /workspace/results/mother.sort.markdup.bam \
  --known-sites /workspace/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.dbsnp_138.b37.success.sorted.vcf.gz \
  --known-sites /workspace/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.Mills_and_1000G_gold_standard.indels.b37.success.sorted.vcf.gz \
  --output /workspace/intermediate/mother_recal_data.table

echo "Applying BQSR for Son..."
gatk ApplyBQSR \
  --reference /workspace/materials/JG.fa \
  --input /workspace/results/son.sort.markdup.bam \
  --bqsr-recal-file /workspace/intermediate/son_recal_data.table \
  --create-output-bam-index true \
  --output /workspace/results/son.bqsr.bam

echo "Applying BQSR for Father..."
gatk ApplyBQSR \
  --reference /workspace/materials/JG.fa \
  --input /workspace/results/father.sort.markdup.bam \
  --bqsr-recal-file /workspace/intermediate/father_recal_data.table \
  --create-output-bam-index true \
  --output /workspace/results/father.bqsr.bam

echo "Applying BQSR for Mother..."
gatk ApplyBQSR \
  --reference /workspace/materials/JG.fa \
  --input /workspace/results/mother.sort.markdup.bam \
  --bqsr-recal-file /workspace/intermediate/mother_recal_data.table \
  --create-output-bam-index true \
  --output /workspace/results/mother.bqsr.bam

echo "=== BQSR completed ==="