#!/bin/bash
# コード12: アライメントの実行

set -e
echo "=== Running alignment ==="

echo "Running BWA-MEM2 alignment for Son..."
bwa-mem2 mem -R "@RG\tID:son\tSM:son\tPL:Illumina\tLB:son" \
  /workspace/materials/JG.fa \
  /workspace/materials/HG005_Son.R1.10X.fastq.gz \
  /workspace/materials/HG005_Son.R2.10X.fastq.gz \
  > /workspace/intermediate/son.sam

echo "Running BWA-MEM2 alignment for Father..."
bwa-mem2 mem -R "@RG\tID:father\tSM:father\tPL:Illumina\tLB:father" \
  /workspace/materials/JG.fa \
  /workspace/materials/HG006_Father.R1.10X.fastq.gz \
  /workspace/materials/HG006_Father.R2.10X.fastq.gz \
  > /workspace/intermediate/father.sam

echo "Running BWA-MEM2 alignment for Mother..."
bwa-mem2 mem -R "@RG\tID:mother\tSM:mother\tPL:Illumina\tLB:mother" \
  /workspace/materials/JG.fa \
  /workspace/materials/HG007_Mother.R1.10X.fastq.gz \
  /workspace/materials/HG007_Mother.R2.10X.fastq.gz \
  > /workspace/intermediate/mother.sam

echo "=== Alignment completed ==="