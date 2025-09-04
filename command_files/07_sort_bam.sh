#!/bin/bash
# コード13: BAMファイルのソートとインデックス作成

set -e
echo "=== Sorting BAM files ==="

echo "Sorting Son BAM..."
samtools sort -o /workspace/intermediate/son.sort.bam /workspace/intermediate/son.bam

echo "Sorting Father BAM..."
samtools sort -o /workspace/intermediate/father.sort.bam /workspace/intermediate/father.bam

echo "Sorting Mother BAM..."
samtools sort -o /workspace/intermediate/mother.sort.bam /workspace/intermediate/mother.bam

echo "Creating indices for sorted BAM files..."
samtools index /workspace/intermediate/son.sort.bam
samtools index /workspace/intermediate/father.sort.bam
samtools index /workspace/intermediate/mother.sort.bam

echo "=== BAM sorting and indexing completed ==="