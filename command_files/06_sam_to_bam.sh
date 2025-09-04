#!/bin/bash
# SAMファイルをBAMファイルへ変換

set -e
echo "=== Converting SAM to BAM ==="

echo "Converting Son SAM to BAM..."
samtools view -bS /workspace/intermediate/son.sam > /workspace/intermediate/son.bam

echo "Converting Father SAM to BAM..."
samtools view -bS /workspace/intermediate/father.sam > /workspace/intermediate/father.bam

echo "Converting Mother SAM to BAM..."
samtools view -bS /workspace/intermediate/mother.sam > /workspace/intermediate/mother.bam

echo "=== SAM to BAM conversion completed ==="