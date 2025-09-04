#!/bin/bash
# コード12: インデックスファイルの作成

set -e
echo "=== Preparing reference genome indices ==="

cd /workspace/materials

echo "Creating BWA-MEM2 index (this may take 10+ minutes and require significant RAM)..."
bwa-mem2 index JG.fa

echo "Creating FASTA index..."
samtools faidx JG.fa

echo "Creating sequence dictionary..."
gatk CreateSequenceDictionary --REFERENCE JG.fa

echo "=== Reference preparation completed ==="