#!/bin/bash
# 重複したリードを除去するMarkDuplicates

set -e
echo "=== Marking duplicates ==="

echo "Marking duplicates for Son..."
gatk MarkDuplicates \
  --INPUT /workspace/intermediate/son.sort.bam \
  --OUTPUT /workspace/results/son.sort.markdup.bam \
  --METRICS_FILE /workspace/results/son.sort.markdup.metrics.txt \
  --CREATE_INDEX true

echo "Marking duplicates for Father..."
gatk MarkDuplicates \
  --INPUT /workspace/intermediate/father.sort.bam \
  --OUTPUT /workspace/results/father.sort.markdup.bam \
  --METRICS_FILE /workspace/results/father.sort.markdup.metrics.txt \
  --CREATE_INDEX true

echo "Marking duplicates for Mother..."
gatk MarkDuplicates \
  --INPUT /workspace/intermediate/mother.sort.bam \
  --OUTPUT /workspace/results/mother.sort.markdup.bam \
  --METRICS_FILE /workspace/results/mother.sort.markdup.metrics.txt \
  --CREATE_INDEX true

echo "=== Mark duplicates completed ==="