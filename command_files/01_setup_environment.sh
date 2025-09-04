#!/bin/bash
# コード1: 環境設定とフォルダ構造の作成

set -e
echo "=== Setting up environment and directory structure ==="

# Create directory structure (already exists in Docker)
echo "Directory structure already created in Docker container"

# Check Java version
echo "Checking Java version..."
java -version

# Check Python version
echo "Checking Python version..."
python3 --version

# Check tools installation
echo "Checking samtools..."
samtools --version

echo "Checking bwa-mem2..."
bwa-mem2 version

echo "Checking GATK..."
gatk --version

echo "=== Environment setup completed ==="