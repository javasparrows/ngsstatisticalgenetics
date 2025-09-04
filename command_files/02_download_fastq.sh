#!/bin/bash
# コード8: FASTQファイルのダウンロード

set -e
echo "=== Downloading FASTQ files ==="

cd /workspace/materials

echo "Downloading 10x FASTQ files (each ~10GB)..."
echo "This may take several hours depending on network speed."

# 10xのFASTQファイルをaria2cで高速ダウンロード
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310196/files/HG005_Son.R1.10X.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310196/files/HG005_Son.R2.10X.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310202/files/HG006_Father.R1.10X.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310202/files/HG006_Father.R2.10X.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310233/files/HG007_Mother.R1.10X.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310233/files/HG007_Mother.R2.10X.fastq.gz

echo "=== FASTQ download completed ==="