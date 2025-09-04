#!/bin/bash
# コード9: JG2.1.0参照ゲノムとリソースファイルのダウンロード

set -e
echo "=== Downloading reference genome and resource files ==="

cd /workspace/materials

echo "Downloading JG2.1.0 FASTA file..."
wget https://jmorp.megabank.tohoku.ac.jp/datasets/tommo-jg2.1.0-20211208/files/jg2.1.0.fa.gz
gzip -cd jg2.1.0.fa.gz > JG.fa && rm jg2.1.0.fa.gz

echo "Downloading GATK Resource Bundle for JG2.1.0..."
wget https://jmorp.megabank.tohoku.ac.jp/datasets/tommo-jg2.1.0-20211208/files/JG2.1.0-ResourceBundle-from-b37.zip
unzip JG2.1.0-ResourceBundle-from-b37.zip && rm JG2.1.0-ResourceBundle-from-b37.zip

echo "=== Reference genome download completed ==="