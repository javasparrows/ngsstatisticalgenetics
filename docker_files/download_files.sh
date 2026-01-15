#!/bin/bash
set -e

# カラー出力の設定
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# ダウンロード進捗を表示する関数
print_header() {
    echo -e "\n${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}\n"
}

print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# ファイルサイズを表示する関数
show_file_size() {
    local url=$1
    local description=$2
    local estimated_size=$3

    echo -e "${YELLOW}File:${NC} $description"
    echo -e "${YELLOW}URL:${NC} $url"
    echo -e "${YELLOW}Estimated Size:${NC} ~${estimated_size}"
    echo ""
}

# ディレクトリ構造の作成
print_header "Creating Directory Structure"
cd /variant_call
mkdir -p tools materials intermediate results logs

print_status "Directory structure created successfully"

# =================================================================
# 1. ツールのダウンロードとインストール
# =================================================================
print_header "Downloading and Installing Tools"

# Samtools
print_status "Downloading Samtools (estimated: ~5 MB compressed)"
cd /variant_call/tools
if [ ! -f "samtools-1.16.1.tar.bz2" ]; then
    show_file_size "https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2" \
                   "Samtools 1.16.1" "5 MB"
    wget -c https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
    print_status "Samtools downloaded"
else
    print_warning "Samtools already exists, skipping download"
fi

# GATK
print_status "Downloading GATK (estimated: ~150 MB)"
cd /variant_call/tools
if [ ! -f "gatk-4.3.0.0.zip" ]; then
    show_file_size "https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip" \
                   "GATK 4.3.0.0" "150 MB"
    aria2c -x 8 -s 8 -c https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip
    unzip -q gatk-4.3.0.0.zip
    print_status "GATK downloaded and extracted"
else
    print_warning "GATK already exists, skipping download"
fi

# Strelka2
print_status "Downloading Strelka2 (estimated: ~40 MB)"
cd /variant_call/tools
if [ ! -f "strelka-2.9.10.centos6_x86_64.tar.bz2" ]; then
    show_file_size "https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2" \
                   "Strelka 2.9.10" "40 MB"
    wget -c https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
    tar -xjf strelka-2.9.10.centos6_x86_64.tar.bz2
    print_status "Strelka2 downloaded and extracted"
else
    print_warning "Strelka2 already exists, skipping download"
fi

# =================================================================
# 2. FASTQファイルのダウンロード（10X版を使用）
# =================================================================
print_header "Downloading FASTQ Files (10X Coverage)"
print_warning "Total estimated size: ~60 GB (6 files × ~10 GB each)"
print_status "This will take a significant amount of time depending on your internet connection"

cd /variant_call/materials

# HG005 Son
if [ ! -f "HG005_Son.R1.10X.fastq.gz" ]; then
    show_file_size "https://zenodo.org/record/7310196/files/HG005_Son.R1.10X.fastq.gz" \
                   "HG005 Son R1 (10X)" "10 GB"
    aria2c -x 8 -s 8 -c https://zenodo.org/record/7310196/files/HG005_Son.R1.10X.fastq.gz
    print_status "HG005 Son R1 downloaded"
else
    print_warning "HG005 Son R1 already exists, skipping"
fi

if [ ! -f "HG005_Son.R2.10X.fastq.gz" ]; then
    show_file_size "https://zenodo.org/record/7310196/files/HG005_Son.R2.10X.fastq.gz" \
                   "HG005 Son R2 (10X)" "10 GB"
    aria2c -x 8 -s 8 -c https://zenodo.org/record/7310196/files/HG005_Son.R2.10X.fastq.gz
    print_status "HG005 Son R2 downloaded"
else
    print_warning "HG005 Son R2 already exists, skipping"
fi

# HG006 Father
if [ ! -f "HG006_Father.R1.10X.fastq.gz" ]; then
    show_file_size "https://zenodo.org/record/7310202/files/HG006_Father.R1.10X.fastq.gz" \
                   "HG006 Father R1 (10X)" "10 GB"
    aria2c -x 8 -s 8 -c https://zenodo.org/record/7310202/files/HG006_Father.R1.10X.fastq.gz
    print_status "HG006 Father R1 downloaded"
else
    print_warning "HG006 Father R1 already exists, skipping"
fi

if [ ! -f "HG006_Father.R2.10X.fastq.gz" ]; then
    show_file_size "https://zenodo.org/record/7310202/files/HG006_Father.R2.10X.fastq.gz" \
                   "HG006 Father R2 (10X)" "10 GB"
    aria2c -x 8 -s 8 -c https://zenodo.org/record/7310202/files/HG006_Father.R2.10X.fastq.gz
    print_status "HG006 Father R2 downloaded"
else
    print_warning "HG006 Father R2 already exists, skipping"
fi

# HG007 Mother
if [ ! -f "HG007_Mother.R1.10X.fastq.gz" ]; then
    show_file_size "https://zenodo.org/record/7310233/files/HG007_Mother.R1.10X.fastq.gz" \
                   "HG007 Mother R1 (10X)" "10 GB"
    aria2c -x 8 -s 8 -c https://zenodo.org/record/7310233/files/HG007_Mother.R1.10X.fastq.gz
    print_status "HG007 Mother R1 downloaded"
else
    print_warning "HG007 Mother R1 already exists, skipping"
fi

if [ ! -f "HG007_Mother.R2.10X.fastq.gz" ]; then
    show_file_size "https://zenodo.org/record/7310233/files/HG007_Mother.R2.10X.fastq.gz" \
                   "HG007 Mother R2 (10X)" "10 GB"
    aria2c -x 8 -s 8 -c https://zenodo.org/record/7310233/files/HG007_Mother.R2.10X.fastq.gz
    print_status "HG007 Mother R2 downloaded"
else
    print_warning "HG007 Mother R2 already exists, skipping"
fi

# =================================================================
# 3. JG2.1.0参照ゲノムとリソースファイルのダウンロード
# =================================================================
print_header "Downloading JG2.1.0 Reference Genome and Resources"

cd /variant_call/materials

# JG2.1.0 FASTAファイル
if [ ! -f "JG.fa" ]; then
    show_file_size "https://jmorp.megabank.tohoku.ac.jp/datasets/tommo-jg2.1.0-20211208/files/jg2.1.0.fa.gz" \
                   "JG2.1.0 Reference Genome" "900 MB compressed, 3 GB uncompressed"
    aria2c -x 8 -s 8 -c https://jmorp.megabank.tohoku.ac.jp/datasets/tommo-jg2.1.0-20211208/files/jg2.1.0.fa.gz
    gzip -cd jg2.1.0.fa.gz > JG.fa && rm jg2.1.0.fa.gz
    print_status "JG2.1.0 reference genome downloaded and extracted"
else
    print_warning "JG.fa already exists, skipping"
fi

# JG2.1.0 Resource Bundle
if [ ! -d "JG2.1.0-ResourceBundle-from-b37" ]; then
    show_file_size "https://jmorp.megabank.tohoku.ac.jp/datasets/tommo-jg2.1.0-20211208/files/JG2.1.0-ResourceBundle-from-b37.zip" \
                   "JG2.1.0 Resource Bundle" "500 MB"
    aria2c -x 8 -s 8 -c https://jmorp.megabank.tohoku.ac.jp/datasets/tommo-jg2.1.0-20211208/files/JG2.1.0-ResourceBundle-from-b37.zip
    unzip -q JG2.1.0-ResourceBundle-from-b37.zip && rm JG2.1.0-ResourceBundle-from-b37.zip
    print_status "JG2.1.0 Resource Bundle downloaded and extracted"
else
    print_warning "JG2.1.0 Resource Bundle already exists, skipping"
fi

# =================================================================
# 4. SnpEff/SnpSiftのダウンロード
# =================================================================
print_header "Downloading SnpEff/SnpSift"

cd /variant_call/tools

if [ ! -d "snpEff" ]; then
    show_file_size "https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip" \
                   "SnpEff Latest Core" "200 MB"
    aria2c -x 8 -s 8 -c https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
    unzip -q snpEff_latest_core.zip
    print_status "SnpEff/SnpSift downloaded and extracted"
else
    print_warning "SnpEff already exists, skipping"
fi

# SnpEff JG2.1.0データベース
cd /variant_call/tools/snpEff
if [ ! -f "data/JG2.1.0/snpEffectPredictor.bin" ]; then
    mkdir -p data/JG2.1.0
    show_file_size "https://zenodo.org/record/7503342/files/snpEffectPredictor.bin" \
                   "SnpEff JG2.1.0 Database" "500 MB"
    cd data/JG2.1.0
    aria2c -x 8 -s 8 -c https://zenodo.org/record/7503342/files/snpEffectPredictor.bin
    cd ../..
    print_status "SnpEff JG2.1.0 database downloaded"
else
    print_warning "SnpEff JG2.1.0 database already exists, skipping"
fi

# =================================================================
# 5. vcfannoのダウンロード
# =================================================================
print_header "Downloading vcfanno"

cd /variant_call/tools

if [ ! -f "vcfanno" ]; then
    show_file_size "https://github.com/brentp/vcfanno/releases/download/v0.3.4/vcfanno_linux64" \
                   "vcfanno v0.3.4" "10 MB"
    wget -c https://github.com/brentp/vcfanno/releases/download/v0.3.4/vcfanno_linux64 -O vcfanno
    chmod +x vcfanno
    print_status "vcfanno downloaded"
else
    print_warning "vcfanno already exists, skipping"
fi

# =================================================================
# 6. slivarのダウンロード
# =================================================================
print_header "Downloading slivar"

cd /variant_call/tools

if [ ! -f "slivar" ]; then
    show_file_size "https://github.com/brentp/slivar/releases/download/v0.2.8/slivar" \
                   "slivar v0.2.8" "15 MB"
    wget -c https://github.com/brentp/slivar/releases/download/v0.2.8/slivar
    chmod +x slivar
    print_status "slivar downloaded"
else
    print_warning "slivar already exists, skipping"
fi

# slivar gnomADデータベース
cd /variant_call/materials
if [ ! -f "gnomad.hg38.genomes.v3.fix.zip" ]; then
    show_file_size "https://slivar.s3.amazonaws.com/gnomad.hg38.genomes.v3.fix.zip" \
                   "gnomAD v3 for slivar" "2 GB"
    aria2c -x 8 -s 8 -c https://slivar.s3.amazonaws.com/gnomad.hg38.genomes.v3.fix.zip
    print_status "gnomAD database downloaded"
else
    print_warning "gnomAD database already exists, skipping"
fi

# =================================================================
# 7. Exomiserのダウンロード
# =================================================================
print_header "Downloading Exomiser"

cd /variant_call/tools

if [ ! -d "exomiser-cli-13.1.0" ]; then
    show_file_size "https://data.monarchinitiative.org/exomiser/13.1.0/exomiser-cli-13.1.0-distribution.zip" \
                   "Exomiser CLI 13.1.0" "100 MB"
    aria2c -x 8 -s 8 -c https://data.monarchinitiative.org/exomiser/13.1.0/exomiser-cli-13.1.0-distribution.zip
    unzip -q exomiser-cli-13.1.0-distribution.zip
    print_status "Exomiser downloaded and extracted"
else
    print_warning "Exomiser already exists, skipping"
fi

# Exomiser データベース (オプション: サイズが大きいため、必要に応じてコメントアウトを解除)
# cd /variant_call/materials
# print_warning "Exomiser databases are very large (>40 GB each). Uncomment if needed."
# if [ ! -f "2209_hg38.zip" ]; then
#     show_file_size "https://data.monarchinitiative.org/exomiser/data/2209_hg38.zip" \
#                    "Exomiser hg38 database" "38.6 GB compressed, 60 GB uncompressed"
#     aria2c -x 16 -s 16 -c https://data.monarchinitiative.org/exomiser/data/2209_hg38.zip
#     unzip -q 2209_hg38.zip
#     print_status "Exomiser hg38 database downloaded"
# fi
#
# if [ ! -f "2209_phenotype.zip" ]; then
#     show_file_size "https://data.monarchinitiative.org/exomiser/data/2209_phenotype.zip" \
#                    "Exomiser phenotype database" "2.6 GB compressed, 13 GB uncompressed"
#     aria2c -x 8 -s 8 -c https://data.monarchinitiative.org/exomiser/data/2209_phenotype.zip
#     unzip -q 2209_phenotype.zip
#     print_status "Exomiser phenotype database downloaded"
# fi

# =================================================================
# ダウンロード完了サマリー
# =================================================================
print_header "Download Summary"

echo -e "${GREEN}All downloads completed successfully!${NC}\n"

echo "Downloaded files summary:"
echo "=========================="
echo "Tools:"
echo "  - Samtools 1.16.1"
echo "  - GATK 4.3.0.0"
echo "  - Strelka 2.9.10"
echo "  - SnpEff/SnpSift"
echo "  - vcfanno"
echo "  - slivar"
echo "  - Exomiser 13.1.0"
echo ""
echo "Data files:"
echo "  - FASTQ files (6 files, ~60 GB total)"
echo "  - JG2.1.0 Reference Genome (~3 GB)"
echo "  - JG2.1.0 Resource Bundle (~500 MB)"
echo "  - SnpEff JG2.1.0 Database (~500 MB)"
echo "  - gnomAD v3 Database (~2 GB)"
echo ""
echo "Total estimated disk usage: ~70-80 GB (excluding Exomiser databases)"
echo ""
print_status "You can now run 'bash run_all.sh' to execute the analysis pipeline"
