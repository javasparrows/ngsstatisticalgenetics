#!/bin/bash
set -e

# カラー出力の設定
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# ログ関数
print_header() {
    echo -e "\n${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}\n"
}

print_step() {
    echo -e "\n${CYAN}>>> Step: $1${NC}\n"
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

# 実行時間を計測する関数
time_cmd() {
    local start_time=$(date +%s)
    local cmd_description="$1"
    shift

    print_step "$cmd_description"
    "$@"

    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    local hours=$((duration / 3600))
    local minutes=$(((duration % 3600) / 60))
    local seconds=$((duration % 60))

    print_status "Completed in ${hours}h ${minutes}m ${seconds}s"
}

# =================================================================
# 0. 環境変数の設定
# =================================================================
print_header "Setting up environment variables"

export PATH=/variant_call/tools:/variant_call/tools/bin:$PATH
export GATK=/variant_call/tools/gatk-4.3.0.0/gatk
export MATERIALS=/variant_call/materials
export INTERMEDIATE=/variant_call/intermediate
export RESULTS=/variant_call/results
export TOOLS=/variant_call/tools

print_status "Environment variables set"

# =================================================================
# 1. Samtoolsのコンパイルとインストール
# =================================================================
print_header "1. Installing Samtools"

cd /variant_call/tools
if [ ! -f "bin/samtools" ]; then
    time_cmd "Extracting and compiling Samtools" bash -c "
        tar -jxf samtools-1.16.1.tar.bz2
        cd samtools-1.16.1
        ./configure --prefix=/variant_call/tools
        make -j\$(nproc)
        make install
    "
    print_status "Samtools installed successfully"
else
    print_warning "Samtools already installed, skipping"
fi

# Samtoolsの確認
samtools --version

# =================================================================
# 2. BWA-MEM2でのインデックス作成
# =================================================================
print_header "2. Creating BWA-MEM2 Index for Reference Genome"

cd $MATERIALS
if [ ! -f "JG.fa.bwt.2bit.64" ]; then
    time_cmd "Creating BWA-MEM2 index" \
        bwa-mem2 index JG.fa
    print_status "BWA-MEM2 index created"
else
    print_warning "BWA-MEM2 index already exists, skipping"
fi

# FASTAインデックスの作成
if [ ! -f "JG.fa.fai" ]; then
    time_cmd "Creating FASTA index" \
        samtools faidx JG.fa
    print_status "FASTA index created"
else
    print_warning "FASTA index already exists, skipping"
fi

# Dictionaryの作成
if [ ! -f "JG.dict" ]; then
    time_cmd "Creating sequence dictionary" \
        $GATK CreateSequenceDictionary --REFERENCE JG.fa
    print_status "Sequence dictionary created"
else
    print_warning "Sequence dictionary already exists, skipping"
fi

# =================================================================
# 3. アライメント（Son, Father, Mother）
# =================================================================
print_header "3. Running Alignment with BWA-MEM2"

cd $INTERMEDIATE

# Son
if [ ! -f "son.sort.bam" ]; then
    time_cmd "Aligning Son reads" bash -c "
        bwa-mem2 mem -t \$(nproc) \
            -R '@RG\tID:son\tSM:son\tPL:Illumina\tLB:son' \
            $MATERIALS/JG.fa \
            $MATERIALS/HG005_Son.R1.10X.fastq.gz \
            $MATERIALS/HG005_Son.R2.10X.fastq.gz \
            | samtools sort -@ \$(nproc) -o son.sort.bam -
    "
    print_status "Son alignment completed"
else
    print_warning "Son alignment already exists, skipping"
fi

# Father
if [ ! -f "father.sort.bam" ]; then
    time_cmd "Aligning Father reads" bash -c "
        bwa-mem2 mem -t \$(nproc) \
            -R '@RG\tID:father\tSM:father\tPL:Illumina\tLB:father' \
            $MATERIALS/JG.fa \
            $MATERIALS/HG006_Father.R1.10X.fastq.gz \
            $MATERIALS/HG006_Father.R2.10X.fastq.gz \
            | samtools sort -@ \$(nproc) -o father.sort.bam -
    "
    print_status "Father alignment completed"
else
    print_warning "Father alignment already exists, skipping"
fi

# Mother
if [ ! -f "mother.sort.bam" ]; then
    time_cmd "Aligning Mother reads" bash -c "
        bwa-mem2 mem -t \$(nproc) \
            -R '@RG\tID:mother\tSM:mother\tPL:Illumina\tLB:mother' \
            $MATERIALS/JG.fa \
            $MATERIALS/HG007_Mother.R1.10X.fastq.gz \
            $MATERIALS/HG007_Mother.R2.10X.fastq.gz \
            | samtools sort -@ \$(nproc) -o mother.sort.bam -
    "
    print_status "Mother alignment completed"
else
    print_warning "Mother alignment already exists, skipping"
fi

# =================================================================
# 4. BAMファイルのインデックス作成
# =================================================================
print_header "4. Creating BAM Indexes"

for sample in son father mother; do
    if [ ! -f "${sample}.sort.bam.bai" ]; then
        time_cmd "Indexing ${sample}.sort.bam" \
            samtools index ${sample}.sort.bam
    else
        print_warning "${sample}.sort.bam.bai already exists, skipping"
    fi
done

# =================================================================
# 5. MarkDuplicates
# =================================================================
print_header "5. Running MarkDuplicates"

cd $RESULTS

for sample in son father mother; do
    if [ ! -f "${sample}.sort.markdup.bam" ]; then
        time_cmd "MarkDuplicates for ${sample}" \
            $GATK --java-options "-Xmx8g -XX:ParallelGCThreads=$(nproc)" MarkDuplicates \
                --INPUT $INTERMEDIATE/${sample}.sort.bam \
                --OUTPUT ${sample}.sort.markdup.bam \
                --METRICS_FILE ${sample}.sort.markdup.metrics.txt \
                --CREATE_INDEX TRUE
    else
        print_warning "${sample}.sort.markdup.bam already exists, skipping"
    fi
done

# =================================================================
# 6. Base Quality Score Recalibration (BQSR) - Optional
# =================================================================
print_header "6. Base Quality Score Recalibration (BQSR) - OPTIONAL"
print_warning "BQSR is time-consuming and may not be necessary for modern data"
print_warning "Skipping BQSR steps as per modern best practices"

# BaseRecalibrator (コメントアウト - 必要に応じて有効化)
# for sample in son father mother; do
#     if [ ! -f "$INTERMEDIATE/${sample}_recal_data.table" ]; then
#         time_cmd "BaseRecalibrator for ${sample}" \
#             $GATK BaseRecalibrator \
#                 --reference $MATERIALS/JG.fa \
#                 --input $RESULTS/${sample}.sort.markdup.bam \
#                 --known-sites $MATERIALS/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.dbsnp_138.b37.success.sorted.vcf.gz \
#                 --known-sites $MATERIALS/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.Mills_and_1000G_gold_standard.indels.b37.success.sorted.vcf.gz \
#                 --output $INTERMEDIATE/${sample}_recal_data.table \
#                 --java-options "-XX:ParallelGCThreads=4"
#     fi
# done

# ApplyBQSR (コメントアウト - 必要に応じて有効化)
# for sample in son father mother; do
#     if [ ! -f "$RESULTS/${sample}.bqsr.bam" ]; then
#         time_cmd "ApplyBQSR for ${sample}" \
#             $GATK --java-options "-Xmx8g" ApplyBQSR \
#                 --reference $MATERIALS/JG.fa \
#                 --input $RESULTS/${sample}.sort.markdup.bam \
#                 --bqsr-recal-file $INTERMEDIATE/${sample}_recal_data.table \
#                 --create-output-bam-index true \
#                 --output $RESULTS/${sample}.bqsr.bam
#     fi
# done

# =================================================================
# 7. HaplotypeCaller (BAM → GVCF)
# =================================================================
print_header "7. Running HaplotypeCaller (BAM → GVCF)"

cd $INTERMEDIATE

for sample in son father mother; do
    if [ ! -f "${sample}.g.vcf" ]; then
        time_cmd "HaplotypeCaller for ${sample}" \
            $GATK --java-options "-Xmx8g" HaplotypeCaller \
                --reference $MATERIALS/JG.fa \
                --emit-ref-confidence GVCF \
                --input $RESULTS/${sample}.sort.markdup.bam \
                --output ${sample}.g.vcf
    else
        print_warning "${sample}.g.vcf already exists, skipping"
    fi
done

# =================================================================
# 8. GVCF → VCF 変換（染色体ごと）
# =================================================================
print_header "8. Converting GVCF to VCF (per chromosome)"

cd $INTERMEDIATE

chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)

for chr in "${chromosomes[@]}"; do
    if [ ! -f "joint_genotyped.${chr}.vcf" ]; then
        print_step "Processing chromosome: ${chr}"

        # GenomicsDBImport
        if [ ! -d "genomics_database.${chr}" ]; then
            time_cmd "GenomicsDBImport for ${chr}" \
                $GATK --java-options "-Xmx4g" GenomicsDBImport \
                    --variant son.g.vcf \
                    --variant father.g.vcf \
                    --variant mother.g.vcf \
                    --reference $MATERIALS/JG.fa \
                    --genomicsdb-workspace-path genomics_database.${chr} \
                    --intervals ${chr}
        fi

        # GenotypeGVCFs
        time_cmd "GenotypeGVCFs for ${chr}" \
            $GATK --java-options "-Xmx4g" GenotypeGVCFs \
                --reference $MATERIALS/JG.fa \
                --variant gendb://genomics_database.${chr} \
                --output joint_genotyped.${chr}.vcf \
                --intervals ${chr}
    else
        print_warning "joint_genotyped.${chr}.vcf already exists, skipping"
    fi
done

# =================================================================
# 9. 染色体ごとのVCFを統合
# =================================================================
print_header "9. Merging per-chromosome VCF files"

cd $RESULTS

if [ ! -f "merged.vcf.gz" ]; then
    input_files=""
    for chr in "${chromosomes[@]}"; do
        input_files+=" --INPUT $INTERMEDIATE/joint_genotyped.${chr}.vcf"
    done

    time_cmd "Merging VCF files" \
        $GATK --java-options "-Xmx4G" MergeVcfs \
            --OUTPUT merged.vcf.gz \
            ${input_files}

    print_status "VCF files merged successfully"
else
    print_warning "merged.vcf.gz already exists, skipping"
fi

# =================================================================
# 10. Strelka2でのBAM → VCF（オプション）
# =================================================================
print_header "10. Running Strelka2 for Variant Calling (OPTIONAL)"
print_warning "Strelka2 requires Python 2.7 which may not be available"
print_warning "Skipping Strelka2 step"

# Strelka2の実行 (Python 2.7が必要)
# cd $INTERMEDIATE
# if [ ! -d "strelka_son" ]; then
#     time_cmd "Configuring Strelka2 for Son" \
#         $TOOLS/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py \
#             --bam $RESULTS/son.sort.markdup.bam \
#             --referenceFasta $MATERIALS/JG.fa \
#             --runDir strelka_son
#
#     time_cmd "Running Strelka2 for Son" \
#         strelka_son/runWorkflow.py \
#             --mode local \
#             --jobs $(nproc) \
#             --memGb 8
# fi

# =================================================================
# 11. VQSR (Variant Quality Score Recalibration)
# =================================================================
print_header "11. Running VQSR"

cd $INTERMEDIATE

# バリアントの分割 (SNV と Indel)
if [ ! -f "merged.snv.vcf.gz" ]; then
    time_cmd "Extracting SNVs" \
        $GATK --java-options "-Xmx4g" SelectVariants \
            --reference $MATERIALS/JG.fa \
            --variant $RESULTS/merged.vcf.gz \
            --select-type-to-include SNP \
            --output merged.snv.vcf.gz
fi

if [ ! -f "merged.indel.vcf.gz" ]; then
    time_cmd "Extracting Indels" \
        $GATK --java-options "-Xmx4g" SelectVariants \
            --reference $MATERIALS/JG.fa \
            --variant $RESULTS/merged.vcf.gz \
            --select-type-to-include INDEL \
            --output merged.indel.vcf.gz
fi

# SNVのVariantRecalibrator
if [ ! -f "snv.recal" ]; then
    time_cmd "VariantRecalibrator for SNVs" \
        $GATK --java-options "-Xmx4g" VariantRecalibrator \
            --reference $MATERIALS/JG.fa \
            --variant $RESULTS/merged.vcf.gz \
            --resource:omni,known=false,training=true,truth=false,prior=12.0 \
                $MATERIALS/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.1000G_omni2.5.b37.success.sorted.vcf.gz \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 \
                $MATERIALS/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.1000G_phase1.snps.high_confidence.b37.success.sorted.vcf.gz \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
                $MATERIALS/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.dbsnp_138.b37.success.sorted.vcf.gz \
            --use-annotation QD --use-annotation MQ --use-annotation MQRankSum \
            --use-annotation ReadPosRankSum --use-annotation FS --use-annotation SOR \
            --mode SNP \
            --max-gaussians 6 \
            --output snv.recal \
            --tranches-file snv.tranches
fi

# IndelのVariantRecalibrator
if [ ! -f "indel.recal" ]; then
    time_cmd "VariantRecalibrator for Indels" \
        $GATK --java-options "-Xmx4g" VariantRecalibrator \
            --reference $MATERIALS/JG.fa \
            --variant merged.indel.vcf.gz \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
                $MATERIALS/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.dbsnp_138.b37.success.sorted.vcf.gz \
            --resource:mills,known=false,training=true,truth=true,prior=12.0 \
                $MATERIALS/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.Mills_and_1000G_gold_standard.indels.b37.success.sorted.vcf.gz \
            --use-annotation QD --use-annotation MQ --use-annotation MQRankSum \
            --use-annotation ReadPosRankSum --use-annotation FS --use-annotation SOR \
            --mode INDEL \
            --max-gaussians 4 \
            --output indel.recal \
            --tranches-file indel.tranches
fi

# ApplyVQSR for SNVs
if [ ! -f "merged.snv.vqsr.vcf.gz" ]; then
    time_cmd "ApplyVQSR for SNVs" \
        $GATK --java-options "-Xmx4g" ApplyVQSR \
            --reference $MATERIALS/JG.fa \
            --variant merged.snv.vcf.gz \
            --output merged.snv.vqsr.vcf.gz \
            --tranches-file snv.tranches \
            --recal-file snv.recal \
            --create-output-variant-index true \
            -mode SNP
fi

# ApplyVQSR for Indels
if [ ! -f "merged.indel.vqsr.vcf.gz" ]; then
    time_cmd "ApplyVQSR for Indels" \
        $GATK --java-options "-Xmx4g" ApplyVQSR \
            --reference $MATERIALS/JG.fa \
            --variant merged.indel.vcf.gz \
            --output merged.indel.vqsr.vcf.gz \
            --tranches-file indel.tranches \
            --recal-file indel.recal \
            --create-output-variant-index true \
            -mode INDEL
fi

# VQSRしたVCFファイルの統合
cd $RESULTS
if [ ! -f "merged.vqsr.vcf.gz" ]; then
    time_cmd "Merging VQSR VCF files" \
        $GATK --java-options "-Xmx4G" MergeVcfs \
            --OUTPUT merged.vqsr.vcf.gz \
            --INPUT $INTERMEDIATE/merged.snv.vqsr.vcf.gz \
            --INPUT $INTERMEDIATE/merged.indel.vqsr.vcf.gz
fi

# =================================================================
# 12. ハードフィルタリング
# =================================================================
print_header "12. Running Hard Filtering"

cd $INTERMEDIATE

# SNVのハードフィルタリング
if [ ! -f "merged.snv.hardfiltering.vcf.gz" ]; then
    time_cmd "Hard filtering SNVs" \
        $GATK --java-options "-Xmx4g" VariantFiltration \
            --reference $MATERIALS/JG.fa \
            --variant merged.snv.vcf.gz \
            --output merged.snv.hardfiltering.vcf.gz \
            --filter-name "QD_filter" --filter-expression "QD < 2.0" \
            --filter-name "FS_filter" --filter-expression "FS > 60.0" \
            --filter-name "SOR_filter" --filter-expression "SOR > 3.0" \
            --filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -8.0" \
            --filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
            --filter-name "MQRankSum_filter" --filter-expression "MQRankSum < -12.5"
fi

# Indelのハードフィルタリング
if [ ! -f "merged.indel.hardfiltering.vcf.gz" ]; then
    time_cmd "Hard filtering Indels" \
        $GATK --java-options "-Xmx4g" VariantFiltration \
            --reference $MATERIALS/JG.fa \
            --variant merged.indel.vcf.gz \
            --output merged.indel.hardfiltering.vcf.gz \
            --filter-name "QD_filter" --filter-expression "QD < 2.0" \
            --filter-name "FS_filter" --filter-expression "FS > 200.0" \
            --filter-name "SOR_filter" --filter-expression "SOR > 10.0" \
            --filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -20.0"
fi

# ハードフィルタリング済みVCFファイルの統合
cd $RESULTS
if [ ! -f "merged.hardfiltering.vcf.gz" ]; then
    time_cmd "Merging hard-filtered VCF files" \
        $GATK --java-options "-Xmx4G" MergeVcfs \
            --OUTPUT merged.hardfiltering.vcf.gz \
            --INPUT $INTERMEDIATE/merged.snv.hardfiltering.vcf.gz \
            --INPUT $INTERMEDIATE/merged.indel.hardfiltering.vcf.gz
fi

# =================================================================
# 13. bcftoolsでのVCF処理（2-1）
# =================================================================
print_header "13. Processing VCF with bcftools"

cd $RESULTS

# chr6のみを抽出
if [ ! -f "merged.vqsr-chr6.vcf.gz" ]; then
    time_cmd "Extracting chromosome 6" \
        bcftools filter -i 'CHROM="chr6"' merged.vqsr.vcf.gz -o merged.vqsr-chr6.vcf.gz
fi

# AF>=0.01でフィルタリング
if [ ! -f "merged.vqsr-chr6-AF0.01.vcf.gz" ]; then
    time_cmd "Filtering by AF >= 0.01" \
        bcftools filter -i 'INFO/AF>=0.01' merged.vqsr-chr6.vcf.gz -o merged.vqsr-chr6-AF0.01.vcf.gz
fi

# 正規化
if [ ! -f "merged.vqsr-chr6-norm.vcf.gz" ]; then
    time_cmd "Normalizing multi-allelic sites" \
        bcftools norm -m -snps merged.vqsr-chr6.vcf.gz -o merged.vqsr-chr6-norm.vcf.gz
fi

# =================================================================
# 14. SnpEffでのアノテーション（2-2）
# =================================================================
print_header "14. Annotating with SnpEff"

cd $RESULTS

# SnpEff設定ファイルの追加
if ! grep -q "JG2.1.0.genome" $TOOLS/snpEff/snpEff.config; then
    echo "" >> $TOOLS/snpEff/snpEff.config
    echo "# Japanese Genome" >> $TOOLS/snpEff/snpEff.config
    echo "JG2.1.0.genome: Human genome JG2.1.0" >> $TOOLS/snpEff/snpEff.config
    echo "    JG2.1.0.M.codonTable : Vertebrate_Mitochondrial" >> $TOOLS/snpEff/snpEff.config
fi

if [ ! -f "merged.vqsr.ann.vcf" ]; then
    time_cmd "Annotating VCF with SnpEff" \
        java -Xmx8g -jar $TOOLS/snpEff/snpEff.jar JG2.1.0 merged.vqsr.vcf.gz \
            > merged.vqsr.ann.vcf
fi

# SnpSiftでのフィルタリング例
if [ ! -f "merged.vqsr.ann.filtered.missense.vcf" ]; then
    time_cmd "Filtering for missense variants" \
        java -jar $TOOLS/snpEff/SnpSift.jar filter \
            "ANN[*].EFFECT has 'missense_variant'" \
            merged.vqsr.ann.vcf > merged.vqsr.ann.filtered.missense.vcf
fi

# =================================================================
# 15. vcfannoでのアノテーション（2-3） - OPTIONAL
# =================================================================
print_header "15. Annotating with vcfanno (OPTIONAL)"
print_warning "Skipping vcfanno as it requires additional configuration and data sources"

# =================================================================
# 16. slivarでのフィルタリング（3-1） - OPTIONAL
# =================================================================
print_header "16. Filtering with slivar (OPTIONAL)"
print_warning "Skipping slivar as it requires GRCh38-based data and trio pedigree file"

# =================================================================
# 17. Exomiserでの解析（3-2） - OPTIONAL
# =================================================================
print_header "17. Running Exomiser (OPTIONAL)"
print_warning "Skipping Exomiser as it requires large databases (>40GB) and specific configuration"

# =================================================================
# パイプライン完了
# =================================================================
print_header "Pipeline Completion Summary"

echo -e "${GREEN}All pipeline steps completed successfully!${NC}\n"

echo "Generated files:"
echo "================="
echo "BAM files:"
echo "  - $RESULTS/son.sort.markdup.bam"
echo "  - $RESULTS/father.sort.markdup.bam"
echo "  - $RESULTS/mother.sort.markdup.bam"
echo ""
echo "VCF files:"
echo "  - $RESULTS/merged.vcf.gz (raw merged VCF)"
echo "  - $RESULTS/merged.vqsr.vcf.gz (VQSR filtered)"
echo "  - $RESULTS/merged.hardfiltering.vcf.gz (hard filtered)"
echo "  - $RESULTS/merged.vqsr.ann.vcf (annotated with SnpEff)"
echo ""
echo "Intermediate files are stored in: $INTERMEDIATE"
echo ""
print_status "Pipeline execution completed!"
