# Variant Calling Commands

## コード1: 環境設定とフォルダ構造の作成

```bash
# ホームフォルダへ移動する
cd

# ファイルやソフトを格納するフォルダを作成する
mkdir variant_call
cd variant_call
mkdir tools
mkdir materials
mkdir intermediate
mkdir results

# ソフトの実行のためにvariant_call内のtoolsにパスを通す．環境や必要性に応じて適当な記述が必要となる
export PATH=~/variant_call/tools:$PATH
```

## コード2: Javaバージョン確認

```bash
# Javaが存在する場合，Java のバージョンが表示される
java -version
```

## コード3: Pythonバージョン確認

```bash
# Pythonが存在する場合，Pythonのバージョンが表示される
python -–version
```

## コード4: Samtoolsのインストール

```bash
# 2023年1月段階での最新
cd ~/variant_call/tools
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2

# ダウンロードしてきたsamtoolsの解凍
tar -jxvf samtools-1.16.1.tar.bz2
cd samtools-1.16.1

# prefix指定では絶対パスのみ有効．
# ${PWD%/*}で~/variant_call/toolsの絶対パスを使用する
./configure --prefix=${PWD%/*}
make
make install

# インストール後，コマンドが認識されない場合は，パスが通っているかを確認する
samtools
# 正しくインストールされていれば，samtoolsの使用法が表示される
```

## コード5: BWA-MEM2のインストール

```bash
# GitHub経由でダウンロードしたソースコードをコンパイルする
# GitHubページにはprecompiled binariesをダウンロードし，インストールする方法も記載されている
brew install brewsci/bio/bwa-mem2
./bwa-mem2
# 正しくインストールされていれば，使用法が表示される
```

## コード6: GATKのインストール

```bash
# 2022年10月段階での最新
cd ~/variant_call/tools
wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip
unzip gatk-4.3.0.0.zip
~/variant_call/tools/gatk-4.3.0.0/gatk
# 正しくインストールされていれば，使用法が表示される。うまくできなかったらPythonの環境を当てはめればいい。
```

## コード7: GATK HaplotypeCallerの使用例

```bash
gatk --java-options "-Xmx4g" HaplotypeCaller \
--reference reference.fasta \
--input input.bam \
-O output.g.vcf.gz \
-ERC GVCF
```

## コード8: FASTQファイルのダウンロード

### 30xのFASTQファイルをダウンロード（各ファイル約30GB）

```bash
cd ~/variant_call/materials
wget https://zenodo.org/record/7306771/files/HG005_Son.R1.fastq.gz
wget https://zenodo.org/record/7307245/files/HG005_Son.R2.fastq.gz
wget https://zenodo.org/record/7307258/files/HG006_Father.R1.fastq.gz
wget https://zenodo.org/record/7307266/files/HG006_Father.R2.fastq.gz
wget https://zenodo.org/record/7307277/files/HG007_Mother.R1.fastq.gz
wget https://zenodo.org/record/7307283/files/HG007_Mother.R2.fastq.gz
```

### aria2cを使用した高速ダウンロード（30x）

```bash
aria2c -x 8 -s 8 -c https://zenodo.org/record/7306771/files/HG005_Son.R1.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7307245/files/HG005_Son.R2.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7307258/files/HG006_Father.R1.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7307266/files/HG006_Father.R2.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7307277/files/HG007_Mother.R1.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7307283/files/HG007_Mother.R2.fastq.gz
```

### 10xのFASTQファイルをダウンロード（各ファイル約10GB）

```bash
cd ~/variant_call/materials
wget https://zenodo.org/record/7310196/files/HG005_Son.R1.10X.fastq.gz
wget https://zenodo.org/record/7310196/files/HG005_Son.R2.10X.fastq.gz
wget https://zenodo.org/record/7310202/files/HG006_Father.R1.10X.fastq.gz
wget https://zenodo.org/record/7310202/files/HG006_Father.R2.10X.fastq.gz
wget https://zenodo.org/record/7310233/files/HG007_Mother.R1.10X.fastq.gz
wget https://zenodo.org/record/7310233/files/HG007_Mother.R2.10X.fastq.gz
```

### aria2cを使用した高速ダウンロード（10x）

```bash
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310196/files/HG005_Son.R1.10X.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310196/files/HG005_Son.R2.10X.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310202/files/HG006_Father.R1.10X.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310202/files/HG006_Father.R2.10X.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310233/files/HG007_Mother.R1.10X.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310233/files/HG007_Mother.R2.10X.fastq.gz
```

## コード9: JG2.1.0参照ゲノムとリソースファイルのダウンロード

```bash
# JG2.1.0のFASTAファイルをダウンロード
cd ~/variant_call/materials
wget https://jmorp.megabank.tohoku.ac.jp/datasets/tommo-jg2.1.0-20211208/files/jg2.1.0.fa.gz
gzip -cd jg2.1.0.fa.gz> JG.fa && rm jg2.1.0.fa.gz

# BQSRおよびVQSR用に必要となるGATK Resource BundleのJG2.1.0対応版をダウンロード
cd ~/variant_call/materials
wget https://jmorp.megabank.tohoku.ac.jp/datasets/tommo-jg2.1.0-20211208/files/JG2.1.0-ResourceBundle-from-b37.zip
unzip JG2.1.0-ResourceBundle-from-b37.zip && rm JG2.1.0-ResourceBundle-fromb37.zip
```

## コード10: GRCh37・GRCh38参照ゲノムのダウンロード

```bash
# GRCh37(hs37d5)のFASTAファイルのダウンロード
aria2c -x 8 -s 8 -c https://ilmn-dragen-giab-samples.s3.amazonaws.com/FASTA/hs37d5.fa

# GRCh38のFASTAファイルのダウンロード
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz
gzip -d hg38.analysisSet.fa.gz
```

## コード11: GATK Resource Bundleのダウンロード

```bash
# GRCh37／GRCh38のGATK Resource Bundleは現在Google Cloud Storageで管理されており，ファイルをダウンロードするにはgsutilが必要
# GRCh37(h19)のGATK Resource Bundle のダウンロード
gsutil -o 'GSUtil:parallel_thread_count=24' \
       -o 'GSUtil:sliced_object_download_max_components=8' \
       -m cp -r \
       "gs://gcp-public-data--broad-references/hg38" \
       ~/variant_call/materials

# GRCh38のGATK Resource Bundleのダウンロード
gsutil -m cp -r "gs://gcp-public-data--broad-references/hg19" ~/variant_call/materials
```

## コード12: アライメント処理

### インデックスファイルの作成

```bash
# FASTAファイルのインデックスの作成
# bwa-mem2 indexでは多くのRAMとストレージを必要とする．https://doi.org/10.5281/zenodo.7567194の全ダウンロードをmaterials内に配置することで実行結果の代用にできる
bwa-mem2 index ~/variant_call/materials/JG.fa
# M4 MacBook Pro 128GBで10分かかった

# FASTAファイルのfasta.faiの作成（関連必要ファイル）
samtools faidx ~/variant_call/materials/JG.fa

# FASTAファイルのdictionaryの作成（関連必要ファイル）
~/variant_call/tools/gatk-4.3.0.0/gatk CreateSequenceDictionary --REFERENCE ~/variant_call/materials/JG.fa
```

### アライメントの実行

```bash
# Son（息子）のアライメント
bwa-mem2 mem -R "@RG\tID:son\tSM:son\tPL:Illumina\tLB:son" ~/variant_call/materials/JG.fa ~/variant_call/materials/HG005_Son.R1.10X.fastq.gz ~/variant_call/materials/HG005_Son.R2.10X.fastq.gz >~/variant_call/intermediate/son.sam
# 完了(2025/09/04)

# Father（父親）のアライメント
bwa-mem2 mem -R "@RG\tID:father\tSM:father\tPL:Illumina\tLB:father" ~/variant_call/materials/JG.fa ~/variant_call/materials/HG006_Father.R1.10X.fastq.gz ~/variant_call/materials/HG006_Father.R2.10X.fastq.gz >~/variant_call/intermediate/father.sam
# 完了(2025/09/12) M4 MacBook Pro 128GBで15時間かかる

## Father（父親）のアライメント（スレッド数を指定して高速化）
bwa-mem2 mem -t 8 -R "@RG\tID:father\tSM:father\tPL:Illumina\tLB:father" \
  ~/variant_call/materials/JG.fa \
  ~/variant_call/materials/HG006_Father.R1.10X.fastq.gz \
  ~/variant_call/materials/HG006_Father.R2.10X.fastq.gz \
  | samtools sort -@ 8 -o ~/variant_call/intermediate/father.sort.bam -

## nohupでバックグラウンド実行
nohup sh -c '
echo "Start time: $(date)"
time bwa-mem2 mem -t 8 -R "@RG\tID:father\tSM:father\tPL:Illumina\tLB:father" \
  ~/variant_call/materials/JG.fa \
  ~/variant_call/materials/HG006_Father.R1.10X.fastq.gz \
  ~/variant_call/materials/HG006_Father.R2.10X.fastq.gz \
  | samtools sort -@ 8 -o ~/variant_call/intermediate/father.sort.bam -
echo "End time: $(date)"
' > ~/variant_call/logs/father_alignment.log 2>&1 &

# Mother（母親）のアライメント
bwa-mem2 mem -R "@RG\tID:mother\tSM:mother\tPL:Illumina\tLB:mother" ~/variant_call/materials/JG.fa ~/variant_call/materials/HG007_Mother.R1.10X.fastq.gz ~/variant_call/materials/HG007_Mother.R2.10X.fastq.gz >~/variant_call/intermediate/mother.sam
# 完了(2025/09/13) M4 MacBook Pro 128GBで15時間かかる

## Mother（母親）のアライメント（スレッド数を指定して高速化）
bwa-mem2 mem -t 8 -R "@RG\tID:mother\tSM:mother\tPL:Illumina\tLB:mother" \
  ~/variant_call/materials/JG.fa \
  ~/variant_call/materials/HG006_Mother.R1.10X.fastq.gz \
  ~/variant_call/materials/HG006_Mother.R2.10X.fastq.gz \
  | samtools sort -@ 8 -o ~/variant_call/intermediate/mother.sort.bam -

## nohupでバックグラウンド実行
nohup sh -c '
echo "Start time: $(date)"
time bwa-mem2 mem -t 8 -R "@RG\tID:mother\tSM:mother\tPL:Illumina\tLB:mother" \
  ~/variant_call/materials/JG.fa \
  ~/variant_call/materials/HG007_Mother.R1.10X.fastq.gz \
  ~/variant_call/materials/HG007_Mother.R2.10X.fastq.gz \
  | samtools sort -@ 8 -o ~/variant_call/intermediate/mother.sort.bam -
echo "End time: $(date)"
' > ~/variant_call/logs/mother_alignment.log 2>&1 &
```

### SAMファイルをBAMファイルへ変換

```bash
# 結果確認のため，headコマンドでSAMファイルの冒頭行を表示する
head -5 ~/variant_call/intermediate/son.sam
head -5 ~/variant_call/intermediate/father.sam
head -5 ~/variant_call/intermediate/mother.sam

# ファイル不足等の際には，以下のコマンドによるBAMファイルの完成後，SAMファイルを削除してもよい
samtools view -bS ~/variant_call/intermediate/son.sam > ~/variant_call/intermediate/son.bam
# 15.6分後に正常終了、以下の標準出力
# 38
samtools view -bS ~/variant_call/intermediate/father.sam > ~/variant_call/intermediate/father.bam
# 5分後に以下エラーで終了
# [E::aux_parse] unrecognized type '\t'
# [W::sam_read1_sam] Parse error at line 73961479
# samtools view: error reading file "/Users/yukik/variant_call/intermediate/father.sam"

# 9/22 801秒で完了。8スレッド並列でやり直したら正常に通った
samtools view -bS ~/variant_call/intermediate/mother.sam > ~/variant_call/intermediate/mother.bam
# 14分後に以下エラーで終了
# 14
# [E::aux_parse] unrecognized type 'S'
# [W::sam_read1_sam] Parse error at line 184200616
# samtools view: error reading file "/Users/yukik/variant_call/intermediate/mother.sam"
```

## コード13: BAMファイルの処理

### BAMファイルのソート

```bash
# ここは-t 8のコマンドでsortできていれば実行しなくて良い
samtools sort -o ~/variant_call/intermediate/son.sort.bam ~/variant_call/intermediate/son.bam
samtools sort -o ~/variant_call/intermediate/father.sort.bam ~/variant_call/intermediate/father.bam
samtools sort -o ~/variant_call/intermediate/mother.sort.bam ~/variant_call/intermediate/mother.bam
```

### ソートしたBAMファイルのインデックス作成

```bash
samtools index ~/variant_call/intermediate/son.sort.bam
# ↑68秒 2025/09/24
samtools index ~/variant_call/intermediate/father.sort.bam
# ↑69秒 2025/09/24
samtools index ~/variant_call/intermediate/mother.sort.bam
# ↑70秒 2025/09/24
```

### 重複したリードを除去するMarkDuplicates (重複箇所にマーキングしているだけ)

```bash
# Son（息子）
~/variant_call/tools/gatk-4.3.0.0/gatk MarkDuplicates \
--INPUT ~/variant_call/intermediate/son.sort.bam \
--OUTPUT ~/variant_call/results/son.sort.markdup.bam \
--METRICS_FILE ~/variant_call/results/son.sort.markdup.metrics.txt \
--CREATE_INDEX TRUE
# 17分経過しても終わらないので以下の並列処理をする 2025/09/24

~/variant_call/tools/gatk-4.3.0.0/gatk MarkDuplicatesSpark \
--input ~/variant_call/intermediate/son.sort.bam \
--output ~/variant_call/results/son.sort.markdup.bam \
--metrics-file ~/variant_call/results/son.sort.markdup.metrics.txt \
--create-output-bam-index true \
--spark-runner LOCAL \
--spark-master local[8]
# ↑Java 21のモジュールシステムがSpark 2.4.5（GATK 4.3.0.0に内蔵）と互換性がなくエラーが出た。/usr/libexec/java_home -V
# ↑Java 11が必要でMacだと実行不可

~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx8g -XX:ParallelGCThreads=8 -XX:ConcGCThreads=2" MarkDuplicates \
--INPUT ~/variant_call/intermediate/son.sort.bam \
--OUTPUT ~/variant_call/results/son.sort.markdup.bam \
--METRICS_FILE ~/variant_call/results/son.sort.markdup.metrics.txt \
--CREATE_INDEX TRUE
# ↑18分 2025/09/25

# Father（父親）
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx8g -XX:ParallelGCThreads=8 -XX:ConcGCThreads=2" MarkDuplicates \
--INPUT ~/variant_call/intermediate/father.sort.bam \
--OUTPUT ~/variant_call/results/father.sort.markdup.bam \
--METRICS_FILE ~/variant_call/results/father.sort.markdup.metrics.txt \
--CREATE_INDEX TRUE
# ↑並列、18分 2025/09/25

# Mother（母親）
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx8g -XX:ParallelGCThreads=8 -XX:ConcGCThreads=2" MarkDuplicates \
--INPUT ~/variant_call/intermediate/mother.sort.bam \
--OUTPUT ~/variant_call/results/mother.sort.markdup.bam \
--METRICS_FILE ~/variant_call/results/mother.sort.markdup.metrics.txt \
--CREATE_INDEX TRUE
# ↑並列、18分 2025/09/25
```

### 結果確認

```bash
# 結果確認のため，headコマンドでBAMファイルの冒頭行を表示する
# なお，header行が必要な場合はsamtools view -h hoge.bamとなる
samtools view ~/variant_call/results/son.sort.markdup.bam | head
samtools view ~/variant_call/results/father.sort.markdup.bam | head
samtools view ~/variant_call/results/mother.sort.markdup.bam | head
# ↑2025/09/25
```

## コード14: BQSR（Base Quality Score Recalibration）

### BaseRecalibrator

```bash
# Son（息子）
~/variant_call/tools/gatk-4.3.0.0/gatk BaseRecalibrator --reference ~/variant_call/materials/JG.fa --input ~/variant_call/results/son.sort.markdup.bam \
  --known-sites ~/variant_call/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.dbsnp_138.b37.success.sorted.vcf.gz \
  --known-sites ~/variant_call/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.Mills_and_1000G_gold_standard.indels.b37.success.sorted.vcf.gz \
  --output ~/variant_call/intermediate/son_recal_data.table
  # ↑ 49分, シングルスレッド, 2025/09/25

~/variant_call/tools/gatk-4.3.0.0/gatk BaseRecalibrator \
  --reference ~/variant_call/materials/JG.fa \
  --input ~/variant_call/results/son.sort.markdup.bam \
  --known-sites ~/variant_call/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.dbsnp_138.b37.success.sorted.vcf.gz \
  --known-sites ~/variant_call/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.Mills_and_1000G_gold_standard.indels.b37.success.sorted.vcf.gz \
  --output ~/variant_call/intermediate/son_recal_data.table \
  --java-options "-XX:ParallelGCThreads=4"

# ↑ 42分, マルチスレッド, 2025/09/27

# Father（父親）
~/variant_call/tools/gatk-4.3.0.0/gatk BaseRecalibrator --reference ~/variant_call/materials/JG.fa --input ~/variant_call/results/father.sort.markdup.bam \
  --known-sites ~/variant_call/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.dbsnp_138.b37.success.sorted.vcf.gz \
  --known-sites ~/variant_call/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.Mills_and_1000G_gold_standard.indels.b37.success.sorted.vcf.gz \
  --output ~/variant_call/intermediate/father_recal_data.table
# ↑シングルスレッド

~/variant_call/tools/gatk-4.3.0.0/gatk BaseRecalibrator \
  --reference ~/variant_call/materials/JG.fa \
  --input ~/variant_call/results/father.sort.markdup.bam \
  --known-sites ~/variant_call/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.dbsnp_138.b37.success.sorted.vcf.gz \
  --known-sites ~/variant_call/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.Mills_and_1000G_gold_standard.indels.b37.success.sorted.vcf.gz \
  --output ~/variant_call/intermediate/father_recal_data.table \
  --java-options "-XX:ParallelGCThreads=4"

# ↑ 分, マルチスレッド, 2025/09/25

# Mother（母親）
~/variant_call/tools/gatk-4.3.0.0/gatk BaseRecalibrator \
  --reference ~/variant_call/materials/JG.fa \
  --input ~/variant_call/results/mother.sort.markdup.bam \
  --known-sites ~/variant_call/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.dbsnp_138.b37.success.sorted.vcf.gz \
  --known-sites ~/variant_call/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.Mills_and_1000G_gold_standard.indels.b37.success.sorted.vcf.gz \
  --output ~/variant_call/intermediate/mother_recal_data.table \
  --java-options "-XX:ParallelGCThreads=4"

# ↑時間かかるばかりで、最近はやる必要はないと言われてきている。高山先生
```

### ApplyBQSR

```bash
# Son（息子）
~/variant_call/tools/gatk-4.3.0.0/gatk ApplyBQSR --reference ~/variant_call/materials/JG.fa \
  --input ~/variant_call/results/son.sort.markdup.bam \
  --bqsr-recal-file ~/variant_call/intermediate/son_recal_data.table \
  --create-output-bam-index true \
  --output ~/variant_call/results/son.bqsr.bam

# ↑シングルスレッド

~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx8g" ApplyBQSR \
  --reference ~/variant_call/materials/JG.fa \
  --input ~/variant_call/results/son.sort.markdup.bam \
  --bqsr-recal-file ~/variant_call/intermediate/son_recal_data.table \
  --create-output-bam-index true \
  --output ~/variant_call/results/son.bqsr.bam

# ↑マルチスレッド

# Father(父親) - マルチスレッド
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx8g" ApplyBQSR \
  --reference ~/variant_call/materials/JG.fa \
  --input ~/variant_call/results/father.sort.markdup.bam \
  --bqsr-recal-file ~/variant_call/intermediate/father_recal_data.table \
  --create-output-bam-index true \
  --output ~/variant_call/results/father.bqsr.bam

# Mother(母親) - マルチスレッド
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx8g" ApplyBQSR \
  --reference ~/variant_call/materials/JG.fa \
  --input ~/variant_call/results/mother.sort.markdup.bam \
  --bqsr-recal-file ~/variant_call/intermediate/mother_recal_data.table \
  --create-output-bam-index true \
  --output ~/variant_call/results/mother.bqsr.bam
```

## コード15: バリアントコーリング（BAM→GVCF）

### HaplotypeCaller実行

```bash
# BAM→GVCFへのバリアントコールを実行（今回はBQSR不実行）
## GATKの代わりになる方法
1. https://github.com/Illumina/strelka → BAMから一気にVCFを作ってくれる
2. https://github.com/luntergroup/octopus
3. https://github.com/google/deepvariant → 高山先生の実験では偽陽性が多かった
# 400人くらいあると精度が良い

# Son（息子）
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx8g" HaplotypeCaller \
  --reference ~/variant_call/materials/JG.fa \
  --emit-ref-confidence GVCF \
  --input ~/variant_call/results/son.sort.markdup.bam \
  --output ~/variant_call/intermediate/son.g.vcf

# Father（父親）
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx8g" HaplotypeCaller \
  --reference ~/variant_call/materials/JG.fa \
  --emit-ref-confidence GVCF \
  --input ~/variant_call/results/father.sort.markdup.bam \
  --output ~/variant_call/intermediate/father.g.vcf

# Mother（母親）
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx8g" HaplotypeCaller \
  --reference ~/variant_call/materials/JG.fa \
  --emit-ref-confidence GVCF \
  --input ~/variant_call/results/mother.sort.markdup.bam \
  --output ~/variant_call/intermediate/mother.g.vcf
```

### GVCF→VCF変換（染色体ごとに分割実行）

```bash
# 通例GVCF→VCFは染色体ごとに分割して実行する. databaseにまとめて，ジェノタイピングを行っている.
# ここの部分はfor文であるため，シェルスクリプトファイル（.sh）にして実行するとよい．

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
do
  ~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx4g" GenomicsDBImport \
    --variant ~/variant_call/intermediate/son.g.vcf \
    --variant ~/variant_call/intermediate/father.g.vcf \
    --variant ~/variant_call/intermediate/mother.g.vcf \
    --reference ~/variant_call/materials/JG.fa \
    --genomicsdb-workspace-path ~/variant_call/intermediate/genomics_database.${chr} \
    --intervals ${chr}
  ~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx4g" GenotypeGVCFs \
    --reference ~/variant_call/materials/JG.fa \
    --variant gendb://~/variant_call/intermediate/genomics_database.${chr} \
    --output ~/variant_call/intermediate/joint_genotyped.${chr}.vcf \
    --intervals ${chr}
done
# ↑2025/11/06 45秒
```
### 染色体ごとのVCFを統合

```bash
# 染色体ごとのVCFを統合する
# ここの部分もfor文であるため，シェルスクリプトファイル（.sh）にして実行するとよい.
input_files=""
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
do
  input_files+=" --INPUT ~/variant_call/intermediate/joint_genotyped.${chr}.vcf"
done

~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx4G" MergeVcfs \
  --OUTPUT ~/variant_call/results/merged.vcf.gz \
  ${input_files}

# ↑2025/11/27 gvcf_vcf.sh, gvcf_vcf_merge.shを作成 60秒で完了
```

# Strelka2でのBAM → VCF
```bash
# Step 1: 解析の設定
~/variant_call/tools/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py \
  --bam ~/variant_call/results/son.sort.markdup.bam \
  --referenceFasta ~/variant_call/materials/JG.fa \
  --runDir ~/variant_call/intermediate/strelka_son

# Step 2  : 実行
~/variant_call/intermediate/strelka_son/runWorkflow.py \
  --mode local \
  --jobs 4 \
  --memGb 8

# ↑2025/11/27 Python2でないと動かないので方法を考え中
# ↑2025/11/27 export PATH="$HOME/.pyenv/versions/2.7.18/bin:$PATH" sh bam_vcf.sh → ダメ
# ↑2025/12/04 bash /Users/yukik/variant_call/run_strelka.sh  44分かかった
```

## コード16: VQSR（Variant Quality Score Recalibration）

### SNV（SNP）のVariantRecalibrator

```bash
# SNV（SNP）についてVariantRecalibratorでrecalibration modelを作成する
# 一般に1塩基の変異はSNVと呼び，その中でも1％以上の頻度で変異が認められた場合をSNP（一塩基多型）と呼ぶ
# HapMapのリソースを利用しない場合
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx4g" VariantRecalibrator \
    --reference ~/variant_call/materials/JG.fa \
    --variant ~/variant_call/results/merged.vcf.gz \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 \
  ~/variant_call/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.1000G_omni2.5.b37.success.sorted.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 \
  ~/variant_call/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.1000G_phase1.snps.high_confidence.b37.success.sorted.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
  ~/variant_call/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.dbsnp_138.b37.success.sorted.vcf.gz \
    --use-annotation QD --use-annotation MQ --use-annotation MQRankSum --use-annotation ReadPosRankSum --use-annotation FS --use-annotation SOR \
    --mode SNP \
    --max-gaussians 6 \
    --output ~/variant_call/intermediate/snv.recal \
    --tranches-file ~/variant_call/intermediate/snv.tranches

# ↑2025/12/15 11秒かかった sh ~/variant_call/tools/run_vqsr.sh
```


### IndelのVariantRecalibrator

```bash
# indelについてVariantRecalibratorでrecalibration modelを作成する
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx4g" VariantRecalibrator \
  --reference ~/variant_call/materials/JG.fa \
  --variant ~/variant_call/intermediate/merged.indel.vcf.gz \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ~/variant_call/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.dbsnp_138.b37.success.sorted.vcf.gz \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 ~/variant_call/materials/JG2.1.0-ResourceBundle-from-b37/human_g1k_v37_to_JG2.1.0.Mills_and_1000G_gold_standard.indels.b37.success.sorted.vcf.gz \
  --use-annotation QD --use-annotation MQ --use-annotation MQRankSum --use-annotation ReadPosRankSum --use-annotation FS --use-annotation SOR \
  --mode INDEL \
  --max-gaussians 4 \
  --output ~/variant_call/intermediate/indel.recal \
  --tranches-file ~/variant_call/intermediate/indel.tranches
```

### ApplyVQSR

```bash
# SNVとindelのそれぞれに構築したrecalibration modelをApplyVQSRで利用する
# SNV用ApplyVQSR
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx4g" ApplyVQSR \
  --reference ~/variant_call/materials/JG.fa \
  --variant ~/variant_call/intermediate/merged.snv.vcf.gz \
  --output ~/variant_call/intermediate/merged.snv.vqsr.vcf.gz \
  --tranches-file ~/variant_call/intermediate/snv.tranches \
  --recal-file ~/variant_call/intermediate/snv.recal \
  --create-output-variant-index true \
  -mode SNP

# Indel用ApplyVQSR
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx4g" ApplyVQSR \
  --reference ~/variant_call/materials/JG.fa \
  --variant ~/variant_call/intermediate/merged.indel.vcf.gz \
  --output ~/variant_call/intermediate/merged.indel.vqsr.vcf.gz \
  --tranches-file ~/variant_call/intermediate/indel.tranches \
  --recal-file ~/variant_call/intermediate/indel.recal \
  --create-output-variant-index true \
  -mode INDEL
```

### VCFファイルの統合

```bash
# SNVとindelのそれぞれでフィルタリングしたVCFをMergeVcfsで1つに統合する
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx4G" MergeVcfs \
  --OUTPUT ~/variant_call/results/merged.vqsr.vcf.gz \
  --INPUT ~/variant_call/intermediate/merged.snv.vqsr.vcf.gz \
  --INPUT ~/variant_call/intermediate/merged.indel.vqsr.vcf.gz
```

## コード17: ハードフィルタリング

### バリアントの分割

```bash
# 検出されたバリアントをフィルタリングするためにSNVとindelに分割する
# SNVの抽出
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx4g" SelectVariants \
  --reference ~/variant_call/materials/JG.fa \
  --variant ~/variant_call/results/merged.vcf.gz \
  --select-type-to-include SNP \
  --output ~/variant_call/intermediate/merged.snv.vcf.gz

# Indelの抽出
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx4g" SelectVariants \
  --reference ~/variant_call/materials/JG.fa \
  --variant ~/variant_call/results/merged.vcf.gz \
  --select-type-to-include INDEL \
  --output ~/variant_call/intermediate/merged.indel.vcf.gz
```

### VariantFiltrationによるフィルタリング

```bash
# VCFファイルの統計量について，事前に設定しておいた値をもとにVariantFiltrationでフィルタリングする
# SNVの設定
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx4g" VariantFiltration \
  --reference ~/variant_call/materials/JG.fa \
  --variant ~/variant_call/intermediate/merged.snv.vcf.gz \
  --output ~/variant_call/intermediate/merged.snv.hardfiltering.vcf.gz \
  --filter-expression "QD < 2.0" \
  --filter-expression "FS > 60.0" \
  --filter-expression "SOR > 3.0" \
  --filter-expression "ReadPosRankSum < -8.0" \
  --filter-expression "MQ < 40.0" \
  --filter-expression "MQRankSum < -12.5"

# indelの設定
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx4g" VariantFiltration \
  --reference ~/variant_call/materials/JG.fa \
  --variant ~/variant_call/intermediate/merged.indel.vcf.gz \
  --output ~/variant_call/intermediate/merged.indel.hardfiltering.vcf.gz \
  --filter-expression "QD < 2.0" \
  --filter-expression "FS > 200.0" \
  --filter-expression "SOR > 10.0" \
  --filter-expression "ReadPosRankSum < -20.0"
```

### 最終的なVCFファイルの統合

```bash
# SNVとINDELのそれぞれでフィルタリングしたVCFをMergeVcfsで1つに統合する
~/variant_call/tools/gatk-4.3.0.0/gatk --java-options "-Xmx4G" MergeVcfs \
  --OUTPUT ~/variant_call/results/merged.hardfiltering.vcf.gz \
  --INPUT ~/variant_call/intermediate/merged.snv.hardfiltering.vcf.gz \
  --INPUT ~/variant_call/intermediate/merged.indel.hardfiltering.vcf.gz
```