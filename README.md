# NGS Statistical Genetics - Variant Calling Pipeline

ゲノム解析のvariant calling（バリアントコーリング）パイプラインを統一されたDocker環境で実行するためのリポジトリです。

## 概要

このプロジェクトは、次世代シーケンシング（NGS）データを用いてゲノム変異を検出するためのパイプラインを提供します。大容量のFASTQファイル（10GB×6個）を使用し、アライメントからバリアントフィルタリングまでの全工程を自動化しています。

### 使用ツール
- **BWA-MEM2**: 高速アライメントツール
- **GATK**: Genome Analysis Toolkit（バリアント検出）
- **samtools**: BAMファイル操作
- **bcftools**: VCFファイル操作

### パイプライン概要
1. 環境設定・ツール確認
2. FASTQファイルダウンロード（HG005/HG006/HG007家系）
3. 参照ゲノム（JG2.1.0）準備
4. リードアライメント
5. BAMファイル処理（ソート、重複除去）
6. Base Quality Score Recalibration（BQSR）
7. バリアントコーリング
8. ジョイントジェノタイピング
9. バリアントフィルタリング

## 🚀 クイックスタート

### 前提条件
- **Docker Engine** 20.10以上
- **利用可能ディスク容量**: 200GB以上
- **メモリ**: 16GB以上推奨（32GB以上が理想）
- **CPU**: 8コア以上推奨
- **ネットワーク**: 高速回線推奨（60GB+のダウンロード）

### 1. コード取得

#### 方法A: GitHubからクローン

```bash
git clone https://github.com/javasparrows/ngsstatisticalgenetics.git
cd ngsstatisticalgenetics
```

#### 方法B: Google DriveからZIPダウンロード

Gitが使用できない環境の場合：

1. **コードのダウンロード**
   - [ngsstatisticalgenetics.zip をダウンロード](https://drive.google.com/drive/folders/19oRvkyUTIScctLUZ9BVYiaYIaxSqgQ37?usp=drive_link)

2. **ファイル展開と配置**
   ```bash
   # 作業用ディレクトリ作成（推奨）
   mkdir -p ~/Projects
   cd ~/Projects
   
   # ダウンロードフォルダからZIPファイルをコピー
   cp ~/Downloads/ngsstatisticalgenetics.zip .
   
   # ZIPファイルを展開
   unzip ngsstatisticalgenetics.zip
   
   # プロジェクトディレクトリに移動
   cd ngsstatisticalgenetics
   
   # 不要になったZIPファイルを削除（任意）
   rm ../ngsstatisticalgenetics.zip
   ```

### 2. 実行権限設定

```bash
chmod +x scripts/data_downloader.sh
chmod +x run_all.sh
chmod +x command_files/*.sh
```

### 3. 大容量データの事前ダウンロード

#### 方法A: 自動ダウンローダー使用

```bash
# 自動ダウンローダーを実行（約60GBのデータをダウンロード）
./scripts/data_downloader.sh

# ダウンロード完了まで待機（数時間〜半日程度）
```

#### 方法B: Google Driveから直接ダウンロード（推奨）

ネットワーク速度が遅い場合や、事前に準備されたデータセットを使用する場合：

1. **Google Driveからダウンロード**
   - [共有データセット（60GB）をダウンロード](https://drive.google.com/drive/folders/19oRvkyUTIScctLUZ9BVYiaYIaxSqgQ37?usp=drive_link)

2. **ファイル配置**
   ```bash
   # データディレクトリ作成
   mkdir -p ~/variant_call_data/{materials,intermediate,results,logs}
   ```

ダウンロードしたファイルを ~/variant_call_data/materials/ の中に移動またはコピーする。

3. **ファイル配置確認**
   ```bash
   # 必要ファイルの確認
   ls ~/variant_call_data/materials/
   # 以下が存在することを確認：
   # - HG005_Son.R1.10X.fastq.gz, HG005_Son.R2.10X.fastq.gz
   # - HG006_Father.R1.10X.fastq.gz, HG006_Father.R2.10X.fastq.gz
   # - HG007_Mother.R1.10X.fastq.gz, HG007_Mother.R2.10X.fastq.gz
   ```

### 4. Docker環境構築

```bash
# 現在のディレクトリ確認（ngsstatisticalgeneticsフォルダ内にいることを確認）
pwd
# 出力例: /Users/username/Projects/ngsstatisticalgenetics

# Dockerfileの存在確認
ls Dockerfile
# Dockerfile が表示されることを確認

# Dockerイメージをビルド（30-60分程度）
docker build -t variant-calling:latest .

# ビルド成功確認
docker images | grep variant-calling
```

### 5. パイプライン実行

```bash
# Dockerコンテナ起動（フル解析）
docker run -it --rm \
  --name variant-calling-pipeline \
  --cpus="8" \
  --memory="32g" \
  -v ~/variant_call_data/materials:/workspace/materials:ro \
  -v ~/variant_call_data/intermediate:/workspace/intermediate \
  -v ~/variant_call_data/results:/workspace/results \
  -v ~/variant_call_data/logs:/workspace/logs \
  variant-calling:latest

# コンテナ内でパイプライン実行
./run_all.sh
```

## 📁 プロジェクト構造

ZIPファイル展開後のディレクトリ構造：

```
~/Projects/ngsstatisticalgenetics/  # ← ここでdocker buildを実行
├── Dockerfile                      # Docker環境定義
├── DOCKER.md                       # 詳細な環境構築ガイド
├── README.md                       # このファイル
├── run_all.sh                      # パイプライン統合実行スクリプト
├── command_files/                  # 個別処理スクリプト（12ステップ）
│   ├── 01_setup_environment.sh
│   ├── 02_download_fastq.sh
│   ├── 03_download_reference.sh
│   ├── 04_prepare_reference.sh
│   ├── 05_alignment.sh
│   ├── 06_sam_to_bam.sh
│   ├── 07_sort_bam.sh
│   ├── 08_mark_duplicates.sh
│   ├── 09_bqsr.sh
│   ├── 10_variant_calling.sh
│   ├── 11_joint_genotyping.sh
│   └── 12_variant_filtering.sh
├── scripts/
│   └── data_downloader.sh          # ホスト側データ準備スクリプト
└── 1-3/
    └── command.md                  # 元のコマンド記録

データ保存用（別途作成される）：
~/variant_call_data/
├── materials/                      # 入力データ（60GB+）
├── intermediate/                   # 中間処理ファイル
├── results/                        # 最終結果
└── logs/                          # ログファイル
```

## 🔧 詳細な使用方法

### データダウンロードの詳細確認

```bash
# ダウンロード状況確認
./scripts/data_downloader.sh --verify-only

# カスタムディレクトリ使用
./scripts/data_downloader.sh --data-dir /path/to/your/data
```

### パイプライン部分実行

```bash
# 特定ステップ範囲実行（例：ステップ5-8のみ）
./run_all.sh --from 5 --to 8

# 実行内容確認（ドライラン）
./run_all.sh --dry-run

# ヘルプ表示
./run_all.sh --help
```

### 個別ステップ実行

```bash
# 特定のステップのみ実行
./command_files/05_alignment.sh
./command_files/10_variant_calling.sh
```

### リソース調整

```bash
# 高性能環境での実行
docker run -it --rm \
  --cpus="16" \
  --memory="64g" \
  --shm-size="8g" \
  -v ~/variant_call_data/materials:/workspace/materials:ro \
  -v ~/variant_call_data/intermediate:/workspace/intermediate \
  -v ~/variant_call_data/results:/workspace/results \
  variant-calling:latest
```

## 📊 出力ファイル

### 主要な結果ファイル

パイプライン完了後、以下のファイルが `~/variant_call_data/results/` に生成されます：

```bash
results/
├── son.sort.markdup.bam        # Son（息子）の前処理済みBAM
├── father.sort.markdup.bam     # Father（父）の前処理済みBAM
├── mother.sort.markdup.bam     # Mother（母）の前処理済みBAM
├── merged.vcf.gz               # 生の変異データ（全染色体統合）
└── merged.hardfiltering.vcf.gz # フィルタリング後の変異データ（最終結果）
```

### ログファイル

```bash
logs/
├── pipeline_YYYYMMDD_HHMMSS.log      # パイプライン全体ログ
├── 01_setup_environment_*.log         # 各ステップの詳細ログ
├── 02_download_fastq_*.log
├── ...
└── 12_variant_filtering_*.log
```

## 🔍 トラブルシューティング

### よくある問題

#### 1. メモリ不足エラー
```bash
# Java heap size調整
export JAVA_OPTS="-Xmx8g"
# または docker run時にメモリ制限調整
docker run --memory="64g" ...
```

#### 2. ディスク容量不足
```bash
# 中間ファイル削除（SAMファイルはBAM作成後に削除可能）
rm -f ~/variant_call_data/intermediate/*.sam
```

#### 3. ダウンロード失敗
```bash
# aria2cのタイムアウト・リトライ調整
aria2c --timeout=300 --retry-wait=30 --max-tries=5 [URL]
```

#### 4. 権限エラー
```bash
# データディレクトリ権限修正
sudo chown -R $(id -u):$(id -g) ~/variant_call_data
```

### ログ確認方法

```bash
# リアルタイムログ監視
tail -f ~/variant_call_data/logs/pipeline_*.log

# エラー検索
grep -i "error" ~/variant_call_data/logs/*.log
```

## 🌐 研究チーム内でのデータ共有

### NFSサーバー使用

```bash
# NFSマウント
sudo mkdir -p /mnt/shared_genomics_data
sudo mount -t nfs server.example.com:/genomics/data /mnt/shared_genomics_data

# 共有データ使用
docker run -it --rm \
  -v /mnt/shared_genomics_data/materials:/workspace/materials:ro \
  -v /mnt/shared_genomics_data/intermediate:/workspace/intermediate \
  variant-calling:latest
```

### クラウドストレージ同期

```bash
# AWS S3から同期
aws s3 sync s3://your-bucket/genomics-data ~/variant_call_data/materials

# Google Cloud Storage
gsutil -m rsync -r gs://your-bucket/genomics-data ~/variant_call_data/materials
```

## 🔬 解析カスタマイズ

### 30Xデータセットの使用

```bash
# command_files/02_download_fastq.sh を編集して30Xデータセット使用
# （各ファイル約30GB、より高精度だが処理時間も増加）
```

### 追加のバリアントフィルタリング

```bash
# カスタムフィルタリング条件
gatk VariantFiltration \
  --filter-expression "QUAL < 30.0" --filter-name "LowQual" \
  --filter-expression "DP < 10" --filter-name "LowDepth"
```

## 📚 参考資料

- **詳細な環境構築**: [DOCKER.md](DOCKER.md)
- **元のコマンド記録**: [1-3/command.md](1-3/command.md)
- **GATK公式ドキュメント**: https://gatk.broadinstitute.org/
- **BWA-MEM2**: https://github.com/bwa-mem2/bwa-mem2

## 🤝 コントリビューション

バグ報告や改善提案は、GitHubのIssuesでお知らせください。

## 📄 ライセンス

このプロジェクトの各ツールは、それぞれのライセンスに従います：
- GATK: BSD 3-Clause License
- BWA-MEM2: MIT License
- samtools/bcftools: MIT/Expat License

## 🏥 引用

このパイプラインを研究で使用する場合は、以下のツールを適切に引用してください：
- GATK: Van der Auwera & O'Connor. (2020). Genomics in the Cloud. O'Reilly.
- BWA-MEM2: Vasimuddin, et al. (2019). BWA-MEM2: fast and memory-efficient genomic alignment. Bioinformatics.
- samtools: Li, et al. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics.
