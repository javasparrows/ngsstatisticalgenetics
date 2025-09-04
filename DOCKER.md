# Docker環境構築方針

## プロジェクト概要
ゲノム解析のバリアントコーリングパイプライン。大容量FASTQファイル（30GB×6個、10GB×6個）を用いたvariant calling解析を実行する。

## Docker構築方式

### 1. 基本方針
- **Dockerfileのみ**: docker-compose.ymlは使用せず、シンプルなDockerfileベースの構築
- **マルチステージビルド**: ツールビルド用とランタイム用を分離してイメージサイズを最適化
- **大容量ファイルは外部マウント**: Docker imageに含めずボリュームマウントで対応

### 2. ファイル分類と配置方針

#### 2.1 Docker image内に含めるもの
- **ツール・ソフトウェア**:
  - samtools
  - BWA-MEM2
  - GATK
  - aria2c
  - bcftools
  - Java（OpenJDK）
  - Python

#### 2.2 外部マウントするもの
- **大容量データファイル**:
  - FASTQファイル（30GB×6個、10GB×6個）
  - 参照ゲノム（JG2.1.0, GRCh37, GRCh38）
  - GATK Resource Bundle
  - 中間ファイル・結果ファイル

### 3. ディレクトリ構造設計

```
/workspace/
├── tools/          # ツールバイナリ（Docker内）
├── materials/      # 外部マウント：入力データ
├── intermediate/   # 外部マウント：中間処理ファイル
└── results/        # 外部マウント：結果出力
```

### 4. Docker構成詳細

#### 4.1 ベースイメージ
- `ubuntu:22.04`: 安定性とパッケージ豊富さから選択

#### 4.2 マルチステージビルド構成

**Stage 1: ツールビルド**
- 開発ツール（gcc, make, cmake, wget）
- ソースコードのダウンロードとコンパイル
- samtools, BWA-MEM2の手動ビルド

**Stage 2: ランタイム**
- 最小限のランタイム依存関係
- 事前ビルドしたバイナリをコピー
- GATK（Java実行形式）

#### 4.3 環境変数
- PATH設定（/workspace/tools）
- JAVA_OPTS設定

### 5. 実行方式

#### 5.1 データ準備
```bash
# ホストマシンでデータ用ディレクトリ作成
mkdir -p ~/variant_call_data/{materials,intermediate,results}

# 大容量ファイルを事前ダウンロード（ホスト側で実行）
cd ~/variant_call_data/materials
# aria2c等でFASTQファイルをダウンロード
```

#### 5.2 Docker実行
```bash
# Docker imageビルド
docker build -t variant-calling:latest .

# 解析実行
docker run -it --rm \
  -v ~/variant_call_data/materials:/workspace/materials:ro \
  -v ~/variant_call_data/intermediate:/workspace/intermediate \
  -v ~/variant_call_data/results:/workspace/results \
  -v ~/variant_call_data/scripts:/workspace/scripts:ro \
  variant-calling:latest bash
```

### 6. 利点とトレードオフ

#### 6.1 利点
- **再現性**: 全研究者が同じ環境で実行可能
- **ポータビリティ**: 異なるOS・ハードウェアでも動作
- **分離性**: ホストシステムを汚染しない
- **バージョン管理**: ツールバージョンを固定

#### 6.2 トレードオフ
- **初回セットアップ時間**: Dockerビルドに30-60分
- **ストレージ要件**: 大容量ファイルは別途管理必要
- **パフォーマンス**: わずかなオーバーヘッド（数%程度）

### 7. 最適化検討事項

#### 7.1 パフォーマンス最適化
- **CPUアフィニティ**: `--cpus`オプションでCPU制限
- **メモリ制限**: `--memory`オプションでメモリ上限設定
- **並列処理**: bwa-mem2, samtoolsの並列オプション活用

#### 7.2 ストレージ最適化
- **ボリュームタイプ**: 高速SSD推奨
- **一時ファイル**: `/tmp`をtmpfsマウント検討
- **圧縮**: 中間ファイルの圧縮保存

### 8. セキュリティ考慮
- **非rootユーザー**: コンテナ内でUID/GIDマッピング
- **読み取り専用マウント**: 入力データは読み取り専用
- **ネットワーク制限**: 解析中のネットワークアクセス制限

### 9. 運用・メンテナンス
- **ログ管理**: 処理ログの永続化
- **バックアップ**: 結果ファイルの定期バックアップ
- **バージョン管理**: Dockerタグでツールバージョン管理
- **ドキュメント**: README.mdでの実行手順明記

この方針により、研究者は大容量データを事前準備した上で、統一されたDocker環境で確実にvariant calling解析を実行できる。

---

## 環境構築手順

### 前提条件
- Docker Engine 20.10以上
- 利用可能ディスク容量: 200GB以上（データ保存用）
- メモリ: 16GB以上推奨（32GB以上が理想）
- CPU: 8コア以上推奨

### 1. リポジトリクローンと事前準備

```bash
# リポジトリクローン
git clone <repository-url>
cd ngsstatisticalgenetics

# 実行権限付与
chmod +x scripts/data_downloader.sh
chmod +x run_all.sh
chmod +x command_files/*.sh
```

### 2. 大容量ファイルの準備

#### 方法A: 自動ダウンローダー使用（推奨）

```bash
# 事前データダウンロード実行
./scripts/data_downloader.sh

# カスタムディレクトリを使用する場合
./scripts/data_downloader.sh --data-dir /path/to/your/data
```

#### 方法B: 手動ダウンロード

```bash
# データディレクトリ作成
mkdir -p ~/variant_call_data/{materials,intermediate,results,logs}
cd ~/variant_call_data/materials

# FASTQ files (10X版、各約10GB)
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310196/files/HG005_Son.R1.10X.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310196/files/HG005_Son.R2.10X.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310202/files/HG006_Father.R1.10X.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310202/files/HG006_Father.R2.10X.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310233/files/HG007_Mother.R1.10X.fastq.gz
aria2c -x 8 -s 8 -c https://zenodo.org/record/7310233/files/HG007_Mother.R2.10X.fastq.gz

# 参照ゲノムとリソースファイル
wget https://jmorp.megabank.tohoku.ac.jp/datasets/tommo-jg2.1.0-20211208/files/jg2.1.0.fa.gz
gzip -cd jg2.1.0.fa.gz > JG.fa && rm jg2.1.0.fa.gz
wget https://jmorp.megabank.tohoku.ac.jp/datasets/tommo-jg2.1.0-20211208/files/JG2.1.0-ResourceBundle-from-b37.zip
unzip JG2.1.0-ResourceBundle-from-b37.zip && rm JG2.1.0-ResourceBundle-from-b37.zip
```

### 3. Docker環境構築

```bash
# Dockerイメージビルド（30-60分程度）
docker build -t variant-calling:latest .

# ビルド成功確認
docker images | grep variant-calling
```

### 4. 解析実行

#### 4.1 フル・パイプライン実行

```bash
# Dockerコンテナ起動
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

#### 4.2 部分実行（特定ステップのみ）

```bash
# ステップ5からステップ8まで実行
./run_all.sh --from 5 --to 8

# ドライラン（実行内容確認）
./run_all.sh --dry-run

# ヘルプ表示
./run_all.sh --help
```

#### 4.3 個別スクリプト実行

```bash
# 特定のスクリプトのみ実行
./command_files/05_alignment.sh
```

### 5. 結果確認

```bash
# 結果ファイル確認
ls -la /workspace/results/

# 主要な出力ファイル
# - son.sort.markdup.bam, father.sort.markdup.bam, mother.sort.markdup.bam
# - merged.vcf.gz (生の変異データ)
# - merged.hardfiltering.vcf.gz (フィルタリング後の変異データ)

# ログファイル確認
ls -la /workspace/logs/
```

## 大容量ファイル共有方法

### 研究チーム内共有方法

#### 方法1: 高速ストレージサーバー（推奨）
```bash
# NFSマウント例
sudo mkdir -p /mnt/shared_genomics_data
sudo mount -t nfs server.example.com:/genomics/variant_call_data /mnt/shared_genomics_data

# 共有データ使用
docker run -it --rm \
  -v /mnt/shared_genomics_data/materials:/workspace/materials:ro \
  -v /mnt/shared_genomics_data/intermediate:/workspace/intermediate \
  -v /mnt/shared_genomics_data/results:/workspace/results \
  variant-calling:latest
```

#### 方法2: オブジェクトストレージ（AWS S3/Google Cloud Storage）
```bash
# AWS S3から同期
aws s3 sync s3://your-genomics-bucket/variant_call_data/materials ~/variant_call_data/materials

# Google Cloud Storageから同期
gsutil -m rsync -r gs://your-genomics-bucket/variant_call_data/materials ~/variant_call_data/materials
```

#### 方法3: 専用転送ツール（Aspera/rclone）
```bash
# rclone設定例
rclone config  # 設定実行
rclone sync remote:genomics_data/materials ~/variant_call_data/materials -P
```

### データ整合性確保

#### チェックサム検証
```bash
# SHA256チェックサム作成
cd ~/variant_call_data/materials
find . -type f -exec sha256sum {} + > checksums.sha256

# チェックサム検証
sha256sum -c checksums.sha256
```

#### ファイル存在確認スクリプト
```bash
# 必要ファイル確認
./scripts/data_downloader.sh --verify-only
```

### セキュリティ考慮事項

- **データ暗号化**: 機密データは保存・転送時に暗号化
- **アクセス制御**: 適切な権限設定（研究者のみアクセス）
- **監査ログ**: データアクセスログの記録・管理

## パフォーマンス最適化

### リソース割り当て調整
```bash
# CPU・メモリ制限調整
docker run -it --rm \
  --cpus="16" \
  --memory="64g" \
  --shm-size="8g" \
  -v ~/variant_call_data/materials:/workspace/materials:ro \
  variant-calling:latest
```

### 高速ストレージ活用
```bash
# SSD上に一時ディレクトリ作成
docker run -it --rm \
  --tmpfs /tmp:rw,size=32g \
  -v /fast/ssd/path:/workspace/intermediate \
  variant-calling:latest
```

## トラブルシューティング

### よくある問題と対処法

#### 1. メモリ不足エラー
```bash
# Java heap sizeを調整
export JAVA_OPTS="-Xmx8g"
```

#### 2. ディスク容量不足
```bash
# 中間ファイル清理
rm -f /workspace/intermediate/*.sam  # SAMファイル削除（BAM作成後）
```

#### 3. ネットワークタイムアウト
```bash
# タイムアウト値調整
aria2c --timeout=300 --retry-wait=30 [URL]
```

#### 4. 権限エラー
```bash
# データディレクトリ権限修正
sudo chown -R $(id -u):$(id -g) ~/variant_call_data
```

### ログファイル確認
```bash
# パイプライン全体ログ
tail -f ~/variant_call_data/logs/pipeline_*.log

# 個別ステップログ
ls ~/variant_call_data/logs/
tail -f ~/variant_call_data/logs/05_alignment_*.log
```

## 継続的改良

### バージョン管理
- Dockerイメージにタグ付け
- ツールバージョン固定化
- 設定ファイルの変更履歴管理

### 性能監視
- 処理時間の記録
- リソース使用量の監視
- ボトルネック特定

この包括的な手順により、研究チーム全体で一貫したvariant calling解析環境を構築・運用できます。