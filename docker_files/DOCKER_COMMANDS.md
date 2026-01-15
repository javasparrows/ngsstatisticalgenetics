# Docker環境 クイックリファレンス

**注意**: 全てのDockerコマンドは `docker_files/` ディレクトリ内で実行してください。

```bash
cd docker_files
```

## よく使うDockerコマンド

### 基本操作

```bash
# 古いコンテナを削除
docker rm -f ngs-variant-calling

# コンテナを起動（バックグラウンド）
docker compose up -d

# コンテナに入る
docker compose exec ngs-pipeline bash

bash download_files.sh
bash run_all.sh

# コンテナから出る（コンテナ内で実行）
exit

# コンテナの状態を確認
docker compose ps

# ログを確認
docker compose logs
docker compose logs -f  # リアルタイム表示

# コンテナを停止
docker compose stop

# コンテナを停止して削除
docker compose down

# コンテナ、ボリューム、ネットワークを全て削除（注意）
docker compose down -v
```

---

## 完全リセット手順

全てをクリーンにして最初からやり直す場合：

```bash
# 1. コンテナとボリュームを全て削除
docker compose down -v

# 2. イメージを削除
docker rmi ngs-statistical-genetics:latest

# 3. システム全体をクリーンアップ（オプション）
docker system prune -a

# 4. ホスト側のデータディレクトリをクリーンアップ（必要に応じて）
rm -rf data/* results/* intermediate/* logs/*

# 5. 再ビルドと起動
docker compose build --no-cache
docker compose up -d

# 6. データのダウンロードから再開
docker compose exec ngs-pipeline bash
bash download_files.sh
bash run_all.sh
```

---

## ファイル修正時の対応コマンド

### 1. download_files.sh や run_all.sh を修正した場合

スクリプトファイルはボリュームマウントされているため、**再ビルド不要**です。

```bash
# コンテナが起動中の場合、そのまま使用可能
docker compose exec ngs-pipeline bash

# コンテナ内で修正されたスクリプトを実行
bash download_files.sh
bash run_all.sh
```

**理由**: docker-compose.ymlで以下のようにマウントされているため、ホスト側の変更が即座にコンテナ内に反映されます。
```yaml
volumes:
  - ./download_files.sh:/variant_call/download_files.sh
  - ./run_all.sh:/variant_call/run_all.sh
```

---

### 2. Dockerfile を修正した場合

Dockerfileの変更は**イメージの再ビルドが必要**です。

```bash
# 方法1: 通常の再ビルド（推奨）
docker compose build

# 方法2: キャッシュを使わずに完全に再ビルド（エラーが出た場合）
docker compose build --no-cache

# 方法3: 完全にクリーンにしてから再ビルド（問題が解決しない場合）
docker compose down -v
docker compose build --no-cache
docker compose up -d
```

**注意**: `-v` オプションを使用すると、名前付きボリューム（ngs-tools）も削除されます。ツールの再ダウンロードが必要になります。

---

### 3. docker-compose.yml を修正した場合

#### 3-1. リソース制限やマウント設定のみを変更した場合

```bash
# コンテナを再作成
docker compose down
docker compose up -d

# または restart（設定の一部は反映されない場合があります）
docker compose restart
```

#### 3-2. build設定（platformsなど）を変更した場合

```bash
# イメージを再ビルド
docker compose down
docker compose build
docker compose up -d
```

---

### イメージ管理

```bash
# イメージをビルド
docker compose build

# キャッシュなしでビルド
docker compose build --no-cache

# イメージの一覧を表示
docker images

# 使用していないイメージを削除
docker image prune

# 全ての停止中のコンテナとイメージを削除（注意）
docker system prune -a
```

### データ管理

```bash
# ボリュームの一覧を表示
docker volume ls

# 特定のボリュームの詳細を確認
docker volume inspect ngs-tools

# 使用していないボリュームを削除
docker volume prune

# 特定のボリュームを削除（例: ngs-tools）
docker volume rm ngs-tools
```

### ファイル操作

```bash
# ホストからコンテナにファイルをコピー
docker cp /path/to/local/file ngs-variant-calling:/variant_call/

# コンテナからホストにファイルをコピー
docker cp ngs-variant-calling:/variant_call/results/merged.vcf.gz ./

# ただし、マウントされたディレクトリを使う方が簡単：
# ./data/ → /variant_call/materials
# ./results/ → /variant_call/results
# ./intermediate/ → /variant_call/intermediate
# ./logs/ → /variant_call/logs
```

### リソース監視

```bash
# コンテナのリソース使用状況をリアルタイム表示
docker stats ngs-variant-calling

# コンテナのプロセスを表示
docker compose top

# コンテナ内のディスク使用量を確認
docker compose exec ngs-pipeline df -h

# Dockerが使用しているディスク容量を確認
docker system df
```

---

## トラブルシューティング

### ビルドが失敗する場合

```bash
# 1. キャッシュをクリアして再ビルド
docker compose build --no-cache

# 2. 完全にクリーンにしてから再ビルド
docker compose down -v
docker system prune -a
docker compose build --no-cache

# 3. Apple Silicon Macでx86_64エミュレーションの確認
grep "platform" Dockerfile docker-compose.yml
# 出力に "linux/amd64" が含まれているか確認
```

### コンテナが起動しない場合

```bash
# ログを確認
docker compose logs

# 詳細なログを確認
docker compose logs -f ngs-pipeline

# コンテナの状態を確認
docker compose ps -a
```

### メモリ不足エラーが出る場合

```bash
# docker-compose.ymlのメモリ制限を調整
# deploy.resources.limits.memory の値を減らす

# Docker Desktop（Mac/Windows）の場合
# Settings → Resources → Memory を増やす
```

### ディスク容量が不足する場合

```bash
# 不要なイメージとコンテナを削除
docker system prune -a

# 中間ファイルを削除（コンテナ内で実行）
docker compose exec ngs-pipeline bash
cd /variant_call/intermediate
rm -f *.sam
rm -f joint_genotyped.chr*.vcf
rm -rf genomics_database.*
```

---

## クイックリファレンス表

| 修正したファイル | 必要なコマンド | 理由 |
|----------------|--------------|------|
| `download_files.sh` | なし（即座に反映） | ボリュームマウント |
| `run_all.sh` | なし（即座に反映） | ボリュームマウント |
| `Dockerfile` | `docker compose build` | イメージ再ビルド必要 |
| `docker-compose.yml`（リソース） | `docker compose down && docker compose up -d` | コンテナ再作成 |
| `docker-compose.yml`（build設定） | `docker compose build && docker compose up -d` | イメージ再ビルド必要 |
| エラーが出た場合 | `docker compose build --no-cache` | キャッシュクリア |
| 完全リセット | `docker compose down -v && docker compose build --no-cache` | 全削除＋再ビルド |

---

**更新日**: 2026-01-15
