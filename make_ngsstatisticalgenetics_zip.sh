#!/bin/bash
# ngsstatisticalgenetics.zip作成スクリプト
# 研究者向けの配布用ZIPファイルを作成する

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_NAME="ngsstatisticalgenetics"
OUTPUT_DIR="$HOME/Desktop"
ZIP_NAME="${PROJECT_NAME}.zip"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" >&2
}

# 作業ディレクトリ作成
TEMP_DIR=$(mktemp -d)
WORK_DIR="$TEMP_DIR/$PROJECT_NAME"

log "=== Creating distribution ZIP for ngsstatisticalgenetics ==="
log "Output: $OUTPUT_DIR/$ZIP_NAME"
log "Working directory: $WORK_DIR"

# 必要なファイルをコピー
log "Copying essential files..."
mkdir -p "$WORK_DIR"

# 必要なファイル・ディレクトリのリスト
essential_files=(
    "Dockerfile"
    "DOCKER.md"
    "README.md"
    "run_all.sh"
    "command_files"
    "scripts"
    "1-3/command.md"
)

# ファイルをコピー
for item in "${essential_files[@]}"; do
    if [[ -e "$SCRIPT_DIR/$item" ]]; then
        if [[ -d "$SCRIPT_DIR/$item" ]]; then
            cp -r "$SCRIPT_DIR/$item" "$WORK_DIR/"
            log "✓ Copied directory: $item"
        else
            cp "$SCRIPT_DIR/$item" "$WORK_DIR/"
            log "✓ Copied file: $item"
        fi
    else
        log_error "Missing file/directory: $item"
        exit 1
    fi
done

# 実行権限を設定
log "Setting executable permissions..."
chmod +x "$WORK_DIR/run_all.sh"
chmod +x "$WORK_DIR/scripts/"*.sh
chmod +x "$WORK_DIR/command_files/"*.sh

# 不要なファイルを除外するための.gitignore相当の処理
log "Cleaning up unnecessary files..."
find "$WORK_DIR" -name ".DS_Store" -delete 2>/dev/null || true
find "$WORK_DIR" -name "*.tmp" -delete 2>/dev/null || true
find "$WORK_DIR" -name "Thumbs.db" -delete 2>/dev/null || true

# ZIPファイル作成
log "Creating ZIP file..."
cd "$TEMP_DIR"
zip -r "$OUTPUT_DIR/$ZIP_NAME" "$PROJECT_NAME" -x "*.git*" "*.serena*" "variant_call/*"

# ファイルサイズ確認
if [[ -f "$OUTPUT_DIR/$ZIP_NAME" ]]; then
    FILE_SIZE=$(du -h "$OUTPUT_DIR/$ZIP_NAME" | cut -f1)
    log "✓ ZIP file created successfully: $OUTPUT_DIR/$ZIP_NAME ($FILE_SIZE)"
else
    log_error "Failed to create ZIP file"
    exit 1
fi

# 作業ディレクトリ削除
rm -rf "$TEMP_DIR"

# 内容確認
log "ZIP file contents:"
unzip -l "$OUTPUT_DIR/$ZIP_NAME" | head -20

log ""
log "=== Distribution ZIP created successfully ==="
log ""
log "Next steps:"
log "1. Upload $OUTPUT_DIR/$ZIP_NAME to Google Drive"
log "2. Update the Google Drive sharing link in README.md"
log "3. Share the Google Drive folder with researchers"
log ""
log "ZIP file location: $OUTPUT_DIR/$ZIP_NAME"
log "File size: $FILE_SIZE"

# 検証用の展開テスト
log ""
log "Performing extraction test..."
TEST_DIR="$HOME/Desktop/test_extraction"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"
unzip -q "$OUTPUT_DIR/$ZIP_NAME"

if [[ -d "$TEST_DIR/$PROJECT_NAME" ]]; then
    log "✓ Extraction test successful"
    
    # 必要なファイルの存在確認
    essential_check=(
        "Dockerfile"
        "README.md"
        "run_all.sh"
        "command_files/01_setup_environment.sh"
        "scripts/data_downloader.sh"
    )
    
    all_present=true
    for file in "${essential_check[@]}"; do
        if [[ -f "$TEST_DIR/$PROJECT_NAME/$file" ]]; then
            log "  ✓ $file"
        else
            log_error "  ✗ Missing: $file"
            all_present=false
        fi
    done
    
    if [[ "$all_present" == "true" ]]; then
        log "✓ All essential files are present in the ZIP"
    else
        log_error "Some essential files are missing"
        exit 1
    fi
else
    log_error "Extraction test failed"
    exit 1
fi

# テスト用ディレクトリ削除
rm -rf "$TEST_DIR"

log "=== Script completed successfully ==="