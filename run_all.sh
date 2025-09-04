#!/bin/bash
# Master script to run the complete variant calling pipeline

set -e

# ログファイルの設定
LOG_DIR="/workspace/logs"
mkdir -p "$LOG_DIR"
MAIN_LOG="$LOG_DIR/pipeline_$(date +%Y%m%d_%H%M%S).log"

# ログ関数
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$MAIN_LOG"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" | tee -a "$MAIN_LOG" >&2
}

# スクリプトの実行関数
run_script() {
    local script_name="$1"
    local script_path="/workspace/command_files/$script_name"
    local log_file="$LOG_DIR/${script_name%.sh}_$(date +%Y%m%d_%H%M%S).log"
    
    log "Starting $script_name..."
    
    if [ ! -f "$script_path" ]; then
        log_error "Script not found: $script_path"
        return 1
    fi
    
    if ! chmod +x "$script_path"; then
        log_error "Failed to make script executable: $script_path"
        return 1
    fi
    
    if ! "$script_path" 2>&1 | tee "$log_file"; then
        log_error "$script_name failed. Check log: $log_file"
        return 1
    fi
    
    log "$script_name completed successfully"
}

# メイン処理
main() {
    log "=== Starting Variant Calling Pipeline ==="
    log "Main log file: $MAIN_LOG"
    
    # 各ステップの実行
    local scripts=(
        "01_setup_environment.sh"
        "02_download_fastq.sh" 
        "03_download_reference.sh"
        "04_prepare_reference.sh"
        "05_alignment.sh"
        "06_sam_to_bam.sh"
        "07_sort_bam.sh"
        "08_mark_duplicates.sh"
        "09_bqsr.sh"
        "10_variant_calling.sh"
        "11_joint_genotyping.sh"
        "12_variant_filtering.sh"
    )
    
    local total_scripts=${#scripts[@]}
    local current_script=0
    
    for script in "${scripts[@]}"; do
        current_script=$((current_script + 1))
        log "=== Step $current_script/$total_scripts: $script ==="
        
        if run_script "$script"; then
            log "Step $current_script/$total_scripts completed successfully"
        else
            log_error "Pipeline failed at step $current_script/$total_scripts: $script"
            log "You can resume from this step by running individual scripts"
            exit 1
        fi
    done
    
    log "=== Variant Calling Pipeline Completed Successfully ==="
    log "Results are available in /workspace/results/"
    log "Main output files:"
    log "  - /workspace/results/merged.vcf.gz (raw variants)"
    log "  - /workspace/results/merged.hardfiltering.vcf.gz (filtered variants)"
}

# 使用方法の表示
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -h, --help     Show this help message"
    echo "  --dry-run      Show what would be executed without running"
    echo "  --from STEP    Start from specific step (1-12)"
    echo "  --to STEP      Stop at specific step (1-12)"
    echo ""
    echo "Steps:"
    echo "  1. Setup environment"
    echo "  2. Download FASTQ files"
    echo "  3. Download reference genome"
    echo "  4. Prepare reference indices"
    echo "  5. Alignment"
    echo "  6. SAM to BAM conversion"
    echo "  7. Sort BAM files"
    echo "  8. Mark duplicates"
    echo "  9. Base Quality Score Recalibration (BQSR)"
    echo "  10. Variant calling"
    echo "  11. Joint genotyping"
    echo "  12. Variant filtering"
}

# コマンドライン引数の処理
DRY_RUN=false
FROM_STEP=1
TO_STEP=12

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            usage
            exit 0
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --from)
            FROM_STEP="$2"
            shift 2
            ;;
        --to)
            TO_STEP="$2"
            shift 2
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# バリデーション
if [[ ! "$FROM_STEP" =~ ^[1-9]$|^1[0-2]$ ]] || [[ ! "$TO_STEP" =~ ^[1-9]$|^1[0-2]$ ]]; then
    log_error "Step numbers must be between 1 and 12"
    exit 1
fi

if [[ $FROM_STEP -gt $TO_STEP ]]; then
    log_error "FROM_STEP cannot be greater than TO_STEP"
    exit 1
fi

if [[ "$DRY_RUN" == "true" ]]; then
    log "DRY RUN - Would execute steps $FROM_STEP to $TO_STEP"
    exit 0
fi

# 部分実行の場合
if [[ $FROM_STEP -ne 1 ]] || [[ $TO_STEP -ne 12 ]]; then
    log "Running partial pipeline: steps $FROM_STEP to $TO_STEP"
    scripts=(
        "01_setup_environment.sh"
        "02_download_fastq.sh" 
        "03_download_reference.sh"
        "04_prepare_reference.sh"
        "05_alignment.sh"
        "06_sam_to_bam.sh"
        "07_sort_bam.sh"
        "08_mark_duplicates.sh"
        "09_bqsr.sh"
        "10_variant_calling.sh"
        "11_joint_genotyping.sh"
        "12_variant_filtering.sh"
    )
    
    for i in $(seq $((FROM_STEP-1)) $((TO_STEP-1))); do
        script="${scripts[$i]}"
        log "=== Step $((i+1)): $script ==="
        if ! run_script "$script"; then
            log_error "Pipeline failed at step $((i+1)): $script"
            exit 1
        fi
    done
    
    log "=== Partial pipeline completed (steps $FROM_STEP to $TO_STEP) ==="
else
    # フル実行
    main
fi