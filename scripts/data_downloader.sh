#!/bin/bash
# Host-side data downloader script
# This script should be run on the host machine before starting the Docker container

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$HOME/variant_call_data"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" >&2
}

# Create data directories
setup_directories() {
    log "Setting up data directories..."
    mkdir -p "$DATA_DIR"/{materials,intermediate,results,logs}
    log "Data directories created at: $DATA_DIR"
}

# Download large files
download_files() {
    cd "$DATA_DIR/materials"
    
    log "Starting download of large data files..."
    log "This will download approximately 60GB+ of data. Ensure sufficient disk space."
    log "Download location: $DATA_DIR/materials"
    
    # Check if aria2c is available
    if ! command -v aria2c &> /dev/null; then
        log_error "aria2c is not installed. Please install it first:"
        log_error "  Ubuntu/Debian: sudo apt-get install aria2"
        log_error "  macOS: brew install aria2"
        log_error "  CentOS/RHEL: sudo yum install aria2"
        exit 1
    fi
    
    # Download 10X FASTQ files (each ~10GB)
    log "Downloading 10X FASTQ files..."
    aria2c -x 8 -s 8 -c https://zenodo.org/record/7310196/files/HG005_Son.R1.10X.fastq.gz
    aria2c -x 8 -s 8 -c https://zenodo.org/record/7310196/files/HG005_Son.R2.10X.fastq.gz
    aria2c -x 8 -s 8 -c https://zenodo.org/record/7310202/files/HG006_Father.R1.10X.fastq.gz
    aria2c -x 8 -s 8 -c https://zenodo.org/record/7310202/files/HG006_Father.R2.10X.fastq.gz
    aria2c -x 8 -s 8 -c https://zenodo.org/record/7310233/files/HG007_Mother.R1.10X.fastq.gz
    aria2c -x 8 -s 8 -c https://zenodo.org/record/7310233/files/HG007_Mother.R2.10X.fastq.gz
    
    # Download reference genome
    log "Downloading JG2.1.0 reference genome..."
    wget https://jmorp.megabank.tohoku.ac.jp/datasets/tommo-jg2.1.0-20211208/files/jg2.1.0.fa.gz
    gzip -cd jg2.1.0.fa.gz > JG.fa && rm jg2.1.0.fa.gz
    
    # Download GATK Resource Bundle
    log "Downloading GATK Resource Bundle..."
    wget https://jmorp.megabank.tohoku.ac.jp/datasets/tommo-jg2.1.0-20211208/files/JG2.1.0-ResourceBundle-from-b37.zip
    unzip JG2.1.0-ResourceBundle-from-b37.zip && rm JG2.1.0-ResourceBundle-from-b37.zip
    
    log "All files downloaded successfully!"
}

# Verify downloaded files
verify_files() {
    log "Verifying downloaded files..."
    
    expected_files=(
        "HG005_Son.R1.10X.fastq.gz"
        "HG005_Son.R2.10X.fastq.gz"
        "HG006_Father.R1.10X.fastq.gz"
        "HG006_Father.R2.10X.fastq.gz"
        "HG007_Mother.R1.10X.fastq.gz"
        "HG007_Mother.R2.10X.fastq.gz"
        "JG.fa"
        "JG2.1.0-ResourceBundle-from-b37"
    )
    
    cd "$DATA_DIR/materials"
    all_present=true
    
    for file in "${expected_files[@]}"; do
        if [[ ! -e "$file" ]]; then
            log_error "Missing file: $file"
            all_present=false
        else
            log "✓ Found: $file"
        fi
    done
    
    if [[ "$all_present" == "true" ]]; then
        log "All expected files are present!"
        log "Total disk usage:"
        du -sh .
    else
        log_error "Some files are missing. Please check the download process."
        exit 1
    fi
}

# Display next steps
show_next_steps() {
    log ""
    log "=== Data preparation completed! ==="
    log ""
    log "Data location: $DATA_DIR"
    log ""
    log "Next steps:"
    log "1. Build the Docker image:"
    log "   cd $(dirname "$SCRIPT_DIR")"
    log "   docker build -t variant-calling:latest ."
    log ""
    log "2. Run the Docker container:"
    log "   docker run -it --rm \\"
    log "     -v $DATA_DIR/materials:/workspace/materials:ro \\"
    log "     -v $DATA_DIR/intermediate:/workspace/intermediate \\"
    log "     -v $DATA_DIR/results:/workspace/results \\"
    log "     -v $DATA_DIR/logs:/workspace/logs \\"
    log "     variant-calling:latest"
    log ""
    log "3. Inside the container, run the pipeline:"
    log "   ./run_all.sh"
    log ""
}

# Main function
main() {
    log "=== Variant Calling Data Downloader ==="
    log "This script will download large genomics data files for the variant calling pipeline."
    
    # Ask for confirmation
    echo ""
    read -p "This will download 60GB+ of data. Do you want to continue? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        log "Download cancelled."
        exit 0
    fi
    
    setup_directories
    download_files
    verify_files
    show_next_steps
}

# Usage message
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "This script downloads large data files required for the variant calling pipeline."
    echo ""
    echo "Options:"
    echo "  -h, --help     Show this help message"
    echo "  --data-dir     Specify custom data directory (default: ~/variant_call_data)"
    echo "  --verify-only  Only verify existing files without downloading"
    echo ""
}

# Command line argument processing
VERIFY_ONLY=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            usage
            exit 0
            ;;
        --data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        --verify-only)
            VERIFY_ONLY=true
            shift
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Execute based on mode
if [[ "$VERIFY_ONLY" == "true" ]]; then
    if [[ -d "$DATA_DIR/materials" ]]; then
        verify_files
    else
        log_error "Data directory not found: $DATA_DIR/materials"
        exit 1
    fi
else
    main
fi