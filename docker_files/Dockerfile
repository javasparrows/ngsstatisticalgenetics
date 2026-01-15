# NGS Statistical Genetics Pipeline Dockerfile
# This Dockerfile sets up a complete environment for variant calling and analysis

# x86_64プラットフォームを明示的に指定（Apple Silicon対応）
FROM --platform=linux/amd64 ubuntu:22.04

# メタデータ
LABEL maintainer="NGS Statistical Genetics Team"
LABEL description="Docker image for NGS variant calling pipeline with GATK, BWA-MEM2, and analysis tools"
LABEL version="1.0"

# 非対話的インストールのための環境変数
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Tokyo

# 基本パッケージとビルドツールのインストール
RUN apt-get update && apt-get install -y \
    # ビルドツール
    build-essential \
    gcc \
    g++ \
    make \
    cmake \
    autoconf \
    automake \
    libtool \
    pkg-config \
    # 圧縮・解凍ツール
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    unzip \
    bzip2 \
    gzip \
    # ダウンロードツール
    wget \
    curl \
    aria2 \
    git \
    # その他のユーティリティ
    vim \
    nano \
    less \
    htop \
    tree \
    # Python関連
    python3 \
    python3-pip \
    python3-dev \
    # Java
    openjdk-17-jdk \
    # ncurses (samtools用)
    libncurses5-dev \
    libncursesw5-dev \
    # bcftools用
    libgsl-dev \
    # その他
    ca-certificates \
    locales \
    && rm -rf /var/lib/apt/lists/*

# ロケールの設定
RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

# Java環境変数
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
ENV PATH=$JAVA_HOME/bin:$PATH

# 作業ディレクトリの作成
RUN mkdir -p /variant_call/tools /variant_call/materials /variant_call/intermediate /variant_call/results /variant_call/logs

# =================================================================
# BWA-MEM2のインストール
# =================================================================
WORKDIR /tmp
RUN git clone --recursive https://github.com/bwa-mem2/bwa-mem2.git && \
    cd bwa-mem2 && \
    make -j$(nproc) && \
    cp bwa-mem2* /usr/local/bin/ && \
    cd .. && \
    rm -rf bwa-mem2

# =================================================================
# htslib, samtools, bcftoolsのインストール
# =================================================================
WORKDIR /tmp

# htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2 && \
    tar -xjf htslib-1.19.1.tar.bz2 && \
    cd htslib-1.19.1 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf htslib-1.19.1*

# samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2 && \
    tar -xjf samtools-1.19.2.tar.bz2 && \
    cd samtools-1.19.2 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf samtools-1.19.2*

# bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 && \
    tar -xjf bcftools-1.19.tar.bz2 && \
    cd bcftools-1.19 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && \
    rm -rf bcftools-1.19*

# =================================================================
# vcfannoのインストール
# =================================================================
RUN wget https://github.com/brentp/vcfanno/releases/download/v0.3.5/vcfanno_linux64 \
    -O /usr/local/bin/vcfanno && \
    chmod +x /usr/local/bin/vcfanno

# =================================================================
# slivarのインストール
# =================================================================
RUN wget https://github.com/brentp/slivar/releases/download/v0.3.1/slivar \
    -O /usr/local/bin/slivar && \
    chmod +x /usr/local/bin/slivar

# =================================================================
# Pythonパッケージのインストール
# =================================================================
RUN pip3 install --no-cache-dir \
    numpy \
    scipy \
    pandas \
    matplotlib \
    seaborn \
    pysam \
    biopython

# =================================================================
# 環境変数の設定
# =================================================================
ENV PATH=/variant_call/tools:/variant_call/tools/bin:/usr/local/bin:$PATH
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

# 作業ディレクトリの設定
WORKDIR /variant_call

# スクリプトのコピー（ビルド時にコピー）
COPY download_files.sh /variant_call/
COPY run_all.sh /variant_call/
RUN chmod +x /variant_call/download_files.sh /variant_call/run_all.sh

# ヘルスチェック用のスクリプト
RUN echo '#!/bin/bash\n\
echo "Checking installed tools..."\n\
echo "Java: $(java -version 2>&1 | head -n 1)"\n\
echo "Python: $(python3 --version)"\n\
echo "BWA-MEM2: $(bwa-mem2 version 2>&1 | head -n 1 || echo "Not found")"\n\
echo "Samtools: $(samtools --version | head -n 1)"\n\
echo "bcftools: $(bcftools --version | head -n 1)"\n\
echo "vcfanno: $(vcfanno 2>&1 | head -n 1 || echo "Installed")"\n\
echo "slivar: $(slivar --version 2>&1 || echo "Installed")"\n\
' > /usr/local/bin/check_tools.sh && chmod +x /usr/local/bin/check_tools.sh

# エントリーポイントの設定（デフォルトはbash）
CMD ["/bin/bash"]

# ビルド情報
RUN echo "Docker image built on $(date)" > /variant_call/build_info.txt
RUN echo "Image version: 1.0" >> /variant_call/build_info.txt
RUN echo "Base: Ubuntu 22.04" >> /variant_call/build_info.txt
RUN echo "Java: OpenJDK 17" >> /variant_call/build_info.txt
