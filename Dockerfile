# Multi-stage build for variant calling pipeline
FROM ubuntu:22.04 AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    wget \
    curl \
    tar \
    bzip2 \
    unzip \
    git \
    autoconf \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Build samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 \
    && tar -jxvf samtools-1.16.1.tar.bz2 \
    && cd samtools-1.16.1 \
    && ./configure --prefix=/build/tools \
    && make \
    && make install

# Build BWA-MEM2
RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 \
    && tar -jxvf bwa-mem2-2.2.1_x64-linux.tar.bz2 \
    && mkdir -p /build/tools/bin \
    && cp bwa-mem2-2.2.1_x64-linux/bwa-mem2* /build/tools/bin/

# Download GATK
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip \
    && unzip gatk-4.3.0.0.zip \
    && mv gatk-4.3.0.0 /build/tools/

# Build bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2 \
    && tar -jxvf bcftools-1.16.tar.bz2 \
    && cd bcftools-1.16 \
    && ./configure --prefix=/build/tools \
    && make \
    && make install

# Runtime stage
FROM ubuntu:22.04

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    openjdk-8-jdk \
    aria2 \
    curl \
    wget \
    bc \
    && rm -rf /var/lib/apt/lists/*

# Set Java environment
ENV JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64

# Create workspace structure
RUN mkdir -p /workspace/{tools,materials,intermediate,results,command_files,scripts}
WORKDIR /workspace

# Copy built tools from builder stage
COPY --from=builder /build/tools /workspace/tools

# Set PATH to include our tools
ENV PATH="/workspace/tools/bin:/workspace/tools:/workspace/tools/gatk-4.3.0.0:${PATH}"

# Create non-root user for security
RUN useradd -m -u 1000 -s /bin/bash researcher
RUN chown -R researcher:researcher /workspace
USER researcher

# Copy command files and scripts
COPY command_files/ /workspace/command_files/
COPY scripts/ /workspace/scripts/
COPY run_all.sh /workspace/

# Make scripts executable
RUN chmod +x /workspace/run_all.sh /workspace/command_files/*.sh /workspace/scripts/*.sh

# Set default command
CMD ["/bin/bash"]