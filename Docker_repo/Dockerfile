# Use Python 3.12.7 base image
FROM python:3.12.7

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    CONDA_DIR=/opt/conda \
    PATH="/tools/blast/bin:/tools/diamond:/tools/mmseqs/bin:/tools/fasttree:/tools/gemoma:$PATH" \
    GEMOMA_JAR=/tools/gemoma/GeMoMa-1.9.jar \
    MPLCONFIGDIR=/tmp/matplotlib \
    FONTCONFIG_PATH=/etc/fonts \
    FONTCONFIG_FILE=/etc/fonts/fonts.conf \
    FONTCONFIG_CACHEDIR=/tmp/fontconfig \
    XDG_CACHE_HOME=/tmp/.cache \
    MPLBACKEND=Agg

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    unzip \
    git \
    openjdk-17-jre-headless \
    curl \
    bzip2 \
    libglib2.0-0 \
    libgomp1 \
    libxext6 \
    libsm6 \
    libxrender1 \
    cpanminus \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh && \
    $CONDA_DIR/bin/conda clean --all --yes

ENV PATH="$CONDA_DIR/bin:$PATH"

# Accept conda ToS
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Install zstandard in base conda env to suppress warnings
RUN conda install -n base -c conda-forge zstandard -y

# Create /tools directory
RUN mkdir -p /tools

# Create dedicated conda env for your Python script + tools
RUN conda create -p /tools/dupy_env python=3.12.7 -y && \
    conda install -p /tools/dupy_env -c conda-forge -c bioconda \
        busco=6.0.0 \
        agat \
        perl-sort-naturally \
        perl-clone \
        perl-http-message \
        perl-libwww-perl \
    -y && \
    conda clean --all --yes

# Copy and run external tool installation script
COPY install_tools.sh /app/install_tools.sh
RUN chmod +x /app/install_tools.sh && /app/install_tools.sh

# Create directories with proper permissions for non-root users
RUN mkdir -p /tmp/matplotlib /tmp/.cache /tmp/.local /tmp/fontconfig && \
    chmod -R 777 /tmp


# Create working directory
WORKDIR /app

# Copy Python dependencies and Python scripts
COPY requirements.txt .
COPY dupylicate_docker.py .
COPY AGAT_wrapper_docker.py .
COPY Fasta_fix_docker.py .


# Install Python dependencies inside dupy_env
RUN /tools/dupy_env/bin/pip install --no-cache-dir -r requirements.txt

# Make /app writable for non-root users (BUSCO needs to write logs)
RUN chmod 777 /app

# Add conda env to PATH (this makes busco/agat/python available globally)
# Add conda env to PATH but keep Diamond priority
ENV PATH="/tools/diamond:/tools/dupy_env/bin:$PATH" \
    PERL5LIB="/tools/dupy_env/lib/perl5/site_perl:/tools/dupy_env/lib/perl5"

# Verify installations
RUN /tools/dupy_env/bin/busco --help > /dev/null 2>&1 && echo "BUSCO installation verified" || echo "BUSCO installation failed"
RUN /tools/dupy_env/bin/agat_convert_sp_gff2gtf.pl --help > /dev/null 2>&1 && echo "AGAT installation verified" || echo "AGAT installation failed"

# Use bash shell by default
SHELL ["/bin/bash", "-c"]

# Default command runs your Python script using the dupy_env Python
ENTRYPOINT ["/tools/dupy_env/bin/python", "dupylicate_docker.py"]

