#!/bin/bash
set -e

# Create /tools directory if not exists
mkdir -p /tools
cd /tools

echo "Installing bioinformatics tools..."

# Install MAFFT
if ! command -v mafft &> /dev/null; then
    echo "Installing MAFFT..."
    wget -q --show-progress https://mafft.cbrc.jp/alignment/software/mafft_7.526-1_amd64.deb
    dpkg -i mafft_7.526-1_amd64.deb
    rm mafft_7.526-1_amd64.deb
    echo "✓ MAFFT installed successfully"
else
    echo "✓ MAFFT already installed"
fi


# Install BLAST
if [ ! -d "blast" ]; then
    echo "Installing BLAST..."
    wget -q --show-progress https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz
    tar -xzf ncbi-blast-2.17.0+-x64-linux.tar.gz
    mv ncbi-blast-2.17.0+ blast
    rm ncbi-blast-2.17.0+-x64-linux.tar.gz
    echo "✓ BLAST installed successfully"
else
    echo "✓ BLAST already installed"
fi

# Install Diamond
if [ ! -f "diamond/diamond" ]; then
    echo "Installing Diamond..."
    mkdir -p diamond
    wget -q --show-progress https://github.com/bbuchfink/diamond/releases/download/v2.1.14/diamond-linux64.tar.gz
    tar -xzf diamond-linux64.tar.gz -C diamond
    chmod +x diamond/diamond
    rm diamond-linux64.tar.gz
    echo "✓ Diamond installed successfully"
else
    echo "✓ Diamond already installed"
fi

# Install MMSeqs2
if [ ! -d "mmseqs" ]; then
    echo "Installing MMSeqs2..."
    wget -q --show-progress https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
    tar -xzf mmseqs-linux-avx2.tar.gz
    chmod +x mmseqs/bin/mmseqs
    rm mmseqs-linux-avx2.tar.gz
    echo "✓ MMSeqs2 installed successfully"
else
    echo "✓ MMSeqs2 already installed"
fi

# Install GeMoMa
if [ ! -d "gemoma" ]; then
    echo "Installing GeMoMa..."
    mkdir -p gemoma
    wget -q --show-progress -O GeMoMa.zip "http://www.jstacs.de/download.php?which=GeMoMa"
    unzip -q GeMoMa.zip -d gemoma
    rm GeMoMa.zip
    echo "✓ GeMoMa installed successfully"
else
    echo "✓ GeMoMa already installed"
fi

# Install FastTree
if [ ! -d "fasttree" ]; then
    echo "Installing FastTree..."
    mkdir -p fasttree
    wget -q --show-progress -O fasttree/FastTree https://morgannprice.github.io/fasttree/FastTree
    chmod +x fasttree/FastTree
    echo "✓ FastTree installed successfully"
else
    echo "✓ FastTree already installed"
fi

echo "All bioinformatics tools installed successfully!"

# Verify installations
echo "Verifying tool installations..."
ls -la /tools/
echo "Installation verification complete."