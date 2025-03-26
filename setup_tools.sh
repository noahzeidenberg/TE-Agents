#!/bin/bash
#SBATCH --account=def-yourpi
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x-%j.out

# Create installation directory in scratch
TOOLS_DIR="/scratch/$USER/tools"
mkdir -p $TOOLS_DIR
cd $TOOLS_DIR

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Load required modules
echo "Loading required modules..."
module purge
module load StdEnv/2020
module load gcc/9.3.0
module load python/3.11
module load blast+/2.14.0
module load hmmer/3.3.2
module load trf/4.09.1
module load perl/5.30.2

# Install rmblast if not present or if rmblastn is not executable
if [ ! -x "rmblast/bin/rmblastn" ]; then
    echo "Installing RMBlast..."
    rm -rf rmblast  # Clean up any partial installation
    wget https://www.repeatmasker.org/rmblast/rmblast-2.14.1+-x64-linux.tar.gz
    tar zxvf rmblast-2.14.1+-x64-linux.tar.gz
    mv rmblast-2.14.1 rmblast
    rm rmblast-2.14.1+-x64-linux.tar.gz
else
    echo "RMBlast already installed, skipping..."
fi

# Function to download with retry
download_with_retry() {
    local url=$1
    local output=$2
    local max_attempts=3
    local attempt=1
    
    while [ $attempt -le $max_attempts ]; do
        echo "Download attempt $attempt of $max_attempts..."
        if wget --no-check-certificate -q -O "$output" "$url"; then
            return 0
        fi
        attempt=$((attempt + 1))
        sleep 5
    done
    return 1
}

# Try to use pre-installed RepeatMasker from cvmfs first
if [ -d "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/repeatmasker" ]; then
    echo "Using system RepeatMasker..."
    ln -sf /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/repeatmasker $TOOLS_DIR/RepeatMasker
else
    # Install RepeatMasker if not present or not configured
    if [ ! -x "RepeatMasker/RepeatMasker" ] || [ ! -f "RepeatMasker/Libraries/Dfam.h5" ]; then
        echo "Setting up RepeatMasker..."
        
        # Clone only if directory doesn't exist
        if [ ! -d "RepeatMasker" ]; then
            # Try multiple mirrors
            if ! git clone https://github.com/rmhubley/RepeatMasker.git; then
                if ! git clone https://gitlab.com/dfam/repeatmasker.git RepeatMasker; then
                    echo "Failed to clone RepeatMasker from any source"
                    exit 1
                fi
            fi
        fi
        
        cd RepeatMasker
        
        # Download Dfam library if needed with multiple fallback URLs
        if [ ! -f "Dfam.h5" ]; then
            echo "Downloading Dfam library..."
            urls=(
                "https://www.dfam.org/releases/Dfam_3.7/families/Dfam.h5.gz"
                "https://dfam.org/releases/Dfam_3.7/families/Dfam.h5.gz"
                "https://www.dfam.org/releases/current/families/Dfam.h5.gz"
            )
            
            success=false
            for url in "${urls[@]}"; do
                if download_with_retry "$url" "Dfam.h5.gz"; then
                    gunzip Dfam.h5.gz
                    success=true
                    break
                fi
            done
            
            if ! $success; then
                echo "Failed to download Dfam library from any source"
                exit 1
            fi
        fi
        
        # Set up Python environment if needed
        if [ ! -d "$TOOLS_DIR/venv" ]; then
            echo "Setting up Python environment..."
            python -m venv $TOOLS_DIR/venv
            source $TOOLS_DIR/venv/bin/activate
            pip install --no-cache-dir numpy h5py
        fi
        
        # Configure RepeatMasker
        echo "Configuring RepeatMasker..."
        source $TOOLS_DIR/venv/bin/activate  # Ensure h5py is available
        perl ./configure \
            -trf_prgm $(which trf) \
            -hmmer_dir $(dirname $(which hmmsearch)) \
            -rmblast_dir $TOOLS_DIR/rmblast/bin \
            -libdir $TOOLS_DIR/RepeatMasker/Libraries \
            -default_search_engine hmmer
        
        cd $TOOLS_DIR
    else
        echo "RepeatMasker already installed and configured, skipping..."
    fi
fi

# Install EDTA and dependencies if not present
if [ ! -x "EDTA/EDTA.pl" ]; then
    echo "Setting up EDTA..."
    
    # Clone EDTA if needed
    if [ ! -d "EDTA" ]; then
        git clone https://github.com/oushujun/EDTA.git
    fi
    
    cd EDTA
    chmod +x EDTA.pl
    if [ -d "util" ]; then
        chmod +x util/*
    fi
    
    # Install cpanm if needed
    if ! command_exists cpanm; then
        echo "Installing cpanminus..."
        curl -L https://cpanmin.us | perl - --self-contained App::cpanminus
        CPANM="$TOOLS_DIR/perl5/bin/cpanm"
    else
        CPANM="cpanm"
    fi
    
    # Install Perl modules
    echo "Installing Perl modules..."
    mkdir -p $TOOLS_DIR/perl5
    PERL_MM_OPT="INSTALL_BASE=$TOOLS_DIR/perl5" \
    PERL5LIB="$TOOLS_DIR/perl5/lib/perl5" \
    $CPANM --local-lib=$TOOLS_DIR/perl5 \
        local::lib File::Which Text::Soundex Hash::Merge LWP::UserAgent \
        Bio::DB::EUtilities Bio::SearchIO Bio::Tools::GFF
    
    cd $TOOLS_DIR
else
    echo "EDTA already installed, skipping..."
fi

# Install GenomeTools if not present
if [ ! -x "genometools/bin/gt" ]; then
    echo "Installing GenomeTools..."
    rm -rf genometools  # Clean up any partial installation
    git clone https://github.com/genometools/genometools.git
    cd genometools
    make threads=yes 64bit=yes
    make install prefix=$TOOLS_DIR/genometools
    cd $TOOLS_DIR
else
    echo "GenomeTools already installed, skipping..."
fi

# Install AnnoSINE2 if not present
if [ ! -f "EDTA/util/AnnoSINE2.pl" ]; then
    echo "Installing AnnoSINE2..."
    if [ ! -d "AnnoSINE2" ]; then
        git clone https://github.com/oushujun/AnnoSINE2.git
    fi
    cd AnnoSINE2
    chmod +x AnnoSINE2.pl
    cd ..
    cp AnnoSINE2/AnnoSINE2.pl EDTA/util/
else
    echo "AnnoSINE2 already installed, skipping..."
fi

# Update configuration file only if needed
CONFIG_FILE="$TOOLS_DIR/tools_config.sh"
if [ ! -f "$CONFIG_FILE" ] || ! grep -q "AnnoSINE2" "$CONFIG_FILE"; then
    echo "Creating/updating tools configuration..."
    cat > "$CONFIG_FILE" << EOF
# Tool paths
export PATH="$TOOLS_DIR/genometools/bin:$TOOLS_DIR/RepeatMasker:$TOOLS_DIR/EDTA:$TOOLS_DIR/rmblast/bin:$TOOLS_DIR/AnnoSINE2:$PATH"
export PERL5LIB="$TOOLS_DIR/perl5/lib/perl5:$PERL5LIB"
# Python virtual environment
source $TOOLS_DIR/venv/bin/activate
EOF
fi

echo "Setup complete!" 