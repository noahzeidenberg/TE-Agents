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

# Load required modules
module purge
module load StdEnv/2020
module load gcc/9.3.0
module load python/3.11
module load blast+/2.14.0
module load hmmer/3.3.2
module load trf/4.09.1
module load perl/5.30.2

# Install rmblast if not present
if [ ! -d "rmblast" ]; then
    echo "Installing RMBlast..."
    wget https://www.repeatmasker.org/rmblast/rmblast-2.14.1+-x64-linux.tar.gz
    tar zxvf rmblast-2.14.1+-x64-linux.tar.gz
    mv rmblast-2.14.1 rmblast
    rm rmblast-2.14.1+-x64-linux.tar.gz
fi

# Install RepeatMasker if not already present
if [ ! -d "RepeatMasker" ]; then
    echo "Cloning RepeatMasker..."
    git clone https://github.com/rmhubley/RepeatMasker.git
else
    echo "RepeatMasker directory already exists, skipping clone..."
fi

cd RepeatMasker

# Download Dfam library if not already present
if [ ! -f "Dfam.h5" ]; then
    echo "Downloading Dfam library..."
    wget https://www.dfam.org/releases/Dfam_3.7/families/Dfam.h5.gz
    gunzip Dfam.h5.gz
else
    echo "Dfam.h5 already exists, skipping download..."
fi

# Install h5py and numpy in a clean virtual environment
echo "Setting up Python environment..."
python -m venv $TOOLS_DIR/venv
source $TOOLS_DIR/venv/bin/activate
pip install --no-cache-dir numpy h5py

# Configure RepeatMasker
echo "Configuring RepeatMasker..."
perl ./configure \
    -trf_prgm $(which trf) \
    -hmmer_dir $(dirname $(which hmmsearch)) \
    -rmblast_dir $TOOLS_DIR/rmblast/bin \
    -libdir $TOOLS_DIR/RepeatMasker/Libraries \
    -default_search_engine hmmer

cd $TOOLS_DIR

# Install EDTA if not already present
if [ ! -d "EDTA" ]; then
    echo "Cloning EDTA..."
    git clone https://github.com/oushujun/EDTA.git
else
    echo "EDTA directory already exists, skipping clone..."
fi

cd EDTA

# Install cpanm if not available
if ! command -v cpanm &> /dev/null; then
    echo "Installing cpanminus..."
    curl -L https://cpanmin.us | perl - --self-contained App::cpanminus
    CPANM="$TOOLS_DIR/perl5/bin/cpanm"
else
    CPANM="cpanm"
fi

# Install required Perl modules
mkdir -p $TOOLS_DIR/perl5
echo "Installing Perl modules..."
PERL_MM_OPT="INSTALL_BASE=$TOOLS_DIR/perl5" \
PERL5LIB="$TOOLS_DIR/perl5/lib/perl5" \
$CPANM --local-lib=$TOOLS_DIR/perl5 local::lib File::Which Text::Soundex Hash::Merge LWP::UserAgent

# Make EDTA scripts executable
if [ -d "util" ]; then
    chmod +x EDTA.pl
    chmod +x util/*
fi

# Install GenomeTools if not present
if [ ! -d "genometools" ]; then
    echo "Installing GenomeTools..."
    git clone https://github.com/genometools/genometools.git
    cd genometools
    make threads=yes 64bit=yes
    make install prefix=$TOOLS_DIR/genometools
    cd $TOOLS_DIR
fi

# Create/update configuration file
echo "Creating/updating tools configuration..."
cat > $TOOLS_DIR/tools_config.sh << EOF
# Tool paths
export PATH="$TOOLS_DIR/genometools/bin:$TOOLS_DIR/RepeatMasker:$TOOLS_DIR/EDTA:$TOOLS_DIR/rmblast/bin:$PATH"
export PERL5LIB="$TOOLS_DIR/perl5/lib/perl5:$PERL5LIB"
# Python virtual environment
source $TOOLS_DIR/venv/bin/activate
EOF

echo "Setup complete!" 