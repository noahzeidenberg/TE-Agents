#!/bin/bash

# Load required modules
module load perl
module load blast+
module load muscle
module load repeatmasker
module load cd-hit
module load trf
module load gt

# Clone EDTA repository
git clone https://github.com/oushujun/EDTA.git

# Install EDTA dependencies using cpanm
curl -L https://cpanmin.us | perl - --sudo App::cpanminus
cpanm File::Which
cpanm Text::Soundex
cpanm Hash::Merge
cpanm LWP::UserAgent

# Make EDTA scripts executable
chmod +x EDTA/EDTA.pl
chmod +x EDTA/util/*

# Test EDTA installation
perl EDTA/EDTA.pl -h 