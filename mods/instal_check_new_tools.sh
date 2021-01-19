#!/bin/bash

green='\033[0;32m'
red='\033[0;31m'
nc='\033[0m' # No Color

# Save original PATH and LD_LIBRARY variables
original_path_variable=$(echo $PATH)
original_ld_library_variable=$(echo $LD_LIBRARY_PATH)


check_files() {
filepath=$1
filename=$(basename $1)

if [ -e $filepath ]
then
        echo -e "${green}The file: $filename is at the right place.${nc}"
else
        echo -e "${red}The file: $filename can not be found at the expected location. Please check and re-download the file.${nc}"
fi
}

check_directories() {
dirpath=$1

if [ -d $dirpath ]
then
        echo -e "${green}The directory: $dirpath exists.${nc}"
else
        echo -e "${red}The directory: $dirpath can not be found at the expected location. Please check and re-run the install script.${nc}"
fi
}

check_tool () {
commandstring=$1
toolname=$2

eval "$commandstring" > /dev/null 2>&1
exit_parameter=$(echo $?)

if [ $exit_parameter == 0 ]
then
        echo -e "${green}$2 seems to be installed correctly.${nc}"
else
        echo -e "${red}Something seems to be wrong with $2. Please check and maybe reinstall again.${nc}"
fi
}

# ERR251006
path="datasets/1000_genomes/ERR251006.filt.fastq"
check_files "$path"
path="datasets/1000_genomes/ERR251006.filt.fa"
check_files "$path"

# ERR016155
path="datasets/1000_genomes/ERR016155.filt.fastq"
check_files "$path"
path="datasets/1000_genomes/ERR016155.filt.fa"
check_files "$path"

# ERR015528
path="datasets/1000_genomes/ERR015528.filt.fastq"
check_files "$path"
path="datasets/1000_genomes/ERR015528.filt.fa"
check_files "$path"

# SRR741411
path="datasets/1000_genomes/SRR741411.filt.fastq"
check_files "$path"
path="datasets/1000_genomes/SRR741411.filt.fa"
check_files "$path"

# DRR001012
path="datasets/ebi/DRR001012.fastq"
check_files "$path"
path="datasets/ebi/DRR001012.fa"
check_files "$path"
path="datasets/BWA/DRR001012/DRR001012.amb"
check_files "$path"
path="datasets/BWA/DRR001012/DRR001012.ann"
check_files "$path"
path="datasets/BWA/DRR001012/DRR001012.bwt"
check_files "$path"
path="datasets/BWA/DRR001012/DRR001012.pac"
check_files "$path"
path="datasets/BWA/DRR001012/DRR001012.sa"
check_files "$path"

# DRR001025
path="datasets/ebi/DRR001025.fastq"
check_files "$path"
path="datasets/ebi/DRR001025.fa"
check_files "$path"
path="datasets/BWA/DRR001025/DRR001025.amb"
check_files "$path"
path="datasets/BWA/DRR001025/DRR001025.ann"
check_files "$path"
path="datasets/BWA/DRR001025/DRR001025.bwt"
check_files "$path"
path="datasets/BWA/DRR001025/DRR001025.pac"
check_files "$path"
path="datasets/BWA/DRR001025/DRR001025.sa"
check_files "$path"

# GRCh38
path="datasets/1000_genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
check_files "$path"
path="datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.amb"
check_files "$path"
path="datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.ann"
check_files "$path"
path="datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.bwt"
check_files "$path"
path="datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.pac"
check_files "$path"
path="datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.sa"
check_files "$path"

# wgs.ANCA.1_200.fsa
path="datasets/clustalOmega/wgs.ANCA.1_200.fsa"
check_files "$path"

# wgs.ANCA.1_400.fsa
path="datasets/clustalOmega/wgs.ANCA.1_400.fsa"
check_files "$path"

# wgs.ANCA.1_500.fsa
path="datasets/clustalOmega/wgs.ANCA.1_500.fsa"
check_files "$path"

# OE-38_R1
path="datasets/SINA/OE-38_R1.fa"
check_files "$path"
path="datasets/SINA/OE-38_R1.fastq"
check_files "$path"

# RefSeq-RDP16S_v2_May2018
path="datasets/SINA/RefSeq-RDP16S_v2_May2018.fa"
check_files "$path"

# GTDB_bac-arc_ssu_r8
path="datasets/SINA/GTDB_bac-arc_ssu_r86.fa"
check_files "$path"

# SILVA_138.1_SSURef_NR99_12_06_20_opt
path="datasets/SINA/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb"
check_files "$path"
path="datasets/SINA/SILVA_138.1_SSURef_NR99_12_06_20_opt.sidx"
check_files "$path"

# nmonchart
path="nmonchart/nmonchart"
check_files "$path"

# cifar-10-batches-bin
path="datasets/tensorflow/cifar-10-batches-bin"
check_directories "$path"

# Tensorflow models
path="datasets/tensorflow/models"
check_directories "$path"

# GROMACS adh_cubic
path="datasets/gromacs/adh_cubic"
check_directories "$path"


# nmonchart
cmd="sh nmonchart/nmonchart -h"
name="nmonchart"
check_tool "$cmd" "$name"

# GCC 7.3.0
cmd="gcc/gcc-installed/bin/gcc --version"
name="GCC 7.3.0"
check_tool "$cmd" "$name"

# BBMap
cmd="BBMap/bbmap/bbmap.sh"
name="BBMap"
check_tool "$cmd" "$name"

# Bowtie2 build
cmd="bowtie2/bowtie2-2.3.4.2/bowtie2-build --version"
name="Bowtie2 build"
check_tool "$cmd" "$name"

# Bowtie2 align
cmd="bowtie2/bowtie2-2.3.4.2/bowtie2 --version"
name="Bowtie2 align"
check_tool "$cmd" "$name"

# BWA
version_line=$(BWA/bwa-0.7.17/bwa 2>&1 | grep Version)
name="BWA"
if [ "$version_line" == "Version: 0.7.17-r1188" ]
then
        echo -e "${green}$name seems to be installed correctly.${nc}"
else
        echo -e "${red}Something seems to be wrong with $name. Please check and maybe reinstall again.${nc}"
fi

# Velveth
cmd="velvet/velveth"
name="Velveth"
check_tool "$cmd" "$name"

# Velvetg
first_line=$(velvet/velvetg | head -n 1)
name="Velvetg"
if [ "$first_line" == "velvetg - de Bruijn graph construction, error removal and repeat resolution" ]
then
        echo -e "${green}$name seems to be installed correctly.${nc}"
else
        echo -e "${red}Something seems to be wrong with $name. Please check and maybe reinstall again.${nc}"
fi

# IDBA
path="datasets/1000_genomes/ERR251006.filt.fa"
name="IDBA"
if [ -e "$path" ]
then
	echo -e "${green}$name seems to be installed correctly.${nc}"
else
	echo -e "${red}Something seems to be wrong with $name. Please check and maybe reinstall again.${nc}"
fi

# Tensorflow
cmd="python datasets/tensorflow/models/tutorials/image/cifar10/cifar10_train.py --help"
name="Tensorflow"
check_tool "$cmd" "$name"

# GROMACS
# Set GCC to 7.3.0
export PATH=$PWD/gcc/gcc-installed/bin:$PATH
export LD_LIBRARY_PATH=$PWD/gcc/gcc-installed/lib64:$LD_LIBRARY_PATH

cmd="/usr/local/gromacs/bin/gmx --version"
name="GROMACS"
check_tool "$cmd" "$name"

# Reset to system compiler
export PATH=$original_path_variable
export LD_LIBRARY_PATH=$original_ld_library_variable

# SPAdes
cmd="python SPAdes/SPAdes-3.12.0-Linux/bin/spades.py --help"
name="SPAdes"
check_tool "$cmd" "$name"

# ClustalOmega
cmd="clustalOmega/clustal-omega-1.2.4/bin/clustalo --version"
name="ClustalOmega"
check_tool "$cmd" "$name"

# MAFFT
version_line=$(MAFFT/mafft-7.475-with-extensions/bin/mafft --help 2>&1 | grep MAFFT)
name="MAFFT"
if [ "$version_line" == "  MAFFT v7.475 (2020/Nov/23)" ]
then
        echo -e "${green}$name seems to be installed correctly.${nc}"
else
        echo -e "${red}Something seems to be wrong with $name. Please check and maybe reinstall again.${nc}"
fi

# SINA
cmd="SINA/sina-1.7.2-linux/sina -V"
name="SINA"
check_tool "$cmd" "$name"

