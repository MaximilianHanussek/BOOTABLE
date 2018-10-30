#!/bin/bash

green='\033[0;32m'
red='\033[0;31m'
nc='\033[0m' # No Color

check_files() {
filepath=$1
filename=$(basename $1)

if [ -e $filepath ]
then
        echo "${green}The file: $filename is at the right place.${nc}"
else
        echo "${red}The file: $filename can not be found at the expected location. Please check and re-download the file.${nc}"
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

# ERR016155
path="File datasets/1000_genomes/ERR016155.filt.fastq"
check_files "$path"







# Bowtie2 build
cmd="bowtie2/bowtie2-2.3.4.2/bowtie2-build --version"
name="Bowtie2 build"
check_tool "$cmd" "$name"

#  Bowtie2 align
cmd="bowtie2/bowtie2-2.3.4.2/bowtie2 --version"
name="Bowtie2 align"
check_tool "$cmd" "$name"

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

# Tensorflow
cmd="python datasets/tensorflow/models/tutorials/image/cifar10/cifar10_train.py --help"
name="Tensorflow"

# GROMACS
cmd="/usr/local/gromacs/bin/gmx --version"
name="GROMACS"
check_tool "$cmd" "$name"

# SPAdes
cmd="python SPAdes/SPAdes-3.12.0-Linux/bin/spades.py --help"
name="SPAdes"
check_tool "$cmd" "$name"

