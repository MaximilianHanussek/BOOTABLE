#!bin/bash

#!/bin/bash

progress_bar() {
  local duration
  local columns
  local space_available
  local fit_to_screen
  local space_reserved

  space_reserved=6   # reserved width for the percentage value
  duration=${1}
  columns=$(tput cols)
  space_available=$(( columns-space_reserved ))
  elapsed=${2}

  if (( duration < space_available )); then
        fit_to_screen=1;
  else
    fit_to_screen=$(( duration / space_available ));
    fit_to_screen=$((fit_to_screen+1));
  fi

  already_done() { for ((done=0; done<(elapsed / fit_to_screen) ; done=done+1 )); do printf "▇"; done }
  remaining() { for (( remain=(elapsed/fit_to_screen) ; remain<(duration/fit_to_screen) ; remain=remain+1 )); do printf " "; done }
  percentage() { echo -n "| $duration% "; }
  clean_line() { printf "\r"; }

  for (( elapsed; elapsed<=duration; elapsed=elapsed+1 )); do
      already_done; remaining; percentage
      clean_line
  done
  clean_line
}

# Create empty directories

progress_bar 0 0

if [ -d log ]
then
        echo "Directory log already exists."
else
        mkdir log
fi

if [ -d benchmark_output ]
then
        echo "Directory benchmark_output already exists."
else
        mkdir benchmark_output
fi

if [ -d benchmark_output/BBMap ]
then
        echo "Directory benchmark_output/BBMap already exists."
else
	mkdir benchmark_output/BBMap
fi

if [ -d benchmark_output/bowtie2 ]
then
        echo "Directory benchmark_output/bowtie2 already exists."
else
	mkdir benchmark_output/bowtie2
fi

if [ -d benchmark_output/BWA ]
then
        echo "Directory benchmark_output/BWA already exists."
else
        mkdir benchmark_output/BWA
fi

if [ -d benchmark_output/clustalOmega ]
then
        echo "Directory benchmark_output/clustalOmega already exists."
else
	mkdir benchmark_output/clustalOmega
fi

if [ -d benchmark_output/gromacs ]
then
        echo "Directory benchmark_output/gromacs already exists."
else
	mkdir benchmark_output/gromacs
fi

if [ -d benchmark_output/IDBA ]
then
        echo "Directory benchmark_output/IDBA already exists."
else
	mkdir benchmark_output/IDBA
fi

if [ -d benchmark_output/MAFFT ]
then
        echo "Directory benchmark_output/MAFFT already exists."
else
        mkdir benchmark_output/MAFFT
fi

if [ -d benchmark_output/SINA ]
then
        echo "Directory benchmark_output/SINA already exists."
else
        mkdir benchmark_output/SINA
fi

if [ -d benchmark_output/SPAdes ]
then
        echo "Directory benchmark_output/SPAdes already exists."
else
	mkdir benchmark_output/SPAdes
fi

if [ -d benchmark_output/tensorflow ]
then
        echo "Directory benchmark_output/tensorflow already exists."
else
	mkdir benchmark_output/tensorflow
fi

if [ -d benchmark_output/velvet ]
then
        echo "Directory benchmark_output/velvet already exists."
else
	mkdir benchmark_output/velvet
fi

if [ -d datasets ]
then
        echo "Directory datasets already exists."
else
	mkdir datasets
fi

if [ -d datasets/1000_genomes ]
then
        echo "Directory datasets/1000_genomes already exists."
else
	mkdir datasets/1000_genomes
fi

if [ -d datasets/BWA ]
then
        echo "Directory datasets/BWA already exists."
else
        mkdir datasets/BWA
fi

if [ -d datasets/BWA/DRR001012/ ]
then
        echo "Directory datasets/BWA/DRR001012/ already exists."
else
        mkdir datasets/BWA/DRR001012/
fi

if [ -d datasets/BWA/DRR001025/ ]
then
        echo "Directory datasets/BWA/DRR001025/ already exists."
else
        mkdir datasets/BWA/DRR001025/
fi

if [ -d datasets/BWA//GRCh38_full_analysis_set_plus_decoy_hla/ ]
then
        echo "Directory datasets/BWA//GRCh38_full_analysis_set_plus_decoy_hla/ already exists."
else
        mkdir datasets/BWA//GRCh38_full_analysis_set_plus_decoy_hla/
fi

if [ -d datasets/clustalOmega ]
then
        echo "Directory datasets/clustalOmega already exists."
else
	mkdir datasets/clustalOmega
fi

if [ -d datasets/gromacs ]
then
        echo "Directory datasets/gromacs already exists."
else
	mkdir datasets/gromacs
fi

if [ -d datasets/SINA ]
then
        echo "Directory datasets/SINA already exists."
else
        mkdir datasets/SINA
fi

if [ -d datasets/tensorflow ]
then
        echo "Directory datasets/tensorflow already exists."
else
	mkdir datasets/tensorflow
fi

if [ -d datasets/ebi ]
then
        echo "Directory datasets/ebi already exists."
else
	mkdir datasets/ebi
fi

if [ -d backed_up_benchmark_results ]
then
        echo "Directory backed_up_benchmark_results already exists."
else
	mkdir backed_up_benchmark_results
fi

if [ -d gromacs ]
then
        echo "Directory gromacs already exists."
else
	mkdir gromacs
fi

if [ -d results ]
then
        echo "Directory results already exists."
else
	mkdir results
fi

if [ -d BBMap ]
then
        echo "Directory BBMap already exists."
else
        mkdir BBMap
fi

if [ -d bowtie2 ]
then
        echo "Directory bowtie2 already exists."
else
	mkdir bowtie2
fi

if [ -d BWA ]
then
        echo "Directory BWA already exists."
else
        mkdir BWA
fi

if [ -d clustalOmega ]
then
        echo "Directory clustalOmega already exists."
else
	mkdir clustalOmega
fi

if [ -d IDBA ]
then
        echo "Directory IDBA already exists."
else
	mkdir IDBA
fi

if [ -d MAFFT ]
then
        echo "Directory MAFFT already exists."
else
        mkdir MAFFT
fi

if [ -d SINA ]
then
        echo "Directory SINA already exists."
else
        mkdir SINA
fi

if [ -d SPAdes ]
then
        echo "Directory SPAdes already exists."
else
	mkdir SPAdes
fi

if [ -d tensorflow ]
then
        echo "Directory tensorflow already exists."
else
	mkdir tensorflow
fi

if [ -d velvet ]
then
        echo "Directory velvet already exists."
else
	mkdir velvet
fi

if [ -d gcc ]
then
        echo "Directory gcc already exists."
else
	mkdir gcc
fi

if [ -d gcc/gcc-build ]
then
        echo "Directory gcc/gcc-build already exists."
else
	mkdir gcc/gcc-build
fi

if [ -d gcc/gcc-installed ]
then
        echo "Directory gcc/gcc-installed already exists."
else
	mkdir gcc/gcc-installed
fi

if [ -d nmon_stats ]
then
        echo "Directory nmon_stats already exists."
else
	mkdir nmon_stats
fi

progress_bar 1 0

# Download benchmark datasets

if [ -e datasets/1000_genomes/ERR251006.filt.fastq ]
then
	echo "File datasets/1000_genomes/ERR251006.filt.fastq already exists."
else
	#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00110/sequence_read/ERR251006.filt.fastq.gz -P datasets/1000_genomes/
	wget https://s3.denbi.uni-tuebingen.de/max/ERR251006.filt.fastq.gz -P datasets/1000_genomes/ >>log/dataset_ERR251006.log 2>&1
	gunzip datasets/1000_genomes/ERR251006.filt.fastq.gz >>log/dataset_ERR251006.log 2>&1
fi

if [ -e datasets/1000_genomes/ERR016155.filt.fastq ]
then
	echo "File datasets/1000_genomes/ERR016155.filt.fastq already exists."
else
	#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00125/sequence_read/ERR016155.filt.fastq.gz -P datasets/1000_genomes/
	wget https://s3.denbi.uni-tuebingen.de/max/ERR016155.filt.fastq.gz -P datasets/1000_genomes/ >>log/dataset_ERR016155.log 2>&1
	gunzip datasets/1000_genomes/ERR016155.filt.fastq.gz >>log/dataset_ERR016155.log 2>&1
fi

if [ -e datasets/1000_genomes/ERR015528.filt.fastq ]
then
        echo "File datasets/1000_genomes/ERR015528.filt.fastq already exists."
else
        #wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00106/sequence_read/ERR015528.filt.fastq.gz -P datasets/1000_genomes/
        wget https://s3.denbi.uni-tuebingen.de/max/ERR015528.filt.fastq.gz -P datasets/1000_genomes/ >>log/dataset_ERR015528.log 2>&1
        gunzip datasets/1000_genomes/ERR015528.filt.fastq.gz >>log/dataset_ERR015528.log 2>&1
fi

if [ -e datasets/1000_genomes/SRR741411.filt.fastq ]
then
        echo "File datasets/1000_genomes/SRR741411.filt.fastq already exists."
else
        #wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00099/sequence_read/SRR741411.filt.fastq.gz -P datasets/1000_genomes/
        wget https://s3.denbi.uni-tuebingen.de/max/SRR741411.filt.fastq.gz -P datasets/1000_genomes/ >>log/dataset_SRR741411.log 2>&1
        gunzip datasets/1000_genomes/SRR741411.filt.fastq.gz >>log/dataset_SRR741411.log 2>&1
fi

if [ -e datasets/ebi/DRR001012.fastq ]
then
        echo "File datasets/ebi/DRR001012.fastq already exists."
else
	#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR001/DRR001012/DRR001012.fastq.gz -P datasets/ebi/
	wget https://s3.denbi.uni-tuebingen.de/max/DRR001012.fastq.gz -P datasets/ebi/ >>log/dataset_DRR001012.log 2>&1
	gunzip datasets/ebi/DRR001012.fastq.gz >>log/dataset_DRR001012.log 2>&1
fi

if [ -e datasets/ebi/DRR001025.fastq ]
then
        echo "File datasets/ebi/DRR001025.fastq already exists."
else
	#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR001/DRR001025/DRR001025.fastq.gz -P datasets/ebi/
	wget https://s3.denbi.uni-tuebingen.de/max/DRR001025.fastq.gz -P datasets/ebi/ >>log/dataset_DRR001025.log 2>&1
	gunzip datasets/ebi/DRR001025.fastq.gz >>log/dataset_DRR001025.log 2>&1
fi

if [ -e datasets/1000_genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa ]
then
        echo "File datasets/1000_genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa already exists."
else
	#wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -P datasets/1000_genomes/
	wget https://s3.denbi.uni-tuebingen.de/max/GRCh38_full_analysis_set_plus_decoy_hla.fa -P datasets/1000_genomes/ >>log/dataset_GRCh38_full_analysis_set_plus_decoy_hla.log 2>&1
fi

if [ -e datasets/BWA/DRR001012/DRR001012.amb ]
then
        echo "File datasets/BWA/DRR001012/DRR001012.amb already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/DRR001012.amb -P datasets/BWA/DRR001012/ >>log/dataset_BWA_DRR001012_amb.log 2>&1
fi

if [ -e datasets/BWA/DRR001012/DRR001012.ann ]
then
        echo "File datasets/BWA/DRR001012/DRR001012.ann already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/DRR001012.ann -P datasets/BWA/DRR001012/ >>log/dataset_BWA_DRR001012_ann.log 2>&1
fi

if [ -e datasets/BWA/DRR001012/DRR001012.bwt ]
then
        echo "File datasets/BWA/DRR001012/DRR001012.bwt already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/DRR001012.bwt -P datasets/BWA/DRR001012/ >>log/dataset_BWA_DRR001012_bwt.log 2>&1
fi

if [ -e datasets/BWA/DRR001012/DRR001012.pac ]
then
        echo "File datasets/BWA/DRR001012/DRR001012.pac already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/DRR001012.pac -P datasets/BWA/DRR001012/ >>log/dataset_BWA_DRR001012_pac.log 2>&1
fi

if [ -e datasets/BWA/DRR001012/DRR001012.sa ]
then
        echo "File datasets/BWA/DRR001012/DRR001012.sa already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/DRR001012.sa -P datasets/BWA/DRR001012/ >>log/dataset_BWA_DRR001012_sa.log 2>&1
fi

if [ -e datasets/BWA/DRR001025/DRR001025.amb ]
then
        echo "File datasets/BWA/DRR001025/DRR001025.amb already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/DRR001025.amb -P datasets/BWA/DRR001025/ >>log/dataset_BWA_DRR001025_amb.log 2>&1
fi

if [ -e datasets/BWA/DRR001025/DRR001025.ann ]
then
        echo "File datasets/BWA/DRR001025/DRR001025.ann already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/DRR001025.ann -P datasets/BWA/DRR001025/ >>log/dataset_BWA_DRR001025_ann.log 2>&1
fi

if [ -e datasets/BWA/DRR001025/DRR001025.bwt ]
then
        echo "File datasets/BWA/DRR001025/DRR001025.bwt already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/DRR001025.bwt -P datasets/BWA/DRR001025/ >>log/dataset_BWA_DRR001025_bwt.log 2>&1
fi

if [ -e datasets/BWA/DRR001025/DRR001025.pac ]
then
        echo "File datasets/BWA/DRR001025/DRR001025.pac already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/DRR001025.pac -P datasets/BWA/DRR001025/ >>log/dataset_BWA_DRR001025_pac.log 2>&1
fi

if [ -e datasets/BWA/DRR001025/DRR001025.sa ]
then
        echo "File datasets/BWA/DRR001025/DRR001025.sa already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/DRR001025.sa -P datasets/BWA/DRR001025/ >>log/dataset_BWA_DRR001025_sa.log 2>&1
fi

if [ -e datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.amb ]
then
        echo "File datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.amb already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/GRCh38_full_analysis_set_plus_decoy_hla.amb -P datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/ >>log/dataset_BWA_GRCh38_full_analysis_set_plus_decoy_hla_amb.log 2>&1
fi

if [ -e datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.ann ]
then
        echo "File datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.ann already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/GRCh38_full_analysis_set_plus_decoy_hla.ann -P datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/ >>log/dataset_BWA_GRCh38_full_analysis_set_plus_decoy_hla_ann.log 2>&1
fi

if [ -e datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.bwt ]
then
        echo "File datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.bwt already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/GRCh38_full_analysis_set_plus_decoy_hla.bwt -P datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/ >>log/dataset_BWA_GRCh38_full_analysis_set_plus_decoy_hla_bwt.log 2>&1
fi

if [ -e datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.pac ]
then
        echo "File datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.pac already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/GRCh38_full_analysis_set_plus_decoy_hla.pac -P datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/ >>log/dataset_BWA_GRCh38_full_analysis_set_plus_decoy_hla_pac.log 2>&1
fi

if [ -e datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.sa ]
then
        echo "File datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.sa already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/GRCh38_full_analysis_set_plus_decoy_hla.sa -P datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/ >>log/dataset_BWA_GRCh38_full_analysis_set_plus_decoy_hla_sa.log 2>&1
fi

if [ -e datasets/clustalOmega/wgs.ANCA.1_200.fsa ]
then
        echo "File datasets/clustalOmega/wgs.ANCA.1_200.fsa already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/wgs.ANCA.1_200.fsa -P datasets/clustalOmega >>log/dataset_wgs.ANCA.1_200.log 2>&1
fi


if [ -e datasets/clustalOmega/wgs.ANCA.1_400.fsa ]
then
        echo "File datasets/clustalOmega/wgs.ANCA.1_400.fsa already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/wgs.ANCA.1_400.fsa -P datasets/clustalOmega >>log/dataset_wgs.ANCA.1_400.log 2>&1
fi


if [ -e datasets/clustalOmega/wgs.ANCA.1_500.fsa ]
then
        echo "File datasets/clustalOmega/wgs.ANCA.1_500.fsa already exists."
else
	wget https://s3.denbi.uni-tuebingen.de/max/wgs.ANCA.1_500.fsa -P datasets/clustalOmega >>log/dataset_wgs.ANCA.1_500.log 2>&1
fi


if [ -d datasets/tensorflow/cifar-10-batches-bin ]
then
	echo "Directory datasets/tensorflow/cifar-10-batches-bin already exists."
else	
	wget https://s3.denbi.uni-tuebingen.de/max/cifar-10-binary.tar.gz -P datasets/tensorflow >>log/dataset_cifar-10-batches-bin.log 2>&1
	tar -xf datasets/tensorflow/cifar-10-binary.tar.gz -C datasets/tensorflow/ >>log/dataset_cifar-10-batches-bin.log 2>&1
fi

if [ -d datasets/tensorflow/models ]
then
        echo "Directory datasets/tensorflow/models already exists."
else    
	wget https://s3.denbi.uni-tuebingen.de/max/models.tar.gz -P datasets/tensorflow >>log/dataset_tensorflow_models.log 2>&1
	tar -xf datasets/tensorflow/models.tar.gz -C datasets/tensorflow/ >>log/dataset_tensorflow_models.log 2>&1
	rm datasets/tensorflow/models.tar.gz >>log/dataset_tensorflow_models.log 2>&1
fi

if [ -e datasets/SINA/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb ]
then
        echo "File datasets/SINA/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb already exists."
else
	#wget https://www.arb-silva.de/fileadmin/arb_web_db/release_138_1/ARB_files/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb.gz
        wget https://s3.denbi.uni-tuebingen.de/max/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb.gz -P datasets/SINA/ >>log/dataset_SINA_SILVA_138.1_SSURef_NR99_12_06_20_opt.arb.log 2>&1
	gunzip datasets/SINA/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb.gz >>log/dataset_SINA_SILVA_138.1_SSURef_NR99_12_06_20_opt.arb.log 2>&1
fi

if [ -e datasets/SINA/SILVA_138.1_SSURef_NR99_12_06_20_opt.sidx ]
then
        echo "File datasets/SINA/SILVA_138.1_SSURef_NR99_12_06_20_opt.sidx already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/SILVA_138.1_SSURef_NR99_12_06_20_opt.sidx -P datasets/SINA/ >>log/dataset_SINA_SILVA_138.1_SSURef_NR99_12_06_20_opt.sidx.log 2>&1
fi

if [ -e datasets/SINA/GTDB_bac-arc_ssu_r86.fa ]
then
        echo "File datasets/SINA/GTDB_bac-arc_ssu_r86.fa already exists."
else 
	#wget https://zenodo.org/record/2541239/files/GTDB_bac-arc_ssu_r86.fa.gz
        wget https://s3.denbi.uni-tuebingen.de/max/GTDB_bac-arc_ssu_r86.fa.gz -P datasets/SINA/ >>log/dataset_SINA_GTDB_bac-arc_ssu_r86.log 2>&1
	gunzip datasets/SINA/GTDB_bac-arc_ssu_r86.fa.gz >>log/dataset_SINA_GTDB_bac-arc_ssu_r86.log 2>&1
fi

if [ -e datasets/SINA/RefSeq-RDP16S_v2_May2018.fa ]
then
        echo "File datasets/SINA/RefSeq-RDP16S_v2_May2018.fa already exists."
else
        #wget https://zenodo.org/record/2541239/files/RefSeq-RDP16S_v2_May2018.fa.gz
        wget https://s3.denbi.uni-tuebingen.de/max/RefSeq-RDP16S_v2_May2018.fa.gz -P datasets/SINA/ >>log/dataset_SINA_RefSeq-RDP16S_v2_May2018.log 2>&1
        gunzip datasets/SINA/RefSeq-RDP16S_v2_May2018.fa.gz >>log/dataset_SINA_RefSeq-RDP16S_v2_May2018.log 2>&1
fi

if [ -e datasets/SINA/OE-38_R1.fastq ]
then
        echo "File datasets/SINA/OE-38_R1.fastq already exists."
else
        #wget https://zenodo.org/record/803376/files/OE-38_R1.fastq.gz
        wget https://s3.denbi.uni-tuebingen.de/max/OE-38_R1.fastq.gz -P datasets/SINA/ >>log/dataset_SINA_OE-38_R1.log 2>&1
        gunzip datasets/SINA/OE-38_R1.fastq.gz >>log/dataset_SINA_OE-38_R1.log 2>&1
fi

if [ -d datasets/gromacs/adh_cubic ]
then
        echo "Directory datasets/gromacs/adh_cubic already exists."
else
	wget https://s3.denbi.uni-tuebingen.de/max/ADH_bench_systems.tar.gz -P datasets/gromacs >>log/dataset_gromacs_adh_cubic.log 2>&1
	tar -xf datasets/gromacs/ADH_bench_systems.tar.gz -C datasets/gromacs/ >>log/dataset_gromacs_adh_cubic.log 2>&1
	rm datasets/gromacs/ADH_bench_systems.tar.gz >>log/dataset_gromacs_adh_cubic.log 2>&1
fi

progress_bar 8 1

# Download GCC compiler version 7.3.0
if [ -d gcc/gcc-7.3.0 ]
then
	echo "Directory gcc/gcc-7.3.0 already exists."
else
	wget https://s3.denbi.uni-tuebingen.de/max/gcc-7.3.0.tar.gz -P gcc >>log/download_gcc.log 2>&1
	tar -xf gcc/gcc-7.3.0.tar.gz -C gcc/ >>log/download_gcc.log 2>&1
	rm gcc/gcc-7.3.0.tar.gz >>log/download_gcc.log 2>&1
fi

# Download nmonchart34
if ! [ -d nmonchart ]
then
	wget https://s3.denbi.uni-tuebingen.de/max/nmonchart34.tar >>log/download_nmonchart34.log 2>&1
	tar -xf nmonchart34.tar >>log/download_nmonchart34.log 2>&1
	rm nmonchart34.tar >>log/download_nmonchart34.log 2>&1
else
	echo "Directory nmonchart already exist."
fi

# Download BBMap binaries
if [ -d BBMap/bbmap ]
then
        echo "Directory BBMap/bbmap already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/BBMap_38.87.tar.gz -P BBMap >>log/download_BBMap.log 2>&1
        tar -xf BBMap/BBMap_38.87.tar.gz -C BBMap/ >>log/download_BBMap.log 2>&1
        rm BBMap/BBMap_38.87.tar.gz >>log/download_BBMap.log 2>&1
fi

# Download Bowtie2 sources
if [ -d bowtie2/bowtie2-2.3.4.2 ]
then
        echo "Directory bowtie2/bowtie2-2.3.4.2 already exists."
else
	wget https://s3.denbi.uni-tuebingen.de/max/bowtie2-2.3.4.2-source.zip -P bowtie2 >>log/download_bowtie2.log 2>&1
	unzip bowtie2/bowtie2-2.3.4.2-source.zip -d bowtie2/ >>log/download_bowtie2.log 2>&1
	rm bowtie2/bowtie2-2.3.4.2-source.zip >>log/download_bowtie2.log 2>&1
fi

# Download BWA sources
if [ -d BWA/bwa-0.7.17 ]
then
        echo "Directory BWA/bwa-0.7.17 already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/bwa-0.7.17.tar.bz2 -P BWA >>log/download_BWA.log 2>&1
        tar -xf BWA/bwa-0.7.17.tar.bz2 -C BWA/ >>log/download_BWA.log 2>&1
        rm BWA/bwa-0.7.17.tar.bz2 >>log/download_BWA.log 2>&1
fi

# Download Velvet github repository
if ! [ -z "$(ls -A "velvet")" ]
then
    echo "Directory velvet/ exists and is not empty."
else
	git clone https://github.com/dzerbino/velvet.git velvet/ >>log/download_velvet.log 2>&1
	rm -rf velvet/.gitignore velvet/.git/ >>log/download_velvet.log 2>&1
fi

# Download IDBA sources
if [ -d IDBA/idba_ud-1.0.9 ]
then
	echo "Directory IDBA/idba_ud-1.0.9 already exists."
else
	wget https://s3.denbi.uni-tuebingen.de/max/idba_ud-1.0.9.tar.gz -P IDBA >>log/download_IDBA.log 2>&1
	tar -xf IDBA/idba_ud-1.0.9.tar.gz -C IDBA/ >>log/download_IDBA.log 2>&1
	rm IDBA/idba_ud-1.0.9.tar.gz >>log/download_IDBA.log 2>&1
fi

# Download GROMACS sources
if [ -d gromacs/gromacs-2018.3 ]
then
	echo "Directory gromacs/gromacs-2018.3 already exists."
else
	#wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-2018.3.tar.gz -P gromacs/
	wget https://s3.denbi.uni-tuebingen.de/max/gromacs-2018.3.tar.gz -P gromacs/ >>log/download_GROMACS.log 2>&1
	tar -xf gromacs/gromacs-2018.3.tar.gz -C gromacs/ >>log/download_GROMACS.log 2>&1
	rm gromacs/gromacs-2018.3.tar.gz >>log/download_GROMACS.log 2>&1
fi

# Download SPAdes binaries
if [ -d SPAdes/SPAdes-3.12.0-Linux ]
then
        echo "Directory SPAdes/SPAdes-3.12.0-Linux already exists."
else
	wget https://s3.denbi.uni-tuebingen.de/max/SPAdes-3.12.0-Linux.tar.gz -P SPAdes >>log/download_SPAdes.log 2>&1
	tar -xf SPAdes/SPAdes-3.12.0-Linux.tar.gz -C SPAdes >>log/download_SPAdes.log 2>&1
	rm SPAdes/SPAdes-3.12.0-Linux.tar.gz >>log/download_SPAdes.log 2>&1
fi

# Download ClustalOmega sources
if [ -d clustalOmega/clustal-omega-1.2.4 ]
then
        echo "Directory clustalOmega/clustal-omega-1.2.4 already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/clustal-omega-1.2.4.tar.gz -P clustalOmega >>log/download_clustalOmega.log 2>&1
        tar -xf clustalOmega/clustal-omega-1.2.4.tar.gz -C clustalOmega/ >>log/download_clustalOmega.log 2>&1
        rm clustalOmega/clustal-omega-1.2.4.tar.gz >>log/download_clustalOmega.log 2>&1
fi

# Download MAFFT sources
if [ -d MAFFT/mafft-7.475-with-extensions-src ]
then
        echo "Directory MAFFT/mafft-7.475-with-extensions-src already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/mafft-7.475-with-extensions-src.tgz -P MAFFT >>log/download_MAFFT.log 2>&1
        tar -xf MAFFT/mafft-7.475-with-extensions-src.tgz -C MAFFT/ >>log/download_MAFFT.log 2>&1
        rm MAFFT/mafft-7.475-with-extensions-src.tgz >>log/download_MAFFT.log 2>&1
fi

# Download SINA binaries
if [ -d SINA/sina-1.7.2-linux ]
then
        echo "Directory SINA/sina-1.7.2-linux already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/sina-1.7.2-linux.tar.gz -P SINA >>log/download_SINA.log 2>&1
        tar -xf SINA/sina-1.7.2-linux.tar.gz -C SINA/ >>log/download_SINA.log 2>&1
        rm SINA/sina-1.7.2-linux.tar.gz >>log/download_SINA.log 2>&1
fi


progress_bar 11 8

# Save original PATH and LD_LIBRARY variables
original_path_variable=$(echo $PATH)
original_ld_library_variable=$(echo $LD_LIBRARY_PATH)

 Compile and install GCC 7.3.0
if [ -e gcc/gcc-installed/bin/gcc ]
then 
	while true; do
        	read -p "GCC 7.3.0 seems already to be installed, do you want to recompile it?" yn
                case $yn in
                        [Yy]* ) rm -rf gcc/gcc-build/* >>log/install_gcc.log 2>&1
				rm -rf gcc/gcc-installed/* >>log/install_gcc.log 2>&1
				cd gcc/gcc-7.3.0
                                ./contrib/download_prerequisites >>../../log/install_gcc.log 2>&1
                                cd ../gcc-build
                                ../gcc-7.3.0/configure --enable-languages=c,c++ --disable-multilib --prefix=$PWD/../gcc-installed >>../../log/install_gcc.log 2>&1
                                make clean >>../../log/install_gcc.log 2>&1
                                make -j$(nproc) >>../../log/install_gcc.log 2>&1
				make install >>../../log/install_gcc.log 2>&1
                                cd ../../
                                break;;
                        [Nn]* ) echo "GCC will not be recompiled"
                                break;;
                        * ) echo "Please answer yes or no."
				;;
                esac
        done
else
        cd gcc/gcc-7.3.0
        ./contrib/download_prerequisites >>../../log/install_gcc.log 2>&1
        cd ../gcc-build
        ../gcc-7.3.0/configure --enable-languages=c,c++ --disable-multilib --prefix=$PWD/../gcc-installed >>../../log/install_gcc.log 2>&1
        make clean >>../../log/install_gcc.log 2>&1
        make -j$(nproc) >>../../log/install_gcc.log 2>&1
	make install >>../../log/install_gcc.log 2>&1
        cd ../../
fi

progress_bar 75 11


# Compile and install BWA
if [ -e BWA/bwa-0.7.17/bwa ]
then
        while true; do
                read -p "BWA seems already to be installed, do you want to recompile it?" yn
                case $yn in
                        [Yy]* ) cd BWA/bwa-0.7.17/
                                make clean >>../../log/install_BWA.log 2>&1
                                make -j$(nproc) >>../../log/install_BWA.log 2>&1
                                cd ../../
                                break;;
                        [Nn]* ) echo "BWA will not be recompiled"
                                break;;
                        * ) echo "Please answer yes or no."
                                ;;
                esac
        done
else
        cd BWA/bwa-0.7.17/
        make clean >>../../log/install_BWA.log 2>&1
        make -j$(nproc) >>../../log/install_bowtie2.log 2>&1
        cd ../../
fi

# Compile and install bowtie2
if [ -e bowtie2/bowtie2-2.3.4.2/bowtie2-align-l ]
then      
        while true; do
                read -p "Bowtie2 seems already to be installed, do you want to recompile it?" yn
                case $yn in
                        [Yy]* )	cd bowtie2/bowtie2-2.3.4.2/
        			make clean >>../../log/install_bowtie2.log 2>&1
        			sudo make static-libs >>../../log/install_bowtie2.log 2>&1
				make -j$(nproc) STATIC_BUILD=1 >>../../log/install_bowtie2.log 2>&1
        			cd ../../
                                break;;
                        [Nn]* ) echo "Bowtie2 will not be recompiled"
                                break;;
                        * ) echo "Please answer yes or no."
                                ;;
                esac
        done
else
	cd bowtie2/bowtie2-2.3.4.2/
	make clean >>../../log/install_bowtie2.log 2>&1
	sudo make static-libs >>../../log/install_bowtie2.log 2>&1
	make -j$(nproc) STATIC_BUILD=1 >>../../log/install_bowtie2.log 2>&1
	cd ../../
fi

progress_bar 80 75

# Set GCC to 7.3.0
export PATH=$PWD/gcc/gcc-installed/bin:$PATH
export LD_LIBRARY_PATH=$PWD/gcc/gcc-installed/lib64:$LD_LIBRARY_PATH

# Compile and install velvet
if [ -e velvet/velvetg ] && [ -e velvet/velveth ]
then
        while true; do
                read -p "Velvet seems already to be installed, do you want to recompile it?" yn
                case $yn in
                        [Yy]* ) cd velvet/
        			make clean >>../log/install_velvet.log 2>&1
        			make -j$(nproc) 'OPENMP=1' >>../log/install_velvet.log 2>&1
        			cd ..
                                break;;
                        [Nn]* ) echo "Velvet will not be recompiled"
                                break;;
                        * ) echo "Please answer yes or no."
                                ;;
                esac
        done
else
	cd velvet/
	make clean >>../log/install_velvet.log 2>&1
	make -j$(nproc) 'OPENMP=1' >>../log/install_velvet.log 2>&1
	cd ..
fi

progress_bar 81 80

# Reset to system compiler as IDBA is not compiling with GCC 7.3.0
export PATH=$original_path_variable
export LD_LIBRARY_PATH=$original_ld_library_variable


# Compile and install MAFFT
if [ -e MAFFT/mafft-7.475-with-extensions/bin/mafft ]
then
        while true; do
                read -p "MAFFT seems already to be installed, do you want to recompile it?" yn
                case $yn in
                        [Yy]* ) sed -i "/PREFIX = /c\PREFIX = $PWD/MAFFT/mafft-7.475-with-extensions" MAFFT/mafft-7.475-with-extensions/core/Makefile
				cd MAFFT/mafft-7.475-with-extensions/core/
                                make clean >>../../../log/install_MAFFT.log 2>&1
                                make -j$(nproc) >>../../../log/install_MAFFT.log 2>&1
				make install >>../../../log/install_MAFFT.log 2>&1
                                cd ../../../
                                break;;
                        [Nn]* ) echo "MAFFT will not be recompiled"
                                break;;
                        * ) echo "Please answer yes or no."
                                ;;
                esac
        done
else
	sed -i "/PREFIX = /c\PREFIX = $PWD/MAFFT/mafft-7.475-with-extensions" MAFFT/mafft-7.475-with-extensions/core/Makefile
        cd MAFFT/mafft-7.475-with-extensions/core/
        make clean >>../../../log/install_MAFFT.log 2>&1
        make -j$(nproc) >>../../../log/install_MAFFT.log 2>&1
	make install >>../../../log/install_MAFFT.log 2>&1
        cd ../../../
fi

# Compile and install IDBA
if [ -e IDBA/idba_ud-1.0.9/bin/idba_ud ]
then
        while true; do
                read -p "IDBA seems already to be installed, do you want to recompile it?" yn
                case $yn in
                        [Yy]* ) cd IDBA/idba_ud-1.0.9/
				./configure >>../../log/install_IDBA.log 2>&1
				make clean >>../../log/install_IDBA.log 2>&1
				make -j$(nproc) >>../../log/install_IDBA.log 2>&1
				cd ../../
                                break;;
                        [Nn]* ) echo "IDBA will not be recompiled"
                                break;;
                        * ) echo "Please answer yes or no."
                                ;;
                esac
        done
else
	cd IDBA/idba_ud-1.0.9/
	./configure >>../../log/install_IDBA.log 2>&1
	make clean >>../../log/install_IDBA.log 2>&1
	make -j$(nproc) >>../../log/install_IDBA.log 2>&1
	cd ../../
fi

progress_bar 83 81


# Compile and install clustalOmega
if [ -e clustalOmega/clustal-omega-1.2.4/bin/clustalo ]
then
	while true; do
                read -p "clustalOmega seems already to be installed, do you want to recompile it?" yn
                case $yn in
                        [Yy]* ) cd clustalOmega/clustal-omega-1.2.4/
                                ./configure --prefix=$PWD >>../../log/install_clustalOmega.log 2>&1
                                make clean >>../../log/install_clustalOmega.log 2>&1
                                make -j$(nproc) >>../../log/install_clustalOmega.log 2>&1
				make install >>../../log/install_clustalOmega.log 2>&1
                                cd ../../
                                break;;
                        [Nn]* ) echo "ClustalOmega will not be recompiled"
                                break;;
                        * ) echo "Please answer yes or no."
                                ;;
                esac
        done
else
	cd clustalOmega/clustal-omega-1.2.4/
        ./configure --prefix=$PWD >>../../log/install_clustalOmega.log 2>&1
        make clean >>../../log/install_clustalOmega.log 2>&1
        make -j$(nproc) >>../../log/install_clustalOmega.log 2>&1
	make install >>../../log/install_clustalOmega.log 2>&1
        cd ../../
fi

progress_bar 86 83


# Compile and install tensorflow (sudo)
#check if pip is in general installed (if 0 then it is installed if 1 needs to be installed)
which pip > /dev/null 2>&1
pip_installed=$(echo $?)
if [ $pip_installed == 0 ]
then
	echo "pip is already installed, nothing to do here." >>log/install_pip.log 2>&1
else
	while true; do
                read -p "pip in general seems not to be installed or is not in your path, do you want to install it via yum?" yn
                case $yn in
                        [Yy]* ) sudo yum install pip -y >>log/install_pip.log 2>&1
                                break;;
                        [Nn]* ) echo "pip will not be installed, please install it on your own, or the Tensorflow benchmarks will not work"
                                break;;
                        * ) echo "Please answer yes or no."
                                ;;
                esac
        done
fi


#check if pip is already installed in Version 9.0.3 (if 1 then it is installed if 0 needs to be installed)
which pip > /dev/null 2>&1
pip_installed=$(echo $?)
pip_version_installed=$(pip list installed  2> /dev/null | grep -c "pip (9.0.3)")
if [ $pip_installed == 0 ] && [ $pip_version_installed == 0 ]
then
	while true; do
                read -p "pip is already installed but not in the required version 9.0.3 do you want to upgrade/downgrade to pip version 9.0.3?" yn
                case $yn in
                        [Yy]* ) sudo pip install --upgrade --force-reinstall pip==9.0.3 >>log/install_pip.log 2>&1
                                break;;
                        [Nn]* ) echo "pip will not be upgraded/downgraded to version 9.0.3, make sure tensorflow version 1.4.0 is working correctly"
                                break;;
                        * ) echo "Please answer yes or no."
                                ;;
                esac
        done
else
	echo "pip is already installed in version 9.0.3, nothing to do here." >>log/install_pip.log 2>&1
fi

#check if tensorflow is installed in version 1.4.0 (if 1 then it is installed if 0 needs to be installed)
pip_version_installed=$(pip list installed  2> /dev/null | grep -c "pip (9.0.3)")
tensorflow_installed=$(pip list installed  2> /dev/null | grep -c "tensorflow (1.4.0)")

if [ $tensorflow_installed == 0 ] && [ $pip_version_installed == 1 ]
then
	sudo pip install tensorflow==1.4.0 >>log/install_tensorflow.log 2>&1
elif [ $tensorflow_installed == 1 ] && [ $pip_version_installed == 1 ]
then 
	while true; do
                read -p "Tensorflow seems already to be installed in version 1.4.0 via pip, do you want to reinstall it?" yn
                case $yn in
                        [Yy]* ) sudo pip install --force-reinstall tensorflow==1.4.0 >>log/install_tensorflow.log 2>&1
                                break;;
                        [Nn]* ) echo "Tensorflow will not be reinstalled"
                                break;;
                        * ) echo "Please answer yes or no."
                                ;;
                esac
        done
else 
	echo "Something seems to be wrong with the pip and tensorflow installation please check if pip is installed in version 9.0.3 and tensorflow in version 1.4.0."
fi


progress_bar 88 86

# Set GCC to 7.3.0 to fasten up GROMACS
export PATH=$PWD/gcc/gcc-installed/bin:$PATH
export LD_LIBRARY_PATH=$PWD/gcc/gcc-installed/lib64:$LD_LIBRARY_PATH

# Compile and install GROMACS (sudo)
if [ -e /usr/local/gromacs/bin/gmx ]
then
        while true; do
                read -p "GROMACS seems already to be installed, do you want to recompile it?" yn
                case $yn in
                        [Yy]* ) rm -rf gromacs/gromacs-2018.3/build/* >>log/install_GROMACS.log 2>&1
				sudo rm -rf /usr/local/gromacs/ >>log/install_GROMACS.log 2>&1
				cd gromacs/gromacs-2018.3/build/
				CC=gcc CXX=g++ cmake3 -DCMAKE_CXX_COMPILER=g++ -DGMX_BUILD_OWN_FFTW=on -DGMX_GPU=off -DGMX_BUILD_MPI=off --build ./ ../../gromacs-2018.3/ >>../../../log/install_GROMACS.log 2>&1
				make clean >>../../../log/install_GROMACS.log 2>&1
				make -j$(nproc) >>../../../log/install_GROMACS.log 2>&1
				sudo make install >>../../../log/install_GROMACS.log 2>&1
				cd ../../../
				break;;
                        [Nn]* ) echo "GROMACS will not be recompiled"
                                break;;
                        * ) echo "Please answer yes or no."
                                ;;
                esac
        done
else
	mkdir gromacs/gromacs-2018.3/build
	cd gromacs/gromacs-2018.3/build/
	CC=gcc CXX=g++ cmake3 -DCMAKE_CXX_COMPILER=g++ -DGMX_BUILD_OWN_FFTW=on -DGMX_GPU=off -DGMX_BUILD_MPI=off --build ./ ../../gromacs-2018.3/ >>../../../log/install_GROMACS.log 2>&1
	make clean >>../../../log/install_GROMACS.log 2>&1
	make -j	$(nproc) >>../../../log/install_GROMACS.log 2>&1
	sudo make install >>../../../log/install_GROMACS.log 2>&1
	cd ../../../
	#/usr/local/gromacs/bin/gmx grompp -f datasets/gromacs/adh_cubic/pme_verlet.mdp -c datasets/gromacs/adh_cubic/conf.gro -p datasets/gromacs/adh_cubic/topol.top -o datasets/gromacs/adh_cubic/topol -po datasets/gromacs/adh_cubic/mdout
#fi

progress_bar 95 88


# Convert fastq files to fasta files
echo "Converting datasets from .fastq to .fa" >>log/file_conversion_ERR016155.log 2>&1

if [ -e datasets/1000_genomes/ERR016155.filt.fa ]
then
	echo "File datasets/1000_genomes/ERR016155.filt.fa already exists."
else
	IDBA/idba_ud-1.0.9/bin/fq2fa datasets/1000_genomes/ERR016155.filt.fastq datasets/1000_genomes/ERR016155.filt.fa >>log/file_conversion_ERR016155.log 2>&1
fi

if [ -e datasets/1000_genomes/ERR251006.filt.fa ]
then
        echo "File datasets/1000_genomes/ERR251006.filt.fa already exists."
else
	IDBA/idba_ud-1.0.9/bin/fq2fa datasets/1000_genomes/ERR251006.filt.fastq datasets/1000_genomes/ERR251006.filt.fa >>log/file_conversion_ERR251006.log 2>&1
fi

if [ -e datasets/1000_genomes/ERR015528.filt.fa ]
then
        echo "File datasets/1000_genomes/ERR015528.filt.fa already exists."
else
        IDBA/idba_ud-1.0.9/bin/fq2fa datasets/1000_genomes/ERR015528.filt.fastq datasets/1000_genomes/ERR015528.filt.fa >>log/file_conversion_ERR015528.log 2>&1
fi

if [ -e datasets/1000_genomes/SRR741411.filt.fa ]
then
        echo "File datasets/1000_genomes/SRR741411.filt.fa already exists."
else
        IDBA/idba_ud-1.0.9/bin/fq2fa datasets/1000_genomes/SRR741411.filt.fastq datasets/1000_genomes/SRR741411.filt.fa >>log/file_conversion_SRR741411.log 2>&1
fi

if [ -e datasets/ebi/DRR001012.fa ]
then
        echo "File datasets/ebi/DRR001012.fa already exists."
else
	IDBA/idba_ud-1.0.9/bin/fq2fa datasets/ebi/DRR001012.fastq datasets/ebi/DRR001012.fa >>log/file_conversion_DRR001012.log 2>&1
fi

if [ -e datasets/ebi/DRR001025.fa ]
then
        echo "File datasets/ebi/DRR001025.fa already exists."
else
	IDBA/idba_ud-1.0.9/bin/fq2fa datasets/ebi/DRR001025.fastq datasets/ebi/DRR001025.fa >>log/file_conversion_DRR001025.log 2>&1
fi

if [ -e datasets/SINA/OE-38_R1.fa ]
then
        echo "File datasets/SINA/OE-38_R1.fa already exists."
else
        IDBA/idba_ud-1.0.9/bin/fq2fa datasets/SINA/OE-38_R1.fastq datasets/SINA/OE-38_R1.fa >>log/file_conversion_DRR001025.log 2>&1
fi


# Reset to system compiler after finishing installation process
export PATH=$original_path_variable
export LD_LIBRARY_PATH=$original_ld_library_variable

progress_bar 100 95
echo "BOOTABLE should be successfully installed, please check with the install_check.sh tool."
