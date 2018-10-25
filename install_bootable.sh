#!bin/bash

# Create empty directories
mkdir benchmark_output
mkdir benchmark_output/bowtie2
mkdir benchmark_output/gromacs
mkdir benchmark_output/IDBA
mkdir benchmark_output/SPAdes
mkdir benchmark_output/tensorflow
mkdir benchmark_output/velvet
mkdir datasets
mkdir datasets/1000_genomes
mkdir datasets/gromacs
mkdir datasets/tensorflow
mkdir dataset/ebi
mkdir gromacs
mkdir results
mkdir bowtie2
mkdir IDBA
mkdir SPAdes
mkdir tensorflow
mkdir velvet
mkdir results
mkdir gcc
mkdir gcc/gcc-build
mkdir gcc/gcc-installed

# Download benchmark datasets
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00110/sequence_read/ERR251006.filt.fastq.gz -P datasets/1000_genomes/

wget https://s3.denbi.uni-tuebingen.de/max/ERR251006.filt.fastq.gz -P datasets/1000_genomes/
gunzip datasets/1000_genomes/ERR251006.filt.fastq.gz

#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00125/sequence_read/ERR016155.filt.fastq.gz -P datasets/1000_genomes/

wget https://s3.denbi.uni-tuebingen.de/max/ERR016155.filt.fastq.gz -P datasets/1000_genomes/
gunzip datasets/1000_genomes/ERR016155.filt.fastq.gz

#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR001/DRR001012/DRR001012.fastq.gz -P datasets/ebi/

wget https://s3.denbi.uni-tuebingen.de/max/DRR001012.fastq.gz -p datasets/ebi/
gunzip datasets/ebi/DRR001012.fastq.gz

#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR001/DRR001025/DRR001025.fastq.gz -P datasets/ebi/

wget https://s3.denbi.uni-tuebingen.de/max/DRR001025.fastq.gz -p datasets/ebi/
gunzip datasets/ebi/DRR001025.fastq.gz

#wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -P datasets/1000_genomes/

wget https://s3.denbi.uni-tuebingen.de/max/GRCh38_full_analysis_set_plus_decoy_hla.fa -P datasets/1000_genomes/

wget https://s3.denbi.uni-tuebingen.de/max/cifar-10-binary.tar.gz -P datasets/tensorflow
tar -xf datasets/tensorflow/cifar-10-binary.tar.gz -C datasets/tensorflow/

wget https://s3.denbi.uni-tuebingen.de/max/models.tar.gz -P datasets/tensorflow
tar -xf datasets/tensorflow/models.tar.gz -C datasets/tensorflow/
rm datasets/tensorflow/models.tar.gz

wget https://s3.denbi.uni-tuebingen.de/max/ADH_bench_systems.tar.gz -P datasets/gromacs
tar -xf datasets/gromacs/ADH_bench_systems.tar.gz -C datasets/gromacs/
rm datasets/gromacs/ADH_bench_systems.tar.gz

# Download GCC compiler version 7.3.0
wget https://s3.denbi.uni-tuebingen.de/max/gcc-7.3.0.tar.gz -P gcc
tar -xf gcc/gcc-7.3.0.tar.gz -C gcc/
rm gcc/gcc-7.3.0.tar.gz

# Download Bowtie2 Sources
wget https://s3.denbi.uni-tuebingen.de/max/bowtie2-2.3.4.2-source.zip -P bowtie2
unzip bowtie2/bowtie2-2.3.4.2-source.zip -d bowtie2/
rm bowtie2/bowtie2-2.3.4.2-source.zip

# Download Velvet github repository
git clone https://github.com/dzerbino/velvet.git velvet/
rm -rf velvet/.gitignore velvet/.git/

# Download IDBA sources
wget https://s3.denbi.uni-tuebingen.de/max/idba_ud-1.0.9.tar.gz -P IDBA
tar -xf IDBA/idba_ud-1.0.9.tar.gz -C IDBA/
rm IDBA/idba_ud-1.0.9.tar.gz

# Download GROMACS sources
#wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-2018.3.tar.gz -P gromacs/

wget https://s3.denbi.uni-tuebingen.de/max/gromacs-2018.3.tar.gz -P gromacs/
tar -xf gromacs/gromacs-2018.3.tar.gz -C gromacs/
rm gromacs/gromacs-2018.3.tar.gz

# Download SPAdes binaries
wget https://s3.denbi.uni-tuebingen.de/max/SPAdes-3.12.0-Linux.tar.gz -P SPAdes
tar -xf SPAdes/SPAdes-3.12.0-Linux.tar.gz -C SPAdes
rm SPAdes/SPAdes-3.12.0-Linux.tar.gz

# Save original PATH and LD_LIBRARY variables
original_path_variable=$(echo $PATH)
original_ld_library_variable=$(echo $LD_LIBRARY_PATH)

# Compile and install GCC 7.3.0
cd gcc/gcc-7.3.0
./contrib/download_prerequisites
cd ../gcc-build
../gcc-7.3.0/configure --enable-languages=c,c++ --disable-multilib --prefix=$PWD/../gcc-installed
make clean
make -j$(nproc) && make install
cd ../../

# Compile and install bowtie2 
cd bowtie2/bowtie2-2.3.4.2/
make clean
sudo make static-libs && make -j$(nproc) STATIC_BUILD=1
cd ../../

# Set GCC to 7.3.0
export PATH=$PWD/gcc/gcc-installed/bin:$PATH
export LD_LIBRARY_PATH=$PWD/gcc/gcc-installed/lib64:$LD_LIBRARY_PATH

# Compile and install velvet
cd velvet/
make clean
make -j$(nproc) 'OPENMP=1'
cd ..

# Reset to system compiler as IDBA is not compiling with GCC 7.3.0
export PATH=$original_path_variable
export LD_LIBRARY_PATH=$original_ld_library_variable

# Compile and install IDBA
cd IDBA/idba_ud-1.0.9/
./configure
make clean
make -j$(nproc)
cd ../../

# Compile and install tensorflow (sudo)
sudo pip install --upgrade --force-reinstall pip==9.0.3
sudo pip install tensorflow==1.4.0

# Set GCC to 7.3.0 to fasten up GROMACS
export PATH=$PWD/gcc/gcc-installed/bin:$PATH
export LD_LIBRARY_PATH=$PWD/gcc/gcc-installed/lib64:$LD_LIBRARY_PATH

# Compile and install GROMACS (sudo)
mkdir gromacs/gromacs-2018.3/build
cd gromacs/gromacs-2018.3/build/
CC=gcc CXX=g++ cmake3 -DCMAKE_CXX_COMPILER=g++ -DGMX_BUILD_OWN_FFTW=on -DGMX_GPU=off -DGMX_BUILD_MPI=off --build ./ ../../gromacs-2018.3/
make clean
make -j$(nproc)
sudo make install
cd ../../../

/usr/local/gromacs/bin/gmx grompp -f datasets/gromacs/adh_cubic/pme_verlet.mdp -c datasets/gromacs/adh_cubic/conf.gro -p datasets/gromacs/adh_cubic/topol.top -o datasets/gromacs/adh_cubic/topol -po datasets/gromacs/adh_cubic/mdout

# Convert fastq files to fasta files
echo "Converting datasets from .fastq to .fa"

IDBA/idba_ud-1.0.9/bin/fq2fa datasets/1000_genomes/ERR016155.filt.fastq datasets/1000_genomes/ERR016155.filt.fa

IDBA/idba_ud-1.0.9/bin/fq2fa datasets/1000_genomes/ERR251006.filt.fastq datasets/1000_genomes/ERR251006.filt.fa

IDBA/idba_ud-1.0.9/bin/fq2fa datasets/ebi/DRR001012.fastq datasets/ebi/DRR001012.fa

IDBA/idba_ud-1.0.9/bin/fq2fa datasets/ebi/DRR001025.fastq datasets/ebi/DRR001025.fa

# Reset to system compiler after finishing installation process
export PATH=$original_path_variable
export LD_LIBRARY_PATH=$original_ld_library_variable

