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

#wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -P datasets/1000_genomes/

wget https://s3.denbi.uni-tuebingen.de/max/GRCh38_full_analysis_set_plus_decoy_hla.fa -P datasets/1000_genomes/

#wget https://s3.denbi.uni-tuebingen.de/fb-test/tf_benchmark/cifar-10-binary.tar.gz -P datasets/tensorflow

wget https://s3.denbi.uni-tuebingen.de/max/cifar-10-binary.tar.gz -P datasets/tensorflow
tar -xf datasets/tensorflow/cifar-10-binary.tar.gz -C datasets/tensorflow/

wget https://s3.denbi.uni-tuebingen.de/max/models.tar.gz -P /datasets/tensorflow
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

# Compile and install GCC 7.3.0
cd gcc/gcc-build
../gcc-7.3.0/configure --enable-languages=c,c++ --disable-multilib --prefix=$PWD/../gcc-installed
make clean
make -j$(nproc) && make install
cd ../../

# Compile and install bowtie2 
cd bowtie2/bowtie2-2.3.4.2/
make clean
make -j$(nproc)
#sudo make static-libs && make STATIC_BUILD=1
cd ../../

# Compile and install velvet
cd velvet/
make clean
make -j$(nproc) 'OPENMP=1'
cd ..

# Compile and install IDBA
cd IDBA/idba_ud-1.0.9/
./configure
make clean
make -j$(nproc)
cd ../../

# Compile and install tensorflow (sudo)
sudo pip install --upgrade --force-reinstall pip==9.0.3
sudo pip install tensorflow==1.4.0

# Compile and install GROMACS (sudo)
mkdir gromacs/gromacs-2018.3/build
cd gromacs/gromacs-2018.3/build/
cmake3 -DGMX_BUILD_OWN_FFTW=on -DGMX_GPU=off -DGMX_BUILD_MPI=off --build ./  ../../gromacs-2018.3/
make clean
make -j$(nproc)
sudo make install
cd ../../../

/usr/local/gromacs/bin/gmx grompp -f datasets/gromacs/adh_cubic/pme_verlet.mdp -c datasets/gromacs/adh_cubic/conf.gro -p datasets/gromacs/adh_cubic/topol.top -o datasets/gromacs/adh_cubic/topol -po datasets/gromacs/adh_cubic/mdout

# Convert fastq files to fasta files
echo "Converting datasets from .fastq to .fa"

IDBA/idba_ud-1.0.9/bin/fq2fa datasets/1000_genomes/ERR016155.filt.fastq datasets/1000_genomes/ERR016155.filt.fa

IDBA/idba_ud-1.0.9/bin/fq2fa datasets/1000_genomes/ERR251006.filt.fastq datasets/1000_genomes/ERR251006.filt.fa
