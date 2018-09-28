#!bin/bash

# Download benchmark datasets
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00110/sequence_read/ERR251006.filt.fastq.gz -P datasets/1000_genomes/
gunzip datasets/1000_genomes/ERR251006.filt.fastq.gz

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00125/sequence_read/ERR016155.filt.fastq.gz -P datasets/1000_genomes/
gunzip datasets/1000_genomes/ERR016155.filt.fastq.gz

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -P datasets/1000_genomes/

wget https://s3.denbi.uni-tuebingen.de/fb-test/tf_benchmark/cifar-10-binary.tar.gz -P datasets/tensorflow
tar xzfv datasets/tensorflow/cifar-10-binary.tar.gz -C /home/centos/datasets/tensorflow/

# Compile and install bowtie2 
cd bowtie2/bowtie2-2.3.4.2/
sudo make static-libs && make STATIC_BUILD=1
cd ../../

# Compile and install velvet
cd velvet/
make 'OPENMP=1'
cd ..

# Compile and install IDBA
cd IDBA/idba_ud-1.0.9/
./configure
make
cd ../../

# Compile and install tensorflow (sudo)
sudo pip install --upgrade --force-reinstall pip==9.0.3
sudo pip install tensorflow==1.4.0

# Download, compile and install GROMACS (sudo)
wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-2018.3.tar.gz -P gromacs/
gzip -dc /home/centos/gromacs/gromacs-2018.3.tar.gz | tar -xvf /home/centos/gromacs/
mkdir /home/centos/gromacs/gromacs-2018.3/build
cd gromacs/gromacs-2018.3/build/
cmake3 -DGMX_BUILD_OWN_FFTW=on -DGMX_GPU=off -DGMX_BUILD_MPI=off --build /home/centos/gromacs/gromacs-2018.3/build/  /home/centos/gromacs/gromacs-2018.3/
make
sudo make install

/usr/local/gromacs/bin/gmx grompp -f /home/centos/datasets/gromacs/adh_cubic/pme_verlet.mdp -c /home/centos/datasets/gromacs/adh_cubic/conf.gro -p /home/centos/datasets/gromacs/adh_cubic/topol.top -o /home/centos/datasets/gromacs/adh_cubic/topol -po /home/centos/datasets/gromacs/adh_cubic/mdout
