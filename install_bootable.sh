#!bin/bash

# Create empty directories
 ../timetest.txt
# Create empty directories
time_start=$(date)
if [ -d benchmark_output ]
then
        echo "Directory benchmark_output already exists."
else
	mkdir benchmark_output
fi

if [ -d benchmark_output/bowtie2 ]
then
        echo "Directory benchmark_output/bowtie2 already exists."
else
	mkdir benchmark_output/bowtie2
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

if [ -d bowtie2 ]
then
        echo "Directory bowtie2 already exists."
else
	mkdir bowtie2
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
time_end=$(date)
mkdir nmon_stats

# Download benchmark datasets

if [ -e datasets/1000_genomes/ERR251006.filt.fastq ]
then
	echo "File datasets/1000_genomes/ERR251006.filt.fastq already exists."
else
	#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00110/sequence_read/ERR251006.filt.fastq.gz -P datasets/1000_genomes/
	wget https://s3.denbi.uni-tuebingen.de/max/ERR251006.filt.fastq.gz -P datasets/1000_genomes/
	gunzip datasets/1000_genomes/ERR251006.filt.fastq.gz
fi

if [ -e datasets/1000_genomes/ERR016155.filt.fastq ]
then
	echo "File datasets/1000_genomes/ERR016155.filt.fastq already exists."
else
	#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00125/sequence_read/ERR016155.filt.fastq.gz -P datasets/1000_genomes/
	wget https://s3.denbi.uni-tuebingen.de/max/ERR016155.filt.fastq.gz -P datasets/1000_genomes/
	gunzip datasets/1000_genomes/ERR016155.filt.fastq.gz
fi

if [ -e datasets/1000_genomes/ERR015528.filt.fastq ]
then
        echo "File datasets/1000_genomes/ERR015528.filt.fastq already exists."
else
        #wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00106/sequence_read/ERR015528.filt.fastq.gz -P datasets/1000_genomes/
        wget https://s3.denbi.uni-tuebingen.de/max/ERR015528.filt.fastq.gz -P datasets/1000_genomes/
        gunzip datasets/1000_genomes/ERR015528.filt.fastq.gz
fi

if [ -e datasets/1000_genomes/SRR741411.filt.fastq ]
then
        echo "File datasets/1000_genomes/SRR741411.filt.fastq already exists."
else
        #wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00099/sequence_read/SRR741411.filt.fastq.gz -P datasets/1000_genomes/
        wget https://s3.denbi.uni-tuebingen.de/max/SRR741411.filt.fastq.gz -P datasets/1000_genomes/
        gunzip datasets/1000_genomes/SRR741411.filt.fastq.gz
fi


if [ -e datasets/ebi/DRR001012.fastq ]
then
        echo "File datasets/ebi/DRR001012.fastq already exists."
else
	#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR001/DRR001012/DRR001012.fastq.gz -P datasets/ebi/
	wget https://s3.denbi.uni-tuebingen.de/max/DRR001012.fastq.gz -P datasets/ebi/
	gunzip datasets/ebi/DRR001012.fastq.gz
fi

if [ -e datasets/ebi/DRR001025.fastq ]
then
        echo "File datasets/ebi/DRR001025.fastq already exists."
else
	#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR001/DRR001025/DRR001025.fastq.gz -P datasets/ebi/
	wget https://s3.denbi.uni-tuebingen.de/max/DRR001025.fastq.gz -P datasets/ebi/
	gunzip datasets/ebi/DRR001025.fastq.gz
fi

if [ -e datasets/1000_genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa ]
then
        echo "File datasets/1000_genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa already exists."
else
	#wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -P datasets/1000_genomes/
	wget https://s3.denbi.uni-tuebingen.de/max/GRCh38_full_analysis_set_plus_decoy_hla.fa -P datasets/1000_genomes/
fi

if [ -e datasets/clustalOmega/wgs.ANCA.1_200.fsa ]
then
        echo "File datasets/clustalOmega/wgs.ANCA.1_200.fsa already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/wgs.ANCA.1_200.fsa -P datasets/clustalOmega
fi


if [ -e datasets/clustalOmega/wgs.ANCA.1_400.fsa ]
then
        echo "File datasets/clustalOmega/wgs.ANCA.1_400.fsa already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/wgs.ANCA.1_400.fsa -P datasets/clustalOmega
fi


if [ -e datasets/clustalOmega/wgs.ANCA.1_500.fsa ]
then
        echo "File datasets/clustalOmega/wgs.ANCA.1_500.fsa already exists."
else
	wget https://s3.denbi.uni-tuebingen.de/max/wgs.ANCA.1_500.fsa -P datasets/clustalOmega
fi


if [ -d datasets/tensorflow/cifar-10-batches-bin ]
then
	echo "Directory datasets/tensorflow/cifar-10-batches-bin already exists."
else	
	wget https://s3.denbi.uni-tuebingen.de/max/cifar-10-binary.tar.gz -P datasets/tensorflow
	tar -xf datasets/tensorflow/cifar-10-binary.tar.gz -C datasets/tensorflow/
fi

if [ -d datasets/tensorflow/models ]
then
        echo "Directory datasets/tensorflow/models already exists."
else    
	wget https://s3.denbi.uni-tuebingen.de/max/models.tar.gz -P datasets/tensorflow
	tar -xf datasets/tensorflow/models.tar.gz -C datasets/tensorflow/
	rm datasets/tensorflow/models.tar.gz
fi

if [ -d datasets/gromacs/adh_cubic ]
then
        echo "Directory datasets/gromacs/adh_cubic already exists."
else
	wget https://s3.denbi.uni-tuebingen.de/max/ADH_bench_systems.tar.gz -P datasets/gromacs
	tar -xf datasets/gromacs/ADH_bench_systems.tar.gz -C datasets/gromacs/
	rm datasets/gromacs/ADH_bench_systems.tar.gz
fi


# Download GCC compiler version 7.3.0
if [ -d gcc/gcc-7.3.0 ]
then
	echo "Directory gcc/gcc-7.3.0 already exists."
else
	wget https://s3.denbi.uni-tuebingen.de/max/gcc-7.3.0.tar.gz -P gcc
	tar -xf gcc/gcc-7.3.0.tar.gz -C gcc/
	rm gcc/gcc-7.3.0.tar.gz
fi

# Download nmonchart34
if ! [ -d nmonchart ]
then
	wget https://s3.denbi.uni-tuebingen.de/max/nmonchart34.tar
	tar -xf nmonchart34.tar 
	rm nmonchart34.tar
else
	echo "Directory nmonchart already exists."
fi

# Download Bowtie2 Sources
if [ -d bowtie2/bowtie2-2.3.4.2 ]
then
        echo "Directory bowtie2/bowtie2-2.3.4.2 already exists."
else
	wget https://s3.denbi.uni-tuebingen.de/max/bowtie2-2.3.4.2-source.zip -P bowtie2
	unzip bowtie2/bowtie2-2.3.4.2-source.zip -d bowtie2/
	rm bowtie2/bowtie2-2.3.4.2-source.zip
fi

# Download Velvet github repository
if ! [ -z "$(ls -A "velvet")" ]
then
    echo "Directory velvet/ exists and is not empty."
else
	git clone https://github.com/dzerbino/velvet.git velvet/
	rm -rf velvet/.gitignore velvet/.git/
fi

# Download IDBA sources
if [ -d IDBA/idba_ud-1.0.9 ]
then
	echo "Directory IDBA/idba_ud-1.0.9 already exists."
else
	wget https://s3.denbi.uni-tuebingen.de/max/idba_ud-1.0.9.tar.gz -P IDBA
	tar -xf IDBA/idba_ud-1.0.9.tar.gz -C IDBA/
	rm IDBA/idba_ud-1.0.9.tar.gz
fi

# Download GROMACS sources
if [ -d gromacs/gromacs-2018.3 ]
then
	echo "Directory gromacs/gromacs-2018.3 already exists."
else
	#wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-2018.3.tar.gz -P gromacs/
	wget https://s3.denbi.uni-tuebingen.de/max/gromacs-2018.3.tar.gz -P gromacs/
	tar -xf gromacs/gromacs-2018.3.tar.gz -C gromacs/
	rm gromacs/gromacs-2018.3.tar.gz
fi

# Download SPAdes binaries
if [ -d SPAdes/SPAdes-3.12.0-Linux ]
then
        echo "Directory SPAdes/SPAdes-3.12.0-Linux already exists."
else
	wget https://s3.denbi.uni-tuebingen.de/max/SPAdes-3.12.0-Linux.tar.gz -P SPAdes
	tar -xf SPAdes/SPAdes-3.12.0-Linux.tar.gz -C SPAdes
	rm SPAdes/SPAdes-3.12.0-Linux.tar.gz
fi

# Download ClustalOmega sources
if [ -d clustalOmega/clustal-omega-1.2.4 ]
then
        echo "Directory clustalOmega/clustal-omega-1.2.4 already exists."
else
        wget https://s3.denbi.uni-tuebingen.de/max/clustal-omega-1.2.4.tar.gz -P clustalOmega
        tar -xf clustalOmega/clustal-omega-1.2.4.tar.gz -C clustalOmega/
        rm clustalOmega/clustal-omega-1.2.4.tar.gz
fi


# Save original PATH and LD_LIBRARY variables
original_path_variable=$(echo $PATH)
original_ld_library_variable=$(echo $LD_LIBRARY_PATH)

# Compile and install GCC 7.3.0
if [ -e gcc/gcc-installed/bin/gcc ]
then 
	while true; do
        	read -p "GCC 7.3.0 seems already to be installed, do you want to recompile it?" yn
                case $yn in
                        [Yy]* ) rm -rf gcc/gcc-build/*
				rm -rf gcc/gcc-installed/*
				cd gcc/gcc-7.3.0
                                ./contrib/download_prerequisites
                                cd ../gcc-build
                                ../gcc-7.3.0/configure --enable-languages=c,c++ --disable-multilib --prefix=$PWD/../gcc-installed
                                make clean
                                make -j$(nproc)
				make install
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
        ./contrib/download_prerequisites
        cd ../gcc-build
        ../gcc-7.3.0/configure --enable-languages=c,c++ --disable-multilib --prefix=$PWD/../gcc-installed
        make clean
        make -j$(nproc)
	make install
        cd ../../
fi


# Compile and install bowtie2
if [ -e bowtie2/bowtie2-2.3.4.2/bowtie2-align-l ]
then      
        while true; do
                read -p "Bowtie2 seems already to be installed, do you want to recompile it?" yn
                case $yn in
                        [Yy]* ) cd bowtie2/bowtie2-2.3.4.2/
        			make clean
        			sudo make static-libs
				make -j$(nproc) STATIC_BUILD=1
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
	make clean
	sudo make static-libs
	make -j$(nproc) STATIC_BUILD=1
	cd ../../
fi

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
        			make clean
        			make -j$(nproc) 'OPENMP=1'
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
	make clean
	make -j$(nproc) 'OPENMP=1'
	cd ..
fi

# Reset to system compiler as IDBA is not compiling with GCC 7.3.0
export PATH=$original_path_variable
export LD_LIBRARY_PATH=$original_ld_library_variable

# Compile and install IDBA
if [ -e IDBA/idba_ud-1.0.9/bin/idba_ud ]
then
        while true; do
                read -p "IDBA seems already to be installed, do you want to recompile it?" yn
                case $yn in
                        [Yy]* ) cd IDBA/idba_ud-1.0.9/
				./configure
				make clean
				make -j$(nproc)
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
	./configure
	make clean
	make -j$(nproc)
	cd ../../
fi

# Compile and install clustalOmega
if [ -e clustalOmega/clustal-omega-1.2.4/bin/clustalo ]
then
	while true; do
                read -p "clustalOmega seems already to be installed, do you want to recompile it?" yn
                case $yn in
                        [Yy]* ) cd clustalOmega/clustal-omega-1.2.4/
                                ./configure --prefix=$PWD
                                make clean
                                make -j$(nproc)
				make install
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
        ./configure --prefix=$PWD
        make clean
        make -j$(nproc)
	make install
        cd ../../
fi

# Compile and install tensorflow (sudo)

#check if pip is in general installed (if 0 then it is installed if 1 needs to be installed)
which pip > /dev/null 2>&1
pip_installed=$(echo $?)
if [ $pip_installed == 0 ]
then
	echo "pip is already installed, nothing to do here."
else
	while true; do
                read -p "Pip in general seems not to be installed or is not in your path, do you want to install it via yum?" yn
                case $yn in
                        [Yy]* ) sudo yum install pip -y
                                break;;
                        [Nn]* ) echo "Pip will not be installed, please install it on your own, or the Tensorflow benchmarks will not work"
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
                        [Yy]* ) sudo pip install --upgrade --force-reinstall pip==9.0.3
                                break;;
                        [Nn]* ) echo "pip will not be upgraded/downgraded to version 9.0.3, make sure tensorflow version 1.4.0 is working correctly"
                                break;;
                        * ) echo "Please answer yes or no."
                                ;;
                esac
        done
else
	echo "Pip is already installed in version 9.0.3, nothing to do here."
fi

#check if tensorflow is installed in version 1.4.0 (if 1 then it is installed if 0 needs to be installed)
pip_version_installed=$(pip list installed  2> /dev/null | grep -c "pip (9.0.3)")
tensorflow_installed=$(pip list installed  2> /dev/null | grep -c "tensorflow (1.4.0)")

if [ $tensorflow_installed == 0 ] && [ $pip_version_installed == 1 ]
then
	sudo pip install tensorflow==1.4.0
elif [ $tensorflow_installed == 1 ] && [ $pip_version_installed == 1 ]
then 
	while true; do
                read -p "Tensorflow seems already to be installed in version 1.4.0 via pip, do you want to reinstall it?" yn
                case $yn in
                        [Yy]* ) sudo pip install --force-reinstall tensorflow==1.4.0
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

# Set GCC to 7.3.0 to fasten up GROMACS
export PATH=$PWD/gcc/gcc-installed/bin:$PATH
export LD_LIBRARY_PATH=$PWD/gcc/gcc-installed/lib64:$LD_LIBRARY_PATH

# Compile and install GROMACS (sudo)
if [ -e /usr/local/gromacs/bin/gmx ]
then
        while true; do
                read -p "GROMACS seems already to be installed, do you want to recompile it?" yn
                case $yn in
                        [Yy]* ) rm -rf gromacs/gromacs-2018.3/build/*
				sudo rm -rf /usr/local/gromacs/ 
				cd gromacs/gromacs-2018.3/build/
				CC=gcc CXX=g++ cmake3 -DCMAKE_CXX_COMPILER=g++ -DGMX_BUILD_OWN_FFTW=on -DGMX_GPU=off -DGMX_BUILD_MPI=off --build ./ ../../gromacs-2018.3/
				make clean
				make -j$(nproc)
				sudo make install
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
	CC=gcc CXX=g++ cmake3 -DCMAKE_CXX_COMPILER=g++ -DGMX_BUILD_OWN_FFTW=on -DGMX_GPU=off -DGMX_BUILD_MPI=off --build ./ ../../gromacs-2018.3/
	make clean
	make -j	$(nproc)
	sudo make install
	cd ../../../
	#/usr/local/gromacs/bin/gmx grompp -f datasets/gromacs/adh_cubic/pme_verlet.mdp -c datasets/gromacs/adh_cubic/conf.gro -p datasets/gromacs/adh_cubic/topol.top -o datasets/gromacs/adh_cubic/topol -po datasets/gromacs/adh_cubic/mdout
fi

# Convert fastq files to fasta files
echo "Converting datasets from .fastq to .fa"

if [ -e datasets/1000_genomes/ERR016155.filt.fa ]
then
	echo "File datasets/1000_genomes/ERR016155.filt.fa already exists."
else
	IDBA/idba_ud-1.0.9/bin/fq2fa datasets/1000_genomes/ERR016155.filt.fastq datasets/1000_genomes/ERR016155.filt.fa
fi

if [ -e datasets/1000_genomes/ERR251006.filt.fa ]
then
        echo "File datasets/1000_genomes/ERR251006.filt.fa already exists."
else
	IDBA/idba_ud-1.0.9/bin/fq2fa datasets/1000_genomes/ERR251006.filt.fastq datasets/1000_genomes/ERR251006.filt.fa
fi

if [ -e datasets/1000_genomes/ERR015528.filt.fa ]
then
        echo "File datasets/1000_genomes/ERR015528.filt.fa already exists."
else
        IDBA/idba_ud-1.0.9/bin/fq2fa datasets/1000_genomes/ERR015528.filt.fastq datasets/1000_genomes/ERR015528.filt.fa
fi

if [ -e datasets/1000_genomes/SRR741411.filt.fa ]
then
        echo "File datasets/1000_genomes/SRR741411.filt.fa already exists."
else
        IDBA/idba_ud-1.0.9/bin/fq2fa datasets/1000_genomes/SRR741411.filt.fastq datasets/1000_genomes/SRR741411.filt.fa
fi

if [ -e datasets/ebi/DRR001012.fa ]
then
        echo "File datasets/ebi/DRR001012.fa already exists."
else
	IDBA/idba_ud-1.0.9/bin/fq2fa datasets/ebi/DRR001012.fastq datasets/ebi/DRR001012.fa
fi

if [ -e datasets/ebi/DRR001025.fa ]
then
        echo "File datasets/ebi/DRR001025.fa already exists."
else
	IDBA/idba_ud-1.0.9/bin/fq2fa datasets/ebi/DRR001025.fastq datasets/ebi/DRR001025.fa
fi

# Reset to system compiler after finishing installation process
export PATH=$original_path_variable
export LD_LIBRARY_PATH=$original_ld_library_variable

