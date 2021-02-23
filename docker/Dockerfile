FROM centos:7.7.1908
RUN yum update -y && yum install epel-release -y && yum group install "Development Tools" -y && yum install -y \
	nano \
	curl \
     	wget \ 
	vim \
	htop \ 
	time \
	git \
	nmon\
	zlib-devel.x86_64 \
	python-pip \
	inxi \ 
	cmake3 \
	tbb-devel.x86_64 \
	argtable-devel \
	util-linux \ 
	hwloc \ 
	hwloc-devel \ 
	gmp-devel \ 
	mpfr-devel \ 
	libmpc-devel \ 
	R

RUN mkdir ~/benchmark_output \
	~/benchmark_output/bowtie2 \
	~/benchmark_output/BBMap \
	~/benchmark_output/BWA \
	~/benchmark_output/clustalOmega \
	~/benchmark_output/MAFFT \
	~/benchmark_output/SINA \ 
	~/benchmark_output/gromacs \
	~/benchmark_output/IDBA \
	~/benchmark_output/SPAdes \
	~/benchmark_output/tensorflow \
	~/benchmark_output/velvet \
	~/datasets \
	~/datasets/1000_genomes \
	~/datasets/clustalOmega \
	~/datasets/gromacs \
	~/datasets/tensorflow \
	~/datasets/ebi \
	~/datasets/BWA \
	~/datasets/BWA/DRR001012 \
	~/datasets/BWA/DRR001025 \
	~/datasets/BWA//GRCh38_full_analysis_set_plus_decoy_hla \
	~/datasets/SINA \
	~/backed_up_benchmark_results \
	~/gromacs \
	~/results \
	~/BBMap \
	~/BWA \
	~/bowtie2 \
	~/IDBA \
	~/SPAdes \
	~/tensorflow \
	~/velvet \
	~/clustalOmega \
	~/MAFFT \
	~/SINA \
	~/gcc \
	~/gcc/gcc-build \
	~/gcc/gcc-installed \
	~/nmon_stats

RUN original_path_variable=$(echo $PATH) \
&& original_ld_library_variable=$(echo $LD_LIBRARY_PATH)

RUN wget https://s3.denbi.uni-tuebingen.de/max/BBMap_38.87.tar.gz -P ~/BBMap \
&& tar -xf ~/BBMap/BBMap_38.87.tar.gz -C ~/BBMap/ \
&& rm ~/BBMap/BBMap_38.87.tar.gz

RUN wget https://s3.denbi.uni-tuebingen.de/max/bowtie2-2.4.2-source.zip -P ~/bowtie2 \
&& unzip ~/bowtie2/bowtie2-2.4.2-source.zip -d ~/bowtie2/ \
&& rm ~/bowtie2/bowtie2-2.4.2-source.zip

RUN wget https://s3.denbi.uni-tuebingen.de/max/bwa-0.7.17.tar.bz2 -P ~/BWA \
&& tar -xf ~/BWA/bwa-0.7.17.tar.bz2 -C ~/BWA/ \
&& rm ~/BWA/bwa-0.7.17.tar.bz2

RUN git clone https://github.com/dzerbino/velvet.git ~/velvet/ \
&& rm -rf ~/velvet/.gitignore ~/velvet/.git/

RUN wget https://s3.denbi.uni-tuebingen.de/max/idba_ud-1.0.9.tar.gz -P ~/IDBA \
&& tar -xf ~/IDBA/idba_ud-1.0.9.tar.gz -C ~/IDBA/ \
&& rm ~/IDBA/idba_ud-1.0.9.tar.gz

RUN wget https://s3.denbi.uni-tuebingen.de/max/gromacs-2018.3.tar.gz -P ~/gromacs/ \
&& tar -xf ~/gromacs/gromacs-2018.3.tar.gz -C ~/gromacs/ \
&& rm ~/gromacs/gromacs-2018.3.tar.gz

RUN wget https://s3.denbi.uni-tuebingen.de/max/SPAdes-3.12.0-Linux.tar.gz -P ~/SPAdes \
&& tar -xf ~/SPAdes/SPAdes-3.12.0-Linux.tar.gz -C ~/SPAdes \
&& rm ~/SPAdes/SPAdes-3.12.0-Linux.tar.gz

RUN wget https://s3.denbi.uni-tuebingen.de/max/clustal-omega-1.2.4.tar.gz -P ~/clustalOmega \
&& tar -xf ~/clustalOmega/clustal-omega-1.2.4.tar.gz -C ~/clustalOmega/ \
&& rm ~/clustalOmega/clustal-omega-1.2.4.tar.gz

RUN wget https://s3.denbi.uni-tuebingen.de/max/mafft-7.475-with-extensions-src.tgz -P ~/MAFFT \
&& tar -xf ~/MAFFT/mafft-7.475-with-extensions-src.tgz -C ~/MAFFT/ \
&& rm ~/MAFFT/mafft-7.475-with-extensions-src.tgz

RUN wget https://s3.denbi.uni-tuebingen.de/max/sina-1.7.2-linux.tar.gz -P ~/SINA \
&& tar -xf ~/SINA/sina-1.7.2-linux.tar.gz -C ~/SINA/ \
&& rm ~/SINA/sina-1.7.2-linux.tar.gz

RUN wget https://s3.denbi.uni-tuebingen.de/max/gcc-7.3.0.tar.gz -P ~/gcc \
&& tar -xf ~/gcc/gcc-7.3.0.tar.gz -C ~/gcc/ \
&& rm ~/gcc/gcc-7.3.0.tar.gz

RUN wget https://s3.denbi.uni-tuebingen.de/max/ERR251006.filt.fastq.gz -P ~/datasets/1000_genomes/ \
&& gunzip ~/datasets/1000_genomes/ERR251006.filt.fastq.gz

RUN wget https://s3.denbi.uni-tuebingen.de/max/ERR016155.filt.fastq.gz -P ~/datasets/1000_genomes/ \
&& gunzip ~/datasets/1000_genomes/ERR016155.filt.fastq.gz

RUN wget https://s3.denbi.uni-tuebingen.de/max/ERR015528.filt.fastq.gz -P ~/datasets/1000_genomes/ \
&& gunzip ~/datasets/1000_genomes/ERR015528.filt.fastq.gz

RUN wget https://s3.denbi.uni-tuebingen.de/max/SRR741411.filt.fastq.gz -P ~/datasets/1000_genomes/ \
&& gunzip ~/datasets/1000_genomes/SRR741411.filt.fastq.gz


RUN wget https://s3.denbi.uni-tuebingen.de/max/GRCh38_full_analysis_set_plus_decoy_hla.fa -P ~/datasets/1000_genomes/

RUN wget https://s3.denbi.uni-tuebingen.de/max/GRCh38_full_analysis_set_plus_decoy_hla.amb -P ~/datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla

RUN wget https://s3.denbi.uni-tuebingen.de/max/GRCh38_full_analysis_set_plus_decoy_hla.ann -P ~/datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla

RUN wget https://s3.denbi.uni-tuebingen.de/max/GRCh38_full_analysis_set_plus_decoy_hla.bwt -P ~/datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla

RUN wget https://s3.denbi.uni-tuebingen.de/max/GRCh38_full_analysis_set_plus_decoy_hla.pac -P ~/datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla

RUN wget https://s3.denbi.uni-tuebingen.de/max/GRCh38_full_analysis_set_plus_decoy_hla.sa -P ~/datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla

RUN wget https://s3.denbi.uni-tuebingen.de/max/wgs.ANCA.1_200.fsa -P ~/datasets/clustalOmega

RUN wget https://s3.denbi.uni-tuebingen.de/max/wgs.ANCA.1_400.fsa -P ~/datasets/clustalOmega

RUN wget https://s3.denbi.uni-tuebingen.de/max/wgs.ANCA.1_500.fsa -P ~/datasets/clustalOmega

RUN  wget https://s3.denbi.uni-tuebingen.de/max/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb.gz -P ~/datasets/SINA \
&& gunzip ~/datasets/SINA/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb.gz

RUN  wget https://s3.denbi.uni-tuebingen.de/max/SILVA_138.1_SSURef_NR99_12_06_20_opt.sidx -P ~/datasets/SINA

RUN  wget https://s3.denbi.uni-tuebingen.de/max/GTDB_bac-arc_ssu_r86.fa.gz -P ~/datasets/SINA \
&& gunzip ~/datasets/SINA/GTDB_bac-arc_ssu_r86.fa.gz

RUN  wget https://s3.denbi.uni-tuebingen.de/max/RefSeq-RDP16S_v2_May2018.fa.gz -P ~/datasets/SINA \
&& gunzip ~/datasets/SINA/RefSeq-RDP16S_v2_May2018.fa.gz

RUN  wget https://s3.denbi.uni-tuebingen.de/max/OE-38_R1.fastq.gz -P ~/datasets/SINA \
&& gunzip ~/datasets/SINA/OE-38_R1.fastq.gz

RUN wget https://s3.denbi.uni-tuebingen.de/max/cifar-10-binary.tar.gz -P ~/datasets/tensorflow \
&& tar xf ~/datasets/tensorflow/cifar-10-binary.tar.gz -C ~/datasets/tensorflow/

RUN wget https://s3.denbi.uni-tuebingen.de/max/models.tar.gz -P ~/datasets/tensorflow \
&& tar -xf ~/datasets/tensorflow/models.tar.gz -C ~/datasets/tensorflow/ \
&& rm ~/datasets/tensorflow/models.tar.gz

RUN wget https://s3.denbi.uni-tuebingen.de/max/DRR001012.fastq.gz -P ~/datasets/ebi/ \
&& gunzip ~/datasets/ebi/DRR001012.fastq.gz

RUN wget https://s3.denbi.uni-tuebingen.de/max/DRR001012.amb -P ~/datasets/BWA/DRR001012

RUN wget https://s3.denbi.uni-tuebingen.de/max/DRR001012.ann -P ~/datasets/BWA/DRR001012

RUN wget https://s3.denbi.uni-tuebingen.de/max/DRR001012.bwt -P ~/datasets/BWA/DRR001012

RUN wget https://s3.denbi.uni-tuebingen.de/max/DRR001012.pac -P ~/datasets/BWA/DRR001012

RUN wget https://s3.denbi.uni-tuebingen.de/max/DRR001012.sa -P ~/datasets/BWA/DRR001012

RUN wget https://s3.denbi.uni-tuebingen.de/max/DRR001025.fastq.gz -P ~/datasets/ebi/ \
&& gunzip ~/datasets/ebi/DRR001025.fastq.gz

RUN wget https://s3.denbi.uni-tuebingen.de/max/DRR001025.amb -P ~/datasets/BWA/DRR001025

RUN wget https://s3.denbi.uni-tuebingen.de/max/DRR001025.ann -P ~/datasets/BWA/DRR001025

RUN wget https://s3.denbi.uni-tuebingen.de/max/DRR001025.bwt -P ~/datasets/BWA/DRR001025

RUN wget https://s3.denbi.uni-tuebingen.de/max/DRR001025.pac -P ~/datasets/BWA/DRR001025

RUN wget https://s3.denbi.uni-tuebingen.de/max/DRR001025.sa -P ~/datasets/BWA/DRR001025

RUN wget https://s3.denbi.uni-tuebingen.de/max/ADH_bench_systems.tar.gz -P ~/datasets/gromacs \
&& tar -xf ~/datasets/gromacs/ADH_bench_systems.tar.gz -C ~/datasets/gromacs/ \
&& rm ~/datasets/gromacs/ADH_bench_systems.tar.gz

RUN wget https://s3.denbi.uni-tuebingen.de/max/nmonchart34.tar -P ~/ \
&& tar -xf ~/nmonchart34.tar -C ~/ \
&& rm ~/nmonchart34.tar

RUN cd ~/gcc/gcc-7.3.0 \
&& ./contrib/download_prerequisites

RUN cd ~/gcc/gcc-build \
&& ../gcc-7.3.0/configure --enable-languages=c,c++ --disable-multilib --prefix=$PWD/../gcc-installed \
&& make clean \
&& make -j$(nproc) \
&& make install

RUN cd ~/BWA/bwa-0.7.17/ \
&& make clean \
&& make -j$(nproc) 

RUN make -C ~/bowtie2/bowtie2-2.4.2/ static-libs \
&& make -C ~/bowtie2/bowtie2-2.4.2/ STATIC_BUILD=1

RUN export PATH=$PWD/gcc/gcc-installed/bin:$PATH \
&& export LD_LIBRARY_PATH=$PWD/gcc/gcc-installed/lib64:$LD_LIBRARY_PATH

RUN make -C ~/velvet/ 'OPENMP=1'

RUN export PATH=$original_path_variable \
&& export LD_LIBRARY_PATH=$original_ld_library_variable

RUN cd ~/IDBA/idba_ud-1.0.9/ \
&& ~/IDBA/idba_ud-1.0.9/configure \
&& make -C ~/IDBA/idba_ud-1.0.9/

RUN cd ~/clustalOmega/clustal-omega-1.2.4/ \
&& ~/clustalOmega/clustal-omega-1.2.4/configure --prefix=$PWD \
&& make -j$(nproc) -C ~/clustalOmega/clustal-omega-1.2.4/ \
&& make -C ~/clustalOmega/clustal-omega-1.2.4/ install

RUN sed -i "/PREFIX = /c\PREFIX = /root/MAFFT/mafft-7.475-with-extensions" ~/MAFFT/mafft-7.475-with-extensions/core/Makefile \
&& cd ~/MAFFT/mafft-7.475-with-extensions/core/ \
&& make clean \
&& make -j$(nproc) \
&& make install

RUN pip install --upgrade --force-reinstall pip==9.0.3
RUN pip install tensorflow==1.4.0

RUN export PATH=$PWD/gcc/gcc-installed/bin:$PATH \
&& export LD_LIBRARY_PATH=$PWD/gcc/gcc-installed/lib64:$LD_LIBRARY_PATH 

RUN mkdir ~/gromacs/gromacs-2018.3/build

RUN cd ~/gromacs/gromacs-2018.3/build/ \
&& CC=gcc CXX=g++ cmake3 -DCMAKE_CXX_COMPILER=g++ -DGMX_BUILD_OWN_FFTW=on -DGMX_GPU=off -DGMX_BUILD_MPI=off --build ./  ~/gromacs/gromacs-2018.3/

RUN make -C ~/gromacs/gromacs-2018.3/build/ \
&& make -C ~/gromacs/gromacs-2018.3/build/ install

RUN export PATH=$original_path_variable \
&& export LD_LIBRARY_PATH=$original_ld_library_variable


#RUN /usr/local/gromacs/bin/gmx grompp -f ~/datasets/gromacs/adh_cubic/pme_verlet.mdp -c ~/datasets/gromacs/adh_cubic/conf.gro -p ~/datasets/gromacs/adh_cubic/topol.top -o ~/datasets/gromacs/adh_cubic/topol -po ~/datasets/gromacs/adh_cubic/mdout
RUN ~/IDBA/idba_ud-1.0.9/bin/fq2fa ~/datasets/1000_genomes/ERR016155.filt.fastq ~/datasets/1000_genomes/ERR016155.filt.fa
RUN ~/IDBA/idba_ud-1.0.9/bin/fq2fa ~/datasets/1000_genomes/ERR251006.filt.fastq ~/datasets/1000_genomes/ERR251006.filt.fa
RUN ~/IDBA/idba_ud-1.0.9/bin/fq2fa ~/datasets/1000_genomes/ERR015528.filt.fastq ~/datasets/1000_genomes/ERR015528.filt.fa
RUN ~/IDBA/idba_ud-1.0.9/bin/fq2fa ~/datasets/1000_genomes/SRR741411.filt.fastq ~/datasets/1000_genomes/SRR741411.filt.fa
RUN ~/IDBA/idba_ud-1.0.9/bin/fq2fa ~/datasets/ebi/DRR001012.fastq ~/datasets/ebi/DRR001012.fa
RUN ~/IDBA/idba_ud-1.0.9/bin/fq2fa ~/datasets/ebi/DRR001025.fastq ~/datasets/ebi/DRR001025.fa
RUN ~/IDBA/idba_ud-1.0.9/bin/fq2fa ~/datasets/SINA/OE-38_R1.fastq ~/datasets/SINA/OE-38_R1.fa


COPY run_benchmarks.sh /root/
COPY BOOTABLE_report_generator.R /root/
COPY install_check.sh /root/
COPY BOOTABLE_collect.sh /root/
