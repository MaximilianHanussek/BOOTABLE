FROM centos:latest
RUN yum update -y
RUN yum install epel-release -y
RUN yum group install "Development Tools" -y
RUN yum install nano curl wget vim htop time git zlib-devel.x86_64 python-pip cmake3 tbb-devel.x86_64 R -y
RUN mkdir ~/benchmark_output
RUN mkdir ~/benchmark_output/bowtie2
RUN mkdir ~/benchmark_output/gromacs
RUN mkdir ~/benchmark_output/IDBA
RUN mkdir ~/benchmark_output/SPAdes
RUN mkdir ~/benchmark_output/tensorflow
RUN mkdir ~/benchmark_output/velvet
RUN mkdir ~/datasets
RUN mkdir ~/datasets/1000_genomes
RUN mkdir ~/datasets/gromacs
RUN mkdir ~/datasets/tensorflow
RUN mkdir ~/gromacs
RUN mkdir ~/results
RUN mkdir ~/bowtie2
RUN mkdir ~/IDBA
RUN mkdir ~/SPAdes
RUN mkdir ~/tensorflow
RUN mkdir ~/velvet
#RUN cd ~/datasets/gromacs/
RUN wget ftp://ftp.gromacs.org/pub/benchmarks/ADH_bench_systems.tar.gz -P ~/datasets/gromacs
RUN tar -xf ~/datasets/gromacs/ADH_bench_systems.tar.gz -C ~/datasets/gromacs/
RUN rm ~/datasets/gromacs/ADH_bench_systems.tar.gz
#RUN cd ~/datasets/tensorflow/
RUN wget https://s3.denbi.uni-tuebingen.de/fb-test/tf_benchmark/models.tar.gz -P ~/datasets/tensorflow
RUN tar -xf ~/datasets/tensorflow/models.tar.gz -C ~/datasets/tensorflow/
RUN rm ~/datasets/tensorflow/models.tar.gz
#RUN cd ~/bowtie2/
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.2/bowtie2-2.3.4.2-source.zip -P ~/bowtie2

RUN unzip ~/bowtie2/bowtie2-2.3.4.2-source.zip -d ~/bowtie2/
RUN rm ~/bowtie2/bowtie2-2.3.4.2-source.zip
#RUN cd ../
RUN git clone https://github.com/dzerbino/velvet.git ~/velvet/
#RUN cd velvet/
RUN rm -rf ~/velvet/.gitignore ~/velvet/.git/
#RUN cd ~/IDBA/
RUN wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/hku-idba/idba_ud-1.0.9.tar.gz -P ~/IDBA
RUN tar -xf ~/IDBA/idba_ud-1.0.9.tar.gz -C ~/IDBA/
RUN rm ~/IDBA/idba_ud-1.0.9.tar.gz
#RUN cd ~/SPAdes/
RUN wget http://cab.spbu.ru/files/release3.12.0/SPAdes-3.12.0-Linux.tar.gz -P ~/SPAdes
RUN tar -xf ~/SPAdes/SPAdes-3.12.0-Linux.tar.gz -C ~/SPAdes
RUN rm ~/SPAdes/SPAdes-3.12.0-Linux.tar.gz
#RUN cd ~/
RUN wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00110/sequence_read/ERR251006.filt.fastq.gz -P ~/datasets/1000_genomes/
RUN gunzip ~/datasets/1000_genomes/ERR251006.filt.fastq.gz
RUN wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00125/sequence_read/ERR016155.filt.fastq.gz -P ~/datasets/1000_genomes/
RUN gunzip ~/datasets/1000_genomes/ERR016155.filt.fastq.gz
#RUN wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -P ~/datasets/1000_genomes/
RUN wget https://s3.denbi.uni-tuebingen.de/max/GRCh38_full_analysis_set_plus_decoy_hla.fa -P ~/datasets/1000_genomes/
RUN wget https://s3.denbi.uni-tuebingen.de/fb-test/tf_benchmark/cifar-10-binary.tar.gz -P ~/datasets/tensorflow
RUN tar xf ~/datasets/tensorflow/cifar-10-binary.tar.gz -C ~/datasets/tensorflow/
#RUN cd ~/bowtie2/bowtie2-2.3.4.2/
RUN make -C ~/bowtie2/bowtie2-2.3.4.2/ static-libs  && make -C ~/bowtie2/bowtie2-2.3.4.2/ STATIC_BUILD=1
#RUN cd ~/velvet/
RUN make -C ~/velvet/ 'OPENMP=1'
RUN cd ~/IDBA/idba_ud-1.0.9/ && ~/IDBA/idba_ud-1.0.9/configure
RUN make -C ~/IDBA/idba_ud-1.0.9/
#RUN cd ~/
RUN pip install --upgrade --force-reinstall pip==9.0.3
RUN pip install tensorflow==1.4.0
RUN wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-2018.3.tar.gz -P ~/gromacs/
RUN tar -xf ~/gromacs/gromacs-2018.3.tar.gz -C ~/gromacs/
RUN rm ~/gromacs/gromacs-2018.3.tar.gz
RUN mkdir ~/gromacs/gromacs-2018.3/build
RUN cd ~/gromacs/gromacs-2018.3/build/ && cmake3 -DGMX_BUILD_OWN_FFTW=on -DGMX_GPU=off -DGMX_BUILD_MPI=off --build ./  ~/gromacs/gromacs-2018.3/
RUN make -C ~/gromacs/gromacs-2018.3/build/ 
RUN make -C ~/gromacs/gromacs-2018.3/build/ install
#RUN cd ~/
RUN /usr/local/gromacs/bin/gmx grompp -f ~/datasets/gromacs/adh_cubic/pme_verlet.mdp -c ~/datasets/gromacs/adh_cubic/conf.gro -p ~/datasets/gromacs/adh_cubic/topol.top -o ~/datasets/gromacs/adh_cubic/topol -po ~/datasets/gromacs/adh_cubic/mdout
RUN ~/IDBA/idba_ud-1.0.9/bin/fq2fa ~/datasets/1000_genomes/ERR016155.filt.fastq ~/datasets/1000_genomes/ERR016155.filt.fa
RUN ~/IDBA/idba_ud-1.0.9/bin/fq2fa ~/datasets/1000_genomes/ERR251006.filt.fastq ~/datasets/1000_genomes/ERR251006.filt.fa
COPY run_benchmarks.sh /root/
COPY threadedBioBenchsuiteStatsGenerator.R /root/
