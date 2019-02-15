# BOOTABLE                                                                                    [![DOI](https://zenodo.org/badge/149121128.svg)](https://zenodo.org/badge/latestdoi/149121128) [![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2274)




BiOinfOrmatics ThreAded Benchmark tooLsuitE (BOOTABLE)
This toolsuite currently consists of the following tools
- Bowtie2
- Velvet
- IDBA
- GROMACS
- Tensowflow
- SPAdes
- Clustal Omega (currently in testing phase)

more are coming soon ...
 
All the tools and datasets included in this bundle are carfeully chosen to 
cover a range of different application cases and make strong use of multithreading 
in order to get some benchmarks of the underlying hardware and about the scalability of the system.
In the following is explained how to install BOOTABLE, how to run it and what kind of options you have. Information about results are stated in the **Results**  section.

## Prerequisites
This current version is only tested for CentOS 7. Ubuntu support will be added soon.  
If you want to install the tools and scripts BOOTABLE use please be sure to have following packages installed via yum (CentOS):
- epel-release
- Development Tools (group install)
- nano
- curl
- wget
- vim
- htop
- time
- git
- nmon
- zlib-devel.x86_64
- python-pip (in Version 9.0.3)
- cmake3
- tbb-devel.x86_64
- argtable-devel
- R
- inxi
- uitl-linux
- hwloc
- hwloc-devel
- gmp-devel 
- mpfr-devel 
- libmpc-devel

<pre>yum update -y
yum install epel-release -y
yum group install "Development Tools" -y
yum install nano curl wget vim htop nmon time git zlib-devel.x86_64 python-pip inxi cmake3 tbb-devel.x86_64 argtable-devel util-linux hwloc hwloc-devel gmp-devel mpfr-devel libmpc-devel R -y</pre>


Further the system needs access to the internet to load the specific datasets and some R packages.

You can also install the required R packages in beforehand which are:
- grid
- gridExtra
- RColorBrewer
- stringr

## Short guide
If you can not wait to test your system or tool and want to read through the whole tutorial here just the commands to start on a bare metal system.

- Clone the github repo
<pre>git clone https://github.com/MaximilianHanussek/BOOTABLE.git</pre>

- Change into the root directory of the cloned repo
<pre>cd BOOTABLE/</pre>

- Install BOOTABLE
<pre>sh install_bootable.sh</pre>

- Check if installation is correct
<pre>sh install_check.sh</pre>

- Run your first benchmark
<pre>sh run_benchmarks.sh</pre>

- After the benchmarking has finished create a report if wanted
<pre>Rscript --vanilla threadedBioBenchsuiteStatsGenerator.R</pre>

## Installation

There are 6 different kinds of installation possibilities.
1. Bare metal installation
2. Using a .qcow2 image with everything already preinstalled
3. Using Docker and the provided Dockerfile to install everything into a Docker Container
4. Using Singularity and the provided Singularity file
5. Transform the Docker container from option 3. into a Singualrity Container
6. Use Ansible and the provided Ansible script to install the tools on a CentOS based system from scratch

### 1. Bare metal installation
In order to compile and install the tools directly on a server without a virtualization technique 
in between, clone this github repo and run the `install_bootable.sh` script like in the following.

Clone the github repo 
<pre>git clone https://github.com/MaximilianHanussek/BOOTABLE.git</pre>

Change into the root directory of the cloned repo
<pre>cd BOOTABLE/</pre>

Run the installation script and all the tools will be compiled and installed in the BOOTABLE directory except some GROMACS tools which will be installed under /usr/local/ at the moment and tensorflow which will be installed via `pip`
<pre>sh install_bootable.sh</pre>

Afterwards please run the `install_ckeck.sh` tool from the BOOTABLE root directory in order to check everything is installed
and working as it should be. 
<pre>sh install_check.sh</pre>

If everything is correct all the output on the screen will appear in green. If not, the problematic tools/files will be marked in red.

If everything is green you can start with running your first benchmark.

### 2. Using a .qcow2 image
This instructions tell you how to use a .qcow2 image with everything preinstalled. It is further assumed that you know how to install an image in an virtual environment or have a virtual environment already running. Make sure you have at least 50GB of disk space. 
The image is deposited in an S3 bucket with the following URL `https://s3.denbi.uni-tuebingen.de/max/benchmark_image.qcow2`.
You can just download it into your current directory via `wget` for example:

<pre>wget https://s3.denbi.uni-tuebingen.de/max/benchmark_image.qcow2</pre>

Just start the image in any virtuel environment you like, for example with `virt-manager` under CentOS.
You can login with the following credentials:
- Username: root
- Password: rootdenbi

After the login please change the password of the user `centos` via the passwd command to any you like.
Afterwards logout and login to the virtual machine as user `centos` with the password you have entered before.
Change into the home directory (`/home/centos/`), where you will find all tools already installed. Start the benchmarks by running the following command:

<pre>sh run_benchmarks</pre>

This will start the calculations of the different tools first with the maximal number of CPU cores, second the half of the maximal number of CPU cores and at the end with one CPU core. Each number of CPU cores is running three times to get the mean value. As a comparison value, on a 28 core machine this will take some days (~80 hours).

After the benchmarks have finished run the following Rscript to generate a brief report in .pdf format for each number of benchmarked CPUs:

<pre>Rscript --vanilla report_generator.R</pre>

If you want to re-run the benchmarks use the -c flag in order to delete the created output files and benchmark results from the run before:

<pre>sh run_benchmarks -c</pre>


### 3. Using a Docker container
This instructions assume that docker is already installed and running. If not you should find most of the information
on the [Docker website](https://www.docker.com/get-started)

The **first** possibility is to pull the already configured BOOTABLE Docker image from [docker hub](https://hub.docker.com/)
with the following docker pull command.

<pre>docker pull maximilianhanussek/bootable</pre>

After the image has been downloaded you can start an interactive Docker container with the following command
<pre>docker run --rm -it maximilianhanussek/bootable</pre>

After you are in the container change into the home directory (root)
<pre>cd /root/</pre>

In the home directory you will find all the tools already installed. Before you start the benchmark you can run the install_check.sh script to be sure everything is fine 
<pre>sh install_check.sh</pre>

In order to start the benchmark just run
<pre>sh run_benchmarks.sh</pre>

or enter some additional parameters via the provided flags.

The **second** option is to build the Docker container by yourself with the provided `Dockerfile`. You can find it in the github repo in the `docker` directory.

To build the Docker container by yourself clone this github repo, change into the `docker` directory and run the following command:

<pre>docker build --tag bootable .</pre>

The build will take a couple of minutes (30-60) as the first step is to update the underlying operating sytem /(centOS7) to the latest version and installing all the required packages. Further some large datasets need to be downloaded and this depends on your network conectivity. Afterwards most of the tools have to be compiled and installed.

After the Docker image has been build you can start the benchmark with the same command like if you had pulled it from Docker Hub:

<pre>docker run --rm -it maximilianhanussek/bootable</pre>
<pre>cd /root/</pre>
<pre>sh run_benchmarks.sh</pre>

### 4. Using a Singularity container
This instructions assume that singularity is already installed and working. If not you should find most of the information
on the [Singularity website](https://www.sylabs.io/guides/2.6/user-guide/installation.html)

The **first** possibility is to pull the already configured and build container from [singularity hub]()
with the following singularity pull command.

<pre>singularity pull shub://MaximilianHanussek/BOOTABLE</pre>

After the image has been downloaded you can start an interactive Singularity container with the following command
<pre>sudo singularity shell bootable.simg</pre>

After you are in the container change to the root directory
<pre>cd /</pre>

In this directory you will find all the tools already installed. Before you start the benchmark you can run the install_check.sh script to be sure everything is fine 
<pre>sh install_check.sh</pre>

In order to start the benchmark just run
<pre>sh run_benchmarks.sh</pre>

or enter some additional parameters via the provided flags.

The **second** option is to pull the already configured and build BOOTABLE Docker container from [docker hub](https://hub.docker.com/) with the following singularity pull command.

<pre>singularity pull docker://maximilianhanussek/bootable</pre>

After the image has been downloaded into your current working directory and converted into a singularity container you can start an interactive Singularity container with the following command
<pre>sudo singularity shell MaximilianHanussek-BOOTABLE-master-latest.simg</pre>

After you are in the container, change to the root home directory /root/
<pre>cd /root/</pre>

In this directory you will find all the tools already installed. In order to start the benchmark just run
<pre>sh run_benchmarks.sh</pre>

or enter some additional parameters via the provided flags.

The **third** option is to build the Singularity container by yourself with the provided `Singularity recipe`. You can find it in the github repo in the `singularity` directory.

To build the Docker container by yourself clone this github repo, change into the `singularity` directory and run the following command:

<pre>sudo singularity build bootable.simg singularity_bootable_recipe</pre>

The build will take a couple of minutes (30-60) as the first step is to update the underlying operating sytem /(centOS7) to the latest version and installing all the required packages. Further some large datasets need to be downloaded and this depends on your network conectivity. Afterwards most of the tools have to be compiled and installed.

After the Singularity image has been build you can start the benchmark with the same command like if you had pulled it from Singularity Hub:

<pre>sudo singularity shell bootable.simg</pre>
<pre>cd /root/</pre>
<pre>sh run_benchmarks.sh</pre>

### 5. Using Ansible

## Start Benchmarks
In order to start a benchmark just stay in the BOOTABLE root directory and run the `run_benchmarks.sh` script.
If you provide no flags the default parameters will be chosen like stated in the following.

- c: Clean option, if you want to run an other new benchmark, the old results will be cleaned up and a backup will be placed in 
in the directory backed_up_benchmark_results if you confirm that. Per default the cleanup will not be done and the old date will be overwritten.

- d: Dataset option, if you want to change the default dataset just set the flag `-d` and one of the corresponding keywords (`large`, `medium`, `small`). The default parameters are the one of the medium group. 

|   Parameter     |        small       |        medium      |                  large                   |
|-----------------|--------------------|--------------------|------------------------------------------|
|Dataset          |ERR016155.filt.fastq|ERR016155.filt.fastq|ERR251006.filt.fastq                      |
|Reference dataset|DRR001012.fa        |DRR001025.fa        |GRCh38_full_analysis_set_plus_decoy_hla.fa|
|Tensorflow steps |1000                |2500                |5000                                      |
|GROMACS steps    |10000               |30000               |50000                                     |

- p: Number of cores you want to use to run the benchmark. You can choose any arbitrary integer number or one of the three keywords `full`, `half` and `one`.
   - full: All cores which are accesible (Default)
   - half: Half of the cores that are accesible
   - one: Only one core (very long runtime)

- r: Number of replicates you want to use. The more replicates the better is the chance to get a trustworthy result and regulate outliers. You can choose any integer value you want. The default value is 3.

- t: Toolgroup you want to use. You can choose between `all` which will use all tools, `genomics` which will use only genomic tools (Bowtie2, Velvet, IDBA, SPAdes), `ml` which will only use Tensorflow and `quant` which will only use GROMACS. The default option is `all`. You can also choose just single tools. Currently these are `bowtie2-build`, `velvet`, `idba`, `tensorflow`, `gromacs`, `SPAdes`, `clustalomega`, which are also the keywords that this flag exceppts.

Some execution examples:

<pre>sh run_benchmarks.sh</pre>
This would be the same if we would run the command:
<pre>sh run_benchmarks.sh -d medium -p full -r 3 -t all</pre>

If you want to run a very long benchmark with maximal CPU usage you could use:
<pre>sh run_benchmarks.sh -d large -p full -r 5 -t all</pre>

If you have already executed a benchmark and want to run a new one use the `-c` flag:
<pre>sh run_benchmarks.sh -c -d medium -p full -r 3 -t all</pre>

If you just want to run some tensorflow benchmarks use for example:
<pre>sh run_benchmarks.sh -d large -p full -r 3 -t ml</pre>

You can also skip flags where you just want to use the default values.

After the benchmarks are finished you can create a Summary file in pdf format running the following command from BOOTABLE root directory:
<pre>Rscript --vanilla threadedBioBenchsuiteStatsGenerator.R</pre>

You will find the pdf file in the BOOTABLE root directory.

## Results
BOOTABLE will produce different result files. Mainly the most interesting file is the `benchmark_summary_\*.txt` file which you will find in the BOOTABLE root directory. The content looks like the following:
<pre>Replica_1 Bowtie2_build with 28 cores on dataset DRR001012
Mi 7. Nov 09:55:58 UTC 2018
real 211.35
user 3097.55
sys 14.13

Replica_1 Bowtie2_align with 28 cores on dataset ERR016155
Mi 7. Nov 09:59:29 UTC 2018
real 6.66
user 66.12
sys 7.77

Replica_1 Velveth with 28 cores on dataset ERR016155
Mi 7. Nov 09:59:36 UTC 2018
real 10.77
user 105.95
sys 6.74

Replica_1 Velvetg with 28 cores on dataset ERR016155
Mi 7. Nov 09:59:47 UTC 2018
real 64.57
user 155.62
sys 1.30

Replica_1 IDBA with 28 cores on dataset ERR016155
Mi 7. Nov 10:00:51 UTC 2018
real 110.99
user 2277.08
sys 10.19

Replica_1 Tensorflow with 28 cores on dataset cifar10 with 1000
Mi 7. Nov 10:02:42 UTC 2018
real 148.45
user 1934.74
sys 310.51

Replica_1 GROMACS with 28 cores on dataset adh_cubic calculating 10000 steps with CPU pinning enabled
Mi 7. Nov 10:05:12 UTC 2018
real 79.11
user 1707.00
sys 465.29

Replica_1 SPAdes with 28 cores on dataset ERR016155
Mi 7. Nov 10:06:31 UTC 2018
real 440.31
user 4642.39
sys 784.05</pre>

Every block gives you information about which replica (Replica_1) of which tool with how many cores and on which dataset or stepsize took how long. The different runtimes are stated in seconds.
- real: The real walltime the tool used from start to end 
- user: The used CPU time taking the number of cores into account
- sys:  The amount of CPU time spent in the kernel for system calls

Further you will find the file `benchmark_summary_\*.pdf` which gives you a more compact and statistical overview about the ran benhcmarks. You will find this file only if you have run the `threadedBioBenchsuiteStatsGenerator.R` script as stated above.
Also in the .pdf file all numbers are stated in seconds. Futher you can find a detailed list, tool by tool, in the results directory. There you can also find the stdout output (`benchmark_<tool>_output_<cores>.txt`) of the single tools and refollow the steps they have done.

Further you will find a summary of the system where the benchmark has been executed, which flags of BOOTABLE has been used and how Bowtie2 and GROMACS has been compiled in the file `bootable_system_info.txt`.

## Uninstall
In order to uninstall BOOTABLE just run the following commands from a directory one level higher than the BOOTABLE directory:
<pre>sudo rm -rf Bootable/
sudo rm -rf /usr/local/gromacs/
sudo pip uninstall tensorflow==1.4.0</pre>







