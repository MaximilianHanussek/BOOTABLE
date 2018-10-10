# BOOTABLE
BiOinfOrmatics ThreAded Benchmark tooLsuitE (BOOTABLE)
This toolsuite currently consists of the follwoing tools
- Bowtie2
- Velvet
- IDBA
- GROMACS
- Tensowflow
- SPAdes

more are coming soon ...
 
All the tools and datasets included in this bundle are carfeully chosen to 
cover a range of different application cases and make strong use of multithreading 
in order to get some benchmarks of the underlying hardware and about the scalability of the system.
In the following is explained how to install BOOTABLE, how to run it and what kind of options you have.

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
- zlib-devel.x86_64
- python-pip (in Version 9.0.3)
- cmake3
- tbb-devel.x86_64
- R

<pre>yum update -y
yum install epel-release -y
yum group install "Development Tools" -y
yum install nano curl wget vim htop time git zlib-devel.x86_64 python-pip cmake3 tbb-devel.x86_64 R -y</pre>


Further the system needs access to the internet to load the specific datasets and some R packages.

You can also install the required R packages in beforehand which are:
- grid
- gridExtra
- RColorBrewer
- stringr

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

Change into the cloned repo
<pre>cd BOOTABLE/</pre>

Run the installation script and all the tools will be compiled and installed in the BOOTABLE directory except some gromacs tools which will be installed under /usr/local/ at the moment and tensorflow which will be installed via `pip`
<pre>sh install_bootable.sh</pre>


### 2. Using a .qcow2 image
This instructions tell you how to use a .qcow2 image with everything preinstalled. It is further assumed that you know how to install an image in an virtuel environment or have a virtuel environment already running. Make sure you have at least 50GB of disk space. 
The image is deposited in an S3 bucket with the following URL `https://s3.denbi.uni-tuebingen.de/max/benchmark_image.qcow2`.
You can just download it into your current directory via `wget` for example:

<pre>wget https://s3.denbi.uni-tuebingen.de/max/benchmark_image.qcow2</pre>

Just start the image in any virtuel environment you like, for example with `virt-manager` under CentOS.
You can login with the following credentials:
- Username: root
- Password: rootdenbi

After the login please change the password of the user `centos` via the passwd command to anyone you like.
Afterwards logout and login to the virtual machine as user `centos` with the password you have entered before.
Change into the home directory (`/home/centos/`), where you will find all tools already installed. Start the benchmarks by running the following command:

<pre>sh run_benchmarks</pre>

This will start the calculations of the different tools first with the maximal number of CPU cores, second the half of the maximal number of CPU cores and at the end with one CPU core. Each number of CPU cores is running three times to get the mean value. As a comparison value, on a 28 core machine this will take some days (~80 hours).

After the benchmarks have finished run the following Rscript to generate a brief report in .pdf format for each number of benchmarked CPUs:

<pre>Rscript --vanilla report_generator.R</pre>

If you want to rerun the benchmarks use the -c flag in order to delete the created output files and benchmark results from the run before:

<pre>sh run_benchmarks -c</pre>



### 3. Using a Docker container
This instructions assume that docker is already installed and running. If not you should find most of the information
on the [Docker website](https://www.docker.com/get-started)

The first possibility is to pull the already configured BOOTABLE Docker image  from [docker hub](https://hub.docker.com/)
with the follwing docker pull command.

<pre>docker pull maximilianhanussek/bootable</pre>

After the image has been downloaded you can just start the benchmarking with the follwoing command

<pre>docker run --rm maximilianhanussek/bootable sh /root/run_benchmarks.sh </pre>

The second option is to build the Docker container by yourself with the provided `Dockerfile`. You can find it in the github repo
in the `docker` directory.

To build the Docker container by yourself clone this github repo, change into the `docker` directory and run the following command:

<pre>docker build --tag bootable .</pre>

The build will take a couple of minutes (30-60) as the first step is to update the underlying operating sytem /(centOS7) to the latest version and installing all the required packages. Further some large datasets need to be downloaded and this depends on your network conectivity. Afterwards most of the tools have to be compiled and installed.

After the Docker image has been build you can start the benchmark with the same command like if you had pulled it from Docker Hub:

<pre>docker run --rm maximilianhanussek/bootable sh /root/run_benchmarks.sh</pre>





### 4. Using a Singularity container

### 5. Transform an existing Docker container into a Singularity container

### 6. Using Ansible





