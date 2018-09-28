# BOOTABLE
BiOinfOrmatics ThreAded Benchmark tooLsuitE (BOOTABLE)
This toolsuite currently consists of the follwoing tools
- Bowtie2
- velvet
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
If you want to install the tools and scripts BOOTABLE use please be sure to have following packages installed via yum:
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
- python-pip
- cmake3
- tbb-devel.x86_64
- R

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
in between, clone this github repo and run the `install_bootable.sh` script

### 2. Using a .qcow2 image

### 3. Using a Docker container

### 4. Using a Singularity container

### 5. Transform an existing Docker container into a Singularity container

### 6. Using Ansible





