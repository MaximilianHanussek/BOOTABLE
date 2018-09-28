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
 
All the tools are included in this bundle and are carfeully chosen to 
cover a range of different application cases and make strong use multithreading 
in order to get some benchmarks of the underlying hardware. In the follwoing is explained 
how to install BOOTABLE, how to run it and what kind of options yopu have.


## Installation

There are different kinds of installation possibilities.
1. Bare metal installation
2. Using a .qcow2 image with everything already preinstalled
3. Using Docker and the provided Dockerfile to install everything into a Docker Container
4. Using Singularity and the provided Singularity file
5. Transform the Docker container from option 3. into a Singualrity Container
6. Use Ansible and the provided Ansible script to install the tools on a CentOS based system from scratch
