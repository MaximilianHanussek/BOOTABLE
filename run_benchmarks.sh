#!/bin/bash
# Run benchmarks

# Get information about the number of cores of the system 
max_cores=$(nproc)				#Get max core number
max_cores_velvet=$(expr $max_cores - 1)		#Get max core number for velvet (max core -1)
half_cores=$( expr $max_cores / 2)		#Calculate half core number
half_cores_velvet=$(expr $half_cores - 1)	#Calculate half core number for velvet 
quarter_cores=$(expr $max_cores / 4)		#Calculate quarter core number
quarter_cores_velvet=$(expr $quarter_cores - 1) #Calculate quarter core number fir velvet
one_core=1					#Set one core variable
clean=0						#Set clean flag initially to 0
dataset="ERR016155"
default_cores=$max_cores
default_cores_velvet=$max_cores_velvet


# Check if the -c flag is set or not for clean up operation
while getopts "cdp:" option; do
	case "${option}" in
		c) 
			clean=1
			;;
		d) 
			dataset=${OPTARG}
			;;
		p)
			default_cores=${OPTARG}
			default_cores_velvet=$(expr $default_cores - 1)
	esac
done

# Define function for benchmarking tools and there parameters
run_benchmark_tools () {
cores=$1					#Input parameter number cores (1,2,3,...)
cores_velvet=$2					#Input parameter number cores velvet (1,2,3,...)
replica=$3					#Input parameter number of replica (1,2,3,..)
dataset=$4					#Input parameter which dataset (medium, large)


# Bowtie2 build index
rm -rf benchmark_output/bowtie2/*		#Clean up bowtie2 output directoy
echo "Running bowtie2 index build benchmark"	
echo "Replica_$replica Bowtie2_build with $cores cores on dataset GRCh38" >> results/benchmark_bowtie_build_time_$cores.txt				 #Create results file with walltime
date >> results/benchmark_bowtie_build_time_$cores.txt	#Add date to walltime file

# Run bowtie2 index builder on dataset GRCh38 
/usr/bin/time -p -a -o results/benchmark_bowtie_build_time_$cores.txt sh -c "bowtie2/bowtie2-2.3.4.2/bowtie2-build --threads $cores --seed 42 datasets/1000_genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa benchmark_output/bowtie2/benchmark" >> results/benchmark_bowtie_build_output_$cores.txt 2>&1
echo "" >> results/benchmark_bowtie_build_time_$cores.txt	#Blank line for clarity and parsing

# Bowtie2 aligner
echo "Running bowtie2 align benchmark"
echo "Replica_$replica Bowtie2_align with $cores cores on dataset $dataset" >> results/benchmark_bowtie_align_time_$cores.txt				 #Create results file with walltime
date >> results/benchmark_bowtie_align_time_$cores.txt	#Add date to walltime file

# Run bowtie2 aligner on dataset $dataset
/usr/bin/time -p -a -o results/benchmark_bowtie_align_time_$cores.txt sh -c "bowtie2/bowtie2-2.3.4.2/bowtie2 --threads $cores -x benchmark_output/bowtie2/benchmark -U datasets/1000_genomes/$dataset.filt.fastq -S benchmark_output/bowtie2/benchmark_$dataset.sam" >> results/benchmark_bowtie_align_output_$cores.txt 2>&1
echo "" >> results/benchmark_bowtie_align_time_$cores.txt	#Blank line for clarity and parsing


# Velvet
echo "Running velvet benchmark"
rm -rf benchmark_output/velvet/*		#Clean up velvet output directory
echo "Replica_$replica Velveth with $cores cores on dataset $dataset" >> results/benchmark_velvet_time_$cores.txt					 #Create results file with walltime
date >> results/benchmark_velvet_time_$cores.txt	#Add date to walltime file

OMP_NUM_THREADS=$cores_velvet			#Set number of threads explicitly with OMP variable

# Run velveth on dataset $dataset
/usr/bin/time -p -a -o results/benchmark_velvet_time_$cores.txt sh -c "velvet/velveth benchmark_output/velvet/ 21 -fastq datasets/1000_genomes/$dataset.filt.fastq" >> results/benchmark_velveth_output_$cores.txt 2>&1
echo "" >> results/benchmark_velvet_time_$cores.txt		#Blank line for clarity and parsing


echo "Replica_$replica Velvetg with $cores cores on dataset $dataset" >> results/benchmark_velvet_time_$cores.txt					 #Create results file with walltime
date >> results/benchmark_velvet_time_$cores.txt #Add date to walltime file


# Run velvetg on dataset $dataset
/usr/bin/time -p -a -o results/benchmark_velvet_time_$cores.txt sh -c "velvet/velvetg benchmark_output/velvet/" >> benchmark_output/velvet/benchmark_velvetg_output.txt 2>&1
echo "" >> results/benchmark_velvet_time_$cores.txt 		#Blank line for clarity and parsing


# IDBA 
echo "Running IDBA benchmark"
rm -rf benchmark_output/IDBA/*
echo "Replica_$replica IDBA with $cores cores on dataset $dataset" >> results/benchmark_idba_time_$cores.txt
date >> results/benchmark_idba_time_$cores.txt

/usr/bin/time -p -a -o results/benchmark_idba_time_$cores.txt sh -c "IDBA/idba_ud-1.0.9/bin/idba_ud -r datasets/1000_genomes/$dataset.filt.fa --num_threads $cores -o benchmark_output/IDBA/" >> results/benchmark_idba_output_$cores.txt 2>&1
echo "" >> results/benchmark_idba_time_$cores.txt


echo "Running Tensorflow benchmark"
rm -rf benchmark_output/tensorflow/*
echo "Replica_$replica Tensorflow with $cores cores on dataset cifar10" >> results/benchmark_tensorflow_time_$cores.txt
date >> results/benchmark_tensorflow_time_$cores.txt

/usr/bin/time -p -a -o results/benchmark_tensorflow_time_$cores.txt sh -c "python datasets/tensorflow/models/tutorials/image/cifar10/cifar10_train.py --data_dir=datasets/tensorflow/ --train_dir=benchmark_output/tensorflow/cifar10_train --max_steps=5000 --threads=$cores" >> results/benchmark_tensorflow_output_$cores.txt 2>&1
echo "" >> results/benchmark_tensorflow_time_$cores.txt



echo "Running GROMACS benchmark"
rm -rf benchmark_output/gromacs/*
echo "Replica_$replica GROMACS with $cores cores on dataset adh_cubic calculating 50000 steps with CPU pinning enabled" >> results/benchmark_gromacs_time_$cores.txt
date >> results/benchmark_gromacs_time_$cores.txt

/usr/bin/time -p -a -o results/benchmark_gromacs_time_$cores.txt sh -c "/usr/local/gromacs/bin/gmx mdrun -v -pin on -nt $cores -s datasets/gromacs/adh_cubic/topol.tpr -o benchmark_output/gromacs/benchmark -cpo benchmark_output/gromacs/benchmark -e benchmark_output/gromacs/benchmark -g benchmark_output/gromacs/benchmark -c benchmark_output/gromacs/benchmark" >> results/benchmark_gromacs_output_$cores.txt 2>&1
echo "" >> results/benchmark_gromacs_time_$cores.txt


echo "Running SPAdes benchmark"
rm -rf benchmark_output/SPAdes/*
echo "Replica_$replica SPAdes with $cores cores on dataset $dataset" >> results/benchmark_SPAdes_time_$cores.txt
date >> results/benchmark_SPAdes_time_$cores.txt

/usr/bin/time -p -a -o results/benchmark_SPAdes_time_$cores.txt sh -c "python SPAdes/SPAdes-3.12.0-Linux/bin/spades.py -s datasets/1000_genomes/$dataset.filt.fastq -o benchmark_output/SPAdes/ -t $cores" >> results/benchmark_SPAdes_output_$cores.txt 2>&1
echo "" >> results/benchmark_SPAdes_time_$cores.txt

touch benchmark_summary_$cores.txt
cat results/benchmark_bowtie_build_time_$cores.txt >> benchmark_summary_$cores.txt

cat results/benchmark_bowtie_align_time_$cores.txt >> benchmark_summary_$cores.txt

cat results/benchmark_velvet_time_$cores.txt >> benchmark_summary_$cores.txt

cat results/benchmark_idba_time_$cores.txt >> benchmark_summary_$cores.txt

cat results/benchmark_tensorflow_time_$cores.txt >> benchmark_summary_$cores.txt

cat results/benchmark_gromacs_time_$cores.txt >> benchmark_summary_$cores.txt

cat results/benchmark_SPAdes_time_$cores.txt >> benchmark_summary_$cores.txt

rm results/benchmark_bowtie_build_time_$cores.txt

rm results/benchmark_bowtie_align_time_$cores.txt

rm results/benchmark_velvet_time_$cores.txt

rm results/benchmark_idba_time_$cores.txt

rm results/benchmark_tensorflow_time_$cores.txt

rm results/benchmark_gromacs_time_$cores.txt

rm results/benchmark_SPAdes_time_$cores.txt

}


if [ $clean == 1 ]
then
	rm benchmark_summary_*
	rm -rf results/*
fi

if [ $dataset == "large" ]
then
	dataset="ERR251006"
elif [ $dataset == "medium" ]
then
	dataset="ERR016155"
elif [ $dataset == "small" ]
then	
	dataset="" 
fi

echo $dataset
echo $default_cores
echo $default_cores_velvet


echo "First run with maximal core number and three replicates"
for replica in {1..3} 
do
	run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset
done

#echo "Second run with half of the maximal core number and three replicates"
#for replica in {1..3} 
#do
#	run_benchmark_tools $half_cores $half_cores_velvet $replica $dataset
#done

#echo "Third run with 1 core and three replicates"
#for replica in {1..3} 
#do
#	run_benchmark_tools $one_core $one_core $replica $dataset
#done

