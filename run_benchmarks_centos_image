#!/bin/bash
#run_benchmarks

max_cores=$(nproc)
max_cores_velvet=$(expr $max_cores - 1)
half_cores=$( expr $max_cores / 2)
half_cores_velvet=$(expr $half_cores - 1)
quarter_cores=$(expr $max_cores / 4)
one_core=1
clean=0

while getopts c option
do
case "${option}"
in
c) clean=1;;
esac
done

run_benchmark_tools () {
cores=$1
cores_velvet=$2
replica=$3

rm -rf /home/centos/benchmark_output/bowtie2/*
echo "Running bowtie2 index build benchmark"
echo "Replica_$replica Bowtie2_build with $cores cores on dataset ERR251006" >> /home/centos/results/benchmark_bowtie_build_time_$cores.txt
date >> /home/centos/results/benchmark_bowtie_build_time_$cores.txt

/usr/bin/time -p -a -o /home/centos/results/benchmark_bowtie_build_time_$cores.txt sh -c "/home/centos/bowtie2/bowtie2-2.3.4.2/bowtie2-build --threads $cores --seed 42 /home/centos/datasets/1000_genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/centos/benchmark_output/bowtie2/benchmark" >> /home/centos/results/benchmark_bowtie_build_output_$cores.txt 2>&1
echo "" >> /home/centos/results/benchmark_bowtie_build_time_$cores.txt


echo "Running bowtie2 align benchmark"
echo "Replica_$replica Bowtie2_align with $cores cores on dataset ERR251006" > /home/centos/results/benchmark_bowtie_align_time_$cores.txt
date >> /home/centos/results/benchmark_bowtie_align_time_$cores.txt

/usr/bin/time -p -a -o /home/centos/results/benchmark_bowtie_align_time_$cores.txt sh -c "/home/centos/bowtie2/bowtie2-2.3.4.2/bowtie2 --threads $cores -x /home/centos/benchmark_output/bowtie2/benchmark -U /home/centos/datasets/1000_genomes/ERR251006.filt.fastq -S /home/centos/benchmark_output/bowtie2/benchmark_ERR251006.sam" >> /home/centos/results/benchmark_bowtie_align_output_$cores.txt 2>&1
echo "" >> /home/centos/results/benchmark_bowtie_align_time_$cores.txt


echo "Running velvet benchmark"
rm -rf /home/centos/benchmark_output/velvet/*
echo "Replica_$replica Velveth with $cores cores on dataset ERR251006" >> /home/centos/results/benchmark_velvet_time_$cores.txt
date >> /home/centos/results/benchmark_velvet_time_$cores.txt

OMP_NUM_THREADS=$cores_velvet

/usr/bin/time -p -a -o /home/centos/results/benchmark_velvet_time_$cores.txt sh -c "/home/centos/velvet/velveth /home/centos/benchmark_output/velvet/ 21 -fastq /home/centos/datasets/1000_genomes/ERR251006.filt.fastq" >> /home/centos/results/benchmark_velveth_output_$cores.txt 2>&1
echo "" >> /home/centos/results/benchmark_velvet_time_$cores.txt


echo "Replica_$replica Velvetg with $cores cores on dataset ERR251006" >> /home/centos/results/benchmark_velvet_time_$cores.txt
date >> /home/centos/results/benchmark_velvet_time_$cores.txt

/usr/bin/time -p -a -o /home/centos/results/benchmark_velvet_time_$cores.txt sh -c "/home/centos/velvet/velvetg /home/centos/benchmark_output/velvet/" >> /home/centos/benchmark_output/velvet/benchmark_velvetg_output.txt 2>&1
echo "" >> /home/centos/results/benchmark_velvet_time_$cores.txt


echo "Running IDBA benchmark"
rm -rf /home/centos/benchmark_output/IDBA/*
echo "Replica_$replica IDBA with $cores cores on dataset ERR251006" >> /home/centos/results/benchmark_idba_time_$cores.txt
date >> /home/centos/results/benchmark_idba_time_$cores.txt

/usr/bin/time -p -a -o /home/centos/results/benchmark_idba_time_$cores.txt sh -c "/home/centos/IDBA/idba_ud-1.0.9/bin/idba_ud -r /home/centos/datasets/1000_genomes/ERR251006.filt.fa --num_threads $cores -o /home/centos/benchmark_output/IDBA/" >> /home/centos/results/benchmark_idba_output_$cores.txt 2>&1
echo "" >> /home/centos/results/benchmark_idba_time_$cores.txt


echo "Running Tensoflow benchmark"
rm -rf /home/centos/benchmark_output/tensorflow/*
echo "Replica_$replica Tensorflow with $cores cores on dataset cifar10" >> /home/centos/results/benchmark_tensorflow_time_$cores.txt
date >> /home/centos/results/benchmark_tensorflow_time_$cores.txt

/usr/bin/time -p -a -o /home/centos/results/benchmark_tensorflow_time_$cores.txt sh -c "python /home/centos/datasets/tensorflow/models/tutorials/image/cifar10/cifar10_train.py --data_dir=/home/centos/datasets/tensorflow/ --train_dir=/home/centos/benchmark_output/tensorflow/cifar10_train --max_steps=5000 --threads=$cores" >> /home/centos/results/benchmark_tensorflow_output_$cores.txt 2>&1
echo "" >> /home/centos/results/benchmark_tensorflow_time_$cores.txt



echo "Running GROMACS benchmark"
rm -rf /home/centos/benchmark_output/gromacs/*
echo "Replica_$replica GROMACS with $cores cores on dataset adh_cubic calculating 50000 steps with CPU pinning enabled" >> /home/centos/results/benchmark_gromacs_time_$cores.txt
date >> /home/centos/results/benchmark_gromacs_time_$cores.txt

/usr/bin/time -p -a -o /home/centos/results/benchmark_gromacs_time_$cores.txt sh -c "/usr/local/gromacs/bin/gmx mdrun -v -pin on -nt $cores -s /home/centos/datasets/gromacs/adh_cubic/topol.tpr -o /home/centos/benchmark_output/gromacs/benchmark -cpo /home/centos/benchmark_output/gromacs/benchmark -e /home/centos/benchmark_output/gromacs/benchmark -g /home/centos/benchmark_output/gromacs/benchmark -c /home/centos/benchmark_output/gromacs/benchmark" >> /home/centos/results/benchmark_gromacs_output_$cores.txt 2>&1
echo "" >> /home/centos/results/benchmark_gromacs_time_$cores.txt


echo "Running SPAdes benchmark"
rm -rf /home/centos/benchmark_output/SPAdes/*
echo "Replica_$replica SPAdes with $cores cores on dataset ERR251006" >> /home/centos/results/benchmark_SPAdes_time_$cores.txt
date >> /home/centos/results/benchmark_SPAdes_time_$cores.txt

/usr/bin/time -p -a -o /home/centos/results/benchmark_SPAdes_time_$cores.txt sh -c "python /home/centos/SPAdes/SPAdes-3.12.0-Linux/bin/spades.py -s /home/centos/datasets/1000_genomes/ERR251006.filt.fastq -o /home/centos/benchmark_output/SPAdes/ -t $cores" >> /home/centos/results/benchmark_SPAdes_output_$cores.txt 2>&1
echo "" >> /home/centos/results/benchmark_SPAdes_time_$cores.txt

touch /home/centos/benchmark_summary_$cores.txt
cat /home/centos/results/benchmark_bowtie_build_time_$cores.txt >> /home/centos/benchmark_summary_$cores.txt

cat /home/centos/results/benchmark_bowtie_align_time_$cores.txt >> /home/centos/benchmark_summary_$cores.txt

cat /home/centos/results/benchmark_velvet_time_$cores.txt >> /home/centos/benchmark_summary_$cores.txt

cat /home/centos/results/benchmark_idba_time_$cores.txt >> /home/centos/benchmark_summary_$cores.txt

cat /home/centos/results/benchmark_tensorflow_time_$cores.txt >> /home/centos/benchmark_summary_$cores.txt

cat /home/centos/results/benchmark_gromacs_time_$cores.txt >> /home/centos/benchmark_summary_$cores.txt

cat /home/centos/results/benchmark_SPAdes_time_$cores.txt >> /home/centos/benchmark_summary_$cores.txt

rm /home/centos/results/benchmark_bowtie_build_time_$cores.txt

rm /home/centos/results/benchmark_bowtie_align_time_$cores.txt

rm /home/centos/results/benchmark_velvet_time_$cores.txt

rm /home/centos/results/benchmark_idba_time_$cores.txt

rm /home/centos/results/benchmark_tensorflow_time_$cores.txt

rm /home/centos/results/benchmark_gromacs_time_$cores.txt

rm /home/centos/results/benchmark_SPAdes_time_$cores.txt

}


if [ $clean == 1 ]
then
	rm /home/centos/benchmark_summary_*
	rm -rf /home/centos/results/*
fi

echo "First run with maximal core number and three replicates"
for replica in {1..3} 
do
	run_benchmark_tools $max_cores $max_cores_velvet $replica
done

echo "Second run with half of the maximal core number and three replicates"
for replica in {1..3} 
do
	run_benchmark_tools $half_cores $half_cores_velvet $replica
done

echo "Third run with 1 core and three replicates"
for replica in {1..3} 
do
	run_benchmark_tools $one_core $one_core $replica
done

