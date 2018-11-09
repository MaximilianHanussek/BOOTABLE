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
clean=0 					#Set clean flag initially to 0
dataset="datasets/1000_genomes/ERR016155.filt.fastq" #Set ERR016155 as default dataset
dataset_idba="datasets/1000_genomes/ERR016155.filt.fa" #Set ERR016155 as default dataset for IDBA
default_reference="datasets/ebi/DRR001025.fa"   #Set DRR001012 as default reference dataset
default_tensorflow_steps=2500			#Set 2500 tensorflow steps as default 
default_gromacs_steps=30000			#Set 30000 gromacs steps as default
default_cores=$max_cores			#Set maximal core number as default
default_cores_velvet=$max_cores_velvet		#Set default core number for velvet tool
default_replicas=3				#Set default number of replicas to 3
default_toolgroup="all"
integer_regex='^[0-9]+$'			#Define regex for integers

# Save original PATH and LD_LIBRARY variables
original_path_variable=$(echo $PATH)
original_ld_library_variable=$(echo $LD_LIBRARY_PATH)

usage="$(basename "$0") [-h] [-cdprt]
where:
    -h  show this help text
    -c  clean up old benchmarks and back them up 
    -d  choose a dataset category (small, medium, large), medium will be default
    -p  number of cores/threads that should be used (one, half, full, or any integer value, full is default)
    -r  number of replica cycles that should be executed (any integer value, 3 is default)  
    -t  toolgroup which should be used for the benchmarks (all, genomics, ml, quant) or a specific tool (bowtie2-build, velvet, idba, tensorflow, gromacs, SPAdes). Default is all"

# Create flag options
while getopts "chd:p:r:t:" option; do
	case "${option}" in
		h)	echo "$usage"
			exit 1
			;;
		c) 
			clean=1
			;;
		d) 
			dataset=${OPTARG}
			;;
		p)
			default_cores=${OPTARG}
			;;
		r)
			default_replicas=${OPTARG}
			;;
		t)
			default_toolgroup=${OPTARG}
			;;
	esac
done

# Define function to check whether a file exists or not
check_results() {
filepath=$1
cores=$2

if [ -e $filepath ]
then
        cat $filepath >> benchmark_summary_$cores.txt
	rm $filepath
fi
}



# Define function for benchmarking tools and there parameters
run_benchmark_tools () {
cores=$1					#Input parameter number cores (1,2,3,...)
cores_velvet=$2					#Input parameter number cores velvet (1,2,3,...)
replica=$3					#Input parameter number of replica (1,2,3,..)
dataset=$4					#Input parameter which dataset (medium, large)
tf_steps=$5					#Input parameter for Tensorflow how much steps
gromacs_steps=$6				#Input parameter for GROMACS how much steps 
reference=$7					#Input parameter for bowtie build index which dataset
reference_name=$(basename $reference | cut -d. -f1)	#Get only the name of the dataset from the filepath
dataset_name=$(basename $dataset | cut -d. -f1) #Get only the name of the dataset from the filepath
dataset_idba=$8					#Input parameter for IDBA (needs .fa file)
toolgroup=$9					#Input parameter which tools should be used

# Save original PATH and LD_LIBRARY variables
original_path_variable=$(echo $PATH)
original_ld_library_variable=$(echo $LD_LIBRARY_PATH)

if [ $toolgroup == "all" ] || [ $toolgroup == "genomics" ] || [ $toolgroup == "bowtie2-build" ]
then
	# Bowtie2 build index
	rm -rf benchmark_output/bowtie2/*		#Clean up bowtie2 output directoy
	echo "Running bowtie2 index build benchmark on dataset $reference_name"	
	echo "Replica_$replica Bowtie2_build with $cores cores on dataset $reference_name" >> results/benchmark_bowtie_build_time_$cores.txt				 #Create results file with walltime
	date >> results/benchmark_bowtie_build_time_$cores.txt	#Add date to walltime file

	# Run bowtie2 index builder on dataset $reference_name
	/usr/bin/time -p -a -o results/benchmark_bowtie_build_time_$cores.txt sh -c "bowtie2/bowtie2-2.3.4.2/bowtie2-build --threads $cores --seed 42 $reference benchmark_output/bowtie2/benchmark" >> results/benchmark_bowtie_build_output_$cores.txt 2>&1
	echo "" >> results/benchmark_bowtie_build_time_$cores.txt	#Blank line for clarity and parsing
else
        echo "Bowtie2 index build will not be started as you did not choose the genomics tools, all tools or the tool itself."
fi

if [ $toolgroup == "all" ] || [ $toolgroup == "genomics" ]
then
	# Bowtie2 aligner
	echo "Running bowtie2 align benchmark on dataset $dataset_name"
	echo "Replica_$replica Bowtie2_align with $cores cores on dataset $dataset_name" >> results/benchmark_bowtie_align_time_$cores.txt				 #Create results file with walltime
	date >> results/benchmark_bowtie_align_time_$cores.txt	#Add date to walltime file

	# Run bowtie2 aligner on dataset $dataset
	/usr/bin/time -p -a -o results/benchmark_bowtie_align_time_$cores.txt sh -c "bowtie2/bowtie2-2.3.4.2/bowtie2 --threads $cores -x benchmark_output/bowtie2/benchmark -U $dataset -S benchmark_output/bowtie2/benchmark_$dataset_name.sam" >> results/benchmark_bowtie_align_output_$cores.txt 2>&1
	echo "" >> results/benchmark_bowtie_align_time_$cores.txt	#Blank line for clarity and parsing
else
	echo "Bowtie2 aligner will not be started as you did not choose the genomics tools, or all tools."
fi

if [ $toolgroup == "all" ] || [ $toolgroup == "genomics" ] || [ $toolgroup == "velvet" ]
then
	# Velvet
	echo "Running velvet benchmark on dataset $dataset_name"
	rm -rf benchmark_output/velvet/*		#Clean up velvet output directory
	echo "Replica_$replica Velveth with $cores cores on dataset $dataset_name" >> results/benchmark_velvet_time_$cores.txt					 #Create results file with walltime
	date >> results/benchmark_velvet_time_$cores.txt	#Add date to walltime file

	OMP_NUM_THREADS=$cores_velvet			#Set number of threads explicitly with OMP variable

	# Run velveth on dataset $dataset_name
	/usr/bin/time -p -a -o results/benchmark_velvet_time_$cores.txt sh -c "velvet/velveth benchmark_output/velvet/ 21 -fastq $dataset" >> results/benchmark_velveth_output_$cores.txt 2>&1
	echo "" >> results/benchmark_velvet_time_$cores.txt		#Blank line for clarity and parsing


	echo "Replica_$replica Velvetg with $cores cores on dataset $dataset_name" >> results/benchmark_velvet_time_$cores.txt					 #Create results file with walltime
	date >> results/benchmark_velvet_time_$cores.txt #Add date to walltime file


	# Run velvetg on dataset $dataset_name
	/usr/bin/time -p -a -o results/benchmark_velvet_time_$cores.txt sh -c "velvet/velvetg benchmark_output/velvet/" >> benchmark_output/velvet/benchmark_velvetg_output.txt 2>&1
	echo "" >> results/benchmark_velvet_time_$cores.txt 		#Blank line for clarity and parsing
else
        echo "Velvet will not be started as you did not choose the genomics tools, all tools or the tool itself."
fi

if [ $toolgroup == "all" ] || [ $toolgroup == "genomics" ] || [ $toolgroup == "idba" ]
then
	# IDBA 
	echo "Running IDBA benchmark on dataset $dataset_name"
	rm -rf benchmark_output/IDBA/*
	echo "Replica_$replica IDBA with $cores cores on dataset $dataset_name" >> results/benchmark_idba_time_$cores.txt
	date >> results/benchmark_idba_time_$cores.txt

	/usr/bin/time -p -a -o results/benchmark_idba_time_$cores.txt sh -c "IDBA/idba_ud-1.0.9/bin/idba_ud -r $dataset_idba --num_threads $cores -o benchmark_output/IDBA/" >> results/benchmark_idba_output_$cores.txt 2>&1
	echo "" >> results/benchmark_idba_time_$cores.txt
else
	echo "IDBA will not be started as you did not choose the genomics tools, all tools or the tool itself."
fi

if [ $toolgroup == "all" ] || [ $toolgroup == "ml" ] || [ $toolgroup == "tensorflow" ]
then
	# Tensorflow
	echo "Running Tensorflow benchmark with $tf_steps steps"
	rm -rf benchmark_output/tensorflow/*
	echo "Replica_$replica Tensorflow with $cores cores on dataset cifar10 with $tf_steps" >> results/benchmark_tensorflow_time_$cores.txt
	date >> results/benchmark_tensorflow_time_$cores.txt

	/usr/bin/time -p -a -o results/benchmark_tensorflow_time_$cores.txt sh -c "python datasets/tensorflow/models/tutorials/image/cifar10/cifar10_train.py --data_dir=datasets/tensorflow/ --train_dir=benchmark_output/tensorflow/cifar10_train --max_steps=$tf_steps --threads=$cores" >> results/benchmark_tensorflow_output_$cores.txt 2>&1
	echo "" >> results/benchmark_tensorflow_time_$cores.txt
else
	echo "Tensorflow will not be started as you did not choose the ml tools, all tools or the tool itself."
fi

if [ $toolgroup == "all" ] || [ $toolgroup == "quant" ] || [ $toolgroup == "gromacs" ]
then
	# GROMACS
	# Load correct compiler paths for GCC 7.3.0
	# Set GCC to 7.3.0 to use GROMACS
	export PATH=$PWD/gcc/gcc-installed/bin:$PATH
	export LD_LIBRARY_PATH=$PWD/gcc/gcc-installed/lib64:$LD_LIBRARY_PATH

	echo "Creating GROMACS test model with $gromacs_steps steps"
	sed -i "s/nsteps.*/nsteps                  = $gromacs_steps/g" datasets/gromacs/adh_cubic/pme_verlet.mdp

	/usr/local/gromacs/bin/gmx grompp -f datasets/gromacs/adh_cubic/pme_verlet.mdp -c datasets/gromacs/adh_cubic/conf.gro -p datasets/gromacs/adh_cubic/topol.top -o datasets/gromacs/adh_cubic/topol -po datasets/gromacs/adh_cubic/mdout >> /dev/null 2>&1

	echo "Running GROMACS benchmark with $gromacs_steps steps"
	rm -rf benchmark_output/gromacs/*
	echo "Replica_$replica GROMACS with $cores cores on dataset adh_cubic calculating $gromacs_steps steps with CPU pinning enabled" >> results/benchmark_gromacs_time_$cores.txt
	date >> results/benchmark_gromacs_time_$cores.txt

	/usr/bin/time -p -a -o results/benchmark_gromacs_time_$cores.txt sh -c "/usr/local/gromacs/bin/gmx mdrun -v -pin on -nt $cores -s datasets/gromacs/adh_cubic/topol.tpr -o benchmark_output/gromacs/benchmark -cpo benchmark_output/gromacs/benchmark -e benchmark_output/gromacs/benchmark -g benchmark_output/gromacs/benchmark -c benchmark_output/gromacs/benchmark" >> results/benchmark_gromacs_output_$cores.txt 2>&1
	echo "" >> results/benchmark_gromacs_time_$cores.txt

	# Reset to system compiler
	export PATH=$original_path_variable
	export LD_LIBRARY_PATH=$original_ld_library_variable
else
	echo "GROMACS will not be started as you did not choose the quant tools, all tools or the tool itself."
fi

if [ $toolgroup == "all" ] || [ $toolgroup == "genomics" ] || [ $toolgroup == "SPAdes" ]
then
	# SPAdes
	echo "Running SPAdes benchmark on dataset $dataset_name"
	rm -rf benchmark_output/SPAdes/*
	echo "Replica_$replica SPAdes with $cores cores on dataset $dataset_name" >> results/benchmark_SPAdes_time_$cores.txt
	date >> results/benchmark_SPAdes_time_$cores.txt

	/usr/bin/time -p -a -o results/benchmark_SPAdes_time_$cores.txt sh -c "python SPAdes/SPAdes-3.12.0-Linux/bin/spades.py -s $dataset -o benchmark_output/SPAdes/ -t $cores" >> results/benchmark_SPAdes_output_$cores.txt 2>&1
	echo "" >> results/benchmark_SPAdes_time_$cores.txt
else
	echo "SPAdes will not be started as you did not choose the genomics tools, all tools or the tool itself."
fi

touch benchmark_summary_$cores.txt

path="results/benchmark_bowtie_build_time_$cores.txt"
check_results $path $cores

path="results/benchmark_bowtie_align_time_$cores.txt"
check_results $path $cores

path="results/benchmark_velvet_time_$cores.txt"
check_results $path $cores

path="results/benchmark_idba_time_$cores.txt"
check_results $path $cores

path="results/benchmark_tensorflow_time_$cores.txt"
check_results $path $cores

path="results/benchmark_gromacs_time_$cores.txt"
check_results $path $cores

path="results/benchmark_SPAdes_time_$cores.txt"
check_results $path $cores
}


if [ $clean == 1 ]
then
	if [ -e benchmark_summary_*.txt ]
	then
		backup_date=$(stat -c %y benchmark_summary_*.txt | cut -d ' ' -f1)
		backup_time=$(stat -c %y benchmark_summary_*.txt | cut -d ' ' -f2 | cut -d. -f1)
		backup_dir_name="$backup_date-$backup_time"
		mkdir backed_up_benchmark_results/$backup_dir_name
		tar -cf backed_up_benchmark_results/$backup_dir_name/results.tar results/*
		cp benchmark_summary_* backed_up_benchmark_results/$backup_dir_name
		cp bootable_system_info.txt backed_up_benchmark_results/$backup_dir_name
		rm benchmark_summary_*
		rm -rf results/*
		rm bootable_system_info.txt
	else
		while true; do
                read -p "There are no files to back them up, do you want to start the benchmark?" yn
                case $yn in
                        [Yy]* ) break;;
                        [Nn]* ) echo "The benchmarks will not be started"
				exit 1
                                ;;
                        * ) echo "Please answer yes or no."
                                ;;
                esac
        	done
	fi

fi

touch bootable_system_info.txt
date=$(date)
if [ $clean == 1 ]
then
        used_command="-c -d $dataset -p $default_cores -r $default_replicas -t $default_toolgroup"
else
        used_command="-d $dataset -p $default_cores -r $default_replicas -t $default_toolgroup"
fi

echo "$date" > bootable_system_info.txt
echo "" >> bootable_system_info.txt
echo "Executed command: run_benchmarks.sh $used_command" >> bootable_system_info.txt
echo "" >> bootable_system_info.txt
echo "System information:" >> bootable_system_info.txt
inxi -C -f -M -m -S -I -D -x >> bootable_system_info.txt

if [ -e /usr/sbin/tuned-adm ]
then
	tuned_out=$(tuned-adm active)
	echo "tuned status: $tuned_out" >> bootable_system_info.txt 
else 
	echo "tuned status: NOT installed." >> bootable_system_info.txt
fi

lscpu -p='Core' | grep -v ^# | sort | uniq -c | awk '{print $1}' | uniq -c | while read -r no_cores threads ;
do
        if [ "$threads" -eq 1 ] ; then
                ht="disabled"
        else
                ht="enabled"
        fi
        echo "Hyperthreading: $ht ($threads thread(s)) per core" >> bootable_system_info.txt
done
echo "" >> bootable_system_info.txt
echo "Bowtie2 compile information:" >> bootable_system_info.txt
bowtie2/bowtie2-2.3.4.2/bowtie2 --version >> bootable_system_info.txt
echo "" >> bootable_system_info.txt
echo "GROMACS compile information:" >> bootable_system_info.txt

# Set GCC to 7.3.0 to use GROMACS
export PATH=$PWD/gcc/gcc-installed/bin:$PATH
export LD_LIBRARY_PATH=$PWD/gcc/gcc-installed/lib64:$LD_LIBRARY_PATH
gromacs_line=$(/usr/local/gromacs/bin/gmx --version | grep -n "GROMACS version:" | cut -d ":" -f 1)
/usr/local/gromacs/bin/gmx --version | tail --lines=+$gromacs_line >> bootable_system_info.txt

# Reset to system compiler
export PATH=$original_path_variable
export LD_LIBRARY_PATH=$original_ld_library_variable

if [ $dataset == "large" ]
then
	dataset="datasets/1000_genomes/ERR251006.filt.fastq"
	dataset_idba="datasets/1000_genomes/ERR251006.filt.fa"
        default_reference="datasets/1000_genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
	default_tensorflow_steps=5000
	default_gromacs_steps=50000

elif [ $dataset == "medium" ]
then
	dataset="datasets/1000_genomes/ERR016155.filt.fastq"
        dataset_idba="datasets/1000_genomes/ERR016155.filt.fa"
	default_reference="datasets/ebi/DRR001025.fa"
	default_tensorflow_steps=2500
	default_gromacs_steps=30000

elif [ $dataset == "small" ]
then	
	dataset="datasets/1000_genomes/ERR016155.filt.fastq"
	dataset_idba="datasets/1000_genomes/ERR016155.filt.fa"
	default_reference="datasets/ebi/DRR001012.fa" 
	default_tensorflow_steps=1000
	default_gromacs_steps=10000

else
	dataset="datasets/1000_genomes/ERR016155.filt.fastq"
        dataset_idba="datasets/1000_genomes/ERR016155.filt.fa"
        default_reference="datasets/ebi/DRR001025.fa"
        default_tensorflow_steps=2500
        default_gromacs_steps=30000

        echo "Parameter is not one of small, medium or large. Please check -d flag again. Default settings will be used (medium)"
fi


if [ $default_cores == "full" ]
then
        default_cores=$max_cores
        default_cores_velvet=$max_cores_velvet

elif [ $default_cores == "half" ]
then
        default_cores=$half_cores
	default_cores_velvet=$half_cores_velvet

elif [ $default_cores == "one" ]
then
        default_cores=$one_core
        default_cores_velvet=$one_core

elif [[ "$default_cores" =~ $integer_regex ]]
then
	default_cores_velvet=$(expr $default_cores - 1)

else
	echo "Parameter is not one of full, half, one or an integer number. Please check -p flag again."
	exit 1
fi

if [[ $default_toolgroup != "all" && $default_toolgroup != "genomics" && $default_toolgroup != "ml" && $default_toolgroup != "quant" && $default_toolgroup != "bowtie2-build" && $default_toolgroup != "velvet" && $default_toolgroup != "idba" && $default_toolgroup != "tensorflow" && $default_toolgroup != "gromacs" && $default_toolgroup != "SPAdes" ]]
then
	echo "Parameter is not one of all, genomics, ml, quant, bowtie2-build, velvet, idba, tensorflow, gromacs or SPAdes. Please check -t flag again."
	exit 1
fi



echo "General genomic dataset for Bowtie2, Velvet, SPAdes: $dataset"
echo "IDBA dataset: $dataset_idba"
echo "Reference dataset: $default_reference"
echo "Number of used cores: $default_cores"
echo "Number of used replicates: $default_replicas"
echo "Number of Tensorflow steps: $default_tensorflow_steps"
echo "Number of GROMACS steps: $default_gromacs_steps"
echo "Toolgroup: $default_toolgroup"

echo "BOOTABLE benchmark run with $default_cores cores and $default_replicas replicates"
for replica in $( seq 1 $default_replicas ) 
do
	run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup
done


