#!/bin/bash -x
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
scaling="none"					#Set scaling initally to "none"
own_tool_path="none"				#Set path for own tools per default to "none"
dataset="datasets/1000_genomes/ERR016155.filt.fastq" #Set ERR016155 as default dataset
dataset_idba="datasets/1000_genomes/ERR015528.filt.fa" #Set ERR015528 as default dataset for IDBA
default_reference="datasets/ebi/DRR001025.fa"   #Set DRR001012 as default reference dataset
default_tensorflow_steps=2500			#Set 2500 tensorflow steps as default 
default_gromacs_steps=30000			#Set 30000 gromacs steps as default
dataset_clustalOmega="datasets/clustalOmega/wgs.ANCA.1_400.fsa"	#Set wgs.ANCA.1_400.fsa as default datatset for ClustalOmega
default_SINA="datasets/SINA/RefSeq-RDP16S_v2_May2018.fa"	#Set RefSeq-RDP16S_v2_May2018.fa as default dataset for SINA
reference_SINA="datasets/SINA/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb"	#Set SILVA_138.1_SSURef_NR99_12_06_20_opt.arb as default dataset for SINA
default_reference_BWA="datasets/BWA/DRR001025/DRR001025"	#Set DRR001025 as default reference dataset for BWA, index files
default_cores=$max_cores			#Set maximal core number as default
default_cores_velvet=$max_cores_velvet		#Set default core number for velvet tool
default_replicas=3				#Set default number of replicas to 3
default_toolgroup="all"
integer_regex='^[0-9]+$'			#Define regex for integers

# Save original PATH and LD_LIBRARY variables
original_path_variable=$(echo $PATH)
original_ld_library_variable=$(echo $LD_LIBRARY_PATH)

usage="$(basename "$0") [-h] [-cdprst]
where:
    -h  show this help text
    -c  clean up old benchmarks and back them up 
    -d  choose a dataset category (small, medium, large), medium will be default
    -o  specify full path to a toolfile to test tools not part of bootable (/home/user/toolfile.btbl) 
    -p  number of cores/threads that should be used (one, half, full, or any integer value, full is default)
    -r  number of replica cycles that should be executed (any integer value, 3 is default)  
    -s  run scaling benchmark (small, medium, large), parameters describe the used datasets
    -t  toolgroup which should be used for the benchmarks (all, genomics, ml, quant) or a specific tool (bbmap, bowtie2-build, bwa, velvet, idba, tensorflow, gromacs, SPAdes, clustalomega, mafft, sina). Default is all"

# Create flag options
while getopts "chd:o:p:r:s:t:" option; do
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
		o)
			own_tool_path=${OPTARG}
			;;
		p)
			default_cores=${OPTARG}
			;;
		r)
			default_replicas=${OPTARG}
			;;
		s)	
			scaling=${OPTARG}
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
dataset=$4					#Input parameter which dataset (small, medium, large)
tf_steps=$5					#Input parameter for Tensorflow how much steps
gromacs_steps=$6				#Input parameter for GROMACS how much steps 
reference=$7					#Input parameter for bowtie build index which dataset
reference_name=$(basename $reference | cut -d. -f1)	#Get only the name of the dataset from the filepath
dataset_name=$(basename $dataset | cut -d. -f1) #Get only the name of the dataset from the filepath
dataset_idba=$8					#Input parameter for IDBA (needs .fa file)
dataset_name_idba=$(basename $dataset_idba | cut -d. -f1) #Get only the name of the dataset for idba from the filepath
toolgroup=$9					#Input parameter which tools should be used
dataset_clustalOmega=${10}                      #Input parameter for Clustal Omega which dataset
dataset_name_clustalOmega=$(basename $dataset_clustalOmega | cut -d. -f-3)
reference_BWA=${11}				#Input parameter for BWA which reference dataset
reference_name_BWA=$(basename $reference_BWA)	#Get only the name of the dataset from the filepath
dataset_SINA=${12}				#Input parameter for SINA which dataset
dataset_name_SINA=$(basename $dataset_SINA | cut -d. -f1) #Get only the name of the dataset from the filepath
reference_SINA=${13}				#Input parameter for SINA which reference dataset
reference_name_SINA=$(basename $reference_SINA | cut -d. -f-2) #Get only the name of the dataset from the filepath
own_tool_path=${14}				#Input parameter path to toolfile (full path)
if [  "$own_tool_path" != "none" ]
then

own_toolname=$(cat "$own_tool_path" | grep "Toolname" | cut -d: -f2)
own_datasetname=$(cat "$own_tool_path" | grep "Dataset" | cut -d: -f2)
own_toolcommand=$(cat "$own_tool_path" | grep "Command" | cut -d: -f2)
substituted_own_toolcommand=$(eval echo $own_toolcommand)

if [ -d benchmark_output/"$own_toolname" ]
then
        echo "Directory benchmark_output for ""$own_toolname"" already exists."
else
        mkdir benchmark_output/"$own_toolname"
fi

#Own tool execution
rm -rf benchmark_output/"$own_toolname"/*
echo "Running ""$own_toolname"" index build benchmark on dataset ""$own_datasetname"" " 
echo "Replica_$replica ""$own_toolname"" with ""$cores"" cores on dataset ""$own_datasetname"" " >> results/benchmark_""$own_toolname""_time_""$cores"".txt #Create results file with walltime
date >> results/benchmark_"$own_toolname"_time_"$cores".txt  #Add date to walltime file

# Start nmon capturing
NMON_FILE_NAME="""$own_toolname""_""$replica""_$(date +"%Y-%m-%d-%H-%M")"
NMON_PID=$(nmon -F "$NMON_FILE_NAME".nmon -m nmon_stats/ -p -s 2 -c 12000000) # -s is the interval between snapshots, -c is the number of snapshots very high to ensure them to the end of the tool run

# Run own tool
/usr/bin/time -p -a -o results/benchmark_"$own_toolname"_time_$cores.txt sh -c "$substituted_own_toolcommand" >> results/benchmark_"$own_toolname"_output_$cores.txt 2>&1

# Stop nmon capturing
kill -USR2 "$NMON_PID"

# Make nmon html graphs
sh nmonchart/nmonchart nmon_stats/"$NMON_FILE_NAME".nmon nmon_stats/"$NMON_FILE_NAME".html

echo "" >> results/benchmark_"$own_toolname"_time_"$cores".txt       #Blank line for clarity and parsing

touch benchmark_summary_"$cores".txt

path="results/benchmark_""$own_toolname""_time_""$cores"".txt"
check_results $path $cores

else
# Save original PATH and LD_LIBRARY variables
original_path_variable=$(echo $PATH)
original_ld_library_variable=$(echo $LD_LIBRARY_PATH)

if [ $toolgroup == "all" ] || [ $toolgroup == "genomics" ] || [ $toolgroup == "bbmap" ]
then
        # BBMap index and mapping
        rm -rf benchmark_output/BBMap/*               #Clean up BBMap output directoy
        echo "Running BBMap benchmark on reference dataset $reference_name and dataset $dataset_name" 
        echo "Replica_$replica BBMap with $cores cores on dataset $reference_name and dataset $dataset_name" >> results/benchmark_bbmap_time_$cores.txt                             #Create results file with walltime
        date >> results/benchmark_bbmap_time_$cores.txt  #Add date to walltime file

        # Start nmon capturing
        NMON_FILE_NAME="bbmap_"$replica"_$(date +"%Y-%m-%d-%H-%M")"
        NMON_PID=$(nmon -F $NMON_FILE_NAME.nmon -m nmon_stats/ -p -s 2 -c 12000000) # -s is the interval between snapshots, -c is the number of snapshots very high to ensure them to the end of the tool run

        # Run BBMap on reference dataset $reference_name and dataset $dataset_name
        /usr/bin/time -p -a -o results/benchmark_bbmap_time_$cores.txt sh -c "BBMap/bbmap/bbmap.sh threads=$cores in=$dataset out=benchmark_output/BBMap/benchmark ref=$reference path=benchmark_output/BBMap/" >> results/benchmark_bbmap_output_$cores.txt 2>&1

        # Stop nmon capturing
        kill -USR2 $NMON_PID

        # Make nmon html graphs
        sh nmonchart/nmonchart nmon_stats/$NMON_FILE_NAME.nmon nmon_stats/$NMON_FILE_NAME.html
        echo "" >> results/benchmark_bbmap_time_$cores.txt       #Blank line for clarity and parsing
else
        echo "BBMap will not be started as you did not choose the genomics tools, all tools or the tool itself."
fi

if [ $toolgroup == "all" ] || [ $toolgroup == "genomics" ] || [ $toolgroup == "bwa" ]
then
        # BWA MEM
        rm -rf benchmark_output/BWA/*               #Clean up BWA output directoy
        echo "Running BWA MEM benchmark on reference dataset $reference_name_BWA and dataset $dataset_name"
        echo "Replica_$replica BWA_MEM with $cores cores on reference dataset $reference_name_BWA and dataset $dataset_name" >> results/benchmark_bwa_mem_time_$cores.txt                             #Create results file with walltime
        date >> results/benchmark_bwa_mem_time_$cores.txt  #Add date to walltime file

        # Start nmon capturing
        NMON_FILE_NAME="bwa_mem_"$replica"_$(date +"%Y-%m-%d-%H-%M")"
        NMON_PID=$(nmon -F $NMON_FILE_NAME.nmon -m nmon_stats/ -p -s 2 -c 12000000) # -s is the interval between snapshots, -c is the number of snapshots very high to ensure them to the end of the tool run

        # Run BWA Mem on dataset $reference_name_BWA
        /usr/bin/time -p -a -o results/benchmark_bwa_mem_time_$cores.txt sh -c "BWA/bwa-0.7.17/bwa mem -t $cores $reference_BWA $dataset > benchmark_output/BWA/benchmark" >> results/benchmark_bwa_mem_output_$cores.txt 2>&1

        # Stop nmon capturing
        kill -USR2 $NMON_PID

        # Make nmon html graphs
        sh nmonchart/nmonchart nmon_stats/$NMON_FILE_NAME.nmon nmon_stats/$NMON_FILE_NAME.html
        echo "" >> results/benchmark_bwa_mem_time_$cores.txt       #Blank line for clarity and parsing
else
        echo "BWA MEM will not be started as you did not choose the genomics tools, all tools or the tool itself."
fi

if [ $toolgroup == "all" ] || [ $toolgroup == "genomics" ] || [ $toolgroup == "bowtie2-build" ]
then
	# Bowtie2 build index
	rm -rf benchmark_output/bowtie2/*		#Clean up bowtie2 output directoy
	echo "Running bowtie2 index build benchmark on dataset $reference_name"	
	echo "Replica_$replica Bowtie2_build with $cores cores on dataset $reference_name" >> results/benchmark_bowtie_build_time_$cores.txt				 #Create results file with walltime
	date >> results/benchmark_bowtie_build_time_$cores.txt	#Add date to walltime file

	# Start nmon capturing
	NMON_FILE_NAME="bowtie2_build_"$replica"_$(date +"%Y-%m-%d-%H-%M")"
	NMON_PID=$(nmon -F $NMON_FILE_NAME.nmon -m nmon_stats/ -p -s 2 -c 12000000) # -s is the interval between snapshots, -c is the number of snapshots very high to ensure them to the end of the tool run

	# Run bowtie2 index builder on dataset $reference_name
	/usr/bin/time -p -a -o results/benchmark_bowtie_build_time_$cores.txt sh -c "bowtie2/bowtie2-2.4.2/bowtie2-build --threads $cores --seed 42 $reference benchmark_output/bowtie2/benchmark" >> results/benchmark_bowtie_build_output_$cores.txt 2>&1

	# Stop nmon capturing
        kill -USR2 $NMON_PID

	# Make nmon html graphs
	sh nmonchart/nmonchart nmon_stats/$NMON_FILE_NAME.nmon nmon_stats/$NMON_FILE_NAME.html 
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

	# Start nmon capturing
	NMON_FILE_NAME="bowtie2_align_"$replica"_$(date +"%Y-%m-%d-%H-%M")"
        NMON_PID=$(nmon -F $NMON_FILE_NAME.nmon -m nmon_stats/ -p -s 2 -c 12000000)

	# Run bowtie2 aligner on dataset $dataset
	/usr/bin/time -p -a -o results/benchmark_bowtie_align_time_$cores.txt sh -c "bowtie2/bowtie2-2.4.2/bowtie2 --threads $cores -x benchmark_output/bowtie2/benchmark -U $dataset -S benchmark_output/bowtie2/benchmark_$dataset_name.sam" >> results/benchmark_bowtie_align_output_$cores.txt 2>&1

	# Stop nmon capturing
        kill -USR2 $NMON_PID

	# Make nmon html graphs
        sh nmonchart/nmonchart nmon_stats/$NMON_FILE_NAME.nmon nmon_stats/$NMON_FILE_NAME.html

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

	# Start nmon capturing
	NMON_FILE_NAME="velveth_"$replica"_$(date +"%Y-%m-%d-%H-%M")"
        NMON_PID=$(nmon -F $NMON_FILE_NAME.nmon -m nmon_stats/ -p -s 2 -c 12000000)

	# Run velveth on dataset $dataset_name
	/usr/bin/time -p -a -o results/benchmark_velvet_time_$cores.txt sh -c "velvet/velveth benchmark_output/velvet/ 21 -fastq $dataset" >> results/benchmark_velveth_output_$cores.txt 2>&1
	
	# Stop nmon capturing
        kill -USR2 $NMON_PID

	# Make nmon html graphs
        sh nmonchart/nmonchart nmon_stats/$NMON_FILE_NAME.nmon nmon_stats/$NMON_FILE_NAME.html
	
	echo "" >> results/benchmark_velvet_time_$cores.txt		#Blank line for clarity and parsing

	echo "Replica_$replica Velvetg with $cores cores on dataset $dataset_name" >> results/benchmark_velvet_time_$cores.txt					 #Create results file with walltime
	date >> results/benchmark_velvet_time_$cores.txt #Add date to walltime file

	# Start nmon capturing
	NMON_FILE_NAME="velvetg_"$replica"_$(date +"%Y-%m-%d-%H-%M")"
        NMON_PID=$(nmon -F $NMON_FILE_NAME.nmon -m nmon_stats/ -p -s 2 -c 12000000)

	# Run velvetg on dataset $dataset_name
	/usr/bin/time -p -a -o results/benchmark_velvet_time_$cores.txt sh -c "velvet/velvetg benchmark_output/velvet/" >> benchmark_output/velvet/benchmark_velvetg_output.txt 2>&1
	
	# Stop nmon capturing
        kill -USR2 $NMON_PID

	# Make nmon html graphs
        sh nmonchart/nmonchart nmon_stats/$NMON_FILE_NAME.nmon nmon_stats/$NMON_FILE_NAME.html

	echo "" >> results/benchmark_velvet_time_$cores.txt 		#Blank line for clarity and parsing
else
        echo "Velvet will not be started as you did not choose the genomics tools, all tools or the tool itself."
fi

if [ $toolgroup == "all" ] || [ $toolgroup == "genomics" ] || [ $toolgroup == "idba" ]
then
	# IDBA 
	echo "Running IDBA benchmark on dataset $dataset_name"
	rm -rf benchmark_output/IDBA/*
	echo "Replica_$replica IDBA with $cores cores on dataset $dataset_name_idba" >> results/benchmark_idba_time_$cores.txt
	date >> results/benchmark_idba_time_$cores.txt

	# Start nmon capturing
	NMON_FILE_NAME="IDBA_"$replica"_$(date +"%Y-%m-%d-%H-%M")"
        NMON_PID=$(nmon -F $NMON_FILE_NAME.nmon -m nmon_stats/ -p -s 2 -c 12000000)

	/usr/bin/time -p -a -o results/benchmark_idba_time_$cores.txt sh -c "IDBA/idba_ud-1.0.9/bin/idba_ud -r $dataset_idba --num_threads $cores -o benchmark_output/IDBA/" >> results/benchmark_idba_output_$cores.txt 2>&1

	# Stop nmon capturing
        kill -USR2 $NMON_PID

	# Make nmon html graphs
        sh nmonchart/nmonchart nmon_stats/$NMON_FILE_NAME.nmon nmon_stats/$NMON_FILE_NAME.html

	echo "" >> results/benchmark_idba_time_$cores.txt
else
	echo "IDBA will not be started as you did not choose the genomics tools, all tools or the tool itself."
fi


if [ $toolgroup == "all" ] || [ $toolgroup == "genomics" ] || [ $toolgroup == "clustalomega" ]
then
	# ClustalOmega
        echo "Running ClustalOmega benchmark on dataset $dataset_name_clustalOmega"
        rm -rf benchmark_output/clustalOmega/*
        echo "Replica_$replica clustalOmega with $cores cores on dataset $dataset_name_clustalOmega" >> results/benchmark_clustalomega_time_$cores.txt
        date >> results/benchmark_clustalomega_time_$cores.txt

        # Start nmon capturing
        NMON_FILE_NAME="ClustalOmega_"$replica"_$(date +"%Y-%m-%d-%H-%M")"
        NMON_PID=$(nmon -F $NMON_FILE_NAME.nmon -m nmon_stats/ -p -s 2 -c 12000000)
	
	/usr/bin/time -p -a -o results/benchmark_clustalomega_time_$cores.txt sh -c "clustalOmega/clustal-omega-1.2.4/bin/clustalo -i $dataset_clustalOmega -o benchmark_output/clustalOmega/$dataset_name_clustalOmega.fa --force --outfmt=fa --threads=$cores" >> results/benchmark_clustalomega_output_$cores.txt 2>&1

	# Stop nmon capturing
        kill -USR2 $NMON_PID

        # Make nmon html graphs
        sh nmonchart/nmonchart nmon_stats/$NMON_FILE_NAME.nmon nmon_stats/$NMON_FILE_NAME.html

        echo "" >> results/benchmark_clustalomega_time_$cores.txt
else
        echo "ClustalOmega will not be started as you did not choose the genomics tools, all tools or the tool itself."
fi



if [ $toolgroup == "all" ] || [ $toolgroup == "ml" ] || [ $toolgroup == "tensorflow" ]
then
	# Tensorflow
	echo "Running Tensorflow benchmark with $tf_steps steps"
	rm -rf benchmark_output/tensorflow/*
	echo "Replica_$replica Tensorflow with $cores cores on dataset cifar10 with $tf_steps" >> results/benchmark_tensorflow_time_$cores.txt
	date >> results/benchmark_tensorflow_time_$cores.txt


	# Start nmon capturing
	NMON_FILE_NAME="tensorflow_"$replica"_$(date +"%Y-%m-%d-%H-%M")"
        NMON_PID=$(nmon -F $NMON_FILE_NAME.nmon -m nmon_stats/ -p -s 2 -c 12000000)

	/usr/bin/time -p -a -o results/benchmark_tensorflow_time_$cores.txt sh -c "python datasets/tensorflow/models/tutorials/image/cifar10/cifar10_train.py --data_dir=datasets/tensorflow/ --train_dir=benchmark_output/tensorflow/cifar10_train --max_steps=$tf_steps --threads=$cores" >> results/benchmark_tensorflow_output_$cores.txt 2>&1
	
	# Stop nmon capturing
        kill -USR2 $NMON_PID

	# Make nmon html graphs
        sh nmonchart/nmonchart nmon_stats/$NMON_FILE_NAME.nmon nmon_stats/$NMON_FILE_NAME.html
	
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


	# Start nmon capturing
	NMON_FILE_NAME="GROMACS_"$replica"_$(date +"%Y-%m-%d-%H-%M")"
        NMON_PID=$(nmon -F $NMON_FILE_NAME.nmon -m nmon_stats/ -p -s 2 -c 12000000)

	/usr/bin/time -p -a -o results/benchmark_gromacs_time_$cores.txt sh -c "/usr/local/gromacs/bin/gmx mdrun -v -pin on -nt $cores -s datasets/gromacs/adh_cubic/topol.tpr -o benchmark_output/gromacs/benchmark -cpo benchmark_output/gromacs/benchmark -e benchmark_output/gromacs/benchmark -g benchmark_output/gromacs/benchmark -c benchmark_output/gromacs/benchmark" >> results/benchmark_gromacs_output_$cores.txt 2>&1

	# Stop nmon capturing
        kill -USR2 $NMON_PID

	# Make nmon html graphs
        sh nmonchart/nmonchart nmon_stats/$NMON_FILE_NAME.nmon nmon_stats/$NMON_FILE_NAME.html
	
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

	# Start nmon capturing
	NMON_FILE_NAME="SPAdes_"$replica"_$(date +"%Y-%m-%d-%H-%M")"
        NMON_PID=$(nmon -F $NMON_FILE_NAME.nmon -m nmon_stats/ -p -s 2 -c 12000000)

	/usr/bin/time -p -a -o results/benchmark_SPAdes_time_$cores.txt sh -c "python SPAdes/SPAdes-3.12.0-Linux/bin/spades.py -s $dataset -o benchmark_output/SPAdes/ -t $cores" >> results/benchmark_SPAdes_output_$cores.txt 2>&1

	# Stop nmon capturing
        kill -USR2 $NMON_PID

	# Make nmon html graphs
        sh nmonchart/nmonchart nmon_stats/$NMON_FILE_NAME.nmon nmon_stats/$NMON_FILE_NAME.html

	echo "" >> results/benchmark_SPAdes_time_$cores.txt
else
	echo "SPAdes will not be started as you did not choose the genomics tools, all tools or the tool itself."
fi

if [ $toolgroup == "all" ] || [ $toolgroup == "genomics" ] || [ $toolgroup == "mafft" ]
then
        # MAFFT
        rm -rf benchmark_output/MAFFT/*               #Clean up MAFFT output directoy
        echo "Running MAFFT benchmark on dataset $dataset_name_clustalOmega" 
        echo "Replica_$replica MAFFT with $cores cores on dataset $dataset_name_clustalOmega" >> results/benchmark_mafft_time_$cores.txt                             #Create results file with walltime
        date >> results/benchmark_mafft_time_$cores.txt  #Add date to walltime file

        # Start nmon capturing
        NMON_FILE_NAME="mafft_"$replica"_$(date +"%Y-%m-%d-%H-%M")"
        NMON_PID=$(nmon -F $NMON_FILE_NAME.nmon -m nmon_stats/ -p -s 2 -c 12000000) # -s is the interval between snapshots, -c is the number of snapshots very high to ensure them to the end of the tool run

        # Run MAFFT on dataset $dataset_name_clustalOmega
        /usr/bin/time -p -a -o results/benchmark_mafft_time_$cores.txt sh -c "MAFFT/mafft-7.475-with-extensions/bin/mafft --auto --thread $cores $dataset_clustalOmega > benchmark_output/MAFFT/benchmark" >> results/benchmark_MAFFT_output_$cores.txt 2>&1

        # Stop nmon capturing
        kill -USR2 $NMON_PID

        # Make nmon html graphs
        sh nmonchart/nmonchart nmon_stats/$NMON_FILE_NAME.nmon nmon_stats/$NMON_FILE_NAME.html
        echo "" >> results/benchmark_mafft_time_$cores.txt       #Blank line for clarity and parsing
else
        echo "MAFFT will not be started as you did not choose the genomics tools, all tools or the tool itself."
fi

if [ $toolgroup == "all" ] || [ $toolgroup == "genomics" ] || [ $toolgroup == "sina" ]
then
        # SINA
        rm -rf benchmark_output/SINA/*               #Clean up SINA output directoy
        echo "Running SINA benchmark on reference datatset $reference_name_SINA and dataset $dataset_name_SINA" 
        echo "Replica_$replica SINA with $cores cores on reference dataset $reference_name_SINA and dataset $dataset_name_SINA" >> results/benchmark_sina_time_$cores.txt                             #Create results file with walltime
        date >> results/benchmark_sina_time_$cores.txt  #Add date to walltime file

        # Start nmon capturing
        NMON_FILE_NAME="sina_"$replica"_$(date +"%Y-%m-%d-%H-%M")"
        NMON_PID=$(nmon -F $NMON_FILE_NAME.nmon -m nmon_stats/ -p -s 2 -c 12000000) # -s is the interval between snapshots, -c is the number of snapshots very high to ensure them to the end of the tool run

        # Run SINA on reference dataset $reference_name_SINA and datatset $dataset_name_SINA
        /usr/bin/time -p -a -o results/benchmark_sina_time_$cores.txt sh -c "SINA/sina-1.7.2-linux/sina --threads $cores -i $dataset_SINA -r $reference_SINA -o benchmark_output/SINA/benchmark" >> results/benchmark_sina_output_$cores.txt 2>&1

        # Stop nmon capturing
        kill -USR2 $NMON_PID

        # Make nmon html graphs
        sh nmonchart/nmonchart nmon_stats/$NMON_FILE_NAME.nmon nmon_stats/$NMON_FILE_NAME.html
        echo "" >> results/benchmark_sina_time_$cores.txt       #Blank line for clarity and parsing
else
        echo "SINA will not be started as you did not choose the genomics tools, all tools or the tool itself."
fi

touch benchmark_summary_$cores.txt

path="results/benchmark_bowtie_build_time_$cores.txt"
check_results $path $cores

path="results/benchmark_bowtie_align_time_$cores.txt"
check_results $path $cores

path="results/benchmark_bbmap_time_$cores.txt"
check_results $path $cores

path="results/benchmark_bwa_mem_time_$cores.txt"
check_results $path $cores

path="results/benchmark_velvet_time_$cores.txt"
check_results $path $cores

path="results/benchmark_idba_time_$cores.txt"
check_results $path $cores

path="results/benchmark_clustalomega_time_$cores.txt"
check_results $path $cores

path="results/benchmark_mafft_time_$cores.txt"
check_results $path $cores

path="results/benchmark_sina_time_$cores.txt"
check_results $path $cores

path="results/benchmark_tensorflow_time_$cores.txt"
check_results $path $cores

path="results/benchmark_gromacs_time_$cores.txt"
check_results $path $cores

path="results/benchmark_SPAdes_time_$cores.txt"
check_results $path $cores

fi
}


if [ "$clean" == 1 ]
then	
	for filepath in benchmark_summary_*.txt; do

		if [ -e "$filepath" ]
		then
			#backup_date=$(stat -c %y benchmark_summary_*.txt | cut -d ' ' -f1)
			#backup_time=$(stat -c %y benchmark_summary_*.txt | cut -d ' ' -f2 | cut -d. -f1)
			#backup_dir_name="$backup_date-$backup_time"
			backup_date=$(stat -c %y $filepath | cut -d ' ' -f1)
                        backup_time=$(stat -c %y $filepath | cut -d ' ' -f2 | cut -d. -f1)
                        backup_dir_name="$backup_date-$backup_time"
			mkdir backed_up_benchmark_results/$backup_dir_name
			tar -cf backed_up_benchmark_results/$backup_dir_name/results.tar results/*
			cp benchmark_summary_* backed_up_benchmark_results/$backup_dir_name
			cp bootable_system_info.txt backed_up_benchmark_results/$backup_dir_name
			cp scaling_plot* backed_up_benchmark_results/$backup_dir_name
			tar -cf backed_up_benchmark_results/$backup_dir_name/nmon_stats.tar nmon_stats/*
			rm benchmark_summary_*
			rm -rf results/*
			rm bootable_system_info.txt
			rm scaling_plot*
			rm -rf nmon_stats/*
			break
		else
			while true; do
                	read -p "No full benchmark could have been detected so there are no files to back them up, all related data from stopped run before will be deleted. Do you want to start the benchmark?" yn
                	case $yn in
                        	[Yy]* )	rm -rf results/*
                			rm -f bootable_system_info.txt
                			rm -rf nmon_stats/*
					break
					;;
                        	[Nn]* ) echo "The benchmarks will not be started"
					exit 1
                                	;;
                        	* ) echo "Please answer yes or no."
                                	;;
                	esac
        		done
			break
		fi
	done
fi

touch bootable_system_info.txt
date=$(date)
if [ $clean == 1 ]
then
        used_command="-c -d $dataset -p $default_cores -r $default_replicas -t $default_toolgroup"
else
	used_command="-d $dataset -p $default_cores -r $default_replicas -t $default_toolgroup -s $scaling"
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
	echo "tuned: NOT installed." >> bootable_system_info.txt
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
bowtie2/bowtie2-2.4.2/bowtie2 --version >> bootable_system_info.txt
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

if [[ $dataset == "large" || "$scaling" == "large" ]]
then
	dataset="datasets/1000_genomes/ERR251006.filt.fastq"
	dataset_idba="datasets/1000_genomes/ERR251006.filt.fa"
	dataset_clustalOmega="datasets/clustalOmega/wgs.ANCA.1_500.fsa"
        default_reference="datasets/1000_genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
	default_SINA="datasets/SINA/GTDB_bac-arc_ssu_r86.fa"
	reference_SINA="datasets/SINA/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb"
	default_reference_BWA="datasets/BWA/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla"
	default_tensorflow_steps=5000
	default_gromacs_steps=50000

elif [[ $dataset == "medium" || "$scaling" == "medium" ]]
then
	dataset="datasets/1000_genomes/ERR016155.filt.fastq"
        dataset_idba="datasets/1000_genomes/ERR015528.filt.fa"
	dataset_clustalOmega="datasets/clustalOmega/wgs.ANCA.1_400.fsa"
	default_reference="datasets/ebi/DRR001025.fa"
	default_SINA="datasets/SINA/RefSeq-RDP16S_v2_May2018.fa"
	reference_SINA="datasets/SINA/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb"
	default_reference_BWA="datasets/BWA/DRR001025/DRR001025"
	default_tensorflow_steps=2500
	default_gromacs_steps=30000

elif [[ $dataset == "small" || "$scaling" == "small" ]]
then	
	dataset="datasets/1000_genomes/ERR016155.filt.fastq"
	dataset_idba="datasets/1000_genomes/SRR741411.filt.fa"
	dataset_clustalOmega="datasets/clustalOmega/wgs.ANCA.1_200.fsa"
	default_reference="datasets/ebi/DRR001012.fa"
	default_SINA="datasets/SINA/OE-38_R1.fa"
        reference_SINA="datasets/SINA/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb"
        default_reference_BWA="datasets/BWA/DRR001012/DRR001012"
	default_tensorflow_steps=1000
	default_gromacs_steps=10000

else
	dataset="datasets/1000_genomes/ERR016155.filt.fastq"
        dataset_idba="datasets/1000_genomes/ERR015528.filt.fa"
	dataset_clustalOmega="datasets/clustalOmega/wgs.ANCA.1_400.fsa"
        default_reference="datasets/ebi/DRR001025.fa"
	default_SINA="datasets/SINA/RefSeq-RDP16S_v2_May2018.fa"
        reference_SINA="datasets/SINA/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb"
        default_reference_BWA="datasets/BWA/DRR001025/DRR001025"
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
	

if [[ $default_toolgroup != "all" && $default_toolgroup != "genomics" && $default_toolgroup != "ml" && $default_toolgroup != "quant" && $default_toolgroup != "bowtie2-build" && $default_toolgroup != "velvet" && $default_toolgroup != "idba" && $default_toolgroup != "tensorflow" && $default_toolgroup != "gromacs" && $default_toolgroup != "SPAdes" && $default_toolgroup != "clustalomega" && $default_toolgroup != "bbmap" && $default_toolgroup != "bwa" && $default_toolgroup != "mafft" && $default_toolgroup != "sina" ]]
then
	echo "Parameter is not one of all, genomics, ml, quant, bowtie2-build, velvet, idba, tensorflow, gromacs, SPAdes, clustalomega, bbmap, bwa, mafft or sina. Please check -t flag again."
	exit 1
fi

if [[ "$scaling" == "none" && "$own_tool_path" == "none" ]]
then
	echo "General genomic dataset for Bowtie2, Velvet, SPAdes, BBMap BWA: $dataset"
	echo "IDBA dataset: $dataset_idba"
	echo "Reference dataset for Bowtie2, BBMap, BWA: $default_reference"
	echo "ClustalOmega and MAFFT dataset: $dataset_clustalOmega"
	echo "SINA reference dataset: $reference_SINA"
	echo "SINA dataset: $default_SINA"
	echo "Number of used cores: $default_cores"
	echo "Number of used replicates: $default_replicas"
	echo "Number of Tensorflow steps: $default_tensorflow_steps"
	echo "Number of GROMACS steps: $default_gromacs_steps"
	echo "Toolgroup: $default_toolgroup"

	echo "BOOTABLE benchmark run with $default_cores cores and $default_replicas replicates"

	for replica in $( seq 1 $default_replicas ) 
	do
		run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
	done

elif [[ "$scaling" != "none" && "$own_tool_path" == "none" ]]
then
	echo "Scaling mode is chosen and scaling benchmark will be conducted."
	echo "BOOTABLE scaling benchmark run with 1 core 1/4 of available cores, 1/2 of available cores, all available cores and $default_replicas replicates each"
	
	for replica in $( seq 1 $default_replicas )
        do
		default_cores=$one_core
		default_cores_velvet=$one_core
		echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=2
                default_cores_velvet=$one_core
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

	for replica in $( seq 1 $default_replicas )
        do
                default_cores=3
                default_cores_velvet=2
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

	for replica in $( seq 1 $default_replicas )
        do
                default_cores=4
                default_cores_velvet=3
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=5
                default_cores_velvet=4
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=6
                default_cores_velvet=5
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=7
                default_cores_velvet=6
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=8
                default_cores_velvet=7
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=9
                default_cores_velvet=8
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=10
                default_cores_velvet=9
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=11
                default_cores_velvet=10
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=12
                default_cores_velvet=11
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=13
                default_cores_velvet=12
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=14
                default_cores_velvet=13
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=15
                default_cores_velvet=14
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=16
                default_cores_velvet=15
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=17
                default_cores_velvet=16
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=18
                default_cores_velvet=17
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=19
                default_cores_velvet=18
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=20
                default_cores_velvet=19
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=21
                default_cores_velvet=20
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=22
                default_cores_velvet=21
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=23
                default_cores_velvet=22
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=24
                default_cores_velvet=23
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=25
                default_cores_velvet=24
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=26
                default_cores_velvet=25
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=27
                default_cores_velvet=26
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=28
                default_cores_velvet=27
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=29
                default_cores_velvet=28
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=30
                default_cores_velvet=29
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=31
                default_cores_velvet=30
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=32
                default_cores_velvet=31
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=33
                default_cores_velvet=32
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=34
                default_cores_velvet=33
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=35
                default_cores_velvet=34
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=36
                default_cores_velvet=35
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

elif [[ "$scaling" == "none" && "$own_tool_path" != "none" ]]
then
	echo "The own tool option is chosen with the toolfile under $own_tool_path"
	
	for replica in $( seq 1 $default_replicas )
        do
		run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
	done

else [[ "$scaling" != "none" && "$own_tool_path" != "none" ]]
	echo "The scaling option is combined with the own tool option and the toolfile under $own_tool_path"
	
	echo "BOOTABLE scaling benchmark run with 1 core 1/4 of available cores, 1/2 of available cores, all available cores and $default_replicas replicates each with your own specified tool"

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=$one_core
                default_cores_velvet=$one_core
                echo "$replica replica of $default_replicas replica with $default_cores core is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=$quarter_cores
                default_cores_velvet=$quarter_cores_velvet
                echo "$replica replica of $default_replicas replica with $default_cores cores is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=$half_cores
                default_cores_velvet=$half_cores_velvet
                echo "$replica replica of $default_replicas replica with $default_cores cores is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done

        for replica in $( seq 1 $default_replicas )
        do
                default_cores=$max_cores
                default_cores_velvet=$max_cores_velvet
                echo "$replica replica of $default_replicas replica with $default_cores cores is running."
                run_benchmark_tools $default_cores $default_cores_velvet $replica $dataset $default_tensorflow_steps $default_gromacs_steps $default_reference $dataset_idba $default_toolgroup $dataset_clustalOmega $default_reference_BWA $default_SINA $reference_SINA $own_tool_path
        done
fi

