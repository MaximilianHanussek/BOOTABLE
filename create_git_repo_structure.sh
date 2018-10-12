mkdir benchmark_output
mkdir bowtie2
mkdir datasets
mkdir gromacs
mkdir gromacs/gromacs-2018.3/build
mkdir IDBA
mkdir results
mkdir SPAdes
mkdir tensorflow

cd benchmark_output

mkdir bowtie2
mkdir gromacs
mkdir IDBA
mkdir SPAdes
mkdir tensorflow
mkdir velvet

cd ../datasets

mkdir 1000_genomes
mkdir gromacs
mkdir tensorflow

cd datasets/gromacs/
#wget ftp://ftp.gromacs.org/pub/benchmarks/ADH_bench_systems.tar.gz
wget https://s3.denbi.uni-tuebingen.de/max/ADH_bench_systems.tar.gz
tar -xf ADH_bench_systems.tar.gz
rm ADH_bench_systems.tar.gz
sed -i 's/nsteps.*/nsteps                  = 50000/' adh_cubic/pme_verlet.mdp

cd ../tensorflow/
#wget https://s3.denbi.uni-tuebingen.de/fb-test/tf_benchmark/models.tar.gz
wget https://s3.denbi.uni-tuebingen.de/max/tf_benchmark/models.tar.gz
tar -xf models.tar.gz
rm models.tar.gz



cd ../../bowtie2/
#wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.2/bowtie2-2.3.4.2-source.zip
wget https://s3.denbi.uni-tuebingen.de/max/bowtie2-2.3.4.2-source.zip
unzip bowtie2-2.3.4.2-source.zip
rm bowtie2-2.3.4.2-source.zip

cd ../../
git clone https://github.com/dzerbino/velvet.git velvet/
cd velvet/
rm -rf .gitignore .git/



cd IDBA/
#wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/hku-idba/idba_ud-1.0.9.tar.gz
wget https://s3.denbi.uni-tuebingen.de/max/idba_ud-1.0.9.tar.gz
tar -xf idba_ud-1.0.9.tar.gz
rm idba_ud-1.0.9.tar.gz

#cd ../gromacs/
#wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-2018.3.tar.gz

cd ../SPAdes
#wget http://cab.spbu.ru/files/release3.12.0/SPAdes-3.12.0-Linux.tar.gz
wget https://s3.denbi.uni-tuebingen.de/max/SPAdes-3.12.0-Linux.tar.gz
tar -xf SPAdes-3.12.0-Linux.tar.gz
rm SPAdes-3.12.0-Linux.tar.gz
