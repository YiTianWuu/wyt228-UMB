# Install Miniconda3 on Linux
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
source miniconda3/bin/activate
conda init
exit

# Create environment for python 2
conda create --name py2env
# Activate environment
conda activate py2env
# Install python 2
conda install python=2.7.16
# Install packages needed
conda install -c conda-forge biopython
conda config --add channels r
conda config --add channels bioconda
conda install pysam
# Install bowtie2 for alignment
conda install bowtie2


# You may exit the python 2 environment
conda deactivate