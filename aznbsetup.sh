# Setup environment for Azure Notebooks
# Steps:
#    - Download, install and setup the latest version of Miniconda
#    - Create a new conda environment 'snakemake-workshop' with everything required for the workshop

# Path to miniconda installation
CONDA_PATH=~/miniconda

# Download the last miniconda installer
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install conda in silent mode
sh Miniconda3-latest-Linux-x86_64.sh -b -p "$CONDA_PATH"

# Setup miniconda for current shell and activate base environment
eval "$("$CONDA_PATH"/bin/conda shell.bash hook)"

# Setup conda in shell permanently
conda init

# Update conda
conda update conda -y

# Create workshop environment
conda env create -f library/workshop.yaml

