# Veszee
Visualization and Analyzing tools for a bayesian forward modelling technique called `eszee`.

## Installing Veszee
Here are some basic instructions for installing `Veszee` to get analytical plots and tools from the uv-space modelling tool `eszee`. <br>
The installation follows a `conda`-based procedure, as it will get everything installed and cross-linked correctly for linux machines. If installing in a mac, it might need some tuning (removing packages) that can't be found.

### conda-powered installation
The simplest way for getting `eszee` ready is to install it using `conda` <br>
Once you have downloaded the repository (by `git`-cloning it on your machine or downloading the zip), you should `cd` in the `Veszee` repository and run:

    conda env create -f environment.yml

This will build a `conda` environment named `Veszee_env` with all the packages required to run `Veszee`. In case you want to change the name of the virtual environment, it will be enough to add `--name yourname` to the command above.<br>
To active this, you can use `conda activate Veszee_env` (or any name you used for the environment).

#### Installing Casatools and Casatasks
Only the installment of casatools and casataks could not be installed automatically. Please install casatools and casataks as follows. Enter the newly created virtual environment with `conda activate Veszee_env` and run:

    pip install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatools==6.0.0.27
    pip install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatasks==6.0.0.27
    
For further information on the installment of casatools and taks, we refer to: https://casa.nrao.edu/casadocs/casa-5.6.0/introduction/casa6-installation-and-usage
