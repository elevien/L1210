
# **Single cell growth dynamics and size control**

## Organization of this repository

The ``src`` folder contains all the code. I've seperated code for processing/simulation and analysis. All the analysis code, which produces the paper figures is in the ``notebooks`` folder. 


This is where all the project specific code is contained. Within this folder you will find the following:
  
The core functionality is in the following files:
* ``gp_pipeline.jl``  Gaussian process pipeline on the experimental data
* ``generate_sums.jl`` Generates the simulated data used for our analysis


Additional code and output is in the folder:
* ``figures`` All figures output here
* ``output`` All the processed data, such as simulated datasets and output from gp pipeline ends up here
* ``notebooks`` Jupyter Notebooks including those which make the paper figures 
* ``gp_models`` Models for the gp pipeline to use (currently only matern32 is used in the pipeline)

###  ``experimental_data``

The raw experimental data as well as the processed files are stored here. In the current working version of the repository these files are ignore by git. 


## How to run the data analysis

1 **setup julia** You need to install julia, then navigate to this repository in terminal, run ``julia``, then type ``] activate .``. This will setup the environment and ensure all the packages are available. 

2 **Generate and process data** Navigate to src folder and run the commands ``julia  --project=./../ gp_pipeline.jl`` and ``julia  --project=./../ sims_pipeline.jl``. 

3 **Run the analysis notebooks** The notebooks which are used to make the paper figures all begin with ``PAPER`` and can be run in any order. 

