
# **Code for "Stochastic fluctuations in instantaneous growth rates drive cell-to-cell variability in leukemia cells"**

## Organization of this repository

The ``src`` folder contains all of the core functionality. No code in this folder needs to be run directly. The ``pipeline`` folder contains the code to preprocess the data, generate simulations and runs the Gaussian process.   All the analysis code, which produces the paper figures is in the ``notebooks`` folder. The notebooks should be run after running the preprocessing and simulation pipeline as described below. The folder ``smr_dta`` contains the mass measurements from the suspended microchannel resonator (SMR). 

## How to run the data analysis

1 **setup julia** You need to install julia, then navigate to this repository in terminal, run ``julia``, then type ``] activate .``. This will setup the environment and ensure all the packages are available. 

2 **Generate and process data** Run the command ``julia  --project=./../ run.jl``.  

3 **Run the analysis notebooks** The notebooks which are used to make the paper figures all begin with ``PAPER`` and can be run in any order. 

