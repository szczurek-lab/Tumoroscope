# Tumoroscope

[![PyPI](https://img.shields.io/pypi/v/mgcpy.svg)](https://pypi.org/project/tumoroscope/)
[![DOI](https://zenodo.org/badge/147731955.svg)](https://doi.org/10.1101/2022.09.22.508914)
[![License](https://img.shields.io/badge/License-GNU%20GPL-blue)](https://opensource.org/licenses/gpl-license)


`tumoroscope` is a Python package for inferring the clonal composition of the cancer tissue across the space using spatial transcriptomics, whole exome sequencing and H&E images.  

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Setting up the development environment](#setting-up-the-development-environment)
- [License](#license)
- [Usage](#usage)

# Overview
``tumoroscope``  The package utilizes a simple class structure to enhance usability while also allowing easy extension of the package for developers. The package can be installed on all major platforms (e.g. BSD, GNU/Linux, OS X, Windows)from Python Package Index (PyPI) and GitHub.

# System Requirements
## Hardware requirements
`tumoroscope` package requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS Requirements
This package is supported for *macOS* and *Linux*. The package has been tested on the following systems:
+ macOS (ventura 13.0 and Catalina 10.15.7)
+ Linux (ubuntu 14.04)

### Python Dependencies
`tumoroscope` mainly depends on the Python scientific stack.

```
numpy
pandas
matplotlib
seaborn
pickle5
scipy
```

# Setting up the development environment:

- Install conda (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- create a new envirnoment with python >=3.7 and activate it:
```
    conda create -n demo python=3.8
    conda activate demo
```

- Install conda on the jupyter notebook as well 
```
    conda install -c anaconda ipykernel
```
- Run the following command for having your created environment in jupyther notebook
```
    python -m ipykernel install --user --name=demo
```

# Installation Guide:

### Install from PyPi
```
pip3 install tumoroscope
```
Installation time is less than one minute.
# License

This project is covered under the **GNU General Public License**.


# Usage


As a demo, you can use demo01.ipynb in demo folder to generate a simulated data and run tumoroscope on that. In this notebook, at the end, you will plot the true values and the inferred values of the fraction of the clones in the spots in the x and y of the plot. Also, you can see the Mean Average Error (MAE) written above the plot as a result. The running time is 1-3 minutes.

Also, for applying on the real data, you can use the following code.
```
from tumoroscope.run_tumoroscope import run_tumoroscope
run_tumoroscope('configs/config_prostate.json')
``` 
config includes adress to 4 files and one folder. 

1. in 'observed' section, st calliing files 
2. in 'C_variation', ssm_data.txt, output of canopy including mutations
3. in 'C_variation', C_tree_1.txt file, transfered output of canopy including genotype
4. in 'n_variation', cell count and annotation

and the folder: "data_tumoroscope_inputs" : "prostate_tumoroscope_input/".


```
{
    "results": {
          "text_result": "inferred_vars", #name of the file which has the important resutls in txt format
          "data_tumoroscope_inputs" : "prostate_tumoroscope_input/generated/", #point to the folder including the input files
          "output" : "Results_output" #point to the folder of the results
    },
    "structure": {
          "optimal_rate": 0.40, #optimal acceptance rate for te MH sampling
          "pi_2D": true, #considering variable pi as a 2D variable 
          "all_spots": false #considering all spots and not only annotated spots
    },
    "observed": { # read count files regarding to the different sections
      "P1.2": "prostate_tumoroscope_input/vardict2_st_calling/vardict2_1.2.ac",
      "P2.4": "prostate_tumoroscope_input/vardict2_st_calling/vardict2_2.4.ac",
      "P3.3": "prostate_tumoroscope_input/vardict2_st_calling/vardict2_3.3.ac"
    },
    "C_variation": {
        "WES_file" : "prostate_tumoroscope_input/ssm_data.txt", #output of the canopy regarding to the somatic mutations 
        "phyloWGS": "prostate_tumoroscope_input/union_mutations", #regarding to the phyloWGS, not supported yet
        "tree": 3350, #regarding to the phyloWGS, not supported yet
        "selected_canopy_tree": "prostate_tumoroscope_input/C_tree_1.txt", #the genotype matrix which is output of the canopy tree
        "method" : "canopy" #the method that we are using for reconstruction of the tree (only canopy supported) 
    },
    "n_variation":{
        "n_sampling" : true,
        "n_file" : "prostate_tumoroscope_input/prostate_cell_count_annotation.txt" #annotation of the spots for finding out the cancerous spots
    },
   "Z_variation":{
         "avarage_clone_in_spot" : 1.5, #the expected value which we predict for the number of clones in the spots
          "threshold": 0.8 #the probability threshhold for putting 0/1 for the presence of the clones
   },
   "sampling":{
      "onLaptop" : false, # if true, the program would multiply the number of iteratioon and all related numbers by 1/10
      "max_iter" : 800, # max iteration - if mcmc is not conoverged yet, we stop the sampling
      "min_iter" : 500, #min iteration - before this min, we do not stop even if mcmc coonverged
      "burn_in" : 100, 
      "batch" : 100,
      "every_n_sample" : 5, # for example we keep one sample every 5 sample
      "changes_batch" : 10
   },
   "Gamma": {
      "phi_gamma" : [0.1, 0.5], #parameters of the gamma distribution over phi 
      "F_epsilon" : [2, 1], 
      "F_fraction": true,
      "F_file": "prostate_tumoroscope_input/F_tree_1.txt"
   },
     "theta": {
    "gamma" : 1,
    "theta_variable" : false
   },
   "criteria": {
        "offset":1,
        "st_read_limit":0,
        "var_calculation" : true
    }
}

```

# Funding

We acknowledge the support from the National Genomics Infrastructure in Stockholm and Uppsala, funded by Science for Life Laboratory, the Knut and Alice Wallenberg Foundation, the Swedish Research Council, and SNIC/Uppsala Multidisciplinary Center for Advanced Computational Science for assistance with massively parallel sequencing and access to the UPPMAX computational infrastructure. This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 766030 and 844712, the Polish National Science Centre PRELUDIUM grant no. 2021/41/N/ST6/03619, the Polish National Science Centre OPUS grant no. 2019/33/B/NZ2/00956, the Polish National Science Centre SONATA BIS grant no. 2020/38/E/NZ2/00305, The Swedish Cancer Society, The Institut Universitaire de France (AC) and the Swedish Research Council Projekt 2018-06217VR.


