# tumoroscope

[![PyPI](https://img.shields.io/pypi/v/mgcpy.svg)](https://pypi.org/project/tumoroscope/)
[![DOI](https://zenodo.org/badge/147731955.svg)](https://doi.org/10.1101/2022.09.22.508914)
[![License](https://img.shields.io/badge/License-GNU%20GPL-blue)](https://opensource.org/licenses/gpl-license)


`tumoroscope` is a Python package for inferring the clonal composition of the cancer tissue across the space using spatial transcriptomics, whole exome sequencing and H&E images.  

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Setting up the development environment](#setting-up-the-development-environment)
- [License](#license)
- [Issues](https://github.com/neurodata/mgcpy/issues)

# Overview
``tumoroscope``  The package utilizes a simple class structure to enhance usability while also allowing easy extension of the package for developers. The package can be installed on all major platforms (e.g. BSD, GNU/Linux, OS X, Windows)from Python Package Index (PyPI) and GitHub.

# System Requirements
## Hardware requirements
`mgcpy` package requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS Requirements
This package is supported for *macOS* and *Linux*. The package has been tested on the following systems:
+ macOS: Mojave (10.14.1)
+ Linux: Ubuntu 16.04

### Python Dependencies
`mgcpy` mainly depends on the Python scientific stack.

```
numpy
scipy
Cython
scikit-learn
pandas
seaborn
```

# Installation Guide:

### Install from PyPi
```
pip3 install mgcpy
```

### Install from Github
```
git clone https://github.com/neurodata/mgcpy
cd mgcpy
python3 setup.py install
```
- `sudo`, if required
- `python3 setup.py build_ext --inplace  # for cython`, if you want to test in-place, first execute this

# Setting up the development environment:
- To build image and run from scratch:
  - Install [docker](https://docs.docker.com/install/)
  - Build the docker image, `docker build -t mgcpy:latest .`
    - This takes 10-15 mins to build
  - Launch the container to go into mgcpy's dev env, `docker run -it --rm --name mgcpy-env mgcpy:latest`
- Pull image from Dockerhub and run:
  - `docker pull tpsatish95/mgcpy:latest` or `docker pull tpsatish95/mgcpy:development`
  - `docker run -it --rm -p 8888:8888 --name mgcpy-env tpsatish95/mgcpy:latest` or `docker run -it --rm -p 8888:8888 --name mgcpy-env tpsatish95/mgcpy:development`


- To run demo notebooks (from within Docker):
  - `cd demos`
  - `jupyter notebook --ip 0.0.0.0 --no-browser --allow-root`
  - Then copy the url it generates, it looks something like this: `http://(0de284ecf0cd or 127.0.0.1):8888/?token=e5a2541812d85e20026b1d04983dc8380055f2d16c28a6ad`
  - Edit this: `(0de284ecf0cd or 127.0.0.1)` to: `127.0.0.1`, in the above link and open it in your browser
  - Then open `mgc.ipynb`

- To mount/load local files into docker container:
  - Do `docker run -it --rm -v <local_dir_path>:/root/workspace/ -p 8888:8888 --name mgcpy-env tpsatish95/mgcpy:latest`, replace `<local_dir_path>` with your local dir path.
  - Do `cd ../workspace` when you are inside the container to view the mounted files. The **mgcpy** package code will be in `/root/code` directory.


## MGC Algorithm's Flow
![MGCPY Flow](https://raw.githubusercontent.com/neurodata/mgcpy/master/MGCPY.png)

## Power Curves
- Recreated Figure 2 in https://arxiv.org/abs/1609.05148, with the addition of MDMR and Fast MGC
![Power Curves](https://raw.githubusercontent.com/neurodata/mgcpy/master/power_curves_dimensions.png)

# License

This project is covered under the **Apache 2.0 License**.




Installation:
```
 pip install tumoroscope
```
Usage:
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


