# Tumoroscope-package

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


