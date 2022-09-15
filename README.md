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


