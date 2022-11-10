# A_fumigatus_GEM
This repository contains the data and scripts (Python and R) that were used in generation of a pan-genome and 252 *A. fumigatus* strain genome-scale metabolic models (GEMs) and their interaction simulations with the lung microbiome described by Mirhakkak et al. (2023).
Running the scripts in A_fumigatus_GEM/src/ will generate the results in A_fumigatus_GEM/res/.

For instance, one can use the following Terminal command to regenerate the strain GEMs.

`$ python A_fumigatus_GEM/src/build_strain_gem.py`

Or the following in interactive python (e.g. ipython3)

`%run A_fumigatus_GEM/src/build_strain_gem.py`

**Dependencies**

* COBRApy 0.17.1

* IBM CPLEX 12.8.0.0



**Python Built-In Modules:**

* glob

* pickle

* numpy

* copy

* pandas




**Citation**

Mirhakkak, M.H., Chen, X., ... , Schäuble, S., Panagiotou, G. A pan-genome resembling genome-scale metabolic model platform of 252 Aspergillus fumigatus strains reveals growth dependencies on the lung microbiome. Nature Communications (2023).
