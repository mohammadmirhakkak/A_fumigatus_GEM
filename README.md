# A_fumigatus_GEM
This repository contains the data and scripts (Python and R) that were used in generation of a pan-genome and 252 *A. fumigatus* strain genome-scale metabolic models (GEMs) and their interaction simulations with the lung microbiome described by Mirhakkak et al. (2023).
Running the scripts in A_fumigatus_GEM/src/ will generate the results in A_fumigatus_GEM/res/.

For instance, one can use the following Terminal command to regenerate the strain GEMs.

`$ python A_fumigatus_GEM/src/build_strain_gem.py`

Or the following in interactive python (e.g. ipython3)

`%run A_fumigatus_GEM/src/build_strain_gem.py`

**Python Dependencies**

* COBRApy 0.17.1

* IBM CPLEX 12.8.0.0



**Python Modules:**

* glob

* pickle

* numpy

* copy

* pandas

* re

* timeit

* itertools


**R packages**

* ggplot2 3.3.6

* ggrepel 0.9.1

* ComplexHeatmap 2.10.0

* circlize 0.4.15

* RColorBrewer 1.1-3

* ggtree 3.2.1

* ape 5.6-2

* tidyverse 1.3.2

* ggstance 0.3.5

* viridis 0.6.2

* dplyr 1.0.9

* ggsignif 0.6.3

* vegan 2.6-2

* tibble 3.1.8

* stringr 1.4.0

* phyloseq 1.34.0

* ggpubr 0.4.0


**Citation**

Mirhakkak, M. H., Chen, X., Ni, Y., Heinekamp, T., Sae-Ong, T., Xu, L. L., ..., Schäuble, S. & Panagiotou, G. (2023). Genome-scale metabolic modeling of Aspergillus fumigatus strains reveals growth dependencies on the lung microbiome. Nature Communications, 14(1), 4369.
