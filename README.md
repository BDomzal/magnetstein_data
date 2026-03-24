# <img src="https://github.com/user-attachments/assets/008fa262-7fa1-458c-a8fc-f4f6e3a2022a" align="left" height="80" width="60"/> Magnetstein data

<img width="2700" height="1200" alt="all_components_variant_4" src="https://github.com/user-attachments/assets/0fe12e68-996b-4c93-bbd8-40ae922ee87c" />

Data for Magnetstein [paper](https://pubs.acs.org/doi/full/10.1021/acs.analchem.3c03594) &amp; [web application](https://bioputer.mimuw.edu.pl/magnetstein). 

For more info about Magnetstein, go [here](https://github.com/BDomzal/magnetstein).

# Contents of this repository

The most important file here is `estimation.py`. This script contains the code for estimation of proportions for all the experiments (1-11).

Folders starting from the word `experiment` contain raw data, scripts used for preprocessing and preprocessed data.

Script `spectra_visualisation.py` contains code for creating visualisations (included in folder `visualisations/`).

Script `AvRE_calculation.py` contains code for computing averaged relative error.

Script `default_parameters_search.py` checks which values of parameters give the least overall error among all the experiments.

In folder `heatmaps/` there are visualisations showing influence of $\kappa_{mixture}$ and $\kappa_{components}$ parameters on the quality of estimation.

# Citing 

If you use data from this repository, please cite:

Domżał, B., Nawrocka, E.K., Gołowicz, D., Ciach, M.A., Miasojedow, B., Kazimierczuk, K., & Gambin, A. (2023). Magnetstein: An Open-Source Tool for Quantitative NMR Mixture Analysis Robust to Low Resolution, Distorted Lineshapes, and Peak Shifts. _Analytical Chemistry_. DOI: [10.1021/acs.analchem.3c03594](https://doi.org/10.1021/acs.analchem.3c03594).


