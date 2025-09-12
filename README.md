# Magnetstein data
Data for Magnetstein [paper](https://pubs.acs.org/doi/full/10.1021/acs.analchem.3c03594) &amp; [web application](https://bioputer.mimuw.edu.pl/magnetstein). 

For more info about Magnetstein, go [here](https://github.com/BDomzal/magnetstein).

# Contents of this repository

The most important file here is `estimation.ipynb`. This notebook contains the code for estimation of proportions for all the experiments (1-11).

Folders starting from the word `experiment` contain raw data, notebooks used for preprocessing and preprocessed data.

Notebook `spectra_visualisation.ipynb` contains code for creating visualisations (included in folder `visualisations/`).

Notebook `AvRE_calculation.ipynb` contains code for computing averaged relative error.

Notebook `default_parameters_search.ipynb` checks which values of parameters give the least overall error among all the experiments.

In folder `heatmaps/` there are visualisations showing influence of $\kappa_{mixture}$ and $\kappa_{components}$ parameters on the quality of estimation.

# Citing 

If you use data from this repository, please cite:

Domżał, B., Nawrocka, E.K., Gołowicz, D., Ciach, M.A., Miasojedow, B., Kazimierczuk, K., & Gambin, A. (2023). Magnetstein: An Open-Source Tool for Quantitative NMR Mixture Analysis Robust to Low Resolution, Distorted Lineshapes, and Peak Shifts. _Analytical Chemistry_. DOI: [10.1021/acs.analchem.3c03594](https://doi.org/10.1021/acs.analchem.3c03594).


