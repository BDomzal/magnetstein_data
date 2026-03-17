
nr_of_experiment = 3
experiment_name = 'experiment_' + str(nr_of_experiment)

variant = 3

default_MTD = 0.25
default_MTD_th = 0.22


import json

with open('config.json', 'r') as f:
    config = json.load(f)

MAGNETSTEIN_PATH = config["magnetstein_path"]
SRC_PATH = config["src_path"]

import sys
sys.path.insert(0, MAGNETSTEIN_PATH)
sys.path.insert(1, SRC_PATH)
from utils import *

COMPONENTS_DICTIONARY = config["components_dictionary"]
PROTONS_DICTIONARY = config["protons_dictionary"]
GROUND_TRUTH_MOLAR_PROPORTIONS = config["ground_truth_molar_proportions"]
if 'variant' in locals():
    GROUND_TRUTH_MOLAR_PROPORTIONS["experiment_10"] = GROUND_TRUTH_MOLAR_PROPORTIONS["experiment_10"][variant]

EXPERIMENTS_FOLDERS = config["experiments_folders"]

XLIMS_LOWER = config["xlims_lower"]
XLIMS_UPPER = config["xlims_upper"] 
YLIMS_LOWER = config["ylims_lower"]
YLIMS_UPPER = config["ylims_upper"]

XLIMS_LOWER_COMP = config["xlims_lower_comp"]
XLIMS_UPPER_COMP = config["xlims_upper_comp"] 
YLIMS_LOWER_COMP = config["ylims_lower_comp"]
YLIMS_UPPER_COMP = config["ylims_upper_comp"]


# ### Loading the data


mix = get_mix_spectrum(nr_of_experiment, EXPERIMENTS_FOLDERS, variant=variant)



spectra, names, how_many_components = get_components_spectra(nr_of_experiment, 
                                                             EXPERIMENTS_FOLDERS, 
                                                             COMPONENTS_DICTIONARY, 
                                                             PROTONS_DICTIONARY, 
                                                             variant=variant)


# ### Estimation with default kappas


default_MTD = 0.25
default_MTD_th = 0.22



start_vis = time.time()
estimation_vis = estimate_proportions(
                                      mix, 
                                      spectra,
                                      MTD=default_MTD, 
                                      MTD_th=default_MTD_th,
                                      verbose=False, 
                                      solver=pulp.GUROBI(msg=False), 
                                      what_to_compare='area'
                                    )
stop_vis = time.time()



print('Proportions (based on area only):')
print(estimation_vis['proportions'])
print('Proportion of noise in components:')
print(estimation_vis['proportion_of_noise_in_components'])
print('Propotion of noise in mixture:')
print(1 - sum(estimation_vis['proportions']))
print('Estimation took '+ str(stop_vis-start_vis) + ' seconds.')


# ### Components without scaling


plot_components_without_scaling(
									spectra, 
									COMPONENTS_DICTIONARY,
									nr_of_experiment,
                                    XLIMS_LOWER,
                                    XLIMS_UPPER,
									variant=0,
									#path_to_save=os.getcwd()+'/visualisations/',
									colors = ['blue', 'orange', 'green', 'red', 'pink']
									)


# ### Mixture vs linear combination of components (noise removed)

# #### Components added in estimated proportions (noise removed)


plot_components_added_in_estimated_proportions(
                                                    spectra,
                                                    mix,
                                                    estimation_vis,
                                                    nr_of_experiment,
                                                    COMPONENTS_DICTIONARY,
                                                    XLIMS_LOWER_COMP,
                                                    XLIMS_UPPER_COMP,
                                                    YLIMS_LOWER_COMP,
                                                    YLIMS_UPPER_COMP
                                                )


# #### Mixture with removed noise



mix_without_noise = get_mix_without_noise(mix, estimation_vis)
#mix_without_noise = get_intensities(mix_without_noise)



plot_mixture(
                mix_without_noise,
                nr_of_experiment,
                COMPONENTS_DICTIONARY,
                XLIMS_LOWER_COMP,
                XLIMS_UPPER_COMP,
                YLIMS_LOWER_COMP,
                YLIMS_UPPER_COMP,
                path_to_save=None
                )


# #### Components added in estimated proportions (noise removed) + mix without noise


plot_components_added_in_estimated_proportions(
                                                    spectra,
                                                    mix_without_noise,
                                                    estimation_vis,
                                                    nr_of_experiment,
                                                    COMPONENTS_DICTIONARY,
                                                    XLIMS_LOWER_COMP,
                                                    XLIMS_UPPER_COMP,
                                                    YLIMS_LOWER_COMP,
                                                    YLIMS_UPPER_COMP,
                                                    include_mixture=True
                                                )

