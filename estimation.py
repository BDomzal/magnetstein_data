# ### Settings & imports

nr_of_experiment = 5
experiment_name = 'experiment_' + str(nr_of_experiment)

variant = 0

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
EXPERIMENTS_FOLDERS = config["experiments_folders"]
BEST_KAPPA_MIX = config["best_kappa_mix"][nr_of_experiment-1]
BEST_KAPPA_MODEL = config["best_kappa_model"][nr_of_experiment-1]
if nr_of_experiment == 10:
    BEST_KAPPA_MIX = config["best_kappa_mix"][nr_of_experiment-1][variant]
    BEST_KAPPA_MODEL = config["best_kappa_model"][nr_of_experiment-1][variant]


# ### Loading the data


mix = get_mix_spectrum(nr_of_experiment, EXPERIMENTS_FOLDERS, variant=variant)


spectra, names, how_many_components = get_components_spectra(nr_of_experiment, 
                                                             EXPERIMENTS_FOLDERS, 
                                                             COMPONENTS_DICTIONARY, 
                                                             PROTONS_DICTIONARY, 
                                                             variant=variant)


# ### Estimation with default parameters

# Proportions of areas (no rescaling based on the number of NMR-active nuclei):


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
print('Proportion of noise in mixture:')
print(1 - sum(estimation_vis['proportions']))
print('Estimation took '+ str(stop_vis-start_vis) + ' seconds.')


# Proportions of concentrations (rescaled based on the numbers of NMR-active nuclei):


start_con = time.time()
estimation_molar = estimate_proportions(
                                        mix, 
                                        spectra,
                                        MTD=default_MTD, 
                                        MTD_th=default_MTD_th,
                                        verbose=False, 
                                        solver=pulp.GUROBI(msg=False), 
                                        what_to_compare='concentration'
                                        )
stop_con = time.time()


print('Proportions (of concentrations):')
print(estimation_molar['proportions'])
print('Estimation took '+ str(stop_con-start_con) + ' seconds.')


# #MTD = 0.25, MTD_th = 0.22

# #1: 'Pinene', 'Benzyl benzoate'
# 0.0937, 0.9063

# #2: 'Pinene', 'Limonene'
# 0.5240, 0.4760

# #3: 'Isopropyl myristate', 'Benzyl benzoate', 'Alpha pinene', 'Limonene'
# 0.7511, 0.1059, 0.0809, 0.0621

# #4: 'Lactate', 'Alanine', 'Creatine', 'Creatinine', 'Choline chloride'
# 0.3165, 0.2137, 0.1341, 0.1960, 0.1397

# #5: 'Lactate', 'Alanine', 'Creatine', 'Creatinine', 'Choline chloride'
# 0.3180, 0.2113, 0.1293, 0.2018, 0.1397

# #6: 'Pinene', 'Benzyl benzoate'
# 0.3931, 0.6069

# #7: 'Benzyl benzoate', 'm Anisaldehyde'
# 0.8428, 0.1572

# #8: 'Benzyl benzoate', 'm Anisaldehyde'
# 0.3762, 0.6238

# #9: Lactate', 'Alanine', 'Creatine', 'Creatinine', 'Choline chloride'
# 0.3180, 0.2125, 0.1333, 0.1967, 0.1395

# #10: 'Leucine', 'Isoleucine', 'Valine'
# # variant 1
# 0.3552, 0.3418, 0.3029

# #10: 'Leucine', 'Isoleucine', 'Valine'
# # variant 2
# 0.2793, 0.2173, 0.5034

# #10: 'Leucine', 'Isoleucine', 'Valine'
# # variant 3
# 0.2761, 0.4778, 0.2461

# #10: 'Leucine', 'Isoleucine', 'Valine'
# # variant 4
# 0.2494, 0.3357, 0.4149

# #11: 'Leucine', 'Isoleucine', 'Valine'
# 0.4960, 0.2519, 0.2521


# ### Estimation with optimal parameters

# Proportions of areas (no rescaling based on the number of NMR-active nuclei):


start_vis_opt = time.time()
estimation_vis_opt = estimate_proportions(
                                          mix, 
                                          spectra,
                                          MTD=BEST_KAPPA_MIX,
                                          MTD_th=BEST_KAPPA_MODEL,
                                          verbose=False, 
                                          solver=pulp.GUROBI(msg=False),
                                          what_to_compare='area'
                                         )
stop_vis_opt = time.time()


print('Proportions (based on area only, best kappa settings):')
print(estimation_vis_opt['proportions'])
print('Proportion of noise in components:')
print(estimation_vis_opt['proportion_of_noise_in_components'])
print('Proportion of noise in mixture:')
print(1 - sum(estimation_vis_opt['proportions']))
print('Estimation took '+ str(stop_vis_opt-start_vis_opt) + ' seconds.')


# Proportions of concentrations (rescaled based on the numbers of NMR-active nuclei):


start_con_opt = time.time()
estimation_molar_opt = estimate_proportions(mix, 
                                            spectra,
                                            MTD=BEST_KAPPA_MIX, 
                                            MTD_th=BEST_KAPPA_MODEL,
                                            verbose=False, 
                                            solver=pulp.GUROBI(msg=False), 
                                            what_to_compare='concentration')
stop_con_opt = time.time()


print('Proportions (of concentrations):')
print(estimation_molar_opt['proportions'])
print('Estimation took '+ str(stop_con_opt-start_con_opt) + ' seconds.')


# ### Saving estimation results

# Default parameters:

save_proportions(
                    estimation_vis, 
                    nr_of_experiment, 
                    COMPONENTS_DICTIONARY, 
                    EXPERIMENTS_FOLDERS, 
                    filename='/results_area_default_parameters', 
                    variant=variant
                )

save_noise_in_components(
                        estimation_vis, 
                        nr_of_experiment, 
                        EXPERIMENTS_FOLDERS, 
                        filename='/proportion_of_noise_in_comp', 
                        variant=variant
                        )

save_noise_in_mixture(
                        estimation_vis, 
                        nr_of_experiment, 
                        EXPERIMENTS_FOLDERS, 
                        filename='/proportion_of_noise_in_mix', 
                        variant=variant
                        )

save_proportions(
                    estimation_molar, 
                    nr_of_experiment, 
                    COMPONENTS_DICTIONARY, 
                    EXPERIMENTS_FOLDERS, 
                    filename='/results_concentration_default_parameters', 
                    variant=variant
                )


# Optimal parameters:


save_proportions(
                    estimation_vis_opt, 
                    nr_of_experiment, 
                    COMPONENTS_DICTIONARY, 
                    EXPERIMENTS_FOLDERS, 
                    filename='/results_area', 
                    variant=variant
                )

save_proportions(
                    estimation_molar_opt, 
                    nr_of_experiment, 
                    COMPONENTS_DICTIONARY, 
                    EXPERIMENTS_FOLDERS, 
                    filename='/results_concentration', 
                    variant=variant
                )

