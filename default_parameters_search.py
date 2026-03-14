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
GROUND_TRUTH_MOLAR_PROPORTIONS = config["ground_truth_molar_proportions"]


list_of_results_for_experiments = get_all_components_results_concentrations(
                                                                             EXPERIMENTS_FOLDERS, 
                                                                             GROUND_TRUTH_MOLAR_PROPORTIONS, 
                                                                             PROTONS_DICTIONARY
                                                                            )


total_error_all_experiments = get_total_error_from_all_experiments(list_of_results_for_experiments)


draw_heatmap(total_error_all_experiments)

a, b = total_error_all_experiments.stack().idxmin()
print('Optimal kappa_mixture:')
print(a)
print('Optimal kappa_components:')
print(b)
