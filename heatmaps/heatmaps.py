
nr_of_experiment = 10
variant = 3


import json

with open('../config.json', 'r') as f:
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
EXPERIMENTS_FOLDERS = {key: '../' + value for key, value in EXPERIMENTS_FOLDERS.items()}
GROUND_TRUTH_MOLAR_PROPORTIONS = config["ground_truth_molar_proportions"]
if 'variant' in locals():
    GROUND_TRUTH_MOLAR_PROPORTIONS["experiment_10"] = GROUND_TRUTH_MOLAR_PROPORTIONS["experiment_10"][variant]



all_components_results_molar = get_results_concentrations(
                                                            nr_of_experiment, 
                                                            EXPERIMENTS_FOLDERS, 
                                                            GROUND_TRUTH_MOLAR_PROPORTIONS, 
                                                            PROTONS_DICTIONARY,
                                                            variant=variant
                                                            )


draw_heatmap_power_norm(
							nr_of_experiment,
							all_components_results_molar,
                            variant=variant
							)

