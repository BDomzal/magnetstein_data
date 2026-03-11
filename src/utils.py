import os
from pathlib import Path

from masserstein import Spectrum, NMRSpectrum, estimate_proportions

import numpy as np
import pandas as pd
import time
from matplotlib import pyplot as plt
import pickle
import seaborn as sns
from textwrap import wrap
import pulp



def get_mix_spectrum(nr_of_experiment, experiments_folders, variant=''):

	if (nr_of_experiment != 9 and nr_of_experiment != 4 and nr_of_experiment != 10):
	    filename = experiments_folders['experiment_' + str(nr_of_experiment)] + '/preprocessed_mix.csv'
	    
	elif nr_of_experiment == 10:
		assert variant != ''
		filename = experiments_folders['experiment_' + str(nr_of_experiment)] + '/preprocessed_mix_variant_'+str(variant+1)+'.csv'

	else:
		filename = experiments_folders['experiment_' + str(nr_of_experiment)] + '/preprocessed_exp'+str(nr_of_experiment)+'_mix.csv'

	mix = np.loadtxt(filename, delimiter=',')
	mix = NMRSpectrum(confs=list(zip(mix[:,0], mix[:,1])))

	return mix