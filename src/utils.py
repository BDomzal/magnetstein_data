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

def get_components_spectra(nr_of_experiment, experiments_folders, components_dictionary, protons_dictionary, variant=''):

	how_many_components = len(components_dictionary['experiment_'+str(nr_of_experiment)])
	names = ['comp' + str(i) for i in range(how_many_components)]

	spectra_np = []

	for i in range(how_many_components):

	    if (nr_of_experiment != 9 and nr_of_experiment != 4 and nr_of_experiment != 10):
	        filename = experiments_folders['experiment_' + str(nr_of_experiment)] + \
	                    '/preprocessed_comp' + str(i) + '.csv'
	    elif nr_of_experiment == 10:
	    	assert variant != ''
	        filename = experiments_folders['experiment_' + str(nr_of_experiment)] + '/preprocessed_variant_'+str(variant+1)+'_comp'+str(i)+'.csv'
	    else:
	        filename = experiments_folders['experiment_' + str(nr_of_experiment)] + \
	                    '/preprocessed_exp'+str(nr_of_experiment)+ '_comp'+ str(i) + '.csv'
	    spectra_np.append(np.loadtxt(filename, delimiter=','))

	spectra = []
	names = []

	for i in range(len(spectra_np)):
	    spectra.append(NMRSpectrum(confs=list(zip(spectra_np[i][:,0], spectra_np[i][:,1])), 
	                                protons=protons_dictionary['experiment_'+str(nr_of_experiment)][i]))
	    names.append('comp'+str(i))

	return spectra, names, how_many_components

def assign_global_names(name, spectra):

	for i, name in enumerate(names):
    	globals()[name] = spectra[i]