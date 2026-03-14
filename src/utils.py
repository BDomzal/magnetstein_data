import os
from pathlib import Path

from masserstein import Spectrum, NMRSpectrum, estimate_proportions

import numpy as np
import pandas as pd
import time
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap, PowerNorm

import pickle
import seaborn as sns
from textwrap import wrap
import pulp


def get_mix_spectrum(nr_of_experiment, experiments_folders, variant=0):

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


def get_components_spectra(nr_of_experiment, experiments_folders, components_dictionary, protons_dictionary, variant=0):

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


def asign_global_names(names, spectra):

	for i, name in enumerate(names):
		globals()[name] = spectra[i]


def save_proportions(estimation, nr_of_experiment, components_dictionary, experiments_folders, filename, variant=''):

	res = pd.DataFrame(np.array(estimation['proportions']).reshape(1,-1), 
	             columns = components_dictionary['experiment_'+str(nr_of_experiment)])

	if (nr_of_experiment != 9 and nr_of_experiment != 4 and nr_of_experiment != 10):
	    res.to_csv(experiments_folders['experiment_' + str(nr_of_experiment)] +
	                   filename + '.csv')
	elif nr_of_experiment ==10:
	    res.to_csv(experiments_folders['experiment_' + str(nr_of_experiment)] +
	                   filename + '_variant_'+str(variant+1) + '.csv')
	else:
	    res.to_csv(experiments_folders['experiment_' + str(nr_of_experiment)] +
	                   filename + '_exp' + str(nr_of_experiment) + '.csv')


def save_noise_in_components(estimation, nr_of_experiment, experiments_folders, filename, variant=''):

	comp_noise = np.array(estimation['proportion_of_noise_in_components']).reshape(-1)

	if (nr_of_experiment != 9 and nr_of_experiment != 4 and nr_of_experiment != 10):

	    np.savetxt(experiments_folders['experiment_' + \
	                                    str(nr_of_experiment)] + \
	                                    filename + \
	                                    '_default_parameters.csv',
	                                      comp_noise)
	elif nr_of_experiment == 10:

	    np.savetxt(experiments_folders['experiment_' + \
	                                    str(nr_of_experiment)] + \
	                                    filename + \
	                                    '_variant' + str(variant+1) + \
	                                    '_default_parameters.csv',
	                                      comp_noise)
	else:

	    np.savetxt(experiments_folders['experiment_' + \
	                                    str(nr_of_experiment)] + \
	                                    filename + \
	                                    '_exp' + str(nr_of_experiment) + \
	                                    '_default_parameters.csv',
	                                      comp_noise)

def save_noise_in_mixture(estimation, nr_of_experiment, experiments_folders, filename, variant=''):

	mix_noise = np.array(1 - sum(estimation['proportions'])).reshape(-1)

	if (nr_of_experiment != 9 and nr_of_experiment != 4 and nr_of_experiment != 10):
	    np.savetxt(experiments_folders['experiment_' + \
	                                    str(nr_of_experiment)] + \
	                                    filename + \
	                                    '_default_parameters.csv',
	                                      mix_noise)
	elif nr_of_experiment == 10:
	    np.savetxt(experiments_folders['experiment_' + \
	                                    str(nr_of_experiment)] + \
	                                    filename + \
	                                    '_variant' + str(variant+1) + \
	                                    '_default_parameters.csv',
	                                      mix_noise)
	else:
	    np.savetxt(experiments_folders['experiment_' + \
	                                    str(nr_of_experiment)] + \
	                                    filename + \
	                                    '_exp' + str(nr_of_experiment) + \
	                                    '_default_parameters.csv',
	                                      mix_noise)

def get_all_componets_results_concentrations(experiments_folders, ground_truth_molar_proportions, protons_dictionary):

	list_of_results_for_experiments = []

	for nr_of_experiment in range(1,10):
	    
	    #MOLAR
	    if (nr_of_experiment != 9 and nr_of_experiment != 4):
	        with open(experiments_folders['experiment_' + str(nr_of_experiment)] + '/results_for_different_kappas.pkl', 'rb') as f:
	            list_of_dataframes_with_results = pd.read_pickle(f)
	    else:
	        with open(experiments_folders['experiment_'+str(nr_of_experiment)] + '/results_for_different_kappas_exp' + str(nr_of_experiment) + '.pkl', 'rb') as f:
	            list_of_dataframes_with_results = pd.read_pickle(f)
	            
	    molar_proportions = ground_truth_molar_proportions['experiment_' + str(nr_of_experiment)]
	    protons = protons_dictionary['experiment_' + str(nr_of_experiment)]
	    
	    temp = []
	    for nr_of_component, df in enumerate(list_of_dataframes_with_results):
	        temp.append(df/protons[nr_of_component])

	    list_of_dataframes_with_results = [df/sum(temp) for df in temp]

	    dataframes_ready_for_heatmap = []

	    for i, df in enumerate(list_of_dataframes_with_results):
	        preprocessed_df = abs(df - molar_proportions[i])
	        dataframes_ready_for_heatmap.append(preprocessed_df)
	        
	    all_components_results_molar = sum(dataframes_ready_for_heatmap)
	    list_of_results_for_experiments.append(all_components_results_molar)

	return list_of_results_for_experiments

	    
def get_all_components_results_areas(experiments_folders, ground_truth_molar_proportions, protons_dictionary):

	list_of_results_for_experiments = []

	for nr_of_experiment in range(1,10):

	    
	    #VISIBLE
	    if (nr_of_experiment != 9 and nr_of_experiment != 4):
	        with open(experiments_folders['experiment_' + str(nr_of_experiment)] + '/results_for_different_kappas.pkl', 'rb') as f:
	            list_of_dataframes_with_results = pd.read_pickle(f)
	    else:
	        with open(experiments_folders['experiment_' + str(nr_of_experiment)] + '/results_for_different_kappas_exp' + str(nr_of_experiment) + '.pkl', 'rb') as f:
	            list_of_dataframes_with_results = pd.read_pickle(f)

	    molar_proportions = ground_truth_molar_proportions['experiment_' + str(nr_of_experiment)]
	    protons = protons_dictionary['experiment_' + str(nr_of_experiment)]
	            
	    real_visible_proportions = [prot*prop for prot, prop in zip(protons, molar_proportions)]
	    real_visible_proportions = [prop/sum(real_visible_proportions) for prop in real_visible_proportions]
	    
	    dataframes_ready_for_heatmap = []
	    for i, df in enumerate(list_of_dataframes_with_results):
	        preprocessed_df = abs(df - real_visible_proportions[i])
	        dataframes_ready_for_heatmap.append(preprocessed_df)
	        
	    all_components_results_vis = sum(dataframes_ready_for_heatmap)

	    list_of_results_for_experiments.append(all_components_results_vis)

	return list_of_results_for_experiments


def get_total_error_from_all_experiments(list_of_results_for_experiments):

	shapes_list = [l.shape for l in list_of_results_for_experiments]

	min_shape = min([sh[0] for sh in shapes_list])

	total_error_all_experiments = sum([df.iloc[:min_shape, :min_shape] for df in list_of_results_for_experiments])

	total_error_all_experiments = total_error_all_experiments.apply(pd.to_numeric, errors = 'coerce', axis=0)

	return total_error_all_experiments


def draw_heatmap(
				total_error_all_experiments,
				colors = ['#2D85C5', '#3CA9EE', '#61BDEE', '#A5D3EB', '#E2E2E2', '#DFA693', '#DC6E55', '#E14B32', '#C33726']
				):

	all_components_results_both = total_error_all_experiments


	vmin = all_components_results_both.min().min()
	vmax = all_components_results_both.max().max()

	my_cmap = ListedColormap(colors)
	my_norm = PowerNorm(0.65, vmin, vmax)

	labels = [round(x,3) for x in all_components_results_both.columns]

	ax = sns.heatmap(all_components_results_both.astype(float), yticklabels=labels, cbar=True,
	                square=True, vmin=vmin, vmax=vmax, xticklabels = labels,
	                cmap=my_cmap, norm=my_norm, cbar_ax=None)

	ax.invert_yaxis()

	for ind, label in enumerate(ax.get_xticklabels()):
	    if (ind+1) % 5 == 0:  # every 10th label is kept
	        label.set_visible(True)
	    else:
	        label.set_visible(False)
	        
	for ind, label in enumerate(ax.get_yticklabels()):
	    if (ind+1) % 5 == 0:  # every 10th label is kept
	        label.set_visible(True)
	    else:
	        label.set_visible(False)

	cbar = ax.collections[0].colorbar
	cbar.set_label('Error of estimation', fontsize=15, labelpad=5)
	# minorticks = [0.1, 0.2]
	# cbar.ax.yaxis.set_ticks(minorticks, minor=True)

	plt.xlabel("Kappa components", fontsize=15, labelpad=5)
	plt.ylabel("Kappa mixture", fontsize=15, labelpad=5)

	plt.tight_layout()
	plt.show()