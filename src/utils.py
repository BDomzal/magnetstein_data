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

def get_ppm(sp):
	return np.array(sp.confs)[:,0]

def get_intensities(sp):
	return np.array(sp.confs)[:,1]

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

def get_results_concentrations(
                                nr_of_experiment, 
                                experiments_folders, 
                                ground_truth_molar_proportions, 
                                protons_dictionary,
                                variant=0
                                ):
	    
	#MOLAR
    if (nr_of_experiment != 9 and nr_of_experiment != 4 and nr_of_experiment != 10):
        with open(experiments_folders['experiment_' + str(nr_of_experiment)] + '/results_for_different_kappas.pkl', 'rb') as f:
            list_of_dataframes_with_results = pd.read_pickle(f)
    elif (nr_of_experiment == 9 or nr_of_experiment == 4):
        with open(experiments_folders['experiment_'+str(nr_of_experiment)] + '/results_for_different_kappas_exp' + str(nr_of_experiment) + '.pkl', 'rb') as f:
            list_of_dataframes_with_results = pd.read_pickle(f)
    else:
        with open(experiments_folders['experiment_'+str(nr_of_experiment)] + '/results_for_different_kappas_variant_' + str(variant+1) + '.pkl', 'rb') as f:
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

    all_components_results_molar = all_components_results_molar.apply(pd.to_numeric, errors = 'coerce', axis=0)

    return all_components_results_molar

def get_all_components_results_concentrations(experiments_folders, ground_truth_molar_proportions, protons_dictionary, variant=0):

	list_of_results_for_experiments = []

	for nr_of_experiment in range(1,10):
	    
	    #MOLAR
	    if (nr_of_experiment != 9 and nr_of_experiment != 4 and nr_of_experiment != 10):
	        with open(experiments_folders['experiment_' + str(nr_of_experiment)] + '/results_for_different_kappas.pkl', 'rb') as f:
	            list_of_dataframes_with_results = pd.read_pickle(f)
	    elif (nr_of_experiment == 9 or nr_of_experiment == 4):
	        with open(experiments_folders['experiment_'+str(nr_of_experiment)] + '/results_for_different_kappas_exp' + str(nr_of_experiment) + '.pkl', 'rb') as f:
	            list_of_dataframes_with_results = pd.read_pickle(f)
	    else:
	        with open(experiments_folders['experiment_'+str(nr_of_experiment)] + '/results_for_different_kappas_variant_' + str(variant+1) + '.pkl', 'rb') as f:
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
	    all_components_results_molar = all_components_results_molar.apply(pd.to_numeric, errors = 'coerce', axis=0)

	    list_of_results_for_experiments.append(all_components_results_molar)

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


def draw_heatmap_power_norm(
							nr_of_experiment,
							all_components_results,
							powers_in_power_norm = [0.45, 1.5, 1.3, 0.65, 0.6, 0.5, 0.525, 0.455, 0.65, [0.75, 0.75, 0.75, 0.75], 0.75],
							colors = ['#2D85C5', '#3CA9EE', '#61BDEE', '#A5D3EB', '#E2E2E2', '#DFA693', '#DC6E55', '#E14B32', '#C33726'],
							remove_edge = True,
							variant=0
							):

	#version with mean of molar and visible proportions
	#powers_in_power_norm = [0.33, 0.8, 0.75, 0.55, 0.39, 0.28, 0.58, 0.342, 0.65]

	#version with only molar proportions
	#powers_in_power_norm = [0.45, 1.5, 1.3, 0.65, 0.6, 0.5, 0.525, 0.455, 0.65, 0.75, 0.75],

	all_components_results_both = all_components_results

	if remove_edge and (nr_of_experiment==3 or nr_of_experiment==8):
		all_components_results_both = all_components_results_both.iloc[:30,:30]

	vmin = all_components_results_both.min().min()
	vmax = all_components_results_both.max().max()

	my_cmap = ListedColormap(colors)

	if nr_of_experiment == 10:
	    my_norm = PowerNorm(powers_in_power_norm[nr_of_experiment-1][variant], vmin, vmax)
	else:
	    my_norm = PowerNorm(powers_in_power_norm[nr_of_experiment-1], vmin, vmax)

	if remove_edge:
	    all_components_results_both = all_components_results_both.iloc[1:, 1:]


	labels = [round(x,3) for x in all_components_results_both.columns]
	if nr_of_experiment==3 or nr_of_experiment==8:
	    labels = labels[:30]
	if remove_edge:
	    labels = labels[1:]
	    

	ax = sns.heatmap(all_components_results_both.astype(float), yticklabels=labels, cbar=True,
	                square=True, vmin=vmin, vmax=vmax, xticklabels = labels,
	                cmap=my_cmap, norm=my_norm, cbar_ax=None)

	ax.invert_yaxis()

	for ind, label in enumerate(ax.get_xticklabels()):
	    if (ind + 3) % 5 == 0:  # every 10th label is kept
	        label.set_visible(True)
	    else:
	        label.set_visible(False)
	        
	for ind, label in enumerate(ax.get_yticklabels()):
	    if (ind + 3) % 5 == 0:  # every 10th label is kept
	        label.set_visible(True)
	    else:
	        label.set_visible(False)

	cbar = ax.collections[0].colorbar
	cbar.set_label('Error of estimation', fontsize=15, labelpad=5)
	cbar.ax.tick_params(labelsize=14)
	# minorticks = [0.1, 0.2]
	# cbar.ax.yaxis.set_ticks(minorticks, minor=True)

	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)

	plt.xlabel("Kappa components", fontsize=15, labelpad=5)
	plt.ylabel("Kappa mixture", fontsize=15, labelpad=5)

	plt.tight_layout()

	if nr_of_experiment != 10:
	    plt.savefig('heatmap_experiment_'+str(nr_of_experiment)+'.png', dpi=300)
	else:
	    plt.savefig('heatmap_experiment_'+str(nr_of_experiment)+'_variant_'+str(variant+1)+'.png', dpi=300)

	plt.show()


def plot_components_without_scaling(
									spectra, 
									components_dictionary,
									nr_of_experiment,
									xlims_lower,
									xlims_upper,
									variant=0,
									path_to_save=None,
									colors = ['blue', 'orange', 'green', 'red', 'pink']
									):


	for i, spectrum in enumerate(spectra):
	    fig, ax = plt.subplots()
	    fig.set_size_inches(9, 4, forward=True)

	    ax.set_xlim(xlims_lower[nr_of_experiment-1], xlims_upper[nr_of_experiment-1])
	    #ax.set_ylim(ylims_lower[nr_of_experiment-1], ylims_upper[nr_of_experiment-1])
	    ax.get_yaxis().set_visible(False)

	    ax.spines['top'].set_visible(False)
	    ax.spines['right'].set_visible(False)
	    ax.spines['bottom'].set_visible(False)
	    ax.spines['left'].set_visible(False)


	    plt.xlabel(chr(0x00b9)+'H, ppm', fontsize=15, labelpad=5)

	    ppm = get_ppm(spectrum)
	    intensity = get_intensities(spectrum)

	    ax.plot(ppm, intensity, color=colors[i],
	           label=components_dictionary['experiment_'+str(nr_of_experiment)][i])
	    ax.fill_between(ppm, 0, intensity, color=colors[i], alpha=1.0)
	    ax.invert_xaxis()

	    ax.legend(prop={'size': 12}, loc='upper left')
	    plt.tight_layout()

	    if path_to_save is not None:
		    if nr_of_experiment != 10:
		        plt.savefig(path_to_save + 'experiment'+str(nr_of_experiment)+'/component'+str(i)+'.png', dpi=300)

		    else:
		        plt.savefig(path_to_save + 'experiment'+str(nr_of_experiment)+'/component'+str(i)+'_variant_'+str(variant+1)+'.png', dpi=300)

	    plt.show()


def get_mix_without_noise(mix, estimation, whether_normalize=True):

    ppm = get_ppm(mix)
    mix_ints = get_intensities(mix)
    mix_without_noise = NMRSpectrum(confs = list(zip(ppm, mix_ints - np.array(estimation['noise']))))
    if whether_normalize:
        mix_without_noise.normalize()
    return mix_without_noise


def get_components_without_noise(components, estimation):

    p = estimation['proportions']
    p = np.array(p).reshape(len(p),1)
    
    components_ints = [get_intensities(comp) for comp in components]
    components = [comp.reshape(-1,1) for comp in components_ints]
    components_no_scaling = np.concatenate(components, axis=1)
    components_scaled = components_no_scaling*p[:,0]

    proportions_point_by_point = (components_scaled/np.sum(components_scaled, axis=1).reshape(-1,1))
    proportions_point_by_point = np.nan_to_num(proportions_point_by_point)
    noise_in_ref = np.array(estimation['noise_in_components']).reshape(-1,1)
    noise_split_for_components = proportions_point_by_point*noise_in_ref

    res = components_scaled - noise_split_for_components
    return res


def get_components_up_to(components_without_noise):

    comp_nr = components_without_noise.shape[1]
    return [np.sum(components_without_noise[:,:(i+1)], axis=1) for i in range(comp_nr)]


def plot_components_added_in_estimated_proportions(
                                                    components,
                                                    mixture,
                                                    estimation,
                                                    nr_of_experiment,
                                                    components_dictionary,
                                                    xlims_lower,
                                                    xlims_upper,
                                                    ylims_lower,
                                                    ylims_upper,
                                                    colors = ['blue', 'orange', 'green', 'red', 'pink'],
                                                    include_mixture=False,
                                                    path_to_save=None
                                                    ):
    
    components_without_noise = get_components_without_noise(components, estimation)
    components_up_to = get_components_up_to(components_without_noise)
    ppm = get_ppm(components[0])
    mixture = get_intensities(mixture)
    
    fig, ax = plt.subplots()
    fig.set_size_inches(9, 4, forward=True)
    
    ax.set_xlim(xlims_lower[nr_of_experiment-1], xlims_upper[nr_of_experiment-1])
    ax.set_ylim(ylims_lower[nr_of_experiment-1], ylims_upper[nr_of_experiment-1])
    ax.get_yaxis().set_visible(False)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    plt.xlabel(chr(0x00b9)+'H, ppm', fontsize=15, labelpad=5)
    #plt.ylabel('Relative intensity', fontsize=15, labelpad=10)
    
    for i, comp_up_to in enumerate(components_up_to):
        ax.plot(ppm, comp_up_to, alpha=1.0, color = colors[i], 
                label=components_dictionary['experiment_'+str(nr_of_experiment)][i])
        if i==0:
            ax.fill_between(ppm, 0, comp_up_to, color=colors[0], alpha=1.0)
        else:
            predecessor = components_up_to[i-1]
            ax.fill_between(ppm, predecessor, comp_up_to, color=colors[i], alpha=1.0)
    if include_mixture:
        ax.plot(ppm, mixture, color='black', alpha=1.0, label='Mixture', linewidth=0.5)
    
    ax.invert_xaxis()
    ax.legend(prop={'size': 12}, loc='upper left')
    plt.tight_layout()
    if path_to_save is not None:
        if nr_of_experiment != 10:
            if include_mixture:
                plt.savefig(path_to_save + 'experiment' + str(nr_of_experiment) + '/mixture_and_all_components.png', dpi=300)
            else:
                plt.savefig(path_to_save + 'experiment' + str(nr_of_experiment) + '/all_components.png', dpi=300)
        else:
            if include_mixture:
                plt.savefig(path_to_save + 'experiment' + str(nr_of_experiment) + '/mixture_and_all_components_variant' + str(variant + 1) + '.png', dpi=300)
            else:
                plt.savefig(path_to_save + 'experiment' + str(nr_of_experiment) + '/all_components_variant_' + str(variant + 1) + '.png', dpi=300)
    plt.show()


def plot_mixture(
                mixture,
                nr_of_experiment,
                components_dictionary,
                xlims_lower,
                xlims_upper,
                ylims_lower,
                ylims_upper,
                path_to_save=None
                ):
    
    ppm = get_ppm(mixture)
    mixture = get_intensities(mixture)

    colors = ['blue', 'orange', 'green', 'red', 'pink']
    fig, ax = plt.subplots()
    #fig.set_size_inches(10, 7, forward=True)
    fig.set_size_inches(9, 4, forward=True)
    
    ax.set_xlim(xlims_lower[nr_of_experiment-1], xlims_upper[nr_of_experiment-1])
    ax.set_ylim(ylims_lower[nr_of_experiment-1], ylims_upper[nr_of_experiment-1])
    ax.get_yaxis().set_visible(False)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    
    plt.xlabel(chr(0x00b9)+'H, ppm', fontsize=15, labelpad=5)
    #plt.ylabel('Relative intensity', fontsize=15, labelpad=10)
    
    ax.plot(ppm, mixture, color='black', alpha=1.0, label='Mixture', linewidth=0.5)
    ax.invert_xaxis()
    
    ax.legend(prop={'size': 12}, loc='upper left')
    
    fig.tight_layout()
    if path_to_save is not None:
        if nr_of_experiment != 10:
            fig.savefig(path_to_save + 'experiment' + str(nr_of_experiment) + '/mixture.png', dpi=300)
        else:
            fig.savefig(path_to_save + 'experiment' + str(nr_of_experiment) + '/mixture_variant_' + str(variant + 1) + '.png', dpi=300)