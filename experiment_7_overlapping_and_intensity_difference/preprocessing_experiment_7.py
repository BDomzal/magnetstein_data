
nr_of_experiment = 7

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
MASS_PROPORTIONS = config["mass_proportions"]
MOLAR_MASSES = config["molar_masses_or_proportions"]

MIX_RAW_FILE = config["mix_raw_files"]["experiment_" + str(nr_of_experiment)]
COMPONENTS_RAW_FILE = config["components_raw_files"]["experiment_" + str(nr_of_experiment)]
official_names = COMPONENTS_DICTIONARY["experiment_" + str(nr_of_experiment)]
protons = PROTONS_DICTIONARY["experiment_" + str(nr_of_experiment)]
mass_proportions = MASS_PROPORTIONS["experiment_" + str(nr_of_experiment)]
molar_proportions = MOLAR_MASSES["experiment_" + str(nr_of_experiment)]


# ### Loading the data


mix = np.loadtxt(MIX_RAW_FILE, skiprows=0, usecols=range(2), 
                      delimiter='\t',dtype=float)
ppm = mix[:,0]
mix = mix[:,1]

comp = np.loadtxt(COMPONENTS_RAW_FILE, skiprows=1, usecols=[0, 2, 4],
                      delimiter='\t',dtype=float)
ppm_components = comp[:,0]
components = comp[:,1:]


components_ints = []
names = []
for i in range(components.shape[1]):
    components_ints.append(components[:,i])
    names.append('comp'+str(i))


spectra = []
for comp_ints in components_ints:
    spectra.append(NMRSpectrum(confs=list(zip(ppm_components, comp_ints))))


mix = NMRSpectrum(confs=list(zip(ppm,mix)))


NMRSpectrum.plot_all([mix]+spectra, profile=True)


# ### Ground truth

molar_proportions = [mass_prop/mol_mass for mass_prop, mol_mass in zip(mass_proportions, molar_proportions)]
molar_proportions = [mol_prop/sum(molar_proportions) for mol_prop in molar_proportions]
print('Molar proportions:')
print(molar_proportions)


real_visible_proportions = [prot*prop for prot, prop in zip(protons, molar_proportions)]
real_visible_proportions = [prop/sum(real_visible_proportions) for prop in real_visible_proportions]
print('Area-based proportions:')
print(real_visible_proportions)


# ### Preprocessing

spectra_and_mixture = spectra + [mix]

preprocessed_spectra = []
for sp in spectra_and_mixture:
    sp2 = make_nonnegative(sp)
    sp2.sort_confs()
    sp2.merge_confs()
    sp2.normalize()
    preprocessed_spectra.append(sp2)


spectra_and_mixture = preprocessed_spectra
spectra = spectra_and_mixture[:-1]
mix = preprocessed_spectra[-1]
del(preprocessed_spectra)


cma = common_ppm_axis(spectra_and_mixture)

preprocessed_spectra = []
for sp in spectra_and_mixture:
    preprocessed_spectra.append(sp.resample(cma))


spectra_and_mixture = preprocessed_spectra
del(preprocessed_spectra)


preprocessed_spectra = []
for sp in spectra_and_mixture:
    sp = make_nonnegative(sp)
    sp.sort_confs()
    sp.merge_confs()
    sp.normalize()
    preprocessed_spectra.append(sp)


spectra_and_mixture = preprocessed_spectra
spectra = spectra_and_mixture[:-1]
mix = preprocessed_spectra[-1]
del(preprocessed_spectra)


for i, sp in enumerate(spectra):
    sp.protons = protons[i]


for i, name in enumerate(names):
    globals()[name] = spectra[i]


components_ints = []
for spectrum in spectra:
    components_ints.append(np.array(spectrum.confs)[:,1])#.reshape(-1,1))


ppm = np.array(mix.confs)[:,0]
mix_ints = np.array(mix.confs)[:,1]


# ### Finding best kappa

lower_limit = 0.01
upper_limit = 0.3
step = 0.01


# list_of_estimation_results = []
# for kappa in np.arange(start=lower_limit, stop=upper_limit, step=step):
#     fixed_kappa = []
#     for kappa_prime in np.arange(start=lower_limit, stop=upper_limit, step=step):
#         start = time.time()
#         estimation = estimate_proportions(mix, [comp0, comp1], 
#                                    MTD=kappa, MTD_th=kappa_prime, verbose=False, solver=pulp.GUROBI(msg=False),
#                                        what_to_compare='area')
#         end = time.time()
#         fixed_kappa.append(estimation['proportions'])
#         print('Estimation for '+str(kappa)+' and '+str(kappa_prime)+' done')
#         print('It took: '+str(end-start)+' seconds.')
#         print(estimation['proportions'])
#     list_of_estimation_results.append(fixed_kappa)


# list_of_dataframes_with_results = []
# for component_number in range(len([comp0, comp1])):
#     results_for_component = pd.DataFrame(columns=np.arange(start=lower_limit, stop=upper_limit, step=step),
#                                    index=np.arange(start=lower_limit, stop=upper_limit, step=step))
#     for i in range(results_for_component.shape[0]):
#         for j in range(results_for_component.shape[1]):
#             results_for_component.iloc[i,j] = list_of_estimation_results[i][j][component_number]
    
#     list_of_dataframes_with_results.append(results_for_component)


# with open('results_for_different_kappas.pkl', 'wb') as f:
#     pickle.dump(list_of_dataframes_with_results, f)


# #### Molar proportions

with open('results_for_different_kappas.pkl', 'rb') as f:
    list_of_dataframes_with_results = pd.read_pickle(f)
#first coordinate: kappa,
#second coordinate: kappa_prime


temp = []
for nr_of_component, df in enumerate(list_of_dataframes_with_results):
    temp.append(df/protons[nr_of_component])
temp2 = []
for df in temp:
    temp2.append(df/sum(temp))
list_of_dataframes_with_results = temp2
del(temp)
del(temp2)


dataframes_ready_for_heatmap = []
for i, df in enumerate(list_of_dataframes_with_results):
    preprocessed_df = abs(df - molar_proportions[i])
    dataframes_ready_for_heatmap.append(preprocessed_df)


all_components_results_molar = sum(dataframes_ready_for_heatmap)


# #### Visible proportions

with open('results_for_different_kappas.pkl', 'rb') as f:
    list_of_dataframes_with_results = pd.read_pickle(f)
#first coordinate: kappa,
#second coordinate: kappa_prime


dataframes_ready_for_heatmap = []
for i, df in enumerate(list_of_dataframes_with_results):
    preprocessed_df = abs(df - real_visible_proportions[i])
    dataframes_ready_for_heatmap.append(preprocessed_df)


all_components_results_vis = sum(dataframes_ready_for_heatmap)


# #### Both together

all_components_results_both = all_components_results_vis + all_components_results_molar


all_components_results_both = all_components_results_both.apply(pd.to_numeric, errors = 'coerce', axis=0)


print('Minimal:')
print(all_components_results_both.min().min())
print(all_components_results_both[0.04][0.30])


# ### Saving preprocessed spectra

for i, sp in enumerate(spectra_and_mixture):
    try:
        np.savetxt('preprocessed_'+str(names[i])+'.csv', np.array(sp.confs), delimiter=',')
    except IndexError:
        np.savetxt('preprocessed_mix.csv', np.array(sp.confs), delimiter=',')

