nr_of_experiment = 5

import numpy as np

#Columns: Ground truth, Magnetstein, ACD/TP, MANIQ

#Exp. 1
#est1 = np.array([[0.0909, 0.0937, 0.0570, 0.1135],[0.9091, 0.9063, 0.9430, 0.8946]])

#after reducing to two significant figures
est1 = np.array([[0.09, 0.09, 0.06, 0.11],[0.91, 0.91, 0.94, 0.89]])


#Exp. 2
# est2 = np.array([[0.4950, 0.4760, 0.4588, 0.4864],
# [0.5050, 0.5240, 0.5412, 0.5136]])

#after reducing to two significant figures
est2 = np.array([[0.50, 0.48, 0.46, 0.49],
[0.51, 0.52, 0.54, 0.51]])


#Exp. 3
# est3 = np.array([[0.1058, 0.1059, 0.1086, 0.1110],
# [0.7265, 0.7511, 0.7133, 0.7578],
# [0.0858, 0.0621, 0, 0.0951],
# [0.0820, 0.0809, 0.1781, 0.0025]])

#after reducing to two significant figures
est3 = np.array([[0.11, 0.11, 0.11, 0.11],
[0.73, 0.75, 0.71, 0.76],
[0.09, 0.06, 0, 0.10],
[0.08, 0.08, 0.18, 0.00]])


#Exp. 4
# est4 = np.array([[0.3022, 0.3165, 0.3488, 0.3311],
# [0.2240, 0.2137, 0.2793, 0.2660],
# [0.1253, 0.1341, 0.1921, 0],
# [0.2028, 0.1960, 0.0109, 0.2717],
# [0.1457, 0.1397, 0.1689, 0.1578]])

#after reducing to two significant figures
est4 = np.array([[0.30, 0.32, 0.35, 0.33],
[0.22, 0.21, 0.28, 0.27],
[0.13, 0.13, 0.19, 0.00],
[0.20, 0.20, 0.01, 0.27],
[0.15, 0.14, 0.17, 0.16]])


#Exp. 5
# est5 = np.array([[0.3022, 0.3180, 0.3052, 0.2780],
# [0.2240, 0.2113, 0.2377, 0.2175],
# [0.1253, 0.1293, 0.0006, 0.1407],
# [0.2028, 0.2018, 0.3156, 0.2255],
# [0.1457, 0.1397, 0.1409, 0.1368]])

#after reducing to two significant figures
est5 = np.array([[0.30, 0.32, 0.31, 0.28],
[0.22, 0.21, 0.24, 0.22],
[0.13, 0.13, 0.00, 0.14],
[0.20, 0.20, 0.32, 0.23],
[0.15, 0.14, 0.14, 0.14]])

all_exps = [est1, est2, est3, est4, est5]

est = all_exps[nr_of_experiment-1]

avre_error = []
for i in range(1, 4):
    error = np.abs(est[:,0]-est[:,i])/est[:,0]
    avre_error.append(np.mean(error))


print('AvRE for Magnetstein; ACD/TP; MANIQ for Experiment ' + str(nr_of_experiment) + ':')
print(avre_error)


