from operator import index
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from glob import glob

#### CONFIG ####

# Use glob to catch all wanted case folders as a list (note: only folders to be included in the grid refinement)
# i.e. [grid_case_1, grid_case_2, grid_case_3]
dirs = glob("./Data/grid*/")
output_path = "./Data/GCI.csv"

################


# Just get the mean value of the Cd after a certain startup time (check where the solutions stabilize)
def calculate_mean_cd(dataframe, start_time):
    
    df = dataframe.loc[start_time:].copy()
    print(df)

    return df["Cd"].mean()

# Prepare a results dataframe
results = pd.DataFrame(index=dirs, columns=["volume", "N", "h", "r", "Cd_mean"])

# Store some relevant case data
for index, row in results.iterrows():
    with open(index + "num_cells.dat", 'r') as f:
        num_cells = float(f.read())
        print(num_cells)

    with open(index + "total_volume.dat", 'r') as f:
        volume = float(f.read())
        print(volume)

    results.loc[index, 'volume'] = volume
    results.loc[index, 'N'] = num_cells

# cell size
results["h"] = (results["volume"]/results["N"])**(1./3.)

# Sort the table based on cell size
results = results.sort_values('h')

# Refinement factor
results["r"] = results.iloc[:, 2].shift(-1)/results.iloc[:, 2]

for index, row in results.iterrows():
    dft = pd.read_csv(index + "export_data.csv", sep="\t", comment="#", index_col=0)
    
    ### CAUTION: Note the start time ###
    cd_mean = calculate_mean_cd(dft, 2.0)

    results.loc[index, 'Cd_mean'] = cd_mean

results["epsilon"] = results.loc[:, "Cd_mean"].shift(-1) - results.loc[:, "Cd_mean"]

e_21 = results.iloc[0, 5]
e_32 = results.iloc[1, 5]
r_21 = results.iloc[0, 3]
r_32 = results.iloc[1, 3]

s = np.sign(e_32/e_21)
results["sgn"] = s

# Find order of the method based on the three grids
def eqn(p):
    return 1.0/np.log(r_21) * np.abs(np.log(np.abs(e_32/e_21))
        + np.log((r_21**p - s)/(r_32**p - s))) - p

p = fsolve(eqn, 1.0)[0]
results["p"] = s
print("Estimated order: ", p)

phi_1 = results.iloc[0, 4]
phi_2 = results.iloc[1, 4]
phi_3 = results.iloc[2, 4]

ext_Cd_21 = (r_21**p * phi_1 - phi_2)/(r_21**p - 1)
ext_Cd_32 = (r_32**p * phi_2 - phi_3)/(r_32**p - 1)
results["Cd_ext"] = [ext_Cd_21, ext_Cd_32, np.nan]

approx_e_21 = np.abs((phi_1 - phi_2)/phi_1)
approx_e_32 = np.abs((phi_2 - phi_3)/phi_2)
results["approx_err"] = [approx_e_21, approx_e_32, np.nan]

ext_e_21 = np.abs((ext_Cd_21 - phi_1)/ext_Cd_21)
ext_e_32 = np.abs((ext_Cd_32 - phi_2)/ext_Cd_32)
results["ext_err"] = [ext_e_21, ext_e_32, np.nan]

safety_factor = 1.25 #Roche (3 grid study)
GCI_f_21 = safety_factor*approx_e_21/(r_21**p - 1)
GCI_f_32 = safety_factor*approx_e_32/(r_32**p - 1)
results["GCI_fine"] = [GCI_f_21, GCI_f_32, np.nan]

pd.set_option("display.precision", 8)
print(results)

results.to_csv(output_path)

# uncomment for plotting

#for index, row in results.iterrows():
#    dft = pd.read_csv(index + "export_data.csv", sep="\t", comment="#", index_col=0)
#
#    dft.plot(kind="line", y="Cd")
#    plt.show()
#
