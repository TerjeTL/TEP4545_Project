from operator import index
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data_folder = "./Data/LT8.0/WD2.0/"
data_path = data_folder + "export_data.csv"

def calculate_strouhal_num(dataframe, start_time):
    
    df = dataframe.loc[start_time:].copy()

    #print(df)

    dt = df.index[1] - df.index[0]
    print("dt =", dt)
    
    Fs = 1.0/dt
    L = len(df)

    df["Y"] = np.fft.fft(df["Cl"].to_numpy())
    df["P2"] = np.abs(df["Y"].to_numpy()/L)
    P1 = df.iloc[:L//2]["P2"].to_numpy()
    P1[1:-2] = 2*P1[1:-2]

    idx_peak = np.argmax(P1)

    idx_list = np.linspace(0, 1, L//2)
    f = Fs*idx_list/L

    shedding_frequency = 1./2. * 1./f[idx_peak]
    #print(1/shedding_frequency)

    print("Shedding frequency: ", shedding_frequency)
    
    St = shedding_frequency * 0.01/1.4776 # Lc = t*LT, U = 1.4776

    plt.plot(f, P1)
    plt.show()
    return St


df = pd.read_csv(data_path, sep="\t", comment="#", index_col=0)

#print(data)
St_num = calculate_strouhal_num(df, 5.0)

print("Strouhal number: ", St_num)

#df.plot(kind="line", y="Cd")

plt.show()
