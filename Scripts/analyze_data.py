from operator import index
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data_folder = "./Data/WD3.0LT5.0/"
data_path = data_folder + "export_data.csv"

def calculate_strouhal_num(dataframe, start_time):
    
    df = dataframe.loc[start_time:].copy()
    print(df)

    dt = df.index[1] - df.index[0]
    print("dt =", dt)
    
    Fs = 1.0/dt
    L = len(df)

    df["Y"] = np.fft.fft(df["Cl"].to_numpy())
    df["P2"] = np.abs(df["Y"].to_numpy()/L)
    P1 = df.iloc[:L//2]["P2"]
    P1[1:-2] = 2*P1[1:-2]

    idx = np.argmax(P1)

    idx_list = np.linspace(0, 1, L//2)
    f = Fs*idx_list/L
    
    print(len(f))
    print(len(P1))
    
    plt.plot(f, P1)
    
    print(df)
    plt.show()
    return dt


df = pd.read_csv(data_path, sep="\t", comment="#", index_col=0)

#print(df)
calculate_strouhal_num(df, 2.0)

df.plot(kind="line", y="Cl")

plt.show()
