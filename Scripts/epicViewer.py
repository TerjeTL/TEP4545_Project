import matplotlib.pyplot as plt
import pandas as pd

case = 'WD1.0LT5.0'
filename = 'Data/{case}/export_data.csv'.format(case=case)
df = pd.read_csv(filename, sep="\t", comment="#")

print(df)

df.plot(kind="line", x="Time", y="Cl", color="blue")
df.plot(kind="line", x="Time", y="Cd", color="red")

plt.show()
