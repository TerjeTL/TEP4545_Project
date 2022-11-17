import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("./../data/some_other_file.csv", sep="\t", comment="#")

print(df)

df.plot(kind="line", x="Time", y="Cl", color="blue")
df.plot(kind="line", x="Time", y="Cd", color="red")

plt.show()
