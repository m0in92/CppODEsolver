import pandas as pd
import matplotlib.pyplot as plt


y = pd.read_csv("cmake-build-debug/results.csv", header=None).to_numpy()

plt.plot(y)
plt.show()