import numpy as np
from scipy.special import zeta
import matplotlib.pyplot as plt

x = np.arange(1.5, 2.5
, step = 0.0001)

zeta_array = np.array([zeta(val) for val in x])

diff = np.abs(zeta_array - 2)
closest_index = diff.argmin()
print(zeta_array[closest_index])
print(x[closest_index])

plt.scatter(x, zeta_array)
plt.hlines(2, 1.5, 2.4, colors="red")
plt.vlines(1.73, 1.4, 2.6, colors="green")
plt.xlabel("")
plt.show()