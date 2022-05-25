import numpy as np
import matplotlib.pyplot as plt

u = np.loadtxt("u.txt")

Nx, Ny = u.shape

yf = np.linspace(0, 1, Ny+1)
yc = (yf[:-1] + yf[1:]) / 2

uc = u[Nx//2, :]

yg = np.array([0, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1])
ug = np.array([0, -0.18109, -0.20196, -0.2222, -0.2973, -0.38289, -0.27805, -0.10648, -0.06080, 0.05702, 0.18719, 0.33304, 0.46604, 0.51117, 0.57492, 0.65928, 1])

plt.plot(yc, uc)
plt.plot(yg, ug, 'ro')
plt.show()
