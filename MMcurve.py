import numpy as np
import matplotlib.pyplot as plt

Cs = np.linspace(0, 100)

k1=1.
k2=1.
k3=2.

KH=(k2+k3)/k1

def MM(Cs):
	return k3*Cs/(Cs+KH)

plt.plot(Cs, MM(Cs))
plt.show()