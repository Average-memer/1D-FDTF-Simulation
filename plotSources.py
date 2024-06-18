import numpy as np
import matplotlib.pyplot as plt

#reads in two files and plots them as 
E = r"C:\Users\lions\source\repos\ichhassegit\EySource.txt"
H = r"C:\Users\lions\source\repos\ichhassegit\HxSource.txt"

Edata = np.genfromtxt(E, delimiter=",")
Hdata = np.genfromtxt(H, delimiter=",")

plt.figure()
plt.plot(Edata, color="blue", label="electric field")
plt.plot(Hdata, color="red", label = "magnetic field")
plt.legend(loc='upper right')
