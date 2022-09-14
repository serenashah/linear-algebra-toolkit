import numpy as np
import matplotlib.pyplot as plt

data_orig = np.loadtxt("dataset1.in",skiprows=1)
data_fit1 = np.loadtxt("dataset-fittings/dataset1_linear.out",skiprows=1)
data_fit2 = np.loadtxt("dataset-fittings/dataset1_exp.out",skiprows=1)
data_fit3 = np.loadtxt("dataset-fittings/dataset1_power.out",skiprows=1)

plt.plot(data_orig[:,0],data_orig[:,1],'o')
plt.plot(data_fit1[:,0],data_fit1[:,1],label="$y_{linear}= -4305.2 + 743.45x$")
plt.plot(data_fit2[:,0],data_fit2[:,1],label="$y_{exp}= 0.25418e^{1.069x}$")
plt.plot(data_fit3[:,0],data_fit3[:,1],label="$y_{power}= 0.000337x^{7.3033}$")
plt.legend()
plt.grid()
plt.savefig("dataset1.pdf")
plt.close()

data_origb = np.loadtxt("dataset2.in",skiprows=1)
data_fit1b = np.loadtxt("dataset-fittings/dataset2_linear.out",skiprows=1)
data_fit2b = np.loadtxt("dataset-fittings/dataset2_exp.out",skiprows=1)
data_fit3b = np.loadtxt("dataset-fittings/dataset2_power.out",skiprows=1)

plt.plot(data_origb[:,0],data_origb[:,1],'o')
plt.plot(data_fit1b[:,0],data_fit1b[:,1],label="$y_{linear}= -4139.6 + 566.14x$")
plt.plot(data_fit2b[:,0],data_fit2b[:,1],label="$y_{exp}= 62.118e^{0.24768x}$")
plt.plot(data_fit3b[:,0],data_fit3b[:,1],label="$y_{power}= 1.1374x^{2.9345}$")
plt.legend()
plt.grid()
plt.savefig("dataset2.pdf")
