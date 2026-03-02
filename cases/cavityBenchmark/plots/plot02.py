#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
import os, sys
from pandas import *

# Benchmark results

csv_benchmark_NFD = read_csv("benchmark/Step02/fiss_dens_AA", skiprows=1, delim_whitespace=True, names=["x", "CNRS", "CNRS2", "PoliMi", "PSI", "TUD", "TUD2"])

benchmark_x = csv_benchmark_NFD['x'].tolist()

benchmark_NFD = csv_benchmark_NFD['PoliMi'].tolist()

# Simulation data

endTime = "20000" # Change if necessary

csv_flux = read_csv("../postProcessing/horizontal/"+endTime+"/horizontal.csv")

x = csv_flux['x'].tolist()

flux0 = csv_flux['flux0'].tolist()
flux1 = csv_flux['flux1'].tolist()
flux2 = csv_flux['flux2'].tolist()
flux3 = csv_flux['flux3'].tolist()
flux4 = csv_flux['flux4'].tolist()
flux5 = csv_flux['flux5'].tolist()

sigma_f = [0.111309, 0.108682, 0.152219, 0.258190, 0.536326, 1.44917]

x = [x_+1 for x_ in x]

NFD = []

for i in range(len(flux0)):
    NFD.append(flux0[i]*sigma_f[0]+flux1[i]*sigma_f[1]+flux2[i]*sigma_f[2]+flux3[i]*sigma_f[3]+flux4[i]*sigma_f[4]+flux5[i]*sigma_f[5])

# Plot

fig = plt.figure("02")
plt.plot(x, NFD, color='tab:green', label='MUSCLE-Foam')
plt.plot(benchmark_x, benchmark_NFD, ':k', label='PoliMi')
plt.xlabel(r'$x [m]$')
plt.ylabel(r'$NFD \ [m^{-3}s^{-1}]$')
plt.legend()
plt.savefig("02_horizontal.png",dpi=1000)