#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
import os, sys
from pandas import *

# Benchmark results

names_AA = ["x", "CNRS", "CNRS2", "PoliMi", "PSI", "TUD", "TUD2"]
names_BB = ["y", "CNRS", "CNRS2", "PoliMi", "PSI", "TUD", "TUD2"]

csv_benchmark_DNP_AA = read_csv("benchmark/Step13/DNP_AA", skiprows=1, sep=r"\s+", names=names_AA)
csv_benchmark_DNP_BB = read_csv("benchmark/Step13/DNP_BB", skiprows=1, sep=r"\s+", names=names_BB)
csv_benchmark_T_AA = read_csv("benchmark/Step13/T_AA", skiprows=1, sep=r"\s+", names=names_AA)
csv_benchmark_T_BB = read_csv("benchmark/Step13/T_BB", skiprows=1, sep=r"\s+", names=names_BB)

benchmark_x = csv_benchmark_DNP_AA['x'].tolist()
benchmark_y = csv_benchmark_DNP_BB['y'].tolist()

benchmark_DNP_AA = csv_benchmark_DNP_AA['PoliMi'].tolist()
benchmark_DNP_BB = csv_benchmark_DNP_BB['PoliMi'].tolist()
benchmark_T_AA = csv_benchmark_T_AA['PoliMi'].tolist()
benchmark_T_BB = csv_benchmark_T_BB['PoliMi'].tolist()

# Simulation data

endTime = "20000" # Change if necessary

csv_horizontal = read_csv("../postProcessing/horizontal/"+endTime+"/horizontal.csv")
csv_vertical = read_csv("../postProcessing/vertical/"+endTime+"/vertical.csv")

# Shift the coordinates from [-1, 1] (mesh) to [0, 2] (benchmark)
x = [x_+1 for x_ in csv_horizontal['x'].tolist()]
y = [z_+1 for z_ in csv_vertical['z'].tolist()]

# DNS<e> is the delayed neutron source of energy group e, i.e. the total
# delayed neutron source multiplied by the (uniform) delayed neutron spectrum.
# Dividing by chiDel0 (see 0.orig/chiDel0) recovers the total source.
chiDel0 = 0.00430325

DNP_AA = [DNS_/chiDel0 for DNS_ in csv_horizontal['DNS0'].tolist()]
DNP_BB = [DNS_/chiDel0 for DNS_ in csv_vertical['DNS0'].tolist()]

T_AA = csv_horizontal['T'].tolist()
T_BB = csv_vertical['T'].tolist()

# Plot

fig = plt.figure("13_DNP_AA")
plt.plot(x, DNP_AA, color='tab:green', label='MUSCLE-Foam')
plt.plot(benchmark_x, benchmark_DNP_AA, ':k', label='PoliMi')
plt.xlabel(r'$x [m]$')
plt.ylabel(r'$DNS \ [m^{-3}s^{-1}]$')
plt.legend()
plt.savefig("13_DNP_horizontal.png",dpi=1000)

fig = plt.figure("13_DNP_BB")
plt.plot(y, DNP_BB, color='tab:green', label='MUSCLE-Foam')
plt.plot(benchmark_y, benchmark_DNP_BB, ':k', label='PoliMi')
plt.xlabel(r'$y [m]$')
plt.ylabel(r'$DNS \ [m^{-3}s^{-1}]$')
plt.legend()
plt.savefig("13_DNP_vertical.png",dpi=1000)

fig = plt.figure("13_T_AA")
plt.plot(x, T_AA, color='tab:green', label='MUSCLE-Foam')
plt.plot(benchmark_x, benchmark_T_AA, ':k', label='PoliMi')
plt.xlabel(r'$x [m]$')
plt.ylabel(r'$T \ [K]$')
plt.legend()
plt.savefig("13_T_horizontal.png",dpi=1000)

fig = plt.figure("13_T_BB")
plt.plot(y, T_BB, color='tab:green', label='MUSCLE-Foam')
plt.plot(benchmark_y, benchmark_T_BB, ':k', label='PoliMi')
plt.xlabel(r'$y [m]$')
plt.ylabel(r'$T \ [K]$')
plt.legend()
plt.savefig("13_T_vertical.png",dpi=1000)
