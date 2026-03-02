#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
import os, sys
from pandas import *
from scipy.integrate import quad
from scipy.special import expn

import figStyle as fs

figSize = 'small'

# Prepare figure

fs.prep(plt, figSize)

# Simulation data

endTime = "10000" # Change if necessary

csv_flux = read_csv("../postProcessing/horizontal/"+endTime+"/horizontal.csv")

x = csv_flux['x'].tolist()

flux0 = csv_flux['flux0'].tolist()
flux1 = csv_flux['flux1'].tolist()

maxFlux0 = max(flux0)
maxFlux1 = max(flux1)

flux0 = [f/maxFlux0 for f in flux0]
flux1 = [f/maxFlux1 for f in flux1]

# Benchmark solution

csv_benchmark0 = read_csv("benchmark/fastFlux.csv")
csv_benchmark1 = read_csv("benchmark/thermalFlux.csv")

x_b0 = csv_benchmark0['x'].tolist()
x_b1 = csv_benchmark1['x'].tolist()

flux0_b = csv_benchmark0['flux0'].tolist()
flux1_b = csv_benchmark1['flux1'].tolist()

maxFlux0_b = max(flux0_b)
maxFlux1_b = max(flux1_b)

flux0_b = [f/maxFlux0_b for f in flux0_b]
flux1_b = [f/maxFlux1_b for f in flux1_b]

# Plot

fig = plt.figure("flux0")
plt.plot(x, flux0, color='tab:green', label='MUSCLE-Foam')
plt.plot(x_b0, flux0_b, ':k', label='Benchmark')
plt.xlabel(r'$x \ [\mathrm{m}]$')
plt.ylabel(r'$\phi_{0} [-]$')
fs.post(fig, figSize, plt.legend())
plt.savefig("flux0.png",dpi=1000)

fig = plt.figure("flux1")
plt.plot(x, flux1, color='tab:green', label='MUSCLE-Foam')
plt.plot(x_b1, flux1_b, ':k', label='Benchmark')
plt.xlabel(r'$x \ [\mathrm{m}]$')
plt.ylabel(r'$\phi_{1} [-]$')
fs.post(fig, figSize, plt.legend())
plt.savefig("flux1.png",dpi=1000)

