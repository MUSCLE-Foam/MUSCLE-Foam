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

endTime = "1000" # Change if necessary

csv_flux = read_csv("../postProcessing/horizontal/"+endTime+"/horizontal.csv")

x = csv_flux['x'].tolist()

flux = csv_flux['totalFlux0'].tolist()

# Analytical solution

sigma_a = 10
psi_in = 100/2

x_an = np.linspace(0, 0.5, 1000)
flux_an = np.array([
        quad(lambda mu: psi_in * np.exp(-sigma_a * xi / mu), 0, 1)[0]
        for xi in x_an
    ])

# Plot

fig = plt.figure("flux")
plt.plot(x, flux, color='tab:green', label='MUSCLE-Foam')
plt.plot(x_an, flux_an, ':k', label='semi-analytical')
plt.xlabel(r'$x \ [\mathrm{m}]$')
plt.ylabel(r'$\phi_{tot}$')
fs.post(fig, figSize, plt.legend())
plt.savefig("flux.png",dpi=1000)

# Angular fluxes

psi_in = 100
neutronFlightDirectionsDict = "../constant/neutronFlightDirectionsDict"

# Extract nDirections from dict
with open(neutronFlightDirectionsDict, 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith('nDirections'):
            start = line.find('nDirections')
            end = line.find(';', start)
            print("nDirections =",line[start+11:end].strip())
            nDirections = int(line[start+11:end].strip())
            break

mu = 0

for i in range(nDirections):
    fig = plt.figure("flux"+str(i))
    phi = csv_flux['flux0_'+str(i)].tolist()
    plt.plot(x, phi, color='tab:green', label='MUSCLE-Foam')

    target_prefix = f'direction{i} '

    # Extract mu from dict
    with open(neutronFlightDirectionsDict, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(target_prefix):
                start = line.find('(')
                end = line.find(')', start)
                if start != -1 and end != -1:
                    vector_str = line[start+1:end]
                    x_str = vector_str.split()[0]
                    try:
                        mu = max(float(x_str),0)
                    except ValueError:
                        raise ValueError(f"Could not convert x-component to float: '{x_str}'")
                break

    phi_an = psi_in * np.exp(-sigma_a * x_an / max(mu,1e-5))

    if mu < 1e-5:
        phi_an = [0 for p in phi_an]

    plt.plot(x_an, phi_an, ':k', label='analytical')

    plt.xlabel(r'$x \ [\mathrm{m}]$')
    plt.ylabel(r'$\varphi(\mu)$')

    fs.post(fig, figSize, plt.legend())

    plt.savefig("phi"+str(i)+".png",dpi=1000)
    print(i, "mu =", mu)

    fig = plt.figure("phi")
    plt.plot(x, phi, color='tab:green')
    plt.plot(x_an, phi_an, ':k')


# dummy lines for legend
plt.plot(0,0, color='tab:green', label='MUSCLE-Foam')
plt.plot(0,0, ':k', label='analytical')

plt.xlabel(r'$x \ [\mathrm{m}]$')
plt.ylabel(r'$\varphi(\mu)$')
fs.post(fig, figSize, plt.legend())
plt.savefig("phi.png",dpi=1000)

