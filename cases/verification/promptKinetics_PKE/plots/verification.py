import numpy as np
import matplotlib.pyplot as plt
import re
import figStyle as fs

figSize = 'small'

# Prepare figure

fs.prep(plt, figSize)

#

times = []
psi = []

current_time = None

with open("../log", "r") as f:
    for line in f:
        if line.startswith("Time ="):
            current_time = float(line.split("=")[1].strip().rstrip("s"))

        elif current_time is not None:
            if line.startswith("Psi:"):
                psi.append(float(line.split()[1]))
                times.append(current_time)

# Convert lists to numpy arrays for calculation
times = np.array(times)
psi = np.array(psi)

# Check if data was parsed successfully
if len(times) == 0:
    print("Error: No data parsed. Check if your log file contains the specified format.")
    exit()

# Analytical solution

reactivity = 1e-5

L = 1e-5

psi_analytical = np.exp(reactivity/L * times)


# Plots

fig = plt.figure("flux")
plt.plot(times, psi, label="MUSCLE-Foam")
plt.plot(times, psi_analytical, ':k', label='analytical')
plt.xlabel("t [s]")
plt.ylabel(r"$\Psi$")
fs.post(fig, figSize, plt.legend())
plt.savefig("psi.png",dpi=1000)
