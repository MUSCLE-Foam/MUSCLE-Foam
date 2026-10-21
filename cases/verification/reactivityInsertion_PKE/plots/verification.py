import numpy as np
import matplotlib.pyplot as plt
import re
import figStyle as fs

figSize = 'small'

# Prepare figure

fs.prep(plt, figSize)

# Fixed kinetic constants matching your solver setup
BETA = 0.007
L_PROMPT = 1e-5
LAMBDA = 0.1

# Read reactivity insertion value

def get_rho_value(file_path):
    with open(file_path, "r") as f:
        content = f.read()

    # Match 'rho' followed by whitespace and a float/int inside reactivityInsertion
    pattern = r"reactivityInsertion\s*\{[^}]*?\brho\s+([0-9\.eE+-]+);"
    match = re.search(pattern, content, re.DOTALL)

    if match:
        return float(match.group(1))
    raise ValueError("rho value not found in reactivityInsertion block")


RHO_INSERTED = get_rho_value("../system/controlDict")

times = []
psi = []
eta = []

current_time = None

with open("../log", "r") as f:
    for line in f:
        if line.startswith("Time ="):
            current_time = float(line.split("=")[1].strip().rstrip("s"))

        elif current_time is not None:
            if line.startswith("Psi:"):
                psi.append(float(line.split()[1]))
                times.append(current_time)
            elif line.startswith("Eta0:"):
                eta.append(float(line.split()[1]))

# Convert lists to numpy arrays for calculation
times = np.array(times)
psi = np.array(psi)
eta = np.array(eta)

# Check if data was parsed successfully
if len(times) == 0:
    print("Error: No data parsed. Check if your log file contains the specified format.")
    exit()

print("========== PKE EXTERNAL VERIFICATION ==========")

# --- TEST 1: PROMPT JUMP APPROXIMATION ---
# Find index right before and right after the t = 0.1s step insertion
idx_before = np.where(times < 0.1)[0][-1]
idx_after  = np.where(times >= 0.11)[0][0] # First time step at or after 0.11s

psi_0 = psi[idx_before]
simulated_jump = psi[idx_after]
expected_jump = psi_0 * (BETA / (BETA - RHO_INSERTED))

jump_error = abs(simulated_jump - expected_jump) / expected_jump * 100

print(f"Prompt Jump Test:")
print(f"  - Expected Flux:  {expected_jump:.5f}")
print(f"  - Simulated Flux: {simulated_jump:.5f}")
print(f"  - Relative Error: {jump_error:.4f}%")
print(f"  --> STATUS: {'PASSED' if abs(jump_error) < 0.5 else 'FAILED'}")
print("-" * 47)

# --- TEST 2: INHOUR REGIME / ASYMPTOTIC GROWTH ---
# Calculate growth rate (omega) at the end of the transient using a linear fit on log(PSI)
tail_idx = -50
omega_sim = np.polyfit(times[tail_idx:], np.log(psi[tail_idx:]), 1)[0]

# Plug simulated growth rate into the analytical Inhour Equation to back-calculate reactivity
rho_back_calculated = (L_PROMPT * omega_sim) + (BETA * omega_sim / (omega_sim + LAMBDA))
inhour_error = abs(rho_back_calculated - RHO_INSERTED) / RHO_INSERTED * 100

print(f"Inhour Asymptotic Test:")
print(f"  - Target Reactivity:       {RHO_INSERTED:.7f}")
print(f"  - Back-calculated Rho:     {rho_back_calculated:.7f}")
print(f"  - Relative Error:          {inhour_error:.4f}%")
print(f"  --> STATUS: {'PASSED' if abs(inhour_error) < 0.5 else 'FAILED'}")
print("===============================================")

# Plots

fig = plt.figure("flux")
plt.plot(times, psi, label="MUSCLE-Foam")
plt.axhline(y=expected_jump, color='black', linestyle=':', label="analytical")
plt.xlabel("t [s]")
plt.ylabel(r"$\Psi$")
fs.post(fig, figSize, plt.legend())
plt.savefig("psi.png",dpi=1000)

fig = plt.figure("eta")
plt.plot(times, eta, label="MUSCLE-Foam")
plt.xlabel("t [s]")
plt.ylabel(r"$\eta$")
fs.post(fig, figSize, plt.legend())
plt.savefig("eta.png",dpi=1000)
