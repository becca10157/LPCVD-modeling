# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 11:31:05 2026

@author: rebec
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

#constants
A1 = 7.076e+06 #back calculated from paper
A2 = 4e20
Ea1 = 152900 #J/mol (1.529e8 J/kmol / 1000)
Ea2 = 186000 #J/mol (1.86e8 J/kmol / 1000)
R_kin = 8.314 #J/molK (used in rate constants only)
R_gas = 8314 #J/kmolK (used in ideal gas law only)
Si_atomic_density = 5e28 #atoms/m^3
Na = 6.022e23 #atoms/mol
N_wafers = 150
P = 26.6 #pascals

#Q at process conditions
P_standard = 101325 #Pa (1 atm)
T_standard = 293 #K (20°C)
Q_sccm = 350 #sccm
T_process = 873 #K
Q = Q_sccm / 60 / 1e6 * (P_standard / P) * (T_process / T_standard)
print(f"Q at process conditions: {Q:.3e} m³/s")

h_target = 1e-7  #m

furnace_diameter = 0.423 #m
furnace_height = 1.25 #m
wafer_diameter = 0.3 #m

wafer_area = N_wafers * np.pi * (wafer_diameter / 2) ** 2
wall_area = np.pi * furnace_diameter * furnace_height
end_area = 2 * np.pi * (furnace_diameter / 2) ** 2
total_area = wafer_area + wall_area + end_area
single_wafer_area = np.pi * (wafer_diameter / 2) ** 2

def rate_constants(T):
    k1 = A1 * np.exp(-Ea1 / (R_kin * T))
    k2 = A2 * np.exp(-Ea2 / (R_kin * T))
    return k1, k2

def silane_concentration(T, x):
    return (P * x) / (R_gas * T) #kmol/m³

def deposition_rate(C, T):
    k1, k2 = rate_constants(T)
    theta_s = np.sqrt(k2) / (np.sqrt(k1 * C) + np.sqrt(k2))
    r1 = k1 * C * theta_s**2
    Si_molar_density = Si_atomic_density / (Na * 1000) #kmol/m³
    R_film = r1 / Si_molar_density #m/s
    return R_film

def concentration_profile(T_wafers, x_input):
    # x_input can be a scalar (fixed x) or array (x varies across stack with dx)
    x_wafers = np.full(len(T_wafers), x_input) if np.isscalar(x_input) else np.asarray(x_input)
    C = silane_concentration(T_wafers[0], x_wafers[0])
    concentrations = []
    Si_molar_density = Si_atomic_density / (Na * 1000)
    for f, T in enumerate(T_wafers):
        concentrations.append(C)
        R_dep = deposition_rate(C, T)
        consumption = R_dep * Si_molar_density * single_wafer_area  # kmol/s per wafer
        C = max(C - consumption / Q, 0)
    return np.array(concentrations)


#C_drop diagnostid
T_test = 873
x_test = 0.5
C_test = silane_concentration(T_test, x_test)
R_dep_test = deposition_rate(C_test, T_test)
Si_molar_density_test = Si_atomic_density / (Na * 1000)
consumption_test = R_dep_test * Si_molar_density_test * single_wafer_area
C_drop_test = consumption_test / Q
print(f"C: {C_test:.3e} kmol/m³")
print(f"R_dep: {R_dep_test:.3e} m/s")
print(f"C_drop: {C_drop_test:.3e} kmol/m³")
print(f"C_drop/C ratio: {C_drop_test/C_test:.3f}")

#Optimization ranges
T_range  = np.linspace(870, 876, 15)
x_range  = np.linspace(0.01, 0.99, 15)
dT_range = np.linspace(0.1, 15.0, 15)
dx_range = np.linspace(0.001, 0.05, 15)

#some times it gets mad without these
best_rate = -1
best_T = None
best_x = None
best_dT = None
best_dx = None
best_uniformity = None

total = len(dT_range) * len(dx_range) * len(T_range) * len(x_range)
count = 0

for dT in dT_range:
    for dx in dx_range:
        for T_avg in T_range:
            T_wafers = np.linspace(T_avg - dT, T_avg + dT, N_wafers)
            for x_avg in x_range:
                C_profile = concentration_profile(T_wafers, x_avg)
                rates = np.array([deposition_rate(C_profile[f], T)
                for f, T in enumerate(T_wafers)])
                R_mean = np.mean(rates)
                uniformity = (np.max(rates) - np.min(rates)) / R_mean
                if uniformity <= 0.03 and R_mean > best_rate:
                    best_rate = R_mean
                    best_T = T_avg
                    best_x = x_avg
                    best_dT = dT
                    best_dx = dx
                    best_uniformity = uniformity
                count += 1
        print(f"\rSearch progress: {100*count/total:.1f}%", end="") #im impatient
print()

#Best achievable uniformity
min_u = float('inf')
min_u_T = min_u_x = min_u_dT = None
for dT in dT_range:
    for T_avg in T_range:
        T_wafers = np.linspace(T_avg - dT, T_avg + dT, N_wafers)
        for x_avg in x_range:
            C_profile = concentration_profile(T_wafers, x_avg)
            rates = np.array([deposition_rate(C_profile[f], T)
            for f, T in enumerate(T_wafers)])
            u = (np.max(rates) - np.min(rates)) / np.mean(rates)
            if u < min_u:
                min_u = u
                min_u_T, min_u_x, min_u_dT = T_avg, x_avg, dT

print(f"Best uniformity: {min_u:.4f} at T={min_u_T:.2f}, x={min_u_x:.5f}, dT={min_u_dT:.3f}")

if min_u > 0.03:
    print(f"Best achievable uniformity ({min_u:.4f}) exceeds 3% target. Adjust ranges.")
    sys.exit()

#Unit convert
rate_m_per_s = best_rate
rate_A_per_min = best_rate * 1e10 * 60

#Contour map
T_grid = np.linspace(843, 893, 50)
x_grid = np.linspace(0.01, 0.99, 50)
rate_map = np.zeros((len(x_grid), len(T_grid)))
uniformity_map = np.zeros((len(x_grid), len(T_grid)))

for i, x_avg in enumerate(x_grid):
    for j, T_avg in enumerate(T_grid):
        T_wafers = np.linspace(T_avg - best_dT, T_avg + best_dT, N_wafers)
        C_profile = concentration_profile(T_wafers, x_avg)
        rates = np.array([deposition_rate(C_profile[f], T)
        for f, T in enumerate(T_wafers)])
        rate_map[i, j] = np.mean(rates) * 1e10 * 60
        uniformity_map[i, j] = (np.max(rates) - np.min(rates)) / np.mean(rates)

plt.figure()
contours = plt.contourf(T_grid, x_grid, rate_map, 20)
plt.colorbar(contours, label="Deposition Rate (Å/min)")
#plt.contour(T_grid, x_grid, uniformity_map, [0.03], colors='red', linewidths=2)
plt.scatter(best_T, best_x, color='white', edgecolor='black', label="Optimal Point")
plt.xlabel("Temperature (K)")
plt.ylabel("Silane Mole Fraction")
plt.title(f"Deposition Rate Map with 3% Uniformity Boundary (dT=±{best_dT:.2f}K)")
plt.legend()
plt.show()

if best_T is None:
    print("No condition found meeting uniformity target.")
else:
    print("Optimal Reactor Conditions:")
    print("--------------------------")
    print(f"Temperature: {best_T:.2f} K  (±{best_dT:.3f} K across stack)")
    print(f"Silane Mole Fraction: {best_x:.5f}  (±{best_dx:.4f} across stack)")
    print(f"Uniformity (R/M): {best_uniformity:.4f}")
    print("\nDeposition Rate")
    print(f"{rate_m_per_s:.3e} m/s")
    print(f"{rate_A_per_min:.2f} Å/min")
    time_to_target = h_target / best_rate
    print(f"\nTime to grow {h_target*1e10:.0f} Å: {time_to_target/60:.2f} minutes")