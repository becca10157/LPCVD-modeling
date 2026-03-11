# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 21:05:52 2026

@author: rebec
"""
from scipy.optimize import fsolve
import numpy as np

T = 873
C_paper = 4.133e-6   # kmol/m³ — pure SiH4 at 30 Pa
R_target = 2.5e-10   # m/s (~150 Å/min from Figure 6)

Ea1 = 152900
Ea2 = 186000
R_kin = 8.314
A2 = 4e20
Na = 6.022e23
Si_atomic_density = 5e28
Si_molar_density = Si_atomic_density / (Na * 1000)

k2 = A2 * np.exp(-Ea2 / (R_kin * T))

def equation(A1):
    k1 = A1 * np.exp(-Ea1 / (R_kin * T))
    theta_s = np.sqrt(k2) / (np.sqrt(k1 * C_paper) + np.sqrt(k2))
    r1 = k1 * C_paper * theta_s**2
    return r1 / Si_molar_density - R_target

A1_solution = fsolve(equation, 1e10)[0]
print(f"A1: {A1_solution:.3e}")

# verify
k1 = A1_solution * np.exp(-Ea1 / (R_kin * T))
theta_s = np.sqrt(k2) / (np.sqrt(k1 * C_paper) + np.sqrt(k2))
r1 = k1 * C_paper * theta_s**2
print(f"R_film at paper conditions: {r1/Si_molar_density:.3e} m/s")
print(f"R_film at paper conditions: {r1/Si_molar_density*1e10*60:.1f} Å/min")

# now check at YOUR conditions
C_yours = 3.665e-7
k1 = A1_solution * np.exp(-Ea1 / (R_kin * T))
theta_s = np.sqrt(k2) / (np.sqrt(k1 * C_yours) + np.sqrt(k2))
r1 = k1 * C_yours * theta_s**2
print(f"\nR_film at your conditions: {r1/Si_molar_density:.3e} m/s")
print(f"R_film at your conditions: {r1/Si_molar_density*1e10*60:.1f} Å/min")