# -*- coding: utf-8 -*-
"""
Vertical LPCVD Polysilicon Furnace Uniformity Model
===================================================

This script models axial and radial thickness uniformity for a vertical
batch LPCVD polysilicon reactor with multiple injectors and heating zones.

The model combines:
- Gas-phase precursor depletion along the wafer stack
- Surface reaction kinetics for silane decomposition
- Radial diffusion limitations at the wafer surface
- Time-integration to a target thickness
- Axial, radial, and combined uniformity metrics (R/M)
- A wafer rotation criterion to justify azimuthal averaging

Primary outputs:
- Axial non-uniformity
- Radial non-uniformity
- Combined non-uniformity (root-sum-square)
- Average product wafer thickness
- Deposition time
- Minimum wafer rotation speed required for azimuthal averaging
- Average silane composition per zone
- Fractional conversion through each zone

Author: Rebecca Rhodes
"""

import numpy as np
import matplotlib.pyplot as plt

#Constants and Geometry

T_inlet = 873.0 #K, 600C
X_inlet = 0.3 #silane mole fraction

A1 = 7.076e6 #Back calculated from  Kimihisa Kaneko et al
A2 = 4.0e20
Ea1 = 152900 #J/mol (1.529e8 J/kmol / 1000)
Ea2 = 186000 #J/mol (1.86e8 J/kmol / 1000)

R_kin = 8.314 #J/molK (used in rate constants only)
R_gas = 8314 #J/kmolK (used in ideal gas law only)

P_standard = 101325 #For flow conversions
T_standard = 293 #for Flow conversions

Si_atomic_density = 5e28 #atoms/m^3
Na = 6.022e23 #atoms/mol
P = 26.6 #Pascals, 200 mTorr

h_target = 1e-7  #m, 1000 angstroms 
uniformity_target = 0.03 #3% R/M
h_min = h_target * (1-uniformity_target)  #m, 970 angstroms 

furnace_diameter = 0.423 #m
furnace_height   = 1.75 #m
wafer_diameter   = 0.300 #m, 300mm

wall_reactivity = 0.2 #guesstimate based on SiC reactiity at 600C

single_wafer_area = 2 * np.pi * (wafer_diameter/2)**2 #m^2
wall_area = 2.07 * np.pi * furnace_diameter * furnace_height #m^2
end_area = 2 * np.pi * (furnace_diameter/2)**2 #m^2

#Wafers
N_dummy_bot = 5
N_dummy_top = 5
N_product = 150

N_tot_wafers = N_product + N_dummy_bot + N_dummy_top
product_slice = slice(N_dummy_bot, N_tot_wafers - N_dummy_top)

wall_area_per_wafer = (wall_area + end_area) / N_tot_wafers

#Radial mass transfer 
    #Represent two effective diffusion path lengths across the wafer surface.
    #These approximate center vs edge depletion without explicit radial fluid dynamics.

delta_r = np.array([1e-4, 5e-4])

##############################################################################
#Functions, functions, and more functions

#Kinetics!

#Two-step surface reaction model for silane decomposition
# SiH4(g) + 2Si(s) ->  2SiH(adsorbed) + Si(bulk)     (rate constant k1)
# 2SiH(adsorbed) ->  2Si(s) + H2(g)  (rate constant k2)
#From  Kimihisa Kaneko et al


#Arrhenius expressions for surface reaction rate constants
    # k = A * exp(-Ea / RT)
def rate_constants(T):
    k1 = A1 * np.exp(-Ea1 / (R_kin * T))
    k2 = A2 * np.exp(-Ea2 / (R_kin * T))
    return k1, k2


#Ideal-gas concentration of silane
    # C = (P * x) / (R * T
def silane_concentration(T, x):
    return (P * x) / (R_gas * T)


#Surface-reaction-limited silicon deposition rate (m/s).
    #Uses a coverage-dependent formulation:
        # theta_s = sqrt(k2) / (sqrt(k1*C) + sqrt(k2))
        # r   = k1 * C * theta_s^2
    #The result is converted from molar to thickness rate.
def deposition_rate(C, T):
    k1, k2 = rate_constants(T)
    theta_s = np.sqrt(k2) / (np.sqrt(k1 * C) + np.sqrt(k2))
    r = k1 * C * theta_s**2
    Si_molar_density = Si_atomic_density / (Na * 1000)
    return r / Si_molar_density


#Diffusion-limited surface deposition rate
    #Solves a 1D steady boundary-layer balance:
        # C_s = C_bulk - (R * rho_Si * delta / D)
        #where C_s is the surface concentration of silane, C_bulk is the bulk
        #concentration silane, R is the deposition rate, rho_Si is density of silane, 
        #delta is the boundary layer thickness, and D is the diffusivity
    #Iterated until surface concentration convergence.
def deposition_rate_surface(C_bulk, T, delta, D):
    Si_molar_density = Si_atomic_density / (Na * 1000)
    C_s = C_bulk
    for _ in range(50):
        R = deposition_rate(C_s, T)
        C_new = max(C_bulk - R * Si_molar_density * delta / D, 0)
        if abs(C_new - C_s) < 1e-15:
            break
        C_s = C_new
    return deposition_rate(C_s, T)


#Wafer rotation should be fast enough to time-average azimuthal variations.
#We require N_rotations full revolutions during deposition:
   # omega ≥ (2pi * N_rotations) / t_process
#Typical vertical LPCVD furnaces rotate ~0.5–1 rpm, there are dimenishing returns
#above 1 rpm, as based on "Throughput Improvement in an LPCVD TEOS Vertical Furnace"
#by Ekbundit & Izzio
def required_rotation_rpm(t_process, N_rotations=60):
    #Minimum wafer rotation speed (RPM) to azimuthally average deposition.
    #N_rotations = number of full rotations during t_process
    
    omega = N_rotations * 2 * np.pi / t_process   #rad/s
    rpm = omega * 60 / (2 * np.pi) #rpm
    return rpm


#Uniformity function: Computes axial, radial, and combined thickness non-uniformity
#for a given furnace configuration. Also produces zonal conversion and average composition.
def uniformity_metrics(Q_sccm_total, N_zones, N_inj_per_zone):

    wafers_per_zone = N_tot_wafers // N_zones
    #Define injector coverage per zone
    injector_bounds = [np.linspace(z*wafers_per_zone, (z+1)*wafers_per_zone, 
                                   N_inj_per_zone + 1,dtype=int)
    
    for z in range(N_zones)]
    #Convert sccm to actual m^3/s at process conditions
    Q_total = (Q_sccm_total / 60e6 * (P_standard / P) * (T_inlet / T_standard))
   
    Q_zone = Q_total / N_zones
    Q_inj  = Q_zone / N_inj_per_zone

    #Binary diffusivity (diffusion coefficient) estimate
    D = 1e-5 * (T_inlet / 300)**1.75 * (101325 / P)

    C_profile = np.zeros(N_tot_wafers)
    C_carry = 0
    Q_carry = 0
    zone_C_sums = np.zeros(N_zones) #sum of C per zone
    zone_counts = np.zeros(N_zones) #wafer count per zone
    zone_consumed = np.zeros(N_zones) #silane consumed by deposition in zone
    zone_supplied = np.zeros(N_zones) #total silane supplied (carry-in + injected)
    
    for z in range(N_zones):
       
        for i in range(N_inj_per_zone):

            w0, w1 = injector_bounds[z][i:i+2]
            C_fresh = silane_concentration(T_inlet, X_inlet)
            Q_after = Q_carry + Q_inj
            C = (C_carry * Q_carry + C_fresh * Q_inj) / Q_after #perfect mixing assumption
            zone_supplied[z] += C_fresh * Q_inj   #fresh silane added this injector
            
            for w in range(w0, w1):
                C_profile[w] = C
                rates = [deposition_rate_surface(C, T_inlet, d, D) for d in delta_r]
              
                zone_C_sums[z] += C
                zone_counts[z] += 1
                R_mean = np.mean(rates)
                Si_molar = Si_atomic_density / (Na * 1000)
                consumption = R_mean * Si_molar * (single_wafer_area + wall_reactivity * wall_area_per_wafer)
                zone_consumed[z] += consumption    #accumulate molar consumption
                
                C = max(C - consumption / Q_after, 0)
                
            C_carry = C
            Q_carry = Q_after
            
    
    #Convert average C per zone back to mole fraction: x = C * R_gas * T / P
    C_avg_per_zone = zone_C_sums / zone_counts
    x_eff_per_zone = C_avg_per_zone * R_gas * T_inlet / P   #effective mole fraction
    
    #Add the carry-in contribution to supplied for zone 1
    #(zones 2+ receive depleted carry gas, already accounted for as leftover from prior zone)
    zone_supplied[0] += silane_concentration(T_inlet, 0) * 0  # zone 1 has no carry-in

    conv_per_zone = np.where(zone_supplied > 0, zone_consumed / zone_supplied, np.nan)
        
    #Radial and axial rates
    rates_radial = np.array([[deposition_rate_surface(C, T_inlet, d, D)
                             for d in delta_r]
                             for C in C_profile])

    rates_mean = rates_radial.mean(axis=1)

    r_mean_prod = rates_mean[product_slice]
    r_radial_prod = rates_radial[product_slice]

    t_process = h_target / np.mean(r_mean_prod) #min, time required to reach desired thickness
    
    rpm_min = required_rotation_rpm(t_process)
    
    #Integrated thickness profiles
    h_axial = r_mean_prod * t_process
    h_radial = r_radial_prod * t_process

    axial_u = (h_axial.max() - h_axial.min()) / h_axial.mean()
    radial_u = ((h_radial.max(axis=1) - h_radial.min(axis=1)) / h_axial).max()

    #Root-Sum-Square combination
    combined_u = np.sqrt(axial_u**2 + radial_u**2)

    h_average_prod = h_axial.mean()   #m, average product thickness
    
    return axial_u, radial_u, combined_u, h_average_prod, t_process, rpm_min, x_eff_per_zone, conv_per_zone

##############################################################################
#Optimization (parametric sweep)

Q_range = [250, 300, 350] #range of total flowrates (sccm)
zone_range = [4, 5, 6] #range of number of zones
N_inj_range  = [10, 20, 30, 40] #range of number of injectors per zone

results = []

for Q in Q_range:
    for Nz in zone_range:
        for Nin in N_inj_range:
            ax, rad, comb, havg, tproc, rpm_min, x_zones, conv = uniformity_metrics(Q, Nz, Nin)
            results.append((Q, Nz, Nin, Nz*Nin, ax, rad, comb, havg, tproc, rpm_min, x_zones, conv))

##############################################################################
#Results!!

results_sorted = sorted(results, key=lambda r: r[6]) #Sort results best uniformity to worst

Qv = np.array([r[0] for r in results_sorted])
Ntot = np.array([r[3] for r in results_sorted])
Comb = np.array([r[6] for r in results_sorted])
Havg = np.array([r[7] for r in results_sorted])
Tproc = np.array([r[8] for r in results_sorted])
RPM = np.array([r[9] for r in results_sorted])

#Valid = good uniformity AND minimumum thickness of 970 angstoms
valid = (Comb <= uniformity_target) & (Havg >= h_min)

print("\nConfigurations:")
for r in results_sorted:
    X_eff = "  ".join(f"Z{i+1}:{v:.5f}" for i, v in enumerate(r[10])) #adjust decimal as desired, often higher zones have similar comps
    zone_conv = "  ".join(f"Z{i+1}:{v:.3f}"  for i, v in enumerate(r[11]))
    print(f"Q={r[0]:4.0f}  zones={r[1]:2d}  inj/zone={r[2]:2d}  "
          f"Ax={r[4]:.5f}  Rad={r[5]:.15f}  Comb={r[6]:.5f}  "
          f"h_avg={r[7]*1e10:.1f} Å  "f"t_proc={r[8]/60:.1f} min  "
          f"RPM = {r[9]:.2f}  Zone compositions: {X_eff}   Zone Conversion: {zone_conv}")
    
plt.figure(figsize=(9,6))

#Plot all configurations
sc = plt.scatter(Qv, Ntot, c=Comb, cmap="plasma", s=90)

#Red circle overlay for valid designs
plt.scatter(Qv[valid], Ntot[valid], facecolors="none", 
            edgecolors="red", s=220, linewidths=2, label="≤3% R/M")


plt.colorbar(sc, label="Combined thickness non-uniformity")
plt.xlabel("Total Flow (sccm)")
plt.ylabel("Total Injectors")
plt.xticks(np.arange(250, 351, 50)) 
plt.yticks(np.arange(40, 245, 30))
plt.title(f"Designs with desired product uniformity highlighted, silane composition = {X_inlet}")
plt.legend(bbox_to_anchor=(0.85, 0.6))
plt.tight_layout()
plt.show()
