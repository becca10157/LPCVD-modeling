# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 21:33:00 2026

@author: rebec
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

#Constants and params
A1 = 7.076e+06 #back calculated
A2 = 4e20
Ea1 = 152900
Ea2 = 186000
R_kin = 8.314 #J/molK  for arrhenius
R_gas = 8314 #J/kmolK for gas law
Si_atomic_density = 5e28
Na = 6.022e23
P = 26.6
""" Capital X = conversion, lowercase x = mol frac """ #should prob pick a different letter for one of these

N_dummy_bottom = 5
N_dummy_top = 5
N_product = 150
N_wafers = N_product + N_dummy_bottom + N_dummy_top #total number wafers
N_zones = 5
N_injectors_per_zone = 30 #150 total
N_injectors_tot = N_zones * N_injectors_per_zone             
wafers_per_zone = N_wafers // N_zones                        

product_slice = slice(N_dummy_bottom, N_wafers - N_dummy_top) #exclude dummies from uniformity calc

wall_reactivity = 0.20 #guessitamte from Kira

#Radial params
delta_center = 8e-3 #boundary layer @ center
delta_edge = 3e-3 #boundary layer @ edge
dT_radial = 1.0 #temp diff between center and edge: guess and check/guesstimate (accurate??)
delta_r = np.array([delta_center, delta_edge])
N_radial = 2 #dont be slow


#Injector boundary indices: injector_boundaries[z][i]::[i+1] = wafer slice
#dont recalc geometry every loop later
injector_boundaries = []
for z in range(N_zones):
    zone_start = z * wafers_per_zone
    zone_end = zone_start + wafers_per_zone
    bounds = np.linspace(zone_start, zone_end, N_injectors_per_zone + 1, dtype=int)
    injector_boundaries.append(bounds)


#Flow at process conditions 
P_standard = 101325
T_standard = 293
Q_sccm_total = 350
T_process = 873

#Flow rate conversions
Q_total    = Q_sccm_total / 60 / 1e6 * (P_standard / P) * (T_process / T_standard)
Q_zone     = Q_total / N_zones
Q_injector = Q_zone / N_injectors_per_zone

D_SiH4 = 1e-5 * (T_process / 300)**1.75 * (101325 / P) #estimate binary diffusivity of silane (Chapman-Enskog scaling)

print(f"Q total          : {Q_total:.3e} m³/s")
print(f"Q per zone       : {Q_zone:.3e}  m³/s")
print(f"Q per injector   : {Q_injector:.3e} m³/s  ({N_injectors_tot} injectors total)")
print(f"D_SiH4           : {D_SiH4:.3e} m²/s")

#Reactor geometry
h_target = 1e-7
furnace_diameter = 0.423
furnace_height = 1.75
wafer_diameter = 0.300

wall_area = np.pi * furnace_diameter * furnace_height
end_area = 2 * np.pi * (furnace_diameter / 2)**2
wall_area_per_slot = (wall_area + end_area) * 1.07 / N_wafers
single_wafer_area = np.pi * (wafer_diameter / 2)**2 #one side, need to be two? does boat cover bottom?

print(f"\nGeometry:")
print(f"  Wafers per zone      : {wafers_per_zone}")
print(f"  Wafers per injector  : ~{wafers_per_zone / N_injectors_per_zone:.1f}")
print(f"  Injector boundaries  (zone 0): {injector_boundaries[0]}")


#Kinetics!
def rate_constants(T):
    k1 = A1 * np.exp(-Ea1 / (R_kin * T))
    k2 = A2 * np.exp(-Ea2 / (R_kin * T))
    return k1, k2

def silane_concentration(T, x):
    return (P * x) / (R_gas * T)

def deposition_rate(C, T):
    k1, k2 = rate_constants(T)
    theta_s = np.sqrt(k2) / (np.sqrt(k1 * C) + np.sqrt(k2))
    r1 = k1 * C * theta_s**2
    Si_molar_density = Si_atomic_density / (Na * 1000)
    return r1 / Si_molar_density

def deposition_rate_surface(C_bulk, T, delta): #coupled mass transport/surface reaction
    Si_molar_density = Si_atomic_density / (Na * 1000)
    C_s = C_bulk
    for _ in range(50):
        R = deposition_rate(C_s, T)
        C_s_new = max(C_bulk - R * Si_molar_density * delta / D_SiH4, 0) #SS flux balance, max so not neg
        if abs(C_s_new - C_s) < 1e-15: #tolerance
            break
        C_s = C_s_new
    return deposition_rate(C_s, T), C_s

def radial_rates(C_bulk, T_bulk, dT_rad=dT_radial):
    T_radial = np.array([T_bulk + dT_rad / 2, T_bulk - dT_rad / 2])
    rates = np.empty(N_radial)
    for i in range(N_radial):
        rates[i], _ = deposition_rate_surface(C_bulk, T_radial[i], delta_r[i])
    return rates

#Multi-zone concentration profile 
def concentration_profile_5zone(T_wafers, x_zones):
    #Each injector adds Q_injector of fresh silane at x_zones[z]
    
    Si_molar_density = Si_atomic_density / (Na * 1000)
    concentrations = np.empty(N_wafers)
    C_zone_inlets = np.empty(N_zones) #for diagnostics

    C_carry = 0
    Q_carry = 0

    for z in range(N_zones):
        bounds = injector_boundaries[z]

        for inj in range(N_injectors_per_zone):
            T_inj   = T_wafers[bounds[inj]]
            C_fresh = silane_concentration(T_inj, x_zones[z])
            Q_after = Q_carry + Q_injector
            C_mix   = (C_carry * Q_carry + C_fresh * Q_injector) / Q_after

            if inj == 0:
                C_zone_inlets[z] = C_mix

            C = C_mix
            for w in range(bounds[inj], bounds[inj + 1]):
                concentrations[w] = C
                r_rad = radial_rates(C, T_wafers[w])
                R_mean_w = np.mean(r_rad)
                consumption = R_mean_w * Si_molar_density * (single_wafer_area + wall_reactivity * wall_area_per_slot)
                C = max(C - consumption / Q_after, 0)

            C_carry = C
            Q_carry = Q_after

    return concentrations, C_zone_inlets

#Uniformity metrics
def compute_all_rates(C_profile, T_wafers): #recompute rates with new conc profile
    rates_radial = np.array([radial_rates(C_profile[w], T_wafers[w])
    for w in range(N_wafers)])
    rates_mean = rates_radial.mean(axis=1)
    return rates_radial, rates_mean

def uniformity_metrics(rates_radial, rates_mean):
    r_mean_prod   = rates_mean[product_slice]
    r_radial_prod = rates_radial[product_slice]

    axial_u = (np.max(r_mean_prod) - np.min(r_mean_prod)) / np.mean(r_mean_prod) #worst case wafer to wafer
    radial_u_per_wafer = ((np.max(r_radial_prod, axis=1) - np.min(r_radial_prod, axis=1)) / r_mean_prod) #calc each wafer
    radial_u = np.max(radial_u_per_wafer) #worst case individual wafer
    combined = max(axial_u, radial_u)
    return axial_u, radial_u, combined

#Fractional conversion diagnostic
def conversion_diagnostic(T_wafers, x_zones):
    Si_molar_density = Si_atomic_density / (Na * 1000)
    C_carry = 0
    Q_carry = 0
    total_consumed_wafers = total_consumed_walls = total_fed = 0

    print("\nFractional Conversion Diagnostic:")
    print(f"  {'Zone':<6} {'Inj':<5} {'C_inlet':>14}  {'C_exit':>14}  {'X_inj':>8}  {'X_cumul':>8}")
    print("  " + "-" * 62)

    for z in range(N_zones):
        bounds = injector_boundaries[z]

        for inj in range(N_injectors_per_zone):
            T_inj = T_wafers[bounds[inj]]
            C_fresh = silane_concentration(T_inj, x_zones[z])
            Q_after = Q_carry + Q_injector
            C_mix = (C_carry * Q_carry + C_fresh * Q_injector) / Q_after
            total_fed += C_fresh * Q_injector

            C = C_mix
            cw_seg = cwl_seg = 0

            for w in range(bounds[inj], bounds[inj + 1]):
                r_rad = radial_rates(C, T_wafers[w])
                R_w = np.mean(r_rad)
                cw = R_w * Si_molar_density * single_wafer_area
                cwl = R_w * Si_molar_density * wall_reactivity * wall_area_per_slot
                cw_seg += cw
                cwl_seg += cwl
                C = max(C - (cw + cwl) / Q_after, 0)

            total_consumed_wafers += cw_seg
            total_consumed_walls  += cwl_seg

            X_inj = (C_mix - C) / C_mix if C_mix > 0 else 0
            X_cumul = 1.0 - (C * Q_after) / total_fed if total_fed > 0 else 0
            print(f"  {z+1:<6} {inj+1:<5} {C_mix:>14.4e}  {C:>14.4e}  {X_inj:>8.4f}  {X_cumul:>8.4f}")

            C_carry = C
            Q_carry = Q_after

        print("  " + "·" * 62)

    total_consumed = total_consumed_wafers + total_consumed_walls
    X_overall = total_consumed / total_fed if total_fed > 0 else 0
    print(f"\n  Overall X          = {X_overall:.4f} ({X_overall*100:.2f}%)")
    print(f"  Consumed on wafers : {total_consumed_wafers/total_consumed*100:.1f}%")
    print(f"  Consumed on walls  : {total_consumed_walls/total_consumed*100:.1f}%")

#Diagnostics
T_diag = np.full(N_wafers, 873.0)
x_diag = np.full(N_zones, 0.5)
C_diag, C_inlets_diag = concentration_profile_5zone(T_diag, x_diag)
rates_rad_diag, rates_mean_diag = compute_all_rates(C_diag, T_diag)
ax_u, rad_u, comb_u = uniformity_metrics(rates_rad_diag, rates_mean_diag)

print(f"\nDiagnostic (x=0.5, T=873 K, uniform):")
for z in range(N_zones):
    print(f"  C inlet zone {z+1}: {C_inlets_diag[z]:.3e} kmol/m³")
print(f"  C exit         : {C_diag[-1]:.3e} kmol/m³")
print(f"\nUniformity (product wafers):")
print(f"  Axial          : {ax_u:.4f}")
print(f"  Radial (worst) : {rad_u:.4f}")
print(f"  Combined       : {comb_u:.4f}")
conversion_diagnostic(T_diag, x_diag)


#Optimizater!
T_range = np.linspace(870, 876, 10)
x_range = np.linspace(0.01, 0.99, 10)
dT_range = np.linspace(0.01, 10.0, 10)
dx_range = np.linspace(0.01, 0.10,  8)

best_rate = -1
best_T = best_x = best_dT = best_dx = best_uniformity = None
best_axial_u = best_radial_u = None

total = len(dT_range) * len(dx_range) * len(T_range) * len(x_range)
count = 0

for dT in dT_range:
    for dx in dx_range:
        for T_avg in T_range:
            T_wafers = np.linspace(T_avg - dT, T_avg + dT, N_wafers)
            for x_avg in x_range:
                x_zones = np.clip(np.linspace(x_avg - dx, x_avg + dx, N_zones), 1e-4, 1.0)
                C_profile, _ = concentration_profile_5zone(T_wafers, x_zones)
                rates_rad, rates_mean = compute_all_rates(C_profile, T_wafers)
                ax_u, rad_u, comb_u = uniformity_metrics(rates_rad, rates_mean)
                R_mean = np.mean(rates_mean[product_slice])

                if comb_u <= 0.03 and R_mean > best_rate:
                    best_rate = R_mean
                    best_T = T_avg
                    best_x = x_avg
                    best_dT = dT
                    best_dx = dx
                    best_uniformity = comb_u
                    best_axial_u = ax_u
                    best_radial_u = rad_u

                count += 1
        print(f"\rSearch progress: {100*count/total:.1f}%", end="") #im impatient but also this isnt accurate anymore, takes a few minutes after 100% reached
print()

#Best achievable uniformity (so far)
min_u = float('inf')
min_u_T = min_u_x = min_u_dT = min_u_dx = None
min_u_axial = min_u_radial = None

for dT in dT_range:
    for dx in dx_range:
        for T_avg in T_range:
            T_wafers = np.linspace(T_avg - dT, T_avg + dT, N_wafers)
            for x_avg in x_range:
                x_zones = np.clip(np.linspace(x_avg - dx, x_avg + dx, N_zones), 1e-4, 1.0)
                C_profile, _ = concentration_profile_5zone(T_wafers, x_zones)
                rates_rad, rates_mean = compute_all_rates(C_profile, T_wafers)
                ax_u, rad_u, comb_u  = uniformity_metrics(rates_rad, rates_mean)
                if comb_u < min_u:
                    min_u = comb_u
                    min_u_T, min_u_x, min_u_dT, min_u_dx = T_avg, x_avg, dT, dx
                    min_u_axial, min_u_radial = ax_u, rad_u

print(f"\nBest combined uniformity : {min_u:.4f}")
print(f"  Axial component        : {min_u_axial:.4f}")
print(f"  Radial component       : {min_u_radial:.4f}")
print(f"  At T={min_u_T:.2f} K, x={min_u_x:.4f}, dT=±{min_u_dT:.2f} K, dx={min_u_dx:.4f}")

if min_u > 0.03:
    print(f"\nBest achievable uniformity ({min_u:.4f}) exceeds 3% target.")
    print("Limiting factor: axial" if min_u_axial > min_u_radial else "Limiting factor: radial")
    sys.exit()

""" This is irrelevent currently...but theoretically helpful later

#Reconstruct optimal profile 
T_wafers_opt = np.linspace(best_T - best_dT, best_T + best_dT, N_wafers)
x_zones_opt = np.clip(np.linspace(best_x - best_dx, best_x + best_dx, N_zones), 1e-4, 1.0)
C_opt, C_inlets_opt = concentration_profile_5zone(T_wafers_opt, x_zones_opt)
rates_rad_opt, rates_mean_opt = compute_all_rates(C_opt, T_wafers_opt)

print("\nConversion diagnostic at optimal conditions:")
conversion_diagnostic(T_wafers_opt, x_zones_opt)

#Contour map
T_grid = np.linspace(843, 893, 40)
x_grid = np.linspace(0.01, 0.99, 40)
rate_map = np.zeros((len(x_grid), len(T_grid)))
uniformity_map = np.zeros((len(x_grid), len(T_grid)))

for i, x_avg in enumerate(x_grid):
    for j, T_avg in enumerate(T_grid):
        T_w = np.linspace(T_avg - best_dT, T_avg + best_dT, N_wafers)
        x_z = np.clip(np.linspace(x_avg - best_dx, x_avg + best_dx, N_zones), 1e-4, 1.0)
        C_p, _ = concentration_profile_5zone(T_w, x_z)
        r_rad, r_mean = compute_all_rates(C_p, T_w)
        _, _, comb = uniformity_metrics(r_rad, r_mean)
        rate_map[i, j] = np.mean(r_mean[product_slice]) * 1e10 * 60
        uniformity_map[i, j] = comb

#Plots 
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

ax = axes[0]
cf = ax.contourf(T_grid, x_grid, rate_map, 20)
fig.colorbar(cf, ax=ax, label="Mean Deposition Rate (Å/min)")
ax.contour(T_grid, x_grid, uniformity_map, [0.03], colors='red',
           linewidths=2, linestyles='--')
if best_T:
    ax.scatter(best_T, best_x, color='white', edgecolor='black', zorder=5, label="Optimal")
ax.set_xlabel("Mean Temperature (K)")
ax.set_ylabel("Mean Silane Mole Fraction")
ax.set_title(f"Deposition Rate Map\n({N_injectors_tot} injectors, 3% boundary in red)")
ax.legend(fontsize=8)

ax2 = axes[1]
w_prod = np.arange(N_dummy_bottom, N_wafers - N_dummy_top)
w_all = np.arange(N_wafers)
ax2.plot(w_prod, rates_rad_opt[product_slice, 0] * 1e10 * 60, color='steelblue', label="Center")
ax2.plot(w_prod, rates_rad_opt[product_slice, 1] * 1e10 * 60, color='tomato', label="Edge")
ax2.plot(w_all[:N_dummy_bottom], rates_rad_opt[:N_dummy_bottom, 0] * 1e10 * 60, color='steelblue', linestyle='--', alpha=0.4)
ax2.plot(w_all[:N_dummy_bottom], rates_rad_opt[:N_dummy_bottom, 1] * 1e10 * 60, color='tomato',    linestyle='--', alpha=0.4)
ax2.plot(w_all[-N_dummy_top:], rates_rad_opt[-N_dummy_top:, 0] * 1e10 * 60, color='steelblue', linestyle='--', alpha=0.4)
ax2.plot(w_all[-N_dummy_top:], rates_rad_opt[-N_dummy_top:, 1] * 1e10 * 60, color='tomato',    linestyle='--', alpha=0.4, label="Dummy")

#Mark injector positions
for z in range(N_zones):
    for inj in range(N_injectors_per_zone):
        w_inj = injector_boundaries[z][inj]
        ax2.axvline(w_inj, color='orange', linewidth=0.7, alpha=0.6, linestyle='--', label="Injector" if (z == 0 and inj == 0) else "")
ax2.axvspan(0, N_dummy_bottom - 0.5, alpha=0.08, color='gray')
ax2.axvspan(N_wafers - N_dummy_top - 0.5, N_wafers-1, alpha=0.08, color='gray')
ax2.set_xlabel("Wafer index (0 = bottom)")
ax2.set_ylabel("Deposition Rate (Å/min)")
ax2.set_title(f"Axial Profile — Center vs Edge\n({N_injectors_tot} injectors, dashed = dummy)")
ax2.legend(fontsize=7)

ax3 = axes[2]
rad_u_per_wafer = ((np.max(rates_rad_opt[product_slice], axis=1) - np.min(rates_rad_opt[product_slice], axis=1)) / rates_mean_opt[product_slice])
ax3.plot(w_prod, rad_u_per_wafer * 100, color='purple')
ax3.axhline(3.0, color='red', linestyle='--', linewidth=1, label="3% target")
ax3.set_xlabel("Wafer index (0 = bottom)")
ax3.set_ylabel("Radial non-uniformity (%)")
ax3.set_title("Per-wafer Radial Non-uniformity\n(center vs edge)")
ax3.legend(fontsize=8)

plt.tight_layout()
plt.show()

#Summary 
if best_T is None:
    print("\nNo condition found meeting combined uniformity target.")
    print("Limiting factor: axial" if min_u_axial > min_u_radial else "Limiting factor: radial")
else:
    rate_A_per_min = best_rate * 1e10 * 60
    print(f"\nOptimal Reactor Conditions ({N_injectors_tot} injectors total):")
    print("--------------------------------------------------------")
    print(f"Mean Temperature     : {best_T:.2f} K  (±{best_dT:.3f} K across stack)")
    print(f"Mean Silane Fraction : {best_x:.5f}  (gradient ±{best_dx:.4f} across zones)")
    print(f"Zone inlet x values  : {np.round(x_zones_opt, 5)}")
    print(f"\nUniformity (product wafers only):")
    print(f"  Axial              : {best_axial_u:.4f}")
    print(f"  Radial (worst)     : {best_radial_u:.4f}")
    print(f"  Combined           : {best_uniformity:.4f}")
    print(f"\nDeposition Rate      : {best_rate:.3e} m/s")
    print(f"                       {rate_A_per_min:.2f} Å/min")
    print(f"Time to grow {h_target*1e10:.0f} Å  : {h_target/best_rate/60:.2f} minutes")
"""