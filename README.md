Vertical LPCVD Polysilicon Furnace Uniformity Model

A simulation for optimizing axial and radial thickness uniformity in a vertical batch LPCVD (Low-Pressure Chemical Vapor Deposition) polysilicon reactor. The model supports multi-injector, multi-zone furnace configurations and performs parametric sweeps to identify process conditions that meet a target film thickness uniformity.

Author: Rebecca Rhodes, Fab Four Buffs Senior Design Team

Debugged using Claude.AI

Made for CHEN 4530 at CU Boulder

Overview:
In a vertical batch LPCVD furnace, silane (SiH₄) decomposes on heated wafer surfaces to deposit a polysilicon thin film. Two competing non-uniformity mechanisms must be controlled:
- Axial non-uniformity: Silane is consumed as gas flows up the wafer stack, depleting the precursor concentration and reducing the deposition rate at wafers further from the injectors.
- Radial non-uniformity: Mass transfer limitations create a concentration gradient from the wafer center to its edge.
This model integrates both effects to compute a combined thickness non-uniformity and identify furnace configurations (flow rate, number of heating zones, injectors per zone) that achieve a target uniformity of ≤ 3% (range/mean).

--- 

Physical Model:

Reactor Geometry: The furnace is modeled as a vertical cylinder containing a stack of 160 wafers (150 product wafers plus 5 dummy wafers at the top and bottom of the stack). The furnace dimensions were calculated separately using known wafer sizes and expected pitch between them. 

Key geometric parameters:
| Parameter | Value |
|---|---|
| Furnace diameter | 0.423 m |
| Furnace height | 1.75 m |
| Wafer diameter | 300 mm |
| Process pressure | 26.6 Pa (200 mTorr) |
| Process temperature | 873 K (600 °C) |
| Inlet silane mole fraction | 0.30 |
| Target film thickness | 1e-7 m (1000 Å) |
| Uniformity target | ≤ 3% range/mean |

The effective wall area per wafer accounts for deposition on the furnace tube walls (scaled by a wall reactivity factor “wall_reactivity = 0.2”, representing reduced SiC reactivity at 600 °C).

---

Surface Reaction Kinetics: Silane decomposition is modeled as a two-step surface reaction:

SiH4(g) + 2Si(s) ->  2SiH(adsorbed) + Si(bulk) (rate constant k1)

2SiH(adsorbed) ->  2Si(s) + H2(g)  (rate constant k2)

Both rate constants follow an Arrhenius expression:

$$k_i = A_i \cdot \exp\!\left(\frac{-E_{a,i}}{RT}\right)$$

| Constant | Pre-exponential factor | Activation energy |
|---|---|---|
| k₁ | 7.076 × 10⁶ | 152,900 J/mol |
| k₂ | 4.0 × 10²⁰ | 186,000 J/mol |

Activation energies and Pre-exponential factor A2 were taken directly from Saxena et al, with A1 being back calculated based on reactor geometry

Function: “rate_constants(T)” — returns k₁ and k₂ at temperature T (K).

---

Silane Concentration: The molar concentration of silane in the gas phase is calculated from the ideal gas law at low pressure:

$$C = \frac{P \cdot x}{R \cdot T}$$

where P is the total pressure (Pa), x is the silane mole fraction, R is the gas constant (8314 J/kmol·K), and T is the temperature (K).

Function: silane_concentration(T, x) — returns concentration in kmol/m³.

---

Deposition Rate (Reaction-Limited): A surface coverage model is used to account for competitive adsorption between silane and hydrogen intermediates. The steady-state surface coverage of reactive sites is:

$$\theta_s = \frac{\sqrt{k_2}}{\sqrt{k_1 C} + \sqrt{k_2}}$$

The molar deposition rate is then:

$$r = k_1 \cdot C \cdot \theta_s^2$$

This is converted from molar flux (kmol/m²·s) to a linear thickness rate (m/s) using the silicon molar density:

$$\rho_{Si,\text{molar}} = \frac{n_{Si}}{N_A \cdot 1000}$$

where n_Si = 5 × 10²⁸ atoms/m³ is the atomic number density of silicon and N_A is Avogadro's number.

Function: “deposition_rate(C, T)” — returns the surface-reaction-limited thickness deposition rate (m/s).

--- 

Deposition Rate (Diffusion-Limited): At the wafer surface, a boundary layer can limit mass transport of silane from the bulk gas to the reactive surface. The model solves a 1D steady-state boundary layer balance iteratively:

$$C_s = C_{\text{bulk}} - \frac{r \cdot \rho_{Si,\text{molar}} \cdot \delta}{D}$$

where C_s is the surface silane concentration, δ is the effective diffusion path length (boundary layer thickness), and D is the binary diffusivity. The surface concentration is iterated to convergence (up to 50 steps, tolerance 10⁻¹⁵ kmol/m³), then the final deposition rate is evaluated at C_s.

Two representative path lengths “delta_r = [1e-4, 5e-4]” m represent center and edge radial positions, providing an estimate of radial non-uniformity without explicit 2D fluid dynamics.

The binary diffusivity D is estimated from a standard Chapman-Enskog-type scaling:

$$D = D_0 \cdot \left(\frac{T}{300}\right)^{1.75} \cdot \frac{101325}{P}$$

with D₀ = 10⁻⁵ m²/s as the reference diffusivity.

Function: “deposition_rate_surface(C_bulk, T, delta, D)” — returns the diffusion-corrected thickness deposition rate (m/s).

--- 

Axial Precursor Depletion: Gas flows upward through the furnace. Within each zone, each injector introduces a fresh stream of silane at the inlet composition. Between injectors, the silane concentration depletes as it is consumed by deposition on wafers and the furnace walls.

For each wafer segment between injector positions, the local mixed concentration after fresh injection is:

$$C_{\text{mixed}} = \frac{C_{\text{carry}} \cdot Q_{\text{carry}} + C_{\text{fresh}} \cdot Q_{\text{inj}}}{Q_{\text{carry}} + Q_{\text{inj}}}$$

This assumes perfect radial mixing within each injector segment — a decent approximation at low pressure where diffusive timescales are short. Silane is then consumed wafer-by-wafer:

$$C_{\text{new}} = \max\!\left(C - \frac{R_{\text{mean}} \cdot \rho_{Si,\text{molar}} \cdot (A_{\text{wafer}} + f_w \cdot A_{\text{wall}})}{Q_{\text{after}}},\ 0\right)$$

where f_w is the wall reactivity factor.

--- 

Uniformity Metrics: Film thickness at each product wafer is computed as:

$$h_i = \bar{r}_i \cdot t_{\text{process}}$$

where the mean rate $\bar{r}_i$ averages the center and edge deposition rates for wafer i, and the process time t_process is set to reach the target thickness at the average rate across all product wafers.

Axial non-uniformity (range/mean across the wafer stack):

$$U_{\text{axial}} = \frac{h_{\max} - h_{\min}}{\bar{h}}$$

Radial non-uniformity (maximum center-to-edge variation, normalized):

$$U_{\text{radial}} = \max_i\left(\frac{h_{i,\text{edge}} - h_{i,\text{center}}}{h_i}\right)$$

Combined non-uniformity (root-sum-square):

$$U_{\text{combined}} = \sqrt{U_{\text{axial}}^2 + U_{\text{radial}}^2}$$

A configuration is marked as valid if both U_combined ≤ 3% and the average product thickness h_avg ≥ 970 Å (i.e., within the 3% lower bound of the 1000 Å target). This is to prevent “valid” configurations that have low uniformity due to low or zero deposition.

Function: “uniformity_metrics(Q_sccm_total, N_zones, N_inj_per_zone) — returns all uniformity metrics, deposition time, minimum rotation speed, zone compositions, and fractional conversion per zone.

--- 

Wafer Rotation Criterion: In a real furnace, wafers rotate to average out azimuthal (angular) non-uniformities caused by gas jets from injectors. For azimuthal averaging to be valid, the wafer must complete a sufficient number of full revolutions during the deposition:

$$\omega_{\min} = \frac{2\pi \cdot N_{\text{rot}}}{t_{\text{process}}}$$

Converting to RPM:

$$\text{RPM}_{\min} = \frac{\omega_{\min} \cdot 60}{2\pi}$$

The default requirement is N_rot = 60 full revolutions during deposition. Based on Ekbundit & Izzio, practical furnaces operate at 0.5–1 RPM, with diminishing returns above 1 RPM. The model reports the minimum RPM required to justify the azimuthal averaging assumption.

Function: “required_rotation_rpm(t_process, N_rotations=60) — returns minimum rotation speed in RPM.

Functions Summary: 

| Function | Description |
|---|---|
| `rate_constants(T)` | Arrhenius rate constants k₁ and k₂ at temperature T |
| `silane_concentration(T, x)` | Ideal-gas silane molar concentration |
| `deposition_rate(C, T)` | Surface-reaction-limited thickness rate |
| `deposition_rate_surface(C_bulk, T, delta, D)` | Diffusion-corrected deposition rate at a given boundary layer thickness |
| `required_rotation_rpm(t_process, N_rotations)` | Minimum wafer rotation speed for azimuthal averaging |
| `uniformity_metrics(Q, N_zones, N_inj)` | Full model evaluation for a given furnace configuration |

---

Parametric Sweep: The script sweeps over three design variables:

| Variable | Values Swept |
|---|---|
| Total flow rate Q (sccm) | 250, 300, 350 |
| Number of heating zones | 4, 5, 6 |
| Injectors per zone | 10, 20, 30, 40 |

This produces 3 × 3 × 4 = 36 configurations, each evaluated for uniformity, thickness, and deposition time. Results are sorted by combined uniformity (best to worst) and printed to the console. Valid designs are highlighted in the output scatter plot.

Outputs:
Console: For each configuration, the script prints:
- Total flow rate, number of zones, injectors per zone
- Axial, radial, and combined non-uniformity
- Average product wafer thickness (Å)
- Deposition time (min)
- Minimum rotation speed (RPM)
- Average silane mole fraction per zone
- Fractional silane conversion per zone

Plot: Scatter plot of total injector count vs. total flow rate, colored by combined non-uniformity. Valid designs (≤ 3% R/M with sufficient thickness) are circled in red.

Literature Citations
(1) Kaneko, K.; Ogino, M.; Shimizu, R.; Koshi, M.; Shimogaki, Y. Film Thickness Prediction of Poly-Silicon LPCVD Process with a Simplified Two-Step Surface Reaction Model. ECS Journal of Solid State Science and Technology 2013, 2 (9), N182–N186. https://doi.org/10.1149/2.024309jss.

(2) Ekbundit, S.; Izzio, B. Characterization of Film Uniformity in LPCVD TEOS Vertical Furnace. 13th Annual IEEE/SEMI Advanced Semiconductor Manufacturing Conference. Advancing the Science and Technology of Semiconductor Manufacturing. ASMC 2002 (Cat. No.02CH37259) 2002, 38–42. https://doi.org/10.1109/asmc.2002.1001570.
