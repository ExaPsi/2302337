# Companion Jupyter Notebooks

## Surface Phenomena & Catalysis (2302337)

Department of Chemistry, Faculty of Science, Chulalongkorn University

---

## Quick Start

### Requirements

| Package | Version | Purpose |
|---------|---------|---------|
| Python | >= 3.9 | Runtime |
| numpy | >= 1.21 | Array operations |
| scipy | >= 1.7 | Optimization, integration, statistics |
| matplotlib | >= 3.5 | Publication-quality plots |
| pandas | >= 1.3 | CSV data loading (NB1, NB2, NB3, NB-S) |

### Installation

```bash
# Create a virtual environment (recommended)
python -m venv .venv

# Activate the virtual environment
# macOS / Linux:
source .venv/bin/activate
# Windows (PowerShell):
.venv\Scripts\Activate.ps1
# Windows (Command Prompt):
.venv\Scripts\activate.bat

# Install packages
pip install numpy scipy matplotlib pandas jupyterlab
```

> **Windows note:** If you see `"jupyter" is not recognized`, make sure you have
> activated the virtual environment first. The `jupyter` command is only available
> inside the environment where it was installed.

### Running the Notebooks

```bash
cd notebooks/
jupyter lab
```

Open any `NB*.ipynb` file and run cells sequentially (Shift+Enter).

> **Alternative (VS Code):** Open any `.ipynb` file directly in VS Code. Select the
> `.venv` Python interpreter via the kernel picker in the top-right corner.

---

## Notebook Overview

Each notebook is the companion for the corresponding lecture note chapter: **NB*N* = Chapter *N***.

| Notebook | Chapter | Topic |
|----------|---------|-------|
| [NB1](#nb1) | Ch 1 | Bridge: CSV handling, Arrhenius fitting, residual analysis |
| [NB2](#nb2) | Ch 2 | Surface Chemistry: Adsorption, kinetics (LH/ER/MVK), temperature dependence |
| [NB3](#nb3) | Ch 3 | Reactor Design & Characterization: CSTR/PFR, BET, TPD, XRD |
| [NB4](#nb4) | Ch 4 | Transport: Thiele modulus, effectiveness factor, Weisz-Prater |
| [NB5](#nb5) | Ch 5 | Selectivity: Parallel/consecutive reactions, yield optimization |
| [NB6](#nb6) | Ch 6 | Zeolites & CNTs: Shape selectivity, micropore/enhanced diffusion |
| [NB-S](#nb-s) | Appendix | Capstone: Suzuki coupling, ee, deactivation, green metrics |

---

## Notebook Details

### NB1
**Chapter 1: Bridge Module — Python Computational Tools for Kinetic Data Analysis**

Your starting point. Introduces the four-step workflow used throughout the course: **Data Input** (CSV with pandas) → **Visualization** (matplotlib) → **Analysis** (linear regression) → **Validation** (residual plots). You will perform a complete Arrhenius analysis, extracting activation energy and pre-exponential factor with uncertainties.

**Key skills:** `pd.read_csv()`, `scipy.stats.linregress()`, publication-quality plots, residual interpretation

---

### NB2
**Chapter 2: Computational Perspectives on Surface Catalysis**

The largest notebook, covering the three pillars of surface catalysis computation:

- **Part I — Adsorption Fundamentals:** Langmuir isotherm, competitive adsorption, coverage maps, TPD analysis with the Redhead equation
- **Part II — Surface Kinetics:** Langmuir-Hinshelwood (LH), Eley-Rideal (ER), and Mars-van Krevelen (MVK) mechanisms; fingerprint diagnostics for model discrimination
- **Part III — Temperature Dependence:** Arrhenius and Eyring analysis, apparent vs intrinsic activation energy, compensation effect

**Key functions:** `langmuir_coverage()`, `lh_rate_bimolecular()`, `er_rate()`, `mvk_rate()`, `arrhenius()`, `eyring()`, `redhead_equation()`

> ☕ Session breaks are suggested between each Part.

---

### NB3
**Chapter 3: Reactor Design and Catalyst Characterization**

Two complementary perspectives on catalytic processes:

- **Part I — Ideal Reactor Design:** Batch, CSTR, PFR design equations; CSTR vs PFR conversion comparison; CSTRs-in-series; TOF and TON calculations
- **Part II — Catalyst Characterization:** BET surface area from N₂ adsorption; TPD analysis with Redhead equation; XRD crystallite sizing with Scherrer equation; metal dispersion and active surface area

**Key functions:** `cstr_conversion()`, `pfr_conversion()`, `bet_analysis()`, `redhead_equation()`, `scherrer_size()`, `dispersion_from_size()`

> ☕ Session break suggested between reactor design and characterization.

---

### NB4
**Chapter 4: Transport Limitations and Effectiveness Factors**

Visualizes the Thiele modulus and effectiveness factor for slab and sphere geometries. Plots internal concentration profiles showing reactant depletion. Performs sensitivity analysis (which parameter best restores η near 1?), implements the Weisz-Prater diagnostic, and integrates diffusion-limited kinetics into a PFR model.

**Key functions:** `effectiveness_slab()`, `effectiveness_sphere()`, `weisz_prater()`, `thiele_modulus()`

---

### NB5
**Chapter 5: Selectivity Engineering**

Models parallel reactions (constant selectivity) and consecutive reactions (yield maximum for the intermediate). Compares PFR vs CSTR for intermediate yield, adds Arrhenius temperature dependence to explore selectivity-temperature tradeoffs, and demonstrates N CSTRs in series approaching PFR behavior.

**Key functions:** `parallel_pfr_analytical()`, `consecutive_pfr_analytical()`, `consecutive_cstr()`, `tau_opt_pfr()`

---

### NB6
**Chapter 6: Zeolites and Carbon Nanotubes in Catalysis**

Two contrasting nanomaterials for catalysis:

- **Part I — Zeolites:** Framework types (FAU, MFI, CHA), Brønsted acid site density, three types of shape selectivity, configurational diffusion, effectiveness factor for microporous crystals, para-xylene selectivity case study
- **Part II — Carbon Nanotubes:** Chiral indices and electronic character, (n,m) classification map, classical Knudsen vs enhanced CNT diffusivity (F = 10–10,000), design optimization for diffusion-limited catalysis

**Key functions:** `acid_site_density()`, `configurational_diffusivity()`, `cnt_diameter()`, `knudsen_diffusivity()`, `effectiveness_cnt()`

> ☕ Session break suggested between zeolites and carbon nanotubes.

---

### NB-S
**Appendix: Catalytic Organic Synthesis — Supplementary Capstone**

Integrates concepts from the entire course. Calculates TON and TOF from Suzuki coupling data at different catalyst loadings. Predicts enantiomeric excess from transition-state energy differences. Fits catalyst deactivation kinetics to extract TOF₀ and k_d. Evaluates green chemistry metrics (atom economy, E-factor). Concludes with a comprehensive catalyst comparison using a radar chart.

**Key functions:** `calculate_ton_tof()`, `ee_from_ddG()`, `deactivation_model()`, `atom_economy()`, `e_factor()`

---

## Data Files

All sample datasets are in the `data/` directory:

| File | Description | Used by |
|------|-------------|---------|
| `arrhenius_data_sample.csv` | Rate constants at 11 temperatures (350–600 K) | NB1, NB2 |
| `kinetic_data_sample.csv` | CO oxidation rates at 5 temperatures and 5 pressures | NB2 |
| `adsorption_isotherm_sample.csv` | N₂ adsorption on activated carbon at 77 K | NB2, NB3 |
| `bet_isotherm_data.csv` | N₂ adsorption on SBA-15 mesoporous silica | NB3 |
| `suzuki_coupling_data.csv` | Suzuki coupling conversion at 3 Pd loadings | NB-S |
| `tof_deactivation_data.csv` | TOF decay over time (with noise) for fitting | NB-S |

---

## Shared Utility Module (`utils.py`)

The `utils.py` module provides 30 reusable functions organized by chapter. Each notebook defines its functions locally for pedagogical clarity, but `utils.py` provides a consolidated reference implementation.

### Function Reference

**Adsorption (Ch 2)**
| Function | Description |
|----------|-------------|
| `langmuir_coverage(P, K)` | Langmuir isotherm: θ = KP/(1+KP) |
| `competitive_langmuir(P_A, K_A, P_I, K_I)` | Competitive adsorption with inhibitor |

**Temperature Dependence (Ch 2)**
| Function | Description |
|----------|-------------|
| `arrhenius(T, A, Ea)` | k = A exp(−Ea/RT) |
| `eyring(T, dH, dS)` | Transition state theory rate constant |
| `K_temperature(T, K0, delta_H)` | van't Hoff equation for K(T) |

**Surface Kinetics (Ch 2)**
| Function | Description |
|----------|-------------|
| `lh_rate(P_A, K_A, k, P_B, K_B)` | Langmuir-Hinshelwood rate |
| `er_rate(theta_A, P_B, k_er)` | Eley-Rideal rate |
| `mvk_rate(P_A, P_O2, k_ox, k_red)` | Mars-van Krevelen rate |

**Reactor Design (Ch 3)**
| Function | Description |
|----------|-------------|
| `cstr_conversion(tau, k, C_A0, order)` | CSTR conversion (1st/2nd order) |
| `pfr_conversion(tau, k, C_A0, order)` | PFR conversion (1st/2nd order) |

**Characterization (Ch 3)**
| Function | Description |
|----------|-------------|
| `bet_transform(P, P0, V_ads)` | BET linearization |
| `bet_surface_area(V_m, sigma, N_A)` | Surface area from monolayer volume |
| `redhead_equation(T_peak, nu, beta)` | Desorption energy from TPD |
| `scherrer_size(K, wavelength, beta_rad, theta_rad)` | Crystallite size from XRD |

**Transport (Ch 4)**
| Function | Description |
|----------|-------------|
| `thiele_modulus(L, k, D_eff)` | Thiele modulus for slab |
| `effectiveness_slab(phi)` | η = tanh(φ)/φ |
| `effectiveness_sphere(phi)` | η = 3(φ coth(φ)−1)/φ² |
| `weisz_prater(phi, eta)` | Weisz-Prater criterion |

**Selectivity (Ch 5)**
| Function | Description |
|----------|-------------|
| `parallel_selectivity(k1, k2, C_A, a1, a2)` | Instantaneous selectivity |
| `consecutive_tau_opt(k1, k2)` | Optimal τ for max intermediate yield |

**Zeolites (Ch 6)**
| Function | Description |
|----------|-------------|
| `configurational_diffusivity(D0, E_diff, T)` | Micropore diffusion |

**Carbon Nanotubes (Ch 6)**
| Function | Description |
|----------|-------------|
| `cnt_diameter(n, m, a)` | Diameter and chirality from (n,m) |
| `knudsen_diffusivity(T, M_gmol, d_pore_nm)` | Classical Knudsen diffusion |
| `effectiveness_cnt(phi_values, F)` | Effectiveness with enhancement factor |

**Catalytic Synthesis (Appendix)**
| Function | Description |
|----------|-------------|
| `calculate_ton_tof(time, conversion, cat_loading)` | TON and TOF from data |
| `ee_from_ddG(ddG, T)` | Enantiomeric excess from energy difference |
| `deactivation_model(time, TOF0, kd)` | First-order deactivation |
| `atom_economy(MW_product, sum_MW_reactants)` | Atom economy percentage |
| `e_factor(mass_waste, mass_product)` | Environmental factor |

### Using `utils.py` in Your Own Scripts

```python
import sys
sys.path.insert(0, '/path/to/notebooks/')
from utils import langmuir_coverage, arrhenius, effectiveness_slab

# Example: coverage at P=1 bar, K=5 bar^-1
theta = langmuir_coverage(1.0, 5.0)
print(f"Coverage = {theta:.3f}")  # 0.833
```

---

## Colorblind-Safe Color Palette

All figures across every notebook use the Wong (2011) palette:

| Name | Hex | Usage |
|------|-----|-------|
| Blue | `#0072B2` | Primary data series |
| Orange | `#E69F00` | Secondary data series |
| Green | `#009E73` | Tertiary data / residuals |
| Vermillion | `#D55E00` | Fitted lines / warnings |
| Sky blue | `#56B4E9` | Backgrounds / fills |
| Purple | `#CC79A7` | Additional series |

---

## Tips for Students

1. **Run cells in order.** Each notebook builds on previous cells. If you get a `NameError`, restart the kernel and run from the top.

2. **Modify parameters.** Sections marked `# ADJUSTABLE PARAMETERS` are designed for exploration. Change values and re-run to build intuition.

3. **Complete the exercises.** Cells with `pass  # Replace with your implementation` are your practice problems. Try before checking solutions.

4. **Check units.** The most common source of errors in catalysis calculations is unit mismatch. Functions document expected units in their docstrings.

5. **Use the reflection questions.** They connect computation to physical understanding — the real goal of this course.

6. **Cross-reference the lecture notes.** Each notebook is the companion for the corresponding chapter: NB*N* = Chapter *N*.

7. **Use session breaks.** Longer notebooks (NB2, NB3, NB6) include ☕ session break markers. Save your work and take a break at those points.

---

## Notebook Dependency Chain

```
NB1 (Ch 1: Bridge)
 |
NB2 (Ch 2: Adsorption + Kinetics + Temperature)
 |
NB3 (Ch 3: Reactors + Characterization)
 |
 +----------+-----------+
 |                      |
NB4 (Ch 4: Transport)  NB5 (Ch 5: Selectivity)
 |                      |
NB6 (Ch 6: Zeolites + CNTs)
 |
NB-S (Appendix: Capstone)
```

---

*Course 2302337 — Surface Phenomena & Catalysis, Chulalongkorn University*
