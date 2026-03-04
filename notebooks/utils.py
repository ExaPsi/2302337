"""
Utility functions for Surface Phenomena & Catalysis course notebooks.

This module provides core functions for:
- Adsorption isotherms (Langmuir single and competitive)
- Temperature dependence (Arrhenius, Eyring, van't Hoff)
- Kinetic models (Langmuir-Hinshelwood, Eley-Rideal, Mars-van Krevelen)
- Reactor design (CSTR and PFR solvers)
- Catalyst characterization (BET, Redhead TPD, Scherrer XRD)
- Transport limitations (Thiele modulus, effectiveness factors)
- Selectivity engineering (parallel/consecutive reactions)
- Zeolite micropore diffusion
- Carbon nanotube enhanced transport
- Catalytic synthesis metrics (TON, TOF, ee, atom economy, E-factor)

Course: 2302337, Chulalongkorn University
"""

import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp

# ==============================================================================
# Constants
# ==============================================================================
R = 8.314  # Gas constant in J/(mol*K)

# Colorblind-safe color palette (Wong, 2011)
COLORS = {
    'blue': '#0072B2',
    'orange': '#E69F00',
    'green': '#009E73',
    'yellow': '#F0E442',
    'skyblue': '#56B4E9',
    'vermillion': '#D55E00',
    'purple': '#CC79A7',
    'black': '#000000'
}


# ==============================================================================
# Plotting Utilities
# ==============================================================================

def setup_plot_style():
    """
    Set up matplotlib style for publication-quality plots.
    Uses colorblind-safe colors.
    """
    import matplotlib.pyplot as plt

    plt.rcParams['figure.figsize'] = (8, 6)
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['legend.fontsize'] = 11
    plt.rcParams['xtick.labelsize'] = 11
    plt.rcParams['ytick.labelsize'] = 11
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.alpha'] = 0.3


# ==============================================================================
# Chapter 2: Adsorption Isotherms
# ==============================================================================

def langmuir_coverage(P, K):
    """Calculate fractional surface coverage using Langmuir isotherm.

    The Langmuir isotherm assumes:
    - Monolayer adsorption
    - Homogeneous surface (all sites equivalent)
    - No lateral interactions between adsorbates

    Parameters
    ----------
    P : float or array-like
        Partial pressure of adsorbate (bar or any consistent unit)
    K : float
        Adsorption equilibrium constant (bar^-1 or inverse of P units)
        K = k_ads / k_des = exp(-deltaG_ads / RT)

    Returns
    -------
    float or ndarray
        Fractional surface coverage theta, range [0, 1]

    Examples
    --------
    >>> langmuir_coverage(0.5, 2.0)  # P=0.5 bar, K=2 bar^-1
    0.5
    >>> langmuir_coverage(1.0, 2.0)
    0.6666666666666666

    Notes
    -----
    theta = K*P / (1 + K*P)

    At low P: theta ~ K*P (linear, Henry's law regime)
    At high P: theta ~ 1 (saturation)
    """
    P = np.asarray(P)
    return K * P / (1 + K * P)


def competitive_langmuir(P_A, K_A, P_I=0, K_I=0):
    """Calculate surface coverages for competitive Langmuir adsorption.

    Two species A and I compete for the same adsorption sites.
    This is the basis for inhibition/poisoning effects in catalysis.

    Parameters
    ----------
    P_A : float or array-like
        Partial pressure of species A (bar)
    K_A : float
        Adsorption equilibrium constant for A (bar^-1)
    P_I : float or array-like, optional
        Partial pressure of inhibitor I (bar), default 0
    K_I : float, optional
        Adsorption equilibrium constant for I (bar^-1), default 0

    Returns
    -------
    theta_A : float or ndarray
        Fractional coverage of species A
    theta_I : float or ndarray
        Fractional coverage of species I
    theta_star : float or ndarray
        Fraction of vacant sites

    Examples
    --------
    >>> theta_A, theta_I, theta_star = competitive_langmuir(1.0, 2.0, 0.01, 100)
    >>> print(f"theta_A={theta_A:.3f}, theta_I={theta_I:.3f}")
    theta_A=0.500, theta_I=0.250

    Notes
    -----
    Site balance: theta_A + theta_I + theta_star = 1

    theta_A = K_A * P_A / (1 + K_A*P_A + K_I*P_I)
    theta_I = K_I * P_I / (1 + K_A*P_A + K_I*P_I)
    """
    P_A = np.asarray(P_A)
    P_I = np.asarray(P_I)

    denom = 1 + K_A * P_A + K_I * P_I
    theta_A = K_A * P_A / denom
    theta_I = K_I * P_I / denom
    theta_star = 1 / denom

    return theta_A, theta_I, theta_star


# ==============================================================================
# Temperature Dependence Functions
# ==============================================================================

def arrhenius(T, A, Ea):
    """Calculate rate constant using Arrhenius equation.

    The Arrhenius equation describes the temperature dependence
    of reaction rate constants.

    Parameters
    ----------
    T : float or array-like
        Temperature in Kelvin
    A : float
        Pre-exponential factor (same units as k)
    Ea : float
        Activation energy in J/mol

    Returns
    -------
    float or ndarray
        Rate constant k (units depend on reaction order)

    Examples
    --------
    >>> arrhenius(500, 1e13, 100000)  # T=500K, A=1e13 s^-1, Ea=100 kJ/mol
    4.939...e-02

    Notes
    -----
    k = A * exp(-Ea / (R*T))

    Linearized form for fitting: ln(k) = ln(A) - Ea/(R*T)
    Plot ln(k) vs 1/T gives slope = -Ea/R
    """
    T = np.asarray(T)
    return A * np.exp(-Ea / (R * T))


def eyring(T, dH, dS, kB=1.380649e-23, h=6.62607e-34):
    """Calculate rate constant using Eyring equation (Transition State Theory).

    Parameters
    ----------
    T : float or array-like
        Temperature in Kelvin
    dH : float
        Activation enthalpy in J/mol
    dS : float
        Activation entropy in J/(mol*K)
    kB : float, optional
        Boltzmann constant in J/K
    h : float, optional
        Planck constant in J*s

    Returns
    -------
    float or ndarray
        Rate constant k in s^-1

    Notes
    -----
    k = (kB*T/h) * exp(dS/R) * exp(-dH/(R*T))

    Linearized form: ln(k/T) = ln(kB/h) + dS/R - dH/(R*T)
    Plot ln(k/T) vs 1/T gives slope = -dH/R
    """
    T = np.asarray(T)
    return (kB * T / h) * np.exp(dS / R) * np.exp(-dH / (R * T))


def K_temperature(T, K0, delta_H):
    """Calculate equilibrium constant K(T) using van't Hoff equation.

    Parameters
    ----------
    T : float or array-like
        Temperature (K)
    K0 : float
        Pre-exponential factor (same units as K)
    delta_H : float
        Enthalpy change (J/mol), negative for exothermic

    Returns
    -------
    K : float or ndarray
        Equilibrium constant at temperature T

    Notes
    -----
    K(T) = K0 * exp(-delta_H / RT)

    For exothermic adsorption, delta_H < 0, so K decreases with increasing T.
    """
    T = np.asarray(T)
    return K0 * np.exp(-delta_H / (R * T))


# ==============================================================================
# Chapter 2: Catalytic Rate Laws (continued)
# ==============================================================================

def lh_rate(P_A, K_A, k, P_B=None, K_B=None):
    """Calculate reaction rate using Langmuir-Hinshelwood kinetics.

    For unimolecular: r = k * theta_A
    For bimolecular: r = k * theta_A * theta_B

    Parameters
    ----------
    P_A : float or array-like
        Partial pressure of reactant A (bar)
    K_A : float
        Adsorption equilibrium constant for A (bar^-1)
    k : float
        Surface reaction rate constant
    P_B : float or array-like, optional
        Partial pressure of reactant B (bar) for bimolecular
    K_B : float, optional
        Adsorption equilibrium constant for B (bar^-1)

    Returns
    -------
    float or ndarray
        Reaction rate

    Notes
    -----
    Unimolecular: r = k * K_A * P_A / (1 + K_A * P_A)
    Bimolecular: r = k * K_A*P_A * K_B*P_B / (1 + K_A*P_A + K_B*P_B)^2
    """
    P_A = np.asarray(P_A)

    if P_B is None:
        # Unimolecular case
        theta_A = K_A * P_A / (1 + K_A * P_A)
        return k * theta_A
    else:
        # Bimolecular case
        P_B = np.asarray(P_B)
        denom = 1 + K_A * P_A + K_B * P_B
        return k * K_A * P_A * K_B * P_B / denom**2


def er_rate(theta_A, P_B, k_er):
    """Calculate Eley-Rideal reaction rate.

    r = k_ER * theta_A * P_B

    Parameters
    ----------
    theta_A : float or array_like
        Coverage of adsorbed species A
    P_B : float or array_like
        Pressure of gas-phase species B (bar)
    k_er : float
        Eley-Rideal rate constant

    Returns
    -------
    r : float or ndarray
        Reaction rate
    """
    return k_er * np.asarray(theta_A) * np.asarray(P_B)


def mvk_rate(P_A, P_O2, k_ox, k_red):
    """Calculate Mars-van Krevelen (MVK) reaction rate.

    The MVK mechanism describes oxidation reactions on metal oxide
    catalysts where lattice oxygen participates in the reaction.
    The surface cycles between oxidized and reduced states.

    Parameters
    ----------
    P_A : float or array-like
        Partial pressure of substrate A (bar)
    P_O2 : float or array-like
        Partial pressure of O2 (bar)
    k_ox : float
        Oxidation rate constant (rate of substrate oxidation by lattice O)
    k_red : float
        Re-oxidation rate constant (rate of lattice re-oxidation by O2)

    Returns
    -------
    float or ndarray
        Reaction rate

    Notes
    -----
    r = (k_ox * P_A * k_red * P_O2) / (k_ox * P_A + k_red * P_O2)

    Limiting cases:
    - k_red * P_O2 >> k_ox * P_A: r ~ k_ox * P_A (oxidation-limited)
    - k_ox * P_A >> k_red * P_O2: r ~ k_red * P_O2 (re-oxidation-limited)
    """
    P_A = np.asarray(P_A, dtype=float)
    P_O2 = np.asarray(P_O2, dtype=float)
    numerator = k_ox * P_A * k_red * P_O2
    denominator = k_ox * P_A + k_red * P_O2
    with np.errstate(divide='ignore', invalid='ignore'):
        rate = np.where(denominator < 1e-30, 0.0, numerator / denominator)
    return rate


# ==============================================================================
# Chapter 3: Reactor Design Equations
# ==============================================================================

def cstr_conversion(tau, k, C_A0, order=1):
    """Calculate CSTR conversion for nth-order reaction.

    Parameters
    ----------
    tau : float or array-like
        Space time (s)
    k : float
        Rate constant (appropriate units for order)
    C_A0 : float
        Inlet concentration (mol/L)
    order : int
        Reaction order (1 or 2)

    Returns
    -------
    X : float or ndarray
        Conversion (dimensionless)

    Notes
    -----
    First order: X = k*tau / (1 + k*tau)
    Second order: X = ((2*Da + 1) - sqrt(4*Da + 1)) / (2*Da), where Da = k*tau*C_A0
    """
    tau = np.asarray(tau)

    if order == 1:
        return k * tau / (1 + k * tau)
    elif order == 2:
        Da = k * tau * C_A0  # Damkohler number
        return ((2 * Da + 1) - np.sqrt(4 * Da + 1)) / (2 * Da)
    else:
        raise ValueError(f"Order {order} not implemented")


def pfr_conversion(tau, k, C_A0, order=1):
    """Calculate PFR conversion for nth-order reaction.

    Parameters
    ----------
    tau : float or array-like
        Space time (s)
    k : float
        Rate constant (appropriate units for order)
    C_A0 : float
        Inlet concentration (mol/L)
    order : int
        Reaction order (1 or 2)

    Returns
    -------
    X : float or ndarray
        Conversion (dimensionless)

    Notes
    -----
    First order: X = 1 - exp(-k*tau)
    Second order: X = k*tau*C_A0 / (1 + k*tau*C_A0)
    """
    tau = np.asarray(tau)

    if order == 1:
        return 1 - np.exp(-k * tau)
    elif order == 2:
        Da = k * tau * C_A0
        return Da / (1 + Da)
    else:
        raise ValueError(f"Order {order} not implemented")


# ==============================================================================
# Chapter 3: Catalyst Characterization (continued)
# ==============================================================================

def bet_transform(P, P0, V_ads):
    """Compute BET linearization: P/(V*(P0-P)) vs P/P0.

    The BET equation in linearized form:
        P / [V_ads * (P0 - P)] = 1/(V_m*C) + (C-1)/(V_m*C) * (P/P0)

    Plotting the left-hand side vs P/P0 gives a straight line
    from which V_m and C can be extracted.

    Parameters
    ----------
    P : float or array-like
        Equilibrium pressure (Pa or any consistent unit)
    P0 : float
        Saturation pressure (same units as P)
    V_ads : float or array-like
        Volume adsorbed at pressure P (cm3(STP)/g or mol/kg)

    Returns
    -------
    x : float or ndarray
        Relative pressure P/P0 (dimensionless)
    y : float or ndarray
        BET transform P / [V_ads * (P0 - P)] (inverse volume units)

    Notes
    -----
    The linear BET region is typically 0.05 < P/P0 < 0.35.
    Outside this range, the BET model breaks down.
    """
    P = np.asarray(P, dtype=float)
    V_ads = np.asarray(V_ads, dtype=float)
    x = P / P0
    y = P / (V_ads * (P0 - P))
    return x, y


def bet_surface_area(V_m, sigma=0.162e-18, N_A=6.022e23):
    """Calculate BET surface area from monolayer volume.

    Parameters
    ----------
    V_m : float
        Monolayer volume in cm3(STP)/g
    sigma : float, optional
        Cross-sectional area of one adsorbate molecule (m2).
        Default 0.162 nm2 = 0.162e-18 m2 for N2 at 77 K.
    N_A : float, optional
        Avogadro's number (mol^-1). Default 6.022e23.

    Returns
    -------
    float
        BET surface area in m2/g

    Notes
    -----
    S_BET = V_m * N_A * sigma / V_molar

    where V_molar = 22414 cm3(STP)/mol is the molar volume at STP.
    """
    V_molar = 22414.0  # cm3(STP)/mol
    return V_m * N_A * sigma / V_molar


def redhead_equation(T_peak, nu=1e13, beta=10.0):
    """Calculate desorption energy from TPD peak temperature using Redhead equation.

    Parameters
    ----------
    T_peak : float
        Temperature at peak maximum (K)
    nu : float, optional
        Pre-exponential factor (s^-1). Default 1e13 s^-1.
    beta : float, optional
        Heating rate (K/s). Default 10.0 K/s.
        Note: Convert K/min to K/s by dividing by 60.

    Returns
    -------
    float
        Desorption activation energy (J/mol)

    Notes
    -----
    E_d = R * T_p * [ln(nu * T_p / beta) - 3.64]

    Rule of thumb: E_d (kJ/mol) ~ 0.25 * T_p (K) for typical parameters.
    Valid when E_d/(R*T_p) is between 20 and 50.
    """
    return R * T_peak * (np.log(nu * T_peak / beta) - 3.64)


def scherrer_size(K, wavelength, beta_rad, theta_rad):
    """Calculate crystallite size from XRD peak broadening using Scherrer equation.

    Parameters
    ----------
    K : float
        Scherrer constant (shape factor), typically 0.89-0.94.
        Use 0.9 for spherical crystallites.
    wavelength : float
        X-ray wavelength (nm). Cu K-alpha = 0.15406 nm.
    beta_rad : float
        Full width at half maximum (FWHM) of the XRD peak (radians).
        Convert from degrees: beta_rad = beta_deg * pi / 180.
    theta_rad : float
        Bragg angle (radians), half of the 2-theta peak position.

    Returns
    -------
    float
        Crystallite size (nm)

    Notes
    -----
    d = K * lambda / (beta * cos(theta))

    The Scherrer equation gives a lower bound on crystallite size.
    Strain broadening is not accounted for.
    """
    return K * wavelength / (beta_rad * np.cos(theta_rad))


# ==============================================================================
# Chapter 4: Transport Limitations
# ==============================================================================

def thiele_modulus(L, k, D_eff, order=1):
    """Calculate Thiele modulus for slab geometry.

    Parameters
    ----------
    L : float
        Half-thickness of slab (m)
    k : float
        Rate constant (s^-1 for first order)
    D_eff : float
        Effective diffusivity (m^2/s)
    order : int
        Reaction order (currently only 1 supported)

    Returns
    -------
    phi : float
        Thiele modulus (dimensionless)

    Notes
    -----
    phi = L * sqrt(k / D_eff)
    """
    if order == 1:
        return L * np.sqrt(k / D_eff)
    else:
        raise ValueError(f"Order {order} not implemented")


def effectiveness_slab(phi):
    """Calculate effectiveness factor for slab (flat plate) geometry.

    Parameters
    ----------
    phi : float or array-like
        Thiele modulus (dimensionless)

    Returns
    -------
    float or ndarray
        Effectiveness factor eta, range (0, 1]

    Examples
    --------
    >>> effectiveness_slab(0.1)  # Kinetic control
    0.9967...
    >>> effectiveness_slab(10)   # Diffusion control
    0.0999...

    Notes
    -----
    eta = tanh(phi) / phi

    Limiting cases:
    - phi << 1: eta ~ 1 (kinetic control)
    - phi >> 1: eta ~ 1/phi (diffusion control)
    """
    phi = np.asarray(phi, dtype=float)

    with np.errstate(divide='ignore', invalid='ignore'):
        eta = np.where(phi < 1e-10, 1.0, np.tanh(phi) / phi)

    return eta


def effectiveness_sphere(phi):
    """Calculate effectiveness factor for spherical geometry.

    Parameters
    ----------
    phi : float or array-like
        Thiele modulus (dimensionless)

    Returns
    -------
    float or ndarray
        Effectiveness factor eta, range (0, 1]

    Notes
    -----
    eta = 3 * (phi * coth(phi) - 1) / phi^2

    Limiting cases:
    - phi << 1: eta ~ 1 (kinetic control)
    - phi >> 1: eta ~ 3/phi (diffusion control)
    """
    phi = np.asarray(phi, dtype=float)

    with np.errstate(divide='ignore', invalid='ignore'):
        coth_phi = 1.0 / np.tanh(phi)
        eta = 3 * (phi * coth_phi - 1) / phi**2
        eta = np.where(phi < 1e-10, 1.0, eta)

    return eta


def weisz_prater(phi, eta):
    """Calculate Weisz-Prater criterion for diffusion limitation diagnosis.

    Parameters
    ----------
    phi : float or array-like
        Thiele modulus
    eta : float or array-like
        Effectiveness factor

    Returns
    -------
    float or ndarray
        Weisz-Prater parameter (phi^2 * eta)

    Notes
    -----
    Interpretation:
    - WP < 0.3: Negligible internal diffusion effects
    - 0.3 < WP < 3: Transition regime
    - WP > 3: Severe diffusion limitation
    """
    return phi**2 * eta


# ==============================================================================
# Chapter 5: Selectivity Engineering
# ==============================================================================

def parallel_selectivity(k1, k2, C_A, alpha1, alpha2):
    """Calculate instantaneous selectivity for parallel reactions.

    For competing reactions A -> B (desired) and A -> C (undesired):
        r1 = k1 * C_A^alpha1
        r2 = k2 * C_A^alpha2

    Parameters
    ----------
    k1 : float
        Rate constant for desired reaction A -> B
    k2 : float
        Rate constant for undesired reaction A -> C
    C_A : float or array-like
        Concentration of reactant A (mol/L)
    alpha1 : float
        Reaction order of desired reaction
    alpha2 : float
        Reaction order of undesired reaction

    Returns
    -------
    float or ndarray
        Instantaneous selectivity S = r1 / (r1 + r2), range [0, 1]

    Notes
    -----
    S = r1 / (r1 + r2) = 1 / (1 + (k2/k1) * C_A^(alpha2 - alpha1))

    When alpha1 > alpha2: high C_A favors desired product (use PFR or batch)
    When alpha1 < alpha2: low C_A favors desired product (use CSTR)
    """
    C_A = np.asarray(C_A, dtype=float)
    r1 = k1 * C_A**alpha1
    r2 = k2 * C_A**alpha2
    with np.errstate(divide='ignore', invalid='ignore'):
        S = np.where((r1 + r2) < 1e-30, 0.0, r1 / (r1 + r2))
    return S


def consecutive_tau_opt(k1, k2):
    """Calculate optimal space time for maximum intermediate yield.

    For consecutive reactions A -> B -> C with first-order kinetics:
        A --k1--> B --k2--> C

    The intermediate B reaches maximum concentration at:
        tau_opt = ln(k2/k1) / (k2 - k1)

    Parameters
    ----------
    k1 : float
        Rate constant for A -> B (s^-1)
    k2 : float
        Rate constant for B -> C (s^-1)

    Returns
    -------
    float
        Optimal space time for maximum B yield (s)

    Notes
    -----
    At tau_opt, the maximum yield of B is:
        Y_B_max = (k1/k2)^(k2/(k2-k1))

    Requires k1 != k2. When k1 = k2, tau_opt = 1/k1.
    """
    if abs(k1 - k2) < 1e-15 * max(abs(k1), abs(k2)):
        # Degenerate case k1 ~ k2
        return 1.0 / k1
    return np.log(k2 / k1) / (k2 - k1)


# ==============================================================================
# Chapter 6: Zeolite Micropore Diffusion
# ==============================================================================

def configurational_diffusivity(D0, E_diff, T):
    """Calculate configurational (micropore) diffusivity in zeolites.

    Diffusion in zeolite micropores is an activated process because
    molecules must squeeze through windows between cavities.

    Parameters
    ----------
    D0 : float
        Pre-exponential factor (m^2/s). Typically 1e-8 to 1e-6 m^2/s.
    E_diff : float
        Activation energy for diffusion (J/mol).
        Typically 10-80 kJ/mol depending on molecule-pore size match.
    T : float or array-like
        Temperature (K)

    Returns
    -------
    float or ndarray
        Configurational diffusivity (m^2/s)

    Notes
    -----
    D_config = D0 * exp(-E_diff / (R*T))

    Configurational diffusion is orders of magnitude slower than
    Knudsen or molecular diffusion (D ~ 1e-12 to 1e-16 m^2/s).
    """
    T = np.asarray(T, dtype=float)
    return D0 * np.exp(-E_diff / (R * T))


# ==============================================================================
# Chapter 6: Carbon Nanotube Transport (continued)
# ==============================================================================

def cnt_diameter(n, m, a=0.246):
    """Calculate CNT diameter and classify chirality.

    Parameters
    ----------
    n : int
        Chiral index n
    m : int
        Chiral index m (0 <= m <= n)
    a : float, optional
        Graphene lattice constant (nm). Default 0.246 nm.

    Returns
    -------
    diameter : float
        CNT diameter (nm)
    chirality : str
        Classification: 'armchair' (n=m), 'zigzag' (m=0), or 'chiral'

    Notes
    -----
    d = a * sqrt(n^2 + nm + m^2) / pi

    The (n,m) indices define the chiral vector on the graphene sheet.
    Armchair (n,n) and zigzag (n,0) are achiral high-symmetry cases.
    """
    diameter = a * np.sqrt(n**2 + n * m + m**2) / np.pi

    if n == m:
        chirality = 'armchair'
    elif m == 0:
        chirality = 'zigzag'
    else:
        chirality = 'chiral'

    return diameter, chirality


def knudsen_diffusivity(T, M_gmol, d_pore_nm):
    """Calculate classical Knudsen diffusivity in a cylindrical pore.

    Parameters
    ----------
    T : float or array-like
        Temperature (K)
    M_gmol : float
        Molar mass of diffusing species (g/mol)
    d_pore_nm : float
        Pore diameter (nm)

    Returns
    -------
    float or ndarray
        Knudsen diffusivity (m^2/s)

    Notes
    -----
    D_K = (d_pore / 3) * sqrt(8 * R * T / (pi * M))

    where d_pore is in meters and M is in kg/mol.
    This is the classical expression; CNTs often show enhancement
    factors F = 10 to 10,000 over this value.
    """
    T = np.asarray(T, dtype=float)
    d_pore_m = d_pore_nm * 1e-9  # Convert nm to m
    M_kgmol = M_gmol * 1e-3  # Convert g/mol to kg/mol
    return (d_pore_m / 3) * np.sqrt(8 * R * T / (np.pi * M_kgmol))


def effectiveness_cnt(phi_values, F=1):
    """Calculate effectiveness factor for CNT-enhanced diffusion.

    Uses slab geometry with an effective Thiele modulus that accounts
    for the CNT enhancement factor F.

    Parameters
    ----------
    phi_values : float or array-like
        Thiele modulus values calculated with classical diffusivity
        (phi = L * sqrt(k / D_classical))
    F : float, optional
        CNT diffusion enhancement factor (D_CNT = F * D_Knudsen).
        Default 1 (no enhancement). Typical values: 10-10,000.

    Returns
    -------
    float or ndarray
        Effectiveness factor eta

    Notes
    -----
    The effective Thiele modulus with enhancement is:
        phi_eff = phi_classical / sqrt(F)

    Enhanced diffusion reduces phi, pushing the system toward
    kinetic control (eta -> 1).
    """
    phi_values = np.asarray(phi_values, dtype=float)
    # Enhanced diffusion reduces the effective Thiele modulus
    phi_eff = phi_values / np.sqrt(F)
    return effectiveness_slab(phi_eff)


# ==============================================================================
# Appendix: Catalytic Synthesis Metrics
# ==============================================================================

def calculate_ton_tof(time, conversion, cat_loading):
    """Calculate turnover number and turnover frequency from conversion data.

    Parameters
    ----------
    time : array-like
        Time points (h)
    conversion : array-like
        Fractional conversion at each time point (0 to 1)
    cat_loading : float
        Catalyst loading as mole fraction (e.g., 0.01 for 1 mol%)

    Returns
    -------
    ton : ndarray
        Turnover number at each time point (dimensionless)
    tof : ndarray
        Turnover frequency at each time point (h^-1)

    Notes
    -----
    TON = conversion / cat_loading
    TOF = d(TON)/dt via centered finite differences (np.gradient)

    At complete conversion (X=1), TON = 1/cat_loading.
    For example, 1 mol% loading gives max TON = 100.
    """
    time = np.asarray(time, dtype=float)
    conversion = np.asarray(conversion, dtype=float)
    ton = conversion / cat_loading
    tof = np.gradient(ton, time)
    return ton, tof


def ee_from_ddG(ddG, T):
    """Calculate enantiomeric excess from transition state energy difference.

    Uses the Boltzmann distribution to predict product enantiomer ratio
    from the energy difference between diastereomeric transition states.

    Parameters
    ----------
    ddG : float or array-like
        Energy difference between diastereomeric transition states (J/mol).
        Positive when the R-enantiomer is favored:
        ddG = DeltaG_S_dagger - DeltaG_R_dagger
    T : float or array-like
        Temperature (K)

    Returns
    -------
    float or ndarray
        Enantiomeric excess as a percentage (0 to 100)

    Notes
    -----
    ratio = exp(ddG / (R*T))
    ee = (ratio - 1) / (ratio + 1) * 100

    Limiting cases:
    - ddG = 0: ee = 0% (racemic)
    - ddG >> RT: ee -> 100% (enantiopure)
    - T -> infinity: ee -> 0% (thermal scrambling)
    """
    ddG = np.asarray(ddG, dtype=float)
    T = np.asarray(T, dtype=float)
    ratio = np.exp(ddG / (R * T))
    ee = (ratio - 1) / (ratio + 1) * 100
    return ee


def deactivation_model(time, TOF0, kd):
    """First-order catalyst deactivation model.

    Models the exponential decay of turnover frequency as the catalyst
    loses active sites over time.

    Parameters
    ----------
    time : float or array-like
        Time points (h)
    TOF0 : float
        Initial turnover frequency (h^-1)
    kd : float
        Deactivation rate constant (h^-1)

    Returns
    -------
    ndarray
        TOF at each time point (h^-1)

    Notes
    -----
    TOF(t) = TOF0 * exp(-kd * t)

    Total TON = TOF0 / kd (integral from 0 to infinity)
    Half-life = ln(2) / kd
    Cumulative TON(t) = (TOF0/kd) * (1 - exp(-kd*t))
    """
    time = np.asarray(time, dtype=float)
    return TOF0 * np.exp(-kd * time)


def atom_economy(MW_product, sum_MW_reactants):
    """Calculate atom economy for a chemical reaction.

    Atom economy measures what fraction of reactant atoms end up
    in the desired product. Higher is better.

    Parameters
    ----------
    MW_product : float
        Molecular weight of the desired product (g/mol)
    sum_MW_reactants : float
        Sum of molecular weights of all reactants (g/mol)

    Returns
    -------
    float
        Atom economy as a percentage (0 to 100)

    Notes
    -----
    AE = (MW_product / sum_MW_reactants) * 100

    Atom economy is a theoretical metric. It does not account for
    solvents, excess reagents, catalysts, or workup waste.
    For a practical measure, use E-factor.

    Example: Suzuki coupling
        PhBr (157) + PhB(OH)2 (122) -> Ph-Ph (154) + B(OH)3 + NaBr
        AE = 154 / (157 + 122) * 100 = 55.2%
    """
    return (MW_product / sum_MW_reactants) * 100


def e_factor(mass_waste, mass_product):
    """Calculate the environmental factor (E-factor).

    E-factor measures total waste per unit product.
    Lower is better: an ideal green process has E-factor = 0.

    Parameters
    ----------
    mass_waste : float
        Total mass of waste produced (kg)
    mass_product : float
        Mass of desired product (kg)

    Returns
    -------
    float
        E-factor (dimensionless, lower is better)

    Notes
    -----
    Typical E-factors by sector:
    - Oil refining: < 0.1
    - Bulk chemicals: 1-5
    - Fine chemicals: 5-50
    - Pharmaceuticals: 25-100+
    """
    return mass_waste / mass_product


# ==============================================================================
# Test
# ==============================================================================

if __name__ == "__main__":
    print("Testing utils.py functions...")

    # Test langmuir_coverage
    theta = langmuir_coverage(0.5, 2.0)
    assert abs(theta - 0.5) < 1e-10, "langmuir_coverage test failed"
    print(f"  langmuir_coverage(0.5, 2.0) = {theta}")

    # Test competitive_langmuir
    theta_A, theta_I, theta_star = competitive_langmuir(1.0, 2.0, 0.01, 100)
    assert abs(theta_A - 0.5) < 1e-10, "competitive_langmuir test failed"
    print(f"  competitive_langmuir: theta_A = {theta_A:.3f}")

    # Test effectiveness_slab
    eta = effectiveness_slab(0.1)
    assert abs(eta - 0.9967) < 0.001, "effectiveness_slab test failed"
    print(f"  effectiveness_slab(0.1) = {eta:.4f}")

    # Test mvk_rate
    rate = mvk_rate(1.0, 0.5, 2.0, 3.0)
    expected = (2.0 * 1.0 * 3.0 * 0.5) / (2.0 * 1.0 + 3.0 * 0.5)
    assert abs(rate - expected) < 1e-10, "mvk_rate test failed"
    print(f"  mvk_rate(1.0, 0.5, 2.0, 3.0) = {rate:.4f}")

    # Test bet_surface_area
    S = bet_surface_area(100.0)
    assert S > 0, "bet_surface_area should return positive value"
    print(f"  bet_surface_area(100) = {S:.1f} m2/g")

    # Test redhead_equation
    E_d = redhead_equation(400, nu=1e13, beta=10.0)
    assert E_d > 0, "redhead_equation should return positive energy"
    print(f"  redhead_equation(400 K) = {E_d/1000:.1f} kJ/mol")

    # Test parallel_selectivity
    S = parallel_selectivity(1.0, 0.5, 1.0, 1, 2)
    assert 0 < S < 1, "parallel_selectivity should be between 0 and 1"
    print(f"  parallel_selectivity(k1=1, k2=0.5, C=1, a1=1, a2=2) = {S:.4f}")

    # Test consecutive_tau_opt
    tau = consecutive_tau_opt(0.5, 0.1)
    assert tau > 0, "consecutive_tau_opt should return positive tau"
    print(f"  consecutive_tau_opt(0.5, 0.1) = {tau:.2f} s")

    # Test cnt_diameter
    d, chirality = cnt_diameter(10, 10)
    assert chirality == 'armchair', "cnt_diameter(10,10) should be armchair"
    print(f"  cnt_diameter(10,10) = {d:.3f} nm ({chirality})")

    # Test ee_from_ddG (lecture notes: 5 kJ/mol at 25C -> ~76%)
    ee = ee_from_ddG(5000, 298.15)
    assert 75 < ee < 78, f"ee_from_ddG(5000, 298.15) should be ~76%, got {ee:.1f}%"
    print(f"  ee_from_ddG(5000, 298.15) = {ee:.1f}%")

    # Test atom_economy (Suzuki: 154/(157+122) = 55.2%)
    ae = atom_economy(154, 157 + 122)
    assert abs(ae - 55.2) < 0.2, "atom_economy Suzuki test failed"
    print(f"  atom_economy(154, 279) = {ae:.1f}%")

    # Test deactivation_model
    tof = deactivation_model(0, 500, 0.1)
    assert abs(tof - 500) < 1e-10, "deactivation_model at t=0 should equal TOF0"
    print(f"  deactivation_model(0, 500, 0.1) = {tof:.0f} h^-1")

    print("\nAll tests passed!")
