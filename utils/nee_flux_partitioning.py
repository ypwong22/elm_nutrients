"""
NEE flux partitioning based on Anthony Walker's, but with the reference temperature
of Q10 set to 15 degC, and night-daytime threshold set to PAR = 50 following Jon Stelling

Reference:
    Walker, A. P. et al. (2017), *J. Geophys. Res. Biogeosci.*, 122, 1078–1097, doi:10.1002/2016JG003711.
"""

import numpy as np
from scipy.optimize import curve_fit
from dataclasses import dataclass
from typing import Optional


# ---------------------------------------------------------------------------
# Model functions
# ---------------------------------------------------------------------------
def q10_function(T, R_ref, Q10, T_ref):
    """Q10 temperature response function.
    
    R(T) = R_ref * Q10^((T - T_ref) / 10)
    
    Parameters:
        T     : temperature (°C or K, consistent with T_ref)
        R_ref : rate at reference temperature
        Q10   : factor by which rate increases per 10°C
        T_ref : reference temperature
    """
    return R_ref * Q10 ** ((T - T_ref) / 10.0)


def q10_wh_function(T, Wh, R_ref, rw, Q10, T_ref):
    """Q10 temperature response function.
    
    R(T) = (R_ref + rw * Wh) * Q10^((T - T_ref) / 10)
    
    Parameters:
        T     : temperature (°C or K, consistent with T_ref)
        R_ref : rate at reference temperature
        Q10   : factor by which rate increases per 10°C
        T_ref : reference temperature
        Wh    : water table depth (m)
    """
    return (R_ref + rw*Wh) * Q10 ** ((T - T_ref) / 10.0)


def rectangular_hyperbola(PAR, alpha, gmax_eff):
    """Rectangular hyperbola (Michaelis-Menten style light response).
    
    Returns: alpha * PAR * gmax_eff / (alpha * PAR + gmax_eff)
    """
    return alpha * PAR * gmax_eff / (alpha * PAR + gmax_eff)


def temp_scalar(T_C, Tmin, Tmax, Topt):
    """Beta temperature response (Yin et al. 1995).
    Exactly zero at Tmin and Tmax, peaks at 1.0 at Topt.
    Typical values: Tmin=0, Topt=20, Tmax=40 (°C)

    Yin, X., Kropff, M.J., McLaren, G., Visperas, R.M. (1995). 
    A nonlinear model for crop development as a function of temperature. 
    Agricultural and Forest Meteorology, 77(1–2), 1–16. 
    https://doi.org/10.1016/0168-1923(95)02236-Q
    """
    T_C = np.asarray(T_C, dtype=float)
    scalar = np.zeros_like(T_C)
    mask = (T_C > Tmin) & (T_C < Tmax)
    T = T_C[mask]
    alpha = (Tmax - Topt) / (Topt - Tmin)
    scalar[mask] = ((T - Tmin) / (Topt - Tmin)) * ((Tmax - T) / (Tmax - Topt)) ** alpha
    return scalar


def gpp_par_only(PAR, alpha, gmax):
    """GPP = rectangular hyperbola of PAR only."""
    return rectangular_hyperbola(PAR, alpha, gmax)


def gpp_par_wh(X, alpha, gmax, gw):
    """GPP = rectangular hyperbola with water table height modifier (Eq. 2).
    
    GPP = alpha * PAR * (gmax + gw * Wh) / (alpha * PAR + (gmax + gw * Wh))
    
    X : 2D array where X[0] = PAR, X[1] = Wh
    """
    PAR, Wh = X[0], X[1]
    gmax_eff = gmax + gw * Wh
    return rectangular_hyperbola(PAR, alpha, gmax_eff)


def gpp_par_tair(X, alpha, gmax, Tmin, Tmax, Topt):
    """GPP = rectangular hyperbola * temperature scalar.
    
    X : 2D array where X[0] = PAR, X[1] = Tair (°C)
    Parameters Ha, Hd in kJ mol⁻¹; Topt in °C (converted internally).
    """
    PAR, Tair = X[0], X[1]
    return rectangular_hyperbola(PAR, alpha, gmax) * temp_scalar(Tair, Tmin, Tmax, Topt)


def gpp_par_wh_tair(X, alpha, gmax, gw, Tmin, Tmax, Topt):
    """GPP = rectangular hyperbola (with Wh) * temperature scalar.
    
    X : 2D array where X[0] = PAR, X[1] = Wh, X[2] = Tair (°C)
    Parameters Ha, Hd in kJ mol⁻¹; Topt in °C (converted internally).
    """
    PAR, Wh, Tair = X[0], X[1], X[2]
    gmax_eff = gmax + gw * Wh
    return rectangular_hyperbola(PAR, alpha, gmax_eff) * temp_scalar(Tair, Tmin, Tmax, Topt)


def fit_q10(T, R, T_ref, Wh=None, p0=None):
    """Fit Q10 model to observed temperature–rate data.
    
    Parameters:
        T     : array-like, temperatures
        R     : array-like, observed rates
        Wh    : array-like, observed water table depth
        T_ref : reference temperature (fixed). If None, fitted as free parameter.
        p0    : initial guesses [R_ref, Q10] or [R_ref, Q10, T_ref]
    
    Returns:
        popt  : optimal parameters
        pcov  : covariance matrix
        label : tuple of parameter names
    """
    T = np.asarray(T, dtype=float)
    R = np.asarray(R, dtype=float)

    # Two-parameter fit (R_ref, Q10); T_ref is fixed
    if Wh is None:
        filt = np.isfinite(T) & np.isfinite(R)
        T = T[filt]
        R = R[filt]

        def model(T, R_ref, Q10):
            return q10_function(T, R_ref, Q10, T_ref)

        if p0 is None:
            p0 = [R.mean(), 2.0]
        popt, pcov = curve_fit(model, T, R, p0=p0, maxfev=10000)
        return popt, pcov, ("R_ref", "Q10"), model
    else:
        Wh = np.asarray(Wh, dtype=float)

        filt = np.isfinite(T) & np.isfinite(R) & np.isfinite(Wh)
        T = T[filt]
        R = R[filt]
        Wh = Wh[filt]

        def model(val, R_ref, rw, Q10):
            return q10_wh_function(val[0,:], val[1,:], R_ref, rw, Q10, T_ref)

        if p0 is None:
            p0 = [R.mean(), 1e-3, 2.0]
        popt, pcov = curve_fit(model, np.vstack([T,Wh]), R, p0=p0, maxfev=10000)
        return popt, pcov, ("R_ref", "rw", "Q10"), model


"""
# ── Example usage ──
if __name__ == "__main__":
    # Synthetic data
    np.random.seed(0)
    T_true = np.linspace(5, 35, 30)
    R_true = 2.0 * 2.3 ** ((T_true - 15.0) / 10.0)
    R_noisy = R_true + np.random.normal(0, 0.3, size=T_true.shape)

    # Fit with fixed T_ref
    popt, pcov, names = fit_q10(T_true, R_noisy, T_ref=15.0)
    perr = np.sqrt(np.diag(pcov))
    print("Fixed T_ref = 15.0")
    for n, v, e in zip(names, popt, perr):
        print(f"  {n} = {v:.4f} ± {e:.4f}")

    # Fit with T_ref as free parameter
    popt3, pcov3, names3 = fit_q10(T_true, R_noisy)
    perr3 = np.sqrt(np.diag(pcov3))
    print("\nFree T_ref")
    for n, v, e in zip(names3, popt3, perr3):
        print(f"  {n} = {v:.4f} ± {e:.4f}")
"""


def fit_gpp(PAR, GPP, Wh=None, Tair=None, p0=None, bounds=None, **kwargs):
    """Fit GPP model to data, selecting model complexity based on inputs.
    
    Parameters:
        PAR   : array, photosynthetically active radiation (µmol m⁻² s⁻¹)
        GPP   : array, observed GPP (µmol m⁻² s⁻¹)
        Wh    : array or None, water table height
        Tair  : array or None, air temperature (°C)
        p0    : initial parameter guesses (optional)
        bounds: parameter bounds for curve_fit (optional)
        **kwargs: passed to scipy.optimize.curve_fit
    
    Returns:
        popt   : fitted parameters
        pcov   : covariance matrix
        names  : tuple of parameter names (units)
        model  : the function that was fitted
    """
    PAR = np.asarray(PAR, dtype=float)
    GPP = np.asarray(GPP, dtype=float)

    use_wh = Wh is not None
    use_tair = Tair is not None

    filt = np.isfinite(PAR) & np.isfinite(GPP)
    if use_wh:
        Wh = np.asarray(Wh, dtype=float)
        filt &= np.isfinite(Wh)
    if use_tair:
        Tair = np.asarray(Tair, dtype=float)
        filt &= np.isfinite(Tair)
    PAR = PAR[filt]
    GPP = GPP[filt]
    if use_wh:
        Wh = Wh[filt]
    if use_tair:
        Tair = Tair[filt]

    # ── PAR only ──
    if not use_wh and not use_tair:
        model = gpp_par_only
        X = PAR
        names = ("alpha", "gmax")
        default_p0 = [0.05, np.max(GPP)]
        default_bounds = ([0, 0], [np.inf, np.inf])

    # ── PAR + Wh ──
    elif use_wh and not use_tair:
        model = gpp_par_wh
        X = np.vstack([PAR, Wh])
        names = ("alpha", "gmax", "gw")
        default_p0 = [0.05, np.max(GPP), 0.0]
        default_bounds = ([0, 0, -np.inf], [np.inf, np.inf, np.inf])

    # ── PAR + Tair ──
    elif not use_wh and use_tair:
        model = gpp_par_tair
        X = np.vstack([PAR, Tair])
        names = ("alpha", "gmax", "Tmin", "Tmax", "Topt") # [°C]
        default_p0 = [0.05, np.max(GPP), 0.0, 45.0, 25.0]
        default_bounds = ([0, 0, -5, 40, 10],
                          [np.inf, np.inf, 10, 60, 40])

    # ── PAR + Wh + Tair ──
    else:
        model = gpp_par_wh_tair
        X = np.vstack([PAR, Wh, Tair])
        names = ("alpha", "gmax", "gw", "Tmin", "Tmax", "Topt") #  [°C]
        default_p0 = [0.05, np.max(GPP), 0.0, 0, 45.0, 25.0]
        default_bounds = ([0, 0, -np.inf, -5, 40, 10],
                          [np.inf, np.inf, np.inf, 10, 60, 40])

    if p0 is None:
        p0 = default_p0
    if bounds is None:
        bounds = default_bounds

    popt, pcov = curve_fit(model, X, GPP, p0=p0, bounds=bounds,
                           maxfev=20000, **kwargs)

    return popt, pcov, names, model


"""
if __name__ == "__main__":
    np.random.seed(42)
    n = 200

    # Synthetic drivers
    PAR = np.random.uniform(0, 2000, n)
    Wh = np.random.uniform(-500, 0, n)
    Tair = np.random.uniform(5, 35, n)

    # True parameters
    alpha_t, gmax_t, gw_t = 0.04, 15.0, 0.01
    Ha_t, Hd_t, Topt_t = 55.0, 195.0, 25.0  # kJ/mol, kJ/mol, °C

    # Synthetic GPP (full model)
    GPP_true = gpp_par_wh_tair(
        np.vstack([PAR, Wh, Tair]),
        alpha_t, gmax_t, gw_t, Ha_t, Hd_t, Topt_t
    )
    GPP_obs = GPP_true + np.random.normal(0, 0.3, n)

    # --- Fit: PAR only ---
    popt, pcov, names, _ = fit_gpp(PAR, GPP_obs)
    print("=== PAR only ===")
    for nm, v, e in zip(names, popt, np.sqrt(np.diag(pcov))):
        print(f"  {nm:>15s} = {v:.4f} ± {e:.4f}")

    # --- Fit: PAR + Wh ---
    popt, pcov, names, _ = fit_gpp(PAR, GPP_obs, Wh=Wh)
    print("\n=== PAR + Wh ===")
    for nm, v, e in zip(names, popt, np.sqrt(np.diag(pcov))):
        print(f"  {nm:>15s} = {v:.4f} ± {e:.4f}")

    # --- Fit: PAR + Tair ---
    popt, pcov, names, _ = fit_gpp(PAR, GPP_obs, Tair=Tair)
    print("\n=== PAR + Tair ===")
    for nm, v, e in zip(names, popt, np.sqrt(np.diag(pcov))):
        print(f"  {nm:>15s} = {v:.4f} ± {e:.4f}")

    # --- Fit: PAR + Wh + Tair (full) ---
    popt, pcov, names, _ = fit_gpp(PAR, GPP_obs, Wh=Wh, Tair=Tair)
    print("\n=== PAR + Wh + Tair ===")
    for nm, v, e in zip(names, popt, np.sqrt(np.diag(pcov))):
        print(f"  {nm:>15s} = {v:.4f} ± {e:.4f}")
"""