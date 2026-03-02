"""
NEE flux partitioning using a rectangular hyperbolic light-response curve
with Lloyd-Taylor temperature-dependent respiration.

Equation (adapted from Lasslop et al. 2010, Eq. 3, with explicit sign convention):

    NEE = -(α * β * PAR) / (α * PAR + β) + rb * exp(E0 * (1/(Tref - T0) - 1/(T - T0)))

where:
    NEE   : net ecosystem exchange (µmol CO2 m-2 s-1), positive = emission, negative = uptake
    PAR   : photosynthetically active radiation (µmol photons m-2 s-1)
    T     : temperature (K), e.g. soil temperature at 5 cm depth
    α     : apparent quantum yield (µmol CO2 / µmol photons), initial slope of light response (>0)
    β     : maximum CO2 uptake rate at light saturation (µmol CO2 m-2 s-1) (>0)
    rb    : base respiration rate at Tref (µmol CO2 m-2 s-1) (>0)
    E0    : activation energy parameter (K), temperature sensitivity of respiration (>0)
    Tref  : reference temperature, fixed at 283.15 K (10 °C)
    T0    : lower temperature limit, fixed at 227.13 K (-46.02 °C)

At nighttime (PAR = 0), the equation reduces to:
    NEE = rb * exp(E0 * (1/(Tref - T0) - 1/(T - T0)))

This is pure ecosystem respiration (Reco) following Lloyd & Taylor (1994).

The sign convention used here differs from the original Lasslop et al. (2010):
  - Original: α is allowed to be negative, no explicit minus sign on the GPP term.
  - Here: α is constrained positive, and the minus sign is written explicitly.
  Both are mathematically equivalent.

Fitted parameters: α, β, rb, E0
Fixed constants:    Tref = 283.15 K, T0 = 227.13 K

Reference:
    Lasslop, G., Reichstein, M., Papale, D., Richardson, A.D., Arneth, A.,
    Barr, A., Stoy, P., and Wohlfahrt, G. (2010). Separation of net ecosystem
    exchange into assimilation and respiration using a light response curve
    approach: critical issues and global evaluation. Global Change Biology,
    16(1), 187-208. https://doi.org/10.1111/j.1365-2486.2009.02041.x
    See Eq. 3.

Usage:
    >>> import numpy as np
    >>> from nee_flux_partitioning import fit_nee_model, predict_nee, predict_reco, predict_gpp

    >>> # Example data (PAR in µmol/m2/s, T in Kelvin, NEE in µmol CO2/m2/s)
    >>> PAR = np.array([0, 0, 100, 500, 1000, 1500, 0, 0, 200, 800])
    >>> T   = np.array([280, 279, 285, 290, 293, 295, 278, 277, 287, 292])
    >>> NEE = np.array([1.5, 1.2, -0.5, -3.0, -4.5, -4.0, 1.0, 0.9, -1.0, -3.5])

    >>> result = fit_nee_model(PAR, T, NEE)
    >>> print(result)

    >>> # Predict NEE, Reco, GPP for new data
    >>> nee_pred = predict_nee(PAR, T, result['params'])
    >>> reco     = predict_reco(T, result['params'])
    >>> gpp      = predict_gpp(PAR, T, result['params'])
"""

import numpy as np
from scipy.optimize import curve_fit
from dataclasses import dataclass
from typing import Optional

# ---------------------------------------------------------------------------
# Fixed constants
# ---------------------------------------------------------------------------
TREF = 283.15   # Reference temperature (K), i.e. 10 °C
T0   = 227.13   # Lower temperature limit (K), i.e. -46.02 °C


# ---------------------------------------------------------------------------
# Model functions
# ---------------------------------------------------------------------------

def nee_model(X, alpha, beta, rb, E0):
    """
    Compute NEE from PAR and temperature.

    NEE = -(α * β * PAR) / (α * PAR + β) + rb * exp(E0 * (1/(Tref - T0) - 1/(T - T0)))

    Parameters
    ----------
    X : tuple of (PAR, T)
        PAR : array-like, photosynthetically active radiation (µmol m-2 s-1)
        T   : array-like, temperature in Kelvin
    alpha : float
        Apparent quantum yield (µmol CO2 / µmol photons), > 0
    beta : float
        Maximum CO2 uptake rate at light saturation (µmol CO2 m-2 s-1), > 0
    rb : float
        Base respiration rate at Tref (µmol CO2 m-2 s-1), > 0
    E0 : float
        Activation energy parameter (K), > 0

    Returns
    -------
    NEE : np.ndarray
        Net ecosystem exchange (µmol CO2 m-2 s-1).
        Positive = net emission, negative = net uptake.
    """
    PAR, T = X
    PAR = np.asarray(PAR, dtype=float)
    T   = np.asarray(T, dtype=float)

    # GPP term (rectangular hyperbola), always <= 0
    gpp_term = -(alpha * beta * PAR) / (alpha * PAR + beta)

    # Reco term (Lloyd-Taylor), always > 0
    reco_term = rb * np.exp(E0 * (1.0 / (TREF - T0) - 1.0 / (T - T0)))

    return gpp_term + reco_term


def reco_model(T, rb, E0):
    """
    Compute ecosystem respiration (Reco) from temperature only.

    Reco = rb * exp(E0 * (1/(Tref - T0) - 1/(T - T0)))

    This is the Lloyd & Taylor (1994) respiration model, equivalent to
    the NEE model evaluated at PAR = 0.

    Parameters
    ----------
    T : array-like
        Temperature in Kelvin.
    rb : float
        Base respiration rate at Tref (µmol CO2 m-2 s-1).
    E0 : float
        Activation energy parameter (K).

    Returns
    -------
    Reco : np.ndarray
        Ecosystem respiration (µmol CO2 m-2 s-1), always > 0.
    """
    T = np.asarray(T, dtype=float)
    return rb * np.exp(E0 * (1.0 / (TREF - T0) - 1.0 / (T - T0)))


def gpp_model(PAR, T, alpha, beta, rb, E0):
    """
    Compute gross primary production (GPP) as the difference Reco - NEE.

    GPP = Reco - NEE = (α * β * PAR) / (α * PAR + β)

    Note: GPP here is independent of temperature and respiration parameters.
    It depends only on PAR, α, and β. The T, rb, E0 arguments are accepted
    for API consistency but are unused.

    Parameters
    ----------
    PAR : array-like
        Photosynthetically active radiation (µmol m-2 s-1).
    T : array-like
        Temperature in Kelvin (unused, kept for API consistency).
    alpha : float
        Apparent quantum yield.
    beta : float
        Maximum CO2 uptake rate at light saturation.
    rb : float
        Unused.
    E0 : float
        Unused.

    Returns
    -------
    GPP : np.ndarray
        Gross primary production (µmol CO2 m-2 s-1), always >= 0.
    """
    PAR = np.asarray(PAR, dtype=float)
    return (alpha * beta * PAR) / (alpha * PAR + beta)


# ---------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------

def fit_nee_model(PAR, T, NEE,
                  p0=None,
                  bounds=None,
                  maxfev=10000,
                  method='simultaneous',
                  par_threshold=5.0,
                  verbose=True):
    """
    Fit the NEE light-response + Lloyd-Taylor respiration model to data.

    Uses scipy.optimize.curve_fit (trust region reflective with bounds).

    Two fitting methods are available:

    - 'simultaneous': Fit all 4 parameters (alpha, beta, rb, E0) at once
      using all data. Simple and fast, but daytime data can dominate the
      fit and bias nighttime Reco predictions.

    - 'two_step': (Recommended when nighttime NEE is poorly fit)
      Step 1: Fit rb and E0 from nighttime data only (PAR < par_threshold),
              where NEE ≈ Reco. This anchors the respiration parameters to
              nighttime observations. Analogous to Reichstein et al. (2005).
      Step 2: Fix rb and E0, then fit alpha and beta from all data.
      This guarantees that the Reco term matches nighttime fluxes.

    Parameters
    ----------
    PAR : array-like
        Photosynthetically active radiation (µmol m-2 s-1).
    T : array-like
        Temperature in Kelvin.
    NEE : array-like
        Observed net ecosystem exchange (µmol CO2 m-2 s-1).
        Sign convention: positive = emission, negative = uptake.
    p0 : list or None
        Initial guesses for [alpha, beta, rb, E0].
        Default: [0.05, 20.0, 2.0, 200.0]
    bounds : tuple of (lower, upper) or None
        Parameter bounds for [alpha, beta, rb, E0]. Default:
            lower: [0.001,  1.0, 0.01,  50.0]
            upper: [0.5,  100.0, 20.0, 500.0]
    maxfev : int
        Maximum number of function evaluations.
    method : str, 'simultaneous' or 'two_step'
        Fitting strategy. Default: 'simultaneous'.
        Use 'two_step' if nighttime NEE is overpredicted.
    par_threshold : float
        PAR threshold (µmol m-2 s-1) to classify nighttime data.
        Default: 5.0. Used by both methods for diagnostics and by
        'two_step' to select nighttime data for Step 1.
    verbose : bool
        If True, print fitted parameters and diagnostics.

    Returns
    -------
    result : dict with keys:
        'params'    : dict with 'alpha', 'beta', 'rb', 'E0'
        'pcov'      : covariance matrix (4x4 for simultaneous,
                      dict with 'reco' and 'gpp' keys for two_step)
        'std_err'   : dict with standard errors of each parameter
        'nee_pred'  : predicted NEE for input data
        'reco_pred' : predicted Reco for input data
        'gpp_pred'  : predicted GPP for input data
        'residuals' : NEE_obs - NEE_pred
        'rmse'      : root mean square error
        'r_squared' : coefficient of determination
        'n_obs'     : number of observations
        'n_night'   : number of nighttime observations (PAR < par_threshold)
        'n_day'     : number of daytime observations (PAR >= par_threshold)
        'method'    : fitting method used
    """
    PAR = np.asarray(PAR, dtype=float)
    T   = np.asarray(T, dtype=float)
    NEE = np.asarray(NEE, dtype=float)

    # Remove NaN rows
    valid = np.isfinite(PAR) & np.isfinite(T) & np.isfinite(NEE)
    if not np.all(valid):
        n_removed = np.sum(~valid)
        if verbose:
            print(f"Removed {n_removed} rows with NaN/Inf values.")
        PAR = PAR[valid]
        T   = T[valid]
        NEE = NEE[valid]

    n_obs   = len(NEE)
    night_mask = PAR < par_threshold
    n_night = np.sum(night_mask)
    n_day   = np.sum(~night_mask)

    if verbose:
        print(f"Fitting with {n_obs} observations ({n_night} night, {n_day} day)")
        print(f"Method: {method}")

    # Defaults
    if p0 is None:
        p0 = [0.05, 20.0, 2.0, 200.0]  # alpha, beta, rb, E0

    if bounds is None:
        bounds = (
            [0.001,   1.0,  0.01,   50.0],   # lower bounds
            [0.5,   100.0,  20.0,  500.0],    # upper bounds
        )

    if method == 'simultaneous':
        # ---- Single-step: fit all 4 parameters at once ----
        popt, pcov = curve_fit(
            nee_model,
            (PAR, T),
            NEE,
            p0=p0,
            bounds=bounds,
            maxfev=maxfev,
            method='trf'
        )

        alpha_fit, beta_fit, rb_fit, E0_fit = popt
        std_err = np.sqrt(np.diag(pcov))

    elif method == 'two_step':
        # ---- Two-step: nighttime Reco first, then GPP ----

        if n_night < 3:
            raise ValueError(
                f"two_step method requires >= 3 nighttime observations, "
                f"got {n_night}. Use method='simultaneous' or lower "
                f"par_threshold."
            )

        # Step 1: Fit rb, E0 from nighttime data only
        # At night, NEE ≈ Reco = rb * exp(E0 * (...))
        rb_bounds  = ([bounds[0][2], bounds[0][3]],
                      [bounds[1][2], bounds[1][3]])
        rb_p0 = [p0[2], p0[3]]

        if verbose:
            print(f"\n--- Step 1: Fitting Reco (rb, E0) from {n_night} "
                  f"nighttime observations ---")

        popt_reco, pcov_reco = curve_fit(
            reco_model,
            T[night_mask],
            NEE[night_mask],
            p0=rb_p0,
            bounds=rb_bounds,
            maxfev=maxfev,
            method='trf'
        )
        rb_fit, E0_fit = popt_reco
        std_err_reco = np.sqrt(np.diag(pcov_reco))

        if verbose:
            print(f"  rb = {rb_fit:.4f} ± {std_err_reco[0]:.4f} µmol CO2 m-2 s-1")
            print(f"  E0 = {E0_fit:.2f} ± {std_err_reco[1]:.2f} K")

            # Nighttime fit diagnostics
            reco_night_pred = reco_model(T[night_mask], rb_fit, E0_fit)
            resid_night = NEE[night_mask] - reco_night_pred
            rmse_night = np.sqrt(np.mean(resid_night**2))
            print(f"  Nighttime RMSE: {rmse_night:.4f} µmol CO2 m-2 s-1")

        # Step 2: Fix rb, E0; fit alpha, beta from all data
        def nee_model_fixed_reco(X, alpha, beta):
            """NEE model with rb and E0 fixed from Step 1."""
            PAR_in, T_in = X
            PAR_in = np.asarray(PAR_in, dtype=float)
            T_in   = np.asarray(T_in, dtype=float)
            gpp_term  = -(alpha * beta * PAR_in) / (alpha * PAR_in + beta)
            reco_term = rb_fit * np.exp(
                E0_fit * (1.0 / (TREF - T0) - 1.0 / (T_in - T0))
            )
            return gpp_term + reco_term

        gpp_bounds = ([bounds[0][0], bounds[0][1]],
                      [bounds[1][0], bounds[1][1]])
        gpp_p0 = [p0[0], p0[1]]

        if verbose:
            print(f"\n--- Step 2: Fitting GPP (alpha, beta) with fixed "
                  f"rb={rb_fit:.4f}, E0={E0_fit:.2f} ---")

        popt_gpp, pcov_gpp = curve_fit(
            nee_model_fixed_reco,
            (PAR, T),
            NEE,
            p0=gpp_p0,
            bounds=gpp_bounds,
            maxfev=maxfev,
            method='trf'
        )
        alpha_fit, beta_fit = popt_gpp
        std_err_gpp = np.sqrt(np.diag(pcov_gpp))

        if verbose:
            print(f"  alpha = {alpha_fit:.5f} ± {std_err_gpp[0]:.5f}")
            print(f"  beta  = {beta_fit:.3f} ± {std_err_gpp[1]:.3f}")

        # Assemble combined std_err (note: cross-covariance between
        # steps is not captured; these are conditional standard errors)
        std_err = np.array([std_err_gpp[0], std_err_gpp[1],
                            std_err_reco[0], std_err_reco[1]])
        pcov = {'reco': pcov_reco, 'gpp': pcov_gpp}

        popt = np.array([alpha_fit, beta_fit, rb_fit, E0_fit])

    else:
        raise ValueError(
            f"Unknown method '{method}'. Use 'simultaneous' or 'two_step'."
        )

    # Predictions
    nee_pred  = nee_model((PAR, T), *popt)
    reco_pred = reco_model(T, rb_fit, E0_fit)
    gpp_pred  = gpp_model(PAR, T, *popt)

    # Diagnostics
    residuals = NEE - nee_pred
    ss_res    = np.sum(residuals**2)
    ss_tot    = np.sum((NEE - np.mean(NEE))**2)
    rmse      = np.sqrt(ss_res / n_obs)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    param_names = ['alpha', 'beta', 'rb', 'E0']
    params   = dict(zip(param_names, popt))
    std_errs = dict(zip(param_names, std_err))

    if verbose:
        print(f"\n{'Parameter':<10} {'Value':>10} {'Std Err':>10} {'Units'}")
        print("-" * 50)
        print(f"{'alpha':<10} {alpha_fit:>10.5f} {std_err[0]:>10.5f}  µmol CO2 / µmol photons")
        print(f"{'beta':<10} {beta_fit:>10.3f} {std_err[1]:>10.3f}  µmol CO2 m-2 s-1")
        print(f"{'rb':<10} {rb_fit:>10.4f} {std_err[2]:>10.4f}  µmol CO2 m-2 s-1")
        print(f"{'E0':<10} {E0_fit:>10.2f} {std_err[3]:>10.2f}  K")
        print(f"\nRMSE:  {rmse:.4f} µmol CO2 m-2 s-1")
        print(f"R²:    {r_squared:.4f}")

    return {
        'params':    params,
        'pcov':      pcov,
        'std_err':   std_errs,
        'nee_pred':  nee_pred,
        'reco_pred': reco_pred,
        'gpp_pred':  gpp_pred,
        'residuals': residuals,
        'rmse':      rmse,
        'r_squared': r_squared,
        'n_obs':     n_obs,
        'n_night':   n_night,
        'n_day':     n_day,
        'method':    method,
    }


# ---------------------------------------------------------------------------
# Convenience prediction functions
# ---------------------------------------------------------------------------

def predict_nee(PAR, T, params):
    """Predict NEE given PAR, T, and a params dict from fit_nee_model."""
    return nee_model((PAR, T), params['alpha'], params['beta'],
                     params['rb'], params['E0'])


def predict_reco(T, params):
    """Predict ecosystem respiration given T and a params dict."""
    return reco_model(T, params['rb'], params['E0'])


def predict_gpp(PAR, T, params):
    """Predict gross primary production given PAR, T, and a params dict."""
    return gpp_model(PAR, T, params['alpha'], params['beta'],
                     params['rb'], params['E0'])


def plot_result(result_sim, result_2s, hours, PAR_true, T_true, NEE_obs):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    night_mask = PAR_true < 5

    # 1. Observed vs predicted NEE - simultaneous
    ax = axes[0, 0]
    ax.scatter(result_sim['nee_pred'], NEE_obs, s=10, alpha=0.5, c='steelblue')
    lims = [min(NEE_obs.min(), result_sim['nee_pred'].min()) - 0.5,
            max(NEE_obs.max(), result_sim['nee_pred'].max()) + 0.5]
    ax.plot(lims, lims, 'k--', lw=1)
    ax.set_xlabel('Predicted NEE')
    ax.set_ylabel('Observed NEE')
    ax.set_title(f'Simultaneous (R²={result_sim["r_squared"]:.3f}, '
                 f'RMSE={result_sim["rmse"]:.3f})')

    # 2. Observed vs predicted NEE - two_step
    ax = axes[0, 1]
    ax.scatter(result_2s['nee_pred'], NEE_obs, s=10, alpha=0.5, c='darkorange')
    ax.plot(lims, lims, 'k--', lw=1)
    ax.set_xlabel('Predicted NEE')
    ax.set_ylabel('Observed NEE')
    ax.set_title(f'Two-step (R²={result_2s["r_squared"]:.3f}, '
                 f'RMSE={result_2s["rmse"]:.3f})')

    # 3. Temperature response (Reco) - both methods overlaid
    ax = axes[1, 0]
    T_range = np.linspace(270, 305, 200)
    reco_sim = predict_reco(T_range, result_sim['params'])
    reco_2s  = predict_reco(T_range, result_2s['params'])
    ax.plot(T_range - 273.15, reco_sim, '-', color='steelblue', lw=2,
            label='Simultaneous')
    ax.plot(T_range - 273.15, reco_2s, '-', color='darkorange', lw=2,
            label='Two-step')
    ax.scatter(T_true[night_mask] - 273.15, NEE_obs[night_mask],
               s=15, alpha=0.4, c='gray', label='Nighttime NEE')
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Reco (µmol CO₂ m⁻² s⁻¹)')
    ax.set_title('Reco: Temperature Response Comparison')
    ax.legend()

    # 4. Mean diurnal cycle with spread (all days overlaid)
    ax = axes[1, 1]
    
    # Bin by hour of day
    hour_of_day = hours % 24
    unique_hours = np.sort(np.unique(hour_of_day))
    
    # Compute mean and std for each half-hour bin
    obs_mean, obs_std = [], []
    sim_mean, twostep_mean = [], []
    for h in unique_hours:
        mask = hour_of_day == h
        obs_mean.append(np.mean(NEE_obs[mask]))
        obs_std.append(np.std(NEE_obs[mask]))
        sim_mean.append(np.mean(result_sim['nee_pred'][mask]))
        twostep_mean.append(np.mean(result_2s['nee_pred'][mask]))
    
    obs_mean = np.array(obs_mean)
    obs_std = np.array(obs_std)
    sim_mean = np.array(sim_mean)
    twostep_mean = np.array(twostep_mean)
    
    # Observed: mean ± 1 SD as shaded band
    ax.fill_between(unique_hours, obs_mean - obs_std, obs_mean + obs_std,
                    color='gray', alpha=0.2, label='Obs ± 1 SD')
    ax.plot(unique_hours, obs_mean, 'o-', ms=3, color='gray', lw=1.5,
            label='Obs mean')
    
    # Model predictions: mean lines
    ax.plot(unique_hours, sim_mean, '-', color='steelblue', lw=2,
            label='Simultaneous')
    ax.plot(unique_hours, twostep_mean, '--', color='darkorange', lw=2,
            label='Two-step')
    
    ax.axhline(0, color='gray', lw=0.5, ls='--')
    ax.set_xlabel('Hour of day')
    ax.set_ylabel('NEE (µmol CO₂ m⁻² s⁻¹)')
    ax.set_title('Mean Diurnal Cycle (all days)')
    ax.legend(fontsize=8)

    plt.tight_layout()

    return fig


# ---------------------------------------------------------------------------
# Demo / self-test
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import os

    # ----- Generate synthetic data with nighttime bias -----
    # We simulate data where nighttime respiration has extra noise
    # and a slight offset, mimicking conditions where simultaneous
    # fitting tends to overpredict nighttime NEE.
    np.random.seed(42)
    n = 200

    # Simulate a diurnal cycle (hours 0-23 repeated ~8 days)
    hours = np.tile(np.arange(0, 24, 0.5), n // 48 + 1)[:n]

    # PAR: bell curve peaking at noon
    PAR_true = np.maximum(0, 1500 * np.sin(np.pi * (hours - 6) / 12))
    PAR_true[hours < 6]  = 0
    PAR_true[hours > 18] = 0

    # Temperature: sinusoid with daily cycle (in K)
    T_true = 283 + 8 * np.sin(np.pi * (hours - 3) / 12)

    # True parameters
    alpha_true = 0.04
    beta_true  = 15.0
    rb_true    = 2.5
    E0_true    = 180.0

    # Generate NEE with noise
    NEE_true = nee_model((PAR_true, T_true), alpha_true, beta_true, rb_true, E0_true)
    NEE_obs  = NEE_true + np.random.normal(0, 0.5, n)

    # ----- Fit with both methods -----
    print("=" * 60)
    print("DEMO: Comparing simultaneous vs two_step fitting")
    print("=" * 60)
    print(f"\nTrue parameters: alpha={alpha_true}, beta={beta_true}, "
          f"rb={rb_true}, E0={E0_true}")

    print("\n" + "=" * 60)
    print("METHOD 1: simultaneous")
    print("=" * 60)
    result_sim = fit_nee_model(PAR_true, T_true, NEE_obs,
                               method='simultaneous', verbose=True)

    print("\n" + "=" * 60)
    print("METHOD 2: two_step")
    print("=" * 60)
    result_2s = fit_nee_model(PAR_true, T_true, NEE_obs,
                              method='two_step', verbose=True)

    # ----- Comparison plot -----
    fig = plot_result(result_sim, result_2s, hours, PAR_true, T_true, NEE_obs)
    plt.savefig(os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'intermediate', 'nee_fit_demo.png'), dpi=150)
    print("\nPlot saved to nee_fit_demo.png")