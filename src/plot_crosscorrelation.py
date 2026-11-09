"""Cross-correlation between atmospheric CO2 and GMST, across the same four windows.

plot_gmst.py, plot_co2.py, plot_holocene.py and plot_industrial.py plot the two
variables against each other's absolute value; plot_derivatives.py plots their
rates of change. Neither says anything about which one moves first. This script
does: for each of the four windows, both records are put on one common,
uniformly discretised timeline, a zero-phase Butterworth high-pass filter
removes the shared secular drift from each so what is left is the faster
oscillation, and the two
residuals are cross-correlated at a range of lags. The lag at which the
correlation peaks is the answer to "does CO2 lead GMST, or the other way round,
and by how much".

SIGN CONVENTION, and it matters: correlate(co2, temp, lag) pairs co2 at time t
with temp at time t + lag. A peak at a POSITIVE lag means the temp series
several steps later looks like the CO2 series looks now, i.e. CO2 LEADS temp by
that lag. A peak at a NEGATIVE lag means the reverse: temp leads, CO2 follows.
This was checked against a synthetic signal with a known shift before trusting
it on real data.

CAVEAT, and it applies more severely here than anywhere else in this set: a
cross-correlation computed at a sampling period S cannot resolve a lag shorter
than S. The best-known result of this kind in the literature - CO2 lagging
Antarctic temperature by some 600-800 years across the last deglaciation
(Caillon et al. 2003) - sits close to the 800 kyr window's 100-year sampling
period's own resolving power; whether it actually shows up here or not, what
that window's cross-correlation mainly speaks to is the phase relationship of
the ~100 kyr glacial cycle as a whole, not any short lag inside a single
termination. Read every lag estimate below against the sampling period it was
computed on, not as a precise date.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt

POLAR_AMPLIFICATION = 2.0  # see plot_co2.py

# ==========================================
# 1. LOAD THE SAME RECORDS AS THE OTHER FIGURES
# ==========================================
df_phan = pd.read_csv("data/PhanDA_GMSTandCO2_percentiles.csv")
phan_age_yr = df_phan["AverageAge"].to_numpy() * 1e6

df_ice_co2 = pd.read_csv(
    "data/antarctica2015co2composite-noaa.txt", comment="#", sep="\t"
).rename(columns={"age_gas_calBP": "Age_BP", "co2_ppm": "co2"})[["Age_BP", "co2"]].dropna()

edc = []
with open("data/edc3deuttemp2007.txt", encoding="utf-8", errors="ignore") as f:
    for line in f:
        tokens = line.split()
        if len(tokens) <= 4:
            continue
        try:
            edc.append([float(tokens[2]), float(tokens[4])])
        except ValueError:
            continue
df_edc_temp = pd.DataFrame(edc, columns=["Age_BP", "temp"])
df_edc_temp["temp"] /= POLAR_AMPLIFICATION

df_holo = pd.read_csv(
    "data/temp12k_allmethods_percentiles.csv", skipinitialspace=True
).rename(columns={"ages": "Age_BP", "global_median": "temp"})[["Age_BP", "temp"]]

df_mlo = pd.read_csv(
    "data/mauna-loa-co2-mm.txt", comment="#", sep=r"\s+", header=None,
    names=["year", "month", "decimal_date", "co2", "deseasonalized", "ndays", "sdev", "unc"],
)
# Monthly, not the annual summary series: 'Time' comes as 'YYYY-MM', converted
# to a decimal year at the middle of the month, matching the convention Mauna
# Loa's own decimal_date already uses.
df_had = pd.read_csv("data/hadcrut5-monthly-global.csv").rename(
    columns={"Anomaly (deg C)": "temp"}
)
_y, _m = df_had["Time"].str.split("-", expand=True).astype(int).to_numpy().T
df_had["year"] = _y + (_m - 0.5) / 12.0
co2_instrumental = pd.concat([
    df_ice_co2.assign(year=1950 - df_ice_co2["Age_BP"]).query("1800 <= year < 1958")[["year", "co2"]],
    df_mlo[["decimal_date", "co2"]].rename(columns={"decimal_date": "year"}),
])


# ==========================================
# 2. COMMON GRID, DETREND, CROSS-CORRELATE
# ==========================================
def common_grid(year1, val1, year2, val2, step_yr):
    """Interpolate two irregularly sampled series onto one shared, uniformly
    discretised timeline of calendar years, over the range where both have data.

    A cross-correlation needs the two series to be indexed the same way; the
    source records here never are (different archives, different sampling).
    Linear interpolation is used, same as elsewhere in this set (see
    resample_uniform in lowpass.py) - it adds no information, it only puts both
    records on a common footing.
    """
    year1, val1 = np.asarray(year1, dtype=float), np.asarray(val1, dtype=float)
    year2, val2 = np.asarray(year2, dtype=float), np.asarray(val2, dtype=float)
    o1, o2 = np.argsort(year1), np.argsort(year2)
    year1, val1 = year1[o1], val1[o1]
    year2, val2 = year2[o2], val2[o2]

    lo = max(year1.min(), year2.min())
    hi = min(year1.max(), year2.max())
    grid = np.arange(lo, hi, step_yr)
    return grid, np.interp(grid, year1, val1), np.interp(grid, year2, val2)


def highpass_residual(values, cutoff_period_yr, samples_per_year):
    """Zero-phase Butterworth high-pass: removes the drift, passes the oscillation.

    A moving average is a poor low-pass filter: its frequency response is a
    sinc, with a first sidelobe only 13 dB down, so it never really stops
    passing anything and leaks secular drift into the residual it is supposed
    to remove. This is the direct alternative: an actual high-pass filter,
    order 4, run forwards and backwards with filtfilt so it introduces no
    phase shift - the same zero-phase construction lowpass.py uses elsewhere
    in this set, applied here as a high-pass instead of a low-pass because it
    is the residual, not the trend, that is wanted.

    Padding uses filtfilt's own default (odd reflection about each end
    point), not the linear-trend extrapolation lowpass.py uses for the
    instrumental figures: that scheme exists to protect a "live" right edge
    from an assumption about what comes after it, which does not apply here -
    both ends of every window in this analysis are equally fixed, equally
    historical, and there is nothing to protect them from.

    cutoff_period_yr is chosen per window to be well above the timescale of
    variability being examined (see the per-window calls below), so this
    removes the shared secular drift without removing the oscillation or
    excursion whose lead-lag relationship is the point of the exercise.
    """
    values = np.asarray(values, dtype=float)
    nyquist_per_year = samples_per_year / 2.0
    normalised = (1.0 / cutoff_period_yr) / nyquist_per_year
    b, a = butter(4, normalised, btype="highpass")
    padlen = min(len(values) - 1, 3 * int(cutoff_period_yr * samples_per_year))
    return filtfilt(b, a, values, padlen=padlen)


def lowpass_band(values, cutoff_period_yr, samples_per_year):
    """Zero-phase Butterworth low-pass, same construction as highpass_residual()
    above but passing periods slower than cutoff_period_yr rather than faster.

    Used to chain onto highpass_residual() and make a band-pass: high-pass at
    a long cutoff to remove the secular trend, then low-pass at a short one to
    also remove the seasonal cycle and anything faster, leaving the smooth
    interannual wiggle in between - see the since-1850 window below, the one
    window fine enough to resolve a seasonal cycle in the first place.
    """
    values = np.asarray(values, dtype=float)
    nyquist_per_year = samples_per_year / 2.0
    normalised = (1.0 / cutoff_period_yr) / nyquist_per_year
    b, a = butter(4, normalised, btype="lowpass")
    padlen = min(len(values) - 1, 3 * int(cutoff_period_yr * samples_per_year))
    return filtfilt(b, a, values, padlen=padlen)


def normalized_xcorr(x, y, max_lag):
    """Pearson correlation between x and y at each lag from -max_lag to +max_lag.

    corr(lag) pairs x[t] with y[t + lag]. A peak at positive lag means y several
    steps later resembles x now - x LEADS y. A peak at negative lag means the
    reverse. Verified against a synthetic shifted signal before use here.
    """
    x = (np.asarray(x, dtype=float) - np.nanmean(x))
    y = (np.asarray(y, dtype=float) - np.nanmean(y))
    x, y = x / np.nanstd(x), y / np.nanstd(y)
    n = len(x)
    lags = np.arange(-max_lag, max_lag + 1)
    corr = np.full(len(lags), np.nan)
    for i, lag in enumerate(lags):
        a = x[:n - lag] if lag > 0 else x[-lag:]
        b = y[lag:] if lag > 0 else y[:n + lag]
        valid = ~(np.isnan(a) | np.isnan(b))
        if valid.sum() > max_lag:
            corr[i] = np.mean(a[valid] * b[valid])
    return lags, corr


# Per window: (name, co2 (year, value), temp (year, value), sampling period
# yr, high-pass cutoff in sampling periods, max lag in sampling periods, lag
# axis divisor+label,
# title)
#
# Each detrend_steps was picked by sweeping a range and keeping the smallest
# one at which the peak lag and its correlation stop moving: 100-200 Myr for
# the Phanerozoic (below 100 Myr the residual still carries enough leftover
# secular drift to throw up a spurious secondary peak - this
# was seen at 50 Myr), 200-800 kyr for the ice cores, and
# 1,000-3,000 yr for the Holocene, all essentially flat once past those
# thresholds. The since-1850 window uses the same 20-year cutoff as the trend
# line drawn throughout the industrial-era figures, rather than a value
# chosen by this sweep, so that "trend" means the same timescale here as
# everywhere else in this set.
WINDOWS = [
    dict(
        name="xcorr_485Ma", step_yr=5_000_000, detrend_steps=30, max_lag=20,
        lag_div=1e6, lag_label="Time shift τ (Myr, positive = CO$_2$ leads GMST)",
        title="485 million years (PhanDA, 5 Myr sampling period)",
        co2=(-phan_age_yr, df_phan["CO2_50"].to_numpy()),
        temp=(-phan_age_yr, df_phan["GMST_50"].to_numpy()),
    ),
    dict(
        name="xcorr_800k", step_yr=100, detrend_steps=4_000, max_lag=500,
        lag_div=1_000, lag_label="Time shift τ (kyr, positive = CO$_2$ leads GMST)",
        title="800,000 years (ice cores, 100-year sampling period)",
        co2=(1950 - df_ice_co2["Age_BP"].to_numpy(), df_ice_co2["co2"].to_numpy()),
        temp=(1950 - df_edc_temp["Age_BP"].to_numpy(), df_edc_temp["temp"].to_numpy()),
    ),
    dict(
        name="xcorr_12k", step_yr=100, detrend_steps=20, max_lag=60,
        lag_div=1_000, lag_label="Time shift τ (kyr, positive = CO$_2$ leads GMST)",
        title="12,000 years (Holocene, 100-year sampling period)",
        co2=(1950 - df_ice_co2["Age_BP"].to_numpy(), df_ice_co2["co2"].to_numpy()),
        temp=(1950 - df_holo["Age_BP"].to_numpy(), df_holo["temp"].to_numpy()),
    ),
    dict(
        name="xcorr_200y", step_yr=1.0 / 12, detrend_steps=240, max_lag=600,
        lag_div=1, lag_label="Time shift τ (years, positive = CO$_2$ leads GMST)",
        title="Since 1850 (instrumental, 1-month sampling period, interannual band)",
        co2=(co2_instrumental["year"].to_numpy(), co2_instrumental["co2"].to_numpy()),
        temp=(df_had["year"].to_numpy(), df_had["temp"].to_numpy()),
        # This window alone is fine enough to resolve a seasonal cycle, so it
        # alone needs a second filter: a low-pass at seasonal_cutoff_yr chained
        # onto the high-pass above turns this into a band-pass, keeping the
        # smooth interannual wiggle and removing the seasonal cycle along with
        # the secular trend. The seasonal cycle itself is not discarded, only
        # moved: it is what the dedicated section below isolates and
        # cross-correlates in its own right, at its own timescale.
        seasonal_cutoff_yr=2,
    ),
]

print(f"{'window':10s} {'samp. period':>12s} {'cutoff':>10s} {'peak lag':>10s} {'peak r':>8s}")
results = []  # collected for the four-panel summary figure in section 3 below
for spec in WINDOWS:
    (co2_yr, co2_val), (temp_yr, temp_val) = spec["co2"], spec["temp"]
    grid, co2_g, temp_g = common_grid(co2_yr, co2_val, temp_yr, temp_val, spec["step_yr"])
    cutoff_period_yr = spec["detrend_steps"] * spec["step_yr"]
    samples_per_year = 1.0 / spec["step_yr"]
    co2_res = highpass_residual(co2_g, cutoff_period_yr, samples_per_year)
    temp_res = highpass_residual(temp_g, cutoff_period_yr, samples_per_year)
    if "seasonal_cutoff_yr" in spec:
        co2_res = lowpass_band(co2_res, spec["seasonal_cutoff_yr"], samples_per_year)
        temp_res = lowpass_band(temp_res, spec["seasonal_cutoff_yr"], samples_per_year)
    lags, corr = normalized_xcorr(co2_res, temp_res, spec["max_lag"])
    lag_years = lags * spec["step_yr"]

    peak_i = np.nanargmax(corr)
    peak_lag_yr, peak_r = lag_years[peak_i], corr[peak_i]
    n_eff = len(grid) - spec["max_lag"]
    floor = 2.0 / np.sqrt(n_eff)
    print(f"{spec['name']:10s} {spec['step_yr']:>10.4g} yr {cutoff_period_yr:>9.4g} yr "
          f"{peak_lag_yr:>+8.0f} yr {peak_r:>8.2f} (floor {floor:.2f}, SNR {peak_r/floor:.1f})")

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.axhline(0.0, color="black", lw=0.8, alpha=0.5)
    ax.axvline(0.0, color="dimgray", lw=0.8, ls="--", alpha=0.7)
    ax.plot(lag_years / spec["lag_div"], corr, color="tab:purple", lw=1.8)
    ax.scatter([peak_lag_yr / spec["lag_div"]], [peak_r], color="tab:purple", s=40, zorder=5)

    # A naive +-2/sqrt(N) noise floor, the textbook large-sample threshold
    # above which two independent white-noise series would rarely correlate
    # by chance; the residuals here are autocorrelated (that is the whole
    # point of them being a real signal), so this understates the true
    # uncertainty and is drawn only as a rough by-eye reference, not a
    # formal test.
    ax.axhline(floor, color="dimgray", lw=0.8, ls=":", alpha=0.6,
               label=f"$\\pm 2/\\sqrt{{n}}$ noise floor ($n \\approx {n_eff}$)")
    ax.axhline(-floor, color="dimgray", lw=0.8, ls=":", alpha=0.6)

    lead = "CO$_2$ leads" if peak_lag_yr > 0 else "GMST leads" if peak_lag_yr < 0 else "in phase"
    # A label offset to the right of the peak runs off the axes if the peak
    # itself already sits well into the right-hand part of the range (as
    # happens for the since-1850 window, where the peak lands close to the
    # +50-year edge); flip the label to the left of the point in that case.
    near_right_edge = (peak_lag_yr / spec["lag_div"]) > 0.4 * (spec["max_lag"] * spec["step_yr"] / spec["lag_div"])
    ha, xoffset = ("right", -10) if near_right_edge else ("left", 10)
    ax.annotate(
        f"peak r={peak_r:.2f} at τ={peak_lag_yr:+.0f} yr\n({lead})",
        xy=(peak_lag_yr / spec["lag_div"], peak_r), xytext=(xoffset, -10),
        textcoords="offset points", ha=ha, va="top",
        fontsize=9, fontweight="bold", color="tab:purple",
    )

    ax.set_xlabel(spec["lag_label"], fontsize=11, fontweight="bold")
    ax.set_ylabel("Cross-correlation, $\\mathrm{CO}_2$ vs GMST residual", fontsize=10, fontweight="bold")
    ax.set_title(spec["title"], fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend(loc="lower right", fontsize=8)

    plt.tight_layout()
    plt.savefig(f"images/co2_vs_temperature_{spec['name']}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote images/co2_vs_temperature_{spec['name']}.png")

    results.append(dict(
        spec=spec, lag_years=lag_years, corr=corr,
        peak_lag_yr=peak_lag_yr, peak_r=peak_r, floor=floor, n_eff=n_eff,
    ))

# ==========================================
# 3. FOUR-PANEL SUMMARY: ALL FOUR WINDOWS AT A GLANCE
# ==========================================
# The four figures above show one window each;
# this one exists to stand alone as
# a lead image - the single frame that shows, without reading a word, that all
# four windows peak close to zero lag and that the peaks get sharper and more
# confident (485 Myr and 200 yr low, 800 kyr and 12 kyr high) towards the two
# middle windows. That contrast is the main result, so it is what a reader
# should see first.
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(11, 10))
for ax, r in zip(axes.flat, results):
    spec = r["spec"]
    ax.axhline(0.0, color="black", lw=0.8, alpha=0.5)
    ax.axvline(0.0, color="dimgray", lw=0.8, ls="--", alpha=0.7)
    ax.plot(r["lag_years"] / spec["lag_div"], r["corr"], color="tab:purple", lw=1.6)
    ax.scatter([r["peak_lag_yr"] / spec["lag_div"]], [r["peak_r"]],
               color="tab:purple", s=35, zorder=5)
    ax.axhline(r["floor"], color="dimgray", lw=0.8, ls=":", alpha=0.6)
    ax.axhline(-r["floor"], color="dimgray", lw=0.8, ls=":", alpha=0.6)

    lead = ("CO$_2$ leads" if r["peak_lag_yr"] > 0
            else "GMST leads" if r["peak_lag_yr"] < 0 else "in phase")
    # Same right-edge check as the single-panel figures above: flip the label
    # to the left of the point rather than let it run off the axes.
    near_right_edge = (r["peak_lag_yr"] / spec["lag_div"]) > 0.4 * (spec["max_lag"] * spec["step_yr"] / spec["lag_div"])
    ha, xoffset = ("right", -8) if near_right_edge else ("left", 8)
    ax.annotate(
        f"r={r['peak_r']:.2f} @ τ={r['peak_lag_yr']:+.0f} yr ({lead})",
        xy=(r["peak_lag_yr"] / spec["lag_div"], r["peak_r"]), xytext=(xoffset, -8),
        textcoords="offset points", ha=ha, va="top",
        fontsize=8, fontweight="bold", color="tab:purple",
    )
    ax.set_xlabel(spec["lag_label"], fontsize=9)
    ax.set_ylabel("Cross-correlation", fontsize=9, fontweight="bold")
    # Wrapped onto two lines: these titles are wider than one subplot column,
    # and at this panel size tight_layout alone cannot stop a long single-line
    # title from overrunning into the neighbouring panel's title.
    ax.set_title(spec["title"].replace(" (", "\n(", 1), fontsize=9.5, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.3)

plt.tight_layout()
fig.subplots_adjust(hspace=0.4, wspace=0.28)
plt.savefig("images/co2_vs_temperature_xcorr_overview.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("wrote images/co2_vs_temperature_xcorr_overview.png")

# ==========================================
# 4. 200-YEAR WINDOW, SEASONAL BAND
# ==========================================
# The interannual xcorr_200y figure above deliberately excludes the seasonal
# cycle - a high-pass at a 20-year cutoff, then a low-pass at 2 years, keeps
# only the smooth wiggle in between. This section is the complement: isolate
# just the seasonal cycle and faster (a high-pass at 2 years, nothing else)
# and cross-correlate that band on its own, over a short lag range where a
# seasonal phase relationship would actually show up. It is a different
# question from the rest of this script - not whether CO2 leads GMST over
# the industrial era, but whether the atmosphere's own seasonal breath and the
# planet's own seasonal temperature cycle are in step - but the same monthly
# sampling period built for the window above answers it for free.
seasonal_spec = next(s for s in WINDOWS if s["name"] == "xcorr_200y")
(co2_yr, co2_val), (temp_yr, temp_val) = seasonal_spec["co2"], seasonal_spec["temp"]
grid, co2_g, temp_g = common_grid(co2_yr, co2_val, temp_yr, temp_val, seasonal_spec["step_yr"])
samples_per_year = 1.0 / seasonal_spec["step_yr"]
SEASONAL_CUTOFF_YR = seasonal_spec["seasonal_cutoff_yr"]
SEASONAL_MAX_LAG_MONTHS = 36  # +-3 years: several seasonal cycles either way

co2_seasonal = highpass_residual(co2_g, SEASONAL_CUTOFF_YR, samples_per_year)
temp_seasonal = highpass_residual(temp_g, SEASONAL_CUTOFF_YR, samples_per_year)
lags, corr = normalized_xcorr(co2_seasonal, temp_seasonal, SEASONAL_MAX_LAG_MONTHS)
lag_years = lags * seasonal_spec["step_yr"]

# This correlation repeats every 12 months by construction (that periodicity
# is the whole point of the figure), so a plain argmax over the full +-36
# month range can just as easily land on a repeat one or two cycles away as
# on the cycle nearest zero lag, and which one wins can flip between runs as
# new monthly data shifts the noise on each near-identical repeat. The
# genuine seasonal phase, read modulo the 12-month period, always falls
# within half a period of zero, so restrict the peak search to +-6 months
# and let only the plot itself show the full range.
near_zero = np.abs(lags) <= 6
peak_i = np.flatnonzero(near_zero)[np.nanargmax(corr[near_zero])]
peak_lag_yr, peak_r = lag_years[peak_i], corr[peak_i]
n_eff = len(grid) - SEASONAL_MAX_LAG_MONTHS
floor = 2.0 / np.sqrt(n_eff)
print(f"\nxcorr_200y_seasonal   {seasonal_spec['step_yr']:.4g} yr   cutoff {SEASONAL_CUTOFF_YR} yr   "
      f"peak {peak_lag_yr*12:+.0f} mo   r {peak_r:.2f}  (floor {floor:.2f}, SNR {peak_r/floor:.1f})")

fig, ax = plt.subplots(figsize=(9, 5))
ax.axhline(0.0, color="black", lw=0.8, alpha=0.5)
ax.axvline(0.0, color="dimgray", lw=0.8, ls="--", alpha=0.7)
ax.plot(lag_years * 12, corr, color="tab:green", lw=1.8)
ax.scatter([peak_lag_yr * 12], [peak_r], color="tab:green", s=40, zorder=5)
ax.axhline(floor, color="dimgray", lw=0.8, ls=":", alpha=0.6,
           label=f"$\\pm 2/\\sqrt{{n}}$ noise floor ($n \\approx {n_eff}$)")
ax.axhline(-floor, color="dimgray", lw=0.8, ls=":", alpha=0.6)

lead = "CO$_2$ leads" if peak_lag_yr > 0 else "GMST leads" if peak_lag_yr < 0 else "in phase"
near_right_edge = (peak_lag_yr * 12) > 0.4 * SEASONAL_MAX_LAG_MONTHS
ha, xoffset = ("right", -10) if near_right_edge else ("left", 10)
ax.annotate(
    f"peak r={peak_r:.2f} at τ={peak_lag_yr*12:+.0f} months\n({lead})",
    xy=(peak_lag_yr * 12, peak_r), xytext=(xoffset, -10),
    textcoords="offset points", ha=ha, va="top",
    fontsize=9, fontweight="bold", color="tab:green",
)
ax.set_xlabel("Time shift τ (months, positive = CO$_2$ leads GMST)", fontsize=11, fontweight="bold")
ax.set_ylabel("Cross-correlation, $\\mathrm{CO}_2$ vs GMST seasonal band", fontsize=10, fontweight="bold")
ax.set_title("Since 1850 (instrumental, 1-month sampling period, seasonal band)", fontsize=11)
ax.grid(True, linestyle="--", alpha=0.3)
ax.legend(loc="lower right", fontsize=8)

plt.tight_layout()
plt.savefig("images/co2_vs_temperature_xcorr_200y_seasonal.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("wrote images/co2_vs_temperature_xcorr_200y_seasonal.png")
