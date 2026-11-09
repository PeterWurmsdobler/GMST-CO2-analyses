"""A single log-time overview of the whole 485-million-year record, rates only.

Companion to plot_overview.py, which stacks each variable's value directly
above its own rate on one log, reversed time axis running from 485 million
years down to 10 years before present - see that script's docstring for the
full reasoning behind the axis and the splice points. This script reuses the
same loading and rate computation, but drops the value panels: this is the
figure to open a discussion of rates of change with, where the value curves
have already been shown once and do not need repeating.

Because a derivative is invariant to a constant offset in its input, this
script does NOT need plot_overview.py's temperature-baseline harmonisation
(PhanDA absolute vs ice-core/Holocene/instrumental anomaly): d(T + c)/dt =
dT/dt for any constant c, so the raw, un-shifted temperature series from each
source is differentiated directly. Only CO2 and temperature VALUES, not
rates, ever needed that harmonisation step.
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lowpass import lowpass, lowpass_rate_per_century, settled  # noqa: E402

REFERENCE_YEAR = 2026  # "present", for the age-before-present log axis
POLAR_AMPLIFICATION = 2.0

# Handoff ages (years before REFERENCE_YEAR), matching the boundaries already
# used throughout the rest of this set's four windows.
AGE_ICECORE_START = 800_000    # PhanDA -> ice-core CO2 / EDC temperature
AGE_HOLOCENE_START = 12_000    # EDC -> Temp12k (CO2 stays on the ice-core composite); Temp12k's own native upper limit
AGE_INSTRUMENTAL_START = 100   # Temp12k / ice-core -> HadCRUT5 / Mauna Loa; HadCRUT5's own native start (1850)
AGE_RIGHT_EDGE = 10            # right edge of the whole overview

WINDOW_PHANEROZOIC_YR = 5_000_000
WINDOW_ICECORE_YR = 100
WINDOW_HOLOCENE_YR = 100
CUTOFF_PERIOD_YR = 20
# Low-pass cutoffs the rate segments below are filtered to before being
# differentiated - see filtered_rate_per_century() - matched to the same
# values plot_derivatives.py's headline rate figures use for these two windows.
# Each series filtered to its own evidence-based cutoff - see
# plot_derivatives.py's module docstring (RATE-FILTER CUTOFFS) for the
# variance-retained and physical-smoothing justification behind each number.
ICECORE_CO2_CUTOFF_YR = 700
ICECORE_TEMP_CUTOFF_YR = 250
HOLOCENE_RATE_CUTOFF_YR = 220


def rate_per_century(age_yr, values, window_yr):
    """Signed d(value)/dt in units per century, a numerical derivative over a sampling period of window_yr.

    Identical to the function of the same name in plot_derivatives.py: ages run
    backwards in time, so the sign is flipped to give a forward-time rate.
    """
    order = np.argsort(age_yr)
    age_yr, values = np.asarray(age_yr)[order], np.asarray(values)[order]
    grid = np.arange(age_yr.min(), age_yr.max(), window_yr)
    gridded = np.interp(grid, age_yr, values)
    rate_per_year = -np.gradient(gridded, grid)
    return grid, rate_per_year * 100.0


def filtered_rate_per_century(age_yr, values, step_yr, cutoff_yr):
    """Zero-phase low-pass, then numerical differentiation, in units per century.

    Identical in intent to the function of the same name in plot_derivatives.py,
    including the fix for the same bug found there: the filter needs several
    samples per cutoff period to be well-conditioned, which step_yr alone does
    not give it (a 220- or 250-year cutoff on a 100-year grid sits at a
    normalised frequency of 0.91 or 0.80, uncomfortably close to that grid's
    own Nyquist limit of 1.0 - confirmed to matter empirically, not just in
    principle, since refining the grid changed one window's fastest-rate
    reading by a factor of two and flipped its sign). So the filter is
    designed and applied on a finer working grid, twenty samples per cutoff
    period, and only the result is resampled back onto step_yr, which is
    still exactly what the text and figures call this window's sampling
    period.
    """
    order = np.argsort(age_yr)
    age_yr, values = np.asarray(age_yr)[order], np.asarray(values)[order]
    fine_step = min(step_yr, cutoff_yr / 20.0)
    fine_grid = np.arange(age_yr.min(), age_yr.max(), fine_step)
    fine_gridded = np.interp(fine_grid, age_yr, values)
    fine_smoothed = lowpass(fine_gridded, cutoff_yr, samples_per_year=1.0 / fine_step)
    fine_rate_per_century = -np.gradient(fine_smoothed, fine_grid) * 100.0
    grid = np.arange(age_yr.min(), age_yr.max(), step_yr)
    rate_per_century = np.interp(grid, fine_grid, fine_rate_per_century)
    return grid, rate_per_century


# ==========================================
# 1. LOAD THE SAME RECORDS AS THE OTHER FIGURES
# ==========================================
df_phan = pd.read_csv("data/PhanDA_GMSTandCO2_percentiles.csv")
phan_year = 1950 - df_phan["AverageAge"].to_numpy() * 1e6
phan_age = REFERENCE_YEAR - phan_year
phan_co2 = df_phan["CO2_50"].to_numpy()
phan_gmst = df_phan["GMST_50"].to_numpy()  # absolute degC - fine, see docstring

df_ice_co2 = pd.read_csv(
    "data/antarctica2015co2composite-noaa.txt", comment="#", sep="\t"
).rename(columns={"age_gas_calBP": "Age_BP", "co2_ppm": "co2"})[["Age_BP", "co2"]].dropna()
ice_year = 1950 - df_ice_co2["Age_BP"].to_numpy()
ice_age = REFERENCE_YEAR - ice_year
ice_co2 = df_ice_co2["co2"].to_numpy()

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
df_edc = pd.DataFrame(edc, columns=["Age_BP", "temp"])
edc_year = 1950 - df_edc["Age_BP"].to_numpy()
edc_age = REFERENCE_YEAR - edc_year
edc_temp = df_edc["temp"].to_numpy() / POLAR_AMPLIFICATION

df_holo = pd.read_csv(
    "data/temp12k_allmethods_percentiles.csv", skipinitialspace=True
).rename(columns={"ages": "Age_BP", "global_median": "temp"})[["Age_BP", "temp"]]
holo_year = 1950 - df_holo["Age_BP"].to_numpy()
holo_age = REFERENCE_YEAR - holo_year
holo_temp = df_holo["temp"].to_numpy()

df_mlo = pd.read_csv(
    "data/mauna-loa-co2-mm.txt", comment="#", sep=r"\s+", header=None,
    names=["year", "month", "decimal_date", "co2", "deseasonalized", "ndays", "sdev", "unc"],
)
mlo_age = REFERENCE_YEAR - df_mlo["decimal_date"].to_numpy()
mlo_co2 = df_mlo["co2"].to_numpy()

# Monthly, not the annual summary series: 'Time' comes as 'YYYY-MM', converted
# to a decimal year at the middle of the month, matching the convention Mauna
# Loa's own decimal_date already uses.
df_had = pd.read_csv("data/hadcrut5-monthly-global.csv").rename(
    columns={"Anomaly (deg C)": "temp"}
)
_y, _m = df_had["Time"].str.split("-", expand=True).astype(int).to_numpy().T
df_had["year"] = _y + (_m - 0.5) / 12.0

# ==========================================
# 2. COMPUTE THE RATES, ONE SEGMENT PER SOURCE, HANDED OFF AT THE SAME
#    BOUNDARIES AS EVERY OTHER FIGURE IN THIS SET
# ==========================================
co2_rate_segments, temp_rate_segments = [], []

grid, rate = rate_per_century(phan_age[phan_age > AGE_ICECORE_START],
                               phan_co2[phan_age > AGE_ICECORE_START], WINDOW_PHANEROZOIC_YR)
co2_rate_segments.append((grid, rate))
grid, rate = rate_per_century(phan_age[phan_age > AGE_ICECORE_START],
                               phan_gmst[phan_age > AGE_ICECORE_START], WINDOW_PHANEROZOIC_YR)
temp_rate_segments.append((grid, rate))

# Each mask below caps only the OLDER edge of its window, not the younger
# one. ice_co2 and edc_temp are each one continuous record running well past
# every one of these windows' younger boundary (ice_co2 to -51 yr BP,
# edc_temp to 38 yr BP) - truncating the input there before filtering would
# manufacture exactly the kind of edge settled() already exists to guard
# against, only at an artificial boundary of this figure's own choosing
# rather than a genuine end of the data. Letting the filter see the whole
# continuous record on its older side too, and applying each window's
# younger boundary only when selecting what to keep, means settled() trims
# relative to where the archive actually stops, not where this figure hands
# off to the next window.
ice_mask = ice_age <= AGE_ICECORE_START
grid, rate = filtered_rate_per_century(ice_age[ice_mask], ice_co2[ice_mask],
                                       WINDOW_ICECORE_YR, ICECORE_CO2_CUTOFF_YR)
keep = settled(grid, ICECORE_CO2_CUTOFF_YR) & (grid > AGE_HOLOCENE_START)
co2_rate_segments.append((grid[keep], rate[keep]))
edc_mask = edc_age <= AGE_ICECORE_START
grid, rate = filtered_rate_per_century(edc_age[edc_mask], edc_temp[edc_mask],
                                       WINDOW_ICECORE_YR, ICECORE_TEMP_CUTOFF_YR)
keep = settled(grid, ICECORE_TEMP_CUTOFF_YR) & (grid > AGE_HOLOCENE_START)
temp_rate_segments.append((grid[keep], rate[keep]))

ice_mask = ice_age <= AGE_HOLOCENE_START
grid, rate = filtered_rate_per_century(ice_age[ice_mask], ice_co2[ice_mask],
                                       WINDOW_HOLOCENE_YR, HOLOCENE_RATE_CUTOFF_YR)
keep = settled(grid, HOLOCENE_RATE_CUTOFF_YR) & (grid > AGE_INSTRUMENTAL_START)
co2_rate_segments.append((grid[keep], rate[keep]))
holo_mask = holo_age <= AGE_HOLOCENE_START
grid, rate = filtered_rate_per_century(holo_age[holo_mask], holo_temp[holo_mask],
                                       WINDOW_HOLOCENE_YR, HOLOCENE_RATE_CUTOFF_YR)
keep = settled(grid, HOLOCENE_RATE_CUTOFF_YR) & (grid > AGE_INSTRUMENTAL_START)
temp_rate_segments.append((grid[keep], rate[keep]))

# Instrumental: the same 20-year zero-phase filtered derivative used as the
# headline curve in plot_derivatives.py, converted from calendar year to age
# and clipped to AGE_INSTRUMENTAL_START - see plot_overview.py's section 4 for
# why (that reasoning is about the filter's settling behaviour, not about
# baselines, so it applies here unchanged).
co2_instrumental = pd.concat([
    df_ice_co2.assign(year=1950 - df_ice_co2["Age_BP"]).query("1800 <= year < 1958")[["year", "co2"]],
    df_mlo[["decimal_date", "co2"]].rename(columns={"decimal_date": "year"}),
])
grid_yr, rate = lowpass_rate_per_century(co2_instrumental["year"], co2_instrumental["co2"], CUTOFF_PERIOD_YR)
keep = settled(grid_yr, CUTOFF_PERIOD_YR)
age = REFERENCE_YEAR - grid_yr[keep]
in_range = age <= AGE_INSTRUMENTAL_START
co2_rate_segments.append((age[in_range], rate[keep][in_range]))

grid_yr, rate = lowpass_rate_per_century(df_had["year"], df_had["temp"], CUTOFF_PERIOD_YR)
keep = settled(grid_yr, CUTOFF_PERIOD_YR)
age = REFERENCE_YEAR - grid_yr[keep]
in_range = age <= AGE_INSTRUMENTAL_START
temp_rate_segments.append((age[in_range], rate[keep][in_range]))


def concat_segments(segments, lo=AGE_RIGHT_EDGE):
    """Join segments end to end for one continuous plotted line - except where
    they are not actually contiguous. Each segment is already trimmed to what
    its own filter and its own record can support (settled(), industrial-era
    exclusion), so a real gap between two segments is age this figure has no
    computed rate for, from either side. Concatenating and sorting alone
    would still connect the dots with a straight line across that gap, which
    draws a slope that was never computed and can be far larger than either
    segment's own genuine variation (most visibly at the Holocene/instrumental
    join, where the gap is a few centuries wide against a y-axis scaled to the
    present-day rate). A NaN inserted at the midpoint of each such gap breaks
    the line there instead.
    """
    ages, vals = [seg[0] for seg in segments], [seg[1] for seg in segments]
    pieces_age, pieces_val = [ages[0]], [vals[0]]
    for prev_age, age_i, val_i in zip(ages, ages[1:], vals[1:]):
        gap_hi, gap_lo = prev_age.min(), age_i.max()
        if gap_hi > gap_lo:
            pieces_age.append(np.array([0.5 * (gap_hi + gap_lo)]))
            pieces_val.append(np.array([np.nan]))
        pieces_age.append(age_i)
        pieces_val.append(val_i)
    age = np.concatenate(pieces_age)
    val = np.concatenate(pieces_val)
    order = np.argsort(age)
    age, val = age[order], val[order]
    keep = age >= lo
    return age[keep], val[keep]


co2_rate_age, co2_rate_val = concat_segments(co2_rate_segments)
temp_rate_age, temp_rate_val = concat_segments(temp_rate_segments)

# ==========================================
# 3. THE RATE-ONLY OVERVIEW
# ==========================================
# CO2 on top, GMST below, matching the panel order used throughout this set's
# value and rate figures alike.
LW = 1.8
fig, (ax_co2, ax_t) = plt.subplots(nrows=2, ncols=1, figsize=(11, 7), sharex=True)
fig.subplots_adjust(hspace=0.12)

ax_co2.axhline(0.0, color="black", lw=1.0, alpha=0.5)
ax_co2.plot(co2_rate_age, co2_rate_val, color="tab:blue", lw=LW)
ax_co2.set_ylabel("$\\mathrm{d}CO_2/\\mathrm{d}t$ (ppm per century)", fontsize=11, fontweight="bold")
ax_co2.grid(True, linestyle="--", alpha=0.3)

ax_t.axhline(0.0, color="black", lw=1.0, alpha=0.5)
ax_t.plot(temp_rate_age, temp_rate_val, color="tab:orange", lw=LW)
ax_t.set_ylabel("$\\mathrm{d}T/\\mathrm{d}t$ (°C per century)", fontsize=11, fontweight="bold")
ax_t.grid(True, linestyle="--", alpha=0.3)

axes = (ax_co2, ax_t)
for age in (AGE_ICECORE_START, AGE_HOLOCENE_START, AGE_INSTRUMENTAL_START):
    for ax in axes:
        ax.axvline(age, color="dimgray", lw=0.8, alpha=0.4, zorder=1)

for ax in axes:
    ax.set_xscale("log")
    ax.set_xlim(5e8, AGE_RIGHT_EDGE)
    ax.set_xticks([5e8, AGE_ICECORE_START, AGE_HOLOCENE_START,
                   AGE_INSTRUMENTAL_START, AGE_RIGHT_EDGE])
ax_t.set_xlabel("Years before present (log scale)", fontsize=11, fontweight="bold")
ax_t.set_xticklabels(["500 Myr", "800 kyr", "12 kyr", "200 yr", "10 yr"])
plt.setp(ax_co2.get_xticklabels(), visible=False)

plt.tight_layout()
plt.savefig("images/co2_vs_temperature_rates_overview.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("wrote images/co2_vs_temperature_rates_overview.png")
