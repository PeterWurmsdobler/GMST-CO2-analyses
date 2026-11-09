"""A single log-time overview of the whole 485-million-year record, values and rates.

The four other pairs of figures in this set (plot_gmst.py/plot_co2.py/
plot_holocene.py/plot_industrial.py for the values, plot_derivatives.py for the
rates) each look at one of four nested windows on a LINEAR time axis, because a
linear axis is what lets the reader judge the size of an event by eye within
that window. This script instead puts all four windows on ONE axis, time
before present on a LOGARITHMIC scale running from 485 million years on the
left down to 10 years on the right. A log axis is the only way to fit windows
fifty to six hundred times apart onto one figure at all, and it has a second,
more useful property here: the last 200 years of the record - two ten-
thousandths of the horizontal span in linear time - end up occupying a
substantial slice of the right-hand side, so the recent acceleration
turns into a near-vertical stroke at the edge of the plot at first glance.
That is the point of this figure.

Two things this script has to do that the other five do not.

First, stitch four sources with no common footing into one continuous curve.
CO2 is unproblematic - it is ppm throughout, PhanDA and the ice-core composite
and Mauna Loa all agree on what zero means. Temperature is not: PhanDA's GMST
is an ABSOLUTE global mean in degC, while the ice-core, Holocene and
instrumental temperatures are all ANOMALIES against the 1850-1900 preindustrial
baseline (see plot_co2.py and plot_holocene.py, which already align EDC and
Temp12k to HadCRUT5 the same way). This figure's temperature panels are drawn
in ABSOLUTE terms throughout - the more natural convention for the
Phanerozoic end of the record, where "anomaly vs 1850-1900" has no physical
meaning - so it is the other three records that are converted, by a single
constant shift back onto PhanDA's convention, not the reverse. That constant
is not chosen to match the ice-core value at the 800 kyr handoff, which would
force a coarse 5-Myr interpolation onto a real ice-age cold spell it cannot
resolve, but anchored at PhanDA's own youngest sample instead, close enough to
the 1850-1900 absolute GMST for the Holocene's independently-known stability
to make that a fair proxy for it. The join is left showing a small, honest
step rather than a forced one.

Second, a log axis cannot reach zero, so "time before present" needs an actual
present. REFERENCE_YEAR fixes that; it also fixes what "10 years" at the right
edge means. The four windows hand off to each other at exactly the boundaries
already used throughout this set - 800 kyr, 12 kyr, 200 yr - so this figure is
a splice of curves the reader has already seen, not a fifth, independently
drawn one.

Third - and this is why the figure plots value and rate together rather than
as two separate images - a log time axis can itself mislead. The ancient
record's wild-looking excursions happen over millions of years; the modern
rise looks comparatively smooth because it happens over decades, close enough
together on this axis to render as one continuous stroke rather than a visibly
jagged one. "Looks jagged" and "changes fast" are not the same thing on a
compressed axis, and a reader unused to log scales could easily read the calm
right-hand edge as the gentler one. Stacking each variable's rate directly
beneath its own value, sharing the same x-axis, closes that gap: the wild
ancient wiggle collapses to a rate indistinguishable from zero directly below
it, and the calm-looking modern rise produces the largest excursion on the
whole figure directly below that.
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lowpass import lowpass, lowpass_rate_per_century, resample_uniform, settled  # noqa: E402

REFERENCE_YEAR = 2026  # "present", for the age-before-present log axis
POLAR_AMPLIFICATION = 2.0
PREINDUSTRIAL_FROM, PREINDUSTRIAL_TO = 1850, 1900

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
phan_gmst_abs = df_phan["GMST_50"].to_numpy()

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
edc_temp_raw = df_edc["temp"].to_numpy() / POLAR_AMPLIFICATION

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
mlo_year = df_mlo["decimal_date"].to_numpy()
mlo_age = REFERENCE_YEAR - mlo_year
mlo_co2 = df_mlo["co2"].to_numpy()

# Monthly, not the annual summary series: 'Time' comes as 'YYYY-MM', converted
# to a decimal year at the middle of the month, matching the convention Mauna
# Loa's own decimal_date already uses.
df_had = pd.read_csv("data/hadcrut5-monthly-global.csv").rename(
    columns={"Anomaly (deg C)": "temp"}
)
_y, _m = df_had["Time"].str.split("-", expand=True).astype(int).to_numpy().T
df_had["year"] = _y + (_m - 0.5) / 12.0
had_year = df_had["year"].to_numpy()
had_age = REFERENCE_YEAR - had_year
had_temp_raw = df_had["temp"].to_numpy()

# ==========================================
# 2. HARMONISE THE TEMPERATURE BASELINES
# ==========================================
# HadCRUT5 -> 1850-1900 preindustrial anomaly first, exactly as every other
# figure computes it; this stage is still needed purely to align EDC onto it
# below, even though the final output here is absolute, not anomaly.
preind = (had_year >= PREINDUSTRIAL_FROM) & (had_year <= PREINDUSTRIAL_TO)
had_baseline = had_temp_raw[preind].mean()
had_temp = had_temp_raw - had_baseline

# EDC -> shifted onto the same zero via their mid-20th-century overlap, exactly
# as plot_co2.py aligns EDC and HadCRUT5.
mid_century = (had_year >= 1940) & (had_year <= 1960)
mid_century_level = had_temp[mid_century].mean()
edc_temp = edc_temp_raw + mid_century_level

# Temp12k needs no shift: it is natively anomalised close enough to 1850-1900
# (see plot_holocene.py) that the other figures already treat it as aligned.

# All three records above are now mutually consistent anomalies vs 1850-1900,
# but the figure this script draws uses ABSOLUTE GMST throughout - the more
# natural convention for the Phanerozoic end of the record, where "anomaly
# vs 1850-1900" is not a physically meaningful baseline (see plot_gmst.py's
# caveat on this) - so PhanDA is left untouched and the other three are
# converted the other way instead, by adding back a single constant.
#
# That constant needs an absolute estimate of the 1850-1900 baseline it is
# an anomaly against. Rather than import an outside number, this reuses
# PhanDA's own YOUNGEST sample, 5,850 years ago (mid-Holocene): at 13.99 degC
# absolute it is close to independent estimates of 1850-1900 GMST (~13.7-13.9
# degC), and the Holocene is independently known (see plot_holocene.py) to have stayed within a few tenths of a degree of that
# baseline until the industrial era - so treating this one PhanDA point as
# "preindustrial" is a physically reasonable, self-consistent anchor, exactly
# as it would be if used to shift PhanDA down to anomaly instead. The
# residual mismatch this leaves at the 800 kyr join (PhanDA's coarse
# interpolation reads warmer than EDC's actual measurement there, because
# PhanDA's five-million-year spacing cannot see the glacial cold spell EDC
# measures at that age) is left visible rather than papered over - it
# belongs to the same resolution story the rest of this set's caveats
# already tell, and is exactly as large in absolute terms as it was in
# anomaly terms, since it is a difference and the shift is a constant.
youngest = np.argmin(phan_age)
absolute_anchor = phan_gmst_abs[youngest]
phan_gmst = phan_gmst_abs
edc_temp = edc_temp + absolute_anchor
holo_temp = holo_temp + absolute_anchor
had_temp = had_temp + absolute_anchor
o = np.argsort(edc_age)
edc_boundary = np.interp(AGE_ICECORE_START, edc_age[o], edc_temp[o])
o = np.argsort(phan_age)
phan_boundary = np.interp(AGE_ICECORE_START, phan_age[o], phan_gmst[o])
print(f"Anomaly-to-absolute anchor: PhanDA's youngest sample "
      f"({phan_age[youngest]:,.0f} yr) = {absolute_anchor:.2f} degC absolute, "
      f"taken as the 1850-1900 baseline")
print(f"Residual mismatch at the {AGE_ICECORE_START:,} yr join: "
      f"PhanDA {phan_boundary:.2f} degC vs EDC {edc_boundary:.2f} degC "
      f"({phan_boundary - edc_boundary:+.2f} degC)")

def concat_segments(segments, lo=AGE_RIGHT_EDGE):
    """Join segments end to end for one continuous plotted line - except where
    they are not actually contiguous. Two segments meeting at a boundary this
    figure itself chose (a splice between archives, or - for the rate
    segments computed further below - a filter's own settled() margin) need
    not leave adjacent samples behind: a real gap between them is age
    neither source covers, not merely age this figure declines to draw.
    Concatenating and sorting alone would still connect the dots with a
    straight line across that gap, implying a value (for a value curve) or a
    slope (for a rate curve) that was never measured or computed, and which
    can be far larger than either segment's own genuine variation - most
    visibly for a rate curve at the Holocene/instrumental join, where the gap
    is a few centuries wide against a y-axis scaled to the present-day rate.
    A NaN inserted at the midpoint of each such gap breaks the line there
    instead, for value curves and rate curves alike.
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


# ==========================================
# 3. SPLICE EACH VARIABLE INTO ONE AGE-ORDERED SERIES
# ==========================================
mlo_start_age = mlo_age.min()  # 1958 in age-before-REFERENCE_YEAR terms

co2_age, co2_val = concat_segments([
    (phan_age[phan_age > AGE_ICECORE_START], phan_co2[phan_age > AGE_ICECORE_START]),
    (ice_age[(ice_age <= AGE_ICECORE_START) & (ice_age > mlo_start_age)],
     ice_co2[(ice_age <= AGE_ICECORE_START) & (ice_age > mlo_start_age)]),
    (mlo_age, mlo_co2),
])

temp_age, temp_val = concat_segments([
    (phan_age[phan_age > AGE_ICECORE_START], phan_gmst[phan_age > AGE_ICECORE_START]),
    (edc_age[(edc_age <= AGE_ICECORE_START) & (edc_age > AGE_HOLOCENE_START)],
     edc_temp[(edc_age <= AGE_ICECORE_START) & (edc_age > AGE_HOLOCENE_START)]),
    (holo_age[(holo_age <= AGE_HOLOCENE_START) & (holo_age > AGE_INSTRUMENTAL_START)],
     holo_temp[(holo_age <= AGE_HOLOCENE_START) & (holo_age > AGE_INSTRUMENTAL_START)]),
    (had_age, had_temp),
])

# ==========================================
# 4. COMPUTE THE RATES (plotting is combined with the values below - see the
#    module docstring for why value and rate are shown stacked together)
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
# headline curve in plot_derivatives.py, converted from calendar year to age.
#
# settled() nominally clears the padding by half a cutoff period (10 years) in
# from the start of the record, 1800 - but that early stretch is sparse Law
# Dome ice-core data (annually to triennially resolved), far coarser than the
# monthly Mauna Loa data the half-cutoff margin was really validated against,
# and the plot_derivatives.py figure that already uses this same computation
# never shows anything before 1850 (its xlim starts there) so this rougher
# patch has never actually been on screen before. Clipping this segment to the
# same AGE_INSTRUMENTAL_START boundary used for the value overview - rather
# than trusting settled()'s own boundary here - keeps only the part of the
# curve that figure has already implicitly vetted, and avoids overlapping the
# ice-core-derived segment just above it, which independently reaches down to
# very nearly the same age.
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


co2_rate_age, co2_rate_val = concat_segments(co2_rate_segments)
temp_rate_age, temp_rate_val = concat_segments(temp_rate_segments)

# ==========================================
# 5. THE COMBINED VALUE-AND-RATE OVERVIEW
# ==========================================
# Each variable's value directly above its own rate, both against the same
# logarithmic, reversed time axis - see the note above for why. Laid out as
# two columns (CO2, GMST) of two rows (value, rate) rather than one column of
# four, so the figure reads closer to square: three of the four panels are a
# flat line for nearly their whole width, and stacking all four in one column
# turned most of the page into white space for no informational gain.
LW = 1.8
# Rows are variables (CO2 on top, matching the convention used throughout
# this set of figures), columns are quantity (value on the left, rate on the
# right) - not variable-per-column as in the first version of this figure,
# which put CO2 and GMST side by side and value-over-rate within each column.
fig, ((ax_co2v, ax_co2r), (ax_tv, ax_tr)) = plt.subplots(
    nrows=2, ncols=2, figsize=(11, 8.25), sharex=True)
fig.subplots_adjust(hspace=0.12, wspace=0.25)

ax_co2v.plot(co2_age, co2_val, color="tab:blue", lw=LW)
# 280 ppm: the Holocene ceiling used throughout this set (see
# plot_industrial.py), i.e. the highest level reached in any interglacial of
# the last 800,000 years - the same pre-industrial reference the GMST panel
# marks below, on the variable this figure puts first.
ax_co2v.axhline(280, color="dimgray", ls="--", lw=1.2,
                 label="Pre-industrial $\\mathrm{CO}_2$ ceiling (280 ppm)")
ax_co2v.set_ylabel("$\\mathrm{CO}_2$ (ppm)", fontsize=11, fontweight="bold")
ax_co2v.grid(True, linestyle="--", alpha=0.3)
ax_co2v.legend(loc="upper left", fontsize=8)

ax_co2r.axhline(0.0, color="black", lw=1.0, alpha=0.5)
ax_co2r.plot(co2_rate_age, co2_rate_val, color="tab:blue", lw=LW)
ax_co2r.set_ylabel("$\\mathrm{d}CO_2/\\mathrm{d}t$ (ppm per century)", fontsize=10, fontweight="bold")
ax_co2r.grid(True, linestyle="--", alpha=0.3)

ax_tv.plot(temp_age, temp_val, color="tab:orange", lw=LW)
ax_tv.axhline(absolute_anchor, color="dimgray", ls="--", lw=1.2,
              label=f"1850-1900 preindustrial baseline (≈{absolute_anchor:.1f}°C)")
ax_tv.set_ylabel("GMST (°C)", fontsize=11, fontweight="bold")
ax_tv.grid(True, linestyle="--", alpha=0.3)
ax_tv.legend(loc="upper left", fontsize=8)

ax_tr.axhline(0.0, color="black", lw=1.0, alpha=0.5)
ax_tr.plot(temp_rate_age, temp_rate_val, color="tab:orange", lw=LW)
ax_tr.set_ylabel("$\\mathrm{d}T/\\mathrm{d}t$ (°C per century)", fontsize=10, fontweight="bold")
ax_tr.grid(True, linestyle="--", alpha=0.3)

axes = (ax_co2v, ax_co2r, ax_tv, ax_tr)
for age in (AGE_ICECORE_START, AGE_HOLOCENE_START, AGE_INSTRUMENTAL_START):
    for ax in axes:
        ax.axvline(age, color="dimgray", lw=0.8, alpha=0.4, zorder=1)

# Fewer ticks than the earlier one-column version: the two-column layout
# halves the pixels available per decade, and "100 Myr"/"10 Myr" were only
# ever decorative - the boundaries that matter (800 kyr, 12 kyr, 200 yr) are
# already marked by the vertical guide lines.
for ax in axes:
    ax.set_xscale("log")
    ax.set_xlim(5e8, AGE_RIGHT_EDGE)
    ax.set_xticks([5e8, AGE_ICECORE_START, AGE_HOLOCENE_START,
                   AGE_INSTRUMENTAL_START, AGE_RIGHT_EDGE])
for ax in (ax_tv, ax_tr):
    ax.set_xlabel("Years before present (log scale)", fontsize=11, fontweight="bold")
    ax.set_xticklabels(["500 Myr", "800 kyr", "12 kyr", "200 yr", "10 yr"])
for ax in (ax_co2v, ax_co2r):
    plt.setp(ax.get_xticklabels(), visible=False)

plt.tight_layout()
plt.savefig("images/co2_vs_temperature_overview.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("wrote images/co2_vs_temperature_overview.png")

