"""Rates of change behind the four absolute-value figures.

plot_gmst.py, plot_co2.py, plot_holocene.py and plot_industrial.py plot CO2 and
global mean surface temperature over 485 Myr, 800 kyr, 12 kyr and 200 yr. This
script differentiates the same records with respect to time and emits one
derivative figure per timescale, so the argument can move from "how much" to
"how fast".

Three deliberate choices make the figures comparable with each other:

  * Every rate is expressed per century, whatever the span of the figure.
  * Every rate is SIGNED and drawn on a LINEAR axis, so a cooling reads as a
    downward excursion and the size of an event can be judged by eye. A log
    |rate| axis fits the whole span onto one figure, but it flattens precisely
    the events - the Ordovician-Silurian cooling, say - that these figures
    exist to show.
  * On every one of the three reconstructed windows (485 Myr, 800 kyr, 12 kyr)
    the present-day rate is comfortably above the fastest natural rate the
    record resolves - four orders of magnitude at 485 Myr, roughly half an
    order to just under two at 800 kyr and 12 kyr (CO2 and temperature are not
    alike here: see the RATE-FILTER CUTOFFS note below for why) - and even the
    smallest of those ratios is enough to flatten the very cycle each figure
    exists to show, if drawn on the same linear axis. All three therefore
    state the present-day rate in the legend rather than draw it. Only the
    200-year figure, built from the instrumental record itself and needing no
    such comparison, plots it directly.

CAVEAT, and it belongs in the caption: a finite difference cannot resolve
anything faster than the record it is taken from. Each record is differentiated
over a window matched to its own sampling (below), so these are rates sustained
over that window, not instantaneous ones. A 150-year excursion of the modern
kind would be invisible in a record sampled every 5 Myr even if it had happened.

RATE-FILTER CUTOFFS (800 kyr and 12 kyr), and why they are not a free
parameter. The 800 kyr and 12 kyr windows are resampled onto a sampling
period fine enough to avoid aliasing the ice cores' own faster content (see
plot_crosscorrelation.py's module docstring for the same argument made about
the cross-correlation figures), which means the sampling period no longer
does any smoothing on its own; filtered_rate_per_century() below puts that
back explicitly, with a zero-phase Butterworth low-pass ahead of the
derivative. That reintroduces the one choice a coarse sampling period used
to make for free: how much smoothing is enough. Too little, and the
reported "fastest natural rate" is
whatever a single noisy sample pair happens to produce - not sustained, and not
a fair thing to compare a multi-decade instrumental trend against. Too much,
and the comparison is rigged the other way: smooth the past hard enough and
any rate looks unprecedented, which is a criticism this set of figures has to
be able to answer, not just deny.

The cutoff for each series is chosen by asking what fraction of that record's
OWN rate-of-change variance - not its value variance, which is dominated by
slow swings regardless of cutoff and so answers a different question - survives
a given cutoff. (Differentiation multiplies the power spectrum by (2*pi*f)**2,
so unlike the value series, the rate series has real energy at short periods,
and how much of it a given cutoff keeps is an empirical question about the
data, not a matter of taste.) The four cutoffs below were chosen as the
shortest that keep at least 80% of that series' own rate-of-change variance:
700 yr for 800 kyr CO2 (90%), 250 yr for 800 kyr EDC temperature (84%), 220 yr
for 12 kyr CO2 (82%) and 220 yr for 12 kyr Temp12k (87%) - see the FIGURES
entries below for the same numbers attached to each window. The two series
within a window need not, and here do not, share one cutoff: CO2 and
temperature are different physical archives with different intrinsic
timescales, and forcing them onto the same filter would flatter one series at
the other's expense.

These are also, independently, close to or finer than what each archive can
physically resolve regardless of any filtering choice. Gas trapped in firn is
smoothed by diffusion before bubble close-off, and that smoothing is a strong
function of accumulation rate: published estimates run from roughly a decade
at a high-accumulation site like Law Dome (the source for the recent portion
of the CO2 composite used here) to several centuries at a low-accumulation,
full-glacial site like Dome C (Bereiter et al. 2015 and the firn-diffusion
literature it draws on). The 700-year CO2 cutoff sits within that range for
the glacial stretches of the 800 kyr window, not below it. EDC's deuterium
record carries its own, depth-dependent diffusive smoothing of the
water-isotope signal, documented at the ice core's oldest, most compacted end
as comparable to the 250-year cutoff used here (Jouzel et al. 2007 and the
EDC3 chronology papers). Temp12k's own authors screened contributing records
for a temporal resolution of at least a few hundred years and built their own
headline analysis around a 200-year smoothing scale (Kaufman et al. 2020),
close to the 220 years used for both series in the 12 kyr window here. None
of the four cutoffs, in other words, asks an archive to hand over more
temporal precision than it is independently documented to have.

One consequence of a single fixed cutoff applied to a record whose own
resolution is not fixed, worth stating plainly since it shows up in the 800
kyr figure: the ice cores are sampled far more densely near the present than
deep in the past. EDC temperature runs a median 47 years between native
samples in the most recent 100 kyr against 621 years in the 700-800 kyr
stretch; the CO2 composite similarly, 44 years against 532 years - both
coarser, at depth, than the 250-year cutoff itself. A 250-year low-pass
genuinely resolves whatever the temperature record varies by on that
timescale where the record is dense enough to support it, and can only
interpolate smoothly across the cutoff where it is not. So the filtered
curve's own oscillation amplitude growing towards the present, and each
glacial termination's own warming looking less abrupt the further back it
lies (peak +dT/dt falls from 0.51 degC/century at Termination I, 11 kyr BP,
to 0.07 degC/century at Termination VIII, 715 kyr BP), are the same effect
seen twice, not two separate ones: declining native sample density with
age, not declining rates of change. It is a property of the archive the
cutoff is applied to, not an artefact of the cutoff itself - tuning the
cutoff per epoch to compensate would trade an honest limitation for a
hidden one. (Note: within a termination's broad age span the record can
also carry a real, separately documented cooling reversal - e.g. the
Antarctic Cold Reversal partway through Termination I - so the largest-
magnitude excursion in that span is not always the termination's own
warming; the numbers above are each termination's peak positive rate,
found deliberately rather than by taking the largest |dT/dt| in the
window.)

Passing the resolution test above is not the same as being a genuine
global-mean swing, and this is worth checking rather than assuming. EDC is
one ice core at one site; the constant global/Antarctic scaling used
throughout this project is calibrated against the slow, orbitally-paced
glacial cycle, not against centennial wiggles riding on it. Published work
on this record reports centennial-scale variability as consistently larger
during glacial periods than interglacials - a direct, checkable prediction.
Testing it against the filtered curve here finds no such contrast: the
Holocene (0-8 kyr BP), the most climatically settled stretch this record
has, gives a rate-of-change standard deviation of 0.20 degC/century,
higher than a distinctly glacial, AIM-event-rich stretch (60-112 kyr BP,
0.16), a full glacial stretch (150-190 kyr BP, 0.14), and comparable to the
Eemian interglacial optimum (120-128 kyr BP, 0.17) - all four eras checked
after excluding filter-settling edges. The Holocene is also the shallowest,
least compacted, most finely sampled part of the core, so resolution (not
climate state) is the better explanation for its variability - which
means some of what a 250-year cutoff lets through even in the best-
resolved stretch is more plausibly single-site noise (local precipitation
intermittency, near-surface diffusion, measurement scatter) than a
sustained global temperature swing. fastest_natural() therefore reports a
defensible UPPER BOUND on natural centennial warming, not a precise
estimate: the present-day rate is *at least* the ratio quoted in the
figures and table, quite possibly more.

fastest_natural() is also a single number, and a single noisy or
locally-sourced sample can dominate a single number. natural_rate_dispersion()
below answers a more robust question instead: not "what is the most extreme
thing this record has ever done" but "how far outside the record's ordinary
spread does today's rate sit", from a mean and standard deviation taken over
several hundred to a thousand samples rather than the one largest of them -
restricted, as above, to the most recent, best-resolved 100 kyr, so noise
that is merely well-resolved is not mistaken for a wide natural spread. Over
that stretch the temperature rate has effectively zero mean (as any rate of
a bounded quantity must, averaged over enough samples) and a standard
deviation of 0.18 degC/century; today's 2.1 degC/century sits roughly twelve
standard deviations from that mean, not one unusually large excursion among
many. This is a descriptive comparison, not a formal significance test - the
underlying rate series is strongly autocorrelated, not independent samples -
but a twelve-sigma-equivalent departure does not depend on trusting any
single peak.
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# These scripts are run from the repo root, so the directory holding them is not
# on sys.path by default.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lowpass import lowpass, lowpass_rate_per_century, resample_uniform, settled  # noqa: E402

# Differencing window per record, matched to its native sampling: PhanDA is
# stage-resolved (median 5.3 Myr), the ice-core composite runs a few hundred to
# a few thousand years, and the instrumental records are annual.
WINDOW_PHANEROZOIC_YR = 5_000_000
WINDOW_ICECORE_YR = 100
WINDOW_HOLOCENE_YR = 100
MODERN_TREND_YEARS = 20  # span used for the present-day reference rate

POLAR_AMPLIFICATION = 2.0  # see plot_co2.py

# The five mass extinctions of the Phanerozoic, in Ma; same list as plot_gmst.py.
MASS_EXTINCTIONS = [
    (445, "Ordovician–Silurian"),
    (372, "Late Devonian"),
    (252, "Permian–Triassic"),
    (201, "Triassic–Jurassic"),
    (66, "Cretaceous–Paleogene"),
]


def rate_per_century(age_yr, values, window_yr):
    """Signed d(value)/dt in units per century, a numerical derivative over a sampling period of window_yr.

    Ages are years before present and so run backwards in time; the sign is
    flipped to give a forward-time rate, positive for rising. The series is
    first placed on a uniformly discretised timeline of sampling period
    window_yr, which both fixes the differencing interval and stops the
    result from being an artefact of uneven proxy sampling.
    """
    order = np.argsort(age_yr)
    age_yr, values = np.asarray(age_yr)[order], np.asarray(values)[order]

    grid = np.arange(age_yr.min(), age_yr.max(), window_yr)
    if len(grid) < 3:
        raise ValueError("window_yr too coarse for the span of this record")
    gridded = np.interp(grid, age_yr, values)

    # Forward in time = decreasing age, hence the negation.
    rate_per_year = -np.gradient(gridded, grid)
    return grid, rate_per_year * 100.0


def filtered_rate_per_century(age_yr, values, step_yr, cutoff_yr):
    """Zero-phase low-pass, then numerical differentiation, in units per century.

    rate_per_century() above resamples onto a sampling period of step_yr and
    differentiates directly; the sampling period itself is then the only
    thing standing between a genuine sustained rate and a two-point
    numerical derivative amplifying whatever noise sits between adjacent
    samples. That is fine while the sampling period is coarse enough to
    average over several samples per differencing interval, but resampling
    finer - to keep faithful to what these records actually resolve rather
    than aliasing it away - removes that implicit smoothing along with the
    aliasing risk. This restores it explicitly: a zero-phase Butterworth
    low-pass at cutoff_yr, run before differentiating rather than left to the
    sampling period to provide by accident, so "sustained" keeps meaning what
    it says regardless of how fine step_yr is.

    The filter itself, though, needs several samples per cutoff period to be
    well-conditioned, and step_yr is not always fine enough to give it that:
    a 220-year cutoff on a 100-year grid is a normalised frequency of 0.91,
    a 250-year cutoff 0.80, both uncomfortably close to that grid's own
    Nyquist limit of 1.0. Confirmed empirically, not just in principle: at
    step_yr = 100 yr this produced a fastest-rate reading for one window that
    changed by a factor of two, and flipped sign, once the grid was refined.
    The fix is standard practice for a filter run near its own Nyquist limit
    - design and apply it on a finer working grid, twenty samples per cutoff
    period, then resample the result back onto step_yr, which stays exactly
    what is described in the text and figures as this window's sampling
    period; only the filter's own internal arithmetic runs finer than that.
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


def rolling_trend_per_century(years, values, span=MODERN_TREND_YEARS):
    """Least-squares trend of the preceding `span` years, year by year.

    Used for the instrumental segment. Defining it as a rolling version of the
    single number quoted as the present-day rate means the red segment ends
    exactly on the dashed reference line rather than near it.

    Returns (age_BP, rate per century), ages descending as elsewhere.
    """
    years = np.asarray(years, dtype=float)
    values = np.asarray(values, dtype=float)
    order = np.argsort(years)
    years, values = years[order], values[order]

    targets = np.arange(np.ceil(years.min()) + span, np.floor(years.max()) + 1)
    slopes = []
    for year in targets:
        window = (years > year - span) & (years <= year)
        slopes.append(np.polyfit(years[window], values[window], 1)[0] * 100.0)
    return 1950.0 - targets, np.array(slopes)


def mark_extinctions(ax, label_axis=False):
    """Draw the five mass extinctions as vertical lines, labelled on one axis."""
    for age, name in MASS_EXTINCTIONS:
        ax.axvline(age, color="dimgray", ls="--", lw=1.0, alpha=0.7, zorder=1)
        if label_axis:
            ax.annotate(
                name, xy=(age, 1.0), xycoords=("data", "axes fraction"),
                xytext=(-4, -6), textcoords="offset points",
                rotation=90, ha="right", va="top",
                fontsize=8, color="dimgray",
            )


# ==========================================
# 1. LOAD THE THREE RECORDS
# ==========================================
# --- Phanerozoic (PhanDA, 485 Myr, ages in Ma) ---
df_phan = pd.read_csv("data/PhanDA_GMSTandCO2_percentiles.csv")
phan_age_yr = df_phan["AverageAge"].to_numpy() * 1e6

# --- Ice core CO2 composite (Bereiter et al. 2015) ---
df_ice_co2 = pd.read_csv(
    "data/antarctica2015co2composite-noaa.txt", comment="#", sep="\t"
).rename(columns={"age_gas_calBP": "Age_BP", "co2_ppm": "co2"})[["Age_BP", "co2"]].dropna()

# --- EPICA Dome C temperature (EDC3), scaled to global-mean equivalent ---
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

# --- Holocene global temperature (Temp12k, Kaufman et al. 2020) ---
df_holo = pd.read_csv(
    "data/temp12k_allmethods_percentiles.csv", skipinitialspace=True
).rename(columns={"ages": "Age_BP", "global_median": "temp"})[["Age_BP", "temp"]]

# --- Modern instrumental ---
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

# Mauna Loa only starts in 1958, too short for a 50-year rolling trend to reach
# back far enough to be drawn as a segment. Carry it back with the ice-core
# composite, whose Holocene section is the high-resolution Law Dome record and
# which overlaps Mauna Loa to 2001 - the same splice plot_co2.py makes.
co2_instrumental = pd.concat([
    df_ice_co2.assign(year=1950 - df_ice_co2["Age_BP"]).query("1800 <= year < 1958")[["year", "co2"]],
    df_mlo[["decimal_date", "co2"]].rename(columns={"decimal_date": "year"}),
])

MODERN_CO2_AGE, MODERN_CO2_SERIES = rolling_trend_per_century(
    co2_instrumental["year"], co2_instrumental["co2"])
MODERN_TEMP_AGE, MODERN_TEMP_SERIES = rolling_trend_per_century(
    df_had["year"], df_had["temp"])

# The single number quoted in the text and the table: the last point of each
# rolling series, i.e. the trend over the most recent MODERN_TREND_YEARS years.
MODERN_CO2_RATE = abs(MODERN_CO2_SERIES[-1])
MODERN_TEMP_RATE = abs(MODERN_TEMP_SERIES[-1])
LATEST_YEAR = int(1950 - min(MODERN_CO2_AGE.min(), MODERN_TEMP_AGE.min()))

# ==========================================
# 2. THE THREE FIGURES
# ==========================================
XLABEL_BP = f"Years Before Present (Present = 1950; instrumental record to {LATEST_YEAR})"

FIGURES = [
    dict(
        name="rates_485Ma",
        window=WINDOW_PHANEROZOIC_YR,
        xdiv=1e6,
        xlabel="Age (Millions of Years Ago)",
        window_label={"co2": "5 Myr", "temp": "5 Myr"},
        co2=(phan_age_yr, df_phan["CO2_50"].to_numpy()),
        temp=(phan_age_yr, df_phan["GMST_50"].to_numpy()),
        xlim=(485, 0),
        # 191 ppm and 2.1 degC per century cannot share a linear axis with a
        # record whose fastest resolved change is four orders of magnitude
        # smaller, so state the present day here rather than drawing it.
        show_modern=False,
        extinctions=True,
    ),
    dict(
        name="rates_800k",
        window=WINDOW_ICECORE_YR,
        # Each series filtered to its own evidence-based cutoff (see the
        # RATE-FILTER CUTOFFS note above), not a shared one: 700 yr keeps
        # 90% of the CO2 record's own rate-of-change variance, 250 yr keeps
        # 84% of the EDC temperature record's - the two series simply do not
        # carry their genuine variability at the same timescale.
        cutoff={"co2": 700, "temp": 250},
        xdiv=1.0,
        xlabel=XLABEL_BP,
        window_label={"co2": "700 yr", "temp": "250 yr"},
        co2=(df_ice_co2["Age_BP"].to_numpy(), df_ice_co2["co2"].to_numpy()),
        temp=(df_edc_temp["Age_BP"].to_numpy(), df_edc_temp["temp"].to_numpy()),
        xlim=(806000, -100),
        # The present-day rate is well over an order of magnitude above the
        # fastest glacial-interglacial rate this window resolves - still enough
        # to flatten the natural cycle to a line along zero if drawn on the
        # same linear axis. State it, don't draw it.
        show_modern=False,
        extinctions=False,
        # Sized to the low-pass filtered curve's own envelope (see
        # compute_rate()), not the raw sampling period's - the cycles this
        # figure exists to show all fit comfortably inside it.
        ylim={"co2": (-5, 6), "temp": (-0.85, 0.85)},
        legend_loc={"co2": "upper left", "temp": "upper left"},
    ),
    dict(
        name="rates_12k",
        window=WINDOW_HOLOCENE_YR,
        # 220 yr keeps 82% (CO2) and 87% (temperature) of each record's own
        # rate-of-change variance - close enough between the two that one
        # shared cutoff serves both without favouring either.
        cutoff={"co2": 220, "temp": 220},
        xdiv=1.0,
        xlabel=XLABEL_BP,
        window_label={"co2": "220 yr", "temp": "220 yr"},
        co2=(df_ice_co2["Age_BP"].to_numpy(), df_ice_co2["co2"].to_numpy()),
        temp=(df_holo["Age_BP"].to_numpy(), df_holo["temp"].to_numpy()),
        xlim=(12000, -100),
        legend_loc={"co2": "upper left"},
        # Same reasoning as 800 kyr above: still enough of a gap to flatten the
        # Holocene's own natural variability if shared with the present-day
        # rate on one linear axis.
        show_modern=False,
        extinctions=False,
    ),
]

PREINDUSTRIAL_AGE_BP = 100  # 1850 CE, HadCRUT5's own native start

# The eight glacial terminations, at the age of peak warming rate each reaches
# in the filtered temperature curve below - see the cooling-reversal note
# above for why this is the peak positive rate in each termination's age
# span, not the largest-magnitude one. Reused as-is in plot_co2.py so the
# markers on the value figure agree with the ones here.
TERMINATIONS_YR_BP = (11438, 134938, 244138, 335538, 425438, 522438, 621438, 715038)


def compute_rate(spec, key):
    """(grid, rate, natural) for one series, low-pass filtered first when spec
    carries a cutoff (the two finer, resampled-for-resolution windows), a plain
    numerical derivative otherwise (485 Myr, whose sampling period did not
    need to change).

    "natural" drops two kinds of edge: points whose differencing or
    filtering window reaches into the industrial era, not merely those centred
    in it, so "natural" means what it says; and, for the filtered windows,
    whichever few points at either end of the record the filter's own padding
    still influences (settled()), so a boundary artefact is never mistaken for
    the fastest rate the record resolves.
    """
    age, val = spec[key]
    if "cutoff" in spec:
        cutoff = spec["cutoff"][key]
        grid, rate = filtered_rate_per_century(age, val, spec["window"], cutoff)
        natural = settled(grid, cutoff) & (grid >= PREINDUSTRIAL_AGE_BP)
    else:
        grid, rate = rate_per_century(age, val, spec["window"])
        natural = grid >= PREINDUSTRIAL_AGE_BP + spec["window"]
    natural = natural & (grid / spec["xdiv"] <= spec["xlim"][0])
    return grid, rate, natural


def fastest_natural(spec, key):
    """Largest rate magnitude a record resolves, with the industrial era removed."""
    grid, rate, natural = compute_rate(spec, key)
    return np.abs(rate[natural]).max()


RECENT_WELL_RESOLVED_YR = 100_000  # the same window the density comparison in
                                    # the sample-density comparison in the docstring uses


def natural_rate_dispersion(spec, key, window_yr=RECENT_WELL_RESOLVED_YR):
    """Mean and standard deviation of the natural rate distribution, restricted
    to the most recent window_yr years - the same "fair", best-resolved
    stretch the sample-density comparison in the module docstring already uses,
    so the two claims are made on the same terms rather than one over the
    whole 800 kyr and the other over a fraction of it.

    fastest_natural() above is one number: the single largest excursion, which
    a single noisy or locally-sourced sample can dominate (see the caveat in
    this module's docstring). A mean and standard deviation taken over
    several hundred to a thousand samples cannot be swung by any one of them,
    so it answers a different, more robust question: not "what is the most
    extreme thing this record has ever done" but "how far outside the
    record's ordinary spread does today's rate sit".
    """
    grid, rate, natural = compute_rate(spec, key)
    mask = natural & (grid <= window_yr)
    return rate[mask].mean(), rate[mask].std()


for spec in FIGURES:
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
    fig.subplots_adjust(hspace=0.15)

    for ax, key, (age, val), colour, ylabel, modern, modern_series, modern_fmt in (
        (ax1, "co2", spec["co2"], "tab:blue", "$\\mathrm{d}CO_2/\\mathrm{d}t$ (ppm per century)",
         MODERN_CO2_RATE, (MODERN_CO2_AGE, MODERN_CO2_SERIES),
         f"{MODERN_CO2_RATE:.0f} ppm/century"),
        (ax2, "temp", spec["temp"], "tab:orange", "$\\mathrm{d}T/\\mathrm{d}t$ (°C per century)",
         MODERN_TEMP_RATE, (MODERN_TEMP_AGE, MODERN_TEMP_SERIES),
         f"{MODERN_TEMP_RATE:.1f} °C/century"),
    ):
        # The same "natural" mask compute_rate() already uses for the reported
        # fastest rate also bounds what gets drawn: the ice-core records run to
        # 806 kyr and would otherwise set the y-scale of the 12 kyr figure, the
        # last stretch before the industrial era is the steepest natural-looking
        # rate a filtered, faithfully-resolved curve can show without actually
        # being industrial (this window's whole reason for excluding it from
        # "fastest natural" rather than just not drawing it), and, for the
        # filtered windows, the filter's own settling margin would otherwise
        # show up as an unexplained excursion right at the same edge.
        grid, rate, visible = compute_rate(spec, key)
        grid, rate = grid[visible], rate[visible]

        ax.axhline(0.0, color="black", lw=0.8, alpha=0.5)
        ax.plot(grid / spec["xdiv"], rate, color=colour, lw=1.8,
                label=f"Reconstruction (Δt = {spec['window_label'][key]})")

        if spec["show_modern"]:
            ax.plot(*modern_series, color="tab:red", lw=1.8, zorder=4,
                    label=f"Instrumental, {MODERN_TREND_YEARS}-year rolling trend")
            ax.axhline(modern, color="tab:red", ls="--", lw=1.0,
                       label=f"Present-day instrumental rate: {modern_fmt}")
        else:
            # The same pre-industrial "fastest natural rate" fastest_natural()
            # reports below -
            # not np.abs(rate).max() over the visible (but not pre-industrial-
            # only) series above, which for the two finer windows would pick up
            # the early industrial uptick already present in the ice-core tail
            # and understate the ratio.
            ratio = modern / fastest_natural(spec, key)
            ax.plot([], [], " ",
                    label=f"Present-day instrumental rate: {modern_fmt} "
                          f"($\\approx${ratio:,.0f}$\\times$ the fastest here, off scale)")

        ax.set_ylabel(ylabel, fontsize=10, fontweight="bold")
        ax.grid(True, linestyle="--", alpha=0.3)
        if spec["name"] == "rates_800k":
            for age in TERMINATIONS_YR_BP:
                ax.axvline(age, color="firebrick", ls="--", lw=1.0, alpha=0.8, zorder=1)
            if ax is ax1:
                ax.plot([], [], color="firebrick", ls="--", lw=1.0, alpha=0.8,
                        label="Glacial terminations")
        if spec["extinctions"]:
            # Headroom at the top for the event labels, and at the bottom for
            # the legend, so neither sits on the curve.
            lo, hi = ax.get_ylim()
            ax.set_ylim(lo - 0.35 * (hi - lo), hi + 0.55 * (hi - lo))
            mark_extinctions(ax, label_axis=(ax is ax1))
        if "ylim" in spec and key in spec["ylim"]:
            ax.set_ylim(*spec["ylim"][key])
        ax.legend(loc=spec.get("legend_loc", {}).get(key, "lower left"), fontsize=8)

    ax2.set_xlabel(spec["xlabel"], fontsize=11, fontweight="bold")
    ax2.set_xlim(*spec["xlim"])
    if spec["xdiv"] == 1.0:
        ax2.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ",")))

    plt.tight_layout()
    plt.savefig(f"images/co2_vs_temperature_{spec['name']}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote images/co2_vs_temperature_{spec['name']}.png  "
          f"(CO2 Δt={spec['window_label']['co2']}, temp Δt={spec['window_label']['temp']})")

# ==========================================
# 3. THE FOURTH WINDOW: 200 YEARS
# ==========================================
# The industrial era is narrower than the width of the line in all three figures
# above, so it gets one of its own. There is no reconstruction to difference
# here, only the instrumental records - and they are the only records in this
# set fast enough to carry anything worth filtering properly.
#
# So the curves here are not rolling trends. Each record is resampled onto a
# uniformly discretised timeline with a monthly sampling period, passed
# through a zero-phase Butterworth low-pass at three cutoff periods, and then
# differentiated (see lowpass.py).
#
# The signal differentiated here is the 20-year curve drawn in
# plot_industrial.py, so this figure shows the slope of a curve the reader can
# actually look at rather than of some other one.
#
# On why it is the 20-year curve and not the 5-year one, which is the better
# looking of the two in the value figure. Differentiation is itself a high-pass,
# +20 dB per decade, so the passband that cleans up a VALUE is nowhere near
# enough to clean up a SLOPE. The 5-year curve looks smooth plotted as a value
# only because a 1 ppm wobble is invisible on an axis 150 ppm tall; as a slope
# the same wobble is enormous. The interannual variability left in it is real -
# ENSO, volcanoes - and it is not what this figure is about:
#
#     Mauna Loa annual growth    1.67 +- 0.71 ppm/yr   -> +-71 ppm/century
#     HadCRUT5 year-to-year      +- 0.117 degC/yr      -> +-12 degC/century
#
# against trends of 240 ppm/century and 2 degC/century. Differentiating the
# 5-year signal duly gives a curve swinging over +-150 ppm/century and
# +-10 degC/century, which tells the reader nothing. Scatter of the filtered
# slope, measured over a quiet stretch of each record:
#
#     cutoff      d(CO2)/dt         dT/dt
#       5 yr   81.3 ppm/cy    5.27 degC/cy
#      10 yr   17.8 ppm/cy    1.80 degC/cy
#      20 yr    5.2 ppm/cy    0.79 degC/cy
#
# Hence 20 years: the shortest cutoff at which the slope is legible, and still
# short enough to resolve the mid-century plateau in the temperature panel.
#
# filtfilt needs data on both sides of every point (see _trend_pad() and
# settled() in lowpass.py), and towards either end of the record it is working
# increasingly on padding rather than data. The padding used there is linear,
# which drives the curvature to zero - so a rate that is in fact still climbing
# eventually appears, close enough to the end, to peak and turn over. Rather
# than pick a single cutoff and either dash its unsupported tail (which still
# ends on the artefact) or paper over the gap with trailing-window endpoints
# (which draws coloured labels for windows that have no curve of their own),
# this figure now plots the same three cutoffs as the value figure, in the same
# colours and relative weights, plus the unfiltered year-on-year change.
#
# settled() trims each cutoff's own padding-affected ends, and a shorter cutoff
# needs less padding: the 20-year curve is settled only to 2016-ish, the
# 10-year to about 2020-2022, the 5-year to about 2022-2024. So the thick,
# smooth 20-year curve carries the headline shape and stops earliest; the
# thinner 10- and 5-year curves are noisier - real interannual variability,
# ENSO and volcanoes, not filter artefacts - but each runs several years closer
# to the present, so between them the last decade is not simply left blank.
# The thinnest line of all is the unfiltered year-on-year change: not the raw
# monthly record (its seasonal cycle would swing the CO2 panel over roughly
# +-2000 ppm/century and swamp everything else on it), but the same records
# reduced to one value a year first. It needs no filter and so has nothing to
# hide: it runs right up to the last measurement, oscillating a good deal, and
# is drawn thin exactly so that oscillation doesn't read as the message of the
# figure.
CUTOFF_PERIODS_YR = (5, 10, 20)
# Kept in sync with the same-named tuples in plot_industrial.py.
CUTOFF_SHADES_CO2 = ("tab:blue", "royalblue", "navy")
CUTOFF_SHADES_TEMP = ("tab:red", "firebrick", "darkred")
CUTOFF_WIDTHS = (1.3, 1.6, 2.0)


def annual_rate_per_century(years, values):
    """Year-on-year change, per century, with no Butterworth filter behind it.

    A DESIGN CHOICE, NOT A DEFAULT. A proper decimation filter - a low-pass
    at the new Nyquist period (2 years, for a 12:1 monthly-to-annual
    decimation) applied before resampling - is the textbook-correct way to
    remove Mauna Loa's seasonal cycle ahead of an annual difference, and
    filtering and differentiating commute (both are linear and
    time-invariant), so it should not matter which comes first. It was
    tried, and rejected here for a specific reason: the CO2 series this
    function is applied to splices sparse, annual-to-triennial Law Dome ice
    core data (pre-1958) onto dense monthly Mauna Loa data, and a 2-year
    filter applied across it produces spurious rate swings of several
    hundred ppm per century scattered through the pre-1958 stretch, not
    confined to the handoff itself - the ice core simply does not resolve
    finely enough to support that cutoff anywhere in it, the same resolution
    limit this set meets everywhere else, here at a splice rather
    than inside one record.

    What is used instead - resample onto a uniformly discretised timeline
    with a strictly annual sampling period, sample_per_year=1, then
    difference - has no such failure mode, and is
    not simply "no filter": every annual sample lands at the same point in
    the seasonal cycle, which removes that cycle (about 6 ppm peak to peak,
    confirmed against NOAA's own deseasonalized series and by FFT of the
    detrended monthly record) much as a notch tuned exactly to its frequency
    would, without being designed as one. Differencing then removes the
    slow secular trend the same way it does for every curve in this figure.
    Together, the two acts like a band-pass: what survives in between is the
    band this curve exists to show, real year-to-year variability (ENSO,
    volcanic eruptions), not an artefact of sampling rate.
    """
    grid, gridded = resample_uniform(years, values, samples_per_year=1)
    return grid, np.gradient(gridded, grid) * 100.0


fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
fig.subplots_adjust(hspace=0.15)

for ax, (years, values), shades, ylabel in (
    (ax1, (co2_instrumental["year"], co2_instrumental["co2"]), CUTOFF_SHADES_CO2,
     "$\\mathrm{d}CO_2/\\mathrm{d}t$ (ppm per century)"),
    (ax2, (df_had["year"], df_had["temp"]), CUTOFF_SHADES_TEMP,
     "$\\mathrm{d}T/\\mathrm{d}t$ (°C per century)"),
):
    ax.axhline(0.0, color="black", lw=0.8, alpha=0.5)

    annual_grid, annual_rate = annual_rate_per_century(years, values)
    ax.plot(annual_grid, annual_rate, color=shades[0], lw=0.7, alpha=0.45,
            label="Year-on-year change (as measured)")

    for cutoff, shade, width in zip(CUTOFF_PERIODS_YR, shades, CUTOFF_WIDTHS):
        grid, rate = lowpass_rate_per_century(years, values, cutoff)
        keep = settled(grid, cutoff)
        ax.plot(grid[keep], rate[keep], color=shade, lw=width, zorder=4,
                label=f"Derivative of the {cutoff}-year filtered record")

    ax.set_ylabel(ylabel, fontsize=10, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)

ax2.set_xlabel("Year", fontsize=11, fontweight="bold")
ax2.set_xlim(1850, 2030)
ax2.set_xticks([1850, 1875, 1900, 1925, 1950, 1975, 2000, 2025])

plt.tight_layout()
plt.savefig("images/co2_vs_temperature_rates_200y.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("wrote images/co2_vs_temperature_rates_200y.png  "
      f"(annual change plus {CUTOFF_PERIODS_YR}-year filtered derivatives)")

print(f"\nPresent-day rates over the last {MODERN_TREND_YEARS} years (to {LATEST_YEAR}):")
print(f"  CO2  {MODERN_CO2_RATE:8.1f} ppm/century")
print(f"  GMST {MODERN_TEMP_RATE:8.2f} °C/century")

# Summary numbers.
print("\nFastest natural rate resolved in each window (pre-industrial only):")
for spec in FIGURES:
    for label, key in (("CO2 ", "co2"), ("GMST", "temp")):
        peak = fastest_natural(spec, key)
        print(f"  {spec['name']:11s} {label}  |max| {peak:10.4f}"
              f"   present day is {(MODERN_CO2_RATE if key == 'co2' else MODERN_TEMP_RATE) / peak:8.0f}x")

print(f"\nNatural rate dispersion, 800 kyr window, most recent {RECENT_WELL_RESOLVED_YR:,} yr only "
      "(pre-industrial, settled):")
rates_800k_spec = next(s for s in FIGURES if s["name"] == "rates_800k")
for label, key, modern in (("CO2 ", "co2", MODERN_CO2_RATE), ("GMST", "temp", MODERN_TEMP_RATE)):
    mean, std = natural_rate_dispersion(rates_800k_spec, key)
    z = (modern - mean) / std
    print(f"  {label}  mean {mean:+8.4f}  std {std:8.4f}   present day is {z:+7.1f} std-devs from the mean")
