"""Zero-phase low-pass filtering for the instrumental records.

Shared by plot_industrial.py (which draws the filtered series alongside the raw
one) and plot_derivatives.py (which differentiates it). Both run from the repo
root, so they add this directory to sys.path before importing.

Why a filter rather than a rolling trend. A rolling least-squares slope over N
years is an FIR differentiator with a rectangular support: its frequency
response is a sinc, so it has a first sidelobe only 13 dB down and it never
really stops passing anything. That is tolerable on a smooth ice core and quite
wrong for the Keeling curve, which carries a large seasonal cycle - about 6 ppm
peak to peak - on top of the trend we actually want. Here the passband is set
explicitly instead: a Butterworth low-pass of order FILTER_ORDER, applied with
filtfilt so that it runs forwards and then backwards and introduces no phase
lag. The effective order is therefore twice FILTER_ORDER and the attenuation
twice the figure quoted below, which is the point of doing it this way - a
causal filter would shift the curve to the right by years and misplace every
feature in it.

At a 5-year cutoff the single-pass response at the 1-year seasonal period is
(1 + 5^8)^-1/2, i.e. about -52 dB, and filtfilt doubles that. Nothing of the
annual cycle survives.
"""

import numpy as np
from scipy.signal import butter, filtfilt

# Monthly, matching the native sampling of Mauna Loa. The annual records are
# interpolated up to it; that adds no information, it only gives every series a
# common sample rate so one filter design serves all of them.
SAMPLES_PER_YEAR = 12
FILTER_ORDER = 4


def resample_uniform(years, values, samples_per_year=SAMPLES_PER_YEAR):
    """Put an irregularly sampled record on a uniformly discretised timeline
    in calendar years.

    filtfilt has no notion of a time axis - it assumes evenly spaced samples -
    so this step is not cosmetic. The spliced CO2 record runs from annual to
    triennial Law Dome samples to monthly Mauna Loa ones, and filtering it as
    delivered would apply a different cutoff to each section.
    """
    years = np.asarray(years, dtype=float)
    values = np.asarray(values, dtype=float)
    order = np.argsort(years)
    years, values = years[order], values[order]

    # np.interp needs a strictly increasing x; the splice can leave a repeated
    # year where two records meet.
    keep = np.append(True, np.diff(years) > 0)
    years, values = years[keep], values[keep]

    grid = np.arange(years[0], years[-1], 1.0 / samples_per_year)
    return grid, np.interp(grid, years, values)


def _trend_pad(values, padlen, fit_len):
    """Extend both ends along a least-squares line fitted to the nearest data.

    scipy's own padding options are all wrong for these records, and the right
    hand edge is exactly the part the figures are about, so it is worth being
    explicit. Constant or even padding forces the slope to zero at the boundary
    and flattens the last decades into a plateau that is not there. The default,
    odd padding, reflects through the final sample - x -> 2*x[-1] - x[n-1-k] -
    which continues a clean trend faithfully but doubles whatever happens to sit
    in the last few samples. HadCRUT5 ends on the 2023-2025 El Nino peak with a
    cool 2022 behind it, so odd padding extrapolates the temperature record at
    something like 13 degC per century and the filtered slope inherits it.

    Fitting a line over the last fit_len samples instead means the extrapolation
    carries the recent trend rather than the most recent wiggle.
    """
    n = len(values)
    fit_len = int(min(fit_len, n))
    index = np.arange(n, dtype=float)

    slope, intercept = np.polyfit(index[:fit_len], values[:fit_len], 1)
    head = intercept + slope * np.arange(-padlen, 0, dtype=float)

    slope, intercept = np.polyfit(index[-fit_len:], values[-fit_len:], 1)
    tail = intercept + slope * np.arange(n, n + padlen, dtype=float)

    return np.concatenate([head, values, tail])


def lowpass(values, cutoff_period_yr, samples_per_year=SAMPLES_PER_YEAR):
    """Zero-phase Butterworth low-pass; cutoff given as a period in years.

    A period is the natural unit here: a 5-year cutoff passes everything slower
    than five years and rejects the seasonal cycle and the year-to-year noise.
    """
    values = np.asarray(values, dtype=float)
    nyquist_per_year = samples_per_year / 2.0
    normalised = (1.0 / cutoff_period_yr) / nyquist_per_year
    if not 0.0 < normalised < 1.0:
        raise ValueError(
            f"a {cutoff_period_yr}-year cutoff is beyond the Nyquist period of "
            f"a record sampled {samples_per_year} times a year"
        )
    b, a = butter(FILTER_ORDER, normalised)

    # scipy's default padlen is 3*max(len(a), len(b)) = 27 samples here, barely
    # two years, while a filter with a 20-year cutoff needs some multiple of
    # twenty years to settle; left at the default the output collapses towards
    # zero at both ends. Pad with six cutoff periods, and fit the extrapolation
    # over one cutoff period so the padding is as smooth as the passband.
    window = int(cutoff_period_yr * samples_per_year)
    padlen = min(len(values) - 1, 6 * window)
    padded = _trend_pad(values, padlen, window)
    return filtfilt(b, a, padded, padlen=0)[padlen:padlen + len(values)]


def settled(grid, cutoff_period_yr):
    """Mask selecting the part of a filtered series that is data, not padding.

    A zero-phase filter at each point uses data from roughly half a cutoff
    period either side of it. Within that distance of the end of the record
    there is nothing on one side, so the output depends on how the padding was
    chosen rather than on what was measured - and since the padding above is
    linear, it drives the curvature to zero and makes a still-accelerating rate
    appear to peak and turn back down a few years before the end. That
    turn-down is an artefact and should not be drawn: a derivative is simply
    not available at the boundary without assuming what happens next.
    """
    half = cutoff_period_yr / 2.0
    return (grid >= grid[0] + half) & (grid <= grid[-1] - half)


def lowpass_rate_per_century(years, values, cutoff_period_yr,
                             samples_per_year=SAMPLES_PER_YEAR):
    """Filter first, then differentiate: d(value)/dt in units per century.

    Filtering and differentiating are both linear and time-invariant, so the
    order makes no difference to the result; doing it this way round simply
    avoids forming the raw monthly difference of the Keeling curve, which is
    dominated by the seasonal cycle and some two orders of magnitude larger than
    the trend it hides.
    """
    grid, gridded = resample_uniform(years, values, samples_per_year)
    smoothed = lowpass(gridded, cutoff_period_yr, samples_per_year)
    return grid, np.gradient(smoothed, grid) * 100.0
