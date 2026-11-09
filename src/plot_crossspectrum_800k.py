"""Cross-spectrum of CO2 and Antarctic temperature, 800,000 years, computed as
the Fourier transform of their cross-correlation function directly (the
classical Blackman-Tukey estimator), rather than by Welch's method.

WHY THIS ROUTE RATHER THAN WELCH'S METHOD. An earlier version of this script
estimated coherence and phase by splitting the record into several overlapping
segments, windowing and detrending each, and averaging their periodograms
(scipy.signal.csd/coherence). That works, but it hides two decisions inside
library defaults: a per-segment linear detrend, and a choice of segment length
that trades frequency resolution against the number of segments averaged. This
version instead does the one thing the question actually needs: resample both
series onto a common, uniformly discretised timeline (a common sampling period
of 100 years, the same one plot_crosscorrelation.py's 800 kyr window uses),
subtract only each series' mean (needed to get a correlation rather than a raw
dot product - not itself a filter), form the cross-correlation function
directly, and Fourier-transform THAT. By the Wiener-Khinchin relation this is
the same cross-spectral density Welch's method estimates, just reached without
segmenting the data or filtering it first.

NO LOW-PASS FILTER IS APPLIED TO EITHER SERIES. A moving average is a low-pass
filter with poor characteristics - a rectangular window function in time is a
sinc in frequency, with a slowly decaying, ~13 dB first sidelobe - and using
one to remove a trend before a spectral estimate risks distorting exactly the
low frequencies the estimate is meant to resolve. Nothing of the kind happens
here: the only thing subtracted from either series is its mean. A Hann window
function IS applied, but to the correlation function, not to the data, and
only to control the leakage that truncating that function abruptly at
+-tau_max would otherwise cause - a standard, well-characterised step in this
method (Blackman & Tukey 1958), not a substitute low-pass filter. This is the
signal-processing sense of "window" (a function multiplied onto a sequence to
control spectral leakage), not the "window" the rest of this set uses
for a nested timescale (485 Myr, 800 kyr, ...); the two are unrelated uses of
the same word.

THE CROSS-CORRELATION FUNCTION. With both series on the common sampling
period Delta t, the cross-correlation function is
    R_xy(tau) = (1/N) * sum_n x[n] * y[n + tau/Delta t]
evaluated over the overlap available at each time shift tau, for tau running
from -tau_max to +tau_max where tau_max is half the 800,000-year window under
investigation (about 400,000 years) - beyond that the overlap is too short to
mean much. This is exactly the function plot_crosscorrelation.py's time-domain
figure already plots for this window (there called against "lag"); the only
difference here is that it is left undetrended and pushed through an FFT
rather than searched for a peak.

ONLY THE LOW FREQUENCIES ARE SHOWN. R_xy(tau)'s Fourier transform is complex
at every frequency up to the Nyquist frequency set by the 100-year sampling
period (one cycle per 200 years), but nothing above the orbital band is of
interest here, and displaying it would only add clutter; the figure and the
values reported below are restricted to periods from 15,000 to 400,000 years,
comfortably spanning the precession, obliquity and eccentricity bands with
some margin either side.

RELIABILITY, WITHOUT WELCH'S COHERENCE. A single correlogram over the whole
record gives one estimate of the phase at each frequency, with no averaging
to say how much to trust it - the averaging Welch's method uses for exactly
that purpose is what this route deliberately does without. In its place, the
same calculation is repeated independently on the first and second halves of
the 800,000-year record (each roughly 400,000 years, tau_max about 200,000
years). Where the phase from the full record, the first half and the second
half all agree in sign and rough size, that agreement is the evidence the
result is not an artefact of one particular stretch of ice; where they do
not agree, it is not.

SIGN CONVENTION: matching plot_crosscorrelation.py, a positive tau means CO2
leads GMST/Antarctic temperature.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

POLAR_AMPLIFICATION = 2.0  # see plot_co2.py

# ==========================================
# 1. LOAD AND RESAMPLE ONTO A COMMON, UNIFORMLY DISCRETISED TIMELINE
# ==========================================
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

co2_yr = 1950 - df_ice_co2["Age_BP"].to_numpy()
co2_val = df_ice_co2["co2"].to_numpy()
temp_yr = 1950 - df_edc_temp["Age_BP"].to_numpy()
temp_val = df_edc_temp["temp"].to_numpy()

SAMPLING_PERIOD_YR = 100.0
o1 = np.argsort(co2_yr)
co2_yr, co2_val = co2_yr[o1], co2_val[o1]
o2 = np.argsort(temp_yr)
temp_yr, temp_val = temp_yr[o2], temp_val[o2]
lo, hi = max(co2_yr.min(), temp_yr.min()), min(co2_yr.max(), temp_yr.max())
t = np.arange(lo, hi, SAMPLING_PERIOD_YR)
co2_t = np.interp(t, co2_yr, co2_val)
temp_t = np.interp(t, temp_yr, temp_val)
co2_t = co2_t - co2_t.mean()
temp_t = temp_t - temp_t.mean()
N = len(t)
print(f"N samples: {N}, sampling period {SAMPLING_PERIOD_YR:.0f} yr, "
      f"span {t.max() - t.min():,.0f} yr")


# ==========================================
# 2. CROSS-CORRELATION FUNCTION, THEN ITS FOURIER TRANSFORM
# ==========================================
def cross_correlation_function(x, y, max_shift_samples):
    """R_xy(tau) for tau = k * Delta t, k = -K..K (biased estimator, /N throughout)."""
    n = len(x)
    ks = np.arange(-max_shift_samples, max_shift_samples + 1)
    R = np.empty(len(ks))
    for i, k in enumerate(ks):
        if k >= 0:
            R[i] = np.sum(x[:n - k] * y[k:]) / n
        else:
            R[i] = np.sum(x[-k:] * y[:n + k]) / n
    return ks, R


def cross_spectrum(x, y, dt, max_shift_samples):
    """FFT of the cross-correlation function, passed through a Hann window
    function first -> complex S_xy(f)."""
    _, R = cross_correlation_function(x, y, max_shift_samples)
    R_windowed = R * np.hanning(len(R))
    S = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(R_windowed)))
    freq = np.fft.fftshift(np.fft.fftfreq(len(R), d=dt))
    return freq, S


TAU_MAX_SAMPLES_FULL = N // 2  # tau up to half the 800 kyr window under investigation
freq_full, S_full = cross_spectrum(co2_t, temp_t, SAMPLING_PERIOD_YR, TAU_MAX_SAMPLES_FULL)

half = N // 2
TAU_MAX_SAMPLES_HALF = half // 2
freq_h1, S_h1 = cross_spectrum(co2_t[:half], temp_t[:half], SAMPLING_PERIOD_YR, TAU_MAX_SAMPLES_HALF)
freq_h2, S_h2 = cross_spectrum(co2_t[half:], temp_t[half:], SAMPLING_PERIOD_YR, TAU_MAX_SAMPLES_HALF)

ORBITAL_BANDS = [(100, "eccentricity"), (41, "obliquity"), (23, "precession")]


def band_report(freq, S, label):
    print(f"\n{label}:")
    out = []
    for nominal_kyr, name in ORBITAL_BANDS:
        idx = np.argmin(np.abs(freq - 1 / (nominal_kyr * 1000)))
        phase = np.angle(S[idx])
        tau_yr = -phase / (2 * np.pi * freq[idx])
        out.append((nominal_kyr, name, 1 / freq[idx] / 1000, np.abs(S[idx]), phase, tau_yr))
        print(f"  {name} (~{nominal_kyr} kyr, bin {1/freq[idx]/1000:.1f} kyr): "
              f"|S|={np.abs(S[idx]):.1f}, phase={np.degrees(phase):+.1f} deg, "
              f"tau={tau_yr:+.0f} yr")
    return out


full_bands = band_report(freq_full, S_full, "Full 800 kyr record")
h1_bands = band_report(freq_h1, S_h1, "First half (~400 kyr)")
h2_bands = band_report(freq_h2, S_h2, "Second half (~400 kyr)")

# ==========================================
# 3. FIGURE: MAGNITUDE AND PHASE-DERIVED TIME SHIFT, LOW FREQUENCIES ONLY
# ==========================================
PERIOD_MIN_KYR, PERIOD_MAX_KYR = 15, 180  # beyond ~180 kyr too few cycles fit the
                                           # 800 kyr record for tau to mean anything


def period_kyr_of(freq):
    with np.errstate(divide="ignore"):
        return np.where(freq > 0, 1.0 / freq / 1000.0, np.inf)


fig, (ax_mag, ax_tau) = plt.subplots(nrows=2, ncols=1, figsize=(9, 7), sharex=True)
fig.subplots_adjust(hspace=0.12)

for freq, S, style, label in (
    (freq_full, S_full, dict(color="tab:purple", lw=2.0, zorder=3), "full record"),
    (freq_h1, S_h1, dict(color="tab:purple", lw=1.0, alpha=0.5, ls="--", zorder=2), "first half"),
    (freq_h2, S_h2, dict(color="tab:purple", lw=1.0, alpha=0.5, ls=":", zorder=2), "second half"),
):
    p = period_kyr_of(freq)
    keep = (p >= PERIOD_MIN_KYR) & (p <= PERIOD_MAX_KYR)
    ax_mag.plot(p[keep], np.abs(S[keep]), label=label, **style)
    phase = np.angle(S)
    tau_yr = np.where(freq != 0, -phase / (2 * np.pi * np.where(freq == 0, np.nan, freq)), 0.0)
    ax_tau.plot(p[keep], tau_yr[keep], **style)

ax_mag.set_ylabel("$|S_{xy}(f)|$ (ppm · °C)", fontsize=11, fontweight="bold")
ax_mag.set_yscale("log")
ax_mag.grid(True, linestyle="--", alpha=0.3)
ax_mag.legend(loc="upper right", fontsize=8)

ax_tau.axhline(0.0, color="black", lw=1.0, alpha=0.5)
ax_tau.set_ylabel("$\\mathrm{CO}_2$-leads-GMST $\\tau$ (yr)", fontsize=11, fontweight="bold")
ax_tau.set_ylim(-4800, 3300)
ax_tau.grid(True, linestyle="--", alpha=0.3)

for nominal_kyr, name in ORBITAL_BANDS:
    idx = np.argmin(np.abs(freq_full - 1 / (nominal_kyr * 1000)))
    p = 1 / freq_full[idx] / 1000
    ax_mag.scatter([p], [np.abs(S_full[idx])], color="tab:red", s=35, zorder=5)
    phase = np.angle(S_full[idx])
    tau_yr = -phase / (2 * np.pi * freq_full[idx])
    ax_tau.scatter([p], [tau_yr], color="tab:red", s=35, zorder=5)

for ax in (ax_mag, ax_tau):
    ax.set_xscale("log")
    for nominal_kyr, name in ORBITAL_BANDS:
        ax.axvline(nominal_kyr, color="dimgray", lw=0.8, alpha=0.5, zorder=1)
ax_mag.set_xticks([15, 23, 41, 100, 180])
ax_mag.set_xticklabels(["15", "23", "41", "100", "180"])
for ax in (ax_mag, ax_tau):
    ax.set_xlim(PERIOD_MIN_KYR, PERIOD_MAX_KYR)
for nominal_kyr, name in ORBITAL_BANDS:
    ax_mag.annotate(name, xy=(nominal_kyr, 1.0), xycoords=("data", "axes fraction"),
                     xytext=(-4, -6), textcoords="offset points",
                     rotation=90, ha="right", va="top", fontsize=8, color="dimgray")
ax_tau.set_xlabel("Period (kyr, log scale)", fontsize=11, fontweight="bold")

plt.tight_layout()
plt.savefig("images/co2_vs_temperature_crossspectrum_800k.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("\nwrote images/co2_vs_temperature_crossspectrum_800k.png")
