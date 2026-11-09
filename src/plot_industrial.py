import os
import sys

import pandas as pd
import matplotlib.pyplot as plt

# Run from the repo root, so the directory holding this script is not on
# sys.path by default.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lowpass import lowpass, resample_uniform  # noqa: E402

# Cutoff periods for the smoothed curves drawn over each raw record. All three
# are the filters plot_derivatives.py differentiates for the companion rate
# figure, in the same colours and relative weights, so that figure is showing
# the slope of curves the reader can see here rather than of some other signal.
#
# 5 years is what it takes to remove the seasonal cycle and leave everything
# else; 20 years is the one whose slope carries the headline numbers,
# for the reason set out at length in plot_derivatives.py; 10 years sits between
# the two, and is drawn mainly so the rate figure has something between "smooth
# but a decade delayed" and "current but noisy" to show.
CUTOFF_PERIODS_YR = (5, 10, 20)
# Lightest/thinnest to darkest/thickest, shortest cutoff to longest - kept in
# sync with the same-named tuples in plot_derivatives.py.
CUTOFF_SHADES_CO2 = ("tab:blue", "royalblue", "navy")
CUTOFF_SHADES_TEMP = ("tab:red", "firebrick", "darkred")
CUTOFF_WIDTHS = (1.3, 1.6, 2.0)

# Fourth and last panel of the series: plot_gmst.py covers 485 Myr, plot_co2.py
# 800 kyr, plot_holocene.py 12 kyr, and this one the 200 years of the industrial
# era - the interval that is narrower than the width of the line in all three of
# the others.
#
# The x-axis is the calendar year rather than years before present: over two
# centuries the 1950 datum of the paleo records is an obstruction rather than a
# convenience.
FIRST_YEAR, LAST_YEAR = 1850, 2030  # HadCRUT5's own native start, not an arbitrary industrial-era date

# The ceiling of the last 800,000 years, never exceeded in any interglacial.
HOLOCENE_CO2_CEILING = 280

PREINDUSTRIAL_FROM, PREINDUSTRIAL_TO = 1850, 1900

# ==========================================
# 1. PARSE MAUNA LOA DATA (Modern Instrumental CO2)
# ==========================================
mlo_path = "data/mauna-loa-co2-mm.txt"
df_mlo = pd.read_csv(
    mlo_path,
    comment='#',
    sep=r'\s+',
    header=None,
    names=['year', 'month', 'decimal_date', 'co2', 'deseasonalized', 'ndays', 'sdev', 'unc']
)

# ==========================================
# 2. PARSE THE ANTARCTIC CO2 COMPOSITE (Bereiter et al. 2015)
# ==========================================
# Over the industrial era the composite is the Law Dome record, which is
# annually to triennially resolved - fine enough to carry this figure on its own
# until Mauna Loa takes over in 1958.
composite_path = "data/antarctica2015co2composite-noaa.txt"
df_ice_co2 = pd.read_csv(composite_path, comment='#', sep='\t').rename(
    columns={'age_gas_calBP': 'Age_BP', 'co2_ppm': 'co2'}
)[['Age_BP', 'co2']].dropna()
df_ice_co2['year'] = 1950 - df_ice_co2['Age_BP']

# ==========================================
# 3. PARSE MODERN INSTRUMENTAL TEMPERATURE (HadCRUT5, monthly)
# ==========================================
# At this scale the uncertainty band is worth showing: it is wide in the 19th
# century, when coverage was sparse, and narrow now.
#
# Monthly, not the annual summary series: 'Time' comes as 'YYYY-MM', converted
# to a decimal year at the middle of the month, matching the convention Mauna
# Loa's own decimal_date already uses.
hadcrut_path = "data/hadcrut5-monthly-global.csv"
df_mod = pd.read_csv(hadcrut_path).rename(
    columns={'Anomaly (deg C)': 'temp_anomaly',
             'Lower confidence limit (2.5%)': 'temp_lo',
             'Upper confidence limit (97.5%)': 'temp_hi'}
)
y, m = df_mod['Time'].str.split('-', expand=True).astype(int).to_numpy().T
df_mod['year'] = y + (m - 0.5) / 12.0
df_mod = df_mod[['year', 'temp_anomaly', 'temp_lo', 'temp_hi']]
# Shift from the native 1961-1990 baseline onto preindustrial, so this figure
# and the three before it agree on what the present-day anomaly is.
baseline = df_mod.loc[
    df_mod['year'].between(PREINDUSTRIAL_FROM, PREINDUSTRIAL_TO), 'temp_anomaly'
].mean()
for column in ('temp_anomaly', 'temp_lo', 'temp_hi'):
    df_mod[column] = df_mod[column] - baseline

# ==========================================
# 4. SPLICE THE CO2 RECORD (one source per interval)
# ==========================================
MLO_START = df_mlo['decimal_date'].min()
df_co2_total = pd.concat([
    df_ice_co2.loc[df_ice_co2['year'] < MLO_START, ['year', 'co2']],
    df_mlo[['decimal_date', 'co2']].rename(columns={'decimal_date': 'year'}),
], ignore_index=True).sort_values(by='year')
df_co2_total = df_co2_total[df_co2_total['year'] >= FIRST_YEAR]

# ==========================================
# 5. GENERATE THE STACKED SUBPLOT LAYOUT
# ==========================================
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
fig.subplots_adjust(hspace=0.15)

# --- TOP PLOT: Carbon Dioxide ---
color_co2 = 'tab:blue'
# The raw monthly record first, in a lighter weight: at this scale the Keeling
# curve's seasonal cycle is visible as a band rather than as an oscillation, and
# it is worth seeing what the filter is being asked to remove.
ax1.plot(df_co2_total["year"], df_co2_total["co2"], color=color_co2, lw=0.8,
         alpha=0.45, label="Law Dome ice core + Mauna Loa (as measured)")
co2_grid, co2_gridded = resample_uniform(df_co2_total["year"], df_co2_total["co2"])
for cutoff, shade, width in zip(CUTOFF_PERIODS_YR, CUTOFF_SHADES_CO2, CUTOFF_WIDTHS):
    ax1.plot(co2_grid, lowpass(co2_gridded, cutoff), color=shade, lw=width,
             zorder=4, label=f"Low-pass filtered, {cutoff}-year cutoff")
ax1.axhline(HOLOCENE_CO2_CEILING, color="dimgray", ls="--", lw=1.0,
            label=f"{HOLOCENE_CO2_CEILING} ppm: the ceiling of the past 800,000 years")
ax1.set_ylabel("Carbon Dioxide (ppm)", fontsize=11, fontweight="bold")
ax1.grid(True, linestyle="--", alpha=0.3)
ax1.legend(loc="upper left", fontsize=8)

current_co2 = df_mlo['co2'].iloc[-1]
current_co2_year = df_mlo['decimal_date'].iloc[-1]
ax1.scatter(current_co2_year, current_co2, color="tab:red", s=30, zorder=5)
ax1.annotate(
    f"{current_co2:.1f} ppm ({int(current_co2_year)})",
    xy=(current_co2_year, current_co2),
    xytext=(1960, current_co2 - 8),
    arrowprops=dict(facecolor='black', arrowstyle='->', lw=0.8),
    fontweight='bold', color='tab:red', fontsize=9, ha='right'
)

# --- BOTTOM PLOT: Global Temperature Anomaly ---
ax2.plot(df_mod["year"], df_mod["temp_anomaly"], color="tab:red", lw=0.8, alpha=0.45,
         label="HadCRUT5 instrumental, global mean (as measured)")
temp_grid, temp_gridded = resample_uniform(df_mod["year"], df_mod["temp_anomaly"])
for cutoff, shade, width in zip(CUTOFF_PERIODS_YR, CUTOFF_SHADES_TEMP, CUTOFF_WIDTHS):
    ax2.plot(temp_grid, lowpass(temp_gridded, cutoff), color=shade, lw=width,
             zorder=4, label=f"Low-pass filtered, {cutoff}-year cutoff")
ax2.fill_between(df_mod["year"], df_mod["temp_lo"], df_mod["temp_hi"],
                 color="tab:red", alpha=0.2, label="2.5th-97.5th percentile")
ax2.axhline(0.0, color="dimgray", ls="--", lw=1.0,
            label="1850-1900 preindustrial baseline")
ax2.set_ylabel("Temperature Anomaly (°C,\nvs. 1850-1900)", fontsize=11, fontweight="bold")
ax2.set_xlabel("Year", fontsize=11, fontweight="bold")
ax2.grid(True, linestyle="--", alpha=0.3)
ax2.legend(loc="upper left", fontsize=8)

current_temp = df_mod['temp_anomaly'].iloc[-1]
current_temp_year = int(df_mod['year'].iloc[-1])
ax2.scatter(current_temp_year, current_temp, color="tab:red", s=30, zorder=5)
ax2.annotate(
    f"+{current_temp:.2f} °C ({current_temp_year})",
    xy=(current_temp_year, current_temp),
    xytext=(1960, current_temp - 0.15),
    arrowprops=dict(facecolor='black', arrowstyle='->', lw=0.8),
    fontweight='bold', color='tab:red', fontsize=9, ha='right'
)

ax2.set_xlim(FIRST_YEAR, LAST_YEAR)
ax2.set_xticks([1850, 1875, 1900, 1925, 1950, 1975, 2000, 2025])

plt.tight_layout()

# Save pristine print assets for your KDP manuscript template
plt.savefig("images/co2_vs_temperature_stacked_200y.png", dpi=300, bbox_inches="tight")
plt.show()
