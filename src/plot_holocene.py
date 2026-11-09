import pandas as pd
import matplotlib.pyplot as plt

# Third panel of the series: plot_gmst.py covers 485 Myr, plot_co2.py 800 kyr,
# this one zooms into the 12,000 years of Holocene stability that agriculture
# and civilisation grew up in - and the industrial departure from it. 12,000
# yr BP is Temp12k's own native upper limit, not an arbitrary round number,
# and sits close to the ICS-official Holocene onset of 11,700 yr BP.
HOLOCENE_START_BP = 12000

# Unlike plot_co2.py, both temperature series here are already GLOBAL means, so
# no polar amplification correction is needed. Both are referenced to the
# preindustrial baseline: Temp12k is anomalised to 1800-1900 CE, and HadCRUT5
# is re-referenced below to 1850-1900 CE, the IPCC convention.
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
df_mlo['Age_BP'] = 1950 - df_mlo['decimal_date']
df_mlo_clean = df_mlo[['Age_BP', 'co2']].copy()

# ==========================================
# 2. PARSE THE ANTARCTIC CO2 COMPOSITE (Bereiter et al. 2015)
# ==========================================
# Over the Holocene the composite is carried by the high-resolution Law Dome
# cores, giving ~330 year spacing - fine enough to resolve the flat Holocene
# plateau and the industrial upturn.
composite_path = "data/antarctica2015co2composite-noaa.txt"
df_ice_co2 = pd.read_csv(composite_path, comment='#', sep='\t').rename(
    columns={'age_gas_calBP': 'Age_BP', 'co2_ppm': 'co2'}
)[['Age_BP', 'co2']].dropna()

# ==========================================
# 3. PARSE HOLOCENE GLOBAL TEMPERATURE (Temp12k, Kaufman et al. 2020)
# ==========================================
# Multi-method global mean surface temperature ensemble, 0-12 kyr BP in 100 yr
# bins, with 5th/95th percentiles. A global reconstruction rather than a single
# Antarctic core, so it needs no amplification scaling.
# https://www.ncei.noaa.gov/pub/data/paleo/reconstructions/kaufman2020/
temp12k_path = "data/temp12k_allmethods_percentiles.csv"
df_holo = pd.read_csv(temp12k_path, skipinitialspace=True).rename(
    columns={'ages': 'Age_BP', 'global_median': 'temp_50',
             'global_5': 'temp_05', 'global_95': 'temp_95'}
)[['Age_BP', 'temp_05', 'temp_50', 'temp_95']]
df_holo = df_holo.sort_values(by="Age_BP", ascending=False)

# ==========================================
# 4. PARSE MODERN INSTRUMENTAL TEMPERATURE (HadCRUT5, monthly)
# ==========================================
# 'Time' comes as 'YYYY-MM', converted to a decimal year at the middle of the
# month, matching the convention Mauna Loa's own decimal_date already uses.
hadcrut_path = "data/hadcrut5-monthly-global.csv"
df_mod = pd.read_csv(hadcrut_path).rename(columns={'Anomaly (deg C)': 'temp_anomaly'})
y, m = df_mod['Time'].str.split('-', expand=True).astype(int).to_numpy().T
df_mod['year'] = y + (m - 0.5) / 12.0
df_mod = df_mod[['year', 'temp_anomaly']]
# Shift from the native 1961-1990 baseline onto preindustrial, so the
# instrumental curve and the Temp12k reconstruction share a zero.
baseline = df_mod.loc[
    df_mod['year'].between(PREINDUSTRIAL_FROM, PREINDUSTRIAL_TO), 'temp_anomaly'
].mean()
df_mod['temp_anomaly'] = df_mod['temp_anomaly'] - baseline
df_mod['Age_BP'] = 1950 - df_mod['year']
df_mod = df_mod.sort_values(by="Age_BP", ascending=False)

# ==========================================
# 5. SPLICE THE CO2 RECORD (one source per age interval)
# ==========================================
MLO_START_AGE = df_mlo_clean["Age_BP"].max()
df_co2_total = pd.concat(
    [df_ice_co2[df_ice_co2["Age_BP"] > MLO_START_AGE], df_mlo_clean],
    ignore_index=True,
)
df_co2_total = df_co2_total.sort_values(by="Age_BP", ascending=False)

# Clip every series to the plotted window. Both source records run far past it
# (the composite to 806 kyr, Temp12k into the deglacial at 12 kyr), and
# matplotlib autoscales the y-axis over all data, not just the visible slice -
# leaving both panels squashed into their top half. Keep a small margin so the
# curves still run to the left-hand edge.
window = HOLOCENE_START_BP * 1.05
df_co2_total = df_co2_total[df_co2_total["Age_BP"] <= window]
df_holo = df_holo[df_holo["Age_BP"] <= window]

# ==========================================
# 6. GENERATE THE STACKED SUBPLOT LAYOUT
# ==========================================
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
fig.subplots_adjust(hspace=0.15)

# --- TOP PLOT: Carbon Dioxide ---
color_co2 = 'tab:blue'
ax1.plot(df_co2_total["Age_BP"], df_co2_total["co2"], color=color_co2, lw=1.8,
         label="Antarctic ice cores + Mauna Loa")
ax1.set_ylabel("Carbon Dioxide (ppm)", fontsize=11, fontweight="bold")
ax1.grid(True, linestyle="--", alpha=0.3)
ax1.legend(loc="upper left", fontsize=8)

current_co2 = df_mlo['co2'].iloc[-1]
current_year = df_mlo['year'].iloc[-1]
ax1.scatter(df_mlo_clean['Age_BP'].min(), current_co2, color="tab:red", s=30, zorder=5)
ax1.annotate(
    f"{current_co2:.1f} ppm ({current_year})",
    xy=(df_mlo_clean['Age_BP'].min(), current_co2),
    xytext=(2200, current_co2 - 15),
    arrowprops=dict(facecolor='black', arrowstyle='->', lw=0.8),
    fontweight='bold', color='tab:red', fontsize=9
)

# --- BOTTOM PLOT: Global Temperature Anomaly ---
color_temp = 'tab:orange'
ax2.plot(df_holo["Age_BP"], df_holo["temp_50"], color=color_temp, lw=1.8,
         label="Temp12k reconstruction (global mean, median)")
ax2.fill_between(df_holo["Age_BP"], df_holo["temp_05"], df_holo["temp_95"],
                 color=color_temp, alpha=0.25, label="5th-95th percentile")
ax2.plot(df_mod["Age_BP"], df_mod["temp_anomaly"], color="tab:red", lw=1.8, zorder=4,
         label="HadCRUT5 instrumental (global mean)")
ax2.set_ylabel("Temperature Anomaly (°C,\nvs. 1850-1900)", fontsize=11, fontweight="bold")
# The records run past 1950, so say how far rather than leaving the reader
# to assume the plot stops at the zero of its own x-axis.
latest_year = int(max(df_mlo["year"].max(), df_mod["year"].max()))
ax2.set_xlabel(f"Years Before Present (Present = 1950; instrumental record to {latest_year})",
               fontsize=11, fontweight="bold")
ax2.grid(True, linestyle="--", alpha=0.3)
ax2.legend(loc="lower left", fontsize=8)

current_temp = df_mod['temp_anomaly'].iloc[-1]
current_temp_year = int(df_mod['year'].iloc[-1])
ax2.scatter(df_mod['Age_BP'].min(), current_temp, color="tab:red", s=30, zorder=5)
ax2.annotate(
    f"+{current_temp:.2f} °C ({current_temp_year})",
    xy=(df_mod['Age_BP'].min(), current_temp),
    xytext=(2200, current_temp - 0.25),
    arrowprops=dict(facecolor='black', arrowstyle='->', lw=0.8),
    fontweight='bold', color='tab:red', fontsize=9
)

# Global Layout Configuration & Geological Constraints
# ax1.set_title("12,000 Years of Holocene Stability - and the Industrial Departure",
#               fontsize=13, fontweight="bold", pad=15)

ax2.set_xlim(HOLOCENE_START_BP, -100)
ax2.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))

plt.tight_layout()

# Save pristine print assets for your KDP manuscript template
plt.savefig("images/co2_vs_temperature_stacked_12k.png", dpi=300, bbox_inches="tight")
plt.show()
