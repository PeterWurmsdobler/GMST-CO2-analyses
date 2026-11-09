import pandas as pd
import matplotlib.pyplot as plt

# The eight glacial terminations, at the age of peak warming rate in the
# low-pass filtered (250-year cutoff) EDC temperature record - the same
# calculation plot_derivatives.py's rates_800k figure uses, so the markers
# here line up with the peaks in that figure rather than an independently
# sourced literature age that need not agree with this ice core's own age
# scale.
TERMINATIONS_YR_BP = (11438, 134938, 244138, 335538, 425438, 522438, 621438, 715038)

# Vostok/EDC anomalies are LOCAL Antarctic temperatures. Antarctic swings run
# roughly twice the global-mean amplitude, so divide the ice-core curve by this
# polar amplification factor to put it on the same scale as the modern
# instrumental record. State the factor in the caption; estimates span ~1.5-2.
POLAR_AMPLIFICATION = 2.0

# Both temperature curves in the bottom panel are put on the IPCC preindustrial
# baseline, the same one plot_holocene.py uses, so the two figures annotate the
# same present-day anomaly instead of differing by the width of that baseline.
PREINDUSTRIAL_FROM, PREINDUSTRIAL_TO = 1850, 1900


def load_noaa_two_column_data(file_path, age_col_idx, data_col_idx):
    """
    Parses NOAA paleoclimate text files line-by-line. Extracts a uniform,
    clean 2-column layout (Age and the targeted variable) by explicitly
    checking row index boundaries before parsing.
    """
    clean_pairs = []
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            # These NOAA files are tab-delimited and may carry empty cells
            # (missing values) plus trailing free-text columns, so keep the
            # column positions intact and only parse the two cells we need.
            tokens = line.rstrip('\n').split('\t') if '\t' in line else line.split()
            max_required_idx = max(age_col_idx, data_col_idx)
            if len(tokens) <= max_required_idx:
                continue
            try:
                age_val = float(tokens[age_col_idx])
                data_val = float(tokens[data_col_idx])
            except ValueError:
                continue
            clean_pairs.append([age_val, data_val])

    df = pd.DataFrame(clean_pairs, columns=['Age_BP', 'Value'])
    return df.dropna()


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
# The single-core files (antarctica2015co2edc/vostok) carry large gaps where
# other cores cover the interval - EDC alone is missing 156-563 kyr BP. The
# published composite splices all Antarctic cores onto one gas chronology and
# runs gap-free from -51 to 806 kyr BP, so use that rather than hand-splicing.
# https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/antarctica2015co2composite.txt
composite_path = "data/antarctica2015co2composite-noaa.txt"
# Format: [Gas_Age, CO2, CO2_1s] -> Gas_Age = 0, CO2 = 1
df_ice_co2 = load_noaa_two_column_data(composite_path, age_col_idx=0, data_col_idx=1)
df_ice_co2 = df_ice_co2.rename(columns={'Value': 'co2'})

# ==========================================
# 3. PARSE EPICA DOME C TEMPERATURE (EDC3, Jouzel et al. 2007)
# ==========================================
# Vostok's deutnat.txt stops at 423 kyr BP; EDC3 reaches 802 kyr BP and so
# matches the span of the CO2 composite. NOTE: gt4nat-noaa.txt is the GT4 age
# scale only (depth, ice age, gas age) and holds no temperature at all.
# https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/epica_domec/edc3deuttemp2007.txt
edc_temp_path = "data/edc3deuttemp2007.txt"
# Format: [Bag, ztop, Age, Deuterium, Temperature] -> Age = 2, Temperature = 4
df_temp = load_noaa_two_column_data(edc_temp_path, age_col_idx=2, data_col_idx=4)
df_temp_clean = df_temp.rename(columns={'Value': 'temp_anomaly'})
df_temp_clean['temp_anomaly'] = df_temp_clean['temp_anomaly'] / POLAR_AMPLIFICATION
df_temp_clean = df_temp_clean.sort_values(by="Age_BP", ascending=False)

# ==========================================
# 4. PARSE MODERN INSTRUMENTAL TEMPERATURE (HadCRUT5, monthly)
# ==========================================
# https://www.metoffice.gov.uk/hadobs/hadcrut5/ - monthly global mean anomaly.
# 'Time' comes as 'YYYY-MM', converted to a decimal year at the middle of the
# month, matching the convention Mauna Loa's own decimal_date already uses.
hadcrut_path = "data/hadcrut5-monthly-global.csv"
df_mod = pd.read_csv(hadcrut_path).rename(columns={'Anomaly (deg C)': 'temp_anomaly'})
y, m = df_mod['Time'].str.split('-', expand=True).astype(int).to_numpy().T
df_mod['year'] = y + (m - 0.5) / 12.0
df_mod = df_mod[['year', 'temp_anomaly']]
# Two re-referencings, both needed to make the panel read correctly.
#
# HadCRUT5 is natively referenced to 1961-1990; shift it onto the IPCC
# preindustrial baseline so this figure and the 12 kyr one agree on what the
# present-day anomaly is.
preindustrial = df_mod.loc[
    df_mod['year'].between(PREINDUSTRIAL_FROM, PREINDUSTRIAL_TO), 'temp_anomaly'
].mean()
df_mod['temp_anomaly'] = df_mod['temp_anomaly'] - preindustrial

# The ice-core anomaly is expressed against present-day Antarctic conditions,
# i.e. its zero sits at roughly the mid-20th century rather than at 1850-1900.
# Shift it by the same interval so both curves in this panel share one zero and
# meet without a step at the join.
mid_century = df_mod.loc[df_mod['year'].between(1940, 1960), 'temp_anomaly'].mean()
df_temp_clean['temp_anomaly'] = df_temp_clean['temp_anomaly'] + mid_century
df_mod['Age_BP'] = 1950 - df_mod['year']
df_mod = df_mod.sort_values(by="Age_BP", ascending=False)

# ==========================================
# 5. SPLICE THE CO2 RECORD (one source per age interval)
# ==========================================
# Mauna Loa from 1958 onwards, the ice-core composite for everything older.
MLO_START_AGE = df_mlo_clean["Age_BP"].max()
df_co2_total = pd.concat(
    [df_ice_co2[df_ice_co2["Age_BP"] > MLO_START_AGE], df_mlo_clean],
    ignore_index=True,
)
df_co2_total = df_co2_total.sort_values(by="Age_BP", ascending=False)

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
for age in TERMINATIONS_YR_BP:
    ax1.axvline(age, color="firebrick", ls="--", lw=1.0, alpha=0.8, zorder=1)
ax1.plot([], [], color="firebrick", ls="--", lw=1.0, alpha=0.8, label="Glacial terminations")
ax1.legend(loc="upper left", fontsize=8)

# Highlight the modern industrial vertical CO2 spike
current_co2 = df_mlo['co2'].iloc[-1]
current_year = df_mlo['year'].iloc[-1]
ax1.scatter(df_mlo_clean['Age_BP'].min(), current_co2, color="tab:red", s=30, zorder=5)
ax1.annotate(
    f"Current Spike: {current_co2:.1f} ppm ({current_year})",
    xy=(df_mlo_clean['Age_BP'].min(), current_co2),
    xytext=(300000, current_co2 - 40),
    arrowprops=dict(facecolor='black', arrowstyle='->', lw=0.8),
    fontweight='bold', color='tab:red', fontsize=9
)

# --- BOTTOM PLOT: Temperature Anomaly ---
color_temp = 'tab:orange'
ax2.plot(df_temp_clean["Age_BP"], df_temp_clean["temp_anomaly"], color=color_temp, lw=1.5, alpha=0.8,
         label=f"EPICA Dome C $\\Delta$T / {POLAR_AMPLIFICATION:g} (global-mean equivalent)")
ax2.plot(df_mod["Age_BP"], df_mod["temp_anomaly"], color="tab:red", lw=1.8, zorder=4,
         label="HadCRUT5 instrumental (global mean)")
ax2.set_ylabel("Temperature Anomaly (°C, global-mean\nequivalent, vs. 1850-1900)", fontsize=11, fontweight="bold")
# The records run past 1950, so say how far rather than leaving the reader
# to assume the plot stops at the zero of its own x-axis.
latest_year = int(max(df_mlo["year"].max(), df_mod["year"].max()))
ax2.set_xlabel(f"Years Before Present (Present = 1950; instrumental record to {latest_year})",
               fontsize=11, fontweight="bold")
ax2.grid(True, linestyle="--", alpha=0.3)
for age in TERMINATIONS_YR_BP:
    ax2.axvline(age, color="firebrick", ls="--", lw=1.0, alpha=0.8, zorder=1)
ax2.legend(loc="lower left", fontsize=8)

# Mirror the CO2 panel: mark the modern warming that is one pixel wide at this scale
current_temp = df_mod['temp_anomaly'].iloc[-1]
current_temp_year = int(df_mod['year'].iloc[-1])
ax2.scatter(df_mod['Age_BP'].min(), current_temp, color="tab:red", s=30, zorder=5)
ax2.annotate(
    f"Current Spike: +{current_temp:.2f} °C ({current_temp_year})",
    xy=(df_mod['Age_BP'].min(), current_temp),
    xytext=(300000, current_temp - 0.8),
    arrowprops=dict(facecolor='black', arrowstyle='->', lw=0.8),
    fontweight='bold', color='tab:red', fontsize=9
)

# Global Layout Configuration & Geological Constraints
# ax1.set_title("800,000-Year Climate Co-Evolution: Atmospheric $CO_2$ vs. Antarctic Temperature",
#               fontsize=13, fontweight="bold", pad=15)

# Span the full available record across both panels
max_historical_age = int(max(df_co2_total["Age_BP"].max(), df_temp_clean["Age_BP"].max()))
ax2.set_xlim(max_historical_age, -500)
ax2.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))

plt.tight_layout()

# Save pristine print assets for your KDP manuscript template
plt.savefig("images/co2_vs_temperature_stacked_800k.png", dpi=300, bbox_inches="tight")
plt.show()
