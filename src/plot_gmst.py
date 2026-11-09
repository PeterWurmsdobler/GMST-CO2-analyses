import pandas as pd
import matplotlib.pyplot as plt

# The five mass extinctions of the Phanerozoic, in Ma. Kept in step with the
# same list in plot_derivatives.py, so the absolute-value figure and the rate
# figure mark the same events at the same ages.
MASS_EXTINCTIONS = [
    (445, "Ordovician–Silurian"),
    (372, "Late Devonian"),
    (252, "Permian–Triassic"),
    (201, "Triassic–Jurassic"),
    (66, "Cretaceous–Paleogene"),
]


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


# 1. Load the data (replace 'phanda_data.csv' with your actual file path)
file_path = "data/PhanDA_GMSTandCO2_percentiles.csv"
df = pd.read_csv(file_path)

# Ensure the data is sorted by age so the fill_between shading maps correctly
df = df.sort_values(by="AverageAge", ascending=False)

# 2. Set up the figure with 2 stacked subplots sharing the x-axis
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 8), sharex=True)
fig.subplots_adjust(hspace=0.1)  # Bring subplots close together

# --- TOP PLOT: Atmospheric CO2, matching the CO2-on-top convention used
#     throughout this set of figures ---
# Median line
ax1.plot(df["AverageAge"], df["CO2_50"], color="blue", label="Median CO2", lw=1.8)
# 16th to 84th percentile uncertainty
ax1.fill_between(df["AverageAge"], df["CO2_16"], df["CO2_84"],
                 color="tab:blue", alpha=0.3, label="16th–84th Percentile")
# 5th to 95th percentile uncertainty
ax1.fill_between(df["AverageAge"], df["CO2_05"], df["CO2_95"],
                 color="tab:blue", alpha=0.1, label="5th–95th Percentile")

ax1.set_ylabel("Atmospheric $CO_2$ (ppm)", fontsize=11, fontweight="bold")
ax1.grid(True, linestyle="--", alpha=0.5)
# Leave headroom for the extinction labels, which hang from the top of the
# axis: a fixed offset would only suit one variable's range, so the margin is
# a fraction of this panel's own autoscaled range instead.
ylo, yhi = ax1.get_ylim()
ax1.set_ylim(top=yhi + 0.35 * (yhi - ylo))
mark_extinctions(ax1, label_axis=True)
# Anchored at a specific data point (age 200 Ma, 2000 ppm) rather than a
# named corner: every corner of this panel is occupied by either the curve
# itself, its shaded uncertainty bands, or a hanging extinction label at some
# point in its history, so a fixed corner keeps drifting back onto one of
# them whenever the data changes. This point sits in the long, low stretch
# between the Permian-Triassic and Triassic-Jurassic extinctions, well above
# the curve and its bands there and well clear of either label.
ax1.legend(loc="center left", bbox_to_anchor=(200, 2300), bbox_transform=ax1.transData)
# ax1.set_title("Phanerozoic Climate Evolution (PhanDA Reconstruction)", fontsize=13, fontweight="bold", pad=15)

# --- BOTTOM PLOT: Global Mean Surface Temperature (GMST) ---
# Median line
ax2.plot(df["AverageAge"], df["GMST_50"], color="black", label="Median GMST", lw=1.8)
# 16th to 84th percentile uncertainty (Darker shading)
ax2.fill_between(df["AverageAge"], df["GMST_16"], df["GMST_84"],
                 color="tab:red", alpha=0.3, label="16th–84th Percentile")
# 5th to 95th percentile uncertainty (Lighter shading)
ax2.fill_between(df["AverageAge"], df["GMST_05"], df["GMST_95"],
                 color="tab:red", alpha=0.1, label="5th–95th Percentile")

ax2.set_ylabel("GMST (°C)", fontsize=11, fontweight="bold")
ax2.set_xlabel("Age (Millions of Years Ago)", fontsize=11, fontweight="bold")
ax2.grid(True, linestyle="--", alpha=0.5)
mark_extinctions(ax2)
ax2.legend(loc="lower left")

# --- GEOLOGICAL CONVENTION: Reverse X-Axis ---
# In geology, charts conventionally read from left (past, 485 Ma) to right (present day, 0 Ma)
ax2.set_xlim(df["AverageAge"].max(), df["AverageAge"].min())

# Save as a print-ready vector file (or 300 DPI PNG) for your KDP manuscript template
plt.savefig("images/co2_vs_temperature_stacked_485Ma.png", dpi=300, bbox_inches="tight")

# Display the plot window
plt.show()

