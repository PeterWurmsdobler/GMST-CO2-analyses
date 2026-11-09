import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# A standalone recreation of Judd et al.'s own GMST-only plot, built from their
# published PhanDA data rather than reproducing their figure directly (a
# licensing, not a scientific, distinction): this is the temperature half of
# the combined CO2/GMST figure plot_gmst.py draws later for the 485-million-
# year window, shown here on its own since the text at this point is making a
# claim about GMST alone.
#
# Judd et al.'s own figure carries a "climate state" header bar, a discrete
# cooler/warmer classification shaded by how confident that classification is.
# The strip below echoes the same idea the text's own "hothouse" wording
# needs, a glance at whether the Earth was in a warm or cold state at any
# point in this record, but renders it differently: a continuous colour
# gradient taken directly from the median GMST curve itself, not a discrete
# category with its own separate confidence shading.
df = pd.read_csv("data/PhanDA_GMSTandCO2_percentiles.csv")
df = df.sort_values(by="AverageAge", ascending=False)
age = df["AverageAge"].to_numpy()
gmst = df["GMST_50"].to_numpy()

fig = plt.figure(figsize=(10, 5.5))
gs = GridSpec(2, 1, height_ratios=[1, 14], hspace=0.06)
ax_strip = fig.add_subplot(gs[0])
ax = fig.add_subplot(gs[1], sharex=ax_strip)

ax_strip.imshow(gmst[np.newaxis, :], aspect="auto", cmap="coolwarm",
                extent=[age.max(), age.min(), 0, 1],
                vmin=gmst.min(), vmax=gmst.max())
ax_strip.set_yticks([])
ax_strip.tick_params(axis="x", which="both", labelbottom=False, length=0)
ax_strip.set_ylabel("Climate\nstate", fontsize=8, rotation=0, ha="right", va="center")
for spine in ax_strip.spines.values():
    spine.set_visible(False)

ax.plot(age, gmst, color="black", lw=1.8, label="Median GMST")
ax.fill_between(age, df["GMST_16"], df["GMST_84"],
                color="tab:red", alpha=0.3, label="16th–84th percentile")
ax.fill_between(age, df["GMST_05"], df["GMST_95"],
                color="tab:red", alpha=0.1, label="5th–95th percentile")

ax.set_ylabel("GMST (°C)", fontsize=11, fontweight="bold")
ax.set_xlabel("Age (Millions of Years Ago)", fontsize=11, fontweight="bold")
ax.grid(True, linestyle="--", alpha=0.5)
# Anchored at a data point, around 300 Ma, rather than a figure corner: every
# corner is crossed by the curve or its bands at some point in this record,
# while the stretch around 300 Ma sits well below the top of the axis.
ax.legend(loc="upper center", bbox_to_anchor=(300, 44), bbox_transform=ax.transData)

# Geological convention: past on the left, present on the right.
ax.set_xlim(age.max(), age.min())

plt.savefig("images/phanda_temperature.png", dpi=300, bbox_inches="tight")
print("wrote images/phanda_temperature.png")
