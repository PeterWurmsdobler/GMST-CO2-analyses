# GMST-CO2-analyses

The calculations and plots behind a series of three essays on the two Earth
variables atmospheric CO<sub>2</sub> and global mean surface temperature (GMST):

1. *A Brief History of Two Earth Variables*: the values themselves;
2. *The Rate of Change of Two Earth Variables*: their derivatives;
3. *The Correlation Between Two Earth Variables*: which one moves first.

Each essay looks at the same two records through four nested windows (485 Myr,
800 kyr, 12 kyr and 200 yr) and, in addition, through one overview on a
logarithmic time axis. Every figure comes from one script in this repository,
with the data it needs checked in alongside.

## Setting up

You need Python 3.10 or newer.

```sh
git clone https://github.com/PeterWurmsdobler/GMST-CO2-analyses.git
cd GMST-CO2-analyses
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

The requirements are numpy, pandas, scipy and matplotlib.

## Running

Run the scripts from the repository root, because the data in `data/` is read
from, and the figures written to `images/`, with paths relative to it. Each script writes its
PNG (300 dpi) into `images/`. The scripts also call `plt.show()` or print
numbers; to run without opening windows, select a non-interactive backend:

```sh
export MPLBACKEND=Agg
python src/plot_gmst.py
```

To make every figure:

```sh
for f in src/plot_*.py; do python "$f"; done
```

All scripts together take about half a minute.

## The plot scripts

| Script | Figures written to `images/` | What it does |
|---|---|---|
| `src/plot_phanda_temperature.py` | `phanda_temperature.png` | GMST alone over the Phanerozoic, redrawn from the PhanDA data: median with 16th--84th and 5th--95th percentile bands, with a colour strip taken from the median curve to show the warm and cold states. |
| `src/plot_gmst.py` | `co2_vs_temperature_stacked_485Ma.png` | CO<sub>2</sub> and GMST from PhanDA over the last 485 Myr, one panel above the other, with the five mass extinctions marked. |
| `src/plot_co2.py` | `co2_vs_temperature_stacked_800k.png` | The last 800 kyr: the EPICA Dome C temperature and the CO<sub>2</sub> ice-core composite, joined to Mauna Loa and HadCRUT5. The Antarctic temperature is divided by a polar amplification factor (2.0) to approach a global mean and is put on the 1850--1900 baseline. The eight glacial terminations are marked. |
| `src/plot_holocene.py` | `co2_vs_temperature_stacked_12k.png` | The last 12 kyr: Temp12k global mean temperature and CO<sub>2</sub>, joined to the instrumental records, both on the 1850--1900 baseline. |
| `src/plot_industrial.py` | `co2_vs_temperature_stacked_200y.png` | The industrial era, 1850 onwards: Mauna Loa and ice-core CO<sub>2</sub> with HadCRUT5, raw and low-pass filtered (5, 10 and 20 year cutoffs). |
| `src/plot_overview.py` | `co2_vs_temperature_overview.png` | All four windows spliced on one logarithmic, reversed time axis from 485 Myr to 10 yr before present, each variable's value above its own rate of change. Temperatures are shifted onto the PhanDA absolute scale by one constant. |
| `src/plot_derivatives.py` | `co2_vs_temperature_rates_485Ma.png`, `_rates_800k.png`, `_rates_12k.png`, `_rates_200y.png` | The signed rate of change per century of both variables in each window. Records are low-pass filtered (zero phase) before differentiating, with cutoffs chosen to keep 80 % of each record's rate variance. Prints the present-day rates, the fastest natural rate each window resolves and how many times the present rate exceeds it, and the rate dispersion over the last 100 kyr. |
| `src/plot_rates_overview.py` | `co2_vs_temperature_rates_overview.png` | The rates only, on the same logarithmic axis as `plot_overview.py`. |
| `src/plot_crosscorrelation.py` | `co2_vs_temperature_xcorr_485Ma.png`, `_xcorr_800k.png`, `_xcorr_12k.png`, `_xcorr_200y.png`, `_xcorr_overview.png`, `_xcorr_200y_seasonal.png` | Cross-correlation of CO<sub>2</sub> and GMST. Both records go onto one uniform timeline and are high-pass filtered to remove the shared drift; the correlation is then computed at a range of lags. A positive lag means CO<sub>2</sub> leads temperature. Prints the lag at the peak for each window; the lag cannot be resolved below the sampling period. |
| `src/plot_crossspectrum_800k.py` | `co2_vs_temperature_crossspectrum_800k.png` | The cross-spectrum of CO<sub>2</sub> and Antarctic temperature over 800 kyr, as the Fourier transform of the cross-correlation function (Blackman--Tukey), for periods of 15 to 400 kyr. The phase is repeated on each half of the record as a check. Prints a report per orbital band. |

`src/lowpass.py` is not run by itself. It holds the zero-phase Butterworth filter,
the resampling onto a uniform timeline and the filtered derivative that the
scripts above share.

The docstring at the top of each script gives its reasoning in full: the
sign convention of the lags, the choice of cutoffs and the caveats.

## Data

The data are in `data/`: PhanDA (485 Myr), the EPICA Dome C temperature and the
Antarctic CO<sub>2</sub> composite (800 kyr), Temp12k (12 kyr), Mauna Loa and
HadCRUT5 (instrumental). Authors, titles and links are in [SOURCES.md](SOURCES.md).

The Mauna Loa and HadCRUT5 files are snapshots, so the present-day rates
printed by `src/plot_derivatives.py` depend on the date the files were
downloaded; replace them with newer downloads to bring the figures up to date.

## Notes

- Matplotlib and numpy versions differ a little in rendering, so the PNGs may
  not match the published figures pixel for pixel.
- The `images/` directory holds the output only; the PNGs are not checked in.
- There is no `pyproject.toml`: the scripts are run, not installed, so
  `requirements.txt` is all they need.
