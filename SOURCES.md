# Sources

The data files in `data/` and the publications that the scripts in `src/` cite.
The entries are taken from the source records of the essays these analyses
belong to. Please cite the original authors when you use the data; the NOAA
files carry their own citation request in their headers.

## Data

### `data/PhanDA_GMSTandCO2_percentiles.csv`

GMST and CO<sub>2</sub> over the last 485 Myr, with percentiles. Used by the
485 Myr figures and the overview.

- Emily J. Judd, Jessica E. Tierney, Daniel J. Lunt, Isabel P. Montañez, Brian
  T. Huber, Scott L. Wing and Paul J. Valdes: PhanDA, <https://github.com/EJJudd/PhanDA>.
- Paper: *A 485-million-year history of Earth's surface temperature*, Science
  385(6716), 1428--1434, 2024. doi:[10.1126/science.adk3705](https://doi.org/10.1126/science.adk3705)

### `data/antarctica2015co2composite-noaa.txt`

Antarctic ice-core CO<sub>2</sub> composite over 800 kyr. Over the Holocene and
the industrial era the composite is carried by the Law Dome record.

- Bernhard Bereiter, Sarah Eggleston, Jochen Schmitt, Christoph Nehrbass-Ahles,
  Thomas F. Stocker, Hubertus Fischer, Sepp Kipfstuhl and Jerome Chappellaz:
  *Antarctic Ice Cores Revised 800KYr CO2 Data (Revision of the EPICA Dome C
  CO2 record from 800 to 600 kyr before present)*, NOAA National Centers for
  Environmental Information, 2015. Paper in Geophysical Research Letters,
  doi:[10.1002/2014GL061957](https://doi.org/10.1002/2014GL061957).
- <https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/antarctica2015co2composite.txt>

### `data/edc3deuttemp2007.txt`

EPICA Dome C deuterium and Antarctic temperature estimate over 800 kyr (EDC3
age scale). Used by the 800 kyr figures.

- Jean Jouzel, Valérie Masson-Delmotte et al.: *EPICA Dome C Ice Core 800KYr
  Deuterium Data and Temperature Estimates (Orbital and Millennial Antarctic
  Climate Variability over the Past 800,000 Years)*, NOAA National Centers for
  Environmental Information. Paper in Science 317(5839), 793--797, 2007.
- <https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/epica_domec/edc3deuttemp2007.txt>

### `data/temp12k_allmethods_percentiles.csv`, `data/readme-kaufman2020gmst.txt`

Holocene global mean surface temperature over the last 12 kyr, from all
reconstruction methods, with percentiles. The readme is NOAA's description of
the files.

- Darrell Kaufman, Nicholas McKay, Cody Routson, Michael Erb, Christoph
  Dätwyler, Philipp Sommer, Oliver Heiri and Basil Davis: *Temp12k: Holocene
  global mean surface temperature, a multi-method reconstruction approach*,
  Scientific Data 7, 2020. doi:[10.1038/s41597-020-0530-7](https://doi.org/10.1038/s41597-020-0530-7)
- <https://www.ncei.noaa.gov/pub/data/paleo/reconstructions/kaufman2020/>

### `data/mauna-loa-co2-mm.txt`

Monthly mean CO<sub>2</sub> at Mauna Loa, from 1958. A snapshot of a file that
NOAA updates monthly.

- Xin Lan, Pieter Tans and Kirk W. Thoning: *Trends in Atmospheric Carbon
  Dioxide: Mauna Loa monthly mean CO2*, NOAA Global Monitoring Laboratory.
- <https://gml.noaa.gov/ccgg/trends/data.html>
- The record was begun by Charles David Keeling, <https://keelingcurve.ucsd.edu/>.

### `data/hadcrut5-monthly-global.csv`

Instrumental global monthly mean near-surface temperature anomaly from 1850.
A snapshot of a file that the Met Office updates monthly.

- Colin P. Morice, John J. Kennedy, Nick A. Rayner et al.: *HadCRUT5 analysis:
  global monthly mean near-surface temperature anomaly (v.5.2.0.0)*, Met Office
  Hadley Centre and Climatic Research Unit, 2021.
- <https://www.metoffice.gov.uk/hadobs/hadcrut5/data/HadCRUT.5.2.0.0/download.html>

## Publication cited in the code

- Nicolas Caillon, Jeffrey P. Severinghaus, Jean Jouzel, Jean-Marc Barnola,
  Jiancheng Kang and Vladimir Y. Lipenkov: *Timing of Atmospheric CO2 and
  Antarctic Temperature Changes Across Termination III*, Science 299(5613),
  1728--1731, 2003. doi:[10.1126/science.1078758](https://doi.org/10.1126/science.1078758).
  Cited in `src/plot_crosscorrelation.py` for the lag of CO<sub>2</sub> behind
  Antarctic temperature across a deglaciation.
