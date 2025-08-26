# Data
The data used in the simulation study and SSI reconstruction analysis are publicly available. We list these data and sources below. To reproduce the results, please put the following data in the directory `siap/data` and name them according to

| Data           | Name it as |
| -------------- | ---------- |
| TSIS-1 SIM SSI | TSIS-1_ssi_tsi_20180314_20230129.nc |
| Synthetic SSI | ssi_synthetic_interpolated.nc |
| CSIM SSI | csim_ssi_L3_latest.nc |

1. [TSIS-1 SIM SSI, v13, accessed Mar 13 2023](https://disc.gsfc.nasa.gov/datacollection/TSIS_SSI_L3_24HR_13.html) (used in SSI reconstruction section)

  Citation to the data:
  - Richard, Erik. 2025. “TSIS SIM Level 3 Solar Spectral Irradiance 24-Hour Means V13.” NASA Goddard Earth Sciences Data and Information Services Center. https://doi.org/10.5067/TSIS/SIM/DATA326.

  Related publications:
  - Coddington, O. M., E. C. Richard, D. Harber, P. Pilewskie, T. N. Woods, M. Snow, K. Chance, X. Liu, and K. Sun. 2023. “Version 2 of the TSIS-1 Hybrid Solar Reference Spectrum and Extension to the Full Spectrum.” Earth and Space Science 10 (3): e2022EA002637. https://doi.org/10.1029/2022EA002637.
  - Richard, Erik, Odele Coddington, Dave Harber, Michael Chambliss, Steven Penton, Keira Brooks, Luke Charbonneau, et al. 2024. “Advancements in Solar Spectral Irradiance Measurements by the TSIS-1 Spectral Irradiance Monitor and Its Role for Long-Term Data Continuity.” Journal of Space Weather and Space Climate 14:10. https://doi.org/10.1051/swsc/2024008.
  - Richard, Erik, Dave Harber, Odele Coddington, Ginger Drake, Joel Rutkowski, Matthew Triplett, Peter Pilewskie, and Tom Woods. 2020. “SI-Traceable Spectral Irradiance Radiometric Characterization and Absolute Calibration of the TSIS-1 Spectral Irradiance Monitor (SIM).” Remote Sensing 12 (11): 1818. https://doi.org/10.3390/rs12111818.

2. Synthetic SSI, v3.2, accessed Apr 18 2023, generated based on [CMIP6, v3.2](https://doi.org/10.22033/ESGF/input4MIPs.1122) record (used in simulation section)

  Citation to the CMIP6 data:
  - Matthes, Katja, Bernd Funke, Tim Kruschke, and Sebastian Wahl. 2017. “Input4mips.Solaris-Heppa.Solar.Cmip.Solaris-HEPPA-3-2.” Earth System Grid Federation. https://doi.org/10.2203/ESGF/input4MIPs.1122.
  
  Related publications:
  - Matthes, Katja, Bernd Funke, Monika E. Andersson, Luke Barnard, Jürg Beer, Paul Charbonneau, Mark A. Clilverd, et al. 2017. “Solar Forcing for CMIP6 (v3.2).” Geoscientific Model Development 10 (6): 2247–2302. https://doi.org/10.5194/gmd-10-2247-2017.

3. [Compact Spectral Irradiance Monitor (CSIM) Level 3 Photodiode SSI, accessed Mar 21 2024](https://lasp.colorado.edu/csim/data-and-ham-radio/) (used in SSI reconstruction section)

  Citation to the data:
  - “Compact Spectral Irradiance Monitor (CSIM) Level 3 Photodiode SSI Data.” 2019. Version 1. Laboratory for Atmospheric and Space Physics (LASP). https://lasp.colorado.edu/csim/data-and-ham-radio/.
  
  Related publications:
  - Richard, Erik, Dave Harber, Ginger Drake, Joel Rutkowsi, Zach Castleman, Matthew Smith, Jacob Sprunck, et al. 2019. “Compact Spectral Irradiance Monitor Flight Demonstration Mission.” In CubeSats and SmallSats for Remote Sensing III, 11131:15–34. SPIE. https://doi.org/10.1117/12.2531268.
