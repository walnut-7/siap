# siap
Code for reproducing the results in "Matrix Factorization-Based Solar Spectral Irradiance Missing Data Imputation with Uncertainty Quantification".

# Restore environment

0. Please check the system requirements in `session_info.txt` to ensure that your system has sufficient memory and storage capacity to run the experiments. 
If the hardware requirements are not met, we also provide intermediate results that allow direct reproduction of the figures and tables.

1. Navigate to the project folder 'your/path/to/siap/' in terminal.

2. Run `bash setup.sh`, which prompt you to install the required tools and, thereafter, prepare the R environment. 

We use [`renv`](https://rstudio.github.io/renv/) to manage and restore the R environment. The .Rprofile file (which sources `renv/activate.R`) is sourced automatically by R when a new R session starts in the project directory.
After this, `renv` would specify a project-level library `renv/library`, where the R packages are/will be installed. 
To install packages globally, you can temporarily deactivate the project and re-activate it after installing the package. 
```
renv::deactivate() # temporarily deactivate
globallib <- .libPaths()
install.packages("mypackage") # install packages globally
renv::activate() # re-activate 
.libPaths(c(.libPaths(), globallib)) # enable R to look in the global library 
```
For more information about `renv`, see this [link](https://rstudio.github.io/renv/).

# Reproduce the results
Set the project root directory `your/path/siap` as the working directory for R.
Running the following commands will generate figures and tables in the `output/simulation` and `output/realdata` directories:
```
cd your/path/siap
Rscript code/output_simulation.R
Rscript code/output_realdata.R
```

The experiments were originally conducted on slurm cluster, managed by R package [`batchtools`](https://github.com/mlr-org/batchtools).
However, the script is also adapted to local machine. 
If you are running the scripts on clusters, please change the account name `account` accordingly in `code/simulation.R` and `code/realdata.R`. 
We also provide scripts (`code/simulation_submit_marss_batches.R` and `code/simulation_submit_marss_batches.sh`) to submit 4000 MARSS jobs in simulation study every 2 hours, since the total number of jobs will likely exceed the maximum job number allowance.

# Data citation and availability
The data used in the simulation study and SSI reconstruction analysis are publicly available. We list these data and sources below.

1. [TSIS-1 SIM SSI](https://disc.gsfc.nasa.gov/datacollection/TSIS_SSI_L3_24HR_13.html) (used in SSI reconstruction section)
- Coddington, O. M., E. C. Richard, D. Harber, P. Pilewskie, T. N. Woods, M. Snow, K. Chance, X. Liu, and K. Sun. 2023. “Version 2 of the TSIS-1 Hybrid Solar Reference Spectrum and Extension to the Full Spectrum.” Earth and Space Science 10 (3): e2022EA002637. https://doi.org/10.1029/2022EA002637.
- Richard, Erik, Odele Coddington, Dave Harber, Michael Chambliss, Steven Penton, Keira Brooks, Luke Charbonneau, et al. 2024. “Advancements in Solar Spectral Irradiance Measurements by the TSIS-1 Spectral Irradiance Monitor and Its Role for Long-Term Data Continuity.” Journal of Space Weather and Space Climate 14:10. https://doi.org/10.1051/swsc/2024008.
- Richard, Erik, Dave Harber, Odele Coddington, Ginger Drake, Joel Rutkowski, Matthew Triplett, Peter Pilewskie, and Tom Woods. 2020. “SI-Traceable Spectral Irradiance Radiometric Characterization and Absolute Calibration of the TSIS-1 Spectral Irradiance Monitor (SIM).” Remote Sensing 12 (11): 1818. https://doi.org/10.3390/rs12111818.

2. Synthetic SSI, v3.2, accessed Apr 18 2023, generated based on [CMIP6, v3.2](https://www.solarisheppa.kit.edu/85.php) record (used in simulation section)
- Matthes, Katja, Bernd Funke, Monika E. Andersson, Luke Barnard, Jürg Beer, Paul Charbonneau, Mark A. Clilverd, et al. 2017. “Solar Forcing for CMIP6 (v3.2).” Geoscientific Model Development 10 (6): 2247–2302. https://doi.org/10.5194/gmd-10-2247-2017.

3. [CSIM SSI](https://lasp.colorado.edu/csim/data-and-ham-radio/) (used in SSI reconstruction section)
- Richard, Erik, Dave Harber, Ginger Drake, Joel Rutkowsi, Zach Castleman, Matthew Smith, Jacob Sprunck, et al. 2019. “Compact Spectral Irradiance Monitor Flight Demonstration Mission.” In CubeSats and SmallSats for Remote Sensing III, 11131:15–34. SPIE. https://doi.org/10.1117/12.2531268.

These data have been made available at the [this URL](https://figshare.com/s/870cc558b0114dc16378).

