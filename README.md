# CONCLUDE

This repository includes the code for running the CONCLUDE model, designed to produce **spatial population projections** across the Mediterranean region at 30 arc seconds spatial resolution.

The **basic model** setup is described in Reimann et al. (2021): Accounting for internal migration in spatial population projections-a gravity-based modeling approach using the Shared Socioeconomic Pathways. In: Environ Res Lett 16, 074025. DOI: https://doi.org/10.1088/1748-9326/ac0b66.

The **extended model** setup accounting for sea-level rise-induced migration and adaptation policies is described in Reimann et al. (2023): Exploring spatial feedbacks between adaptation policies and internal migration patterns due to sea-level rise. In: Nature Communications 14, 2630. DOI: https://doi.org/10.1038/s41467-023-38278-y.


### Overview of model scripts

**00_preprocessing**
-	Calculate distance matrix for country / countries ('cntr') to be studied 
-	Create a matrix with the distances of the each point (i.e. center of a raster cell) to the next points (i.e. raster cells) up to the maximum distance of the gravity window 
-	Save as .mtx.gz 
-	Note: needed to be done on a hpc because a lot of RAM needed (at least 64 GB); allow up to 18h per country (for TUR) depending on the number of cells
-	Data needed: 
    -   Country 3-letter ISO code and numeric ID
    -   Spatial mask / other cntr dataset (to derive coordinates of each raster cell)
-	**Needs to be done only once and can be used in each following step**
 
**01_beta_calibration**
-	Run beta optimization function
-	Find the optimal beta value representing the observed changes in population patterns best (i.e. the beta that results in the lowest error) 
-	Load data needed for beta opt (per cntr) 
    -   distance matrix 
    -   total pop of base year
    -   values of spatial mask
    -   population change from base year to next time step in urban versus rural locations (i.e. 1990-2000; 2000-2010) 
    -   urb/rur pop of base year and base year + t
-	Run beta opt function for tolerance ('tol') = 1e-3, both calibration periods, urb versus rur
-	Write betas in csv file / save entire workspace

**02_alpha_calculation**
-	Calculate alphas per cell
-	calculate 'correction factor' per cell (i.e. where patterns produced by using the calibrated betas do not represent the observed values) 
-	load data needed 
    -   population per t
    -   spatial mask of 2000
    -   coastal mask
    -   established betas
-	Calculate alpha 
-	Cap per cntr:  use IQR = upper and lower 25 % removed
-	Scale per cntr: scale to values that can range -100-100
-	Convert into rasters + merge all cntr rasters 

**03_validation**
-	validate model based on calibrated betas and alphas
-	Run projections script for time steps of available observed data:  
    -   Use calibration of 1990-2000 and project to 2010 
    -   Use calibration of 2000-2010 and project to 2015 (less robust)
-	Load respective betas and alphas (use unscaled alphas)
-	Calculate error terms (raster-based and for entire cntr/region) 

**04_pop_proj_nSLR** (Reimann et al. 2021)
-	Produce population projections from 2020 to 2100 for each cntr and each ssp
-	Use population numbers (urb/rur/tot) per cntr from the SSPs 

**04_pop_proj_wSLR** (Reimann et al. 2023)
-	Produce population projections from 2020 to 2100 for three SSP-RCP combinations (i.e. SSP1-RCP2.6, SSP3-RCP4.5, SSP5-RCP8.5) and two adaptation policy scenarios (i.e. nA, wA).
-	No adaptation (nA) 
    -   Load inundation layers 
    -   'Pick up' population in submerged areas + set spatial mask to 0 
    -   Redistribute population picked up along with population allocation 
-	With adaptation (wA)
    -   Load inundation layers (accounting for protected land + retreat (SSP1 only)) 
    -   'Pick up' population in submerged areas + set spatial mask to 0 
    -   Use new spatial mask including setback zones (in SSPs 1 + 5) 
    -   Redistribute population picked up along with population allocation 

