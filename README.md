# Drivers of site fidelity in ungulates
Morrison TA, Merkle J, Hopcraft JGC, Aikens EO, Beck JL, Boone RB, Courtemanch A, Dwinnell SP, Fairbanks WS, Griffith B, Middleton AD, Monteith KL, Oates B, Riotte-Lambert L, Sawyer H, Smith K, Stabach JA, Taylor KL, Kauffman MJ. 2021. Drivers of site fidelity in ungulates. Journal of Animal Ecology, 90, 4, 955-966.
https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.13425

R code used in the analysis of site fidelity

Data files only include: (1) example raw GPS telemetry dataset from Serengeti wildebeest and zebra as an example, (2) derived site fidelity file

There are 2 main 'run' scripts: 
(1) sitefidelity.R calculates inter-year site fidelity from raw telemetry data and compiles all inidivduals into an aggregated dataframe, which is saved as ‘fidelity_toshare.csv’. It also calculates a number of covariates from MODIS ndvi data and from the spatial attributes of the telemetry data. Note this script takes quite a while to run. 
(2) analysis_toShare.R runs the statistical analyses and plotting of data that are included in the paper, using the aggregated dataset  ‘fidelity_toshare.csv’. Each row of this dataset is a unique individual, with site fidelity summarized across year1 and year2. 

If interested in the actual calculation of site fidelity, it is written as a function in the 'functions/iyd.R' file. 

If interested in predictability metrics, they are written as function in 'functions/CTime.R' and 'functions/Cspace.R' files. 


Contact: <thomas.morrison@glasgow.ac.uk>
