# Drivers of site fidelity in ungulates

Morrison TA, Merkle J, Hopcraft JGC, Aikens EO, Beck JL, Boone RB, Courtemanch A, Dwinnell SP, Fairbanks WS, Griffith B, Middleton AD, Monteith KL, Oates B, Riotte-Lambert L, Sawyer H, Smith K, Stabach JA, Taylor KL, Kauffman MJ. 2021. Drivers of site fidelity in ungulates. Journal of Animal Ecology, 90, 4, 955-966.
https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.13425

There are 2 main 'run' scripts: 
(1) sitefidelity.R calculates inter-year site fidelity from raw telemetry data and compiles all inidivduals into an aggregated dataframe, which is saved as ‘fidelity_toshare.csv’. It also calculates a number of covariates from MODIS ndvi data and from the spatial attributes of the telemetry data. Note this script takes quite a while to run. 
(2) analysis_toShare.R runs the statistical analyses and plotting of data that are included in the paper, using the aggregated dataset  ‘fidelity_toshare.csv’. Each row of this dataset is a unique individual, with site fidelity summarized across year1 and year2. 

If interested in the actual calculation of site fidelity, it is written as a function in the [functions/iyd_toShare.R][iyd] file. 

If interested in predictability metrics, they are written as function in [functions/CTime.R][ctime] and [functions/Cspace.R][cspace] files. 

Contact: <thomas.morrison@glasgow.ac.uk>

[iyd]: https://github.com/morrisontom/sitefidelity/blob/91b649241475bf0204b8ea879f19458ec199d475/functions/iyd_toshare.R
[ctime]: https://github.com/morrisontom/sitefidelity/blob/91b649241475bf0204b8ea879f19458ec199d475/functions/Ctime.R
[cspace]: https://github.com/morrisontom/sitefidelity/blob/80f206c2d66ffd2aa7fad3eebd6d55611bf9be58/functions/Cspace.R
