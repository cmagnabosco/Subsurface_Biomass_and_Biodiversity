# Contents of DiversityAnalyses

This directory contains code and data files to perform the diversity analyses of "The biomass and biodiversity of the continental subsurface". Codes are adapted from https://github.com/LennonLab/ScalingMicroBiodiversity/ to include diversity data from the continental subsurface.

## Code

* **Estimate_S_via_lognormal.ipynb**: Estimates the number of species in the continental subsurface using the lognormal model of biodiversity and various scalings of Nmax. This code is adapted from Fig3.py of https://github.com/LennonLab/ScalingMicroBiodiversity/.
    + **Ncells_Dominance_fit.pickle**
    + **Nreads_Dominance_fit.pickle**
* **Fig3_ABCD_CM.ipynb**: Generates Figure 3 of "The Biomass and Biodiversity of the Continental Subsurface". This code is adapted from Fig1.py of https://github.com/LennonLab/ScalingMicroBiodiversity/.
* **FigS23_predicted_vs_observed_S.ipynb**: Generates Figure S23 of "The Biomass and Biodiversity of the Continental Subsurface".
* **FigS24_Species_Volume.ipynb**: Generates Figure S24 of "The Biomass and Biodiversity of the Continental Subsurface".
* **FigS25_Distance_Decay.ipynb**: Generates Figure S25 of "The Biomass and Biodiversity of the Continental Subsurface".

## Data 
All codes are linked to data stored in the "data" directory of this project. The data directory holds several datasets:
* **data/SubsurfaceDistanceDecay**
* **data/SubsurfaceVolumes**
* **data/nsubsurface**
* **data/subsurface**

**NOTE:** a data/micro directory is needed to run all within this project. Download micro data via the set-up instructions below.

## Set-up
* Codes are linked to a "data/micro" directory. This directory can be dowloaded from  https://github.com/LennonLab/ScalingMicroBiodiversity/tree/master/data/micro and should be added to the "data" directory prior to executing the code found here.

## References
Locey, Kenneth J., and Jay T. Lennon. "Scaling laws predict global microbial diversity." Proceedings of the National Academy of Sciences 113.21 (2016): 5970-5975.

Magnabosco, Cara et al. "The Biomass and Biodiversity of the Continental Subsurface" (2018)
