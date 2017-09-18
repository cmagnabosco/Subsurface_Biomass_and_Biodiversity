# Subsurface_Biomass_and_Biodiversity
This repository contains data and code used in the preparation of "The biomass and biodiversity of the continental subsurface"

# Data File Description
- **1000_indices_for_bootstrap.csv**: indices for bootstrap analysis
- **cores_with_PCR.csv**: subsurface cell count data derived from core samples
- **metadata_by_grid.csv**: metadata (columns) associated with the 4303 map-pixels (rows) used for the integration and visualiation of continental subsurface biomass
- **subsurface_database_used_in_shiny.csv**: dataset associated with interactive visualization https://caramagnabosco.shinyapps.io/subsufacebiomass/

# Code File Description
- C**rust_Specific_Fits.R**: code to estimate the subsurface biomass based on crust-specific power fits with depth
- **Depth_and_Temperature_Fits.R**: code to estimate the subsurface biomass based a power fits with depth and logarithmic fit with temperature
- **Depth_and_Temperature_GLM.R**: code to estimate the subsurface biomass based on a generalized linear model using depth and temperature
