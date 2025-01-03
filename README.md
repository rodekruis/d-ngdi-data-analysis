# New Generation Drought Index

## Introduction
Deltares shared two datasets with us to check their potential to use them in an early warning system for droughts. 
The case study is Ethiopia. Some functions were made to quickly inspect the shared datasets

- Project: project name (link to project in CRM)
- Contact: [Björn Bolhuis](https://github.com/BjornBolhuis)

## Input and output
- Inputs:
  - 510 - Anticipatory Action - Next Generation Drought Index\20240914_2300_wflow_sbm_ethiopia_20240823_forecast_seas5_drought_indicators.nc (2 GB)
  - 510 - Anticipatory Action - Next Generation Drought Index\wflow_sbm_ethiopia_20240823_climatology.nc (7.5 GB)
- Output: no file outputs

## Code
- Structure: The functions are found in deltares_drought_functions.py. A Jupyter Notebook is added to quickly read the data and visualize the results.
- Setup/requirements:
  - os
  - pandas
  - numpy
  - geopandas
  - shapely
  - netCDF4
  - matplotlib
  - time

## How to run the code:
The notebook runs top to bottom, no changes are required. Explanations are added in markdown

See example [here](https://github.com/rodekruis/GloFAS-river-depth-analysis)
