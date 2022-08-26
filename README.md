# Swainson's Warbler Pressure Analysis

This analysis is part of a multi-state collaboration which deployed multi-sensor geolocators on Swainson's Warblers (Limnothlypis swainsonii) in 2021 and recovered tags in 2022. 17 geolocators recorded light and atmospheric pressure. The [`GeoPressureR`](https://raphaelnussbaumer.com/GeoPressureR/) package will be used to refine geolocation probabilities previously analyzed by GeoLight and SGAT.

Here we are implementing the [`GeoPressureTemplate`](https://github.com/Rafnuss/GeoPressureTemplate). Therefore, files and analysis are built upon the original structure for reproducible purposes.

##Template Structure:
```
GeoPressureTemplate
├── DESCRIPTION          		                # project metadata and dependencies
├── README.md            		                # top-level description of content and guide to users
├── GeoPressureTemplate.Rproj               # R project file
├── data                                    # Folder structured by order of use
│   ├── 0_PAM                               # Folder with raw geolocator data grouped by gdl_id
│   │   ├── 18LX
│   │   │   ├── 18LX_20180725.acceleration
│   │   │   ├── 18LX_20180725.glf
│   │   │   ├── 18LX_20180725.pressure 
│   │   │   └── ...
│   │   └── 22BT
│   │       └── ...
│   ├── 1_pressure                          # Data generated with analyis/1-pressure.R
│   │   ├── 18LX_pressure_prob.Rdata
│   │   └── labels
│   │       ├── 18LX_act_pres-labeled.csv
│   │       ├── 18LX_act_pres.csv
│   │       └── ...                    
│   ├── 2_light                             # Data generated with analyis/2-light.R
│   │   ├── 18LX_light_prob.Rdata
│   │   └── labels
│   │       ├── 18LX_light-labeled.csv
│   │       ├── 18LX_light.csv
│   │       └── ...    
│   ├── 3_static                            # Data generated with analyis/3-static.R
│   │   ├── 18LX_static_prob.Rdata
│   │   └── ...
│   ├── 4_basic_graph                       # Data generated with analyis/3-basic_graph.R
│   │   ├── 18LX_basic_graph.Rdata
│   │   └── ...
│   ├── 5_wind_graph
│   │   └── ERA5_wind
│   │       ├──
│   │       └── ...
│   └── gpr_settings.xlsx
├── analysis                                # R code used to analyse your data. Follow the order
│   ├── 1-pressure.R
│   ├── 2-light.R
│   ├── 3-static.R
│   ├── 4-basic-graph.R
│   ├── 5-1-wind-graph_request.R
│   ├── 5-2-wind-graph_transfer.R
│   ├── 5-3-wind-graph_create.R
│   ├── 5-4-wind-graph_analyse.R
│   └── 99-combined.R
└── reports                                 # Generate HTML report to be shared (see below for details)
│   ├── _basic_trajectory.Rmd
│   ├── _site.yml
│   ├── _technical_details.Rmd
│   ├── basic_trajectory
│   │   └── 18LX.html
│   ├── technical_details
│   │   └── 18LX.html
│   ├── index.Rmd
│   └── make_reports.R
└── docs                                      # Folder where your reports will be served as a website on Github Page
    └── ...
```
</details>


Now that you are set-up, it's time to start the serious work. :grimacing: Follow the order of the `.R` code in the `analysis/` folder. They follow the same order as the vignettes (but with different numerotation).

|  GeoPressureTemplate analysis |  GeoPressureR vignettes  |
|---|---|
|  `1-pressure.R`  |  [Creating probability maps from pressure data](https://raphaelnussbaumer.com/GeoPressureR/articles/pressure-map.html) |
|  `2-light.R` |  [Creating probability maps from light data](https://raphaelnussbaumer.com/GeoPressureR/articles/light-map.html) |
|  `3-static.R` | [Preparing data for trajectory modelling](https://raphaelnussbaumer.com/GeoPressureR/articles/preparing-data.html)  |
|  `4-basic-graph.R` |  [Modeling trajectory with a graph](https://raphaelnussbaumer.com/GeoPressureR/articles/basic-graph.html) |
|  `5-1-wind-graph_request.R` |  [Improving the graph with wind - Request wind data on ERA5](https://raphaelnussbaumer.com/GeoPressureR/articles/wind-graph.html#download-wind-data) |
|  `5-2-wind-graph_transfer.R` |  [Improving the graph with wind - Download wind data on ERA5](https://raphaelnussbaumer.com/GeoPressureR/articles/wind-graph.html#download-wind-data) |
|  `5-3-wind-graph_create.R` |  [Improving the graph with wind - Create](https://raphaelnussbaumer.com/GeoPressureR/articles/wind-graph.html#add-wind-to-graph) |
|  `5-4-wind-graph_analyse.R` |  [Improving the graph with wind - Outputs](https://raphaelnussbaumer.com/GeoPressureR/articles/wind-graph.html#output-1-shortest-path-with-wind) |
|  `99-combined.R` |  Run all steps for multiple tracks. |


