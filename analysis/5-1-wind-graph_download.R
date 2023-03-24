library(ecmwfr)
library(GeoPressureR)

# Define which track to work with
gdl <- "CB627"

# Load
load(paste0("data/3_static/", gdl, "_static_prob.Rdata"))

# Set credential
Sys.setenv( cds.key="162f79de-0442-42b9-abed-8ad75f09f17d")
Sys.setenv( cds.user="151327")
cds.key <- Sys.getenv("cds.key")
cds.user <- Sys.getenv("cds.user")
wf_set_key(user = cds.user, key = cds.key, service = "cds")

graph_download_wind(pam,
                    area = static_prob,
                    cds_key= cds.key,
                    cds_user= cds.user
)


# Check request at https://cds.climate.copernicus.eu/cdsapp#!/yourrequests
