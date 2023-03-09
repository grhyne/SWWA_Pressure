# This script generate the reports and the website
library(rmarkdown)
library(readxl)

#Storage: "CB594","CB596","CB598", "CB599","CB601","CB608","CB611",

# Define the track to consider
gdl_list <- c("CB619","CB620","CB621",
              "CB623","CB624","CB627")
# set <- read_excel("data/gdl_settings.xlsx")
# gdl_list = set$gdl_id

# Define the report to produce
report_list <- c("basic_trajectory")

# Generate the following report for each tracks
# This will produce an HTML file for each report
for (report in report_list) {
  dir.create(paste0("reports/", report), showWarnings = F)
  for (gdl in gdl_list) {
    render(paste0("reports/_", report, ".Rmd"),
      output_file = paste0(report, "/", gdl, ".html"),
      params = list(
        gdl_id = gdl
      )
    )
  }
}

# The index page (List of Reports) and all the reports are then combined into a website and moved to
# the docs/ folder to be serve on Github easily.
# Change _site._yml for the overall look of the website and index.Rmd for the content of the homepage

# Update _site._yml with str:
str <- ""
for (report in report_list) {
  str <- paste0(str, '- text: "', report, '"\n  menu:\n')
  for (gdl in gdl_list) {
    str <- paste0(str, '   - text: "', gdl, '"\n')
    str <- paste0(str, '     href: "/SWWA_Pressure/', report, "/", gdl, ".html\n")
  }
}
writeLines(str)

# Render site
render_site("./reports")
