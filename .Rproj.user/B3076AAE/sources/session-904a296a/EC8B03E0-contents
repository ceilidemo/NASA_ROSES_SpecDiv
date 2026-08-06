## 00_DownloadFiles: Download data needed for plant div and RTM##

library(neonUtilities)
library(dplyr)

# DP1.10058.001 = plant presence percent cover
# DP1.10098.001 = vegetation structure
# DP1.10026.001 = foliar traits

setwd("C:/Users/ceilidemarais/Projects/NASA_ROSES_SpecDiv/NEON_PlantDiv")

# Load datasets 
cover_all <- loadByProduct(
  dpID = "DP1.10058.001", 
  site = "all", 
  enddate = "2024-12", 
  package = "basic", 
  release = "current",
  check.size = FALSE,
  include.provisional = TRUE,
  token = "eyJhbGciOiJFUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnLyIsImF1ZCI6Imh0dHBzOi8vZGF0YS5uZW9uc2NpZW5jZS5vcmcvYXBpL3YwLyIsImlhdCI6MTc4MzQ0NDQyOSwiZXhwIjoxODE0OTgwNDI5LCJzdWIiOiJkZW1hcmFpc2NlaWxpQGdtYWlsLmNvbSIsImVtYWlsIjoiZGVtYXJhaXNjZWlsaUBnbWFpbC5jb20iLCJzY29wZSI6InJhdGU6cHVibGljIn0.cPxc8CjWr0Vbpwxj3Li6uN0Thawq7SvN8tyZXnqvGiSlbnOWGR3zG1s6ctt4tadUCGuEV2oLftj0Hqb-60-osw"
)

structure_all <- loadByProduct(
  dpID = "DP1.10098.001",
  site = "all",
  enddate = "2024-12",
  check.size = FALSE,
  include.provisional = TRUE,
  token = "eyJhbGciOiJFUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnLyIsImF1ZCI6Imh0dHBzOi8vZGF0YS5uZW9uc2NpZW5jZS5vcmcvYXBpL3YwLyIsImlhdCI6MTc4MzQ0NDQyOSwiZXhwIjoxODE0OTgwNDI5LCJzdWIiOiJkZW1hcmFpc2NlaWxpQGdtYWlsLmNvbSIsImVtYWlsIjoiZGVtYXJhaXNjZWlsaUBnbWFpbC5jb20iLCJzY29wZSI6InJhdGU6cHVibGljIn0.cPxc8CjWr0Vbpwxj3Li6uN0Thawq7SvN8tyZXnqvGiSlbnOWGR3zG1s6ctt4tadUCGuEV2oLftj0Hqb-60-osw"
)

#save it all
save(structure_all, cover_all, file = "data_in/veg_data_2024.RData")


## Foliar traits

# sites that i have for the other stuff to limit the functional trait data

my_sites <- c("ABBY", "BART", "BLAN", "BONA", "CLBJ", "CPER", "DEJU", "DELA",
              "DSNY", "GRSM", "GUAN", "HARV", "HEAL", "JERC", "JORN", "KONA", 
              "KONZ", "OEAS", "STER", "DCFS", "STER",
              "LAJA", "LENO", "MLBS", "MOAB", "NIWO", "ONAQ", "ORNL", "OSBS",
              "PUUM", "RMNP", "SCBI", "SERC", "SJER", "SOAP", "SRER", "STEI",
              "TALL", "TEAK", "TREE", "UKFS", "UNDE", "WOOD", "WREF", "YELL")

dir.create("C:/tmp", showWarnings = FALSE) # r being annoying with file length

zipsByProduct(
  dpID                = "DP1.10026.001",
  site                = "all",
  package             = "expanded",
  check.size          = FALSE,
  include.provisional = TRUE,
  savepath            = "C:/tmp",
  token = "eyJhbGciOiJFUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnLyIsImF1ZCI6Imh0dHBzOi8vZGF0YS5uZW9uc2NpZW5jZS5vcmcvYXBpL3YwLyIsImlhdCI6MTc4MzQ0NDQyOSwiZXhwIjoxODE0OTgwNDI5LCJzdWIiOiJkZW1hcmFpc2NlaWxpQGdtYWlsLmNvbSIsImVtYWlsIjoiZGVtYXJhaXNjZWlsaUBnbWFpbC5jb20iLCJzY29wZSI6InJhdGU6cHVibGljIn0.cPxc8CjWr0Vbpwxj3Li6uN0Thawq7SvN8tyZXnqvGiSlbnOWGR3zG1s6ctt4tadUCGuEV2oLftj0Hqb-60-osw"
)

neonUtilities::stackByTable(filepath = "C:/tmp/filesToStack10026/", folder = TRUE)

stacked_path <- "C:/tmp/filesToStack10026/stackedFiles/"

foliar_all <- list(
  cfc_fieldData      = read.csv(file.path(stacked_path, "cfc_fieldData.csv")),
  cfc_LMA            = read.csv(file.path(stacked_path, "cfc_LMA.csv")),
  cfc_carbonNitrogen = read.csv(file.path(stacked_path, "cfc_carbonNitrogen.csv")),
  cfc_chlorophyll    = read.csv(file.path(stacked_path, "cfc_chlorophyll.csv")),
  cfc_lignin         = read.csv(file.path(stacked_path, "cfc_lignin.csv")),
  cfc_elements       = read.csv(file.path(stacked_path, "cfc_elements.csv"))
)

save(foliar_all, file = "data_in/foliar_data.RData")

unlink("C:/tmp", recursive = TRUE)
