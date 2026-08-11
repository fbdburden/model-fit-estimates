### FERG2/EDTF/DIARRHEA
### 11/04/2025

## load incidence files
diarrhea_inc <- readRDS("../diarrhea-inc-ihme.rds")
diarrhea_mrt <- readRDS("../diarrhea-mrt-who.rds")

## check order of years/countries
all(names(diarrhea_inc) == names(diarrhea_mrt))
all(names(diarrhea_inc[[1]]) == names(diarrhea_mrt[[1]]))
all(sapply(diarrhea_inc[[1]][[1]], function(x) x$AGE) ==
      sapply(diarrhea_mrt[[1]][[1]], function(x) x$AGE))

# str(diarrhea_mrt$`2021`$BEL)
sum(sapply(diarrhea_mrt$`2021`$BEL, function(x) x$POP))
mean(rowSums(sapply(diarrhea_mrt$`2021`$BEL, function(x) x$MRT_NR)))
sum(sapply(diarrhea_inc$`2021`$BEL, function(x) x$POP))
mean(rowSums(sapply(diarrhea_inc$`2021`$BEL, function(x) x$INC_NR)))

## add YLL
for (year in as.character(all_yrs)) {
  for (country in names(diarrhea_mrt[[year]])) {
    for (agesex in seq_along(diarrhea_mrt[[year]][[country]])) {
      diarrhea_mrt[[year]][[country]][[agesex]]$YLL_NR <-
        diarrhea_mrt[[year]][[country]][[agesex]]$MRT_NR *
        dalymod:::get_rle(
          rle = rle,
          age = diarrhea_mrt[[year]][[country]][[agesex]]$AGE,
          sex = diarrhea_mrt[[year]][[country]][[agesex]]$SEX)
    }
  }
}
# str(diarrhea_mrt[[1]][[1]])

## define country groupings
countries_lmic <-
  c("AFG", "AGO", "ALB", "ARG", "ARM", "AZE", "BDI", "BEN", "BFA", 
    "BGD", "BGR", "BIH", "BLR", "BLZ", "BOL", "BRA", "BTN", "BWA", 
    "CAF", "CHN", "CIV", "CMR", "COD", "COG", "COL", "COM", "CPV", 
    "CRI", "CUB", "DJI", "DMA", "DOM", "DZA", "ECU", "EGY", "ERI", 
    "ETH", "FJI", "FSM", "GAB", "GEO", "GHA", "GIN", "GMB", "GNB", 
    "GNQ", "GRD", "GTM", "HND", "HTI", "IDN", "IND", "IRN", "IRQ", 
    "JAM", "JOR", "KAZ", "KEN", "KGZ", "KHM", "KIR", "LAO", "LBN", 
    "LBR", "LBY", "LCA", "LKA", "LSO", "MAR", "MDA", "MDG", "MDV", 
    "MEX", "MHL", "MKD", "MLI", "MMR", "MNE", "MNG", "MOZ", "MRT", 
    "MUS", "MWI", "MYS", "NAM", "NER", "NGA", "NIC", "NPL", "PAK", 
    "PER", "PHL", "PLW", "PNG", "PRK", "PRY", "RUS", "RWA", "SDN", 
    "SEN", "SLB", "SLE", "SLV", "SOM", "SRB", "SSD", "STP", "SUR", 
    "SWZ", "SYC", "SYR", "TCD", "TGO", "THA", "TJK", "TKM", "TLS", 
    "TON", "TUN", "TUR", "TUV", "TZA", "UGA", "UKR", "UZB", "VCT", 
    "VEN", "VNM", "VUT", "WSM", "YEM", "ZAF", "ZMB", "ZWE")
countries_hic <-
  c("AND", "ARE", "ATG", "AUS", "AUT", "BEL", "BHR", "BHS", "BRB", 
    "BRN", "CAN", "CHE", "CHL", "COK", "CYP", "CZE", "DEU", "DNK", 
    "ESP", "EST", "FIN", "FRA", "GBR", "GRC", "GUY", "HRV", "HUN", 
    "IRL", "ISL", "ISR", "ITA", "JPN", "KNA", "KOR", "KWT", "LTU", 
    "LUX", "LVA", "MCO", "MLT", "NIU", "NLD", "NOR", "NRU", "NZL", 
    "OMN", "PAN", "POL", "PRT", "QAT", "ROU", "SAU", "SGP", "SMR", 
    "SVK", "SVN", "SWE", "TTO", "URY", "USA")
