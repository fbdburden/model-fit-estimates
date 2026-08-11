#' ---
#' title: Diarrhea Etiology Mortality • Data summary
#' author:
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 2
#'     html_preview: true
#' ---

#' # Settings
## required packages ----
library(bd)
library(FERG2)
library(ggplot2)
library(metafor)
library(readxl)
library(rmarkdown)
library(sf)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(knitr)

## global options ----
knitr::opts_chunk$set(fig.width = 10)

#' # Data

## import data ----
all_dta <- list()
all_dta$dta <-
  read_xlsx("ferg2-edtf-diarrheal-highincome-20251208.xlsx", "Mortality") %>%
  filter(!is.na(SOURCE_ID)) 

## data management ----
# harmonize pathogen names
all_dta$dta <- all_dta$dta %>%
  mutate(pathogen = case_when(
    str_detect(OPT_DISEASE, "Campylobacter") ~"CAMP",
    str_detect(OPT_DISEASE, "Crypstosporidium") | str_detect(OPT_DISEASE, "Cryptosporidium") ~ "CRYP",
    str_detect(OPT_DISEASE, "Cyclospora") ~ "CYCL",
    str_detect(OPT_DISEASE, "EAEC") ~ "EAEC",
    str_detect(OPT_DISEASE, "EPEC") ~ "EPEC",
    OPT_DISEASE == "ETEC" | OPT_DISEASE == "ETEC, foodborne" ~ "ETEC",
    OPT_DISEASE == "STEC" | OPT_DISEASE == "EHEC (STEC)" | str_detect(OPT_DISEASE, "VTEC") | str_detect(OPT_DISEASE, "STEC") ~ "STEC",
    str_detect(OPT_DISEASE, "Giardia") ~ "GIAR",
    str_detect(OPT_DISEASE, "Salmonella") | OPT_DISEASE == "S. enterica" ~ "SALM",
    str_detect(OPT_DISEASE, "Norovirus") | str_detect(OPT_DISEASE, "Norwalk") ~ "NORO",
    str_detect(OPT_DISEASE, "Shigella") ~ "SHIG",
    str_detect(OPT_DISEASE, "Vibrio") ~ "VIBR",
    str_detect(OPT_DISEASE, "Rotavirus") ~ "ROTA"
  ))
tbl <- table(all_dta$dta$pathogen, useNA ="always") 
kable(tbl, 
      col.names = c("Pathogen","Number of datapoints"))

#  For STEC some studies need to be summed
STEC <- all_dta$dta %>%
  filter(SOURCE_ID %in% c(3,22,42,49) & OPT_DISEASE %in% c("STEC", "E coli")) %>% 
  mutate(pathogen = "STEC")
i <- c(26:40)  
STEC[, i] <- apply(STEC[, i], 2, function(x) as.numeric(as.character(x)))
col_STEC <- colnames(STEC)
STEC <- aggregate( STEC[,26:40], STEC[,c(1:15,41)], FUN = sum )
col_STEC <- setdiff(col_STEC, colnames(STEC))
for (c in col_STEC){
  STEC[[c]] <- NA
}

all_dta$dta <- all_dta$dta %>%
  filter(!(SOURCE_ID %in% c(3,22,42,49) & OPT_DISEASE %in% c("STEC", "E coli")))
all_dta$dta <- rbind(all_dta$dta, STEC)

# For study 15, 37 STEC only O-157 is reported, this should be removed
all_dta$dta$FLAG <- 0 
all_dta$dta$FLAG <- if_else(all_dta$dta$SOURCE_ID %in% c(15,37) & all_dta$dta$pathogen ==  "STEC",
                            5,
                            all_dta$dta$FLAG)
all_dta$dta$FLAG <- if_else(all_dta$dta$SOURCE_ID %in% c(17,22,37,41) & is.na(all_dta$dta$pathogen),
                            5,
                            all_dta$dta$FLAG)

# remove Salmonella Typhi
all_dta$dta$FLAG <- if_else(all_dta$dta$OPT_SEROTYPE != "Typhi" | is.na(all_dta$dta$OPT_SEROTYPE),
                            all_dta$dta$FLAG,
                            5)
all_dta$dta$FLAG <- if_else(all_dta$dta$OPT_SEROTYPE != "typhi" | is.na(all_dta$dta$OPT_SEROTYPE),
                            all_dta$dta$FLAG,
                            5)

# remove 'Diarrheagenic E. coli other than STEC and ETEC'
all_dta$dta$FLAG[
  grepl("other than STEC and ETEC", all_dta$dta$OPT_DISEASE)] <- 5

# CC: for NLD study was added, only DP from Shigella need to be kept
all_dta$dta$FLAG[
  all_dta$dta$SOURCE_ID == "CC_48" &
    all_dta$dta$OPT_DISEASE != "Shigella"] <-5

Territories <- read_xlsx("Territories_R_20250221.xlsx")
Flag_territory <- unlist(Territories)

all_dta$dta$FLAG_REF_LOCATION <- as.integer(apply(sapply(Flag_territory, function(x) grepl(x, all_dta$dta$REF_LOCATION, ignore.case = TRUE)), 1, any))
all_dta$dta$FLAG_REF_NOTES <- as.integer(apply(sapply(Flag_territory, function(x) grepl(x, all_dta$dta$REF_NOTES, ignore.case = TRUE)), 1, any))
all_dta$dta$FLAG_SOURCE_TITLE <- as.integer(apply(sapply(Flag_territory, function(x) grepl(x, all_dta$dta$SOURCE_TITLE, ignore.case = TRUE)), 1, any))
all_dta$dta$FLAG_TERRITORY <- if_else(all_dta$dta$FLAG_REF_LOCATION + all_dta$dta$FLAG_REF_NOTES + all_dta$dta$FLAG_SOURCE_TITLE >=1 , 1, 0)

all_dta$dta$FLAG <-
  if_else(all_dta$dta$FLAG_TERRITORY == 1 & all_dta$dta$FLAG == 0, 
          1, all_dta$dta$FLAG)

all_dta$dta$COUNTRY <- all_dta$dta$REF_LOCATION_ISO3

all_dta$dta$REF_LOCATION_ISO3 <-
  if_else(is.na(all_dta$dta$REF_LOCATION_ISO3) & all_dta$dta$REF_LOCATION == "New Zealand",
          "NZL", all_dta$dta$REF_LOCATION_ISO3)
sum(is.na(all_dta$dta$REF_LOCATION_ISO3))
all_dta$dta$REF_LOCATION_ISO3 <-
  if_else(is.na(all_dta$dta$REF_LOCATION_ISO3) & all_dta$dta$REF_LOCATION == "Switzerland",
          "CHE", all_dta$dta$REF_LOCATION_ISO3)

all_dta$dta$ISO3 <- all_dta$dta$REF_LOCATION_ISO3
all_dta$dta$ID <- all_dta$dta$SOURCE_ID
all_dta$dta$REG2 <-
  FERG2:::countries$REG2[match(all_dta$dta$ISO3, FERG2:::countries$ISO3)]
all_dta$dta$SUB2 <-
  FERG2:::countries$SUB2[match(all_dta$dta$ISO3, FERG2:::countries$ISO3)]

all_dta$dta <- all_dta$dta %>% 
  mutate(REF_AGE_START = 0, 
         REF_AGE_END = 125, 
         REF_SEX = case_when(is.na(REF_SEX) ~ "All sexes", 
                             TRUE ~ REF_SEX))

all_dta$dta$YEAR <- rowMeans(cbind(all_dta$dta$REF_YEAR_START, all_dta$dta$REF_YEAR_END))

all_dta$dta$REF_YEAR_START2 <- all_dta$dta$REF_YEAR_START
all_dta$dta$REF_YEAR_END2 <- all_dta$dta$REF_YEAR_END

all_dta$dta$REF_YEAR_START <- round(all_dta$dta$YEAR)
all_dta$dta$REF_YEAR_END <- round(all_dta$dta$YEAR)

if (sum(is.na(all_dta$dta$REF_SAMPLE_SIZE))>0) {
  all_dta$dta <- add_pop(all_dta$dta)
} else {
  all_dta$dta$POP <- NA_real_
}

all_dta$dta$REF_YEAR_START <- all_dta$dta$REF_YEAR_START2
all_dta$dta$REF_YEAR_END <- all_dta$dta$REF_YEAR_END2

all_dta$dta$REF_YEAR_START2 <- NULL
all_dta$dta$REF_YEAR_END2 <- NULL

all_dta$dta <- all_dta$dta %>%
  mutate(REF_SAMPLE_SIZE = case_when(is.na(REF_SAMPLE_SIZE) 
                                     & REF_LOC_LEVEL=="National"
                                     ~ POP, 
                                     TRUE ~ as.numeric(REF_SAMPLE_SIZE)))

#  All are cases except study 31
all_dta$dta <- all_dta$dta %>% 
  mutate(VALUE_X = case_when(
    !is.na(VALUE_MEDIAN) & is.na(VALUE_X) ~ VALUE_MEDIAN, # Study 7
    !is.na(VALUE_MEAN) & is.na(VALUE_X) & is.na(VALUE_DENOM) ~ VALUE_MEAN, # Studies 24, 39, 42 & 43
    !is.na(VALUE_MEAN) & is.na(VALUE_X) & !is.na(VALUE_DENOM) ~ REF_SAMPLE_SIZE * VALUE_MEAN / VALUE_DENOM, # Study 31
    .default = VALUE_X))

sum(is.na(all_dta$dta$VALUE_X)) 

all_dta$dta$PERSONYEARS100 <- all_dta$dta$REF_SAMPLE_SIZE / 1e5
all_dta$dta$PERSONYEARS100[all_dta$dta$PERSONYEARS100 == 0] <- NA

## filter data ----

# Studies excluded: from countries which don't belong to WHO 
sum(is.na(all_dta$dta$REG2))
# Studies excluded: before 1990 
length(which(all_dta$dta$YEAR<1990))

data.frame(subset(all_dta$dta, FLAG != 0))

tbl <- table(all_dta$dta$pathogen, useNA ="always") 
kable(tbl, 
      col.names = c("Pathogen","Number of datapoints"))

## estimate mortality ----

pathogens <- c("CAMP","CRYP","CYCL","ETEC","GIAR",
               "NORO","ROTA","SALM","SHIG","STEC","VIBR")
es <- list()

for (p in pathogens) {
  print(p)
  all_dta[[p]] <- all_dta$dta %>% filter(pathogen == p)
  
  ## .. "IRLN" for the log-transformed mortality rate
  es[[p]] <- escalc(xi = VALUE_X, ti = PERSONYEARS100,
                    measure = "IRLN", data = all_dta[[p]])
  es[[p]]$sei <- sqrt(es[[p]]$vi)
  es[[p]]$RAW_MRT <-1e5*es[[p]]$VALUE_X/es[[p]]$REF_SAMPLE_SIZE
  es[["all"]] <- rbind(es[["all"]],es[[p]])
  es[[p]] <- es[[p]] %>% filter(!is.na(yi))
}

all_dta$Exclude <- subset(all_dta$dta, is.na(pathogen))
all_dta$Exclude$yi <- NA
all_dta$Exclude$vi <- NA
all_dta$Exclude$sei <- NA
all_dta$Exclude$RAW_MRT <- NA

es$all <- rbind(es$all, all_dta$Exclude) # es$all should be used to list studies for countries

sum(is.na(es$all$yi))
es$all$FLAG <-
  if_else(is.na(es$all$yi) & es$all$FLAG == 0,
          4, es$all$FLAG)
es[["all_yi"]] <- es$all %>% filter(!is.na(yi))

tbl <- table(subset(es$all_yi, FLAG == 0)$pathogen, useNA ="always") 
kable(tbl, 
      col.names = c("Pathogen","Number of datapoints"))

#' # Retained data points
cols <- c("SOURCE_ID", "SOURCE_AUTHOR", "SOURCE_YEAR", "REF_LOCATION_ISO3",
          "OPT_DISEASE", "OPT_SEROTYPE", "RAW_MRT", "yi", "sei")

cat(sum(es[["CAMP"]]$FLAG == 0), "data points\n")
knitr::kable(es[["CAMP"]] |> subset(FLAG == 0, select = cols), row.names = FALSE, digits = 3)

cat(sum(es[["CRYP"]]$FLAG == 0), "data points\n")
knitr::kable(es[["CRYP"]] |> subset(FLAG == 0, select = cols), row.names = FALSE, digits = 3)

cat(sum(es[["CYCL"]]$FLAG == 0), "data points\n")
knitr::kable(es[["CYCL"]] |> subset(FLAG == 0, select = cols), row.names = FALSE, digits = 3)

cat(sum(es[["ETEC"]]$FLAG == 0), "data points\n")
knitr::kable(es[["ETEC"]] |> subset(FLAG == 0, select = cols), row.names = FALSE, digits = 3)

cat(sum(es[["GIAR"]]$FLAG == 0), "data points\n")
knitr::kable(es[["GIAR"]] |> subset(FLAG == 0, select = cols), row.names = FALSE, digits = 3)

cat(sum(es[["NORO"]]$FLAG == 0), "data points\n")
knitr::kable(es[["NORO"]] |> subset(FLAG == 0, select = cols), row.names = FALSE, digits = 3)

cat(sum(es[["ROTA"]]$FLAG == 0), "data points\n")
knitr::kable(es[["ROTA"]] |> subset(FLAG == 0, select = cols), row.names = FALSE, digits = 3)

cat(sum(es[["SALM"]]$FLAG == 0), "data points\n")
knitr::kable(es[["SALM"]] |> subset(FLAG == 0, select = cols), row.names = FALSE, digits = 3)

cat(sum(es[["SHIG"]]$FLAG == 0), "data points\n")
knitr::kable(es[["SHIG"]] |> subset(FLAG == 0, select = cols), row.names = FALSE, digits = 3)

cat(sum(es[["STEC"]]$FLAG == 0), "data points\n")
knitr::kable(es[["STEC"]] |> subset(FLAG == 0, select = cols), row.names = FALSE, digits = 3)

cat(sum(es[["VIBR"]]$FLAG == 0), "data points\n")
knitr::kable(es[["VIBR"]] |> subset(FLAG == 0, select = cols), row.names = FALSE, digits = 3)

#' # Graphical representations 
#' ## Campylobacter: data availability
p <- "CAMP"

#Data availability 
#+ fig.height=3
plot_data(es[[p]], by = "REG2")
#+ fig.height=4
plot_data(es[[p]], by = "SUB2")
#+ fig.height=4
plot_world_imputation(es[[p]], sub = "SUB2")
#+ fig.height=4
plot_world_data(es[[p]])
#+ fig.height=4
ggplot(es[[p]], aes(x = REG2, y = RAW_MRT, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw mortality per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_MRT)),
                     breaks = pretty(es[[p]]$RAW_MRT, n = 20))
#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_MRT, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported mortality per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$REG2), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("WHO region", "Mean", "Lower", "Upper"))

#' ## Cryptosporidium: data availability
p <- "CRYP"

#Data availability 
#+ fig.height=3
plot_data(es[[p]], by = "REG2")
#+ fig.height=4
plot_data(es[[p]], by = "SUB2")
#+ fig.height=4
plot_world_imputation(es[[p]], sub = "SUB2")
#+ fig.height=4
plot_world_data(es[[p]])
#+ fig.height=4
ggplot(es[[p]], aes(x = REG2, y = RAW_MRT, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw mortality per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_MRT)),
                     breaks = pretty(es[[p]]$RAW_MRT, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_MRT, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported mortality per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$REG2), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("WHO region", "Mean", "Lower", "Upper"))

#' ## Cyclospora : data availability and quality 
p <- "CYCL"

#Data availability 
#+ fig.height=3
plot_data(es[[p]], by = "REG2")
#+ fig.height=4
plot_data(es[[p]], by = "SUB2")
#+ fig.height=4
plot_world_imputation(es[[p]], sub = "SUB2")
#+ fig.height=4
plot_world_data(es[[p]])
#+ fig.height=4
ggplot(es[[p]], aes(x = REG2, y = RAW_MRT, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw mortality per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_MRT)),
                     breaks = pretty(es[[p]]$RAW_MRT, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_MRT, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported mortality per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$REG2), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("WHO region", "Mean", "Lower", "Upper"))


#' ## Enteroaggregative E. coli: data availability and quality
p <- "EAEC"

#Data availability 
print("No data available for Enteroaggregative E. coli (EAEC)")

#' ## Entamoeba histolytica: data availability and quality
p <- "ENTA"

#Data availability 
print("No data available for Entamoeba histolytica (ENTA)")


#'## Enteropathogenic E. coli: data availability and quality 
p <- "EPEC"

#Data availability 
print("No data available for Enteropathogenic E. coli (EPEC)")

#' ## Enterotoxigenic E. coli : data availability and quality
p <- "ETEC"

#Data availability 
#+ fig.height=3
plot_data(es[[p]], by = "REG2")
#+ fig.height=4
plot_data(es[[p]], by = "SUB2")
#+ fig.height=4
plot_world_imputation(es[[p]], sub = "SUB2")
#+ fig.height=4
plot_world_data(es[[p]])
#+ fig.height=4
ggplot(es[[p]], aes(x = REG2, y = RAW_MRT, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw mortality per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_MRT)),
                     breaks = pretty(es[[p]]$RAW_MRT, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_MRT, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported mortality per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$REG2), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("WHO region", "Mean", "Lower", "Upper"))


#'## Giardia: data availability and quality 
p <- "GIAR"

#Data availability 
#+ fig.height=3
plot_data(es[[p]], by = "REG2")
#+ fig.height=4
plot_data(es[[p]], by = "SUB2")
#+ fig.height=4
plot_world_imputation(es[[p]], sub = "SUB2")
#+ fig.height=4
plot_world_data(es[[p]])
#+ fig.height=4
ggplot(es[[p]], aes(x = REG2, y = RAW_MRT, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw mortality per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_MRT)),
                     breaks = pretty(es[[p]]$RAW_MRT, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_MRT, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported mortality per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$REG2), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("WHO region", "Mean", "Lower", "Upper"))


#'## Norovirus: data availability and quality
p <- "NORO"

#Data availability 
#+ fig.height=3
plot_data(es[[p]], by = "REG2")
#+ fig.height=4
plot_data(es[[p]], by = "SUB2")
#+ fig.height=4
plot_world_imputation(es[[p]], sub = "SUB2")
#+ fig.height=4
plot_world_data(es[[p]])
#+ fig.height=4
ggplot(es[[p]], aes(x = REG2, y = RAW_MRT, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw mortality per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_MRT)),
                     breaks = pretty(es[[p]]$RAW_MRT, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_MRT, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported mortality per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$REG2), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("WHO region", "Mean", "Lower", "Upper"))


#'## Rotavirus: data availability and quality
p <- "ROTA"

#Data availability 
#+ fig.height=3
plot_data(es[[p]], by = "REG2")
#+ fig.height=4
plot_data(es[[p]], by = "SUB2")
#+ fig.height=4
plot_world_imputation(es[[p]], sub = "SUB2")
#+ fig.height=4
plot_world_data(es[[p]])
#+ fig.height=4
ggplot(es[[p]], aes(x = REG2, y = RAW_MRT, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw mortality per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_MRT)),
                     breaks = pretty(es[[p]]$RAW_MRT, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_MRT, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported mortality per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$REG2), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("WHO region", "Mean", "Lower", "Upper"))


#'## Salmonella: data availability and qualit
p <- "SALM"

#Data availability 
#+ fig.height=3
plot_data(es[[p]], by = "REG2")
#+ fig.height=4
plot_data(es[[p]], by = "SUB2")
#+ fig.height=4
plot_world_imputation(es[[p]], sub = "SUB2")
#+ fig.height=4
plot_world_data(es[[p]])
#+ fig.height=4
ggplot(es[[p]], aes(x = REG2, y = RAW_MRT, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw mortality per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_MRT)),
                     breaks = pretty(es[[p]]$RAW_MRT, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_MRT, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported mortality per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$REG2), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("WHO region", "Mean", "Lower", "Upper"))


#'## Shigella: data availability and quality
p <- "SHIG"

#Data availability 
#+ fig.height=3
plot_data(es[[p]], by = "REG2")
#+ fig.height=4
plot_data(es[[p]], by = "SUB2")
#+ fig.height=4
plot_world_imputation(es[[p]], sub = "SUB2")
#+ fig.height=4
plot_world_data(es[[p]])
#+ fig.height=4
ggplot(es[[p]], aes(x = REG2, y = RAW_MRT, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw mortality per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_MRT)),
                     breaks = pretty(es[[p]]$RAW_MRT, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_MRT, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported mortality per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$REG2), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("WHO region", "Mean", "Lower", "Upper"))


#'## Shiga toxin-producing E. coli: data availability and quality
p <- "STEC"

#Data availability 
#+ fig.height=3
plot_data(es[[p]], by = "REG2")
#+ fig.height=4
plot_data(es[[p]], by = "SUB2")
#+ fig.height=4
plot_world_imputation(es[[p]], sub = "SUB2")
#+ fig.height=4
plot_world_data(es[[p]])
#+ fig.height=4
ggplot(es[[p]], aes(x = REG2, y = RAW_MRT, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw mortality per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_MRT)),
                     breaks = pretty(es[[p]]$RAW_MRT, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_MRT, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported mortality per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$REG2), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("WHO region", "Mean", "Lower", "Upper"))


#'## Vibrio cholerae:  data availability & quality
p <- "VIBR"

#Data availability 
#+ fig.height=3
plot_data(es[[p]], by = "REG2")
#+ fig.height=4
plot_data(es[[p]], by = "SUB2")
#+ fig.height=4
plot_world_imputation(es[[p]], sub = "SUB2")
#+ fig.height=4
plot_world_data(es[[p]])
#+ fig.height=4
ggplot(es[[p]], aes(x = REG2, y = RAW_MRT, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw mortality per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_MRT)),
                     breaks = pretty(es[[p]]$RAW_MRT, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_MRT, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported mortality per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_MRT, list(es[[p]]$REG2), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("WHO region", "Mean", "Lower", "Upper"))

#' # Session info
Date <- format(Sys.Date(), "%Y%m%d")
saveRDS(es, paste0("es_HIC_MRT_", Date, ".rds"))
sessioninfo::session_info()

##rmarkdown::render("01-data-hic-mrt-v20251208.R")
