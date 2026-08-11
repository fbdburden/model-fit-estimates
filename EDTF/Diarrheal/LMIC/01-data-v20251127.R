#' ---
#' title: Diarrhea Etiology • Data summary • v20250819
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
library(kableExtra)

## global options ----
knitr::opts_chunk$set(fig.width = 10)

#' # Data

## import data ----
all_dta <- list()
check <- list()
sources <- read_xlsx("../2025-11-27_FERG_results_UVA-SOM.xlsx", "SOURCE_index")
sites <- read_xlsx("../2025-11-27_FERG_results_UVA-SOM.xlsx", "SITE_index")
pathogens <-
  c("CAMP", "CRYP", "CYCL", "EAEC", "ENTA", "EPTP", "EPAP", "ETLT", "ETST",
    "GIAR", "NORO", "ROTA", "SALM", "SHIG", "STEC", "VIBR")

zero_cases <-
  read_xlsx("endemic_countries.xlsx") %>%
  select(ISO3, Country, edtf_diarrheal) %>% 
  rename(COUNTRY = ISO3, COUNTRY_LABEL = Country) 

Territories <-
  read_xlsx("Territories_R_20250221.xlsx")
Flag_territory <- unlist(Territories)

for (i in seq_along(pathogens)) {
  
  dta <- read_xlsx("../2025-11-27_FERG_results_UVA-SOM.xlsx", pathogens[i]) %>%
         left_join(sites) %>% left_join(sources) %>%
         unique()
  dta <- dta[rowSums(is.na(dta)) < ncol(dta) - 4, ]
  
  dta <- dta[, c(25,1,2,42,44,45,47,48,35,36,28,27,30,29,5,9,10,11,4,12:16)]
  
  cols <- c("SOURCE_ID","SITE_ID","EST_ID","SOURCE_AUTHOR","SOURCE_YEAR","SOURCE_TITLE", "SOURCE_DOI", "SOURCE-DESIGN",
            "REF_DATE_START", "REF_DATE_END","REF_LOC_LEVEL","COUNTRY_INCOME", 
            "REF_LOCATION", "ISO3", "AGE_GROUP", "REF_AGE_START", "REF_AGE_END","REF_SEX",  
            "SYNDROM", "DX","STRAIN","SUBJECTS","REF_SAMPLE_SIZE","VALUE_X")
  colnames(dta) <- cols
  
  class(dta) <- "data.frame"
  str(dta)
  
## data management ----

dta$COUNTRY <- dta$ISO3
dta$ID <- dta$SOURCE_ID
dta$REG2 <-
  FERG2:::countries$REG2[match(dta$ISO3, FERG2:::countries$ISO3)]
dta$SUB2 <-
  FERG2:::countries$SUB2[match(dta$ISO3, FERG2:::countries$ISO3)]

dta<-dta %>% 
     mutate( REF_SEX = case_when(REF_SEX == "Both" ~ "All sexes", 
                                 TRUE ~ REF_SEX),
            pathogen = paste(pathogens[i]),
            month_start = month(REF_DATE_START),
            year_start = year(REF_DATE_START), 
            month_end = month(REF_DATE_END),
            year_end = year(REF_DATE_END), 
            REF_YEAR_START = case_when( 
                                        (month_start == month_end |
                                         month_start == month_end + 1) 
                                         & (year_start == year_end |
                                            year_end == year_start + 1) 
                                        & month_start>6 ~ year_end,
                                        TRUE ~ year_start),
                REF_YEAR_END=case_when( 
                                        (month_start == month_end | 
                                         month_start == month_end + 1) 
                                         & (year_start== year_end | 
                                            year_end == year_start + 1) 
                                         & month_start<=6 ~ year_start,
                                        TRUE ~ year_end))
                
dta$YEAR <- rowMeans(cbind(dta$REF_YEAR_START, dta$REF_YEAR_END))

dta <- dta %>% mutate(REF_AGE_START = case_when(is.na(REF_AGE_START) ~ 0, 
                                                TRUE ~ REF_AGE_START), 
                      REF_AGE_END = case_when(is.na(REF_AGE_END) ~ 125, 
                                              TRUE ~ REF_AGE_END), 
                      REF_SEX = case_when(is.na(REF_SEX) ~ "All sexes", 
                                          TRUE ~ REF_SEX))

if (sum(is.na(dta$REF_SAMPLE_SIZE))>0) {
  dta <- add_pop(dta)
  } else {dta$POP <- NA_real_
}

dta <- dta %>%
       mutate(REF_SAMPLE_SIZE = case_when(is.na(REF_SAMPLE_SIZE) 
                                                 & REF_LOC_LEVEL=="National"
                                                ~ POP, 
                                                TRUE ~ REF_SAMPLE_SIZE))
# EXCLUDE ASYMPTOMATIC
dta <- dta %>%
  mutate(FLAG = case_when(
    SYNDROM == "Asymptomatic" ~ 7,
    .default = 0))

dta$FLAG_REF_LOCATION <- as.integer(apply(sapply(Flag_territory, function(x) grepl(x, dta$REF_LOCATION, ignore.case = TRUE)), 1, any))
dta$FLAG_SOURCE_TITLE <- as.integer(apply(sapply(Flag_territory, function(x) grepl(x, dta$SOURCE_TITLE, ignore.case = TRUE)), 1, any))
dta$FLAG_TERRITORY <- if_else(dta$FLAG_REF_LOCATION + dta$FLAG_SOURCE_TITLE >=1 , 1, 0)

dta$FLAG <- if_else(dta$FLAG_TERRITORY == 1 & dta$FLAG == 0,
                    1, 
                    dta$FLAG)

## estimate incidence

## .. "PLO" for the logit transformed proportion 

es <- escalc(xi = VALUE_X, ni = REF_SAMPLE_SIZE,
             measure = "PLO", data = dta)
es$sei <- sqrt(es$vi)

es <- es[!is.na(es$REG2) & es$YEAR>=1990, ]
es <- es %>% filter(!is.na(YEAR)) # filter out studies without year
# Exclude studies from HIC countries
es <- merge(es, zero_cases[,c("COUNTRY", "edtf_diarrheal")], by.x="ISO3", by.y="COUNTRY")
es <- es %>% 
  mutate(FLAG = case_when(
    FLAG == 0 & edtf_diarrheal == 0 ~ 5,
    .default = FLAG))
es$edtf_diarrheal <- NULL

## Classification of age and syndrom variables
es <- es %>%
  mutate(AGE = case_when(
    AGE_GROUP == "Adolescents" | AGE_GROUP == "Adults" | AGE_GROUP == "School age children" ~ 1, 
    AGE_GROUP == "Pre-school age children" ~ 2,
    AGE_GROUP == "Combined ages" ~ 3),
    SYNDROMTYPE = case_when(
      SYNDROM == "Medically attended diarrhea - inpatient" | SYNDROM == "Mortality" ~ 1, 
      SYNDROM == "Medically attended diarrhea - outpatient" | SYNDROM == "Community detected diarrhea" ~ 2,
      SYNDROM == "Asymptomatic" ~ 3),
    DIAGNOSIS = case_when(
      DX == "Culture" ~ 1,
      DX == "Immunoassay" ~ 3,
      DX == "Microscopy" ~ 4,
      DX == "Molecular" ~ 2,
      DX == "Other" ~ 5,
      DX == "Unclear/unspecified" ~ 5,
      DX == "Unclear/Unspecified" ~ 5),
    REFERENCE = case_when(
      DIAGNOSIS == 2 ~ 1, 
      .default = 2))

es$AGE <- factor(es$AGE, c(1,2,3), c("Age above or equal 5", "Age below 5", "Mixed ages"))
es$SYNDROMTYPE <- factor(es$SYNDROMTYPE, c(3,1,2), c("Asymptomatic", "Inpatient", "Outpatient"))
es$DIAGNOSIS <- factor(es$DIAGNOSIS, c(1,2,3,4,5), c("Culture", "Molecular", "Immunoassay", "Microscopy","Other/UNK"))
es$REFERENCE <- factor(es$REFERENCE, c(1,2), c("Reference", "Other"))

es <- es %>% 
  mutate(FLAG = case_when(
    FLAG == 0 & is.na(REG2) ~ 2, 
    FLAG == 0 & YEAR < 1990 ~ 3, 
    FLAG == 0 & is.na(yi) ~ 4,
    .default=FLAG))

es$RAW_PROP_100000 <- es$VALUE_X/es$REF_SAMPLE_SIZE

all_dta[[i]] <- es
check[[i]] <- dta

}

## filter data ----

check <- bind_rows(check) 
# Studies excluded: from countries which don't belong to WHO 
sum(is.na(check$REG2))
id <- is.na(check$REG2)
check[id, c(29,1,2,4,5,13) ]
# Studies excluded: before 1990 
length(which(check$YEAR<1990))
id <- (check$YEAR<1990)
check[id, c(29,1,4,5,13) ]
sum(is.na(es$yi))
# Exclude studies from HIC countries
check <- merge(check, zero_cases[,c("COUNTRY", "edtf_diarrheal")], by.x="ISO3", by.y="COUNTRY")
addmargins(table(check$pathogen, check$edtf_diarrheal), 2)
sum(check$edtf_diarrheal)
id <- (check$edtf_diarrheal == 0)
check[id, c(29,1,4,5,13) ]

#' # Graphical representations 

#' ## Campylobacter: data availability and quality

#Data availability 
i <- 1
#+ fig.height=3
plot_data(subset(all_dta[[i]], FLAG == 0), by = "REG2")
#+ fig.height=4
plot_data(subset(all_dta[[i]], FLAG == 0), by = "SUB2")
#+ fig.height=4
plot_world_imputation(subset(all_dta[[i]], FLAG == 0), sub = "SUB2")
#+ fig.height=4
plot_world_data(subset(all_dta[[i]], FLAG == 0))

#Data quality

#+ warning=FALSE
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi, group = REG2)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "loess") +
  facet_wrap(~REG2) +
  theme_bw()

#+ fig.height=8
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, aes(col = REG2)) +
  facet_wrap(SYNDROMTYPE ~ AGE) +
  theme_bw()

ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = REG2, y = RAW_PROP_100000, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw proportion (VALUE_X/SAMPLE_SIZE) reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000)),
                     breaks = pretty(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000, n = 20))

#' ## Cryptosporidium: data availability and quality
i <- 2
# Data availability
#+ fig.height=3
plot_data(subset(all_dta[[i]], FLAG == 0), by = "REG2")
#+ fig.height=4
plot_data(subset(all_dta[[i]], FLAG == 0), by = "SUB2")
#+ fig.height=4
plot_world_imputation(subset(all_dta[[i]], FLAG == 0), sub = "SUB2")
#+ fig.height=4
plot_world_data(subset(all_dta[[i]], FLAG == 0))

# Data quality

#+ warning=FALSE
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi, group = REG2)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "loess") +
  facet_wrap(~REG2) +
  theme_bw()

#+ fig.height=8
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, aes(col = REG2)) +
  facet_wrap(SYNDROMTYPE ~ AGE) +
  theme_bw()

ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = REG2, y = RAW_PROP_100000, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw proportion (VALUE_X/SAMPLE_SIZE) reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000)),
                     breaks = pretty(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000, n = 20))

#' ## Cyclospora : data availability and quality 
i <- 3
# Data availability
#+ fig.height=3
plot_data(subset(all_dta[[i]], FLAG == 0), by = "REG2")
#+ fig.height=4
plot_data(subset(all_dta[[i]], FLAG == 0), by = "SUB2")
#+ fig.height=4
plot_world_imputation(subset(all_dta[[i]], FLAG == 0), sub = "SUB2")
#+ fig.height=4
plot_world_data(subset(all_dta[[i]], FLAG == 0))

# Data quality

#+ warning=FALSE
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi, group = REG2)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "loess") +
  facet_wrap(~REG2) +
  theme_bw()

#+ fig.height=8
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, aes(col = REG2)) +
  facet_wrap(SYNDROMTYPE ~ AGE) +
  theme_bw()

ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = REG2, y = RAW_PROP_100000, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw proportion (VALUE_X/SAMPLE_SIZE) reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000)),
                     breaks = pretty(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000, n = 20))

#' ## Enteroaggregative E. coli: data availability and quality
i <- 4
# Data availability
#+ fig.height=3
plot_data(subset(all_dta[[i]], FLAG == 0), by = "REG2")
#+ fig.height=4
plot_data(subset(all_dta[[i]], FLAG == 0), by = "SUB2")
#+ fig.height=4
plot_world_imputation(subset(all_dta[[i]], FLAG == 0), sub = "SUB2")
#+ fig.height=4
plot_world_data(subset(all_dta[[i]], FLAG == 0))

# Data quality

#+ warning=FALSE
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi, group = REG2)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "loess") +
  facet_wrap(~REG2) +
  theme_bw()

#+ fig.height=8
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, aes(col = REG2)) +
  facet_wrap(SYNDROMTYPE ~ AGE) +
  theme_bw()

ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = REG2, y = RAW_PROP_100000, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw proportion (VALUE_X/SAMPLE_SIZE) reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000)),
                     breaks = pretty(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000, n = 20))
  
#' ## Entamoeba histolytica: data availability and quality
i <- 5
# Data availability

#+ fig.height=3
plot_data(subset(all_dta[[i]], FLAG == 0), by = "REG2")
#+ fig.height=4
plot_data(subset(all_dta[[i]], FLAG == 0), by = "SUB2")
#+ fig.height=4
plot_world_imputation(subset(all_dta[[i]], FLAG == 0), sub = "SUB2")
#+ fig.height=4
plot_world_data(subset(all_dta[[i]], FLAG == 0))

# Data quality

#+ warning=FALSE
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi, group = REG2)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "loess") +
  facet_wrap(~REG2) +
  theme_bw()

#+ fig.height=8
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, aes(col = REG2)) +
  facet_wrap(SYNDROMTYPE ~ AGE) +
  theme_bw()

ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = REG2, y = RAW_PROP_100000, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw proportion (VALUE_X/SAMPLE_SIZE) reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000)),
                     breaks = pretty(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000, n = 20))
  
#'## EPEC-typical: data availability and quality 

# Data availability
i <- 6
#+ fig.height=3
plot_data(subset(all_dta[[i]], FLAG == 0), by = "REG2")
#+ fig.height=4
plot_data(subset(all_dta[[i]], FLAG == 0), by = "SUB2")
#+ fig.height=4
plot_world_imputation(subset(all_dta[[i]], FLAG == 0), sub = "SUB2")
#+ fig.height=4
plot_world_data(subset(all_dta[[i]], FLAG == 0))

# Data quality

#+ warning=FALSE
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi, group = REG2)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "loess") +
  facet_wrap(~REG2) +
  theme_bw()

#+ fig.height=8
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, aes(col = REG2)) +
  facet_wrap(SYNDROMTYPE ~ AGE) +
  theme_bw()

ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = REG2, y = RAW_PROP_100000, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw proportion (VALUE_X/SAMPLE_SIZE) reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000)),
                     breaks = pretty(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000, n = 20))

#'## EPEC-atypical: data availability and quality 

# Data availability
i <- 7
#+ fig.height=3
plot_data(subset(all_dta[[i]], FLAG == 0), by = "REG2")
#+ fig.height=4
plot_data(subset(all_dta[[i]], FLAG == 0), by = "SUB2")
#+ fig.height=4
plot_world_imputation(subset(all_dta[[i]], FLAG == 0), sub = "SUB2")
#+ fig.height=4
plot_world_data(subset(all_dta[[i]], FLAG == 0))

# Data quality

#+ warning=FALSE
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi, group = REG2)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "loess") +
  facet_wrap(~REG2) +
  theme_bw()

#+ fig.height=8
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, aes(col = REG2)) +
  facet_wrap(SYNDROMTYPE ~ AGE) +
  theme_bw()

ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = REG2, y = RAW_PROP_100000, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw proportion (VALUE_X/SAMPLE_SIZE) reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000)),
                     breaks = pretty(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000, n = 20))

#' ## ETEC-LT: data availability and quality
i <- 8
#+ fig.height=3
plot_data(subset(all_dta[[i]], FLAG == 0), by = "REG2")
#+ fig.height=4
plot_data(subset(all_dta[[i]], FLAG == 0), by = "SUB2")
#+ fig.height=4
plot_world_imputation(subset(all_dta[[i]], FLAG == 0), sub = "SUB2")
#+ fig.height=4
plot_world_data(subset(all_dta[[i]], FLAG == 0))

# ETEC-LT Data quality

#+ warning=FALSE
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi, group = REG2)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "loess") +
  facet_wrap(~REG2) +
  theme_bw()

#+ fig.height=8
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, aes(col = REG2)) +
  facet_wrap(SYNDROMTYPE ~ AGE) +
  theme_bw()

ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = REG2, y = RAW_PROP_100000, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw proportion (VALUE_X/SAMPLE_SIZE) reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000)),
                     breaks = pretty(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000, n = 20))

#' ## ETEC-ST: data availability and quality
i <- 9
#+ fig.height=3
plot_data(subset(all_dta[[i]], FLAG == 0), by = "REG2")
#+ fig.height=4
plot_data(subset(all_dta[[i]], FLAG == 0), by = "SUB2")
#+ fig.height=4
plot_world_imputation(subset(all_dta[[i]], FLAG == 0), sub = "SUB2")
#+ fig.height=4
plot_world_data(subset(all_dta[[i]], FLAG == 0))

# ETEC-ST Data quality

#+ warning=FALSE
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi, group = REG2)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "loess") +
  facet_wrap(~REG2) +
  theme_bw()

#+ fig.height=8
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, aes(col = REG2)) +
  facet_wrap(SYNDROMTYPE ~ AGE) +
  theme_bw()

ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = REG2, y = RAW_PROP_100000, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw proportion (VALUE_X/SAMPLE_SIZE) reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000)),
                     breaks = pretty(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000, n = 20))

#'## Giardia: data availability and quality 
i <- 10
# Data availability

#+ fig.height=3
plot_data(subset(all_dta[[i]], FLAG == 0), by = "REG2")
#+ fig.height=4
plot_data(subset(all_dta[[i]], FLAG == 0), by = "SUB2")
#+ fig.height=4
plot_world_imputation(subset(all_dta[[i]], FLAG == 0), sub = "SUB2")
#+ fig.height=4
plot_world_data(subset(all_dta[[i]], FLAG == 0))

# Data quality

#+ warning=FALSE
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi, group = REG2)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "loess") +
  facet_wrap(~REG2) +
  theme_bw()

#+ fig.height=8
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, aes(col = REG2)) +
  facet_wrap(SYNDROMTYPE ~ AGE) +
  theme_bw()

ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = REG2, y = RAW_PROP_100000, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw proportion (VALUE_X/SAMPLE_SIZE) reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000)),
                     breaks = pretty(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000, n = 20))
  
#'## Norovirus: data availability and quality

# Data availability 
i <- 11
#+ fig.height=3
plot_data(subset(all_dta[[i]], FLAG == 0), by = "REG2")
#+ fig.height=4
plot_data(subset(all_dta[[i]], FLAG == 0), by = "SUB2")
#+ fig.height=4
plot_world_imputation(subset(all_dta[[i]], FLAG == 0), sub = "SUB2")
#+ fig.height=4
plot_world_data(subset(all_dta[[i]], FLAG == 0))

# Data quality

#+ warning=FALSE
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi, group = REG2)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "loess") +
  facet_wrap(~REG2) +
  theme_bw()

#+ fig.height=8
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, aes(col = REG2)) +
  facet_wrap(SYNDROMTYPE ~ AGE) +
  theme_bw()

ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = REG2, y = RAW_PROP_100000, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw proportion (VALUE_X/SAMPLE_SIZE) reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000)),
                     breaks = pretty(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000, n = 20))
  
#'## Rotavirus: data availability and quality
i <- 12
# Data availability 
#+ fig.height=3
plot_data(subset(all_dta[[i]], FLAG == 0), by = "REG2")
#+ fig.height=4
plot_data(subset(all_dta[[i]], FLAG == 0), by = "SUB2")
#+ fig.height=4
plot_world_imputation(subset(all_dta[[i]], FLAG == 0), sub = "SUB2")
#+ fig.height=4
plot_world_data(subset(all_dta[[i]], FLAG == 0))

# Data quality

#+ warning=FALSE
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi, group = REG2)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "loess") +
  facet_wrap(~REG2) +
  theme_bw()

#+ fig.height=8
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, aes(col = REG2)) +
  facet_wrap(SYNDROMTYPE ~ AGE) +
  theme_bw()

ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = REG2, y = RAW_PROP_100000, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw proportion (VALUE_X/SAMPLE_SIZE) reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000)),
                     breaks = pretty(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000, n = 20))
  
#'## Salmonella: data availability and quality

# Data availability
i <- 13
#+ fig.height=3
plot_data(subset(all_dta[[i]], FLAG == 0), by = "REG2")
#+ fig.height=4
plot_data(subset(all_dta[[i]], FLAG == 0), by = "SUB2")
#+ fig.height=4
plot_world_imputation(subset(all_dta[[i]], FLAG == 0), sub = "SUB2")
#+ fig.height=4
plot_world_data(subset(all_dta[[i]], FLAG == 0))

# Data quality

#+ warning=FALSE
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi, group = REG2)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "loess") +
  facet_wrap(~REG2) +
  theme_bw()

#+ fig.height=8
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, aes(col = REG2)) +
  facet_wrap(SYNDROMTYPE ~ AGE) +
  theme_bw()

ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = REG2, y = RAW_PROP_100000, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw proportion (VALUE_X/SAMPLE_SIZE) reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000)),
                     breaks = pretty(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000, n = 20))
  
#'## Shigella: data availability and quality

# Data availability
i <- 14
#+ fig.height=3
plot_data(subset(all_dta[[i]], FLAG == 0), by = "REG2")
#+ fig.height=4
plot_data(subset(all_dta[[i]], FLAG == 0), by = "SUB2")
#+ fig.height=4
plot_world_imputation(subset(all_dta[[i]], FLAG == 0), sub = "SUB2")
#+ fig.height=4
plot_world_data(subset(all_dta[[i]], FLAG == 0))

# Data quality

#+ warning=FALSE
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi, group = REG2)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "loess") +
  facet_wrap(~REG2) +
  theme_bw()

#+ fig.height=8
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, aes(col = REG2)) +
  facet_wrap(SYNDROMTYPE ~ AGE) +
  theme_bw()

ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = REG2, y = RAW_PROP_100000, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw proportion (VALUE_X/SAMPLE_SIZE) reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000)),
                     breaks = pretty(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000, n = 20))
  
#'## Shiga toxin-producing E. coli: data availability and quality

# Data availability
i <- 15
#+ fig.height=3
plot_data(subset(all_dta[[i]], FLAG == 0), by = "REG2")
#+ fig.height=4
plot_data(subset(all_dta[[i]], FLAG == 0), by = "SUB2")
#+ fig.height=4
plot_world_imputation(subset(all_dta[[i]], FLAG == 0), sub = "SUB2")
#+ fig.height=4
plot_world_data(subset(all_dta[[i]], FLAG == 0))

# Data quality

#+ warning=FALSE
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi, group = REG2)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "loess") +
  facet_wrap(~REG2) +
  theme_bw()

#+ fig.height=8
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, aes(col = REG2)) +
  facet_wrap(SYNDROMTYPE ~ AGE) +
  theme_bw()

ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = REG2, y = RAW_PROP_100000, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw proportion (VALUE_X/SAMPLE_SIZE) reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000)),
                     breaks = pretty(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000, n = 20))
  
#'## Vibrio cholerae:  data availability & quality
  
# Data availability
i <- 16
#+ fig.height=3
plot_data(subset(all_dta[[i]], FLAG == 0), by = "REG2")
#+ fig.height=4
plot_data(subset(all_dta[[i]], FLAG == 0), by = "SUB2")
#+ fig.height=4
plot_world_imputation(subset(all_dta[[i]], FLAG == 0), sub = "SUB2")
#+ fig.height=4
plot_world_data(subset(all_dta[[i]], FLAG == 0))

# Data quality

#+ warning=FALSE
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi, group = REG2)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "loess") +
  facet_wrap(~REG2) +
  theme_bw()

#+ fig.height=8
ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = YEAR, y = yi)) +
  geom_point(aes(col = REG2)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, aes(col = REG2)) +
  facet_wrap(SYNDROMTYPE ~ AGE) +
  theme_bw()

ggplot(subset(all_dta[[i]], FLAG == 0), aes(x = REG2, y = RAW_PROP_100000, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw proportion (VALUE_X/SAMPLE_SIZE) reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000)),
                     breaks = pretty(subset(all_dta[[i]], FLAG == 0)$RAW_PROP_100000, n = 20))
  
#'# Age and syndrome classification 
#' ## General
dta_cat <- bind_rows(all_dta) %>% filter(FLAG == 0)
kable(addmargins(table(dta_cat$pathogen, dta_cat$AGE_GROUP),FUN = list(list(Total = sum)),c(1,2)))
kable(addmargins(table(dta_cat$pathogen, dta_cat$AGE),FUN = list(list(Total = sum)),c(1,2)))

kable(addmargins(table(dta_cat$pathogen, dta_cat$SYNDROM),FUN = list(list(Total = sum)),c(1,2)))
kable(addmargins(table(dta_cat$pathogen, dta_cat$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#' ## Campylobacter
kable(addmargins(table(subset(all_dta[[1]], FLAG== 0)$AGE, subset(all_dta[[1]], FLAG== 0)$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#' ## Cryptosporidium
kable(addmargins(table(subset(all_dta[[2]], FLAG== 0)$AGE, subset(all_dta[[2]], FLAG== 0)$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#' ## Cyclospora
kable(addmargins(table(subset(all_dta[[3]], FLAG== 0)$AGE, subset(all_dta[[3]], FLAG== 0)$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#' ## Enteroaggregative E coli
kable(addmargins(table(subset(all_dta[[4]], FLAG== 0)$AGE, subset(all_dta[[4]], FLAG== 0)$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#' ## Entamoeba histolytica
kable(addmargins(table(subset(all_dta[[5]], FLAG== 0)$AGE, subset(all_dta[[5]], FLAG== 0)$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#' ## Enteropathogenic E coli - Typical
kable(addmargins(table(subset(all_dta[[6]], FLAG== 0)$AGE, subset(all_dta[[6]], FLAG== 0)$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#' ## Enteropathogenic E coli - Atypical
kable(addmargins(table(subset(all_dta[[7]], FLAG== 0)$AGE, subset(all_dta[[7]], FLAG== 0)$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#' ## Enterotoxigenic E coli - LT
kable(addmargins(table(subset(all_dta[[8]], FLAG== 0)$AGE, subset(all_dta[[8]], FLAG== 0)$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#' ## Enterotoxigenic E coli - ST
kable(addmargins(table(subset(all_dta[[9]], FLAG== 0)$AGE, subset(all_dta[[9]], FLAG== 0)$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#' ## Giardia
kable(addmargins(table(subset(all_dta[[10]], FLAG== 0)$AGE, subset(all_dta[[10]], FLAG== 0)$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#' ## Norovirus
kable(addmargins(table(subset(all_dta[[11]], FLAG== 0)$AGE, subset(all_dta[[11]], FLAG== 0)$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#' ## Rotavirus
kable(addmargins(table(subset(all_dta[[12]], FLAG== 0)$AGE, subset(all_dta[[12]], FLAG== 0)$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#' ## Salmonella
kable(addmargins(table(subset(all_dta[[13]], FLAG== 0)$AGE, subset(all_dta[[13]], FLAG== 0)$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#' ## Shigella
kable(addmargins(table(subset(all_dta[[14]], FLAG== 0)$AGE, subset(all_dta[[14]], FLAG== 0)$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#' ## Shiga toxin-producing E coli
kable(addmargins(table(subset(all_dta[[15]], FLAG== 0)$AGE, subset(all_dta[[15]], FLAG== 0)$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#' ## Vibrio cholerae
kable(addmargins(table(subset(all_dta[[16]], FLAG== 0)$AGE, subset(all_dta[[16]], FLAG== 0)$SYNDROMTYPE),FUN = list(list(Total = sum)),c(1,2)))

#'# Frequency tables of additional variables ---
#' ## DX: diagnosis method by pathogen
dta_cat <- bind_rows(all_dta) %>% filter(FLAG == 0)
kable(addmargins(table(dta_cat$pathogen, dta_cat$DX),FUN = list(list(Total = sum)),c(1,2)))

Date <- format(Sys.Date(), "%Y%m%d")
es <- all_dta
names(es) <- pathogens
saveRDS(es, paste0("es_LMIC_",Date,".RDS"))

#'## presence of HI countries in data ---
# zero_cases<- read_xlsx("//sciensano.be/fs/11401_LifeChron_FERG_DWH/Estimates/00_ADMIN/Endemic_countries.xlsx")%>%
#   select(REG2, SUB2, ISO3, Country, edtf_diarrheal) %>% 
#   rename(COUNTRY=ISO3, COUNTRY_LABEL = Country, DISEASEFREE = edtf_diarrheal)
# 
# HI_countries <- zero_cases %>%
#   filter(DISEASEFREE == 0) %>%
#   select(COUNTRY) %>%
#   unlist()
# 
# dta_all_HI_zero <- dta_cat %>%
#   filter(ISO3 %in% HI_countries) %>%
#   filter(VALUE_X == 0)
# 
# dta_all_HI_nonzero <- dta_cat %>%
#   filter(ISO3 %in% HI_countries) %>%
#   filter(VALUE_X != 0)
# 
# kable(addmargins(table(dta_all_HI_zero$ISO3, dta_all_HI_zero$pathogen),FUN = list(list(Total = sum)),c(1,2)))
# kable(addmargins(table(dta_all_HI_nonzero$ISO3, dta_all_HI_nonzero$pathogen),FUN = list(list(Total = sum)),c(1,2)))

#' ## Study types present in data
kable(table(dta_cat$SOURCE.DESIGN))

#' # Session info
sessioninfo::session_info()

##bd::render_today("01-data-v20250819.R")
