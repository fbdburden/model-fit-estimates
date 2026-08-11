#' ---
#' title: Diarrhea Etiology Incidence • Data summary
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
  read_xlsx("ferg2-edtf-diarrheal-highincome-20251011.xlsx", "Incidence") %>%
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
    OPT_DISEASE == "STEC" | OPT_DISEASE == "EHEC (STEC)" | str_detect(OPT_DISEASE, "VTEC") ~ "STEC",
    str_detect(OPT_DISEASE, "Giardia") ~ "GIAR",
    str_detect(OPT_DISEASE, "Salmonella") | OPT_DISEASE == "S. enterica" ~ "SALM",
    str_detect(OPT_DISEASE, "Norovirus") | str_detect(OPT_DISEASE, "Norwalk") ~ "NORO",
    str_detect(OPT_DISEASE, "Shigella") ~ "SHIG",
    str_detect(OPT_DISEASE, "Vibrio") ~ "VIBR",
    str_detect(OPT_DISEASE, "Rotavirus") ~ "ROTA",
    str_detect(OPT_DISEASE, "E coli") & str_detect(Comments, "Classifty as EPEC.") ~ "EPEC",
    str_detect(OPT_DISEASE, "E Coli") & str_detect(Comments, "Classifty as EPEC.") ~ "EPEC",
    str_detect(OPT_DISEASE, "O157") ~ "STEC" #> BD added
  ))
subset(all_dta$dta, pathogen == "STEC")$OPT_DISEASE
# View(all_dta$dta %>% select(OPT_DISEASE, pathogen) %>% distinct)
# View(all_dta$dta %>% select(OPT_DISEASE, pathogen) %>% distinct %>% filter(is.na(pathogen)))
tbl <- table(all_dta$dta$pathogen, useNA = "always") 
kable(tbl, 
      col.names = c("Pathogen", "Number of datapoints"))

# For STEC some studies need to be summed
STEC <- all_dta$dta %>%
  filter(SOURCE_ID %in% c(3,18,22,29) & pathogen == "STEC" & Comments != "Exclude")
i <- c(26:41)  
STEC[, i] <- apply(STEC[, i], 2, function(x) as.numeric(as.character(x)))
col_STEC <- colnames(STEC)
STEC <- aggregate( STEC[,26:41], STEC[,c(1:15,42)], FUN = sum )
col_STEC <- setdiff(col_STEC, colnames(STEC))
for (c in col_STEC){
  STEC[[c]] <- NA
}

all_dta$dta <- all_dta$dta %>%
  filter(!(SOURCE_ID %in% c(3,18,22,29) & pathogen == "STEC" & Comments != "Exclude"))
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

# EDTF agreed that this sub national study could be used as national study
all_dta$dta <- all_dta$dta %>%
  mutate(
    REF_LOC_LEVEL = case_when(
      SOURCE_ID %in% c(1,34) & REF_LOCATION == "Japan" ~ "National",  
      .default = REF_LOC_LEVEL))

# Study 33, the ref sample size should be calculated using the cases and incidence rate, this information needs to be added to the data
Study33 <- subset(all_dta$dta, SOURCE_ID ==33)
Study33 <- subset(Study33, !(OPT_DISEASE == "Rotavirus" & is.na(OPT_SEROTYPE)))
Study33 <- Study33 %>% 
  mutate(VALUE_X = case_when(
    OPT_DISEASE == "Campylobacter" ~ 32,
    OPT_DISEASE == "Cryptosporidium" ~ 3,
    OPT_DISEASE == "EAEC" ~ 18,
    OPT_DISEASE == "EPEC" ~ 1,
    OPT_DISEASE == "ETEC" ~ 10,
    OPT_DISEASE == "Giardia" ~ 2,
    OPT_DISEASE == "Rotavirus" & OPT_SEROTYPE == "Group A" ~ 26,
    OPT_DISEASE == "Rotavirus" & OPT_SEROTYPE == "Group C" ~ 2,
    OPT_DISEASE == "Salmonella" ~ 8,
    OPT_DISEASE == "Shigella" ~ 1,
    OPT_DISEASE == "VTEC" ~ 3))
Study33$REF_SAMPLE_SIZE <- Study33$VALUE_X*Study33$VALUE_DENOM/Study33$VALUE_MEAN

# next to this the two observations on rotavirus need to be summed
Study33_rota <- subset(Study33, OPT_DISEASE == "Rotavirus")
Study33_rota$VALUE_MEAN <- NA
col_Study33_rota <- colnames(Study33_rota)
Study33_rota <- aggregate( Study33_rota[,26:41], Study33_rota[,c(1:15,42)], FUN = sum )
col_Study33_rota <- setdiff(col_Study33_rota, colnames(Study33_rota))
for (c in col_Study33_rota){
  Study33_rota[[c]] <- NA
}
Study33_rota$FLAG <- 0

Study33 <- subset(Study33, OPT_DISEASE != "Rotavirus")
all_dta$dta <- subset(all_dta$dta, SOURCE_ID != 33)
all_dta$dta <- rbind(all_dta$dta, Study33, Study33_rota)

# CC: for NLD study was added, only DP from Shigella need to be kept
all_dta$dta$FLAG[
  all_dta$dta$SOURCE_ID == "CC_48" &
    all_dta$dta$OPT_DISEASE != "Shigella"] <-5

all_dta$dta$COUNTRY <- all_dta$dta$REF_LOCATION_ISO3
sum(is.na(all_dta$dta$REF_LOCATION_ISO3))
all_dta$dta$REF_LOCATION_ISO3 <- if_else(is.na(all_dta$dta$REF_LOCATION_ISO3) & all_dta$dta$REF_LOCATION == "Switzerland",
                                         "CHE",
                                         all_dta$dta$REF_LOCATION_ISO3)
all_dta$dta$ISO3 <- all_dta$dta$REF_LOCATION_ISO3
all_dta$dta$ID <- all_dta$dta$SOURCE_ID
all_dta$dta$REG2 <-
  FERG2:::countries$REG2[match(all_dta$dta$ISO3, FERG2:::countries$ISO3)]
all_dta$dta$SUB2 <-
  FERG2:::countries$SUB2[match(all_dta$dta$ISO3, FERG2:::countries$ISO3)]

all_dta$dta$YEAR <- rowMeans(cbind(all_dta$dta$REF_YEAR_START, all_dta$dta$REF_YEAR_END))

all_dta$dta <- all_dta$dta %>% 
  mutate(REF_AGE_START = 0, 
         REF_AGE_END = 125, 
         REF_SEX = case_when(is.na(REF_SEX) ~ "All sexes", 
                             TRUE ~ REF_SEX))

all_dta$dta$REF_YEAR_START2 <- all_dta$dta$REF_YEAR_START
all_dta$dta$REF_YEAR_END2 <- all_dta$dta$REF_YEAR_END

all_dta$dta$REF_YEAR_START <- round(all_dta$dta$YEAR)
all_dta$dta$REF_YEAR_END <- round(all_dta$dta$YEAR)

if (sum(is.na(all_dta$dta$REF_SAMPLE_SIZE))>0) {
  all_dta$dta <- add_pop(all_dta$dta)
} else {all_dta$dta$POP <- NA_real_
}
# 1 missing with flag = 5 so not important

all_dta$dta$REF_YEAR_START <- all_dta$dta$REF_YEAR_START2
all_dta$dta$REF_YEAR_END <- all_dta$dta$REF_YEAR_END2

all_dta$dta$REF_YEAR_START2 <- NULL
all_dta$dta$REF_YEAR_END2 <- NULL

(id <- grep("to", all_dta$dta$VALUE_X))
# cbind(all_dta$dta$VALUE_X, as.numeric(all_dta$dta$VALUE_X))
all_dta$dta <- all_dta$dta %>%
  mutate(REF_SAMPLE_SIZE = case_when(is.na(REF_SAMPLE_SIZE) 
                                     & REF_LOC_LEVEL=="National"
                                     ~ POP, 
                                     TRUE ~ as.numeric(REF_SAMPLE_SIZE)),
         VALUE_X = case_when(
           VALUE_X == "75271 to 214386" ~ mean(c(75271, 214386)),
           VALUE_X == "21572 to 103005" ~ mean(c(21572, 103005)),
           VALUE_X == "283052 to 599650" ~ mean(c(283052, 599650)),
           .default = as.numeric(VALUE_X)))
all_dta$dta$VALUE_X[id]

#  E-mail Elaine 05/12 saying that if denom value is present it is always a rate, while if denom is not present it is the number of cases
all_dta$dta <- all_dta$dta %>% 
  mutate(VALUE_XMEAN = case_when(
    SOURCE_ID == 33 ~ VALUE_X,
    !is.na(VALUE_X) & is.na(VALUE_MEAN) & is.na(VALUE_MEDIAN) ~ VALUE_X,
    !is.na(VALUE_MEAN) & is.na(VALUE_X) & is.na(VALUE_MEDIAN) ~ VALUE_MEAN,
    is.na(VALUE_X) & is.na(VALUE_MEAN) & !is.na(VALUE_MEDIAN) ~ VALUE_MEDIAN
  )) %>%
  mutate(VALUE_X = case_when(
    is.na(VALUE_DENOM) | SOURCE_ID == 33 ~ VALUE_XMEAN, 
    .default=NA),
    VALUE_MEAN = case_when(
      !is.na(VALUE_DENOM) ~ VALUE_XMEAN, 
      .default=NA))

sum(is.na(all_dta$dta$VALUE_X))
all_dta$dta <- all_dta$dta %>% 
  mutate(VALUE_X = case_when(
    is.na(VALUE_X) & !is.na(VALUE_MEAN) &
      !is.na(VALUE_DENOM) & !is.na(REF_SAMPLE_SIZE)
    ~ REF_SAMPLE_SIZE * VALUE_MEAN / VALUE_DENOM,
    TRUE ~ VALUE_X
  ))
sum(is.na(all_dta$dta$VALUE_X))

Territories <- read_xlsx("Territories_R_20250221.xlsx")
Flag_territory <- unlist(Territories)

all_dta$dta$FLAG_REF_LOCATION <- as.integer(apply(sapply(Flag_territory, function(x) grepl(x, all_dta$dta$REF_LOCATION, ignore.case = TRUE)), 1, any))
all_dta$dta$FLAG_REF_NOTES <- as.integer(apply(sapply(Flag_territory, function(x) grepl(x, all_dta$dta$REF_NOTES, ignore.case = TRUE)), 1, any))
all_dta$dta$FLAG_SOURCE_TITLE <- as.integer(apply(sapply(Flag_territory, function(x) grepl(x, all_dta$dta$SOURCE_TITLE, ignore.case = TRUE)), 1, any))
all_dta$dta$FLAG_TERRITORY <- if_else(all_dta$dta$FLAG_REF_LOCATION + all_dta$dta$FLAG_REF_NOTES + all_dta$dta$FLAG_SOURCE_TITLE >=1 , 1, 0)

all_dta$dta$FLAG <- if_else(all_dta$dta$FLAG_TERRITORY == 1 & all_dta$dta$FLAG == 0, 
                    1, 
                    all_dta$dta$FLAG)

all_dta$dta$PERSONYEARS100 <- all_dta$dta$REF_SAMPLE_SIZE / 1e5
all_dta$dta$PERSONYEARS100[all_dta$dta$PERSONYEARS100 == 0] <- NA


## filter data ----
# Studies excluded: from countries which don't belong to WHO 
sum(is.na(all_dta$dta$REG2))
# Studies excluded: before 1990 
length(which(all_dta$dta$YEAR<1990))

data.frame(subset(all_dta$dta, FLAG != 0))

tbl <- table(subset(all_dta$dta, FLAG == 0)$pathogen, useNA ="always") 
kable(tbl, 
      col.names = c("Pathogen","Number of datapoints"))

## estimate incidence ----
pathogens <- c("CAMP","CRYP","CYCL","EAEC","EPEC","ETEC","GIAR",
               "NORO","ROTA","SALM","SHIG","STEC","VIBR")
es <- list()

for (p in pathogens) {
  print(p)
  all_dta[[p]] <- all_dta$dta %>% filter(pathogen == p)
  
  ## .. "IRLN" for the log-transformed incidence rate
  es[[p]] <- escalc(xi = VALUE_X, ti = PERSONYEARS100,
                    measure = "IRLN", data = all_dta[[p]])
  es[[p]]$sei <- sqrt(es[[p]]$vi)
  es[[p]]$RAW_INC <-1e5 * es[[p]]$VALUE_X / es[[p]]$REF_SAMPLE_SIZE
  es[["all"]] <- rbind(es[["all"]], es[[p]])
  es[[p]] <- es[[p]] %>% filter(!is.na(yi))
}

all_dta$Exclude <- subset(all_dta$dta, is.na(pathogen))
all_dta$Exclude$yi <- NA
all_dta$Exclude$vi <- NA
all_dta$Exclude$sei <- NA
all_dta$Exclude$RAW_INC <- NA

es$all <- rbind(es$all, all_dta$Exclude) # es$all should be used to list studies for countries

sum(is.na(es$all$yi))
es[["all_yi"]] <- es$all %>% filter(!is.na(yi))

tbl <- table(subset(es$all_yi, FLAG == 0)$pathogen, useNA ="always") 
kable(tbl, 
      col.names = c("Pathogen","Number of datapoints"))

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
ggplot(es[[p]], aes(x = REG2, y = RAW_INC, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw incidence per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_INC)),
                     breaks = pretty(es[[p]]$RAW_INC, n = 20))
#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_INC, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported incidence per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$REG2), FUN=mean_ci) 
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
ggplot(es[[p]], aes(x = REG2, y = RAW_INC, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw incidence per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_INC)),
                     breaks = pretty(es[[p]]$RAW_INC, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_INC, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported incidence per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$REG2), FUN=mean_ci) 
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
ggplot(es[[p]], aes(x = REG2, y = RAW_INC, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw incidence per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_INC)),
                     breaks = pretty(es[[p]]$RAW_INC, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_INC, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported incidence per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$REG2), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("WHO region", "Mean", "Lower", "Upper"))


#' ## Enteroaggregative E. coli: data availability and quality
p <- "EAEC"

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
ggplot(es[[p]], aes(x = REG2, y = RAW_INC, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw incidence per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_INC)),
                     breaks = pretty(es[[p]]$RAW_INC, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_INC, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported incidence per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$REG2), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("WHO region", "Mean", "Lower", "Upper"))


#' ## Entamoeba histolytica: data availability and quality
p <- "ENTA"

#Data availability 
print("No data available for Entamoeba histolytica (ENTA)")

#'## Enteropathogenic E. coli: data availability and quality 
p <- "EPEC"

#Data availability 
print("No data available for Enteropathogenic E. coli (EPEC)")
#+ fig.height=3
plot_data(es[[p]], by = "REG2")
#+ fig.height=4
plot_data(es[[p]], by = "SUB2")
#+ fig.height=4
plot_world_imputation(es[[p]], sub = "SUB2")
#+ fig.height=4
plot_world_data(es[[p]])
#+ fig.height=4
ggplot(es[[p]], aes(x = REG2, y = RAW_INC, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw incidence per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_INC)),
                     breaks = pretty(es[[p]]$RAW_INC, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_INC, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported incidence per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$REG2), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("WHO region", "Mean", "Lower", "Upper"))


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
ggplot(es[[p]], aes(x = REG2, y = RAW_INC, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw incidence per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_INC)),
                     breaks = pretty(es[[p]]$RAW_INC, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_INC, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported incidence per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$REG2), FUN=mean_ci) 
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
ggplot(es[[p]], aes(x = REG2, y = RAW_INC, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw incidence per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_INC)),
                     breaks = pretty(es[[p]]$RAW_INC, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_INC, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported incidence per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$REG2), FUN=mean_ci) 
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
ggplot(es[[p]], aes(x = REG2, y = RAW_INC, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw incidence per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_INC)),
                     breaks = pretty(es[[p]]$RAW_INC, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_INC, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported incidence per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$REG2), FUN=mean_ci) 
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
ggplot(es[[p]], aes(x = REG2, y = RAW_INC, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw incidence per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_INC)),
                     breaks = pretty(es[[p]]$RAW_INC, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_INC, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported incidence per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$REG2), FUN=mean_ci) 
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
ggplot(es[[p]], aes(x = REG2, y = RAW_INC, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw incidence per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_INC)),
                     breaks = pretty(es[[p]]$RAW_INC, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_INC, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported incidence per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$REG2), FUN=mean_ci) 
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
ggplot(es[[p]], aes(x = REG2, y = RAW_INC, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw incidence per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_INC)),
                     breaks = pretty(es[[p]]$RAW_INC, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_INC, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported incidence per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$REG2), FUN=mean_ci) 
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
ggplot(es[[p]], aes(x = REG2, y = RAW_INC, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw incidence per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_INC)),
                     breaks = pretty(es[[p]]$RAW_INC, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_INC, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported incidence per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$REG2), FUN=mean_ci) 
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
ggplot(es[[p]], aes(x = REG2, y = RAW_INC, group = REG2, color=REG2)) +
  geom_jitter(width = 0.2, height = 0, ) +  # Adds random noise horizontally for better distribution
  stat_summary(fun.min = function(z) { quantile(z, 0.25) },
               fun.max = function(z) { quantile(z, 0.75) },
               fun = median, geom = "pointrange", color = "black", size = 0.5) +
  labs(y = "Raw incidence per 100000 reported in the studies included", x = "WHO Region") +
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(limits = c(0, max(es[[p]]$RAW_INC)),
                     breaks = pretty(es[[p]]$RAW_INC, n = 20))

#+ fig.height=4
ggplot(es[[p]], aes(x=YEAR, y=RAW_INC, color=REG2)) + 
  geom_point() + 
  labs(y="Raw reported incidence per 100000", x="Year")

#Raw results
T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$YEAR), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("Year", "Mean", "Lower", "Upper"))

T <- aggregate(es[[p]]$RAW_INC, list(es[[p]]$REG2), FUN=mean_ci) 
T1 <- T$Group.1
T2 <- unlist(T$x) %>% as.data.frame()
T <- cbind(T1, T2)
kable(T, col.names=c("WHO region", "Mean", "Lower", "Upper"))

#' # Session info
Date <- format(Sys.Date(), "%Y%m%d")
saveRDS(es, paste0("es_HIC_INC_", Date, ".rds"))
sessioninfo::session_info()

##bd::render_today("01-data-hic-inc.R")
