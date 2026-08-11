#' ---
#' title: Diarrheal Etiology • Normalization • MIX+ROTA+STRAIN+ZERO+AVG
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

## global options ----
knitr::opts_chunk$set(fig.width = 10)

## initialization ----
dir_avg_v4 <- "../ESTIMATES_AVG_ROTA_V4/01-average-output/"
dir_asy_v4 <- "../ESTIMATES_ASYMPT_AF_ROTA_V4/02-models-lmic/"
dir_avg_v5 <- "../ESTIMATES_AVG_ROTA_V5/01-average-output/"
dir_asy_v5 <- "../ESTIMATES_ASYMPT_AF_ROTA_V5/02-models-lmic/"

hazards <-
data.frame(
  pathogens = c("CAMP", "CRYP", "CYCL", "EAEC", "ENTA", "EPTP", "EPAP", "ETLT", "ETST",
                "GIAR", "NORO", "ROTA", "SALM", "SHIG", "STEC", "VIBR"), 
  paths = c(paste0(dir_avg_v4, "01-campy"),
            paste0(dir_asy_v5, "02-crypto"),
            paste0(dir_asy_v5, "03-cyclo"),
            paste0(dir_asy_v4, "04-eaec"),
            paste0(dir_asy_v5, "05-enta"),
            paste0(dir_asy_v4, "06-eptp"),
            paste0(dir_asy_v4, "07-epap"),
            paste0(dir_asy_v4, "08-etlt"),
            paste0(dir_asy_v4, "09-etst"),
            paste0(dir_avg_v5, "10-giardia"),
            paste0(dir_asy_v4, "11-noro"),
            paste0(dir_asy_v4, "12-rota"),
            paste0(dir_asy_v4, "13-salmo"),
            paste0(dir_asy_v4, "14-shig"),
            paste0(dir_avg_v4, "15-stec"),
            paste0(dir_asy_v4, "16-vibrio")))

hazards_op <-
  subset(hazards, !(pathogens %in% "EPAP"))
nrow(hazards_op)

hazards_ip <-
  subset(hazards, !(pathogens %in% c("EPAP", "ETLT", "GIAR")))
nrow(hazards_ip)

#' # Simulations

## load simulations ----

# outpatients
rm(list = ls(pattern = "sim")); gc()
sim_all_outpat <- list()

for (i in seq(nrow(hazards_op))) {
  sim_path3 <- readRDS(dir(hazards_op$paths[i], pattern = "sim_all_Outpatient_Older5", full.names = TRUE))
  sim_path4 <- readRDS(dir(hazards_op$paths[i], pattern = "sim_all_Outpatient_Younger5", full.names = TRUE))
  sim_path  <- rbind(sim_path3, sim_path4)
  sim_path$hazard <- hazards_op$pathogens[i]
  h <- hazards_op$pathogens[i]
  sim_all_outpat[[h]] <- sim_path
}

sim_all_outpat <-
lapply(sim_all_outpat,
  function(x) {names(x)[names(x) == "edtf_cholera"] <- "edtf_diarrheal"; return(x)})
sim_all_outpat <- lapply(sim_all_outpat, function(x) subset(x, edtf_diarrheal == 1))

qs::qsave(sim_all_outpat, "sim_all_outpat_lmic.qs")

## inpatients
rm(list = ls(pattern = "sim")); gc()
sim_all_inpat <- list()

for (i in seq(nrow(hazards_ip))) {
  sim_path1 <- readRDS(dir(hazards_ip$paths[i], pattern = "sim_all_Inpatient_Older5", full.names = TRUE))
  sim_path2 <- readRDS(dir(hazards_ip$paths[i], pattern = "sim_all_Inpatient_Younger5", full.names = TRUE))
  sim_path  <- rbind(sim_path1, sim_path2)
  sim_path$hazard <- hazards_ip$pathogens[i]
  h <- hazards_ip$pathogens[i]
  sim_all_inpat[[h]] <- sim_path
}

sim_all_inpat <-
  lapply(sim_all_inpat,
         function(x) {names(x)[names(x) == "edtf_cholera"] <- "edtf_diarrheal"; return(x)})
sim_all_inpat <- lapply(sim_all_inpat, function(x) subset(x, edtf_diarrheal == 1))

qs::qsave(sim_all_inpat, "sim_all_inpat_lmic.qs")

#' # Proportions: build array

cols <-
  c("REG2", "SUB2", "COUNTRY", "YEAR", "SYNDROMTYPE", "AGE", "COUNTRY_LABEL", "edtf_diarrheal")

## outpatients
rm(list = ls(pattern = "sim")); gc()
type <- "outpat"
sim_all <- qs::qread(sprintf("sim_all_%s_lmic.qs", type))

## cnt simulations ----
all_cnt_prop <- list()

for (h in names(sim_all)){
  all_cnt_prop[[h]] <- t(sapply(sim_all[[h]]$PROP, mean_ci))
  all_cnt_prop[[h]] <- data.frame(all_cnt_prop[[h]])
  colnames(all_cnt_prop[[h]]) <- c("VAL_MEAN", "VAL_LWR", "VAL_UPR")
  all_cnt_prop[[h]] <- cbind(sim_all[[h]][, cols], all_cnt_prop[[h]])
  all_cnt_prop[[h]]$LOCATION <- "Country"
  all_cnt_prop[[h]]$LOCATION_NAME <- all_cnt_prop[[h]]$COUNTRY_LABEL
  all_cnt_prop[[h]]$COUNTRY_LABEL <- NULL
  all_cnt_prop[[h]]$METRIC <- "Proportion"
  str(all_cnt_prop[[h]])
}

qs::qsave(all_cnt_prop, paste0("prop_", type, "_lmic.qs"))

## inpatients
rm(list = ls(pattern = "sim")); gc()
type <- "inpat"
sim_all <- qs::qread(sprintf("sim_all_%s_lmic.qs", type))

## cnt simulations ----
all_cnt_prop <- list()

for (h in names(sim_all)){
  all_cnt_prop[[h]] <- t(sapply(sim_all[[h]]$PROP, mean_ci))
  all_cnt_prop[[h]] <- data.frame(all_cnt_prop[[h]])
  colnames(all_cnt_prop[[h]]) <- c("VAL_MEAN", "VAL_LWR", "VAL_UPR")
  all_cnt_prop[[h]] <- cbind(sim_all[[h]][, cols], all_cnt_prop[[h]])
  all_cnt_prop[[h]]$LOCATION <- "Country"
  all_cnt_prop[[h]]$LOCATION_NAME <- all_cnt_prop[[h]]$COUNTRY_LABEL
  all_cnt_prop[[h]]$COUNTRY_LABEL <- NULL
  all_cnt_prop[[h]]$METRIC <- "Proportion"
  str(all_cnt_prop[[h]])
}

qs::qsave(all_cnt_prop, paste0("prop_", type, "_lmic.qs"))

#' # Normalization: by patient type

#' ## Outpatients

rm(list = ls(pattern = "sim|Diarrheal|sum")); gc()
type <- "outpat"

## import simulations
sim_all <- qs::qread(sprintf("sim_all_%s_lmic.qs", type))

## compile proportions
Diarrheal_prop <-
  rbind(sim_all$CAMP[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$CRYP[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$CYCL[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$EAEC[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$ENTA[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$EPTP[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$EPAP[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$ETLT[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$ETST[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$GIAR[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$NORO[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$ROTA[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$SALM[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$SHIG[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$STEC[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$VIBR[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")])
nrow(Diarrheal_prop)
rm(sim_all); gc()

# function to normalize a list of proportions
normalize_list_outpat <-
  function(prop_list) {
    y <- do.call("cbind", prop_list)
    y <- y / rowSums(y)
    y <- y * (1-0.063) # OTHUNK
    y <- apply(y, 2, c, simplify = FALSE)
    return(y)
  }

# normalize proportions by COUNTRY, YEAR, and AGE
Diarrheal_norm <-
Diarrheal_prop %>%
  group_by(COUNTRY, YEAR, AGE) %>%
  mutate(SAMPLES = normalize_list_outpat(PROP)) %>%
  ungroup()
Diarrheal_norm$PROP <- NULL
Diarrheal_norm$AGE <-
  factor(x = Diarrheal_norm$AGE,
         levels = c("Age below 5", "Age above or equal 5"),
         labels = c("ch", "ad"))

# check
head(Diarrheal_norm)

# split into list
sim_all_norm <- split(data.frame(Diarrheal_norm), Diarrheal_norm$hazard)
sim_all_norm <- lapply(sim_all_norm, function(df) df[, !(names(df) == "hazard")])
length(sim_all_norm)

qs::qsave(
  sim_all_norm,
  file = sprintf("sim_all_%s_norm_lmic.qs", type))

sum_before <-
  Diarrheal_prop %>%
  group_by(YEAR, COUNTRY, AGE) %>%
  dplyr::summarize(PROP = dalymod:::list_sum(PROP)) %>%
  ungroup()
sum_before
summary(sum_before$PROP)
hist(sum_before$PROP)

sum_after <-
  Diarrheal_norm %>%
  group_by(YEAR, COUNTRY, AGE) %>%
  dplyr::summarize(PROP = dalymod:::list_sum(SAMPLES)) %>%
  ungroup()
sum_after
summary(sum_after$PROP)
hist(sum_after$PROP)

## normalisation summmary

Diarrheal_norm_mean <-
  t(sapply(Diarrheal_norm$SAMPLES, mean_ci))
Diarrheal_norm_mean <-
  cbind(Diarrheal_norm[, c("hazard", "COUNTRY", "YEAR", "AGE")],
        Diarrheal_norm_mean)
str(Diarrheal_norm_mean)

## list with normalized cnt ----
low_cnt_norm <- list()
for (i in seq_along(hazards_op$pathogens)) {
  h <- hazards_op$pathogens[i]
  df_norm <- subset(Diarrheal_norm_mean, hazard == h)
  colnames(df_norm) <-
    c("HAZARD", "COUNTRY", "YEAR", "AGE", "VAL_MEAN", "VAL_LWR", "VAL_UPR")
  df_norm$REG2 <-
    FERG2:::countries$REG2[match(df_norm$COUNTRY, FERG2:::countries$ISO3)]
  df_norm$SUB2 <-
    FERG2:::countries$SUB2[match(df_norm$COUNTRY, FERG2:::countries$ISO3)]
  df_norm$SYNDROMTYPE <- "Outpatient" #> IN/OUT
  df_norm$edtf_diarrheal <- 1
  df_norm$LOCATION <- "Country"
  df_norm$LOCATION_NAME <-
    FERG2:::countries$COUNTRY[match(df_norm$COUNTRY, FERG2:::countries$ISO3)]
  df_norm$METRIC <- "Proportion"
  df_norm$HAZARD <- NULL
  low_cnt_norm[[h]] <- df_norm
  # str(low_cnt_norm[[h]])
}
low_cnt_norm$OTHUNK <- low_cnt_norm[[1]]
low_cnt_norm$OTHUNK$VAL_MEAN <-
  low_cnt_norm$OTHUNK$VAL_LWR <-
  low_cnt_norm$OTHUNK$VAL_UPR <- 0.063 #> IN/OUT

qs::qsave(low_cnt_norm, paste0("norm_", type, "_lmic.qs"))

#' ## Inpatients

rm(list = ls(pattern = "sim|Diarrheal|sum")); gc()
type <- "inpat"

## import simulations
sim_all <- qs::qread(sprintf("sim_all_%s_lmic.qs", type))

## compile proportions
Diarrheal_prop <-
  rbind(sim_all$CAMP[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$CRYP[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$CYCL[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$EAEC[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$ENTA[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$EPTP[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$EPAP[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$ETLT[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$ETST[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$GIAR[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$NORO[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$ROTA[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$SALM[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$SHIG[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$STEC[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")],
        sim_all$VIBR[, c("hazard", "COUNTRY", "YEAR", "AGE", "PROP")])
nrow(Diarrheal_prop) / 16
rm(sim_all); gc()

# function to normalize a list of proportions
normalize_list_inpat <-
  function(prop_list) {
    y <- do.call("cbind", prop_list)
    y <- y / rowSums(y)
    y <- y * (1-0.036) # OTHUNK
    y <- apply(y, 2, c, simplify = FALSE)
    return(y)
  }

# normalize proportions by COUNTRY, YEAR, and AGE
Diarrheal_norm <-
Diarrheal_prop %>%
  group_by(COUNTRY, YEAR, AGE) %>%
  mutate(SAMPLES = normalize_list_inpat(PROP)) %>%
  ungroup()
Diarrheal_norm$PROP <- NULL
Diarrheal_norm$AGE <-
  factor(x = Diarrheal_norm$AGE,
         levels = c("Age below 5", "Age above or equal 5"),
         labels = c("ch", "ad"))

# check
head(Diarrheal_norm)

# split into list
sim_all_norm <- split(data.frame(Diarrheal_norm), Diarrheal_norm$hazard)
sim_all_norm <- lapply(sim_all_norm, function(df) df[, !(names(df) == "hazard")])
length(sim_all_norm)

qs::qsave(
  sim_all_norm,
  file = sprintf("sim_all_%s_norm_lmic.qs", type))

sum_before <-
Diarrheal_prop %>%
  group_by(YEAR, COUNTRY, AGE) %>%
  dplyr::summarize(PROP = dalymod:::list_sum(PROP)) %>%
  ungroup()
sum_before
summary(sum_before$PROP)
hist(sum_before$PROP)

sum_after <-
Diarrheal_norm %>%
  group_by(YEAR, COUNTRY, AGE) %>%
  dplyr::summarize(PROP = dalymod:::list_sum(SAMPLES)) %>%
  ungroup()
sum_after
summary(sum_after$PROP)
hist(sum_after$PROP)

## normalisation summmary

Diarrheal_norm_mean <-
  t(sapply(Diarrheal_norm$SAMPLES, mean_ci))
Diarrheal_norm_mean <-
  cbind(Diarrheal_norm[, c("hazard", "COUNTRY", "YEAR", "AGE")],
        Diarrheal_norm_mean)
str(Diarrheal_norm_mean)

## list with normalized cnt ----
low_cnt_norm <- list()
for (i in seq_along(hazards_ip$pathogens)) {
  h <- hazards_ip$pathogens[i]
  df_norm <- subset(Diarrheal_norm_mean, hazard == h)
  colnames(df_norm) <-
    c("HAZARD", "COUNTRY", "YEAR", "AGE", "VAL_MEAN", "VAL_LWR", "VAL_UPR")
  df_norm$REG2 <-
    FERG2:::countries$REG2[match(df_norm$COUNTRY, FERG2:::countries$ISO3)]
  df_norm$SUB2 <-
    FERG2:::countries$SUB2[match(df_norm$COUNTRY, FERG2:::countries$ISO3)]
  df_norm$SYNDROMTYPE <- "Inpatient" #> IN/OUT
  df_norm$edtf_diarrheal <- 1
  df_norm$LOCATION <- "Country"
  df_norm$LOCATION_NAME <-
    FERG2:::countries$COUNTRY[match(df_norm$COUNTRY, FERG2:::countries$ISO3)]
  df_norm$METRIC <- "Proportion"
  df_norm$HAZARD <- NULL
  low_cnt_norm[[h]] <- df_norm
  # str(low_cnt_norm[[h]])
}
low_cnt_norm$OTHUNK <- low_cnt_norm[[1]]
low_cnt_norm$OTHUNK$VAL_MEAN <-
  low_cnt_norm$OTHUNK$VAL_LWR <-
  low_cnt_norm$OTHUNK$VAL_UPR <- 0.036 #> IN/OUT

qs::qsave(low_cnt_norm, paste0("norm_", type, "_lmic.qs"))

#' # Session info
sessioninfo::session_info()

#bd::render_today("04-normalization-mix-rota-zero-avg.R")#~40min