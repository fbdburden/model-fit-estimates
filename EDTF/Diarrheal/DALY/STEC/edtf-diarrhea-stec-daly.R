### FERG2/EDTF/DIARRHEA/STEC/DALY

## SETTINGS ----

## required packages
library(dalymod)
library(FERG2)

## source generic FERG2 settings
source("ferg-settings.R")
set.seed(264)

## hazard specific settings
hazard <- "edtf-diarrhea-stec"
hazard_sa <- "STEC"
n_samples <- 500

## DATA ----

## source diarrhea estimates
source("../diarrhea-data.R")

## HIC
hic_inc <-
  dalymod:::get_samples(
    n_samples,
    "../../ESTIMATES_HIC/13-stec/sim_all_inc_reg_ct_20251011.rds",
    transformation = "log",
    denominator = 1e5)
summarize(hic_inc)

hic_mrt <-
  dalymod:::get_samples(
    n_samples,
    "../../ESTIMATES_HIC/13-stec/sim_all_mrt_reg_ct_20251208.rds",
    transformation = "log",
    denominator = 1e5)
summarize(hic_mrt)

## HIC / age distribution

# define age distribution
hic_age <-
  list(`2000` = c("<5" = 0.200, "5-9" = 0.128, "10-19" = 0.235, "20-64" = 0.326, "65+" = 0.110),
       `2010` = c("<5" = 0.292, "5-9" = 0.126, "10-19" = 0.169, "20-64" = 0.323, "65+" = 0.089),
       `2021` = c("<5" = 0.206, "5-9" = 0.058, "10-19" = 0.126, "20-64" = 0.458, "65+" = 0.152))
hic_age <-
  lapply(hic_age, function(x) x/sum(x))

# proportion <5, by country and year
p_ch_hic <- lapply(countries_hic, get_prob_ch_country, x = hic_age)
p_ch_hic <- do.call("rbind", p_ch_hic)
head(p_ch_hic)

# average RLE by <5 and 5+, by country and year
rle_ch_hic_w <- get_weighted_rle("<5")

rle_ad_hic <-
  rbind(cbind(AGE = "5-9", get_weighted_rle("5-9")),
        cbind(AGE = "10-19", get_weighted_rle("10-19")),
        cbind(AGE = "20-64", get_weighted_rle("20-64")),
        cbind(AGE = "65+", get_weighted_rle("65+")))

p_all_hic <- lapply(countries_hic, get_prob_all_country, x = hic_age)
p_all_hic <- do.call("rbind", p_all_hic)
head(p_all_hic)
sum(p_all_hic$PROB)

rle_ad_hic_w <-
  merge(rle_ad_hic, p_all_hic,
        by.x = c("AGE", "ISO3", "YEAR"),
        by.y = c("AGEGRP", "COUNTRY", "YEAR"))
rle_ad_hic_w <-
  group_by(rle_ad_hic_w, ISO3, YEAR) |>
  reframe(RLE = weighted.mean(RLE, PROB)) |>
  as.data.frame()
head(rle_ad_hic_w)
tail(rle_ad_hic_w)

# average LLE by <5 and 5+, by country and year
lle_ch_w <- get_weighted_lle("<5")

lle_ad <-
  merge(
    rbind(
      data.frame(AGEGRP = "5-9", get_weighted_lle("5-9")),
      data.frame(AGEGRP = "10-19", get_weighted_lle("10-19")),
      data.frame(AGEGRP = "20-64", get_weighted_lle("20-64")),
      data.frame(AGEGRP = "65+", get_weighted_lle("65+"))),
    p_all_hic,
    by.x = c("AGEGRP", "YEAR", "ISO3"),
    by.y = c("AGEGRP", "YEAR", "COUNTRY"))
head(lle_ad)
tail(lle_ad)

lle_ad_w <-
  group_by(lle_ad, YEAR, ISO3) |>
  reframe(LLE = weighted.mean(LLE, PROB))
head(lle_ad_w)

## LMIC / age distribution

# proportion <5, by country and year
p_ch_lmic <- lapply(countries_lmic, get_prob_ch_country, x = hic_age)
p_ch_lmic <- do.call("rbind", p_ch_lmic)
head(p_ch_lmic)

# average RLE by <5 and 5+, by country and year
rle_ch_lmic_w <- get_weighted_rle("<5")

rle_ad_lmic <-
  rbind(cbind(AGE = "5-9", get_weighted_rle("5-9")),
        cbind(AGE = "10-19", get_weighted_rle("10-19")),
        cbind(AGE = "20-64", get_weighted_rle("20-64")),
        cbind(AGE = "65+", get_weighted_rle("65+")))

p_all_lmic <- lapply(countries_lmic, get_prob_all_country, x = hic_age)
p_all_lmic <- do.call("rbind", p_all_lmic)
head(p_all_lmic)
sum(p_all_lmic$PROB)

rle_ad_lmic_w <-
  merge(rle_ad_lmic, p_all_lmic,
        by.x = c("AGE", "ISO3", "YEAR"),
        by.y = c("AGEGRP", "COUNTRY", "YEAR"))
rle_ad_lmic_w <-
  group_by(rle_ad_lmic_w, ISO3, YEAR) |>
  reframe(RLE = weighted.mean(RLE, PROB)) |>
  as.data.frame()
head(rle_ad_lmic_w)
tail(rle_ad_lmic_w)

## LMIC
set.seed(264)
id <- sample(seq(10000), n_samples)

lmic_af_op <-
  qs::qread("../../ESTIMATES_MIX_ROTA_V4/sim_all_outpat_lmic.qs")
lmic_af_op <- lmic_af_op$STEC
lmic_af_op$SAMPLES <- lapply(lmic_af_op$PROP, function(x) x[id])
lmic_af_op$AGE <-
  factor(lmic_af_op$AGE,
         c("Age below 5", "Age above or equal 5"),
         c("ch", "ad"))
summarize(subset(lmic_af_op, AGE == "ch"))
summarize(subset(lmic_af_op, AGE == "ad"))

lmic_af_ip <-
  qs::qread("../../ESTIMATES_MIX_ROTA_V4/sim_all_inpat_lmic.qs")
lmic_af_ip <- lmic_af_ip$STEC
lmic_af_ip$SAMPLES <- lapply(lmic_af_ip$PROP, function(x) x[id])
lmic_af_ip$AGE <-
  factor(lmic_af_ip$AGE,
         c("Age below 5", "Age above or equal 5"),
         c("ch", "ad"))
summarize(subset(lmic_af_ip, AGE == "ch"))
summarize(subset(lmic_af_ip, AGE == "ad"))

## YLD settings / diarrhea
dur_ch <- dalymod:::sim_pert(n_samples, 0.013835616, 0.011780822, 0.023013699)
mean_ci(dur_ch)
dur_ad <- 0.007671233
dsw <- dalymod:::sim_mean_ci(n_samples, 0.098, 0.074, 0.125, "PROB")
mean_ci(dsw)

## YLD settings / O157-HUS-ESRD
sub2 <- c("AMRA", "AMRB", "EURA", "EURB", "EURC", "WPRA")
p_O157 <-
rbind(
  data.frame(
    ISO3 = subset(FERG2:::countries, SUB2 %in% sub2)$ISO3,
    PROB = 0.311),
  data.frame(
    ISO3 = subset(FERG2:::countries, !(SUB2 %in% sub2))$ISO3,
    PROB = 0.100))

p_hus <- # estimate P(HUS|STEC) ~ type, age
  dalymod:::get_samples(
    n_samples,
    "../../../edtf-stec/ESTIMATES/edtf-stec/sim_all_20250924.rds",
    transformation = "logit",
    denominator = 1)
knitr::kable(unique(t(sapply(p_hus$SAMPLES, mean_ci))))

p_hus_esrd <- # estimate P(ESRD|HUS,STEC) ~ REG
  dalymod:::get_samples(
    n_samples,
    "../../../edtf-stec/ESTIMATES/edtf-stechus/sim_all_20250221.rds",
    transformation = "logit",
    denominator = 1)
knitr::kable(unique(t(sapply(p_hus_esrd$SAMPLES, mean_ci))))

cfr_hus_ch <- 0.02
cfr_hus_ad <- 0.04

# note - these are for HIC only
mult_p_hus_O157 <- 1.8 / 32
mult_p_hus_nO157 <- 1.8 / 106.8

dur_hus <-
  dalymod:::sim_pert(n_samples, 0.076712329, 0.038356164, 0.115068493)
mean_ci(dur_hus)
dsw_hus <-
  dalymod:::sim_mean_ci(n_samples, 0.247, 0.164, 0.348, "PROB")
mean_ci(dsw_hus)
dsw_esrd <-
  dalymod:::sim_mean_ci(n_samples, 0.571, 0.398, 0.725, "PROB")
mean_ci(dsw_esrd)

cfr_esrd_node <-
  dalymod:::import_node(
    "../disease_models/ferg2-daly-stec.xlsx",
    "PROB_O157_HUS_ERSD_CFR#<5")
cfr_esrd_node$val$dnmn <- 1
cfr_esrd <-
  dalymod:::pre_sample_node(cfr_esrd_node, n_samples)
str(head(cfr_esrd$val$samp))

# hus/esrd age groups
p_ch_hus <- 0.407
p_ad_hus <- c(0.185, 0.259, 0.148)

lle_ch_hus <- lle_ad_hus <- get_weighted_lle("<5")
lle2_hus <- get_weighted_lle("5-14")
lle3_hus <- get_weighted_lle("15-54")
lle4_hus <- get_weighted_lle("55+")
lle_ad_hus$LLE <-
  apply(cbind(lle2_hus$LLE, lle3_hus$LLE, lle4_hus$LLE),
        1, weighted.mean, p_ad_hus)

rle_ch_hus <- rle_ad_hus <- get_weighted_rle("<5")
rle2_hus <- get_weighted_rle("5-14")
rle3_hus <- get_weighted_rle("15-54")
rle4_hus <- get_weighted_rle("55+")
rle_ad_hus$RLE <-
  apply(cbind(rle2_hus$RLE, rle3_hus$RLE, rle4_hus$RLE),
        1, weighted.mean, p_ad_hus)

## get source attribution estimates
sa <- get_sa(hazard_sa, n_samples)
knitr::kable(t(sapply(sa, mean_ci)))

## GENERIC SCRIPT ----

## calculate DALYs per year
out <- NULL

for (year in as.character(all_yrs)) {
  
  # overall list - all nodes
  dalycalc <- list(diarrhea_inc[[year]], diarrhea_mrt[[year]])
  class(dalycalc) <- "dalycalc"
  dalycalc <- 
    dalycalc_aggregate_agesex(
      dalycalc_aggregate_nodes(dalycalc), age_agg, sex_agg)
  # str(dalycalc)
  
  # adjust measures / HIC
  for (country in countries_hic) {

    p_ch_country <-
      subset(p_ch_hic, COUNTRY == country & YEAR == year)$PROB
    
    # INC_NR
    inc_all <-
      subset(hic_inc, COUNTRY == country & YEAR == year)$SAMPLES[[1]] *
      (dalycalc[[country]][[1]]$POP + dalycalc[[country]][[2]]$POP)
    dalycalc[[country]][[1]]$INC_NR <-
      inc_all * p_ch_country
    dalycalc[[country]][[2]]$INC_NR <-
      inc_all * (1 - p_ch_country)
    
    # MRT_NR
    mrt_all <-
      subset(hic_mrt, COUNTRY == country & YEAR == year)$SAMPLES[[1]] *
      (dalycalc[[country]][[1]]$POP + dalycalc[[country]][[2]]$POP)
    dalycalc[[country]][[1]]$MRT_NR <-
      mrt_all * p_ch_country
    dalycalc[[country]][[2]]$MRT_NR <-
      mrt_all * (1 - p_ch_country)
    
    # YLD_NR
    dalycalc[[country]][[1]]$YLD_NR <-
      dalycalc[[country]][[1]]$INC_NR * dur_ch * dsw
    dalycalc[[country]][[2]]$YLD_NR <-
      dalycalc[[country]][[2]]$INC_NR * dur_ad * dsw 
    
    # YLL_NR
    dalycalc[[country]][[1]]$YLL_NR <-
      dalycalc[[country]][[1]]$MRT_NR *
      subset(rle_ch_hic_w, ISO3 == country & YEAR == year)$RLE
    dalycalc[[country]][[2]]$YLL_NR <-
      dalycalc[[country]][[2]]$MRT_NR *
      subset(rle_ad_hic_w, ISO3 == country & YEAR == year)$RLE
    
    # # DALY_NR
    # dalycalc[[country]][[1]]$DALY_NR <-
    #   dalycalc[[country]][[1]]$YLD_NR +
    #   dalycalc[[country]][[1]]$YLL_NR
    # dalycalc[[country]][[2]]$DALY_NR <-
    #   dalycalc[[country]][[2]]$YLD_NR +
    #   dalycalc[[country]][[2]]$YLL_NR
  }
  
  # inc_hic_ch <-
  #   1e5 * rowSums(sapply(countries_hic, function(x) dalycalc[[x]][[1]]$INC_NR)) /
  #   sum(sapply(countries_hic, function(x) dalycalc[[x]][[1]]$POP))
  # mean_ci(inc_hic_ch)  
  # inc_hic_ad <-
  #   1e5 * rowSums(sapply(countries_hic, function(x) dalycalc[[x]][[2]]$INC_NR)) /
  #   sum(sapply(countries_hic, function(x) dalycalc[[x]][[2]]$POP))
  # mean_ci(inc_hic_ad) 
  # inc_hic_all <-
  #   1e5 *
  #   (rowSums(sapply(countries_hic, function(x) dalycalc[[x]][[1]]$INC_NR)) +
  #      rowSums(sapply(countries_hic, function(x) dalycalc[[x]][[2]]$INC_NR))) /
  #   (sum(sapply(countries_hic, function(x) dalycalc[[x]][[1]]$POP)) +
  #      sum(sapply(countries_hic, function(x) dalycalc[[x]][[2]]$POP)))
  # mean_ci(inc_hic_all)
  
  # mrt_hic_ch <-
  #   1e5 * rowSums(sapply(countries_hic, function(x) dalycalc[[x]][[1]]$MRT_NR)) /
  #   sum(sapply(countries_hic, function(x) dalycalc[[x]][[1]]$POP))
  # mean_ci(mrt_hic_ch)  
  # mrt_hic_ad <-
  #   1e5 * rowSums(sapply(countries_hic, function(x) dalycalc[[x]][[2]]$MRT_NR)) /
  #   sum(sapply(countries_hic, function(x) dalycalc[[x]][[2]]$POP))
  # mean_ci(mrt_hic_ad) 
  # mrt_hic_all <-
  #   1e5 *
  #   (rowSums(sapply(countries_hic, function(x) dalycalc[[x]][[1]]$MRT_NR)) +
  #      rowSums(sapply(countries_hic, function(x) dalycalc[[x]][[2]]$MRT_NR))) /
  #   (sum(sapply(countries_hic, function(x) dalycalc[[x]][[1]]$POP)) +
  #      sum(sapply(countries_hic, function(x) dalycalc[[x]][[2]]$POP)))
  # mean_ci(mrt_hic_all)
  
  # adjust measures / LMIC
  for (country in countries_lmic) {

    # INC_NR
    # dalycalc[[country]][[1]]$INC_NR <-
    #   inc_hic_ch * dalycalc[[country]][[1]]$POP / 1e5
    # dalycalc[[country]][[2]]$INC_NR <-
    #   inc_hic_ad * dalycalc[[country]][[2]]$POP / 1e5
    # dalycalc[[country]][[1]]$INC_NR <- rep(0, n_samples)
    # dalycalc[[country]][[2]]$INC_NR <- rep(0, n_samples)
    dalycalc[[country]][[1]]$INC_NR <-
      dalycalc[[country]][[1]]$INC_NR *
        subset(lmic_af_op,
               AGE == "ch" & COUNTRY == country & YEAR == year)$SAMPLES[[1]]
    dalycalc[[country]][[2]]$INC_NR <-
      dalycalc[[country]][[2]]$INC_NR *
        subset(lmic_af_op,
               AGE == "ad" & COUNTRY == country & YEAR == year)$SAMPLES[[1]]
    
    # MRT_NR
    # dalycalc[[country]][[1]]$MRT_NR <-
    #   mrt_hic_ch * dalycalc[[country]][[1]]$POP / 1e5
    # dalycalc[[country]][[2]]$MRT_NR <-
    #   mrt_hic_ad * dalycalc[[country]][[2]]$POP / 1e5
    # dalycalc[[country]][[1]]$MRT_NR <- rep(0, n_samples)
    # dalycalc[[country]][[2]]$MRT_NR <- rep(0, n_samples)
    dalycalc[[country]][[1]]$MRT_NR <-
      dalycalc[[country]][[1]]$MRT_NR *
        subset(lmic_af_ip,
               AGE == "ch" & COUNTRY == country & YEAR == year)$SAMPLES[[1]]
    dalycalc[[country]][[2]]$MRT_NR <-
      dalycalc[[country]][[2]]$MRT_NR *
        subset(lmic_af_ip,
               AGE == "ad" & COUNTRY == country & YEAR == year)$SAMPLES[[1]]
    
    # YLD_NR
    dalycalc[[country]][[1]]$YLD_NR <-
      dalycalc[[country]][[1]]$INC_NR * dur_ch * dsw
    dalycalc[[country]][[2]]$YLD_NR <-
      dalycalc[[country]][[2]]$INC_NR * dur_ad * dsw 
    
    # YLL_NR
    # dalycalc[[country]][[1]]$YLL_NR <-
    #   dalycalc[[country]][[1]]$MRT_NR *
    #   subset(rle_ch_lmic_w, ISO3 == country & YEAR == year)$RLE
    # dalycalc[[country]][[2]]$YLL_NR <-
    #   dalycalc[[country]][[2]]$MRT_NR *
    #   subset(rle_ch_lmic_w, ISO3 == country & YEAR == year)$RLE
    # dalycalc[[country]][[1]]$YLL_NR <- rep(0, n_samples)
    # dalycalc[[country]][[2]]$YLL_NR <- rep(0, n_samples)
    dalycalc[[country]][[1]]$YLL_NR <-
      dalycalc[[country]][[1]]$YLL_NR *
        subset(lmic_af_ip,
               AGE == "ch" & COUNTRY == country & YEAR == year)$SAMPLES[[1]]
    dalycalc[[country]][[2]]$YLL_NR <-
      dalycalc[[country]][[2]]$YLL_NR *
        subset(lmic_af_ip,
               AGE == "ad" & COUNTRY == country & YEAR == year)$SAMPLES[[1]]
    
    # # DALY_NR
    # dalycalc[[country]][[1]]$DALY_NR <-
    #   dalycalc[[country]][[1]]$YLD_NR +
    #   dalycalc[[country]][[1]]$YLL_NR
    # dalycalc[[country]][[2]]$DALY_NR <-
    #   dalycalc[[country]][[2]]$YLD_NR +
    #   dalycalc[[country]][[2]]$YLL_NR
  }
  
  # add HUS/ESRD burden
  for (country in names(dalycalc)) {
    
    # HUS_INC_NR
    inc_hus_ch_O157 <-
      dalycalc[[country]][[1]]$INC_NR *
        subset(p_O157, ISO3 == country)$PROB *
        p_hus$SAMPLES$V_O157_Age19 *
        mult_p_hus_O157
    inc_hus_ad_O157 <-
      dalycalc[[country]][[2]]$INC_NR *
        subset(p_O157, ISO3 == country)$PROB *
        p_hus$SAMPLES$V_O157_Age125 *
        mult_p_hus_O157
    inc_hus_ch_nO157 <-
      dalycalc[[country]][[1]]$INC_NR *
        (1 - subset(p_O157, ISO3 == country)$PROB) *
        p_hus$SAMPLES$V_nonO157_Age19 *
        mult_p_hus_nO157
    inc_hus_ad_nO157 <-
      dalycalc[[country]][[2]]$INC_NR *
        (1 - subset(p_O157, ISO3 == country)$PROB) *
        p_hus$SAMPLES$V_nonO157_Age125 *
        mult_p_hus_nO157
    inc_hus_ch <- inc_hus_ch_O157 + inc_hus_ch_nO157
    inc_hus_ad <- inc_hus_ad_O157 + inc_hus_ad_nO157
    
    # HUS_MRT_NR
    mrt_hus_ch <- inc_hus_ch * cfr_hus_ch
    mrt_hus_ad <- inc_hus_ad * cfr_hus_ad
    
    # ESRD_INC_NR
    inc_esrd_ch <-
      inc_hus_ch * subset(p_hus_esrd, COUNTRY == country)$SAMPLES[[1]]
    inc_esrd_ad <-
      inc_hus_ad * subset(p_hus_esrd, COUNTRY == country)$SAMPLES[[1]]
    
    inc_esrd_ch_surv <-
      inc_esrd_ch *
        (1 - subset(cfr_esrd$val$samp, COUNTRY == country)$SAMPLES[[1]])
    inc_esrd_ad_surv <-
      inc_esrd_ad *
        (1 - subset(cfr_esrd$val$samp, COUNTRY == country)$SAMPLES[[1]])
    
    # ESRD_MRT_NR
    mrt_esrd_ch <-
      inc_esrd_ch *
        subset(cfr_esrd$val$samp, COUNTRY == country)$SAMPLES[[1]]
    mrt_esrd_ad <-
      inc_esrd_ad *
        subset(cfr_esrd$val$samp, COUNTRY == country)$SAMPLES[[1]]
        
    # YLD_NR
    dalycalc[[country]][[1]]$YLD_NR <-
      dalycalc[[country]][[1]]$YLD_NR +
        inc_hus_ch * dsw_hus * dur_hus +
        inc_esrd_ch_surv * dsw_esrd * subset(lle_ch_hus, ISO3 == country & YEAR == year)$LLE
    dalycalc[[country]][[2]]$YLD_NR <-
      dalycalc[[country]][[2]]$YLD_NR +
        inc_hus_ad * dsw_hus * dur_hus +
        inc_esrd_ad_surv * dsw_esrd * subset(lle_ad_hus, ISO3 == country & YEAR == year)$LLE
    
    # MRT_NR
    dalycalc[[country]][[1]]$MRT_NR <-
      dalycalc[[country]][[1]]$MRT_NR + mrt_hus_ch + mrt_esrd_ch
    
    dalycalc[[country]][[2]]$MRT_NR <-
      dalycalc[[country]][[2]]$MRT_NR + mrt_hus_ad + mrt_esrd_ad
    
    # YLL_NR
    dalycalc[[country]][[1]]$YLL_NR <-
      dalycalc[[country]][[1]]$YLL_NR +
        mrt_hus_ch * subset(rle_ch_hus, ISO3 == country & YEAR == year)$RLE +
        mrt_esrd_ch * subset(rle_ch_hus, ISO3 == country & YEAR == year)$RLE
      
    dalycalc[[country]][[2]]$YLL_NR <-
      dalycalc[[country]][[2]]$YLL_NR +
        mrt_hus_ad * subset(rle_ad_hus, ISO3 == country & YEAR == year)$RLE +
        mrt_esrd_ad * subset(rle_ad_hus, ISO3 == country & YEAR == year)$RLE
    
    # DALY_NR
    dalycalc[[country]][[1]]$DALY_NR <-
      dalycalc[[country]][[1]]$YLD_NR +
      dalycalc[[country]][[1]]$YLL_NR
    dalycalc[[country]][[2]]$DALY_NR <-
      dalycalc[[country]][[2]]$YLD_NR +
      dalycalc[[country]][[2]]$YLL_NR
  }
  
  ## check
  length(dalycalc)           # countries
  length(dalycalc[[1]])      # age-sex
  length(dalycalc[[1]][[1]]) # values

  ## summaries
  dalycalc_summary(
    dalycalc_aggregate_country(dalycalc, glb))
  
  dalycalc_summary(
    dalycalc_aggregate_country(
      dalycalc_aggregate_agesex(dalycalc, age_agg_all), glb))
  
  ## apply source attribution
  dalycalc_food <- dalycalc
  for (i in seq_along(dalycalc_food)) { # per country
    reg2 <- subset(FERG2:::countries, ISO3 == names(dalycalc)[i])$SUB2
    dalycalc_food[i] <-
      dalycalc_mult(dalycalc_food[i], sa[[reg2]])
  }
  
  ## create summaries
  out_yr_all <- dalyout_agg(dalycalc)
  out_yr_all$YEAR <- year
  out_yr_all$SOURCE <- "ALL"
  out_yr_all <- out_yr_all[, c(11, 1:2, 12, 3:10)]
  # str(out_yr_all)
  
  out_yr_food <- dalyout_agg(dalycalc_food)
  out_yr_food$YEAR <- year
  out_yr_food$SOURCE <- "FOOD"
  out_yr_food <- out_yr_food[, c(11, 1:2, 12, 3:10)]
  # str(out_yr_food)
  
  out_yr <- rbind(out_yr_all, out_yr_food)
  out <- rbind(out, out_yr)
  # str(out)

  ## save simulations - aggregates age-sex
  saveRDS(
    dalycalc,
    sprintf("SIM/%s-dalycalc-%s.rds", hazard, year))
  saveRDS(
    dalycalc_food,
    sprintf("SIM_FOOD/%s-dalycalc-food-%s.rds", hazard, year))
}

out$MEASURE <- factor(out$MEASURE)
out$METRIC <- factor(out$METRIC)
out$COUNTRY <- factor(out$COUNTRY)
saveRDS(out, sprintf("%s-dalyout-%s.rds", hazard, today()))
writexl::write_xlsx(out, sprintf("%s-dalyout-%s.xlsx", hazard, today()))
# str(out)

DT::datatable(
  rownames = FALSE,
  subset(out, AGE == "ALL AGES" & COUNTRY == "GLOBAL" & METRIC == "NR")) |>
  DT::formatRound(columns = 8:12, digits = c(0, rep(3, 4)))

##bd::render_today("edtf-diarrhea-stec-daly.R")