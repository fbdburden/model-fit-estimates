#' ---
#' params:
#'   gbs: TRUE
#' ---
#' 

### FERG2/EDTF/DIARRHEA/CAMPYLOBACTER/DALY

## SETTINGS ----

## required packages
library(dalymod)
library(FERG2)

## source generic FERG2 settings
source("ferg-settings.R")

## hazard specific settings
hazard <- "edtf-diarrhea-campy"
hazard_sa <- "CAMPY"
n_samples <- 500

## DATA ----

## source diarrhea estimates
source("../diarrhea-data.R")

## HIC
hic_inc <-
  dalymod:::get_samples(
    n_samples,
    "../../ESTIMATES_HIC/01-campy/sim_all_inc_reg_yr_20251011.rds",
    transformation = "log",
    denominator = 1e5)
summarize(hic_inc)

hic_mrt <-
  dalymod:::get_samples(
    n_samples,
    "../../ESTIMATES_HIC/01-campy/sim_all_mrt_reg_ct_20251208.rds",
    transformation = "log",
    denominator = 1e5)
summarize(hic_mrt)

## HIC / age distribution

# define age distribution
hic_age <-
  list(`2000` = c("<5" = 0.132, "5-9" = 0.051, "10-19" = 0.082, "20-64" = 0.647, "65+" = 0.087),
       `2010` = c("<5" = 0.123, "5-9" = 0.052, "10-19" = 0.099, "20-64" = 0.605, "65+" = 0.121),
       `2021` = c("<5" = 0.111, "5-9" = 0.042, "10-19" = 0.069, "20-64" = 0.569, "65+" = 0.209))
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
lmic_af_op <- lmic_af_op$CAMP
lmic_af_op$SAMPLES <- lapply(lmic_af_op$PROP, function(x) x[id])
lmic_af_op$AGE <-
  factor(lmic_af_op$AGE,
         c("Age below 5", "Age above or equal 5"),
         c("ch", "ad"))
summarize(subset(lmic_af_op, AGE == "ch"))
summarize(subset(lmic_af_op, AGE == "ad"))

lmic_af_ip <-
  qs::qread("../../ESTIMATES_MIX_ROTA_V4/sim_all_inpat_lmic.qs")
lmic_af_ip <- lmic_af_ip$CAMP
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
dsw <- dalymod:::sim_mean_ci(n_samples, 0.106, 0.081,	0.134, "PROB")
mean_ci(dsw)

## YLD settings / GBS

## load prevalence files
gbs_prev_all <- readRDS("../gbs-prev-ihme.rds")
gbs_prev_cvd <- readRDS("../gbs-covid-prev-ihme.rds")

## check order of years/countries
all(names(gbs_prev_all) == names(gbs_prev_cvd))
all(names(gbs_prev_all) == names(diarrhea_mrt))
all(names(gbs_prev_all[[1]]) == names(gbs_prev_cvd[[1]]))
all(names(gbs_prev_all[[1]]) == names(diarrhea_mrt[[1]]))
all.equal(
  sapply(gbs_prev_all[[1]][[1]], function(x) x$POP),
  sapply(diarrhea_mrt[[1]][[1]], function(x) x$POP))

## disease model settings
af_gbs_campy <- dalymod:::sim_pert(n_samples, 0.31, 0.11, 0.45)
mean_ci(af_gbs_campy)
dsw_gbs <- dalymod:::sim_mean_ci(n_samples, 0.296, 0.198, 0.414, "PROB")
mean_ci(dsw_gbs)
cfr_gbs_campy <- dalymod:::sim_pert(n_samples, 0.041, 0.024, 0.060)
mean_ci(cfr_gbs_campy)

# calculate weighted LLE/RLE

# age dist
#     <5	5-14	15-24	25-64	 65+
#   0.11	0.08	 0.10	 0.57	0.14

gbs_ch <- 0.11
gbs_ad <- c(0.08, 0.10, 0.57, 0.14)
lle_ch <- lle_ad <- get_weighted_lle("<5")
lle2 <- get_weighted_lle("5-14")
lle3 <- get_weighted_lle("15-24")
lle4 <- get_weighted_lle("25-64")
lle5 <- get_weighted_lle("65+")
lle_ad$LLE <-
  apply(cbind(lle2$LLE, lle3$LLE, lle4$LLE, lle5$LLE),
        1, weighted.mean, gbs_ad)

rle_ch <- rle_ad <- get_weighted_rle("<5")
rle2 <- get_weighted_rle("5-14")
rle3 <- get_weighted_rle("15-24")
rle4 <- get_weighted_rle("25-64")
rle5 <- get_weighted_rle("65+")
rle_ad$RLE <-
  apply(cbind(rle2$RLE, rle3$RLE, rle4$RLE, rle5$RLE),
        1, weighted.mean, gbs_ad)

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
    # dalycalc[[country]][[1]]$YLL_NR <- 0
    # dalycalc[[country]][[2]]$YLL_NR <- 0
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
  
  # add GBS burden
  
  if (params$gbs) {
    
    # compile across ages
    # note: 'INC_NR' is actually 'PREV_NR'
    dalycalc_gbs_all <- gbs_prev_all[[year]]
    dalycalc_gbs_cvd <- gbs_prev_cvd[[year]]
    class(dalycalc_gbs_all) <- class(dalycalc_gbs_cvd) <- "dalycalc_agg"
    dalycalc_gbs_all <-
      dalycalc_aggregate_agesex(dalycalc_gbs_all, age_agg, sex_agg)
    dalycalc_gbs_cvd <-
      dalycalc_aggregate_agesex(dalycalc_gbs_cvd, age_agg, sex_agg)
    # str(dalycalc_gbs_all)
    # str(dalycalc_gbs_cvd)
    
    for (country in names(dalycalc)) {
      
      # MRT_NR
      dalycalc[[country]][[1]]$MRT_NR <-
        dalycalc[[country]][[1]]$MRT_NR +
        (((dalycalc_gbs_all[[country]][[1]]$INC_NR -
             dalycalc_gbs_cvd[[country]][[1]]$INC_NR) /
            subset(lle_ch, ISO3 == country & YEAR == year)$LLE) * 
           af_gbs_campy *
           cfr_gbs_campy)
      dalycalc[[country]][[2]]$MRT_NR <-
        dalycalc[[country]][[2]]$MRT_NR +
        (((dalycalc_gbs_all[[country]][[2]]$INC_NR -
             dalycalc_gbs_cvd[[country]][[2]]$INC_NR) /
            subset(lle_ad, ISO3 == country & YEAR == year)$LLE) * 
           af_gbs_campy *
           cfr_gbs_campy)
      
      # YLD_NR
      dalycalc[[country]][[1]]$YLD_NR <-
        dalycalc[[country]][[1]]$YLD_NR +
        ((dalycalc_gbs_all[[country]][[1]]$INC_NR -
            dalycalc_gbs_cvd[[country]][[1]]$INC_NR) * 
           af_gbs_campy *
           dsw_gbs)
      dalycalc[[country]][[2]]$YLD_NR <-
        dalycalc[[country]][[2]]$YLD_NR +
        ((dalycalc_gbs_all[[country]][[2]]$INC_NR -
            dalycalc_gbs_cvd[[country]][[2]]$INC_NR) * 
           af_gbs_campy *
           dsw_gbs)
      
      # YLL_NR
      dalycalc[[country]][[1]]$YLL_NR <-
        dalycalc[[country]][[1]]$YLL_NR +
        (((dalycalc_gbs_all[[country]][[1]]$INC_NR -
             dalycalc_gbs_cvd[[country]][[1]]$INC_NR) /
            subset(lle_ch, ISO3 == country & YEAR == year)$LLE) * 
           af_gbs_campy *
           cfr_gbs_campy *
           subset(rle_ch, ISO3 == country & YEAR == year)$RLE)
      dalycalc[[country]][[2]]$YLL_NR <-
        dalycalc[[country]][[2]]$YLL_NR +
        (((dalycalc_gbs_all[[country]][[2]]$INC_NR -
             dalycalc_gbs_cvd[[country]][[2]]$INC_NR) /
            subset(lle_ad, ISO3 == country & YEAR == year)$LLE) * 
           af_gbs_campy *
           cfr_gbs_campy *
           subset(rle_ad, ISO3 == country & YEAR == year)$RLE)
      
      # DALY_NR
      dalycalc[[country]][[1]]$DALY_NR <-
        dalycalc[[country]][[1]]$YLD_NR +
        dalycalc[[country]][[1]]$YLL_NR
      dalycalc[[country]][[2]]$DALY_NR <-
        dalycalc[[country]][[2]]$YLD_NR +
        dalycalc[[country]][[2]]$YLL_NR
    }
    
  } else {
    for (country in names(dalycalc)) {
      # DALY_NR
      dalycalc[[country]][[1]]$DALY_NR <-
        dalycalc[[country]][[1]]$YLD_NR +
        dalycalc[[country]][[1]]$YLL_NR
      dalycalc[[country]][[2]]$DALY_NR <-
        dalycalc[[country]][[2]]$YLD_NR +
        dalycalc[[country]][[2]]$YLL_NR
    }
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
    sprintf("%sSIM/%s-dalycalc-%s.rds",
            ifelse(params$gbs, "", "NOGBS/"), hazard, year))
  saveRDS(
    dalycalc_food,
    sprintf("%sSIM_FOOD/%s-dalycalc-food-%s.rds",
            ifelse(params$gbs, "", "NOGBS/"), hazard, year))
}

out$MEASURE <- factor(out$MEASURE)
out$METRIC <- factor(out$METRIC)
out$COUNTRY <- factor(out$COUNTRY)
saveRDS(
  out,
  sprintf("%s%s-dalyout-%s.rds",
          ifelse(params$gbs, "", "NOGBS/"), hazard, today()))
writexl::write_xlsx(
  out,
  sprintf("%s%s-dalyout-%s.xlsx",
          ifelse(params$gbs, "", "NOGBS/"), hazard, today()))
# str(out)

DT::datatable(
  rownames = FALSE,
  subset(out, AGE == "ALL AGES" & COUNTRY == "GLOBAL" & METRIC == "NR")) |>
  DT::formatRound(columns = 8:12, digits = c(0, rep(3, 4)))

##bd::render_today("edtf-diarrhea-campy-daly.R", params = list(gbs = TRUE))
##bd::render_today("edtf-diarrhea-campy-daly.R", params = list(gbs = FALSE))