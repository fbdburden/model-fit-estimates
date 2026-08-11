### FERG2/EDTF/DIARRHEA/ENTAMOEBA/DALY

## SETTINGS ----

## required packages
library(dalymod)
library(FERG2)

## source generic FERG2 settings
source("ferg-settings.R")

## hazard specific settings
hazard <- "edtf-diarrhea-enta"
hazard_sa <- "ENTA"
n_samples <- 500

## DATA ----

## source diarrhea estimates
source("../diarrhea-data.R")

## HIC / INC only - hard-coded
set.seed(264)
hic_inc <- dalymod:::sim_gamma(n_samples, 1, 2340)
mean_ci(hic_inc) # 1/2340 = 0.0004

## HIC / age distribution

# define age distribution
hic_age <-
  list(`2000` = c("<20" = 0.004, "20-29" = 0.396, "30-39" = 0.329, "40-49" = 0.118, "50-59" = 0.067, "60+" = 0.086),
       `2010` = c("<20" = 0.004, "20-29" = 0.396, "30-39" = 0.329, "40-49" = 0.118, "50-59" = 0.067, "60+" = 0.086),
       `2021` = c("<20" = 0.016, "20-29" = 0.300, "30-39" = 0.424, "40-49" = 0.128, "50-59" = 0.048, "60+" = 0.084))
hic_age <-
  lapply(hic_age, function(x) x/sum(x))

# proportion <5, by country and year
p_ch_hic <- lapply(countries_hic, get_prob_ch_country, x = hic_age)
p_ch_hic <- do.call("rbind", p_ch_hic)
head(p_ch_hic)

# average RLE by <5 and 5+, by country and year
# given very small proportion '<5', we do not correct '<20' RLE
rle_ch_hic_w <- get_weighted_rle("<5")

rle_ad_hic <-
  rbind(cbind(AGE = "<20", get_weighted_rle("<20")),
        cbind(AGE = "20-29", get_weighted_rle("20-29")),
        cbind(AGE = "30-39", get_weighted_rle("30-39")),
        cbind(AGE = "40-49", get_weighted_rle("40-49")),
        cbind(AGE = "50-59", get_weighted_rle("50-59")),
        cbind(AGE = "60+", get_weighted_rle("60+")))

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

## LMIC
set.seed(264)
id <- sample(seq(10000), n_samples)

lmic_af_op <-
  qs::qread("../../ESTIMATES_MIX_ROTA_V5/sim_all_outpat_lmic.qs")
lmic_af_op <- lmic_af_op$ENTA
lmic_af_op$SAMPLES <- lapply(lmic_af_op$PROP, function(x) x[id])
lmic_af_op$AGE <-
  factor(lmic_af_op$AGE,
         c("Age below 5", "Age above or equal 5"),
         c("ch", "ad"))
summarize(subset(lmic_af_op, AGE == "ch"))
summarize(subset(lmic_af_op, AGE == "ad"))

lmic_af_ip <-
  qs::qread("../../ESTIMATES_MIX_ROTA_V5/sim_all_inpat_lmic.qs")
lmic_af_ip <- lmic_af_ip$ENTA
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
dsw <- dalymod:::sim_mean_ci(n_samples, 0.085, 0.061, 0.113, "PROB")
mean_ci(dsw)

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
    # inc_all <-
    #   subset(hic_inc, COUNTRY == country & YEAR == year)$SAMPLES[[1]] *
    #   (dalycalc[[country]][[1]]$POP + dalycalc[[country]][[2]]$POP)
    inc_all <-
      hic_inc *
      (dalycalc[[country]][[1]]$POP + dalycalc[[country]][[2]]$POP)
    dalycalc[[country]][[1]]$INC_NR <-
      inc_all * p_ch_country
    dalycalc[[country]][[2]]$INC_NR <-
      inc_all * (1 - p_ch_country)
    
    # MRT_NR
    dalycalc[[country]][[1]]$MRT_NR <- rep(0, n_samples)
    dalycalc[[country]][[2]]$MRT_NR <- rep(0, n_samples)
    # mrt_all <-
    #   subset(hic_mrt, COUNTRY == country & YEAR == year)$SAMPLES[[1]] *
    #   (dalycalc[[country]][[1]]$POP + dalycalc[[country]][[2]]$POP)
    # dalycalc[[country]][[1]]$MRT_NR <-
    #   mrt_all * p_ch_country
    # dalycalc[[country]][[2]]$MRT_NR <-
    #   mrt_all * (1 - p_ch_country)
    
    # YLD_NR
    dalycalc[[country]][[1]]$YLD_NR <-
      dalycalc[[country]][[1]]$INC_NR * dur_ch * dsw
    dalycalc[[country]][[2]]$YLD_NR <-
      dalycalc[[country]][[2]]$INC_NR * dur_ad * dsw 
    
    # YLL_NR
    dalycalc[[country]][[1]]$YLL_NR <- rep(0, n_samples)
    dalycalc[[country]][[2]]$YLL_NR <- rep(0, n_samples)
    # dalycalc[[country]][[1]]$YLL_NR <-
    #   dalycalc[[country]][[1]]$MRT_NR *
    #   subset(rle_ch_hic_w, ISO3 == country & YEAR == year)$RLE
    # dalycalc[[country]][[2]]$YLL_NR <-
    #   dalycalc[[country]][[2]]$MRT_NR *
    #   subset(rle_ad_hic_w, ISO3 == country & YEAR == year)$RLE
    
    # DALY_NR
    dalycalc[[country]][[1]]$DALY_NR <-
      dalycalc[[country]][[1]]$YLD_NR +
      dalycalc[[country]][[1]]$YLL_NR
    dalycalc[[country]][[2]]$DALY_NR <-
      dalycalc[[country]][[2]]$YLD_NR +
      dalycalc[[country]][[2]]$YLL_NR
  }
  
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
    dalycalc[[country]][[1]]$YLL_NR <-
      dalycalc[[country]][[1]]$YLL_NR *
        subset(lmic_af_ip,
               AGE == "ch" & COUNTRY == country & YEAR == year)$SAMPLES[[1]]
    dalycalc[[country]][[2]]$YLL_NR <-
      dalycalc[[country]][[2]]$YLL_NR *
        subset(lmic_af_ip,
               AGE == "ad" & COUNTRY == country & YEAR == year)$SAMPLES[[1]]
    
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

##bd::render_today("edtf-diarrhea-enta-daly.R")