## get GAM thresholds

library(gam)
library(tidyverse)
library(sf)
library(tidylog)

## function to find value in curve
load("code/functions/root_interpolation_function.Rdata")

## function to find closest values to zero
load("code/functions/find_closest_values.RData")
## load DF 

load(file = "ignore/04_quantGams_smooths_predictions_updatedSites.RData")
head(DF)

## data
load(file = "final_data/01_bugs_algae_flow_joined_by_masterid.RData")
AllDataLongx <- AllData
head(AllDataLongx)

AllDataLongx %>% group_by(Flow.Metric.Name) %>% summarise(max(na.omit(deltah_final)), min(na.omit(deltah_final)))

# format for thresholds ----------------------------------------------------------

## full names for labels
labels <- read.csv("input_data/ffm_names.csv")
labels <- labels[1:24, ]
labels <- labels %>% rename(hydro.endpoints = Flow.Metric.Code)
labels[25, 1] <- "Peak Flow Magnitude (Q99, cfs)"
labels[25, 2] <- "Q99"
labels[25, 3] <- "Peak Flow Magnitude"
labels

DF

## FIX NAMES TO MATCH LABELS AND LIMITS 
dfx <- DF %>% 
  mutate(hydro.endpoints = case_when(Variable == "d_ds_mag_50" ~ "DS_Mag_50",            
                                     Variable == "d_fa_mag" ~ "FA_Mag",
                                     Variable == "d_peak_10" ~ "Peak_10",
                                     Variable == "d_peak_2" ~ "Peak_2",
                                     Variable == "d_peak_5" ~ "Peak_5",
                                     Variable == "d_sp_mag" ~ "SP_Mag",
                                     Variable == "d_wet_bfl_mag_10" ~ "Wet_BFL_Mag_10",
                                     Variable == "d_wet_bfl_mag_50" ~ "Wet_BFL_Mag_50", 
                                     Variable == "delta_q99" ~ "Q99")) %>%
  mutate(comb_code_type = paste0(Metric, "_", hydro.endpoints, "_", Smooths, "_", Quant)) %>%
  rename(HydroValue = deltah_final, BioValue = predictedVals)

dfx

all_data <- left_join(dfx, labels, by ="hydro.endpoints")

head(all_data)


# calculate thresholds ----------------------------------------------------

## create df
df <- as.data.frame(matrix(ncol=20))
colnames(df) <- c("metric", "ThresholdNAT_Lower", "ThresholdNAT_Upper", "ThresholdNATMed_Lower", "ThresholdNATMed_Upper",
                  "ThresholdNATLow_Lower", "ThresholdNATLow_Upper",
                  "ThresholdNATHigh_Lower", "ThresholdNATHigh_Upper", "ThresholdHB_Lower", "ThresholdHB_Upper",
                  "ThresholdSB0_Lower", "ThresholdSB0_Upper", "ThresholdSB2_Lower", "ThresholdSB2_Upper",
                  "n", "Index", "Hydro_endpoint","SmoothingFunc", "Quantile")
df

## define metrics
metrics <- unique(all_data$comb_code_type)
metrics

i=18
i
## loop through metrics
for(i in 1: length(metrics)) {
  
  ## define metric
  met <- metrics[i]
  
  ## filter by metric
  hydroxx <- all_data %>%
    filter(comb_code_type == met)
  
  head(hydroxx)
  
  ## define thresholds 
    
    ## get curves values at different probabilities, get closest value to zero for both -ve and +ve
    threshNAT <- find_closest_values(RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.79))
    threshNAT
    # threshNAT <- ifelse(length(threshNAT) == 0, NA, threshNAT)
    
    threshNATMed <- find_closest_values(RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.79)) ## sep ref threshold for overall tally
    threshNATMed
    threshNATLow <- find_closest_values(RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.63))


    threshNATHigh <- find_closest_values(RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.92))
    
    threshHB <- find_closest_values(RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.67))
    
    threshSB0 <- find_closest_values(RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.78))
    # threshSB0 <- ifelse(length(threshSB0) == 0, NA, threshSB0)
    
    threshSB2 <- find_closest_values(RootLinearInterpolant(hydroxx$HydroValue, hydroxx$BioValue, 0.75))
    # threshSB2 <- ifelse(length(threshSB2) == 0, NA, threshSB2)

  ## add info to df
  df[i, 1] <- met
  df[i, 2] <- threshNAT$closest_negative
  df[i, 3] <- threshNAT$closest_positive
  df[i, 4] <- threshNATMed$closest_negative
  df[i, 5] <- threshNATMed$closest_positive
  df[i, 6] <- threshNATLow$closest_negative
  df[i, 7] <- threshNATLow$closest_positive
  df[i, 8] <- threshNATHigh$closest_negative
  df[i, 9] <- threshNATHigh$closest_positive
  df[i, 10] <- threshHB$closest_negative
  df[i, 11] <- threshHB$closest_positive
  df[i, 12] <- threshSB0$closest_negative
  df[i, 13] <- threshSB0$closest_positive
  df[i, 14] <- threshSB2$closest_negative
  df[i, 15] <- threshSB2$closest_positive
  df[i, 16] <- length(hydroxx$BioValue)
  df[i ,17] <- hydroxx$Metric[1]
  df[i, 18] <- hydroxx$hydro.endpoints[1]
  df[i, 19] <- hydroxx$Smooths[1]
  df[i, 20] <- paste(hydroxx$Quant[1])

  
}


df ## NAs - where curve doesn't reach threshold
write.csv(df, "final_data/05_delta_thresholds_GAMs.csv")

# df <- read.csv("final_data/05_delta_thresholds_GAMs.csv")

head(df)
unique(df$Hydro_endpoint)

### metrics and smooths

# csci
## DS_Mag_50 - 6
## FA_Mag - 6
## Peak_2 - 3
## Peak_5 - 3
## Peak_10 - 3
## wet_bfl_mag_10 - 6
## wet_bfl_mag_10 - 6

sm6csci <- c("DS_Mag_50", "FA_Mag", "Peak_10", "Peak_2", "Peak_5", "Wet_BFL_Mag_50", "Wet_BFL_Mag_10")
sm6csci

## Q99 - 3
## Sp_mag - 3


sm3csci <- c( "SP_Mag",  "Q99")
sm3csci


unique(df$Hydro_endpoint)

## subset models/thresholds to only the ones we're using
##  quants 0.3, 0.9 & 0.5 for all metrics
dfsub <- df %>%
  # filter(Quantile == c(0.9, 0.5)) %>%
  mutate(ChosenSmth = case_when((Index == "csci" & Hydro_endpoint %in% sm6csci) ~ 6,
                               (Index == "csci" & Hydro_endpoint %in% sm3csci) ~ 3)) %>%
  mutate(ModelstoUse = ifelse(SmoothingFunc == ChosenSmth, "Yes", "No")) %>%
  filter(ModelstoUse == "Yes") %>%
  filter(Quantile %in% c(0.90, 0.50))

unique(dfsub$Hydro_endpoint)
  
write.csv(dfsub, "final_data/05_chosen_models_thresholds.csv")
dfsub


# format to apply data to categories----------------------------------------------------

df <- read.csv("final_data/05_chosen_models_thresholds.csv")
# add best observed thresholds
## apply one ffm - raw data
## how many modified channels are within thresholds

head(df) ## thresholds

## format 
limits <- df %>%
  dplyr::select( -n, -X) %>%
  mutate(metric = paste0(Index, "_", Hydro_endpoint)) %>%  # make code 
  pivot_longer(ThresholdNAT_Lower:ThresholdSB2_Upper, names_to = "Threshold", values_to = "DeltaH") %>% # make threshold names longer
  separate(Threshold, into=c("Threshold", "Type")) %>%
  mutate(Threshold = gsub("Threshold", "", Threshold)) %>%
  mutate(BioThresh = case_when((Index == "csci" & Threshold == "NAT") ~ 0.79,
                               (Index == "csci" & Threshold == "NATMed") ~ 0.79,
                               (Index == "csci" & Threshold == "NATLow") ~ 0.63,
                               (Index == "csci" & Threshold == "NATHigh") ~ 0.92,
                               (Index == "csci" & Threshold == "HB") ~ 0.67,
                               (Index == "csci" & Threshold == "SB0") ~ 0.78,
                               (Index == "csci" & Threshold == "SB2") ~ 0.75))

limits

## make wider with type - lower/upper
limits <- limits %>%
  pivot_wider(names_from = Type, values_from = DeltaH)  %>%
  mutate(MetricCurve = paste0(metric,"_", SmoothingFunc, "_", Quantile))

unique(limits$Threshold)

## organise limits by type of curve: AboveCurve, BelowCurve
# ## 0.5s should be below, and 0.9s above - there are some exceptions in 0.5 - we'll adjust manually below
AboveCurve <- c(0.9)
BelowCurve <- c(0.5)

## assign labels to know which curves are above or below thresholds
## some thresholds hit curve, so make sure it's NA too
limits <- limits %>%
  mutate(CurveStatus = ifelse(is.na(Upper) & Quantile %in% AboveCurve, "Above", "Passes")) %>%
  mutate(CurveStatus = ifelse(is.na(Lower) & Quantile  %in% BelowCurve, "Below", CurveStatus)) 

## add in minimum delta H values, where no limits can be found
## here: alteration isn't an issue as curve always within limits

MinDel <- min(na.omit(AllDataLongx$deltah_final)) ## -51409.35
MaxDel <- max(na.omit(AllDataLongx$deltah_final)) ## 8059.891

## if a curve is below the threshold, value should be zero - should result in "outside limits"
## Q99 and peak 5 have funky 0.3 curves, should have lower and upper limits of 0 too
limits2 <- limits %>%
  mutate(Lower2 = ifelse(is.na(Lower) & CurveStatus == "Below", 0, Lower),
         Upper2 = ifelse(is.na(Upper) & CurveStatus == "Below", 0, Upper)) #%>%
  # mutate(Upper2 = ifelse(Hydro_endpoint == "Q99" & Quantile == 0.3 & CurveStatus == "Below" & Upper2 > 0, 0, Upper2),
  #        Upper2 = ifelse(Hydro_endpoint == "Peak_5" & Quantile == 0.3 & CurveStatus == "Below"& Upper2 > 0, 0, Upper2))

## if a curve is above the threshold, max value should be added 
limits3 <- limits2 %>%
  mutate(Lower2 = ifelse(is.na(Lower2) & CurveStatus == "Above", MinDel, Lower2),
         Upper2 = ifelse(is.na(Upper2) & CurveStatus == "Above", MaxDel, Upper2)) 
names(limits3)

## exceptions to the above/below rule - 
## Peak 2, 0.3,  has a lower limit for 0.63, upper limit is also negative - change to 0
## Peak 5 0.5 has a lower limit for 0.63, but rest of curve does not hit - change to max
## Q99 shouldn't have any positive values, change to 0
## wet bfl 10 basically peaks at the 0.78 threshold, so upper value will be negative - change to 0

limits4 <- limits3 %>%
  mutate(Upper2 = ifelse(is.na(Upper2) & MetricCurve == "csci_Peak_5_6_0.5", MaxDel, Upper2)) %>%
  mutate(Upper2 = ifelse(is.na(Upper2) & Hydro_endpoint %in% c("Q99", "Peak_2"), 0, Upper2)) %>%
## all remaining NAs are 0.9 quantiles, so are above curve - change to min value
  mutate(Lower2 = ifelse(is.na(Lower2), MinDel, Lower2)) 
  # select(-X.1)

unique(limits4$Hydro_endpoint)

## save out
write.csv(limits4, "final_data/05_QGAM_FFM_Ranges_With_Adjustments.csv")

limits4 <- read.csv("final_data/05_QGAM_FFM_Ranges_With_Adjustments.csv")

# Add to predicted value --------------------------------------------------

names(AllDataLongx)

## take only csci and rename to match columns
AllDataLong2 <- AllDataLongx %>%
  select( -Class2) %>%
  filter(Metric %in% c("csci")) %>%
  drop_na(deltah_final) %>%
  rename(Hydro_endpoint = hydro.endpoints,
         Index = Metric,
         Threshold = channel_engineering_class,
         IndexValue = MetricValue)

head(AllDataLong2)

### sep df for overall thresholds - nat low/high

AllDataLongOverall <- AllDataLong2 %>%
  mutate(ThresholdNatLow = "NATLow", ThresholdNatMed = "NATMed", ThresholdNatHigh = "NATHigh") %>%
  select(-c(Threshold)) %>%
  pivot_longer(ThresholdNatLow:ThresholdNatHigh, names_to = "Check", values_to = "Threshold") %>%
  select(-Check) 

AllDataLongOverall
## join categorised and overall dfs together

AllDF <- bind_rows(AllDataLongOverall, AllDataLong2) %>% distinct()
class(AllDF)
names(AllDF)

## count site per class with FFM
tallyFFM <- AllDF %>%
  group_by(Threshold) %>%
  select(masterid, comid) %>%
  distinct() %>%
  drop_na(Threshold) %>%
  tally()

tallyFFM

## join limits to data
names(AllDF)
names(limits4)


allLims <- full_join(limits4, AllDF, by = c("Index", "Hydro_endpoint", "Threshold"), relationship = "many-to-many")
names(allLims)

unique(allLims$Hydro_endpoint)

# Get categories ----------------------------------------------------------

## some masterids repeated - take highest if same year, take most recent if different years

allLimsx <- allLims %>%
  group_by(masterid, Index, sampleyear) %>% ## group variables
  drop_na(IndexValue) %>%
  mutate(IndexValueMax = max(na.omit(IndexValue))) %>% ## take max
  dplyr::select(-IndexValue) %>% ## remove indexvalue
  ungroup(sampleyear) %>% ## remove year from group
  mutate(RecentYear = max(sampleyear)) %>% ## get most recernt year
  mutate(YearKeep = ifelse(sampleyear == RecentYear, "Yes", "No")) %>% ## add an indicator
  filter(YearKeep == "Yes") %>% ## remove less recent years
  dplyr::select(-c(RecentYear, YearKeep, sampleyear)) %>% ## remove extra columns
  distinct() %>% ## remove duplicates
  drop_na(Lower2) ## NAs places with no bio data

names(allLimsx)
  
## define if hydro is within mod limits & bio above threshold
imps <- allLimsx %>%
  group_by(metric, Index, Hydro_endpoint, Threshold, BioThresh, masterid, comid, Quantile) %>%
  mutate(WithinHydroLimits = ifelse(deltah_final <= Upper2 & deltah_final >= Lower2, "Within", "NotWithin")) %>%  ## within hydro limits
  mutate(WithinBioLimits = ifelse(IndexValueMax >= BioThresh, "Within", "NotWithin"))  %>% ## above bio thresholds
  ungroup() %>%
  dplyr::select( -c(metric, ChosenSmth, SmoothingFunc, CurveStatus))

unique(imps$Hydro_endpoint)

## check Nas

ind <- which(is.na(imps$WithinHydroLimits))
test <- imps[ind,] # 0

## save table in wide format with FFM ranges
write.csv(imps, "ignore/05_impact_ffm_bio_ffm_limits.csv")

imps <- read.csv("ignore/05_impact_ffm_bio_ffm_limits.csv")
imps

names(imps)

## get ranges and save

ranges <- imps %>%
  # rename(Q0.9 = "0.9", Q0.5 = "0.5") %>%
  select(Index, masterid, Threshold, BioThresh, Flow.Metric.Name, Quantile, Lower2, Upper2)

write.csv(ranges, "final_data/05_ffm_ranges.csv")


# categorize the sites ----------------------------------------------------

impsx <- imps %>%
  dplyr::select(-c(Lower2, Upper2,X.1, X, Lower, Upper,MetricCurve) ) %>%
  pivot_wider(names_from = Quantile, values_from = WithinHydroLimits) %>%
  rename(Q0.9 = "0.9", Q0.5 = "0.5") %>%
  mutate(Result = case_when((WithinBioLimits == "Within" & Q0.9 == "Within" & Q0.5 == "Within") ~ "HBUS", ## healthy and unlikely to be stressed
                            # ( WithinBioLimits == "NotWithin" & Q0.9 == "Within" & Q0.5 == "Within" & Q0.3 == "Within") ~ "UHVUF", ## Unhealthy and very unlikely to be stressed
                            ( WithinBioLimits == "NotWithin" & Q0.9 == "Within" & Q0.5 == "Within") ~ "UBUS", ## Unhealthy and unlikely to be stressed 
                            ( WithinBioLimits == "Within"& Q0.9 == "Within" & Q0.5 == "NotWithin") ~ "HBLS", ## healthy,  likely to be stressed
                            ( WithinBioLimits == "NotWithin" & Q0.9 == "Within" & Q0.5 == "NotWithin") ~ "UBLS", ## unhealthy, mod likely to be stressed
                            ( WithinBioLimits == "Within" & Q0.9 == "NotWithin" & Q0.5 == "NotWithin") ~ "HBVLS", ## healthy, very likely to be stressed
                            ( WithinBioLimits == "NotWithin" & Q0.9 == "NotWithin" & Q0.5 == "NotWithin") ~ "UBVLS")) %>% ## unhealthy, very likely to be stressed
  mutate(Result = factor(Result, levels = c( "HBUS", "UBUS", "HBLS", "UBLS", "HBVLS", "UBVLS"))) # %>%
  
unique(imps$Hydro_endpoint)
  ## check NAs
  # 
  # ind <- which(is.na(impsx$Result))
  # test <- impsx[ind,]
  # # # 
  # # # ## get masterids from na sites
  # # # 
  # checkids <- unique(test$masterid)
  # # # 
  # test2 <- allLimsx %>%
  #   filter(masterid %in% checkids, Hydro_endpoint %in% c("DS_Mag_50", "FA_Mag")) %>%
  #   select(Hydro_endpoint, Quantile,Lower2, Upper2, deltah_final)
  ## was a bug above, changed to 0 instead of min value

names(impsx)
write.csv(impsx, "ignore/05_impact_ffm_bio.csv")

unique(imps$Threshold)

imps <- read.csv("ignore/05_impact_ffm_bio.csv")

## count NAs in result
sum(is.na(imps$Result)) ## 0

## tally of impact per ffm

tallyImpact <- imps %>%
  group_by(Index, Hydro_endpoint, Flow.Metric.Name, Threshold, Result) %>%
  distinct() %>%
  tally() #%>%
  # drop_na(Result) %>%
  # mutate(PercChans = (n/sum(n)*100))

tallyImpact

write.csv(tallyImpact, "final_data/05_count_impact.csv")

# tallyImpact <- read.csv("output_data/05_count_impact.csv")
unique(tallyImpact$Hydro_endpoint)

# Number of strikes -------------------------------------------------------

head(impsx)

unique(impsx$masterid)

## count number of FFM that are in UHLF per site

strikes <- imps %>%
  select(masterid, Threshold, Index, Hydro_endpoint, Result) %>%
  group_by(masterid, Threshold, Index, Result) %>%
  distinct() %>%
  tally()

strikes

write.csv(strikes, "final_data/05_Number_ffm_per_result.csv")

unique(impsx$Hydro_endpoint)

