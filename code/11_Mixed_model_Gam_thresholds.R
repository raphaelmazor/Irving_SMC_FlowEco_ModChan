## get GAM thresholds

library(gam)
library(tidyverse)
library(sf)
library(sf)
library(mgcv)
library(itsadug)
library(tidylog)

## function to find value in curve
# load("code/functions/root_interpolation_function.Rdata")

## function to get delta H at csci score
load(file = "code/functions/10_find_deltaH_exact_values.Rdata")
## load DF 

load(file = "ignore/09_mixed_effects_model_predictions.RData")
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


## smoothing functions per FFM
sm6csci <- c("DS_Mag_50", "FA_Mag", "Peak_10", "Peak_2", "Peak_5", "Wet_BFL_Mag_50", "Wet_BFL_Mag_10")

sm3csci <- c( "SP_Mag",  "Q99")

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
  mutate(comb_code_type = paste0(Metric, "_", hydro.endpoints, "_", Smooths)) %>% ## make metric
  rename(HydroValue = deltah_final, BioValue = SlopeAndIntercept) %>% ## change names
  select(-Intercept, -Slope) %>% ## remove other models
  mutate(ChosenSmth = case_when((Metric == "csci" & hydro.endpoints %in% sm6csci) ~ 6,
                                (Metric == "csci" & hydro.endpoints %in% sm3csci) ~ 3)) %>% ## define correct smooths
  mutate(ModelstoUse = ifelse(Smooths == ChosenSmth, "Yes", "No")) %>%
  filter(ModelstoUse == "Yes")  ## take only smooths needed
  

head(dfx)

## join with labels
all_data <- left_join(dfx, labels, by ="hydro.endpoints")
head(all_data)

###Gams

load(file = "final_data/09_mixed_models.RData")
gam_lme

## look up table 
bio_h_summary <- read.csv("final_data/09_coefs_mixed_effects_model.csv")
bio_h_summary


## add row number as a column for an index and change FFM names

bio_h_summary <- bio_h_summary %>%
  mutate(Index = row.names(bio_h_summary)) %>%
  rename(flow.endpoints = Variable) %>%
  mutate(Hydro_endpoint = case_when(flow.endpoints == "d_ds_mag_50" ~ "DS_Mag_50",            
                                    flow.endpoints == "d_fa_mag" ~ "FA_Mag",
                                    flow.endpoints == "d_peak_10" ~ "Peak_10",
                                    flow.endpoints == "d_peak_2" ~ "Peak_2",
                                    flow.endpoints == "d_peak_5" ~ "Peak_5",
                                    flow.endpoints == "d_sp_mag" ~ "SP_Mag",
                                    flow.endpoints == "d_wet_bfl_mag_10" ~ "Wet_BFL_Mag_10",
                                    flow.endpoints == "d_wet_bfl_mag_50" ~ "Wet_BFL_Mag_50", 
                                    flow.endpoints == "delta_q99" ~ "Q99"))
bio_h_summary
## filter to slope/intercept and smoothing functions per ffm
Med_GAM_models <- bio_h_summary %>%
  mutate(ChosenSmth = case_when((Metric == "csci" & Hydro_endpoint %in% sm6csci) ~ 6,
                                (Metric == "csci" & Hydro_endpoint %in% sm3csci) ~ 3)) %>%
  mutate(ModelstoUse = ifelse(Smooths == ChosenSmth, "Yes", "No")) %>%
  filter(ModelstoUse == "Yes", Model == "gamm_int_slope")


# calculate thresholds ----------------------------------------------------

## create df
df <- as.data.frame(matrix(ncol=21))
colnames(df) <- c("metric", "ThresholdNAT_Lower", "ThresholdNAT_Upper", "ThresholdNATMed_Lower", "ThresholdNATMed_Upper",
                  "ThresholdNATLow_Lower", "ThresholdNATLow_Upper",
                  "ThresholdNATHigh_Lower", "ThresholdNATHigh_Upper", "ThresholdHB_Lower", "ThresholdHB_Upper",
                  "ThresholdSB0_Lower", "ThresholdSB0_Upper","ThresholdSB1_Lower", "ThresholdSB1_Upper",
                  "ThresholdSB2_Lower", "ThresholdSB2_Upper",
                  "n", "Index", "Hydro_endpoint","SmoothingFunc")
df

head(all_data)
## define FFMs

ffms <- unique(AllData$flow_metric)
ffms

i=1
i
## loop through metrics
for(i in 1:length(ffms)) {
  
  ## define metric
  met <- ffms[i]
  
  ## filter by metric
  hydroxx <- all_data %>%
    filter(Variable == met)
  
  ## get ffm from bio summary
  
  ffmIndex <- Med_GAM_models %>%
    filter(flow.endpoints == ffms[i])
  
  ## get index for model
  Ind <- as.numeric(ffmIndex$Index)
  ## extract model for ffm
  mod <- gam_lme[[Ind]]
  
  ## find ranges of deltah from model

  ## standard
    threshNATMed <- find_delta_h_exact(hydroxx$HydroValue, mod,0.79, "NAT") ## sep ref threshold for overall tally
    
    ## low
    threshNATLow <- find_delta_h_exact(hydroxx$HydroValue, mod, 0.63,  "NAT")
   
    # high
    threshNATHigh <- find_delta_h_exact(hydroxx$HydroValue, mod, 0.92,  "NAT")
    
    ## get curves values at different scores get closest value to zero for both -ve and +ve
    threshNAT <- find_delta_h_exact(hydroxx$HydroValue, mod, 0.79,   "NAT")
    
    ## hard bottom
    threshHB <- find_delta_h_exact(hydroxx$HydroValue, mod,  0.67,"NAT")
    
    ## soft bottom 0
    threshSB0 <- find_delta_h_exact(hydroxx$HydroValue, mod, 0.78, "NAT")
    threshSB0
    ## soft bottom 2
    threshSB2 <- find_delta_h_exact(hydroxx$HydroValue, mod, 0.75, "NAT")
    # threshSB2 <- ifelse(length(threshSB2) == 0, NA, threshSB2)
    
    threshSB1 <- find_delta_h_exact(hydroxx$HydroValue, mod, 0.79, "NAT")

  ## add info to df
  df[i, 1] <- met
  df[i, 2] <- threshNAT$delta_h_negative
  df[i, 3] <- threshNAT$delta_h_positive
  df[i, 4] <- threshNATMed$delta_h_negative
  df[i, 5] <- threshNATMed$delta_h_positive
  df[i, 6] <- threshNATLow$delta_h_negative
  df[i, 7] <- threshNATLow$delta_h_positive
  df[i, 8] <- threshNATHigh$delta_h_negative
  df[i, 9] <- threshNATHigh$delta_h_positive
  df[i, 10] <- threshHB$delta_h_negative
  df[i, 11] <- threshHB$delta_h_positive
  df[i, 12] <- threshSB0$delta_h_negative
  df[i, 13] <- threshSB0$delta_h_positive
  df[i, 14] <- threshSB1$delta_h_negative
  df[i, 15] <- threshSB1$delta_h_positive
  df[i, 16] <- threshSB2$delta_h_negative
  df[i, 17] <- threshSB2$delta_h_positive
  df[i, 18] <- length(hydroxx$BioValue) ## not needed here but keep in for coding 
  df[i ,19] <- hydroxx$Metric[1]
  df[i, 20] <- hydroxx$hydro.endpoints[1]
  df[i, 21] <- hydroxx$Smooths[1]
  # df[i, 20] <- paste(hydroxx$Quant[1])

  
}

df
df ## NAs - where curve doesn't reach threshold
write.csv(df, "final_data/11_delta_thresholds_Mixed_Effects_GAMs.csv")

# df <- read.csv("final_data/11_delta_thresholds_GAMs.csv")

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
  filter(ModelstoUse == "Yes") 

unique(dfsub$Hydro_endpoint)
  
write.csv(dfsub, "final_data/11_chosen_models_thresholds_Mixed_Effects.csv")
dfsub


# format to apply data to categories----------------------------------------------------

df <- read.csv("final_data/11_chosen_models_thresholds_Mixed_Effects.csv")
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
                               (Index == "csci" & Threshold == "SB2") ~ 0.75,
                               (Index == "csci" & Threshold == "SB1") ~ 0.79)) 

limits

unique(limits$Threshold)
## make wider with type - lower/upper
limits <- limits %>%
  pivot_wider(names_from = Type, values_from = DeltaH)  %>%
  mutate(MetricCurve = paste0(metric,"_", SmoothingFunc)) %>%
  filter(!Threshold %in% c("NATMed",  "NATLow" , "NATHigh" )) ## remove standard thresholds as don't need them

## organise limits by type of curve: AboveCurve, BelowCurve
# ## 0.5s should be below, and 0.9s above - there are some exceptions in 0.5 - we'll adjust manually below
# AboveCurve <- c("NAT", "SB1")
# BelowCurve <- c("HB", "SB0", "SB2", "NATLow")
# 
# ## assign labels to know which curves are above or below thresholds
# ## some thresholds hit curve, so make sure it's NA too
# limits <- limits %>%
#   mutate(CurveStatus = ifelse(is.na(Upper) & Threshold %in% AboveCurve, "Above", "Passes")) %>%
#   mutate(CurveStatus = ifelse(is.na(Lower) & Threshold  %in% BelowCurve, "Below", CurveStatus)) 

## Q99 doesn't reach beyond zero, add zero as upper limit
# limits <- limits %>%
#   mutate(Upper = ifelse(CurveStatus == "Above" & Hydro_endpoint == "Q99", 0, Upper))


## an NA means the NAT curve is above the threshold for HB or SBO
## add in min and max in these instances

MinDel <- min(na.omit(AllDataLongx$deltah_final)) ## -51409.35
MaxDel <- max(na.omit(AllDataLongx$deltah_final)) ## 8059.891

## if a curve is below the threshold, value should be zero - should result in "outside limits"
## Q99 and peak 5 have funky 0.3 curves, should have lower and upper limits of 0 too
limits2 <- limits %>%
  mutate(Lower = ifelse(is.na(Lower), MinDel, Lower),
         Upper = ifelse(is.na(Upper), MaxDel, Upper)) #%>%

limits2

limits4 <- limits2
## save out
write.csv(limits4, "final_data/11__Mixed_Effects_GAMs_FFM_Ranges_With_Adjustments.csv")

limits4 <- read.csv("final_data/11__Mixed_Effects_GAMs_FFM_Ranges_With_Adjustments.csv")

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
  distinct() #%>% ## remove duplicates
  # drop_na(Lower2) ## NAs places with no bio data

names(allLimsx)
  
## define if hydro is within mod limits & bio above threshold
imps <- allLimsx %>%
  group_by(metric, Index, Hydro_endpoint, Threshold, BioThresh, masterid, comid) %>%
  mutate(WithinHydroLimits = ifelse(deltah_final <= Upper & deltah_final >= Lower, "Within", "NotWithin")) %>%  ## within hydro limits
  mutate(WithinBioLimits = ifelse(IndexValueMax >= BioThresh, "Within", "NotWithin"))  %>% ## above bio thresholds
  ungroup() %>%
  dplyr::select( -c(metric, ChosenSmth, SmoothingFunc)) %>%
  filter(!Threshold %in% c("NATHigh", "SB1") ) %>% ## remove high threshold
  drop_na()

unique(imps$Hydro_endpoint)

## check Nas

# ind <- which(is.na(imps$WithinHydroLimits))
# test <- imps[ind,] # 0
# test ## nathighs and SB1s can be removed

## save table in wide format with FFM ranges
write.csv(imps, "ignore/11_impact_ffm_bio_ffm_limits_Mixed_Effects.csv")

imps <- read.csv("ignore/11_impact_ffm_bio_ffm_limits_Mixed_Effects.csv")
imps

names(imps)

## get ranges and save
ranges
ranges <- imps %>%
  # rename(Q0.9 = "0.9", Q0.5 = "0.5") %>%
  select(Index, masterid, Threshold, BioThresh, Flow.Metric.Name, Lower, Upper)

write.csv(ranges, "final_data/11_ffm_ranges_Mixed_Effects.csv")

