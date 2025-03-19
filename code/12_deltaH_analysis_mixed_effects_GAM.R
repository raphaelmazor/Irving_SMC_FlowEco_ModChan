## using different curves: 50 GAM, Quantile GAM, Mixed Model GAM
## find delta H to get to: 0.79, 0.67, incremental improvement from obs CSCI

## packages
library(qgam); library(MASS)
library(tidyverse)
library(sf)
library(mgcv)
library(itsadug)
# library(mapview)
library(tidylog)
library(RColorBrewer)

## labels for FFMs
labels <- read.csv("input_data/ffm_names.csv")
labels <- labels[1:24, ]
labels <- labels %>% rename(hydro.endpoints = Flow.Metric.Code)
labels[25, 1] <- "Peak Flow Magnitude (Q99, cfs)"
labels[25, 2] <- "Q99"
labels[25, 3] <- "Peak Flow Magnitude"
labels

## load in root linear interpolation function
load("code/functions/root_interpolation_function.Rdata")

## function to find closest values to zero
load("code/functions/find_closest_values.RData")

## function to get delta H at csci score
load(file = "code/functions/10_find_deltaH_mixed_model.Rdata")

## directory for figures
out.dir <- "final_figures/"

# Upload GAM Model -----------------------------------------------------------

### Median Gam

load(file = "ignore/09_mixed_models.RData")
gam_lme

## look up table 
bio_h_summary <- read.csv("final_data/09_coefs_mixed_effects_model.csv")
bio_h_summary
## smoothing functions per FFM
sm6csci <- c("DS_Mag_50", "FA_Mag", "Peak_10", "Peak_2", "Peak_5", "Wet_BFL_Mag_50", "Wet_BFL_Mag_10")

sm3csci <- c( "SP_Mag",  "Q99")


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
head(bio_h_summary)
## filter to slope/intercept and smoothing functions per ffm
Med_GAM_models <- bio_h_summary %>%
  mutate(ChosenSmth = case_when((Metric == "csci" & Hydro_endpoint %in% sm6csci) ~ 6,
                                (Metric == "csci" & Hydro_endpoint %in% sm3csci) ~ 3)) %>%
  mutate(ModelstoUse = ifelse(Smooths == ChosenSmth, "Yes", "No")) %>%
  filter(ModelstoUse == "Yes", Model == "gamm_int_slope")
  

Med_GAM_models

## upload thresholds

df <- read.csv("final_data/11__Mixed_Effects_GAMs_FFM_Ranges_With_Adjustments.csv") %>%
  # filter( Threshold %in%c("NAT", "SB2","SB0", "HB")) %>%
  select(Hydro_endpoint, Threshold, BioThresh, Lower, Upper) 
head(df)


# Upload Site data ---------------------------------------------------

## data
load(file = "final_data/01_bugs_algae_flow_joined_by_masterid.RData")  

AllData <- AllData %>%
  drop_na(channel_engineering_class)

head(AllData)

## take the ffm, get a site, see where it hits the curve in delta H and score
##
## how far to get it to the delta H limits?
## how far for an incremental increase?

### loop around FFMs
## empty df
ChangeDx <- NULL
f=1
## define FFMs

ffms <- unique(AllData$flow_metric)
ffms
for(f in 1:length(ffms)) {
  
  ## filter to FFM
  FFMDf <- AllData %>%
    filter(flow_metric == ffms[f]) %>%
    distinct() %>%
    select(masterid, sampleyear, MetricValue, deltah_final, hydro.endpoints, channel_engineering_class) %>%
    drop_na(MetricValue)
  
  ## get ffm from bio summary
  
  ffmIndex <- Med_GAM_models %>%
    filter(flow.endpoints == ffms[f])
  
  ## get index for model
  Ind <- as.numeric(ffmIndex$Index)
  
  ## extract model for ffm
  mod <- gam_lme[[Ind]]
  
  ## findgam_lme## find DH for incremental value - 0.1 increase of CSCI
  
  ## define increment values
  FFMDf <- FFMDf %>%
    mutate(PredCSCI = predict(mod, FFMDf)) %>% ## predict csci 
    mutate(csci_Incr = PredCSCI+0.1) ## increment

  ## define sites/rownames - adjust here to remove duplicate sites
  sites <- rownames(FFMDf)
  head(FFMDf)
  # loop over sites

  for(s in 1:length(sites)) {

    siteDelta <- FFMDf[s,"deltah_final"] ## get delta
    siteIncr <- FFMDf[s,"csci_Incr"] ## get adjusted csci
    ChanEng <- FFMDf[s,"channel_engineering_class"] ## get adjusted csci
    target_csci_score_INC <-  siteIncr  #  target csci_score - predicted csci increased by 0.1
    
    ## use function to predict delta H for both neg and pos values
    predicted_delta_h <- find_delta_h_optimize_mm(target_csci_score, mod, min(FFMDf$deltah_final), max(FFMDf$deltah_final), siteDelta, ChanEng)
   
    ## add to main df## aFFMDfdd to main df
    FFMDf$newDeltaINC[s] <- ifelse(siteDelta > 0, predicted_delta_h$delta_h_positive, predicted_delta_h$delta_h_negative)
    
  }
  

  FFMDf <- FFMDf %>%
    mutate(newDeltaINC = round(newDeltaINC, digits = 2)) %>%
    mutate(ChangeInDeltaINC = deltah_final - newDeltaINC) %>% ## absolute change in delta for inc 
    mutate(ChangeInDeltaINC = round(ChangeInDeltaINC, digits = 2)) %>%
    mutate(RelChangeInDeltaINC = (abs(ChangeInDeltaINC)/abs(deltah_final))*100) %>% ## relative change
    mutate(RelChangeInDeltaINC = round(RelChangeInDeltaINC, digits = 2))

  ## delta H will either be within the limits or out.
  ## if it is within, then issue is likely not flow
  ## if it is outside then can calculate change in flow needed for an improved score
  
  ## join site to ffm limits data by ffm and threshold
  
  FFMDFLims <- right_join(FFMDf, df, by =c("hydro.endpoints" = "Hydro_endpoint")) %>%
    rename(ChannelType = channel_engineering_class) %>%
    drop_na(deltah_final)
  
  head(FFMDFLims) ## missing sites are SB0, these can be changed to SB2 for now. 

  
  ## add column of whether within limits or not
  FFMDFLims <- FFMDFLims %>%
    # mutate(WithinLims = ifelse(Threshold == "NAT" & deltah_final >= Lower & deltah_final <= Upper, "Yes", "No")) %>% ## no flow problem
    mutate(HealthyCSCI = ifelse(MetricValue >= 0.79, "Yes", "No")) %>% ## add column for healthy csci (remove in figure later)
    mutate(HealthyCSCIBO = ifelse(MetricValue < 0.79 & MetricValue >= BioThresh, "Yes", "No")) #%>% ##  flow within BO range
    # mutate(WithinLimsBO = ifelse(deltah_final >= Lower & deltah_final <= Upper, "Yes", "No"))
  
  # Delta H to thresholds ---------------------------------------------------
  
  
  ## if delta H is positive it must be above the upper limit
  ## calculate difference of delta at site to upper limit
  ## if negative then calculate difference to lower limit
  
  ChangeD <- FFMDFLims %>%
    mutate(Keep = ifelse(ChannelType == Threshold | Threshold == "NAT", "Yes", "No")) %>% ## define channel that match with thresholds
    filter(Keep == "Yes") %>% ## remove those that don't
    mutate(ChangeInDelta = ifelse(deltah_final >= 0, deltah_final- Upper, Lower - deltah_final)) %>% ## absolute change
    # mutate(ChangeInDelta = ifelse(WithinLims == "Yes" | WithinLimsBO == "Yes", NA, ChangeInDelta)) %>% ## keep as NA if DH within limits
    mutate(newDeltaBO = ifelse(deltah_final < 0, deltah_final + ChangeInDelta, deltah_final - ChangeInDelta)) %>% ## DH target
    mutate(RelChangeInDelta = (abs(ChangeInDelta)/abs(deltah_final))*100) %>% ## relative change
    mutate(DeltaSD = sd(deltah_final)) ## get sd of delta 

 # head(ChangeD)
  
  
  ChangeDx <- rbind(ChangeDx,ChangeD)
  
}

# 412M08642
# 408M03005
head(ChangeDx)

## save
write.csv(ChangeDx, "ignore/12_Change_in_delta_to_hit_thresholds_Mixed_Model.csv")

## plot 

# p1 <- ggplot(data = ChangeDx, aes(x=RelChangeInDelta, y= MetricValue, group = ChannelType, col = ChannelType)) +
#   geom_smooth()+
#   facet_wrap(~hydro.endpoints, scales = "free")
# 
# p1 ## plot of rel change vs csci score by channel type and ffm facet
# 
# p2 <- ggplot(data = ChangeDx, aes(x=RelChangeInDelta, y= MetricValue, group = hydro.endpoints, col = hydro.endpoints)) +
#   geom_smooth()+
#   facet_wrap(~ChannelType, scales = "free")
# 
# p2 ## plot of rel change vs csci score per ffm and channel type facet


# Calculate categories ----------------------------------------------------

ChangeDx <- read.csv("ignore/12_Change_in_delta_to_hit_thresholds_Mixed_Model.csv")
head(ChangeDx)

## also filter dataset, only need the rows where channel type matches threshold
## relative change
ChangeDxWide <- ChangeDx %>%
  select(-c( X)) %>% ##
  mutate(Match = ifelse(Threshold %in% c("NAT") , "Yes", "No")) %>% ## keep all ref threshold for Type B 
  mutate(Match = ifelse(Threshold %in% c("NAT") & ChannelType %in% c("SB1"), "Yes", Match)) %>% ## for SB1s Type B
  mutate(Match = ifelse(Threshold %in% c("SB0", "SB2") & ChannelType %in% c("SB0", "SB2"), "Yes", Match)) %>% ## for SBs for type c
  mutate(Match = ifelse(Threshold == "HB" & ChannelType == "HB", "Yes", Match)) %>% ## all other thresholds that match channel types for type C
  filter(Match == "Yes") %>%
  distinct() %>%
  select(-c( Keep, Match))


unique(ChangeDxWide$Threshold)

## save 
write.csv(ChangeDxWide, "ignore/12_matching_channel_typw_threshold_Mixed_Model.csv")
## remember to remove duplicates/take max or mean of same year samples!!!!!!!

## categories

# Biology	  Hydrology	  Action
# Above ref	Within ref	No action/protect
#           Outside ref	No action/monitor
# 
# 
# Above BOs	Within ref	Investigate nonflow
#           Outside ref	Identify improvement needed to get ref (<10, 10 to 50, >50)
#                       Identify improvement needed to get +0.1 CSCI (<10, 10 to 50, >50)
# 
# Below BO	Within ref	Investigate nonflow
#           Outside ref	Identify improvement needed to get ref (<10, 10 to 50, >50)
#                       Identify improvement needed to get bo (<10, 10 to 50, >50)
#                       Identify improvement needed to get +0.1 CSCI (<10, 10 to 50, >50)


## function to get the categories
head(ChangeDxWide)

## add column for the bio criteria
BioChange <- ChangeDxWide %>%
  group_by(masterid, hydro.endpoints, sampleyear) %>% ## group data
  mutate(HealthyCSCIBO2 = last(HealthyCSCIBO)) %>% ## take the last value of group and fill = "Yes"
  mutate(BioStatus = case_when(
    HealthyCSCI == "Yes" ~ "Above Reference",
    HealthyCSCIBO2 == "Yes" ~ "Above Best Observed",
    HealthyCSCI == "No" & HealthyCSCIBO == "No" ~ "Below Best Observed"
  ))

## add column for flow criteria
FlowChange <- BioChange %>%
  mutate(FlowStatus = case_when( ## define within/outside for NAT threshold
    Threshold == "NAT" & deltah_final >= Lower & deltah_final <= Upper ~ "Within Reference",
    Threshold == "NAT" & (deltah_final < Lower | deltah_final > Upper) ~ "Outside Reference"
  )) %>% ## define within/outside for NAT threshold but mod channel types
  group_by(masterid, hydro.endpoints, sampleyear) %>%
  mutate(FlowStatus2 = first(FlowStatus)) %>% ## take the status for the masterid and ffm
  select(-FlowStatus) %>%
  mutate(BOFlows = case_when( ## define within/outside for BO threshold - to add sites within BO as NAs later
    Threshold != "NAT" & deltah_final >= Lower & deltah_final <= Upper ~ "Within BO",
    Threshold != "NAT" & (deltah_final < Lower | deltah_final > Upper) ~ "Outside BO"
  )) %>%
  mutate(BOFlows = last(BOFlows)) #%>% ## take the status for the masterid and ffm
  

class(FlowChange)

## Add action column
AllChange <- FlowChange %>%
  mutate(Action = case_when(
    BioStatus == "Above Reference" & FlowStatus2 == "Within Reference" ~ "No Action/Protect",
    BioStatus == "Above Reference" & FlowStatus2 == "Outside Reference" ~ "No Action/Monitor",
    BioStatus == "Above Best Observed" & FlowStatus2 == "Within Reference" ~ "Non-flow related",
    BioStatus == "Above Best Observed" & FlowStatus2 == "Outside Reference" ~ "Identify improvement needed",
    BioStatus == "Below Best Observed" & FlowStatus2 == "Within Reference" ~ "Non-flow related",
    BioStatus == "Below Best Observed" & FlowStatus2 == "Outside Reference" ~ "Identify improvement needed",
  ))
names(AllChange)

## Identify improvement needed for each scenario
ImpNeeded <- AllChange %>%
  select(masterid, sampleyear,ChannelType, PredCSCI, csci_Incr, MetricValue, RelChangeInDeltaINC, Threshold, RelChangeInDelta, BioStatus:Action) %>%
  mutate(Threshold2 = ifelse(Threshold == "NAT", "Standard", "Modified")) %>% ## put all mod thresholds together as they don't overlap 
  select(-Threshold) %>%
  pivot_wider(names_from = "Threshold2", values_from = RelChangeInDelta) %>% ## pivot longer so can take from column
  mutate(RefImprovement = case_when( ## improvement needed to get ref flows
    BioStatus == "Above Reference" & Action == "No Action/Protect" ~ "No Action", ## no action
    BioStatus == "Above Reference" & Action == "No Action/Monitor"~ "No Action",  ## no action
    ## above BO
    BioStatus == "Above Best Observed" & Action == "Non-flow related" ~ "Non-flow related", ## other stressors
    BioStatus == "Above Best Observed" & Action == "Identify improvement needed" & Standard <= 10 ~ "Less than 10%",
    BioStatus == "Above Best Observed" & Action == "Identify improvement needed" & Standard > 10 & Standard <= 50 ~ "10 - 50%",
    BioStatus == "Above Best Observed" & Action == "Identify improvement needed" & Standard > 50 ~ "Over 50%",
    ## Below BO
    BioStatus == "Below Best Observed" & Action == "Non-flow related" ~ "Non-flow related", ## other stressors
    BioStatus == "Below Best Observed" & Action == "Identify improvement needed" & Standard <= 10 ~ "Less than 10%",
    BioStatus == "Below Best Observed" & Action == "Identify improvement needed" & Standard > 10 & Standard <= 50 ~ "10 - 50%",
    BioStatus == "Below Best Observed" & Action == "Identify improvement needed" & Standard > 50 ~ "Over 50%"
  )) %>%  ## improvement needed to get INC flows
  mutate(INCImprovement = case_when( ## improvement needed 
    BioStatus == "Above Reference" & Action == "No Action/Protect" ~ "No Action", ## no action
    BioStatus == "Above Reference" & Action == "No Action/Monitor"~ "No Action",  ## no action
    ## above BO
    BioStatus == "Above Best Observed" & Action == "Non-flow related" ~ "Non-flow related", ## other stressors
    BioStatus == "Above Best Observed" & Action == "Identify improvement needed" & RelChangeInDeltaINC <= 10 ~ "Less than 10%",
    BioStatus == "Above Best Observed" & Action == "Identify improvement needed" & RelChangeInDeltaINC > 10 & RelChangeInDeltaINC <= 50 ~ "10 - 50%",
    BioStatus == "Above Best Observed" & Action == "Identify improvement needed" & RelChangeInDeltaINC > 50 ~ "Over 50%",
    ## Below BO
    BioStatus == "Below Best Observed" & Action == "Non-flow related" ~ "Non-flow related", ## other stressors
    BioStatus == "Below Best Observed" & Action == "Identify improvement needed" & RelChangeInDeltaINC <= 10 ~ "Less than 10%",
    BioStatus == "Below Best Observed" & Action == "Identify improvement needed" & RelChangeInDeltaINC > 10 & RelChangeInDeltaINC <= 50 ~ "10 - 50%",
    BioStatus == "Below Best Observed" & Action == "Identify improvement needed" & RelChangeInDeltaINC > 50 ~ "Over 50%"
  )) %>%  ## improvement needed to get BO flows
  mutate(BOImprovement = case_when( ## improvement needed to get
    BioStatus == "Above Reference" & Action == "No Action/Protect" ~ "No Action", ## no action
    BioStatus == "Above Reference" & Action == "No Action/Monitor"~ "No Action",  ## no action
    ## above BO
    BioStatus == "Above Best Observed" & Action == "Non-flow related" ~ "Non-flow related", ## other stressors
    BioStatus == "Above Best Observed" & Action == "Identify improvement needed"  ~ NA,
    BioStatus == "Above Best Observed" & Action == "Identify improvement needed" ~ NA,
    BioStatus == "Above Best Observed" & Action == "Identify improvement needed" ~ NA,
    ## Below BO
    ## if sites are within BO flows then NA (as only ref and INC practical)
    BioStatus == "Below Best Observed" & Action == "Identify improvement needed" & BOFlows == "Within BO" ~ NA,
    BioStatus == "Below Best Observed" & Action == "Non-flow related" ~ "Non-flow related", ## other stressors
    BioStatus == "Below Best Observed" & Action == "Identify improvement needed" & Modified <= 10 ~ "Less than 10%",
    BioStatus == "Below Best Observed" & Action == "Identify improvement needed" & Modified > 10 & Modified <= 50 ~ "10 - 50%",
    BioStatus == "Below Best Observed" & Action == "Identify improvement needed" & Modified > 50 ~ "Over 50%"
  ))

## note that there just aren't many sites above BO and below NAT
names(ImpNeeded)
## can add col here for "best" for each site???
## save
write.csv(ImpNeeded, "ignore/12_improvements_in_flow.csv")
## tidy up DF with only columns needed

categories <- ImpNeeded 

# view(categories) 

## save 
write.csv(categories, "final_data/12_categories_Mixed_Model.csv")
head(categories)
str(ChangeDxWide)

## remove values from same site - if from same year max(), if from different years then most recent - can change
categoriesx <- categories %>%
  group_by(masterid, sampleyear) %>% ## group by site and year and ffm
  mutate(MaxScore = max(MetricValue)) %>%
  mutate(Match = ifelse(MaxScore == MetricValue, "Yes", "No")) %>%
  filter(Match == "Yes" ) %>%
  select(-Match) %>%
  group_by(masterid, hydro.endpoints, ChannelType) %>% ## group by site and year and ffm
  mutate(sampleyear2 = max(sampleyear)) %>% ## get max year i.e., most recent - change here for lowest csci etc
  mutate(Match = ifelse(sampleyear2 == sampleyear, "Yes", "No")) %>%
  filter(Match == "Yes") %>%
  dplyr::select(-sampleyear2, -Match) %>% ## remove original score
  distinct()

write.csv(categoriesx, "final_data/12_categories_Mixed_Model.csv")

## count catergories per channel type
names(categoriesx)
length(unique(categoriesx$masterid)) ## 396

## tally per channel type and target 
tallyCats1 <- categoriesx %>%
  select(-c(RelChangeInDeltaINC:Modified)) %>%
  pivot_longer(RefImprovement:BOImprovement, names_to = "Target", values_to = "Categories") %>% ## make long
  drop_na(Categories) %>% ## remove NAs - sites that are below BO csci but within BO flows - focus on INC and ref improvement
  mutate(ChannelType = ifelse(ChannelType %in% c("SB0", "SB2"), "SB", ChannelType)) %>% ## add SB2 & 0 together
  mutate(ChannelType = ifelse(ChannelType %in% c("NAT", "SB1"), "NAT", ChannelType)) %>% ## add SB1 &NAT together
  group_by(ChannelType, hydro.endpoints, Target) %>% ## group
  mutate(TotalSites = length(unique(masterid))) %>% ## get total n sites
  group_by(ChannelType, hydro.endpoints,Target, Categories) %>% # add categories gto group
  mutate(Tally = length(unique(masterid))) %>% # count sites per category
  mutate(PercentSites = Tally/TotalSites*100) %>% ## get percentages
  select(Tally, TotalSites, PercentSites) %>% ## remove cols
  distinct() %>%
  drop_na(ChannelType)

## get percentages for overall column
tallyCatsOV <- categoriesx %>%
  select(-c(RelChangeInDeltaINC:Modified)) %>%
  pivot_longer(RefImprovement:BOImprovement, names_to = "Target", values_to = "Categories") %>%
  drop_na(Categories) %>% ## remove NAs - sites that are below BO csci but within BO flows - focus on INC and ref improvement
  group_by(hydro.endpoints, Target) %>%
  mutate(TotalSites = length(unique(masterid))) %>% ## get total sites for overall column
  group_by(hydro.endpoints, Target, Categories) %>% ## group by ffm and cats
  mutate(Tally = length(unique(masterid))) %>% ## count sites for each ffm and cat
  mutate(PercentSites = Tally/TotalSites*100) %>% ## percentage per ffm and cat
  select(Tally, TotalSites, PercentSites) %>%
  mutate(ChannelType = "All") %>%
  distinct() %>%
  drop_na(ChannelType)

## add channel tally and overall tally together

tallyCats2 <- bind_rows(tallyCats1, tallyCatsOV)
unique(tallyCats2$Categories)  
## 801M16916
## SMC02270

## make factor for figure 
tallyCats <- tallyCats2 %>%
  mutate(ChannelType = factor(ChannelType, levels = c("NAT", "SB", "HB", "All"),
                              labels = c("NAT & SB1", "SB0 & SB2", "HB", "All"))) %>%
  inner_join(labels, by = "hydro.endpoints") %>%
  mutate(Categories = factor(Categories, levels = c("No Action","Non-flow related","Less than 10%","10 - 50%","Over 50%")))

## 
unique(tallyCats$ChannelType)
unique(tallyCats$Categories)
range(tallyCats$PercentSites)

## save
write.csv(tallyCats, "final_data/12_tally_categories_per_target_for_figures_V2.csv")
unique(categoriesx1$ChannelType)
unique(categoriesx1$Flow.Metric.Name)

## tally per channel type and target - sites can be in more than category so use length(masterid) not unique
## OR take best target per category
head(categoriesxRed)
unique(categoriesxRed$ChannelType)

# categoriesxRed <- categoriesx %>%
#   select(-c(RelChangeInDeltaINC, BioStatus:Modified, MaxScore)) %>%
#   pivot_longer(RefImprovement:BOImprovement, names_to = "Target", values_to = "Categories") %>% ## make long
#   drop_na(Categories) %>% ## remove NAs - sites that are below or above BO csci but within BO flows - focus on INC and ref improvement
#   mutate(Target = factor(Target, levels = c("RefImprovement", "BOImprovement", "INCImprovement"),
#                        labels = c("Reference Condition", "Best Observed Condition", "Incremental Improvement"))) %>%
#   filter(!Categories %in% c("No Action", "Non-flow related")) %>% ## keep only categoried sites for now
#   group_by(masterid, hydro.endpoints, Categories) %>% # group
#   mutate(Target2 = first(Target))# %>% ## some INc taken instead of BO cond
#   mutate(BioThresh = case_when( ## add threshold value
#     ChannelType == "HB" ~ 0.67,
#     ChannelType == "NAT" ~ 0.79,
#     ChannelType == "SB0" ~ 0.78,
#     ChannelType == "SB1" ~ 0.79,
#     ChannelType == "SB2" ~ 0.75
#   )) %>% ## take last when incremental increase is less than BO
#   mutate(Target3 = ifelse(Target2 != "Reference Condition" & csci_Incr < BioThresh, last(Target), Target2)) %>%
#   mutate(Target2 = case_when(
#     Target3 == 1 ~ "Reference Condition",
#     Target3 == 2 ~ "Best Observed Condition",
#     Target3 == 3 ~ "Incremental Improvement",
#   )) %>%
#   select(-Target, -Target3)
 
## count sites per category
tallyCats1c <- categoriesx %>%
  pivot_longer(RefImprovement:BOImprovement, names_to = "Target", values_to = "Categories") %>% ## make long
  drop_na(Categories) %>% ## remove NAs - sites that are below BO csci but within BO flows - focus on INC and ref improvement
  mutate(ChannelType = ifelse(ChannelType %in% c("SB0", "SB2"), "SB", ChannelType)) %>% ## add SB2 & 0 together
  mutate(ChannelType = ifelse(ChannelType %in% c("NAT", "SB1"), "NAT", ChannelType)) %>% ## add SB1 &NAT together
  group_by(ChannelType, hydro.endpoints, Categories) %>% ## group
  mutate(TotalSites = length(masterid)) %>% ## get total n sites
  group_by(ChannelType, hydro.endpoints,Target, Categories) %>% # add categories gto group
  mutate(Tally = length(masterid)) %>% # count sites per category
  mutate(PercentSites = Tally/TotalSites*100) %>% ## get percentages
  select(Tally, TotalSites, PercentSites) %>% ## remove cols
  distinct() %>%
  drop_na(ChannelType)

tallyCats1c
## get percentages for overall column
tallyCatsOVc <- categoriesx %>%
  pivot_longer(RefImprovement:BOImprovement, names_to = "Target", values_to = "Categories") %>%
  drop_na(Categories) %>% ## remove NAs - sites that are below BO csci but within BO flows - focus on INC and ref improvement
  group_by(hydro.endpoints, Categories) %>%
  mutate(TotalSites = length(masterid)) %>% ## get total sites for overall column
  group_by(hydro.endpoints, Target, Categories) %>% ## group by ffm and cats
  mutate(Tally = length(masterid)) %>% ## count sites for each ffm and cat
  mutate(PercentSites = Tally/TotalSites*100) %>% ## percentage per ffm and cat
  select(Tally, TotalSites, PercentSites) %>%
  mutate(ChannelType = "All") %>%
  distinct() %>%
  drop_na(ChannelType)

## add channel tally and overall tally together

tallyCats2c <- bind_rows(tallyCats1c, tallyCatsOVc)
unique(tallyCats2c$Target)  
unique(tallyCats2c$ChannelType)

## data for categories figure 
tallyCats2c <- tallyCats2c %>%
  mutate(ChannelType = factor(ChannelType, levels = c( "NAT", "SB", "HB", "All"),
                              labels = c("NAT & SB1", "SB0 & SB2", "HB", "All"))) %>%
  inner_join(labels, by = "hydro.endpoints") %>%
  mutate(Categories = factor(Categories, levels = c("No Action", "Non-flow related","Less than 10%","10 - 50%","Over 50%"))) %>%
  mutate(Target = factor(Target, levels = c("RefImprovement", "BOImprovement", "INCImprovement"),
                        labels = c("Reference Condition", "Best Observed Condition", "Incremental Improvement")))
tallyCats2c
## save
write.csv(tallyCats2c, "final_data/12_tally_categories_per_category_for_figures_V2.csv")

# Figures -----------------------------------------------------------------

head(tallyCats)
## plot all FFM per Target - bars rep improvement categories

tars <- unique(tallyCats$Target)

for(t in 1:length(tars)) {
  
  ## filter to target
  tallyCatsx <- tallyCats %>%
    filter(Target == tars[t])
  
  ## plot
  a1 <- ggplot(tallyCatsx, aes(fill=Categories, y=PercentSites, x=ChannelType)) + 
    geom_bar(position="stack", stat="identity") +
    # scale_fill_manual(values=c("chartreuse4", "dodgerblue1", "orange","pink", "firebrick3"))+ ## colour of points
    scale_fill_manual(values=hcl.colors(n=5, palette = "Zissou 1"))+ ## colour of points
    facet_wrap(~Flow.Metric.Name) +
    # scale_fill_manual(values=catPal)+
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Sites (%)")
  
  a1

  file.name1 <- paste0(out.dir, "12_", tars[t], "_target_mixed_mod_percent_sites.jpg")
  ggsave(a1, filename=file.name1, dpi=600, height=7, width=10)
  
}

## plot per categories - bars rep target achievable
cats <- unique(tallyCats2c$Categories)

for(c in 1:length(cats)) {
  
  ## filter to target
  data <- tallyCats2c %>%
    filter(Categories == cats[c])

  ## plot
  a1 <- ggplot(data, aes(fill=Target, y=PercentSites, x=ChannelType)) + 
    geom_bar(position="stack", stat="identity") +
    # scale_fill_manual(values=c("chartreuse4", "dodgerblue1", "orange","pink", "firebrick3"))+ ## colour of points
    scale_fill_manual(values=hcl.colors(n=3, palette = "Zissou 1"))+ ## colour of points
    facet_wrap(~Flow.Metric.Name) +
    # scale_fill_manual(values=catPal)+
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Sites (%)")
  a1
  
  file.name1 <- paste0(out.dir, "12_", cats[c], "_category_mixed_mod_percent_sites.jpg")
  ggsave(a1, filename=file.name1, dpi=600, height=7, width=12)
  
}


## spatial figures

head(categories2)

## get coordinates
head(AllData) ## full df

coords <- AllData %>%
  select(masterid, longitude, latitude) %>%
  distinct() ## select only site info 

### add to categories df

categoriesSP <- categories2 %>%
  inner_join(coords, by = "masterid") %>%
  st_as_sf(coords=c( "longitude", "latitude"), crs=4326, remove=F) %>% ## make spatial
  filter(!Categories %in% c("Type X", "Type X2")) %>% ## remove healthy sites
  drop_na(ChannelType) %>% ## remove NAs
  mutate(ChannelType = factor(ChannelType, levels = c("NAT", "SB1", "SB0", "SB2", "HB"),
                              labels = c("NAT & SB1", "NAT & SB1", "SB0 & SB2","SB0 & SB2", "HB"))) %>% ## make channel type a factor
  inner_join(labels, by = "hydro.endpoints") ## join with formal names
## upload ca boundary

ca_sf <- st_read("ignore/SpatialData/California/Ca_State_poly.shp")
## upload ca counties

counties_sf <- st_read("ignore/SpatialData/Counties/Counties.shp")

## get only socal
counties_socal_sf<-counties_sf %>%
  filter(NAME_PCASE %in% c("Ventura","Los Angeles", "Orange","San Bernardino", "Riverside", "San Diego"))

## upload watersheds

sheds_sf<- st_read("ignore/SpatialData/SMCSheds2009/SMCSheds2009.shp")


## create beige palet

beige_pal<-c("#f2dbb7","#eed9c4","#fff0db","#e4d5b7","#d9b99b",
             "#d9c2ba","#9c8481","#e2cbb0","#a69279","#f2dbb7",
             "#f6e6bf","#a69279","#f2dbb7","#9c8481","#e4d5b7")


# Plot per flow metric and result

metrics <- unique(categoriesSP$hydro.endpoints)
metrics
m=1

head(categoriesSP)
## loop
for(m in 1:length(metrics)) {
  
  ## extract metric
  CatsFFM <- categoriesSP %>%
    filter(hydro.endpoints == metrics[m])
  
  ## get metric fancy name
  metName <- CatsFFM$Flow.Metric.Name[1]
  
  CatsFFM

    ## plot
    m1 <- ggplot()+
      geom_sf(data=ca_sf)+ ## california
      geom_sf(data=sheds_sf, aes(fill=SMC_Name), color=NA)+ ## watersheds
      scale_fill_manual(values=beige_pal, guide= "none")+ ## beige colour for watersheds
      geom_sf(data=counties_socal_sf, fill=NA)+ ## counties
      geom_sf(data=CatsFFM,aes(colour = Categories), size = 2) + ## results
      # scale_colour_manual(values=c("chartreuse4", "dodgerblue2", "mediumpurple2", "firebrick3"))+
      scale_colour_manual(values=hcl.colors(n=6, palette = "Zissou 1"))+ ## colour of points
      coord_sf(xlim=c(-119.41,-116.4),
               ylim=c(32.5, 34.8),
               crs=4326) +
      theme_bw()+
      labs(title = paste(metName)) +
      theme(legend.title = element_blank(), 
            legend.position = "bottom",
            legend.text=element_text(size=15),
            plot.title = element_text(size = 15),
            axis.text = element_text(size = 8),
            strip.text.x = element_text(size = 15),
            strip.text.y = element_text(size = 15)) +
      # facet_grid(rows = vars(BioResult), cols = vars(FlowResult)) +
      facet_wrap(~ChannelType) +
      guides(col= guide_legend(title= ""))
    
    
    file.name1 <- paste0(out.dir, "12_", metrics[m], "_mixed_mod_maps_per_result_New.jpg")
    
    ggsave(m1, filename=file.name1, dpi=500, height=10, width=17)
    
  }
  

m1

# Comparing to Natural Flows ----------------------------------------------

## load change in delta
data <- read.csv("ignore/12_Change_in_delta_to_hit_thresholds_Mixed_Model.csv")
head(data)

# Type X Healthy
# Type A Unhealthy, flow not altered
# Type B Can reach ref DH within 10% 
# Type C Can reach BO with 10%
# Type D Can reach Inc (0.1) with 10% 
# Type E Can't achieve inc improvement with 10%

## get categories
categories <- data %>%
  mutate(
    Categories = case_when(
      HealthyCSCI == "Yes" ~ "Type X",
      NoFlowProblem == "Yes" & HealthyCSCI == "No" ~ "Type A",
      Threshold == "NAT" & RelChangeInDelta <= 10 ~ "Type B",
      ChannelType != "NAT" & RelChangeInDelta <= 10 ~ "Type C",
      RelChangeInDeltaINC <= 10 ~ "Type D",
      TRUE ~ "Type E"  # Assign NA if none of the conditions are met
    )
  )

# view(categories)
head(categories)

## if type X or A add deltah_final to change in delta - but maybe just NA as these don't matter in analysis
## if type B or C - add deltah_final to change in delta
## if type D - new delta
## if type E - NA
## then add on 50th percentile of natural flows
## and see if its between the 10th and 90th percentile

categoriesDF <- categories %>%
  select(masterid:ChannelType, Threshold, BioThresh, Categories, newDelta, ChangeInDelta)   %>%
  mutate(AdjustedDelta = ifelse(Categories == "Type D", newDelta, (deltah_final+ChangeInDelta))) %>%
  mutate(AdjustedDelta = ifelse(Categories %in% c("Type A", "Type X"), deltah_final, AdjustedDelta)) %>%
  # mutate(AdjustedDelta = ifelse(Categories == "Type E", NA, AdjustedDelta)) %>%
  select(-newDelta)

categoriesDF

write.csv(categoriesDF, "final_data/12_change_in_delta_current_delta_Mixed_Model.csv")

## upload natural flows from Megan

natFlows <- read.csv("input_data/change_in_delta_range.csv") %>%
  select(masterid, sampleyear, hydro.endpoints, comid:p90)

names(natFlows)
head(natFlows)

head(categoriesDF)

## join new categories to natural flows

catsNat <- full_join(categoriesDF, natFlows, by = c("masterid","sampleyear", "hydro.endpoints"), relationship = "many-to-many")
## the ones that did not match are Q99 - no nat flows so OK
# ind <- which(is.na(catsNat$comid))
# 
# test <- catsNat[ind,]
# length(unique(test$masterid))

## get current flow by adding delta to nat flow
## find whether it's in range

catsNat2 <- catsNat %>% 
  mutate(CurrentCFS = p50 + AdjustedDelta) %>% 
  mutate(WithinNaturalRange = if_else(CurrentCFS < p90 & CurrentCFS > p10, "yes", "no"))
head(catsNat2)

unique(counts$ChannelType)

write_csv(catsNat2, "ignore/12_currentCSF_nat_flows_categories.csv")

## count number sites within range for each Type
## then each type and wyt

counts <- catsNat2 %>%
  drop_na(ChannelType) %>%
  mutate(ChannelType = ifelse(ChannelType %in% c("SB0", "SB2"), "SB", ChannelType)) %>% ## add SB2 & 0 toegther
  mutate(ChannelType = ifelse(ChannelType %in% c("NAT", "SB1"), "NAT", ChannelType)) %>% ## add SB1 &NAT together
  filter(wyt == "all") %>% ## use only "all" wyt
  select(masterid, hydro.endpoints, ChannelType, Categories, WithinNaturalRange) %>% ## remove columns
  distinct() %>% ## remove duplicates
  group_by(hydro.endpoints, ChannelType, Categories, WithinNaturalRange) %>% tally() %>% ## count n sites per category
  ungroup(WithinNaturalRange) %>% 
  mutate(TotalPerCategory = sum(n)) %>% ## calculate total sites per category
  mutate(PercPerCategory = n/TotalPerCategory*100) %>% ## calculate percentage
  drop_na(Categories, WithinNaturalRange) %>% ## drop NAs Q99 no nat flows
  mutate(Categories = factor(Categories, levels = c("Type X", "Type A", "Type B", "Type C", "Type D", "Type E")))  %>%
  mutate(ChannelType = factor(ChannelType, levels = c("NAT", "SB", "HB"),
                            labels = c("NAT & SB1",  "SB0 & SB2", "HB"))) ## make channel type a factor
  
head(counts)

# sum(is.na(counts$WithinNaturalRange)) # 16
# ind <- which(is.na(counts$WithinNaturalRange))
# counts[ind,] ## all Q99 so no ref flows

## save
write.csv(counts, "final_data/12_n_sites_per_category_within_ref_flows.csv")

## plot stacked bar plot

ffms <- unique(counts$hydro.endpoints)

for(f in 1:length(ffms)) {
  
  countsx <- counts %>%
    filter(hydro.endpoints == ffms[f]) %>%
    filter(!Categories == "Type X")
  
  ## plot percentages
  p1 <- ggplot(countsx, aes(fill=WithinNaturalRange, y=n, x=Categories)) + 
    geom_bar(position="fill", stat="identity") +
    facet_wrap(~ChannelType, scales = "free_y") +
    theme_classic() +
    scale_y_continuous(name = "Number of Sites (%)") +
    scale_x_discrete(name = "") +
    scale_fill_discrete(name = "Within Natural Range?") +
    theme(legend.title = element_text(size=15), 
          legend.position = "bottom",
          legend.text=element_text(size=15),
          axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 15),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15))
  p1
  file.name1 <- paste0(out.dir, "12_", ffms[f], "_stacked_bar_n_sites_nat_flows_percent.jpg")
  ggsave(p1, filename=file.name1, dpi=300, height=7, width=9)
  
  ### plot N
  p2 <- ggplot(countsx, aes(fill=WithinNaturalRange, y=n, x=Categories)) + 
    geom_bar(position="stack", stat="identity") +
    facet_wrap(~ChannelType, scales = "free_y") +
    theme_classic() +
    scale_y_continuous(name = "Number of Sites") +
    scale_x_discrete(name = "") +
    scale_fill_discrete(name = "Within Natural Range?") +
    theme(legend.title = element_text(size=15), 
          legend.position = "bottom",
          legend.text=element_text(size=15),
          axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 15),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15))
  
  file.name1 <- paste0(out.dir, "12_", ffms[f], "_stacked_bar_n_sites_nat_flows_count.jpg")
  ggsave(p2, filename=file.name1, dpi=300, height=7, width=9)
  
  
}


## number of sites that need to go to zero flow. 

head(catsNat2)

## subset zero flow sites (based on <1 cfs 90th%)

zeroFlow <- catsNat2 %>%
  filter(p50 == 0) %>%
  mutate(Categories = factor(Categories, levels = c("Type X", "Type A", "Type B", "Type C", "Type D", "Type E")),
         ChannelType = factor(ChannelType, levels = c("NAT", "SB0", "SB1", "SB2", "HB"))) %>%
  select(masterid,comid,hydro.endpoints, deltah_final, ChannelType, Categories, p50) %>%
  distinct()
zeroFlow

## save for appendix
write.csv(zeroFlow, "Tables/12_appendix_sites_ref_zero.csv")

unique(zeroFlow$hydro.endpoints) ## "DS_Mag_50"      "Wet_BFL_Mag_10"
length(unique(zeroFlow$masterid)) ## 26
unique(zeroFlow$Categories) ## "Type X" "Type E" "Type C" "Type A"

## count number of sites 

cFlows <- zeroFlow %>%
  select(masterid, hydro.endpoints, ChannelType, Categories) %>% ## remove columns
  distinct() %>% ## remove duplicates
  group_by(hydro.endpoints, ChannelType, Categories) %>% tally() %>%
  inner_join(labels, by = "hydro.endpoints")

cFlows  

write.csv(cFlows, "Tables/12_num_sites_to_zero.csv")


# count sites within SD of delta 0 ----------------------------------------

## load change in delta
data <- read.csv("ignore/12_Change_in_delta_to_hit_thresholds_Mixed_Model.csv")
head(data)

## get categories
categories <- ChangeDxWide %>%
  mutate(
    Categories = case_when(
      HealthyCSCI == "Yes" ~ "Type X",
      HealthyCSCIBO == "Yes" ~ "Type X2",
      WithinLims == "Yes" & HealthyCSCI == "No" ~ "Type A",
      WithinLimsBO == "Yes" & HealthyCSCIBO == "No" ~ "Type A2",
      Threshold == "NAT" & RelChangeInDelta <= 10 ~ "Type B",
      !(ChannelType %in% c("NAT", "SB1")) & RelChangeInDelta <= 10 ~ "Type C",
      RelChangeInDeltaINC <= 10 ~ "Type D",
      TRUE ~ "Type E"  # Assign NA if none of the conditions are met
    )
  )


# view(categories)
head(categories)


## add change in delta to delta H final

absDelta <- categories %>%
  mutate(NewDeltaAll = ifelse(Categories == "Type E", newDelta, newDeltaBO))
  

head(absDelta)
## check that only channel types are counted on their own thresholds
