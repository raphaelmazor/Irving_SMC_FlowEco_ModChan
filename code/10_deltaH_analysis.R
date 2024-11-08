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

## directory for figures
out.dir <- "final_figures/"

# Upload GAM Model -----------------------------------------------------------

### Median Gam

load(file = "ignore/04_csci_asci_quantile_gam_updatedSites.RData")
gam.lm

## look up table 
load( file="final_data/04_csci_asci_quantile_gam_rsqds_updatedSites.RData")
head(bio_h_summary)

## smoothing functions per FFM
sm6csci <- c("DS_Mag_50", "FA_Mag", "Peak_10", "Peak_2", "Peak_5", "Wet_BFL_Mag_50", "Wet_BFL_Mag_10")

sm3csci <- c( "SP_Mag",  "Q99")


## add row number as a column for an index and change FFM names

bio_h_summary <- bio_h_summary %>%
  mutate(Index = row.names(bio_h_summary)) %>%
  mutate(Hydro_endpoint = case_when(flow.endpoints == "d_ds_mag_50" ~ "DS_Mag_50",            
                                    flow.endpoints == "d_fa_mag" ~ "FA_Mag",
                                    flow.endpoints == "d_peak_10" ~ "Peak_10",
                                    flow.endpoints == "d_peak_2" ~ "Peak_2",
                                    flow.endpoints == "d_peak_5" ~ "Peak_5",
                                    flow.endpoints == "d_sp_mag" ~ "SP_Mag",
                                    flow.endpoints == "d_wet_bfl_mag_10" ~ "Wet_BFL_Mag_10",
                                    flow.endpoints == "d_wet_bfl_mag_50" ~ "Wet_BFL_Mag_50", 
                                    flow.endpoints == "delta_q99" ~ "Q99"))

## filter to median (0.5) and smoothing functions per ffm
Med_GAM_models <- bio_h_summary %>%
  mutate(ChosenSmth = case_when((biol.endpoints == "csci" & Hydro_endpoint %in% sm6csci) ~ 6,
                                (biol.endpoints == "csci" & Hydro_endpoint %in% sm3csci) ~ 3)) %>%
  mutate(ModelstoUse = ifelse(smooth_funcs == ChosenSmth, "Yes", "No")) %>%
  filter(ModelstoUse == "Yes") %>%
  filter(quants %in% c( 0.50))

Med_GAM_models

## upload thresholds

df <- read.csv("final_data/05_QGAM_FFM_Ranges_With_Adjustments.csv") %>%
  filter(Quantile == 0.5, Threshold %in%c("NAT", "HB", "SB2", "SB0")) %>%
  select(Hydro_endpoint, Threshold, BioThresh, Lower2, Upper2)
df

## function to get delta H at csci score
find_delta_h_optimize <- function(target_csci_score, model, lower, upper, observed_delta_h) {
  # Define the objective function with a secondary criterion
  objective_function <- function(deltah_final) {
    predicted_csci <- predict(model, newdata = data.frame(deltah_final = deltah_final))
    diff_csci <- abs(predicted_csci - target_csci_score)  # Primary objective: absolute difference from target
    
    # Secondary criterion: distance from observed_delta_h
    diff_reference <- abs(deltah_final -observed_delta_h)
    
    # Combine objectives: prioritize csci match, then proximity to reference point
    diff_csci + 0.001 * diff_reference  # Small weight on reference distance
  }
  
  # Optimize for negative delta_h range
  result_negative <- optimize(objective_function, lower = lower, upper = 0)
  delta_h_negative <- result_negative$minimum
  
  # Split the positive range and optimize each part
  midpoint <- (upper + 0) / 2
  result_positive_1 <- optimize(objective_function, lower = 0, upper = midpoint)
  result_positive_2 <- optimize(objective_function, lower = midpoint, upper = upper)
  
  # Choose the best result for positive delta_h based on objective values
  if (result_positive_1$objective < result_positive_2$objective) {
    delta_h_positive <- result_positive_1$minimum
  } else {
    delta_h_positive <- result_positive_2$minimum
  }
  
  # Return both results
  list(delta_h_negative = delta_h_negative, delta_h_positive = delta_h_positive)
}

# Example usage with dynamic lower and upper
target_csci_score <- 0.5  # target csci_score

predicted_delta_h_values <- find_delta_h_optimize(
  target_csci_score, mod, 
  min(FFMDf$deltah_final), max(FFMDf$deltah_final), 
  observed_delta_h = 0  # observed delta H
)


print(predicted_delta_h_values)
head(FFMDf)


## save 
save(find_delta_h_optimize, file = "code/functions/10_find_deltaH.Rdata")

# Upload Site data ---------------------------------------------------

## data
load(file = "final_data/01_bugs_algae_flow_joined_by_masterid.RData") 

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
  mod <- gam.lm[[Ind]]
  
  ## find DH for incremental value - 0.1 increase of CSCI
  
  ## define increment values
  FFMDf <- FFMDf %>%
    mutate(csci_Incr = MetricValue+0.1)

  ## define sites/rownames - adjust here to remove duplicate sites
  sites <- rownames(FFMDf)

  # loop over sites
  
  for(s in 1:length(sites)) {
    
    siteDelta <- FFMDf[s,"deltah_final"] ## get delta
    siteIncr <- FFMDf[s,"csci_Incr"] ## get adjusted csci

    target_csci_score <-  siteIncr  #  target csci_score - observed csci increased by 0.1
    
    ## use function to predict delta H for both neg and pos values
    predicted_delta_h <- find_delta_h_optimize(target_csci_score, mod, min(FFMDf$deltah_final), max(FFMDf$deltah_final), siteDelta)
    
    ## add to main df
    FFMDf$newDelta[s] <- ifelse(siteDelta > 0, predicted_delta_h$delta_h_positive, predicted_delta_h$delta_h_negative)

  }
  
  FFMDf <- FFMDf %>%
    mutate(ChangeInDeltaINC = deltah_final - newDelta) %>% ## absolute change in delta for inc 
    mutate(ChangeInDeltaINC = round(ChangeInDeltaINC, digits = 2)) %>%
    mutate(RelChangeInDeltaINC = (abs(ChangeInDeltaINC)/abs(deltah_final))*100) %>% ## relative change
    mutate(RelChangeInDeltaINC = round(RelChangeInDeltaINC, digits = 2))

  ## delta H will either be within the limits or out.
  ## if it is within, then issue is likely not flow
  ## if it is outside then can calculate change in flow needed for an improved score
  
  ## first find out whether DH is within or outside limits
  
  head(df)
  head(FFMDf)
  ## join site to ffm limits data by ffm and threshold
  
  FFMDFLims <- right_join(FFMDf, df, by =c("hydro.endpoints" = "Hydro_endpoint")) %>%
    rename(ChannelType = channel_engineering_class)  
  
  
  head(FFMDFLims) ## missing sites are SB0, these can be changed to SB2 for now. 
  
  ## add column of whether within limits or not 
  
  FFMDFLims <- FFMDFLims %>%
    mutate(NoFlowProblem = ifelse(deltah_final >= Lower2 & deltah_final <= Upper2, "Yes", "No" )) %>%
    mutate(HealthyCSCI = ifelse(MetricValue >= 0.79, "Yes", "No")) ## add column for healthy csci (remove in figure later)
  
  
  # Delta H to thresholds ---------------------------------------------------
  
  head(FFMDFLims)
  
  ## remove sites within limits, to join back later
  # WithinLims <- FFMDFLims %>%
  #   filter(WithinLimits == "Yes") %>%
  #   mutate(ChangeInDelta = NA, RelChangeInDelta = NA) ## add change columns to match following data
  
  ## remove them from other df
  OutsideLims <- FFMDFLims #%>%
    # filter(WithinLimits == "No")
  
  head(OutsideLims)
  
  ## if delta H is positive it must be above the upper limit
  ## calculate difference of delta at site to upper limit
  ## if negative then calculate difference to lower limit
  
  OutsideLims2 <- OutsideLims %>%
    mutate(ChangeInDelta = ifelse(deltah_final >= 0, deltah_final- Upper2, Lower2 - deltah_final )) %>% ## absolute change
    mutate(ChangeInDelta = ifelse(NoFlowProblem == "Yes", NA, ChangeInDelta)) %>%
    mutate(RelChangeInDelta = (abs(ChangeInDelta)/abs(deltah_final))*100) ## relative change
  
  # names(WithinLims)
  # names(OutsideLims2)
  
  ChangeD <- OutsideLims2
  # head(ChangeD)
  
  
  ChangeDx <- rbind(ChangeDx,ChangeD)
  
}

head(ChangeDx)

## save
write.csv(ChangeDx, "final_data/10_Change_in_delta_to_hit_thresholds.csv")

## plot 

p1 <- ggplot(data = ChangeDx, aes(x=RelChangeInDelta, y= MetricValue, group = ChannelType, col = ChannelType)) +
  geom_smooth()+
  facet_wrap(~hydro.endpoints, scales = "free")

p1 ## plot of rel change vs csci score by channel type and ffm facet

p2 <- ggplot(data = ChangeDx, aes(x=RelChangeInDelta, y= MetricValue, group = hydro.endpoints, col = hydro.endpoints)) +
  geom_smooth()+
  facet_wrap(~ChannelType, scales = "free")

p2 ## plot of rel change vs csci score per ffm and channel type facet


# Calculate categories ----------------------------------------------------

head(ChangeDx)
## also filter dataset, only need the rows where channel type matches threshold
## relative change
ChangeDxWide <- ChangeDx %>%
  select(-c(Lower2, Upper2, deltah_final, ChangeInDelta,ChangeInDeltaINC,newDelta, BioThresh)) %>% ##  
  mutate(Match = ifelse(ChannelType == Threshold, "Yes", "No")) %>%
  filter(Match == "Yes") %>%
  distinct() #%>%
  # pivot_wider(names_from = Threshold, values_from = RelChangeInDelta)

## remember to remove duplicates/take max or mean of same year samples!!!!!!!

## from excel file
# Type A	Can reach ref threshold
# Type B	Can't reach ref but can reach BO
# Type C	Can't reach BO but can get inc imrpveoment
# Type D	Can't get inc improvement

## modified for now

# Type X Healthy
# Type A Unhealthy, flow not altered
# Type B Can reach ref DH within 10% 
# Type C Can reach BO with 10%
# Type D Can reach Inc (0.1) with 10% 
# Type E Can't achieve inc improvement with 10%


## 
## unhealthy

## can reach defined as 10% relative change - can change later!!! discuss with Rafi

## function to get the categories
head(ChangeDxWide)
## Healthy CSCI

categories <- ChangeDxWide %>%
  mutate(
    Categories = case_when(
      HealthyCSCI == "Yes" ~ "Type X",
      NoFlowProblem == "Yes" & HealthyCSCI == "No" ~ "Type A",
      ChannelType == "NAT" & RelChangeInDelta <= 10 ~ "Type B",
      ChannelType != "NAT" & RelChangeInDelta <= 10 ~ "Type C",
      RelChangeInDeltaINC <= 10 ~ "Type D",
      TRUE ~ "Type E"  # Assign NA if none of the conditions are met
    )
  )

view(categoriesx)
head(categories)
str(ChangeDxWide)

## 	405M10598
## 	801CYC121

## remove values from same site - if from same year max(), if from different years then most recent? - can change
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

sum(is.na(categories$Categories))

## how many csci scores per masterid
# categories %>% group_by(masterid, hydro.endpoints) %>% tally(MetricValue)

## count catergories per channel type
names(categories)
length(unique(categories$masterid)) ## 376

tallyCats <- categoriesx %>%
  filter(!Categories == "Type X") %>%
  group_by(ChannelType, hydro.endpoints) %>%
  mutate(TotalSites = length(unique(masterid))) %>%
  group_by(ChannelType, hydro.endpoints, Categories) %>%
  mutate(Tally = length(unique(masterid))) %>%
  mutate(PercentSites = Tally/TotalSites*100) %>%
  select(Tally, TotalSites, PercentSites) %>%
  distinct() %>%
  drop_na(ChannelType)


## ## remove healthy sites for figure
tallyCats <- tallyCats %>%
  mutate(ChannelType = factor(ChannelType, levels = c("NAT", "SB0", "SB1", "SB2", "HB"))) %>%
  inner_join(labels, by = "hydro.endpoints")
  # filter(!Categories == "Type X")

## plot

a1 <- ggplot(tallyCats, aes(fill=Categories, y=PercentSites, x=ChannelType)) + 
  geom_bar(position="stack", stat="identity") +
  # scale_fill_manual(values=c("chartreuse4", "dodgerblue1", "orange","pink", "firebrick3"))+ ## colour of points
  scale_fill_manual(values=hcl.colors(n=5, palette = "Zissou 1"))+ ## colour of points
  facet_wrap(~Flow.Metric.Name) +
  # scale_fill_manual(values=catPal)+
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Sites (%)")

a1


file.name1 <- paste0(out.dir, "10_percent_sites_per_cat_bar.jpg")
ggsave(a1, filename=file.name1, dpi=600, height=7, width=10)


# Upload Mixed GAM Model --------------------------------------------------

## data
load(file = "final_data/01_bugs_algae_flow_joined_by_masterid.RData") 

## remove SB1
AllDataLongx <- AllData %>%
  dplyr::filter(!channel_engineering_class == "SB1")
head(AllDataLongx)




