## combine table for data packet
## sites code, scores, ffm value and stress levels, number of strikes

library(tidyverse)
library(sf)
library(tidylog)

# Upload data -------------------------------------------------------------

## within limits doc
imps <- read.csv("ignore/05_impact_ffm_bio.csv") %>%
  select( -X)

## for now remove natlow,med,high
removes <- unique(imps$Threshold)[c(1,5,7)]
imps <- imps %>%
  filter(!Threshold %in% removes)

imps

# Ranges ------------------------------------------------------------------

## upload, remove asci and stanard thresholds
ranges <- read.csv("final_data/05_ffm_ranges.csv") %>%
  filter(!Threshold %in% removes, Index == "csci") %>%
  dplyr::select(-c(X, masterid, Index)) %>%
  distinct() 

## change min/max sub values to associated min/max of metric

# Flow.Metric.Name                  `max(deltah_final)` `min(deltah_final)`
# <chr>                                           <dbl>               <dbl>
#   1 10-year flood magnitude                        -1.55             -19120. 
# 2 2-year flood magnitude                         -5.01              -7944. 
# 3 5-year flood magnitude                         -1.54             -18874. 
# 4 Dry-season median baseflow                     94.6                 -26.9
# 5 Fall pulse magnitude                          922.                 -227. 
# 6 Magnitude of largest annual storm              -0.701            -44934. 
# 7 Spring recession magnitude                   1700.                -1050. 
# 8 Wet-season low baseflow                        74.6                 -89.2
# 9 Wet-season median baseflow                    205.                 -368. 

## define values we want to swap
minVal <- (ranges$Lower2)[1]
maxVal <- (ranges$Upper2)[1]

ranges2 <- ranges %>%
  mutate(Lower = case_when(Lower2 == minVal & Flow.Metric.Name == "Spring recession magnitude" ~ -1050.,
                            Lower2 == minVal & Flow.Metric.Name == "10-year flood magnitude" ~ -19120. ,
                            Lower2 == minVal & Flow.Metric.Name == "2-year flood magnitude" ~ -7944.,
                            Lower2 == minVal & Flow.Metric.Name == "5-year flood magnitude" ~ -18874,
                            Lower2 == minVal & Flow.Metric.Name == "Dry-season median baseflow" ~ -26.9,
                            Lower2 == minVal & Flow.Metric.Name == "Fall pulse magnitude" ~ -227.,
                            # Lower2 == minVal & Flow.Metric.Name == "Magnitude of largest annual storm" ~ 1700.,
                            Lower2 == minVal & Flow.Metric.Name == "Wet-season low baseflow" ~ -89.2,
                            Lower2 == minVal & Flow.Metric.Name == "Wet-season median baseflow" ~  -368.)) %>% 
  mutate(Lower = ifelse(is.na(Lower), Lower2, Lower)) %>%
  mutate(Upper = case_when(Upper2 == maxVal & Flow.Metric.Name == "Spring recession magnitude" ~ 1700,
                            Upper2 == maxVal & Flow.Metric.Name == "10-year flood magnitude" ~ 0 ,
                            Upper2 == maxVal & Flow.Metric.Name == "2-year flood magnitude" ~ 0,
                            Upper2 == maxVal & Flow.Metric.Name == "5-year flood magnitude" ~ 0,
                            Upper2 == maxVal & Flow.Metric.Name == "Dry-season median baseflow" ~ 94.6,
                            Upper2 == maxVal & Flow.Metric.Name == "Fall pulse magnitude" ~ 922.,
                            # Upper2 == maxVal & Flow.Metric.Name == "Magnitude of largest annual storm" ~ 1700.,
                            Upper2 == maxVal & Flow.Metric.Name == "Wet-season low baseflow" ~ 74.6,
                            Upper2 == maxVal & Flow.Metric.Name == "Wet-season median baseflow" ~  205)) %>%
  mutate(Upper = ifelse(is.na(Upper), Upper2, Upper)) %>%
  dplyr::select(-Lower2, - Upper2)

ranges
## save ranges 
write.csv(ranges, "final_data/07_ffm_ranges_modified_csci.csv")

## combine upper lower columns and convert to text. make wide based on quantiles
ranges_wide <- ranges2 %>%
  mutate(Lower = round(Lower, digits = 0)) %>% ## round range values
  mutate(Upper = round(Upper, digits = 0)) %>%
  mutate(FFMRange = paste0(Lower, " to ", Upper)) %>% ## paste together 
  dplyr::select(-c(Lower, Upper)) %>% ## remove non-combined limits
  pivot_wider(names_from= Quantile, values_from = FFMRange) %>% ## make wide
  dplyr::select(Threshold:Flow.Metric.Name, "0.5", "0.9") %>% ## change order to be more logical
  rename(ModifiedChannelType = Threshold, CSCI.Threshold = BioThresh) %>% ## rename 
  mutate(ModifiedChannelType = factor(ModifiedChannelType, levels = c("NAT", "SB0", "SB2", "HB")))

write.csv(ranges_wide, "final_data/07_ffm_ranges_modified_csci_wide.csv")


# Format ------------------------------------------------------------------

## format names
impsx <- imps %>% 
  dplyr::select(Index, Hydro_endpoint, Threshold, BioThresh,  masterid, COMID, Flow.Metric.Name, Flow.Component, Result, longitude, latitude)  %>%
  mutate(Threshold = factor(Threshold, levels = c("NAT", "SB0", "SB2", "HB"), 
                            labels = c("Natural", "Soft Bottom (0)" , "Soft Bottom (2)", "Hard Bottom"))) %>%
  mutate(BioResult = case_when(Result %in% c("HUF", "HMF","HLF") ~ "Healthy Biology",
                               Result %in% c("UHUF", "UHMF", "UHLF") ~ "Unhealthy Biology")) %>%
  mutate(FlowResult = case_when(Result %in% c("HUF", "UHUF") ~ "Unlikely Stressed",
                                Result %in% c("HMF", "UHMF") ~ "Likely Stressed",
                                Result %in% c("HLF", "UHLF") ~ "Likely Very Stressed")) %>%
  mutate(FlowResult = factor(FlowResult, levels = c("Unlikely Stressed", "Likely Stressed", "Likely Very Stressed"))) %>%
  mutate(Result = factor(Result, levels = c("HUF", "UHUF", "HMF", "UHMF", "HLF", "UHLF"),
                         labels = c("Healthy Biology, Unlikely Stressed",
                                    "Unhealthy Biology, Unlikely Stressed",
                                    "Healthy Biology, Likely Stressed",
                                    "Unhealthy Biology, Likely Stressed",
                                    "Healthy Biology, Likely Very Stressed",
                                    "Unhealthy Biology, Likely Very Stressed"))) 



# Number of strikes -------------------------------------------------------
removes
## upload number of ffm strikes 
strikes <- read.csv("final_data/05_Number_ffm_per_result.csv") %>%
  filter(!Threshold %in% removes) %>%
  mutate(Threshold = factor(Threshold, levels = c("NAT", "SB0", "SB2", "HB"), 
                            labels = c("Natural", "Soft Bottom (0)" , "Soft Bottom (2)", "Hard Bottom")), 
         NumStrikes = as.factor(n)) %>%
  mutate(BioResult = case_when(Result %in% c("HUF", "HMF","HLF") ~ "Healthy Biology",
                               Result %in% c("UHUF", "UHMF", "UHLF") ~ "Unhealthy Biology")) %>%
  mutate(FlowResult = case_when(Result %in% c("HUF", "UHUF") ~ "Unlikely Stressed",
                                Result %in% c("HMF", "UHMF") ~ "Likely Stressed",
                                Result %in% c("HLF", "UHLF") ~ "Very Likely Stressed")) %>%
  mutate(FlowResult = factor(FlowResult, levels = c("Unlikely Stressed", "Likely Stressed", "Very Likely Stressed"))) %>%
  mutate(Result = factor(Result, levels = c("HUF", "UHUF", "HMF", "UHMF", "HLF", "UHLF"),
                         labels = c("Healthy Biology, Unlikely Stressed",
                                    "Unhealthy Biology, Unlikely Stressed",
                                    "Healthy Biology, Likely Stressed",
                                    "Unhealthy Biology, Likely Stressed",
                                    "Healthy Biology,  Very Likely Stressed",
                                    "Unhealthy Biology,  Very Likely Stressed"))) 

# FFm and bio values ------------------------------------------------------

load(file = "final_data/03_bugs_algae_flow_joined_by_masterid_noPeakAugs.RData") 
AllDataLongx

## get sites with flow data
AllData <- AllDataLongx %>%
  dplyr::select(-c( longitude, latitude, COMID, comid, flow_metric, Class2, hydro.endpoints, Flow.Component)) %>% ## remove redundant columns
  distinct() %>%
  drop_na(deltah_final) %>% ## remove sites with no FFM
  rename(Metric.Score = MetricValue, Modified.Class = channel_engineering_class, 
         Flow.Metric.Value = deltah_final) %>% ## humanise names
  pivot_wider(names_from = Flow.Metric.Name, values_from = Flow.Metric.Value) %>% ## format wide
  filter(!Metric == "asci") %>% ## remove asci
  drop_na(Modified.Class) ## remove sites with no classification
  
## take max of scores with more than one value
AllData2 <- AllData %>%
  group_by(masterid, sampleyear) %>% ## group by site and year
  mutate(Metric.Score2 = max(Metric.Score)) %>% ## get max score
  dplyr::select(-Metric.Score) %>% ## remove orriginal score
  distinct() %>% ## distinct with new scores
  dplyr::select(masterid:sampleyear, Metric.Score2, Modified.Class, 
                "Dry-season median baseflow":"Magnitude of largest annual storm") %>% ## change order of columns
  rename(CSCI.Score = Metric.Score2) #%>% ## rename metric score
  # pivot_longer("Dry-season median baseflow":"Magnitude of largest annual storm", names_to = "Flow.Metric.Name", values_to = "Flow.Metric.Value")


head(AllData2)
names(AllData2)
## get list of sites

sites <- unique(AllData2$masterid)
sites

## upload bio sites shape to get spatial info

bio_sf <- st_read("ignore/01_bio_sites_all.shp")
bio_sf

## subset all sites to sites in model

bio_sf_sub <- bio_sf %>%
  filter(masterid %in% sites) %>%
  as.data.frame() %>% ## make into df
  dplyr::select(-geometry)


## join with smc sites

ffm_bio_tab <- full_join(bio_sf_sub, AllData2) %>%
  mutate(Threshold = factor(Modified.Class, levels = c("NAT", "SB0", "SB2", "HB"), 
                            labels = c("Natural", "Soft Bottom (0)" , "Soft Bottom (2)", "Hard Bottom"))) %>% ## change level names of modified class
  dplyr::select(-Modified.Class) %>%
  pivot_longer("Dry-season median baseflow":"Magnitude of largest annual storm", names_to = "Flow.Metric.Name", values_to = "Change_in_Flow.Metric_Value_(cfs)")

head(ffm_bio_tab)


# Add flow alteration stress ----------------------------------------------

head(impsx)
## format data for join
imps1 <- impsx %>%
  filter(Index == "csci") %>%
  dplyr::select(-c(Index, Hydro_endpoint, COMID, Flow.Component, longitude, latitude, Result)) %>% ## remove redundant columns
  # mutate(Flow.Metric.Name = paste0(Flow.Metric.Name, " Stress Level")) %>%
  distinct() #%>%
  # pivot_wider(names_from = Flow.Metric.Name, values_from = FlowResult)

names(imps1)
head(imps1)
head(ffm_bio_tab)

## join with ffm and bio data

ffm_bio_stress_tab <- full_join(ffm_bio_tab, imps1, by = c("Threshold", "masterid", "Flow.Metric.Name"))

sum(is.na(ffm_bio_stress_tab$masterid)) ## 0 - NAs are from sites we dpn't have mod channel data, can remove for now

# Format column order -----------------------------------------------------

names(ffm_bio_stress_tab)

infoTable <- ffm_bio_stress_tab %>%
  # dplyr::select(masterid:CSCI.Score, Threshold:BioThresh, BioResult,F) %>%
  distinct() %>%
  rename(Modified.Class = Threshold, CSCI.Threshold = BioThresh, CSCI.Level = BioResult, Flow.Alteration.Stress = FlowResult) %>%
  group_by(masterid) %>%
  mutate(RecentYear = max(sampleyear)) %>%
  mutate(YearKeep = ifelse(sampleyear == RecentYear, "Yes", "No")) %>%
  filter(YearKeep == "Yes") %>%
  dplyr::select(-c(RecentYear, YearKeep, COMID, COMID1)) %>%
  drop_na(Modified.Class, "Change_in_Flow.Metric_Value_(cfs)")

names(infoTable)
head(infoTable)

write.csv(infoTable, "final_data/07_ffm_bio_stresLevel_data_packet_all_counties.csv")

length(unique(infoTable$masterid)) ## 132

unique(infoTable$county)

# Add in Ranges -----------------------------------------------------------

## format to match
ranges_wide <- ranges_wide %>%
  mutate(Modified.Class = factor(ModifiedChannelType, levels = c("NAT", "SB0", "SB2", "HB"), 
                            labels = c("Natural", "Soft Bottom (0)" , "Soft Bottom (2)", "Hard Bottom"))) %>%
  dplyr::select(-ModifiedChannelType)

names(ranges_wide)
## join to clean table

infoTablex <- full_join(infoTable, ranges_wide, by = c("CSCI.Threshold", "Flow.Metric.Name", "Modified.Class")) %>%
  dplyr::select(masterid:CSCI.Level, "0.5", "0.9", "Flow.Alteration.Stress")


names(infoTablex)

## save
write.csv(infoTablex, "final_data/07_ffm_bio_stresLevel_data_packet_all_counties_with_ranges.csv")

# Separate by county ------------------------------------------------------

## san diego
infoTableSD <- infoTablex %>%
  filter(county == "San Diego")

write.csv(infoTableSD, "/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/BioIntegrity in Modified Channels Statewide - SMC/Flow Ecology in Modified Channels/DataPackets/SanDiego_flowEcology_modChannels_stress.csv")

## ventura

infoTableVen <- infoTablex %>%
  filter(county == "Ventura")

write.csv(infoTableVen, "/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/BioIntegrity in Modified Channels Statewide - SMC/Flow Ecology in Modified Channels/DataPackets/Ventura_flowEcology_modChannels_stress.csv")

## LA

infoTableLA <- infoTablex %>%
  filter(county == "Los Angeles")

write.csv(infoTableLA, "/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/BioIntegrity in Modified Channels Statewide - SMC/Flow Ecology in Modified Channels/DataPackets/Los_Angeles_flowEcology_modChannels_stress.csv")

## San Bernardino

infoTableSB <- infoTablex %>%
  filter(county == "San Bernardino")

write.csv(infoTableSB, "/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/BioIntegrity in Modified Channels Statewide - SMC/Flow Ecology in Modified Channels/DataPackets/San_Bernardino_flowEcology_modChannels_stress.csv")

## Ornage
infoTableOC <- infoTablex %>%
  filter(county == "Orange")

write.csv(infoTableOC, "/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/BioIntegrity in Modified Channels Statewide - SMC/Flow Ecology in Modified Channels/DataPackets/Orange_flowEcology_modChannels_stress.csv")

## Riverside

infoTableRS <- infoTablex %>%
  filter(county == "Riverside")

write.csv(infoTableRS, "/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/BioIntegrity in Modified Channels Statewide - SMC/Flow Ecology in Modified Channels/DataPackets/Riverside_flowEcology_modChannels_stress.csv")


# format number of strikes ---------------------------------------------------
unique(strikes$Result)

## strikes have to be separate due to the formatting of the data
strikes1 <- strikes %>%
  filter(Index == "csci") %>% ## remove asci
  dplyr::select(masterid, Threshold, Result, NumStrikes) %>% ## remove redundant columns
  mutate(NumStrikes = as.numeric(NumStrikes)) %>%
  pivot_wider(names_from = Result, values_from = NumStrikes) %>%
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>%
  mutate("Total Stressed FFM" = `Unhealthy Biology, Likely Stressed`+ `Unhealthy Biology,  Very Likely Stressed`) %>%
  rename(Modified.Class = Threshold) %>%
  dplyr::select(masterid, Modified.Class, "Total Stressed FFM")

## get all site info from above
infoTableSub <- infoTablex %>%
  dplyr::select(masterid:Modified.Class)

## join with strikes
strikes2 <- inner_join(infoTableSub, strikes1, by = c("masterid", "Modified.Class"))

strikes2 <- strikes2 %>%
  group_by(masterid) %>%
  mutate(RecentYear = max(sampleyear)) %>%
  mutate(YearKeep = ifelse(sampleyear == RecentYear, "Yes", "No")) %>%
  filter(YearKeep == "Yes") %>%
  dplyr::select(-c(RecentYear, YearKeep))


## save to repo
write.csv(strikes2, "final_data/07_sum_strikes_all_sites_unhealthy_bio.csv")  

## save to teams
write.csv(strikes2, "/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/BioIntegrity in Modified Channels Statewide - SMC/Flow Ecology in Modified Channels/DataPackets/Number_FFM_per_Site.csv")

