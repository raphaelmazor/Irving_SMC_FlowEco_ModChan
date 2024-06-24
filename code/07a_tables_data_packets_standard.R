## combine table for data packet
## sites code, scores, ffm value and stress levels, number of strikes

library(tidyverse)
library(sf)
library(tidylog)

# Upload data -------------------------------------------------------------

## within limits doc
imps <- read.csv("ignore/05_impact_ffm_bio.csv") %>%
  select( -X)

## remove all but nat med

imps <- imps %>%
  filter(Threshold == "NATMed")

## add eng data
## upload
BioEng <- read.csv("ignore/02_chan_eng.csv") %>% ## upload data
  select(-c(X,channel_engineering_classification_date, channel_engineering_personnel, channel_engineering_comments, Class2)) 

## join to results and format to match
imps2<- left_join(imps, BioEng, by = "masterid") %>%
  dplyr::select(-Threshold) %>%
  rename(Threshold = channel_engineering_class)

imps2

# Format ------------------------------------------------------------------

## format names
impsx <- imps2 %>% 
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

unique(impsx$masterid) %in% sites

# Number of strikes -------------------------------------------------------
## upload number of ffm strikes 
strikes <- read.csv("output_data/05_Number_ffm_per_result.csv") %>%
  filter(Threshold == "NATMed", Index =="csci")  %>%
  dplyr::select(-X)

## add eng data 

## upload
BioEng <- read.csv("ignore/02_chan_eng.csv") %>% ## upload data
  select(-c(X,channel_engineering_classification_date, channel_engineering_personnel, channel_engineering_comments, Class2)) 

## join to results and format to match
strikesx<- left_join(strikes, BioEng, by = "masterid") %>%
  dplyr::select(-Threshold) %>%
  rename(Threshold = channel_engineering_class) %>%
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

load(file = "final_data/01_bugs_algae_flow_joined_by_masterid.RData") 

## get sites with flow data
AllData <- AllData %>%
  dplyr::select(-c(X, longitude, latitude, COMID, comid, flow_metric, Class2, hydro.endpoints, Flow.Component)) %>% ## remove redundant columns
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
  dplyr::select(-Modified.Class)

head(ffm_bio_tab)


# Add flow alteration stress ----------------------------------------------

head(imps2)

## format data for join
imps1 <- impsx %>%
  filter(Index == "csci") %>%
  dplyr::select(-c(Index, Hydro_endpoint, COMID, Flow.Component, longitude, latitude, Result)) %>% ## remove redundant columns
  mutate(Flow.Metric.Name = paste0(Flow.Metric.Name, " Stress Level")) %>%
  distinct() %>%
  pivot_wider(names_from = Flow.Metric.Name, values_from = FlowResult)

names(imps1)
names(ffm_bio_tab)

unique(imps2$masterid) %in% unique(ffm_bio_tab)

## join with ffm and bio data

ffm_bio_stress_tab <- full_join(ffm_bio_tab, imps1, by = c("Threshold", "masterid"))


# Format column order -----------------------------------------------------

infoTable <- ffm_bio_stress_tab %>%
  dplyr::select(masterid:CSCI.Score, Threshold:BioThresh, BioResult,
                "Spring recession magnitude", "Spring recession magnitude Stress Level",
                "Dry-season median baseflow", "Dry-season median baseflow Stress Level",
                "Fall pulse magnitude", "Fall pulse magnitude Stress Level",
                "Wet-season low baseflow", "Wet-season low baseflow Stress Level",
                "Wet-season median baseflow", "Wet-season median baseflow Stress Level",
                "Magnitude of largest annual storm", "Magnitude of largest annual storm Stress Level",
                "2-year flood magnitude", "2-year flood magnitude Stress Level",
                "5-year flood magnitude", "5-year flood magnitude Stress Level",
                "10-year flood magnitude", "10-year flood magnitude Stress Level") %>%
  distinct() %>%
  rename(Modified.Class = Threshold, CSCI.Threshold = BioThresh, CSCI.Level = BioResult) %>%
  group_by(masterid) %>%
  mutate(RecentYear = max(sampleyear)) %>%
  mutate(YearKeep = ifelse(sampleyear == RecentYear, "Yes", "No")) %>%
  filter(YearKeep == "Yes") %>%
  dplyr::select(-c(RecentYear, YearKeep, COMID, COMID1)) %>%
  drop_na(Modified.Class)

names(infoTable)

write.csv(infoTable, "final_data/07_ffm_bio_stresLevel_data_packet_all_counties_standard.csv")

length(unique(infoTable$masterid)) ## 132
unique(infoTable$county)

# Separate by county ------------------------------------------------------

## san diego
infoTableSD <- infoTable %>%
  filter(county == "San Diego")

write.csv(infoTableSD, "/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/BioIntegrity in Modified Channels Statewide - SMC/Flow Ecology in Modified Channels/DataPackets/SanDiego_flowEcology_modChannels_stress_standard_thresholds.csv")

## ventura

infoTableVen <- infoTable %>%
  filter(county == "Ventura")

write.csv(infoTableVen, "/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/BioIntegrity in Modified Channels Statewide - SMC/Flow Ecology in Modified Channels/DataPackets/Ventura_flowEcology_modChannels_stress_standard_thresholds.csv")

## LA

infoTableLA <- infoTable %>%
  filter(county == "Los Angeles")

write.csv(infoTableLA, "/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/BioIntegrity in Modified Channels Statewide - SMC/Flow Ecology in Modified Channels/DataPackets/Los_Angeles_flowEcology_modChannels_stress_standard_thresholds.csv")

## San Bernardino

infoTableSB <- infoTable %>%
  filter(county == "San Bernardino")

write.csv(infoTableSB, "/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/BioIntegrity in Modified Channels Statewide - SMC/Flow Ecology in Modified Channels/DataPackets/San_Bernardino_flowEcology_modChannels_stress_standard_thresholds.csv")

## Ornage
infoTableOC <- infoTable %>%
  filter(county == "Orange")

write.csv(infoTableOC, "/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/BioIntegrity in Modified Channels Statewide - SMC/Flow Ecology in Modified Channels/DataPackets/Orange_flowEcology_modChannels_stress_standard_thresholds.csv")

## Riverside

infoTableRS <- infoTable %>%
  filter(county == "Riverside")

write.csv(infoTableRS, "/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/BioIntegrity in Modified Channels Statewide - SMC/Flow Ecology in Modified Channels/DataPackets/Riverside_flowEcology_modChannels_stress_standard_thresholds.csv")


# format number of strikes ---------------------------------------------------
unique(strikes$Result)

## strikes have to be separate due to the formatting of the data
strikes1 <- strikesx %>%
  filter(Index == "csci") %>% ## remove asci
  dplyr::select(masterid, Threshold, Result, NumStrikes) %>% ## remove redundant columns
  mutate(NumStrikes = as.numeric(NumStrikes)) %>%
  drop_na(Threshold) %>%
  pivot_wider(names_from = Result, values_from = NumStrikes) %>%
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>%
  mutate("Total Stressed FFM" = `Unhealthy Biology, Likely Stressed`+ `Unhealthy Biology,  Very Likely Stressed`) %>%
  rename(Modified.Class = Threshold)

## get all site info from above
infoTableSub <- infoTable %>%
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
write.csv(strikes2, "final_data/07_sum_strikes_all_sites_unhealthy_bio_standard_threshold.csv")  

## save to teams
write.csv(strikes2, "/Users/katieirving/Library/CloudStorage/OneDrive-SharedLibraries-SCCWRP/BioIntegrity in Modified Channels Statewide - SMC/Flow Ecology in Modified Channels/DataPackets/Number_FFM_per_Site_standard_thresholds.csv")

