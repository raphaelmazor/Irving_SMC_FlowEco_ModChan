
library(tidyverse)
library(sf)
library(tidylog)
# library(mapview)

out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/SMC_Modified_Channels/final_figures/"

getwd()


# Flow data ---------------------------------------------------------------

## upload data
dh_data <- read.csv("ignore/2024-07-26_RFpred_output_alldata_chan.engCOMIDs_med_dlt_FFM_test12_test2.csv")
head(dh_data)

## pivot longer
dh_median <- dh_data %>%
  pivot_longer(d_ds_mag_50:delta_q99, names_to = "flow_metric", values_to = "deltah_final") %>%
  mutate(comid = as.character(comid))

dim(dh_median) # 13329
head(dh_median)
str(dh_median)

## full names for FFM labels
labels <- read.csv("input_data/ffm_names.csv")
labels <- labels[1:24, ]
labels <- labels %>% rename(hydro.endpoints = Flow.Metric.Code) ## rename to match
## add Q99
labels[25, 1] <- "Magnitude of largest annual storm"
labels[25, 2] <- "Q99"
labels[25, 3] <- "Peak flow"

## change all names to match dh_data
labels <- labels %>%
  mutate(flow_metric = case_when(hydro.endpoints == "DS_Mag_50" ~ "d_ds_mag_50",
                                 hydro.endpoints == "FA_Mag" ~ "d_fa_mag",
                                 hydro.endpoints == "Peak_10" ~ "d_peak_10",
                                 hydro.endpoints == "Peak_2" ~ "d_peak_2",
                                 hydro.endpoints == "Peak_5" ~ "d_peak_5",
                                 hydro.endpoints == "SP_Mag" ~ "d_sp_mag",
                                 hydro.endpoints == "Wet_BFL_Mag_10" ~ "d_wet_bfl_mag_10",
                                 hydro.endpoints == "Wet_BFL_Mag_50" ~ "d_wet_bfl_mag_50",
                                 hydro.endpoints == "Q99" ~ "delta_q99")) %>%
  drop_na(flow_metric)

## count sites with flow data

flowSites <- unique(dh_median$masterid)

length(flowSites) ## 1149

## channel engineering data

BioEng <- read.csv("ignore/Chan_eng_all_SMC.csv") %>% ## upload data
  # select(-c(X,channel_engineering_classification_date, channel_engineering_personnel, channel_engineering_comments)) %>% ## remove redundant columns
  mutate(Class2 = ifelse(channel_engineering_class =="NAT", "Natural", "Modified")) %>% ## add overall modification class
  mutate(comid = as.character(comid))

head(BioEng)

## count number of sites per class
tallyFFMClass <- BioEng %>%
  group_by(channel_engineering_class) %>%
  distinct() %>%
  drop_na(channel_engineering_class) %>%
  tally()

tallyFFMClass


# Bio data ----------------------------------------------------------------

##  upload all sites

## bugs data
load(file = "ignore/CSCI_CA_Aug2024.RData")
dim(csci)
bug_tax_ca <- csci

head(bug_tax_ca)

# Join bio and flow -------------------------------------------------------

## how many flow in bug sites
sum(flowSites %in% bug_tax_ca$masterid) ## 1073


## filter bug data using masterid ### remove reps format date
csciScores <- bug_tax_ca %>%
  filter(fieldreplicate == 1 ) %>%
  mutate(Metric = "csci", csci = as.numeric(csci)) %>%
  rename(MetricValue = csci) %>%
  select(masterid, sampleyear, Metric, MetricValue, longitude, latitude, comid)

length(unique(csciScores$masterid)) ## 5729 sites (all ca)
str(csciScores)

# Join bio sites to flow data ---------------------------------------------

## join channel class to bio sites
engsites <- left_join(csciScores, BioEng, by = c("masterid", "comid"), relationship = "many-to-many") #%>%
  # mutate(comid = as.integer(comid))
  
engsites

length(unique(engsites$masterid)) ## 5729
sum(is.na(engsites$Metric))

## count asci csci by class 

tallyclass <- engsites %>%
  group_by(Metric, channel_engineering_class) %>%
  select(-sampleyear, - MetricValue) %>% distinct() %>%
  tally()

tallyclass

## define sites with FFM
flowsites <- dh_median %>%
  select(masterid) %>%
  mutate(HasFFM = "Yes")

flowsites

## join to channel class
AllData1 <- full_join(engsites, flowsites, by = c("masterid"), relationship = "many-to-many") %>%
  mutate(HasFFM = replace_na(HasFFM, "No")) ## replace NAs with No (no ffm)
head(AllData1)

## count FFM per channel class

tallyFFM <- AllData1 %>%
  group_by(channel_engineering_class, HasFFM) %>%
  select(-sampleyear, - MetricValue, -Metric) %>% distinct() %>%
  # drop_na(channel_engineering_class) %>%
  tally()

tallyFFM

## join all data

names(dh_median) %in% names(engsites)
dim(dh_median)

names(engsites)

AllData <- right_join(engsites, dh_median, by = c("masterid", "comid", "channel_engineering_class", "huc", "county", "smcshed"), relationship = "many-to-many") %>%
          left_join( labels, by = "flow_metric") %>%
  drop_na(deltah_final)

AllData


## save out
save(AllData, file = "final_data/01_bugs_algae_flow_joined_by_masterid.RData")


