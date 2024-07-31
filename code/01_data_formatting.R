library(tidylog)
library(tidyverse)
library(sf)
# library(mapview)

out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/SMC_Modified_Channels/final_figures/"

getwd()


# Flow data ---------------------------------------------------------------

## upload data
dh_data <- read.csv("ignore/2024-07-26_RFpred_output_alldata_chan.engCOMIDs_med_dlt_FFM_test12_test2.csv")
head(dh_data)

## pivot longer
dh_median <- dh_data %>%
  pivot_longer(d_ds_mag_50:delta_q99, names_to = "flow_metric", values_to = "deltah_final") 

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

labels
## count sites with flow data

flowSites <- unique(dh_median$masterid)

length(flowSites) ## 1149

## channel engineering data

BioEng <- read.csv("ignore/02_chan_eng.csv") %>% ## upload data
  select(-c(X,channel_engineering_classification_date, channel_engineering_personnel, channel_engineering_comments)) %>% ## remove redundant columns
  mutate(Class2 = ifelse(channel_engineering_class =="NAT", "Natural", "Modified")) ## add overall modification class

BioEng

## count number of sites per class
tallyFFMClass <- BioEng %>%
  group_by(channel_engineering_class) %>%
  distinct() %>%
  drop_na(channel_engineering_class) %>%
  tally()

tallyFFMClass


# Bio data ----------------------------------------------------------------

##  upload all sites

bioSites <- st_read("ignore/01_bio_sites_all.shp")
head(bioSites)
dim(bioSites)

## bugs data
load(file = "ignore/SMC_csci_cali_Mar2023_v2.RData")

# csciScores <- read.csv("ignore/01_csci_comp_mets_comid_socal.csv")
head(bug_tax_ca)

### algae scores
load(file="ignore/SMC_asci_cali_Mar2023_v3.RData")
head(alg_tax_ca)
# asciScores <- read.csv("ignore/01_asci_comp_mets_comid_socal.csv")

# Join bio and flow -------------------------------------------------------

## how many flow in bug sites
sum(flowSites %in% bug_tax_ca$masterid) ## 982

sum(flowSites %in% alg_tax_ca$masterid) ## 722

## filter bug data using masterid ### remove reps format date
csciScores <- bug_tax_ca %>%
  filter(fieldreplicate == 1 ) %>%
  mutate(Metric = "csci", csci = as.numeric(csci)) %>%
  rename(MetricValue = csci) %>%
  select(masterid, sampleyear, Metric, MetricValue, longitude, latitude, comid)

length(unique(csciScores$masterid)) ## 4413 sites (all ca)


## filter bug data using masterid ### remove reps, format date and data
asciScores <- alg_tax_ca %>%
  filter(replicate == 1 ) %>%
  separate(sampledate, into = c("sampledate", "Time"), sep= " ", remove = F) %>%
  separate(sampledate, into = c("sampleyear", "Month", "Day"), sep= "-", remove = F) %>%
  mutate(sampleyear = as.numeric(sampleyear))  %>%
  mutate(Metric = "asci") %>%
  rename(MetricValue = result) %>%
  mutate(MetricValue = as.numeric(MetricValue)) %>%
  select(masterid, sampleyear, Metric, MetricValue, longitude, latitude,  comid)


str(asciScores)
length(unique(asciScores$masterid)) ## 2306 sites 

## check names to join
names(asciScores)
names(csciScores)

# Join bio sites to flow data ---------------------------------------------

## join asci and csci
scoresSites <- bind_rows(asciScores, csciScores)
sum(is.na(scoresSites$Metric))

## join channel class to bio sites
engsites <- left_join(scoresSites, BioEng, by = "masterid") %>%
  mutate(comid = as.integer(comid))
  
engsites

length(unique(engsites$masterid)) ## 4427
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
AllData1 

## count FFM per channel class

tallyFFM <- AllData1 %>%
  group_by(channel_engineering_class, HasFFM) %>%
  select(-sampleyear, - MetricValue, -Metric) %>% distinct() %>%
  # drop_na(channel_engineering_class) %>%
  tally()

tallyFFM

## join all data

names(dh_median)
names(engsites)

AllData <- right_join(engsites, dh_median, by = c("masterid", "comid", "channel_engineering_class"), relationship = "many-to-many") %>%
          left_join( labels, by = "flow_metric")

AllData


## save out
save(AllData, file = "final_data/01_bugs_algae_flow_joined_by_masterid.RData")


