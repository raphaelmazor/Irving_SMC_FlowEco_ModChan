## level 2: site specific analysis,  how many of mod chans have issues with flow metric, keep 4 outcomes 

# library(gam)
library(tidyverse)
library(sf)
library(tidylog)

## directory for figures
out.dir <- "final_figures/"

## workflow
## use 6 outcome categories to count how many standard channels have issues with flow metrics
## which flow metrics - visualise
## which flow metrics also issue with bio - visualise


# Upload Data -------------------------------------------------------------

## within limits doc
imps <- read.csv("ignore/05_impact_ffm_bio.csv") %>%
  select( -X)

unique(imps$Result)
## remove all but nat med

imps <- imps %>%
  filter(Threshold == "NATMed")
  view(imps)
## add eng data
## upload
BioEng <- read.csv("ignore/02_chan_eng.csv") %>% ## upload data
  select(-c(X,channel_engineering_classification_date, channel_engineering_personnel, channel_engineering_comments, Class2)) 

## join to results and format to match
imps2<- left_join(imps, BioEng, by = "masterid") %>%
  dplyr::select(-Threshold) %>%
  rename(Threshold = channel_engineering_class)

## tally using eng classes

talx <- imps2 %>%
  group_by(Index, Hydro_endpoint, Flow.Metric.Name, Threshold, Result) %>%
  distinct() %>%
  tally() 

talx

## get percentages
tal <- talx %>%
    group_by(Index, Hydro_endpoint, Flow.Metric.Name, Threshold) %>%
    # filter(Threshold == "NATMed") %>%
    mutate(TotalPerCat = sum(n)) %>%
    mutate(PercChans = (n/TotalPerCat)*100)

tal

write.csv(tal, "final_data/06_count_impact_standard.csv")


# Clean Tables ------------------------------------------------------------


## change order of cols, rename cols, change names of standard
FinalTablex <- tal %>%
  ungroup() %>%
  select(Index, Flow.Metric.Name,  Threshold, Result,  n, PercChans) %>% # 
  rename(FlowMetric = Flow.Metric.Name, ModifiedClass = Threshold, Impact = Result, PercentageOfSites = PercChans, NumberOfSites = n,) %>% #  
  mutate(ModifiedClass = case_when((ModifiedClass == "HB") ~ "Hard Bottom",
                                   (ModifiedClass == "NAT") ~ "Natural",
                                   (ModifiedClass == "NATHigh") ~ "Overall 1st",
                                   (ModifiedClass == "NATLow") ~ "Overall 30th",
                                   (ModifiedClass == "NATMed") ~ "Overall 10th",
                                   (ModifiedClass == "SB0") ~ "Soft Bottom (no hard sides)",
                                   (ModifiedClass == "SB2") ~ "Soft Bottom (two hard sides)"))
FinalTablex    
unique(FinalTablex$Impact)
## get total sites 
sums <- FinalTablex %>%
  group_by(Index, FlowMetric, ModifiedClass) %>%
  summarise(NumberSitesPerClass = sum(NumberOfSites))

sums

write.csv(sums, "final_data/06_number_of_sites_standard_class_per_FFM.csv")

  
## pivot impact results wider
FinalTable <- FinalTablex %>%  
  select(-NumberOfSites) %>%
  pivot_wider(names_from = Impact, values_from = PercentageOfSites) %>%
  drop_na(ModifiedClass) %>%
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>%
  select(Index:ModifiedClass, HUF, UHUF, HMF, UHMF, HLF, UHLF)

## join number of sites

FinalTable <- full_join(FinalTable, sums, by = c("Index", "FlowMetric", "ModifiedClass"))
FinalTable

## save out

write.csv(FinalTable, "final_data/06_percent_impacts_standard_Class.csv")

## table with columns as channel class

FinalTableWide <- FinalTablex %>%
  select(-NumberOfSites) %>%
  drop_na(ModifiedClass) %>%
  pivot_wider(names_from = ModifiedClass, values_from = PercentageOfSites) #%>%
  # mutate(across(everything(), .fns = ~replace_na(.,0))) 
  # select(-)
  str(FinalTableWide)

  FinalTableWide
  
  write.csv(FinalTableWide, "final_data/06_percent_impacts_standard_Class_wide.csv")
  


# Map results -------------------------------------------------------------
  
  library(viridisLite)
  library(tmaptools)
  
  ## upload ca boundary
  
  ca_sf <- st_read("ignore/SpatialData/California/Ca_State_poly.shp")
  ## upload ca counties
  
  counties_sf <- st_read("ignore/SpatialData/Counties/Counties.shp")
  
  ## get only socal
  counties_socal_sf<-counties_sf %>%
    filter(NAME_PCASE %in% c("Ventura","Los Angeles", "Orange","San Bernardino", "Riverside", "San Diego"))
  
  ## upload watersheds
  
  sheds_sf<- st_read("ignore/SpatialData/SMCSheds2009/SMCSheds2009.shp")
  
## format names
impsx <- imps2 %>% 
  dplyr::select(Index, Hydro_endpoint, Threshold, BioThresh,  masterid, comid, Flow.Metric.Name, Flow.Component, Result, longitude, latitude)  %>%
  mutate(Threshold = factor(Threshold, levels = c("NAT", "SB0", "SB2", "HB"), 
                            labels = c("Natural", "Soft Bottom (0)" , "Soft Bottom (2)", "Hard Bottom"))) %>%
  mutate(BioResult = case_when(Result %in% c("HUF", "HMF","HLF", "HVUF") ~ "Healthy Biology",
                               Result %in% c("UHUF", "UHMF", "UHLF", "UHVUF") ~ "Unhealthy Biology")) %>%
  mutate(FlowResult = case_when(Result %in% c("HUF", "UHUF") ~ "Unlikely Stressed",
                                Result %in% c("HVUF", "UHVUF") ~ "Very Unlikely Stressed",
                                Result %in% c("HMF", "UHMF") ~ "Likely Stressed",
                                Result %in% c("HLF", "UHLF") ~ "Very Likely Stressed")) %>%
  mutate(FlowResult = factor(FlowResult, levels = c("Very Unlikely Stressed", "Unlikely Stressed", "Likely Stressed", "Very Likely Stressed"))) %>%
  mutate(Result = factor(Result, levels = c("HVUF", "UHVUF", "HUF", "UHUF", "HMF", "UHMF", "HLF", "UHLF"),
                         labels = c("Healthy Biology, Very Unlikely Stressed",
                                    "Unhealthy Biology, Very Unlikely Stressed",
                                    "Healthy Biology, Unlikely Stressed",
                                    "Unhealthy Biology, Unlikely Stressed",
                                    "Healthy Biology, Likely Stressed",
                                    "Unhealthy Biology, Likely Stressed",
                                    "Healthy Biology,  Very Likely Stressed",
                                    "Unhealthy Biology,  Very Likely Stressed")))  


## make spatial
imps_sf <- impsx %>%
  st_as_sf(coords=c( "longitude", "latitude"), crs=4326, remove=F) 

## set bounding box 
st_bbox(imps_sf)

## create beige palet

beige_pal<-c("#f2dbb7","#eed9c4","#fff0db","#e4d5b7","#d9b99b",
             "#d9c2ba","#9c8481","#e2cbb0","#a69279","#f2dbb7",
             "#f6e6bf","#a69279","#f2dbb7","#9c8481","#e4d5b7")


# Plot per flow metric and result

metrics <- unique(imps_sf$Hydro_endpoint)
metrics

## define index

ind <- unique(imps_sf$Index)
ind

## loop over metrics and index
i=2
m=8
imps_sf
  
  for(m in 1:length(metrics)) {
    
    ## extract metric
    imps_sf1 <- imps_sf %>%
      filter(Hydro_endpoint == metrics[m])
    
    ## get metric fancy name
    metName <- imps_sf1$Flow.Metric.Name[1]

    for(i in 1:length(ind)) {
      
      ## extract - index
      imps_sfx <- imps_sf1 %>%
        filter(Index == ind[i])
      
      # imps_sfx
      
      ## get bio index and change to upper case
      IndName <- str_to_upper(imps_sfx$Index[1])
      
    m1 <- ggplot()+
      geom_sf(data=ca_sf)+ ## california
      geom_sf(data=sheds_sf, aes(fill=SMC_Name), color=NA)+ ## watersheds
      scale_fill_manual(values=beige_pal, guide= "none")+ ## beige colour for watersheds
      geom_sf(data=counties_socal_sf, fill=NA)+ ## counties
      geom_sf(data=imps_sfx,aes(colour = Threshold), size = 2) + ## results
      scale_colour_manual(values=c("chartreuse4", "dodgerblue2", "mediumpurple2", "firebrick3"))+
      coord_sf(xlim=c(-119.41,-116.4),
               ylim=c(32.5, 34.8),
               crs=4326) +
      theme_bw()+
      theme(axis.title = element_blank(),
            panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank()) +
      # guides(size = "none") +
      ggtitle(paste0(IndName,": ", metName)) +
      facet_grid(rows = vars(BioResult), cols = vars(FlowResult))
    
    m1

    head(imps_sfx)
    file.name1 <- paste0(out.dir, "06_", ind[i], "_", metrics[m], "_standard_map_6_cats_per_bioHealth_likelihood_flow_impact_per_result.jpg")

    ggsave(m1, filename=file.name1, dpi=500, height=10, width=15)
    
  }
  
}

# Plot per flow metric and channel type

## define impact acategory
type <- unique(imps_sf$Threshold)

for(m in 1:length(metrics)) {
  
  ## extract metric
  imps_sf1 <- imps_sf %>%
    filter(Hydro_endpoint == metrics[m])
  
  ## get metric fancy name
  metName <- imps_sf1$Flow.Metric.Name[1]
  
  for(i in 1:length(ind)) {
    
    ## extract - index
    imps_sfx <- imps_sf1 %>%
      filter(Index == ind[i])
    
    ## get bio index and change to upper case
    IndName <- str_to_upper(imps_sfx$Index[1])
    
    # map
    m1 <- ggplot()+
      geom_sf(data=ca_sf)+ ## california
      geom_sf(data=sheds_sf, aes(fill=SMC_Name), color=NA)+ ## watersheds
      scale_fill_manual(values=beige_pal, guide= "none")+ ## beige colour for watersheds
      geom_sf(data=counties_socal_sf, fill=NA)+ ## counties
      geom_sf(data=imps_sfx,aes(colour = Result), size = 2) + ## results
      coord_sf(xlim=c(-119.41,-116.4),
               ylim=c(32.5, 34.8),
               crs=4326) +
      theme_bw()+
      theme(axis.title = element_blank(),
            panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank()) +
      # guides(size = "none") +
      ggtitle(paste0(IndName,": ", metName)) +
      facet_wrap(~Threshold)
    
    m1
    
    file.name1 <- paste0(out.dir, "06_", ind[i], "_", metrics[m], "_standard__map_6_cats_per_bioHealth_likelihood_flow_impact_per_channel_type.jpg")

    ggsave(m1, filename=file.name1, dpi=500, height=10, width=15)
    
  }
  
}



# Plot number of strikes --------------------------------------------------

## upload number of ffm strikes 
strikes <- read.csv("final_data/05_Number_ffm_per_result.csv") %>%
  filter(Threshold == "NATMed", Index =="csci")  %>%
  dplyr::select(-X)

## add eng data 

## upload
BioEng <- read.csv("ignore/02_chan_eng.csv") %>% ## upload data
  select(-c(X,channel_engineering_classification_date, channel_engineering_personnel, channel_engineering_comments, Class2)) 

## join to results and format to match
strikes2<- left_join(strikes, BioEng, by = "masterid") %>%
  dplyr::select(-Threshold) %>%
  rename(Threshold = channel_engineering_class) %>%
  mutate(Threshold = factor(Threshold, levels = c("NAT", "SB0", "SB2", "HB"), 
                            labels = c("Natural", "Soft Bottom (0)" , "Soft Bottom (2)", "Hard Bottom")), 
         NumStrikes = as.factor(n)) %>%
  mutate(BioResult = case_when(Result %in% c("HUF", "HMF","HLF") ~ "Healthy Biology",
                               Result %in% c("UHUF", "UHMF", "UHLF", "UHVUF") ~ "Unhealthy Biology")) %>%
  mutate(FlowResult = case_when(Result %in% c("HUF", "UHUF") ~ "Unlikely Stressed",
                                Result %in% c("HVUF", "UHVUF") ~ "Very Unlikely Stressed",
                                Result %in% c("HMF", "UHMF") ~ "Likely Stressed",
                                Result %in% c("HLF", "UHLF") ~ "Very Likely Stressed")) %>%
  mutate(FlowResult = factor(FlowResult, levels = c("Very Unlikely Stressed", "Unlikely Stressed", "Likely Stressed", "Very Likely Stressed"))) %>%
  mutate(Result = factor(Result, levels = c("HVUF", "UHVUF", "HUF", "UHUF", "HMF", "UHMF", "HLF", "UHLF"),
                         labels = c("Healthy Biology, Very Unlikely Stressed",
                                    "Unhealthy Biology, Very Unlikely Stressed",
                                    "Healthy Biology, Unlikely Stressed",
                                    "Unhealthy Biology, Unlikely Stressed",
                                    "Healthy Biology, Likely Stressed",
                                    "Unhealthy Biology, Likely Stressed",
                                    "Healthy Biology,  Very Likely Stressed",
                                    "Unhealthy Biology,  Very Likely Stressed")))

head(strikes2)

## get strikes for likely or very likely stress
boxData2 <- strikes2 %>%
  filter(Index == "csci", Result %in% c("Unhealthy Biology, Likely Stressed", "Unhealthy Biology,  Very Likely Stressed")) %>% ## remove
  dplyr::select( -NumStrikes, - BioResult, -FlowResult) %>% ## remove redundant columns
  drop_na(Threshold) %>%
  pivot_wider(names_from = Result, values_from = n) %>% ## make wide
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>% ## NAs to 0
  mutate(AllStress = `Unhealthy Biology, Likely Stressed` + `Unhealthy Biology,  Very Likely Stressed`) ## combine

head(boxData)

T1 <- (ggplot(boxData2,  aes(x=Threshold, y=AllStress, fill = Threshold)) +
         geom_boxplot() +
         scale_x_discrete(name = "") +
         scale_fill_manual(values=c("chartreuse4", "dodgerblue2", "mediumpurple2", "firebrick3"))+ ## colour of boxes
         scale_y_continuous(name = "Number of FFM", breaks = scales::pretty_breaks(9)))


T1

file.name1 <- paste0(out.dir, "/06_stikes_boxplot_standard.jpg")
ggsave(T1, filename=file.name1, dpi=500, height=10, width=15)



#  Bar plots --------------------------------------------------------------

tal
tallyImpactx <- tal %>%
  dplyr::select(-n) %>%
  filter(Index == "csci") %>%
  mutate(ModifiedClass = factor(Threshold, 
                                levels = c("NAT", "SB0", "SB2", "HB"), 
                                labels = c("Natural", "Soft Bottom (0)" , "Soft Bottom (2)", "Hard Bottom"))) %>%
  drop_na(ModifiedClass) %>%
  # pivot_longer(NoImpact:BothImpact, names_to = "Result", values_to = "PercChans") %>%
    mutate(Result = factor(Result, levels = c("HUF", "UHUF", "HMF", "UHMF", "HLF", "UHLF"), 
                           labels = c("Healthy Biology, Unlikely Stressed","Unhealthy Biology, Unlikely Stressed",
                                    "Healthy Biology, Likely Stressed", "Unhealthy Biology, Likely Stressed",
                                      "Healthy Biology, Very Likely Stressed","Unhealthy Biology, Very Likely Stressed")))


catPal <- c("lightblue3", "lightpink3","dodgerblue", "red1",  "darkred") 


a1 <- ggplot(tallyImpactx, aes(fill=Result, y=PercChans, x=ModifiedClass)) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~Flow.Metric.Name) +
  scale_fill_manual(values=catPal)+
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Sites (%)")

a1

file.name1 <- paste0(out.dir, "06_all_hydro_stacked_perc_standard.jpg")
ggsave(a1, filename=file.name1, dpi=600, height=7, width=10)

unique(tallyImpactx$Flow.Metric.Name)
### get just dry season baseflow

dryTally <- tallyImpactx %>%
  filter(Flow.Metric.Name == "Dry-season median baseflow")

dryTally

a2 <- ggplot(dryTally, aes(fill=Result, y=PercChans, x=ModifiedClass)) + 
  geom_bar(position="stack", stat="identity") +
  # facet_wrap(~Flow.Metric.Name) +
  scale_fill_manual(values=catPal)+
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Sites (%)")

a2

file.name1 <- paste0(out.dir, "06_DS_baseflow_stacked_perc_standard.jpg")
ggsave(a2, filename=file.name1, dpi=600, height=7, width=10)

### get just wet season baseflow

wetTally <- tallyImpactx %>%
  filter(Flow.Metric.Name == "Wet-season median baseflow")

wetTally

a3 <- ggplot(wetTally, aes(fill=Result, y=PercChans, x=ModifiedClass)) + 
  geom_bar(position="stack", stat="identity") +
  # facet_wrap(~Flow.Metric.Name) +
  scale_fill_manual(values=catPal)+
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Sites (%)")

a3

file.name1 <- paste0(out.dir, "06_WS_baseflow_stacked_perc_standard.jpg")
ggsave(a3, filename=file.name1, dpi=600, height=7, width=10)

### get just peaks

peakTally <- tallyImpactx %>%
  filter(Flow.Metric.Name %in% c("10-year flood magnitude", "2-year flood magnitude",
                                 "5-year flood magnitude", "Magnitude of largest annual storm") )

peakTally

a4 <- ggplot(peakTally, aes(fill=Result, y=PercChans, x=ModifiedClass)) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~Flow.Metric.Name) +
  scale_fill_manual(values=catPal)+
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Sites (%)")

a4

file.name1 <- paste0(out.dir, "06_peak_flows_stacked_perc_standard.jpg")
ggsave(a4, filename=file.name1, dpi=600, height=7, width=11)
