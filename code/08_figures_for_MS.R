## all figures for manuscript 

library(tidyverse)
library(sf)
library(tidylog)


# FFM boxplots ------------------------------------------------------------

## directory for figures
out.dir <- "final_figures/"

getwd()

# Ranges of Delta H -------------------------------------------------------

## upload data
load(file = "final_data/01_bugs_algae_flow_joined_by_masterid.RData")
head(AllData)

AllDatax <- AllData %>%
  drop_na(channel_engineering_class) %>% ## remove sites with no channel class
  mutate(channel_engineering_class = recode_factor(channel_engineering_class, NAT = "Natural", SB0 = "Soft Bottom \n (no hard sides)",
                                                   SB1 = "Soft Bottom (1 hard side)", SB2 = "Soft Bottom \n (2 hard sides)", 
                                                   HB = "Hard Bottom")) %>% ## change names
  filter(!channel_engineering_class == "Soft Bottom (1 hard side)") ## remove 0 hard sides
## check names
unique(AllDatax$channel_engineering_class)


## box plot of Ranges
m=1

## define metrics
mets <- unique(na.omit(AllDatax$Flow.Metric.Name))
mets

## loop through metrics one by one
for(m in 1:length(mets)) {
  
  ## boxplot
  T1 <- (ggplot(subset(AllDatax, Flow.Metric.Name == mets[m]),  aes(x=channel_engineering_class, y=deltah_final, fill = channel_engineering_class)) +
           geom_boxplot() +
           scale_fill_manual(values=c("chartreuse4", "dodgerblue2", "mediumpurple2", "firebrick3"))+
           scale_x_discrete(name=paste("")) +
           scale_y_continuous(name = paste0("Delta: ", mets[m])) +
           theme_classic() +
           geom_hline(yintercept = 0, linetype="dashed", linewidth=0.5, color = "grey50") +
           theme(legend.position="none")) +
            theme(legend.title = element_blank(), 
        # legend.position = "bottom",
            legend.text=element_text(size=15),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            strip.text.x = element_text(size = 15),
            strip.text.y = element_text(size = 15))
  
  
  T1
  
  file.name1 <- paste0(out.dir, "08_", mets[m], "_boxplot_delta_range.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=7, width=9)
}


# Quantile Gams -----------------------------------------------------------

## upload prediction data
load(file = "ignore/04_quantGams_smooths_predictions_updatedSites.RData")

## upload raw data for points
load(file = "final_data/01_bugs_algae_flow_joined_by_masterid.RData")
AllDataLongx <- AllData

### predicted figures
bio <- unique(DF$Metric)
mets <- unique(DF$Variable)

  
  for(m in 1:length(mets)) {
    
    ptsbiox <- AllDataLongx %>%
      filter(flow_metric == mets[m]) %>%
      mutate(channel_engineering_class = as.factor(channel_engineering_class))
    
    ptsbiox
    
    T1 <- ggplot() +
      geom_smooth(data = subset(DF1, Variable == mets[m]), aes(y=predictedVals, x=deltah_final, col = Quant), linewidth = 1)+
      geom_point(data = ptsbiox, aes(x=deltah_final, y = MetricValue,
                                     col = channel_engineering_class)) +
      # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
      geom_hline(yintercept = 0.79,  linetype="dashed", linewidth=0.5, color = "grey50") +
      geom_hline(yintercept = 0.67,  linetype="dashed", linewidth=0.5, color = "red") +
      geom_vline(xintercept = 0) +
      facet_wrap(~Smooths, scales = "free") +
      scale_x_continuous(name=paste(mets[m])) +
      scale_y_continuous(name = paste0("CSCI Score")) +
      theme(legend.title = element_blank(), 
            # legend.position = "bottom",
            legend.text=element_text(size=15),
            axis.text.x = element_text(size = 15, angle = 20, vjust = 0.5,hjust = 0.3),
            axis.text.y = element_text(size = 15),
            axis.title = element_text(size = 15),
            strip.text.x = element_text(size = 15),
            strip.text.y = element_text(size = 15))
      
    
  
    # theme(legend.position = "none"))
    
    T1
    
    file.name1 <- paste0(out.dir, "08_", bio[b], "_", mets[m], "_flow_response_predicted_gam_combined_all_quants_updatedSites.jpg")
    ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
    
  }
  
  



# Figures for smoothing function ------------------------------------------

### predicted figures
bio <- unique(DF$Metric)
mets <- unique(DF$Variable)

names(DF)

DFx <- DF %>%
  filter(Quant == 0.5, Smooths %in% c(3,6,9,12))

## plot curve with points for each quant without points

  
  for(m in 1:length(mets)) {
    
    
    T1 <- ggplot() +
      geom_smooth(data = subset(DFx, Variable == mets[m]), aes(y=predictedVals, x=deltah_final, col = Quant), linewidth = 1)+
      # geom_point(data = ptsbiox, aes(x=deltah_final, y = MetricValue,
      #                                col = channel_engineering_class)) +
      # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
      # geom_hline(yintercept = refT,  linetype="dashed", linewidth=0.5, color = "grey50") +
      # geom_hline(yintercept = hbT,  linetype="dashed", linewidth=0.5, color = "red") +
      geom_vline(xintercept = 0) +
      facet_wrap(~Smooths) +
      scale_x_continuous(name=paste(mets[m])) +
      scale_y_continuous(name = paste0("CSCI Score"))
    
    
    # theme(legend.position = "none"))
    
    T1
    
    file.name1 <- paste0(out.dir, "08_", bio[b], "_", mets[m], "_response_smoothing_functions_updatedSites.jpg")
    ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
    
  }
  
## same plot with more smoothings

DFx <- DF %>%
  filter(Quant == 0.5, Smooths %in% c(3,6,9,12, 15))

## plot curve with points for each quant without points


for(m in 1:length(mets)) {
  
  
  T1 <- ggplot() +
    geom_smooth(data = subset(DFx, Variable == mets[m]), aes(y=predictedVals, x=deltah_final, col = Quant), linewidth = 1)+
    # geom_point(data = ptsbiox, aes(x=deltah_final, y = MetricValue,
    #                                col = channel_engineering_class)) +
    # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
    # geom_hline(yintercept = refT,  linetype="dashed", linewidth=0.5, color = "grey50") +
    # geom_hline(yintercept = hbT,  linetype="dashed", linewidth=0.5, color = "red") +
    geom_vline(xintercept = 0) +
    facet_wrap(~Smooths) +
    scale_x_continuous(name=paste(mets[m])) +
    scale_y_continuous(name = paste0("CSCI Score"))
  
  
  # theme(legend.position = "none"))
  
  T1
  
  file.name1 <- paste0(out.dir, "08_", bio[b], "_", mets[m], "_response_smoothing_functions_ALL_updatedSites.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
  
}

### figures with Quantiles

DFx <- DF %>%
  filter( Smooths %in% c(6))

## plot curve with points for each quant without points


for(m in 1:length(mets)) {
  
  
  T1 <- ggplot() +
    geom_smooth(data = subset(DFx, Variable == mets[m]), aes(y=predictedVals, x=deltah_final, col = Quant), linewidth = 1)+
    # geom_point(data = ptsbiox, aes(x=deltah_final, y = MetricValue,
    #                                col = channel_engineering_class)) +
    # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
    # geom_hline(yintercept = refT,  linetype="dashed", linewidth=0.5, color = "grey50") +
    # geom_hline(yintercept = hbT,  linetype="dashed", linewidth=0.5, color = "red") +
    geom_vline(xintercept = 0) +
    facet_wrap(~Smooths) +
    scale_x_continuous(name=paste(mets[m])) +
    scale_y_continuous(name = paste0("CSCI Score"))
  
  
  # theme(legend.position = "none"))
  
  T1
  
  file.name1 <- paste0(out.dir, "08_", bio[b], "_", mets[m], "_response_quantiles_updatedSites.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
  
}



# Maps per category - Standard -------------------------------------------------------

library(viridisLite)
library(tmaptools)

## within limits doc
imps <- read.csv("ignore/05_impact_ffm_bio.csv") %>%
  select( -X)

unique(imps$Result)
## remove all but nat med

imps <- imps %>%
  filter(Threshold == "NATMed")

## add eng data
## upload
BioEng <- read.csv("ignore/Chan_eng_all_SMC.csv") %>% ## upload data
  mutate(Class2 = ifelse(channel_engineering_class =="NAT", "Natural", "Modified")) #%>% ## add overall modification class
  # mutate(comid = as.character(comid)) 

## join to results and format to match
imps2<- left_join(imps, BioEng, by = c("masterid", "comid", "huc", "county", "smcshed")) %>%
  dplyr::select(-Threshold) %>%
  rename(Threshold = channel_engineering_class)

head(imps2)
head(BioEng)
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
                            labels = c("Natural", "Soft Bottom \n (no hard sides)" , "Soft Bottom \n (2 hard sides)", "Hard Bottom"))) %>%
  mutate(Threshold = recode_factor(channel_engineering_class, NAT = "Natural", SB0 = "Soft Bottom \n (no hard sides)",
                                                   SB1 = "Soft Bottom (1 hard side)", SB2 = "Soft Bottom \n (2 hard sides)", 
                                                   HB = "Hard Bottom")) %>% ## change names
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
                                    "Unhealthy Biology,  Very Likely Stressed"))) %>%
  drop_na(Threshold)


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
i=1
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
      # ggtitle(paste0(IndName,": ", metName)) +
      theme(legend.title = element_text("Channel Type"), 
            legend.position = "bottom",
            legend.text=element_text(size=15),
            plot.title = element_text(size = 15),
            axis.text = element_text(size = 12),
            strip.text.x = element_text(size = 15),
            strip.text.y = element_text(size = 15)) +
      facet_grid(rows = vars(BioResult), cols = vars(FlowResult)) +
      guides(col= guide_legend(title= ""))
    
    m1
 
    file.name1 <- paste0(out.dir, "08_", ind[i], "_", metrics[m], "_maps_per_result_Standard.jpg")
    
    ggsave(m1, filename=file.name1, dpi=500, height=10, width=17)
    
  }
  
}


# Mapped results - Modified -----------------------------------------------

## within limits doc
imps <- read.csv("ignore/05_impact_ffm_bio.csv") %>%
  select( -X)

## for now remove natlow,med,high
removes <- unique(imps$Threshold)[c(2,3,4)]

imps <- imps %>%
  filter(!Threshold %in% removes)

## format names
impsx <- imps %>% 
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
                                    "Unhealthy Biology,  Very Likely Stressed"))) %>%
  drop_na(Threshold)

unique(impsx$FlowResult)

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
i=1
m=2
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
      filter(Index == "csci")
    
    ## get bio index and change to upper case
    IndName <- str_to_upper(imps_sfx$Index[1])
    
    m1 <- ggplot()+
      geom_sf(data=ca_sf)+ ## california
      geom_sf(data=sheds_sf, aes(fill=SMC_Name), color=NA)+ ## watersheds
      scale_fill_manual(values=beige_pal, guide= "none")+ ## beige colour for watersheds
      geom_sf(data=counties_socal_sf, fill=NA)+ ## counties
      geom_sf(data=imps_sfx,aes(colour = Threshold), size = 2) + ## results
      scale_colour_manual(values=c("chartreuse4", "dodgerblue2", "mediumpurple2", "firebrick3"))+ ## colour of points
      coord_sf(xlim=c(-119.41,-116.4), ## axis limits
               ylim=c(32.5, 34.8),
               crs=4326) +
      theme_bw()+
      # ggtitle(paste0(IndName,": ", metName)) +
      theme(legend.title = element_text("Channel Type"), 
            legend.position = "bottom",
            legend.text=element_text(size=15),
            plot.title = element_text(size = 15),
            axis.text = element_text(size = 12),
            strip.text.x = element_text(size = 15),
            strip.text.y = element_text(size = 15)) +
      facet_grid(rows = vars(BioResult), cols = vars(FlowResult)) +
      guides(col= guide_legend(title= ""))
    
    m1
    
    ## create directory 
    file.name1 <- paste0(out.dir, "/08_", ind[i], "_", metrics[m], "_maps_per_result_Modified.jpg")
    ## save
    ggsave(m1, filename=file.name1, dpi=500, height=10, width=15)
    
  }
  
}


# Column chart ------------------------------------------------------------

## modified

## tally of categories
tal <- read.csv("final_data/05_count_impact.csv")
head(tal)

## get percentages

talM <- tal %>%
  group_by(Index, Hydro_endpoint, Flow.Metric.Name, Threshold) %>%
  filter(!Threshold %in% removes) %>%
  mutate(TotalPerCat = sum(n)) %>%
  mutate(PercChans = (n/TotalPerCat)*100) %>%
  mutate(Type = "Modified")

talM
unique(talM$Hydro_endpoint)

### standard

talS <- read.csv("final_data/06_count_impact_standard.csv") %>%
  mutate(Type = "Standard")

names(talS)
names(talM)

## combine

tal <- bind_rows(talM, talS)
head(tal)

sum(is.na(tal$Threshold))

## format 
tallyImpactx <- tal %>%
  dplyr::select(-n, -X) %>%
  filter(Index == "csci") %>%
  mutate(ModifiedClass = factor(Threshold, 
                                levels = c("NAT", "SB0", "SB2", "HB"), 
                                labels = c("Natural", "Soft Bottom (0)" , "Soft Bottom (2)", "Hard Bottom"))) %>%
  # pivot_longer(NoImpact:BothImpact, names_to = "Result", values_to = "PercChans") %>%
  mutate(Result = factor(Result, levels = c("HVUF", "UHVUF", "HUF", "UHUF", "HMF", "UHMF", "HLF", "UHLF"))) %>%
                         # labels = c("Healthy Biology, Very Unlikely Stressed",
                         #            "Unhealthy Biology, Very Unlikely Stressed",
                         #            "Healthy Biology, Unlikely Stressed",
                         #            "Unhealthy Biology, Unlikely Stressed",
                         #            "Healthy Biology, Likely Stressed",
                         #            "Unhealthy Biology, Likely Stressed",
                         #            "Healthy Biology,  Very Likely Stressed",
                         #            "Unhealthy Biology,  Very Likely Stressed"))) %>%
  mutate(Type = factor(Type, levels = c("Standard", "Modified"))) %>%
  drop_na(Threshold, ModifiedClass)

length(unique(tallyImpactx$Result))
length(catPal)

sum(is.na(tallyImpactx$Threshold))
# view(tallyImpactx)

## define colours
catPal <- c("lightcyan1", "mistyrose1" ,"lightblue3", "lightpink3","dodgerblue", "red1", "blue", "darkred") 

## format colours to factor levels
names(catPal) <- levels(tallyImpactx$Result)
# colScale <- scale_fill_manual(name = "Result",values = catPal)

## define metrics
mets <- unique(tallyImpactx$Hydro_endpoint)
m=1

for(m in 1:length(mets)) {
  
  tallyImpactx1 <- tallyImpactx %>%
    filter(Hydro_endpoint == mets[m])
  
  a1 <- ggplot(tallyImpactx1, aes(fill=Result, y=PercChans, x=Type)) + 
    geom_bar(position="stack", stat="identity") +
    facet_wrap(~ModifiedClass) +
    scale_fill_manual(name = "Result",values = catPal)+
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Sites (%)") +
    theme(legend.title = element_blank(), 
          # legend.position = "bottom",
          legend.text=element_text(size=15),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15))


  file.name1 <- paste0(out.dir, "08_", mets[m],  "_percentage_sites_per_cat_bar_ALL.jpg")
  ggsave(a1, filename=file.name1, dpi=600, height=7, width=10)
  
}


# Strikes Figure -----------------------------------------------------------

### standard 
## upload number of ffm strikes 
strikes <- read.csv("final_data/05_Number_ffm_per_result.csv") %>%
  filter(Threshold == "NATMed", Index =="csci")  %>%
  dplyr::select(-X)

## add eng data 

## upload
BioEng <- read.csv("ignore/Chan_eng_all_SMC.csv") %>% 
  mutate(Class2 = ifelse(channel_engineering_class =="NAT", "Natural", "Modified")) %>% ## add overall modification class
  mutate(comid = as.character(comid))


## join to results and format to match
strikes2<- left_join(strikes, BioEng, by = "masterid") %>%
  dplyr::select(-Threshold) %>%
  rename(Threshold = channel_engineering_class) %>%
  mutate(Threshold = factor(Threshold, levels = c("NAT", "SB0", "SB2", "HB"), 
                            labels = c("Natural",  "Soft Bottom \n (no hard sides)" , "Soft Bottom \n (2 hard sides)", "Hard Bottom")), 
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
                                    "Unhealthy Biology,  Very Likely Stressed"))) %>%
  mutate(Type = "Standard")

head(strikes2)

## get strikes for likely or very likely stress
boxData2 <- strikes2 %>%
  filter(Index == "csci", Result %in% c("Unhealthy Biology, Likely Stressed", "Unhealthy Biology,  Very Likely Stressed")) %>% ## remove
  dplyr::select( -NumStrikes, - BioResult, -FlowResult, -huc, -smcshed, - Class2, -county, -comid) %>% ## remove redundant columns
  drop_na(Threshold) %>%
  distinct() %>%
  pivot_wider(names_from = Result, values_from = n) %>% ## make wide
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>% ## NAs to 0
  mutate(AllStress = `Unhealthy Biology, Likely Stressed` + `Unhealthy Biology,  Very Likely Stressed`) ## combine
# view(boxData2)

## modified

## within limits doc
imps <- read.csv("ignore/05_impact_ffm_bio.csv") %>%
  select( -X)

##  remove natlow,med,high
removes <- unique(imps$Threshold)[c(2,3,4)]
removes

## upload number of ffm strikes - some sites missing - FIX!!!!!!!
strikes <- read.csv("final_data/05_Number_ffm_per_result.csv") %>%
  filter(!Threshold %in% removes) %>%
  mutate(Threshold = factor(Threshold, levels = c("NAT", "SB0", "SB2", "HB"), 
                            labels = c("Natural", "Soft Bottom \n (no hard sides)" , "Soft Bottom \n (2 hard sides)", "Hard Bottom")), 
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
                                    "Unhealthy Biology,  Very Likely Stressed"))) %>%
  mutate(Type = "Modified")

## get strikes for likely or very likely stress
boxData <- strikes %>%
  filter(Index == "csci", Result %in% c("Unhealthy Biology, Likely Stressed", "Unhealthy Biology,  Very Likely Stressed")) %>% ## remove
  dplyr::select(-X, -NumStrikes, - BioResult, -FlowResult) %>% ## remove redundant columns
  pivot_wider(names_from = Result, values_from = n) %>% ## make wide
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>% ## NAs to 0
  mutate(AllStress = `Unhealthy Biology, Likely Stressed` + `Unhealthy Biology,  Very Likely Stressed`) ## combine


## join box data together
names(boxData)
names(boxData2)

boxDataAll <- bind_rows(boxData, boxData2) %>%
  mutate(Type = factor(Type, levels = c("Standard", "Modified")))
  

## plot

T1 <- (ggplot(boxDataAll,  aes(x=Threshold, y=AllStress, fill = Threshold)) +
         geom_boxplot() +
         scale_x_discrete(name = "") +
         scale_fill_manual(values=c("chartreuse4", "dodgerblue2", "mediumpurple2", "firebrick3"))+ ## colour of boxes
         scale_y_continuous(name = "Number of FFM (out of 9)", breaks = scales::pretty_breaks(9), limits = c(0,9)) +
         facet_wrap(~Type) +
         theme(legend.position = "none",
         strip.text.x = element_text(size = 15),
         strip.text.y = element_text(size = 15),
         axis.text = element_text(size = 13),
         axis.title = element_text(size = 15)))

T1

file.name1 <- paste0(out.dir, "/08_stikes_boxplot_both.jpg")
ggsave(T1, filename=file.name1, dpi=500, height=10, width=15)


# All sites map -----------------------------------------------------------

## data with sites and other info
imps <- read.csv("ignore/05_impact_ffm_bio.csv") %>%
  select( -X)

## get only columns needed
names(imps)

impsx <- imps %>%
  select(masterid:latitude, comid, Threshold) %>%
  filter(!Threshold %in% c("NATMed", "NATLow", "NATHigh")) %>%
  distinct()

unique(impsx$Threshold)
## make spatial and format channel type
imps_sf <- impsx %>%
  st_as_sf(coords=c( "longitude", "latitude"), crs=4326, remove=F) %>%
  mutate(ModifiedClass = factor(Threshold, 
                                levels = c("NAT", "SB0", "SB2", "HB"), 
                                labels = c("Natural", "Soft Bottom (0)" , "Soft Bottom (2)", "Hard Bottom"))) 

## map
m1 <- ggplot()+
  geom_sf(data=ca_sf)+ ## california
  geom_sf(data=sheds_sf, aes(fill=SMC_Name), color=NA)+ ## watersheds
  scale_fill_manual(values=beige_pal, guide= "none")+ ## beige colour for watersheds
  geom_sf(data=counties_socal_sf, fill=NA)+ ## counties
  geom_sf(data=imps_sf,aes(colour = ModifiedClass), size = 2) + ## results
  scale_colour_manual(values=c("chartreuse4", "dodgerblue2", "mediumpurple2", "firebrick3"))+ ## colour of points
  coord_sf(xlim=c(-119.41,-116.4), ## axis limits
           ylim=c(32.5, 34.8),
           crs=4326) +
  theme_bw()+
  # ggtitle(paste0(IndName,": ", metName)) +
  theme(legend.text=element_text(size=15),
        plot.title = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15)) +
  guides(col= guide_legend(title= ""))
  # facet_grid(rows = vars(BioResult), cols = vars(FlowResult)) ## facet by categories

m1

file.name1 <- paste0(out.dir, "08_all_sites_map.jpg")
ggsave(m1, filename=file.name1, dpi=600, height=7, width=10)


# Conceptual gam figure ---------------------------------------------------

## data 
load(file = "ignore/04_quantGams_smooths_predictions_updatedSites.RData")

### predicted figures
bio <- unique(DF$Metric)
mets <- unique(DF$Variable)

### take the assigned K from script 05, just csci and DS baseflow, only 0.5 & 0.9 

## filter to index and only 
DF1 <- DF %>%
  filter(Metric == bio[b], Quant %in% c( 0.3, 0.5, 0.9), Smooths == 6)


T1 <- ggplot() +
  geom_smooth(data = subset(DF1, Variable == mets[m]), aes(y=predictedVals, x=deltah_final, col = Quant), linewidth = 1)+
  # geom_point(data = ptsbiox, aes(x=deltah_final, y = MetricValue,
  #                                col = channel_engineering_class)) +
  # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
  geom_hline(yintercept = 0.79,  linetype="dashed", linewidth=0.5, color = "grey50") +
  geom_hline(yintercept = 0.67,  linetype="dotted", linewidth=0.5, color = "grey50") +
  geom_vline(xintercept = 0) +
  guides(col=guide_legend(title="Flow Alteration Level")) +
  scale_x_continuous(name= "Dry Season Baseflow (cfs)", breaks=seq(-25,100,25)) +
  # expand_limits(x = -25) +
  scale_y_continuous(name = paste0("CSCI Score"))


# theme(legend.position = "none"))

T1

file.name1 <- paste0(out.dir, "08_DS_baseflow_example_updatedSites.jpg")
ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
