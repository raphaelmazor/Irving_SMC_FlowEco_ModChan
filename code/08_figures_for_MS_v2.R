## Figures for MS New analysis

library(tidyverse)
library(sf)
library(tidylog)

## labels for FFMs
labels <- read.csv("input_data/ffm_names.csv")
labels <- labels[1:24, ]
labels <- labels %>% rename(hydro.endpoints = Flow.Metric.Code)
labels[25, 1] <- "Peak Flow Magnitude (Q99, cfs)"
labels[25, 2] <- "Q99"
labels[25, 3] <- "Peak Flow Magnitude"
labels

## raw data
load(file = "final_data/01_bugs_algae_flow_joined_by_masterid.RData") 

## directory for figures
out.dir <- "MS_Figures/"

getwd()

# All sites map -----------------------------------------------------------

## ca state
ca_sf <- st_read("ignore/SpatialData/California/Ca_State_poly.shp")
## upload ca counties

counties_sf <- st_read("ignore/SpatialData/Counties/Counties.shp")

## get only socal
counties_socal_sf<-counties_sf %>%
  filter(NAME_PCASE %in% c("Ventura","Los Angeles", "Orange","San Bernardino", "Riverside", "San Diego"))

## upload watersheds

sheds_sf<- st_read("ignore/SpatialData/SMCSheds2009/SMCSheds2009.shp")

## set bounding box 


## create beige palet

beige_pal<-c("#f2dbb7","#eed9c4","#fff0db","#e4d5b7","#d9b99b",
             "#d9c2ba","#9c8481","#e2cbb0","#a69279","#f2dbb7",
             "#f6e6bf","#a69279","#f2dbb7","#9c8481","#e4d5b7")


## data
load(file = "final_data/01_bugs_algae_flow_joined_by_masterid.RData") 

## get only columns needed
names(AllData)

impsx <- AllData %>%
  dplyr::select(masterid, longitude, latitude, comid, channel_engineering_class) %>%
  # filter(!Threshold %in% c("NATMed", "NATLow", "NATHigh")) %>%
  drop_na() %>%
  distinct(masterid, .keep_all = T )

length(unique(impsx$masterid)) ## 396
length(unique(impsx$comid)) ## 271

impsx %>% group_by(channel_engineering_class) %>% tally()


## make spatial and format channel type
imps_sf <- impsx %>%
  st_as_sf(coords=c( "longitude", "latitude"), crs=4326, remove=F) %>%
  mutate(ModifiedClass = factor(channel_engineering_class, 
                                levels = c("NAT", "SB0", "SB1", "SB2", "HB"), 
                                labels = c("Natural", "Soft Bottom (0)" ,"Soft Bottom (1)", "Soft Bottom (2)", "Hard Bottom"))) 

## map
m1 <- ggplot()+
  geom_sf(data=ca_sf)+ ## california
  geom_sf(data=sheds_sf, aes(fill=SMC_Name), color=NA)+ ## watersheds
  scale_fill_manual(values=beige_pal, guide= "none")+ ## beige colour for watersheds
  geom_sf(data=counties_socal_sf, fill=NA)+ ## counties
  geom_sf(data=imps_sf,aes(colour = ModifiedClass), size = 2) + ## results
  scale_colour_manual(values=c("chartreuse4", "dodgerblue2", "darkblue", "mediumpurple2", "firebrick3"))+ ## colour of points
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

file.name1 <- paste0(out.dir, "08_all_sites_map_v2.jpg")
ggsave(m1, filename=file.name1, dpi=600, height=7, width=10)

# Ranges of Delta H -------------------------------------------------------

## upload data
load(file = "final_data/01_bugs_algae_flow_joined_by_masterid.RData")
head(AllData)

AllDatax <- AllData %>%
  drop_na(channel_engineering_class) %>% ## remove sites with no channel class
  mutate(ChannelType = factor(channel_engineering_class, levels = c("NAT", "SB1", "SB0", "SB2", "HB"),
                              labels = c("NAT & SB1", "NAT & SB1", "SB0 & SB2","SB0 & SB2", "HB"))) ## make channel type a factor

unique(AllDatax$channel_engineering_class)

## box plot of Ranges
m=1
m 
## define metrics
mets <- unique(na.omit(AllDatax$Flow.Metric.Name))
mets


# AllDatax <- AllData %>%
#   drop_na(channel_engineering_class) %>% ## remove sites with no channel class
#   mutate(channel_engineering_class = recode_factor(channel_engineering_class, NAT = "Natural", SB0 = "Soft Bottom \n (no hard sides)",
#                                                    SB1 = "Soft Bottom (1 hard side)", SB2 = "Soft Bottom \n (2 hard sides)", 
#                                                    HB = "Hard Bottom")) %>% ## change names
#   filter(!channel_engineering_class == "Soft Bottom (1 hard side)") ## remove 0 hard sides
# ## check names
# unique(AllDatax$channel_engineering_class)
# 
# 
# ## box plot of Ranges
# m=1
# 
# ## define metrics
# mets <- unique(na.omit(AllDatax$Flow.Metric.Name))
# mets

## loop through metrics one by one
for(m in 1:length(mets)) {
  
  ## boxplot
  T1 <- (ggplot(subset(AllDatax, Flow.Metric.Name == mets[m]),  aes(x=ChannelType, y=deltah_final, fill = ChannelType)) +
           geom_boxplot() +
           scale_fill_manual(values=c("chartreuse4", "dodgerblue2", "firebrick3"))+
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

## all on one plot

## boxplot
T2 <- (ggplot(AllDatax,  aes(x=ChannelType, y=deltah_final, fill = ChannelType)) +
         geom_boxplot() +
         scale_fill_manual(values=c("chartreuse4", "dodgerblue2", "firebrick3"))+
         scale_x_discrete(name=paste("")) +
         scale_y_continuous(name = paste0("Delta FFM")) +
         theme_classic() +
         geom_hline(yintercept = 0, linetype="dashed", linewidth=0.5, color = "grey50") +
         theme(legend.position="none")) +
  facet_wrap(~Flow.Metric.Name, scales = "free_y") +
  theme(legend.title = element_blank(), 
        # legend.position = "bottom",
        legend.text=element_text(size=25),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 25))


T2

file.name1 <- paste0(out.dir, "08_ALL_boxplot_delta_range.jpg")
ggsave(T2, filename=file.name1, dpi=300, height=20, width=22)


# Bar plots of categories -------------------------------------------------

## data
tallyCats <- read.csv("ignore/12_tally_categories_per_target_for_figures_V2.csv")
head(tallyCats)

## plot all FFM

a1 <- ggplot(tallyCats, aes(fill=Categories, y=PercentSites, x=ChannelType)) + 
  geom_bar(position="stack", stat="identity") +
  # scale_fill_manual(values=c("chartreuse4", "dodgerblue1", "orange","pink", "firebrick3"))+ ## colour of points
  scale_fill_manual(values=hcl.colors(n=5, palette = "Zissou 1"))+ ## colour of points
  facet_wrap(~Flow.Metric.Name) +
  # scale_fill_manual(values=catPal)+
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Sites (%)")

a1


file.name1 <- paste0(out.dir, "08_mixed_mod_percent_sites_per_cat_bar.jpg")
ggsave(a1, filename=file.name1, dpi=600, height=7, width=10)

## loop for FFM separately

## define ffms
ffms <- unique(tallyCats$hydro.endpoints)
s=1
## loop
for(s in 1:length(ffms)) {
  
  ## filter to ffm
  ffmCats <- tallyCats %>%
    filter(hydro.endpoints == ffms[s])
  
  ## define formal name
  ffmNam <- ffmCats$Flow.Metric.Name[s]
  
  ## plot
  a2 <- ggplot(ffmCats, aes(fill=Categories, y=PercentSites, x=ChannelType)) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(values=hcl.colors(n=5, palette = "Zissou 1"))+ ## colour of points
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.text=element_text(size=15),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          title = element_text(size = 15)) +
    labs(title = paste(ffmNam)) +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Sites (%)") 
  
  a2
  ## save
  file.name1 <- paste0(out.dir, "08_", ffms[s], "_mixed_mod_percent_sites_per_cat.jpg")
  ggsave(a2, filename=file.name1, dpi=600, height=7, width=10)
  
} ## end loop


# Maps of categories ------------------------------------------------------

## data 
categories2 <- read.csv("ignore/12_categories_rel_changes_all_values.csv")

head(categories2)

## get coordinates
head(AllData) ## full df

coords <- AllData %>%
  select(masterid, longitude, latitude) %>%
  distinct() ## select only site info 

head(coords)

### add to categories df

categoriesSP <- categories2 %>%
  inner_join(coords, by = "masterid") %>%
  st_as_sf(coords=c( "longitude", "latitude"), crs=4326, remove=F) %>% ## make spatial
  filter(!Categories == "Type X") %>% ## remove healthy sites
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
    scale_colour_manual(values=hcl.colors(n=5, palette = "Zissou 1"))+ ## colour of points
    coord_sf(xlim=c(-119.41,-116.4),
             ylim=c(32.5, 34.8),
             crs=4326) +
    theme_bw()+
    labs(title = paste(metName)) +
    theme(legend.title = element_blank(), 
          legend.position = "bottom",
          legend.text=element_text(size=18),
          plot.title = element_text(size = 15),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10, hjust = 0.05),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)) +
    # facet_grid(rows = vars(BioResult), cols = vars(FlowResult)) +
    facet_wrap(~ChannelType, nrow= 2) +
    guides(col= guide_legend(title= ""))
  
  
  file.name1 <- paste0(out.dir, "08_", metrics[m], "_mixed_mod_maps_per_result_New.jpg")
  ggsave(m1, filename=file.name1, dpi=500, height=10, width=10)
  
}


# Gams plots --------------------------------------------------------------

load(file= "ignore/09_mixed_effects_model_predictions.RData")
head(DF)

## smoothing functions per FFM
sm6csci <- c("DS_Mag_50", "FA_Mag", "Peak_10", "Peak_2", "Peak_5", "Wet_BFL_Mag_50", "Wet_BFL_Mag_10")

sm3csci <- c( "SP_Mag",  "Q99")

## format data
DF <- DF %>% ## format names
  # select(-X) %>%
  mutate(hydro.endpoints = case_when(Variable == "d_ds_mag_50" ~ "DS_Mag_50",            
                                     Variable == "d_fa_mag" ~ "FA_Mag",
                                     Variable == "d_peak_10" ~ "Peak_10",
                                     Variable == "d_peak_2" ~ "Peak_2",
                                     Variable == "d_peak_5" ~ "Peak_5",
                                     Variable == "d_sp_mag" ~ "SP_Mag",
                                     Variable == "d_wet_bfl_mag_10" ~ "Wet_BFL_Mag_10",
                                     Variable == "d_wet_bfl_mag_50" ~ "Wet_BFL_Mag_50", 
                                     Variable == "delta_q99" ~ "Q99")) %>%
  mutate(ChosenSmth = case_when((Metric == "csci" & hydro.endpoints %in% sm6csci) ~ 6,
                                (Metric == "csci" & hydro.endpoints %in% sm3csci) ~ 3)) %>% ## take smooths for FFMs
  mutate(ModelstoUse = ifelse(Smooths == ChosenSmth, "Yes", "No")) %>%
  filter(ModelstoUse == "Yes") %>% ## take only smooths needed
  left_join(labels, by ="hydro.endpoints")

head(DF)

## define FFMs
mets <- unique(DF$hydro.endpoints)
mets
m=1

## start loop
for(m in 1:length(mets)) {
  
  DFX <- DF %>%
    filter(hydro.endpoints == mets[m])
  
  ## format names
  DFX1 <- DFX %>%
    pivot_longer(Intercept:SlopeAndIntercept, names_to = "ModelType", values_to = "predictedVals") %>%
    mutate(ModelType = factor(ModelType, levels = c("Intercept", "Slope", "SlopeAndIntercept"), 
                              labels = c("Intercept", "Slope", "Slope and Intercept"))) %>%
    mutate(channel_engineering_class = factor(channel_engineering_class, 
                                              levels = c("NAT", "SB0", "SB1", "SB2", "HB"), 
                                              labels = c("NAT", "SB0", "SB1", "SB2", "HB"))) 
  
  ## plot
  T2 <- ggplot() +
    geom_smooth(data = DFX1, aes(y=predictedVals, x=deltah_final, col = channel_engineering_class), linewidth = 1)+
    geom_hline(yintercept = 0.79,  linetype="dashed", linewidth=0.5, color = "grey50") +
    geom_hline(yintercept = 0.67,  linetype="dotted", linewidth=0.5, color = "grey50") +
    geom_vline(xintercept = 0) +
    scale_colour_manual(values=c("chartreuse4", "dodgerblue2", "darkblue", "mediumpurple2", "firebrick3"))+ ## colour of points
    scale_x_continuous(name="Delta (CFS)") +
    facet_wrap(~ModelType) +
    scale_y_continuous(name = paste0("CSCI Score"), limits = c(0, 1)) +
    theme_classic() +
    theme(legend.title = element_blank(), 
          # legend.position = "bottom",
          legend.text=element_text(size=12),
          axis.text.x = element_text(size = 12, angle = 20, vjust = 0.5,hjust = 0.3),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12))
  
  
  T2
  # theme(legend.position = "none"))
  

  file.name1 <- paste0(out.dir, "08_", mets[m],  "_mixed_effects_GAM_predicted.jpg")
  ggsave(T2, filename=file.name1, dpi=300, height=5, width=7.5)
  
}

## peak flow ranges

## above 0
DFU <- DF %>%
  filter(hydro.endpoints == mets[5],
         deltah_final > 0) 

## below 0 
DFB <- DF %>%
  filter(hydro.endpoints == mets[5],
         deltah_final < 0) 

head(DFB)

range(DFB$deltah_final) ## [1] -25260.564408     -1.718226
round(range(DFU$deltah_final)) # 0 8059


