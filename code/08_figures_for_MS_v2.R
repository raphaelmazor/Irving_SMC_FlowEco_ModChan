## Figures for MS New analysis

library(tidyverse)
library(sf)
library(tidylog)

# FFM boxplots ------------------------------------------------------------

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
st_bbox(imps_sf)

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
