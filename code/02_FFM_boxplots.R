## boxplots of FFM

library(tidyverse)
library(sf)
# library(mapview)
library(tidylog)

## directory for figures
out.dir <- "final_figures/"

getwd()

# Ranges of Delta H -------------------------------------------------------

## upload data
load(file = "final_data/01_bugs_algae_flow_joined_by_masterid.RData")
head(AllData)

AllDatax <- AllData %>%
  drop_na(channel_engineering_class) %>% ## remove sites with no channel class
  mutate(channel_engineering_class = recode_factor(channel_engineering_class, NAT = "Natural", SB0 = "Soft Bottom (no hard sides)",
                                                   SB1 = "Soft Bottom (1 hard side)", SB2 = "Soft Bottom (2 hard sides)", 
                                                   HB = "Hard Bottom")) ## change names
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
  T1 <- (ggplot(subset(AllDatax, Flow.Metric.Name == mets[m]),  aes(x=channel_engineering_class, y=deltah_final)) +
           geom_boxplot() +
           scale_x_discrete(name=paste("")) +
           scale_y_continuous(name = paste0("Delta: ", mets[m])))
  
  
  T1
  
  file.name1 <- paste0(out.dir, "02_", mets[m], "_boxplot_delta_range.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
}

## jitter plot: height to be 0, unchanges, changed below a therehold, changed above a threshold. delta h limits below critical 


# Number of sites per ffm -------------------------------------------------

perFFM <- AllData %>%
  drop_na(deltah_final) %>%
  select(masterid, Metric, Flow.Metric.Name) %>%
  filter(Metric == "csci") %>%
  distinct() %>%
  group_by(Metric, Flow.Metric.Name) %>%
  summarise(NSites = length(unique(masterid)))

perFFM 


perEng <- AllData %>%
  drop_na(deltah_final) %>%
  select(masterid, Metric, channel_engineering_class) %>%
  filter(Metric == "csci") %>%
  distinct() %>%
  group_by(Metric, channel_engineering_class) %>%
  summarise(NSites = length(unique(masterid)))

perEng 
