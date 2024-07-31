## peak metrics show unexpected relationship with bio, i.e., positive
## remove from analysis, can cite RB9 tech report

library(tidyverse)
library(tidylog)

## directory for figures
out.dir <- "final_figures/"

getwd()

# Upload & Format Data -------------------------------------------------------------

load(file = "final_data/01_bugs_algae_flow_joined_by_masterid.RData")
head(AllData)

## how many bio sites 
length(unique(AllData$masterid)) ## 479


## take only acsi h and csci
AllDataLong2 <- AllData %>%
  filter(Metric %in% c("csci", "asci"))

## check number of sites
length(unique(AllDataLong2$masterid)) ## 473

# Remove peak flow augmentation -------------------------------------------

head(AllDataLong2)

## define peak mets
peakmets <- unique(na.omit(AllDataLong2$flow_metric))[c(3,4,5,9)]
peakmets

## assign which to be removed
AllDataLongA <- AllDataLong2 %>%
  # dplyr::select(-X) %>%
  distinct() %>%
  drop_na(deltah_final) %>%
  mutate(PeakAug = ifelse(flow_metric %in% peakmets & deltah_final > 0, "Yes", "No"))

## how many bio sites 
length(unique(AllDataLongA$masterid)) ## 187, 364

## count number of sites with augmnented peaks
talTest <- AllDataLongA %>%
  filter(flow_metric %in% peakmets) %>%
  group_by(Metric, flow_metric) %>%
  summarise(NoSitesAugmented = sum(na.omit(PeakAug == "Yes")))

talTest

## plot CSCI data points

## subset to only csci
cscitest <- AllDataLongA %>%
  filter(Metric == "csci")

names(cscitest)

P1 <- ggplot() +
  geom_point(data = cscitest, aes(x=deltah_final, y = MetricValue,
                                  col = channel_engineering_class)) + ## plot points
  geom_vline(xintercept = 0) + ## add line at 0
  facet_wrap(~Flow.Metric.Name, scales = "free")## facet by flow metric

P1

## save out
file.name1 <- paste0(out.dir, "03_raw_data_plot_CSCI.jpg")
ggsave(P1, filename=file.name1, dpi=300, height=5, width=7.5)


## plot ASCI

## subset to asci
ascitest <- AllDataLongA %>%
  filter(Metric == "asci")


P2 <- ggplot() +
  geom_point(data = ascitest, aes(x=deltah_final, y = MetricValue,
                                  col = channel_engineering_class)) + ## plot points
  geom_vline(xintercept = 0) + ## line at 0
  facet_wrap(~Flow.Metric.Name, scales = "free")## facet by flow metric


P2

file.name1 <- paste0(out.dir, "03_raw_data_plot_ASCI.jpg")
ggsave(P2, filename=file.name1, dpi=300, height=5, width=7.5)


## remove augmented peak flows and save
AllDataLongx <- AllDataLongA %>%
  filter(!PeakAug == "Yes")

## check number of sites
unique(AllDataLongx$masterid) ## 187
unique(AllDataLongA$masterid) ## 187

## save
save(AllDataLongx, file = "final_data/03_bugs_algae_flow_joined_by_masterid_noPeakAugs.RData")
