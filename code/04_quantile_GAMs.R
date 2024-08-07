### run the quantile GAMs
library(qgam); library(MASS)
library(tidyverse)
library(sf)
library(gam)
# library(mapview)
library(tidylog)

out.dir <- "final_figures/"

getwd()
# Upload data ---------------------------------------------------

load(file = "final_data/01_bugs_algae_flow_joined_by_masterid.RData")
AllDataLongx <- AllData
head(AllDataLongx)

# Create df of model configurations ---------------------------------------

## bio 
biol.endpoints<-"csci"
biol.endpoints

## flow
flow.endpoints<- unique(na.omit(AllDataLongx$flow_metric))
flow.endpoints

## smoothing functions - testing smoothness
smooth_funcs <- c(3,6,9,12,15,18,20)

## quantiles - testing each quantile
quants <- c(0.99, 0.90, 0.70, 0.50, 0.1, 0.3) ## adding in lower quantile s for Chad

### make grid of all configurations
bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints, 
                             quants = quants, smooth_funcs = smooth_funcs, stringsAsFactors = F)

dim(bio_h_summary) ## 378, 4


# Run the Models ----------------------------------------------------------

i=1 ## to test

## model of each configuration - takes a little while - time for coffee!
gam.lm <-lapply(1:nrow(bio_h_summary), function(i)
{
  ## take each model element
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  qmet<-as.numeric(bio_h_summary[i,"quants"])
  smet<-as.numeric(bio_h_summary[i,"smooth_funcs"])
  
  ## subset from data
  mydat<-AllDataLongx %>%
    filter(Metric == bmet, ## bio metric
           flow_metric == tmet) %>% ## flow metric
    dplyr::select(Metric, MetricValue, deltah_final, channel_engineering_class) %>% ## remove redendant columns
    drop_na( deltah_final) %>% ## drop na from flow
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct()
  
  ## format data
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  
  ## run qgams
  qgam(MetricValue~s(deltah_final, k=smet), ## formula with smoothing
       data = mydat, ## data
       qu = qmet) ## quantiles
  
})

## save models
save(gam.lm, file = "ignore/04_csci_asci_quantile_gam_updatedSites.RData")
gam.lm


# Model coefficients ------------------------------------------------------

## load models
load(file = "ignore/04_csci_asci_quantile_gam_updatedSites.RData")

### get rsqds and pvals

for(i in 1:length(gam.lm)) {
  
  ## define model
  mod <- summary(gam.lm[[i]])
  
  
  bio_h_summary$DevExplained[i] <- mod$dev.expl ##1-mod$deviance/mod$null.deviance ## mcfaddens r2
  bio_h_summary$PValue[i] <- mod$s.table[4] ## p values
  bio_h_summary$R2[i] <- mod$r.sq ## r sq
  bio_h_summary$n[i] <- mod$n ## number of data points
  bio_h_summary$EDF[i] <- mod$edf ## edf
  
}

## save configs and r sqds
save(bio_h_summary, file="final_data/04_csci_asci_quantile_gam_rsqds_updatedSites.RData")


# Plot model performance --------------------------------------------------

## plot quantiles with r2s 

head(bio_h_summary)

P1 <- ggplot(data = bio_h_summary, aes(y=R2, x=quants, group = flow.endpoints, colour = flow.endpoints)) +
  geom_line( linewidth = 1)+
  facet_wrap(~biol.endpoints, scales = "free") +
  facet_grid(rows = vars(smooth_funcs), cols = vars(biol.endpoints), scales = "free") +
  scale_x_continuous(name= "Quantiles") +
  scale_y_continuous(name = "McFaddens (?) Rsq") 

P1

file.name1 <- paste0(out.dir, "04_quants_rsq_smooths_updatedSites.jpg")
ggsave(P1, filename=file.name1, dpi=300, height=5, width=7.5)

## plot k vals with AIC

P2 <- ggplot(data = bio_h_summary, aes(y=DevExplained, x=quants, group = flow.endpoints, colour = flow.endpoints)) +
  geom_line( linewidth = 1)+
  facet_wrap(~biol.endpoints, scales = "free") +
  facet_grid(rows = vars(smooth_funcs), cols = vars(biol.endpoints), scales = "free") +
  scale_x_continuous(name= "Quantiles") +
  scale_y_continuous(name = "Deviance Explained") 

P2

file.name1 <- paste0(out.dir, "04_Quantiles_DevExpl_smooths_updatedSites.jpg")
ggsave(P2, filename=file.name1, dpi=300, height=5, width=7.5)


# Predictions -------------------------------------------------------------

## make df of predicted values to predict on - need to be different for each flow metric
### this takes a loooong time

## blank df
DF <- NULL
DF <- as.data.frame(DF)
## reduce to only peaks etc to test augmentations
# bio_h_summary_sub <- bio_h_summary %>%
#   mutate(Index = rownames(bio_h_summary)) %>%
#   filter(flow.endpoints %in% c("d_peak_10", "d_peak_2", "d_peak_5", "delta_q99"),
#          smooth_funcs %in% c(3,6),
#          quants %in% c(0.5, 0.9, 0.1))
# 
# 
# bio_h_summary_sub$Index
# 
# gam.lm[[4]]
# 
# modInd <- as.numeric(bio_h_summary_sub$Index)
# modInd
# gam.lm_sub <- gam.lm[modInd]
# gam.lm_sub
# flow.endpoints <- flow.endpoints[6:9]
# biol.endpoints <- "csci"
# 
# ## smoothing functions - testing smoothness
# smooth_funcs <- c(3,6)
# 
# ## quantiles - testing each quantile
# quants <- c( 0.90, 0.50, 0.1) ## adding in lower quantile s for Chad

### get predictions and fitted values
for(i in 1:length(gam.lm)) {
  
  print(paste0("Model ", i))
  
  tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
  bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
  qmet<-as.numeric(bio_h_summary[i,"quants"])
  smet<-as.numeric(bio_h_summary[i,"smooth_funcs"])
  
  ## get input data
  mydat<-AllDataLongx %>%
    filter(Metric == bmet,
           flow_metric == tmet) %>%
    dplyr::select(Metric, MetricValue, deltah_final, Class2, channel_engineering_class) %>% ## only metrics needed
    drop_na(deltah_final, MetricValue) %>%
    # filter(deltah_final > 0) %>%
    filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
    distinct()
  
  incr <- mean(mydat$deltah_final)/500
  
  ## get new data from range of iunput data
  newdata <- as.data.frame(seq(min(mydat$deltah_final), max(mydat$deltah_final), abs(incr)))
  colnames(newdata) <- "deltah_final" ## change name to match original
  
  ## get model, predict, extract all data and categories
  mod <- gam.lm[[i]]
  
  predictedVals <- predict(mod, newdata, type = "response") ## add newdata if needed
  
  ## join
  DFX <- cbind(newdata, as.data.frame(predictedVals)) %>%
    # rename(Value = "s(deltah_final)") %>%
    mutate(Variable = tmet,  Metric = bmet, Quant = qmet, Smooths = smet) %>%
    distinct()
  
  # head(DFX)
  
  DF <- bind_rows(DF, DFX)
  
}

head(DF)

## change back to numeric
DF <- DF %>%
  mutate(deltah_final = as.numeric(deltah_final),
         predictedVals = as.numeric(predictedVals),
         Quant = as.factor(Quant))

save(DF, file = "ignore/04_quantGams_smooths_predictions_updatedSites.RData")
str(DF)

load(file = "ignore/04_quantGams_smooths_predictions_updatedSites.RData")

### predicted figures
bio <- unique(DF$Metric)
mets <- unique(DF$Variable)

## plot curve with points for each quant

for(b in 1:length(bio)) {
  
  ## filter to index
  DF1 <- DF %>%
    filter(Metric == bio[b])
  
  ## get thresholds for ref and HB channels
  if(bio[b] == "asci") {
    
    refT <- 0.86
    hbT <- 0.87
    
  } else {
    
    refT <- 0.79
    hbT <- 0.67
  }
  
  # names(AllDataLong2)
  # head(DF1)
  
  ### get observed data points for overlay
  ptsbio <- AllDataLongx %>%
    filter(Metric == bio[b])
  
  # str(ptsbio)
  
  for(m in 1:length(mets)) {
    
    ptsbiox <- ptsbio %>%
      filter(flow_metric == mets[m]) %>%
      mutate(channel_engineering_class = as.factor(channel_engineering_class))
    
    ptsbiox
    
    T1 <- ggplot() +
      geom_smooth(data = subset(DF1, Variable == mets[m]), aes(y=predictedVals, x=deltah_final, col = Quant), linewidth = 1)+
      geom_point(data = ptsbiox, aes(x=deltah_final, y = MetricValue,
                                     col = channel_engineering_class,
                                     size=0.5)) +
      # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
      geom_hline(yintercept = refT,  linetype="dashed", linewidth=0.5, color = "grey50") +
      geom_hline(yintercept = hbT,  linetype="dashed", linewidth=0.5, color = "red") +
      geom_vline(xintercept = 0) +
      facet_wrap(~Smooths, scales = "free") +
      scale_x_continuous(name=paste(mets[m])) +
      scale_y_continuous(name = paste0("Index Score  (", bio[b], ")"))
    
    
    # theme(legend.position = "none"))
    
    T1
    
    file.name1 <- paste0(out.dir, "04_", bio[b], "_", mets[m], "_flow_response_predicted_gam_combined_all_quants_updatedSites.jpg")
    ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
    
  }
  
  
  
  
}


# Simplified Plot ---------------------------------------------------------

## plot curve with NO points 

### predicted figures
bio <- unique(DF$Metric)
mets <- unique(DF$Variable)

## thresholds for ref and HB channels
refT <- 0.79
hbT <- 0.67

DF1 <- DF %>%
  filter(Smooths %in% c(3,6),
         Quant %in% c(0.9,0.5,0.3))

head(DF1)

  for(m in 1:length(mets)) {

    T1 <- ggplot() +
      geom_smooth(data = subset(DF1, Variable == mets[m]), aes(y=predictedVals, x=deltah_final, col = Quant), linewidth = 1)+
      # geom_point(data = ptsbiox, aes(x=deltah_final, y = MetricValue,
      #                                col = channel_engineering_class)) +
      # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
      geom_hline(yintercept = refT,  linetype="dashed", linewidth=0.5, color = "grey50") +
      geom_hline(yintercept = hbT,  linetype="dashed", linewidth=0.5, color = "red") +
      geom_vline(xintercept = 0) +
      facet_wrap(~Smooths, scales = "free") +
      scale_x_continuous(name=paste(mets[m])) +
      scale_y_continuous(name = paste0("Index Score (CSCI)"))
    
    
    # theme(legend.position = "none"))
    
    T1
    
    file.name1 <- paste0(out.dir, "04_CSCI_", mets[m], "_flow_response_predicted_gam_combined_all_quants_updatedSites_simple.jpg")
    ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
    
  }
  
  

# load(file = "ignore/04_quantGams_smooths_predictions_noPeakAugs.RData")
str(DF)

### predicted figures
bio <- unique(DF$Metric)
mets <- unique(DF$Variable)
mets
bio

b=1
m=1
### take the assigned K from script 05, just csci and DS baseflow, only 0.5 & 0.9 

## filter to index and only 
DF1 <- DF %>%
  filter(Metric == bio[b], Quant %in% c( 0.3, 0.5, 0.9), Smooths == 6)

## thresholds
refT <- 0.79
hbT <- 0.67


## get metric
# ptsbiox <- ptsbio %>%
#       filter(flow_metric == mets[m]) %>%
#       mutate(channel_engineering_class = as.factor(channel_engineering_class))

# ptsbiox

T1 <- ggplot() +
  geom_smooth(data = subset(DF1, Variable == mets[m]), aes(y=predictedVals, x=deltah_final, col = Quant), linewidth = 1)+
  # geom_point(data = ptsbiox, aes(x=deltah_final, y = MetricValue,
  #                                col = channel_engineering_class)) +
  # stat_smooth(method = "lm", formula = y ~ x + I(x^2), linewidth = 1)+
  geom_hline(yintercept = refT,  linetype="dashed", linewidth=0.5, color = "grey50") +
  geom_hline(yintercept = hbT,  linetype="dashed", linewidth=0.5, color = "red") +
  geom_vline(xintercept = 0) +
  guides(col=guide_legend(title="Flow Alteration Level")) +
  scale_x_continuous(name= "Dry Season Baseflow") +
  scale_y_continuous(name = paste0("Index Score  (", bio[b], ")"))


# theme(legend.position = "none"))

T1

file.name1 <- paste0(out.dir, "04_", bio[b], "_", mets[m], "_DS_baseflow_example_updatedSites.jpg")
ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)








