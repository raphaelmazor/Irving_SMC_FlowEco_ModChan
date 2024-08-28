## trying mixed model gam
# install.packages("itsadug")
### run the quantile GAMs
library(qgam); library(MASS)
library(tidyverse)
library(sf)
library(mgcv)
library(itsadug)
# library(mapview)
library(tidylog)


## data
load(file = "final_data/01_bugs_algae_flow_joined_by_masterid.RData") 

## remove SB1
AllDataLongx <- AllData %>%
  dplyr::filter(!channel_engineering_class == "SB1")
head(AllDataLongx)

## out directory
out.dir <- "final_figures/"

# Create df of model configurations ---------------------------------------

## bio 
biol.endpoints<-"csci"
biol.endpoints

## flow
flow.endpoints<- unique(na.omit(AllDataLongx$flow_metric))
flow.endpoints

## smoothing functions - testing smoothness
smooth_funcs <- c(3,6)

## quantiles - testing each quantile
# quants <- c(0.99, 0.90, 0.70, 0.50) ## adding in lower quantile s for Chad

### make grid of all configurations
bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints, 
                            smooth_funcs = smooth_funcs, stringsAsFactors = F)

dim(bio_h_summary) ## 18, 4

## blank df
DF <- NULL
DF <- as.data.frame(DF)

## blank coefs 
coefs <- NULL
coefs <- as.data.frame(coefs)

i=1


# Loop over FFMs: Mixed Model Gams ----------------------------------------------------------


for(i in 1:nrow(bio_h_summary)) {
  
## take each model element
tmet<-as.character(bio_h_summary[i,"flow.endpoints"])
bmet<-as.character(bio_h_summary[i,"biol.endpoints"])
# qmet<-as.numeric(bio_h_summary[i,"quants"])
smet<-as.numeric(bio_h_summary[i,"smooth_funcs"])

## subset from data
mydat<-AllDataLongx %>%
  filter(Metric == bmet, ## bio metric
         flow_metric == tmet) %>% ## flow metric
  dplyr::select(Metric, MetricValue, deltah_final, channel_engineering_class) %>% ## remove redendant columns
  drop_na( deltah_final) %>% ## drop na from flow
  filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
  distinct() %>%
  mutate(channel_engineering_class = as.factor(channel_engineering_class)) %>%
  drop_na()

## format data
mydat<-mydat[order(mydat$MetricValue),] ## order by csci value

## run models

### random intercept
gamm_intercept <- mgcv::gam(MetricValue ~ s(deltah_final, k=smet) + 
                        s(channel_engineering_class, bs = "re"), 
                      data = mydat,
                      method = "REML")

### random slopes
gamm_slope <- gam(MetricValue ~ s(deltah_final, k=smet) + 
                    s(deltah_final, channel_engineering_class, bs = "re"), 
                  data = mydat,
                  method = "REML")

### random slope & intercept
gamm_int_slope <- gam(MetricValue ~ s(deltah_final, k=smet) + 
                        s(channel_engineering_class, bs = "re") + 
                        s(channel_engineering_class, deltah_final, bs = "re"), 
                        data = mydat, 
                        method = "REML")



## coefs
coefsx <- AIC(gamm_intercept, gamm_slope, gamm_int_slope) %>% ## get AIC for comparison
  mutate(DevianceExplained = c(summary(gamm_intercept)$dev.expl, summary(gamm_slope)$dev.expl, summary(gamm_int_slope)$dev.expl)) %>% ## deviance explained
  mutate(RSquared = c(summary(gamm_intercept)$r.sq, summary(gamm_slope)$r.sq, summary(gamm_int_slope)$r.sq)) %>%
  mutate(Variable = tmet,  Metric = bmet, Smooths = smet)


## predict for figure
mydat_long <- mydat %>%
  mutate(Intercept = predict(gamm_intercept),
         Slope = predict(gamm_slope),
         SlopeAndIntercept = predict(gamm_int_slope)) %>%
  pivot_longer(Intercept:SlopeAndIntercept, names_to = "ModelType", values_to = "predictedVals")


T1 <- ggplot(data = mydat_long, aes(y=MetricValue, x=deltah_final, group = channel_engineering_class, col = channel_engineering_class)) +
  geom_smooth( aes(y=predictedVals, x=deltah_final), linewidth = 1)+
  # geom_point(data = mydat, aes(x=deltah_final, y = MetricValue,
  #                                col = channel_engineering_class)) +
  geom_hline(yintercept = 0.79,  linetype="dashed", linewidth=0.5, color = "grey50") +
  geom_hline(yintercept = 0.67,  linetype="dotted", linewidth=0.5, color = "grey50") +
  geom_vline(xintercept = 0) +
  facet_wrap(~ModelType, scales = "free") +
  scale_x_continuous(name=paste(tmet)) +
  scale_y_continuous(name = paste0("CSCI Score")) +
  theme(legend.title = element_blank(), 
        # legend.position = "bottom",
        legend.text=element_text(size=15),
        axis.text.x = element_text(size = 15, angle = 20, vjust = 0.5,hjust = 0.3),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

T1

file.name1 <- paste0(out.dir, "09_", tmet, "_", smet,  "_mixed_effects_GAM_raw.jpg")
ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)


## increments for newdata
incr <- mean(mydat$deltah_final)/500

## get new data from range of input data, incluyding channel class for random effects
newdata <- as.data.frame(seq(min(mydat$deltah_final), max(mydat$deltah_final), abs(incr))) %>%
  rename(deltah_final = 1) %>% ## change name to match original
  mutate(channel_engineering_class1 = "NAT", ## add channel class as columns
         channel_engineering_class2 = "SB0",
         channel_engineering_class4 = "SB2",
         channel_engineering_class5 = "HB") %>%
  pivot_longer(channel_engineering_class1:channel_engineering_class5,  ## make channel class longer
               names_to = "deletethiscolumn", values_to = "channel_engineering_class") %>%
  select(-deletethiscolumn) ## remove unwanted column

# head(newdata)

## predict on all 3 models
predictedValsI <- as.data.frame(predict(gamm_intercept, newdata, type = "response")) #%>% ## add newdata 
predictedValsS = as.data.frame(predict(gamm_slope, newdata, type = "response"))
predictedValsIS  = as.data.frame(predict(gamm_int_slope, newdata, type = "response"))

## join all predictions
predictedVals <- cbind(predictedValsI,predictedValsS,predictedValsIS)

## change names
colnames(predictedVals) = c("Intercept", "Slope", "SlopeAndIntercept")

# head(predictedVals)
## join 
DFX <- cbind(newdata, as.data.frame(predictedVals)) %>%
  # rename(Value = "s(deltah_final)") %>%
  mutate(Variable = tmet,  Metric = bmet, Smooths = smet) %>%
  distinct()

head(DFX)

DFX1 <- DFX %>%
  pivot_longer(Intercept:SlopeAndIntercept, names_to = "ModelType", values_to = "predictedVals")


T2 <- ggplot() +
  geom_smooth(data = DFX1, aes(y=predictedVals, x=deltah_final, col = channel_engineering_class), linewidth = 1)+
  geom_hline(yintercept = 0.79,  linetype="dashed", linewidth=0.5, color = "grey50") +
  geom_hline(yintercept = 0.67,  linetype="dotted", linewidth=0.5, color = "grey50") +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 0) +
  facet_wrap(~ModelType, scales = "free") +
  scale_x_continuous(name=paste(tmet)) +
  scale_y_continuous(name = paste0("CSCI Score"))


# theme(legend.position = "none"))

T2

file.name1 <- paste0(out.dir, "09_", tmet, "_", smet,  "_mixed_effects_GAM_predicted.jpg")
ggsave(T2, filename=file.name1, dpi=300, height=5, width=7.5)


# head(DFX)
## combine 
DF <- bind_rows(DF, DFX)
coefs <- bind_rows(coefs, coefsx)


}

## save coefs
write.csv(coefs, "final_data/09_coefs_mixed_effects_model.csv")

## save predictions

save(DF, file= "ignore/09_mixed_effects_model_predictions.RData")


# Loop over FFMs: Quantile mixed model Gams -----------------------------------------------------------

## quantiles - testing each quantile
quants <- c( 0.90, 0.50) 

### make grid of all configurations
bio_h_summary<-  expand.grid(biol.endpoints=biol.endpoints,flow.endpoints=flow.endpoints, 
                             quants = quants, smooth_funcs = smooth_funcs, stringsAsFactors = F)

dim(bio_h_summary) ## 36, 4

## blank df
DF <- NULL
DF <- as.data.frame(DF)

## blank coefs 
coefs <- NULL
coefs <- as.data.frame(coefs)

i=28

for(i in 1:nrow(bio_h_summary)) {
  
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
    distinct() %>%
    mutate(channel_engineering_class = as.factor(channel_engineering_class)) %>%
    drop_na()
  
  ## format data
  mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
  
  ## run models
  
  ### random intercept
  gamm_intercept <- qgam(MetricValue ~ s(deltah_final, k=smet) + 
                                s(channel_engineering_class, bs = "re"), 
                              data = mydat,
                              qu = qmet)
  
  ### random slopes
  gamm_slope <- qgam(MetricValue ~ s(deltah_final, k=smet) + 
                      s(deltah_final, channel_engineering_class, bs = "re"), 
                    data = mydat,
                    qu = qmet)
  
  ### random slope & intercept
  gamm_int_slope <- qgam(MetricValue ~ s(deltah_final, k=smet) + 
                          s(channel_engineering_class, bs = "re") + 
                          s(channel_engineering_class, deltah_final, bs = "re"), 
                        data = mydat,
                        qu = qmet)
  
  
  
  ## coefs
  coefsx <- AIC(gamm_intercept, gamm_slope, gamm_int_slope) %>% ## get AIC for comparison
    mutate(DevianceExplained = c(summary(gamm_intercept)$dev.expl, summary(gamm_slope)$dev.expl, summary(gamm_int_slope)$dev.expl)) %>% ## deviance explained
    mutate(RSquared = c(summary(gamm_intercept)$r.sq, summary(gamm_slope)$r.sq, summary(gamm_int_slope)$r.sq)) %>%
    mutate(Variable = tmet,  Metric = bmet, Smooths = smet, Quants = qmet)
  
  
  ## predict for figure
  mydat_long <- mydat %>%
    mutate(Intercept = predict(gamm_intercept),
           Slope = predict(gamm_slope),
           SlopeAndIntercept = predict(gamm_int_slope),
           Quantile = qmet) %>%
    pivot_longer(Intercept:SlopeAndIntercept, names_to = "ModelType", values_to = "predictedVals")
  
  head(mydat_long)
  
  
  T1 <- ggplot(data = mydat_long, aes(y=MetricValue, x=deltah_final, group = channel_engineering_class, col = channel_engineering_class)) +
    geom_smooth( aes(y=predictedVals, x=deltah_final), linewidth = 1)+
    # geom_point(data = mydat, aes(x=deltah_final, y = MetricValue,
    #                                col = channel_engineering_class)) +
    geom_hline(yintercept = 0.79,  linetype="dashed", linewidth=0.5, color = "grey50") +
    geom_hline(yintercept = 0.67,  linetype="dotted", linewidth=0.5, color = "grey50") +
    geom_vline(xintercept = 0) +
    facet_wrap(~ModelType, scales = "free") +
    scale_x_continuous(name=paste(tmet)) +
    scale_y_continuous(name = paste0("CSCI Score")) +
    theme(legend.title = element_blank(), 
          # legend.position = "bottom",
          legend.text=element_text(size=15),
          axis.text.x = element_text(size = 15, angle = 20, vjust = 0.5,hjust = 0.3),
          axis.text.y = element_text(size = 15),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15))
  
  T1
  
  file.name1 <- paste0(out.dir, "09_", tmet, "_", smet, "_", qmet,  "_mixed_effects_Quantile_GAM_raw.jpg")
  ggsave(T1, filename=file.name1, dpi=300, height=5, width=7.5)
  
  ## increments for newdata
  incr <- mean(mydat$deltah_final)/500
  
  ## get new data from range of input data, incluyding channel class for random effects
  newdata <- as.data.frame(seq(min(mydat$deltah_final), max(mydat$deltah_final), abs(incr))) %>%
    rename(deltah_final = 1) %>% ## change name to match original
    mutate(channel_engineering_class1 = "NAT", ## add channel class as columns
           channel_engineering_class2 = "SB0",
           channel_engineering_class4 = "SB2",
           channel_engineering_class5 = "HB") %>%
    pivot_longer(channel_engineering_class1:channel_engineering_class5,  ## make channel class longer
                 names_to = "deletethiscolumn", values_to = "channel_engineering_class") %>%
    select(-deletethiscolumn) ## remove unwanted column
  
  # head(newdata)
  
  ## predict on all 3 models
  predictedValsI <- as.data.frame(predict(gamm_intercept, newdata, type = "response")) #%>% ## add newdata 
  predictedValsS = as.data.frame(predict(gamm_slope, newdata, type = "response"))
  predictedValsIS  = as.data.frame(predict(gamm_int_slope, newdata, type = "response"))
  
  ## join all predictions
  predictedVals <- cbind(predictedValsI,predictedValsS,predictedValsIS)
  
  ## change names
  colnames(predictedVals) = c("Intercept", "Slope", "SlopeAndIntercept")
  
  # head(predictedVals)
  ## join 
  DFX <- cbind(newdata, as.data.frame(predictedVals)) %>%
    # rename(Value = "s(deltah_final)") %>%
    mutate(Variable = tmet,  Metric = bmet, Smooths = smet, Quantile = qmet) %>%
    distinct()
  
  head(DFX)
  
  DFX1 <- DFX %>%
    pivot_longer(Intercept:SlopeAndIntercept, names_to = "ModelType", values_to = "predictedVals")
  
  
  T2 <- ggplot() +
    geom_smooth(data = DFX1, aes(y=predictedVals, x=deltah_final, col = channel_engineering_class), linewidth = 1)+
    geom_hline(yintercept = 0.79,  linetype="dashed", linewidth=0.5, color = "grey50") +
    geom_hline(yintercept = 0.67,  linetype="dotted", linewidth=0.5, color = "grey50") +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = 0) +
    facet_wrap(~ModelType, scales = "free") +
    scale_x_continuous(name=paste(tmet)) +
    scale_y_continuous(name = paste0("CSCI Score"))
  
  
  # theme(legend.position = "none"))
  
  T2
  
  file.name1 <- paste0(out.dir, "09_", tmet, "_", smet, "_", qmet,  "_mixed_effects_GAM_predicted.jpg")
  ggsave(T2, filename=file.name1, dpi=300, height=5, width=7.5)
  
  
  # head(DFX)
  ## combine 
  DF <- bind_rows(DF, DFX)
  coefs <- bind_rows(coefs, coefsx)
  
}

## save coefs
write.csv(coefs, "final_data/09_coefs_mixed_effects_model_Quantile.csv")

## save predictions

save(DF, file= "ignore/09_mixed_effects_model_predictions_Quantile.RData")


# Quantile Figures --------------------------------------------------------

load(file= "ignore/09_mixed_effects_model_predictions_Quantile.RData")

head(DF)

## filter out smoothing function 6 for now
DFx <- DF %>%
  pivot_longer(Intercept:SlopeAndIntercept, names_to = "ModelType", values_to = "predictedVals") %>% ## make longer
  filter(Smooths %in% c(6), ModelType == "Intercept") %>%
  mutate(Quantile = as.factor(Quantile)) %>% ## change quantile to factor
  mutate(channel_engineering_class = factor(channel_engineering_class, levels = c("NAT", "SB0", "SB2", "HB"))) ## define levels of channel type

str(DFx)
  
## define ffms
mets <- unique(DF$Variable)

## plot curve for each quant per channel type

m=1
for(m in 1:length(mets)) {
  
  
  T3 <- ggplot() +
    geom_smooth(data = subset(DFx, Variable == mets[m]), aes(y=predictedVals, x=deltah_final, col = Quantile), linewidth = 1)+
    geom_hline(yintercept = 0.79,  linetype="dashed", linewidth=0.5, color = "grey50") +
    geom_hline(yintercept = 0.67,  linetype="dotted", linewidth=0.5, color = "grey50") +
    geom_vline(xintercept = 0) +
    facet_wrap(~channel_engineering_class) +
    scale_x_continuous(name=paste(mets[m])) +
    scale_y_continuous(name = paste0("CSCI Score"))
  
  
  T3
  
  file.name1 <- paste0(out.dir, "09_", mets[m], "_mixed_model_Quants_predicted.jpg")
  ggsave(T3, filename=file.name1, dpi=300, height=5, width=7.5)
  
}


## plot channel types, only intercept and 0.5 quant

DFx1 <- DFx %>% 
  filter(Quantile == 0.5) ## filter only 0.5

for(m in 1:length(mets)) {
  
  
  T4 <- ggplot() +
    geom_smooth(data = subset(DFx1, Variable == mets[m]), aes(y=predictedVals, x=deltah_final, col = channel_engineering_class), linewidth = 1)+
    geom_hline(yintercept = 0.79,  linetype="dashed", linewidth=0.5, color = "grey50") +
    geom_hline(yintercept = 0.67,  linetype="dotted", linewidth=0.5, color = "grey50") +
    geom_vline(xintercept = 0) +
    # facet_wrap(~channel_engineering_class) +
    scale_x_continuous(name=paste(mets[m])) +
    scale_y_continuous(name = paste0("CSCI Score"))
  
  
  T4
  
  file.name1 <- paste0(out.dir, "09_", mets[m], "_mixed_model_intercept_predicted.jpg")
  ggsave(T4, filename=file.name1, dpi=300, height=5, width=7.5)
  
}

