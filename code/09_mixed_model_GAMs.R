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
AllDataLongx <- AllData #%>%
  # dplyr::filter(!channel_engineering_class == "SB1")
head(AllDataLongx)

unique(AllDataLongx$channel_engineering_class)

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

dim(bio_h_summary) ## 18, 3

## blank df
DF <- NULL
DF <- as.data.frame(DF)

## blank coefs 
coefs <- NULL
coefs <- as.data.frame(coefs)

i=1
## blank mods
gam_lme <- NULL

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
  drop_na(deltah_final) %>% ## drop na from flow
  filter_all(all_vars(!is.infinite(.))) %>% ## remove all missing values
  distinct() %>%
  mutate(channel_engineering_class = as.factor(channel_engineering_class)) %>%
  drop_na()

# unique(mydat$channel_engineering_class)

## format data
mydat<-mydat[order(mydat$MetricValue),] ## order by csci value
head(mydat)
## run models

### random intercept
gamm_intercept <- mgcv::gam(MetricValue ~ s(deltah_final, k=smet) + 
                        s(channel_engineering_class, bs = "re"), 
                      data = mydat,
                      method = "REML")

### random slopes
gamm_slope <- mgcv::gam(MetricValue ~ s(deltah_final, k=smet) + 
                    s(deltah_final, channel_engineering_class, bs = "re"), 
                  data = mydat,
                  method = "REML")

# Extract standard deviations
# random_std_dev <- vcomp["s(deltah_final,channel_engineering_class)", "std.dev"]  # Random effect for channel_engineering_class
# residual_std_dev <- vcomp["scale", "std.dev"]  # Residual standard deviation

# Compute variances
# random_variance <- random_std_dev^2
# residual_variance <- residual_std_dev^2
# 
# # Compute ICC
# icc <- random_variance / (random_variance + residual_variance)
# print(paste("ICC:", round(icc, digits = 2)))

### random slope & intercept
gamm_int_slope <- mgcv::gam(MetricValue ~ s(deltah_final, k=smet) + 
                        s(channel_engineering_class, bs = "re") + ## intercept
                        s(channel_engineering_class, deltah_final, bs = "re"), ## slope
                        data = mydat, 
                        method = "REML")


# vcomp <- gam.vcomp(gamm_slope)
# 
# icc <- vcomp[["fac"]] / (vcomp[["fac"]] + vcomp[["residual"]])
# 
# print(paste("ICC:", icc))
# print(vcomp)

## mods 

# modsx <- c(gamm_intercept, gamm_slope, gamm_int_slope)
## coefs
coefsx <- AIC(gamm_intercept, gamm_slope, gamm_int_slope) %>% ## get AIC for comparison
  mutate(DevianceExplained = c(summary(gamm_intercept)$dev.expl, summary(gamm_slope)$dev.expl, summary(gamm_int_slope)$dev.expl)) %>% ## deviance explained
  mutate(RSquared = c(summary(gamm_intercept)$r.sq, summary(gamm_slope)$r.sq, summary(gamm_int_slope)$r.sq)) %>%
  mutate(Variable = tmet,  Metric = bmet, Smooths = smet) %>%
  mutate(Model = c("gamm_intercept", "gamm_slope", "gamm_int_slope"))

## predict for figure
mydat_long <- mydat %>%
  mutate(Intercept = predict(gamm_intercept),
         Slope = predict(gamm_slope),
         SlopeAndIntercept = predict(gamm_int_slope)) %>%
  pivot_longer(Intercept:SlopeAndIntercept, names_to = "ModelType", values_to = "predictedVals") %>%
  mutate(ModelType = factor(ModelType, levels = c("Intercept", "Slope", "SlopeAndIntercept"), 
                            labels = c("Intercept", "Slope", "Slope and Intercept"))) %>%
  mutate(channel_engineering_class = factor(channel_engineering_class, 
                                levels = c("NAT", "SB0", "SB1", "SB2", "HB"), 
                                labels = c("NAT", "SB0", "SB1", "SB2", "HB"))) 

head(mydat_long)

T1 <- ggplot(data = mydat_long, aes(y=MetricValue, x=deltah_final, group = channel_engineering_class, col = channel_engineering_class)) +
  geom_smooth( aes(y=predictedVals, x=deltah_final), linewidth = 1)+
  # geom_point(data = mydat, aes(x=deltah_final, y = MetricValue,
  #                                col = channel_engineering_class)) +
  geom_hline(yintercept = 0.79,  linetype="dashed", linewidth=0.5, color = "grey50") +
  geom_hline(yintercept = 0.67,  linetype="dotted", linewidth=0.5, color = "grey50") +
  geom_vline(xintercept = 0) +
  facet_wrap(~ModelType, scales = "free") +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue2", "darkblue", "mediumpurple2", "firebrick3"))+ ## colour of points
  scale_x_continuous(name="Delta (CFS)") +
  scale_y_continuous(name = paste0("CSCI Score")) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        # legend.position = "bottom",
        legend.text=element_text(size=12),
        axis.text.x = element_text(size = 12, angle = 20, vjust = 0.5,hjust = 0.3),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))

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
         channel_engineering_class4 = "SB1",
         channel_engineering_class5 = "SB2",
         channel_engineering_class6 = "HB") %>%
  pivot_longer(channel_engineering_class1:channel_engineering_class6,  ## make channel class longer
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
  

head(DFX1)

DFX1 <- DFX %>%
  pivot_longer(Intercept:SlopeAndIntercept, names_to = "ModelType", values_to = "predictedVals") %>%
  mutate(ModelType = factor(ModelType, levels = c("Intercept", "Slope", "SlopeAndIntercept"), 
                            labels = c("Intercept", "Slope", "Slope and Intercept"))) %>%
  mutate(channel_engineering_class = factor(channel_engineering_class, 
                                            levels = c("NAT", "SB0", "SB1", "SB2", "HB"), 
                                            labels = c("NAT", "SB0", "SB1", "SB2", "HB"))) 



T2 <- ggplot() +
  geom_smooth(data = DFX1, aes(y=predictedVals, x=deltah_final, col = channel_engineering_class), linewidth = 1)+
  geom_hline(yintercept = 0.79,  linetype="dashed", linewidth=0.5, color = "grey50") +
  geom_hline(yintercept = 0.67,  linetype="dotted", linewidth=0.5, color = "grey50") +
  geom_vline(xintercept = 0) +
  scale_colour_manual(values=c("chartreuse4", "dodgerblue2", "darkblue", "mediumpurple2", "firebrick3"))+ ## colour of points
  scale_x_continuous(name="Delta (CFS)") +
  facet_wrap(~ModelType, scales = "free") +
  scale_y_continuous(name = paste0("CSCI Score")) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        # legend.position = "bottom",
        legend.text=element_text(size=12),
        axis.text.x = element_text(size = 12, angle = 20, vjust = 0.5,hjust = 0.3),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))



# theme(legend.position = "none"))

T2

file.name1 <- paste0(out.dir, "09_", tmet, "_", smet,  "_mixed_effects_GAM_predicted.jpg")
ggsave(T2, filename=file.name1, dpi=300, height=5, width=7.5)


# head(DFX)
## combine 
DF <- bind_rows(DF, DFX)
coefs <- bind_rows(coefs, coefsx)
gam_lme <- c(gam_lme, list(gamm_intercept, gamm_slope, gamm_int_slope))
# sapply(gam_lme, class)
}

gam_lme
## remove blank model
# gam_lme <- gam_lme[-1]

## save models
save(gam_lme, file = "final_data/09_mixed_models.RData")

## save coefs
write.csv(coefs, "final_data/09_coefs_mixed_effects_model.csv")
## without SB1 r2 ~ 0.4-0.45
## save predictions

save(DF, file= "ignore/09_mixed_effects_model_predictions.RData")
load(file= "ignore/09_mixed_effects_model_predictions.RData")
head(DF)


# Format coefs for MS -----------------------------------------------------

## full names for labels
labels <- read.csv("input_data/ffm_names.csv")
labels <- labels[1:24, ]
labels <- labels %>% rename(hydro.endpoints = Flow.Metric.Code)
labels[25, 1] <- "Peak Flow Magnitude (Q99, cfs)"
labels[25, 2] <- "Q99"
labels[25, 3] <- "Peak Flow Magnitude"
labels

## upload coefs
coefs <- read.csv("final_data/09_coefs_mixed_effects_model.csv")
head(coefs)

## smoothing functions per FFM
sm6csci <- c("DS_Mag_50", "FA_Mag", "Peak_10", "Peak_2", "Peak_5", "Wet_BFL_Mag_50", "Wet_BFL_Mag_10")

sm3csci <- c( "SP_Mag",  "Q99")

## format data
coefsx <- coefs %>% ## format names
  select(-X) %>%
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

## organise into tables
head(coefsx)

CoefsTab <- coefsx %>%
  select(-c(Metric, hydro.endpoints:ModelstoUse)) %>%
  mutate(Model = case_when(Model == "gamm_intercept" ~ "Intercept Only",            
                           Model == "gamm_slope" ~ "Slope Only",
                           Model == "gamm_int_slope" ~ "Slope and Intercept")) %>%
  select(Flow.Metric.Name, Model, DevianceExplained, RSquared, df, AIC)
  
CoefsTab

write.csv(CoefsTab, "Tables/09_mixed_model_gam_coefs.csv")

# Loop over FFMs: Quantile mixed model Gams -----------------------------------------------------------

## quantiles - testing each quantile
quants <- c(0.90, 0.50) 

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
           channel_engineering_class4 = "SB1",
           channel_engineering_class5 = "SB2",
           channel_engineering_class6 = "HB") %>%
    pivot_longer(channel_engineering_class1:channel_engineering_class6,  ## make channel class longer
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

load(file= "ignore/09_mixed_effects_model_predictions.RData")

head(DF)

## change ffm names 
DF <- DF %>%
  rename(flow.endpoints = Variable) %>%
  mutate(Hydro_endpoint = case_when(flow.endpoints == "d_ds_mag_50" ~ "DS_Mag_50",            
                                  flow.endpoints == "d_fa_mag" ~ "FA_Mag",
                                  flow.endpoints == "d_peak_10" ~ "Peak_10",
                                  flow.endpoints == "d_peak_2" ~ "Peak_2",
                                  flow.endpoints == "d_peak_5" ~ "Peak_5",
                                  flow.endpoints == "d_sp_mag" ~ "SP_Mag",
                                  flow.endpoints == "d_wet_bfl_mag_10" ~ "Wet_BFL_Mag_10",
                                  flow.endpoints == "d_wet_bfl_mag_50" ~ "Wet_BFL_Mag_50", 
                                  flow.endpoints == "delta_q99" ~ "Q99"))

## smoothing functions per FFM
sm6csci <- c("DS_Mag_50", "FA_Mag", "Peak_10", "Peak_2", "Peak_5", "Wet_BFL_Mag_50", "Wet_BFL_Mag_10")

sm3csci <- c( "SP_Mag",  "Q99")

## raw data for points
head(AllDataLongx)

## select columns, rename to match, make channel a factor
pts <- AllDataLongx %>%
  select(channel_engineering_class, hydro.endpoints, Flow.Metric.Name, deltah_final, MetricValue) %>%
  rename(Hydro_endpoint = hydro.endpoints) %>%
  mutate(channel_engineering_class = factor(channel_engineering_class, levels = c("NAT", "SB0","SB1", "SB2", "HB"))) ## define levels of channel type


### plot slope&intercept for all FFMs
## uses only 0.5 quantile


## filter out smoothing functions
DFx <- DF %>%
  pivot_longer(Intercept:SlopeAndIntercept, names_to = "ModelType", values_to = "predictedVals") %>% ## make longer
  mutate(ChosenSmth = case_when((Metric == "csci" & Hydro_endpoint %in% sm6csci) ~ 6,
                                (Metric == "csci" & Hydro_endpoint %in% sm3csci) ~ 3)) %>%
  mutate(ModelstoUse = ifelse(Smooths == ChosenSmth, "Yes", "No")) %>%
  filter(ModelstoUse == "Yes") %>%
  # filter(Quantile %in% c( 0.50)) %>% #
  filter(ModelType == "SlopeAndIntercept") %>%
  # mutate(Quantile = as.factor(Quantile)) %>% ## change quantile to factor
  mutate(channel_engineering_class = factor(channel_engineering_class, levels = c("NAT", "SB0", "SB1", "SB2", "HB"))) ## define levels of channel type

unique(DFx$channel_engineering_class)
## define colours


cols = c("NAT" = "#7CAE00", "SB0" = "#00BFC4", "SB1" = "#00C19A",  "SB2" = "#C77CFF", "HB" = "#F8766D" )
### plot slope and intercept 
T5 <- ggplot() +
  geom_smooth(data = DFx, aes(y=predictedVals, x=deltah_final, col = channel_engineering_class), linewidth = 1)+
  geom_point(data = pts, aes(y = MetricValue, x = deltah_final, col = channel_engineering_class)) +
  geom_hline(yintercept = 0.79,  linetype="dashed", linewidth=0.5, color = "grey50") +
  geom_hline(yintercept = 0.67,  linetype="dotted", linewidth=0.5, color = "grey50") +
  geom_vline(xintercept = 0) +
  facet_wrap(~Hydro_endpoint, scales = "free") +
  scale_color_manual(values = cols) +
  # scale_colour_manual(values=hcl.colors(n=4, palette = "Zissou 1"))+
  # scale_x_continuous(name=paste(mets[m])) +
  scale_y_continuous(name = paste0("CSCI Score"))


T5

file.name1 <- paste0(out.dir, "09_mixed_model_slopeintercept_all_FFM_predicted.jpg")
ggsave(T5, filename=file.name1, dpi=300, height=5, width=7.5)

### plot slope and intercept, no dots
T6 <- ggplot() +
  geom_smooth(data = DFx, aes(y=predictedVals, x=deltah_final, col = channel_engineering_class), linewidth = 1)+
  # geom_point(data = pts, aes(y = MetricValue, x = deltah_final, col = channel_engineering_class)) +
  geom_hline(yintercept = 0.79,  linetype="dashed", linewidth=0.5, color = "grey50") +
  geom_hline(yintercept = 0.67,  linetype="dotted", linewidth=0.5, color = "grey50") +
  geom_vline(xintercept = 0) +
  facet_wrap(~Hydro_endpoint, scales = "free") +
  scale_color_manual(values = cols) +
  # scale_colour_manual(values=hcl.colors(n=4, palette = "Zissou 1"))+
  # scale_x_continuous(name=paste(mets[m])) +
  scale_y_continuous(name = paste0("CSCI Score"))


T6

file.name1 <- paste0(out.dir, "09_mixed_model_slopeintercept_all_FFM_predicted_noDots.jpg")
ggsave(T6, filename=file.name1, dpi=300, height=5, width=7.5)

### plot only dots
T7 <- ggplot() +
  # geom_smooth(data = DFx, aes(y=predictedVals, x=deltah_final, col = channel_engineering_class), linewidth = 1)+
  geom_point(data = pts, aes(y = MetricValue, x = deltah_final, col = channel_engineering_class)) +
  geom_hline(yintercept = 0.79,  linetype="dashed", linewidth=0.5, color = "grey50") +
  geom_hline(yintercept = 0.67,  linetype="dotted", linewidth=0.5, color = "grey50") +
  geom_vline(xintercept = 0) +
  facet_wrap(~Hydro_endpoint, scales = "free") +
  scale_color_manual(values = cols) +
  # scale_colour_manual(values=hcl.colors(n=4, palette = "Zissou 1"))+
  scale_x_continuous(name=paste(mets[m])) +
  scale_y_continuous(name = paste0("CSCI Score"))


T7

file.name1 <- paste0(out.dir, "09_observed_data_dots.jpg")
ggsave(T7, filename=file.name1, dpi=300, height=5, width=7.5)


#### other plots
## filter out smoothing functions
DFx <- DF %>%
  pivot_longer(Intercept:SlopeAndIntercept, names_to = "ModelType", values_to = "predictedVals") %>% ## make longer
  mutate(ChosenSmth = case_when((Metric == "csci" & Hydro_endpoint %in% sm6csci) ~ 6,
                                (Metric == "csci" & Hydro_endpoint %in% sm3csci) ~ 3)) %>%
  mutate(ModelstoUse = ifelse(smooth_funcs == ChosenSmth, "Yes", "No")) %>%
  filter(ModelstoUse == "Yes") #%>%
# filter(quants %in% c( 0.50))
filter(ModelType == "Intercept") %>%
  mutate(Quantile = as.factor(Quantile)) %>% ## change quantile to factor
  mutate(channel_engineering_class = factor(channel_engineering_class, levels = c("NAT", "SB0", "SB2", "HB"))) ## define levels of channel type

str(DFx)
head(DFx)
  
## define ffms
mets <- unique(DFx$Hydro_endpoint)

## plot curve for each quant per channel type

m=1
for(m in 1:length(mets)) {
  
  
  T3 <- ggplot() +
    geom_smooth(data = subset(DFx, Hydro_endpoint == mets[m]), aes(y=predictedVals, x=deltah_final, col = Quantile), linewidth = 1)+
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
    geom_smooth(data = subset(DFx1, Hydro_endpoint == mets[m]), aes(y=predictedVals, x=deltah_final, col = channel_engineering_class), linewidth = 1)+
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




