# Irving_SMC_FlowEco_ModChan
Flow ecology analysis in modified channels of southern California

Scripts

## 01_data_formatting.R

combines delta FFM (flow data) with bio data. uses both csci and asci

produces file - 01_bugs_algae_flow_joined_by_masterid.RData

## 02_FFM_boxplots.R

creates boxplots of all FFMs per modified stream type
Example naming cnvention - 02_Fall pulse magnitude_boxplot_delta_range

## 03_remove_augmented_peak_flows.R

removes all delta H augmented delta h values for peak flows
After exploratory analysis the relationship between augmented peak flows and csci tended
to be in an unexpected direction

uses data from script 01 - 01_bugs_algae_flow_joined_by_masterid.RData
produces file - 03_bugs_algae_flow_joined_by_masterid_noPeakAugs.RData
produces figures with raw data points - csci as a function of delta ffm, includes all data
03_raw_data_plot_CSCI.jpg

## 04_quantile_GAMs.R

runs the flow ecology analysis
tests various quantiles and smoothing Ks for the GAM
predicts csci scores on the observed data range

produces files - 
04_csci_asci_quantile_gam_noPeakAugs.RData - models
04_csci_asci_quantile_gam_rsqds_noPeakAugs.RData - model coeficients
04_quantGams_smooths_predictions_noPeakAugs.RData - predictions

produces figures
04_quants_rsq_smooths_noPeakAugs.jpg - rsqs with different smoothing ks and quantiles
04_Quantiles_DevExpl_smooths_noPeakAugs.jpg - deviance explained with different smoothing ks and quantiles
"index & FFM" _flow_response_predicted_gam_combined_all_quants_noPeakAugs.jpg - predictions with quantiles 0.5 & 0.9
04_csci_DS_baseflow_example_noPeakAugs.jpg - simplified pot for illustrations

## 05_Gam_thresholds.R

calculates flow ranges from modifed csci thresholds using the predictions from script 04 - 04_quantGams_smooths_predictions_noPeakAugs.RData

uses root_interpolation_function.Rdata function to calculate which takes a csci score amnd finds the corresponding ffm where it intersects the curve for both the 0.5 & 0.9 quantiles

takes a subset of model chosen through visualisation and expert judgement
alot of formatting and data wrangling happens in this script to make sure we have the correct limits in the correct place. take note of script notes

catergorises each site into - healthy/unhealthy biology and unlikely/likely/very likely stress by flow alteration
tallies up the results

also counts number the ffm per site that are likely or very likely stressed by flow alteration

produces files - 
05_delta_thresholds_GAMs.csv - all models
05_chosen_models_thresholds.csv - chosen models
05_QGAM_FFM_Ranges_With_Adjustments.csv - chosen models and adjusted ranges after formatting
05_impact_ffm_bio_ffm_limits.csv - ranges of ffm for each mod channel type
05_ffm_ranges.csv - ranges simple df used later
05_impact_ffm_bio.csv -tidy table of ranges
05_count_impact.csv - tally of categories oer mod channel type
05_Number_ffm_per_result.csv - number of stressed ffm per site

## 06_modified_site_specific_analysis.R

all based on modifed analysis
uses categories from script 05 and counts amount and percentages per ffm
creates maps of mod channel types in each category

produces files - 
06_number_of_sites_each_class_per_FFM.csv - number of sites in each category
06_percent_impacts_each_Class.csv - percent sites in each category
06_percent_impacts_each_Class_wide.csv - percent sites in wide format

produces maps - 
"index_ffm_" _map_6_cats_per_bioHealth_likelihood_flow_impact_per_result - facets by category, colour coded by channel type
"index_ffm_" _map_6_cats_per_bioHealth_likelihood_flow_impact_per_channel_type - facets by channel type, colours by category

produces plots - 
06_stikes_boxplot.jpg - numner of ffm per site er category
06_all_hydro_stacked_perc.jpg - stacked columns with % site per category
06_DS_baseflow_stacked_perc.jpg - stacked column dry season
06_WS_baseflow_stacked_perc.jpg - stacked column wet season
06_peak_flows_stacked_perc.jpg - stacked column peak flows

## 06a_standard_site_specific_analysis.R

same as 06 but using only standard thresholds - 0.79 - for comparison
some data wrangling to get the thresholds is included
need to upload some earlier data to extract the info


## 07_tables_data_packets.R

Tidy and intuitive tables to send out to SMC members for modified analysis


## 07a_tables_data_packets_standard.R

same as 07 but for standard 0.79 analysis