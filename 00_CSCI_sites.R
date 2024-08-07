############################
#
# Purpose: Find bioassessment data for sites in channel engineering table, for Katie.
#          Add masterid, comid & channel class to channel engineering data for export.
#
#
# July 2024
############################

library(tidyverse)

#####-- Get data --#####

con <- DBI::dbConnect(RPostgreSQL::PostgreSQL(), host = 'geobiology.cottkh4djef2.us-west-2.rds.amazonaws.com', port = 5432, dbname = 'smc', user = 'smcread', password = '1969$Harbor')

lustations.1

lustations_query = "select * from sde.lu_stations"           # Station information
tbl_lustations   = tbl(con, sql(lustations_query))
lustations.1     = as.data.frame(tbl_lustations)
rm(lustations_query, tbl_lustations)

csci_query = "select * from analysis_csci_core"              # CSCI data
tbl_csci   = tbl(con, sql(csci_query))
smc.csci     = as.data.frame(tbl_csci)
rm(csci_query, tbl_csci)

asci_query = "select * from analysis_asci"                   # asci data
tbl_asci   = tbl(con, sql(asci_query))
smc.asci     = as.data.frame(tbl_asci)
rm(asci_query, tbl_asci)

ce_query = "select * from unified_channelengineering"        # Channel Engineering data
tbl_ce   = tbl(con, sql(ce_query))
smc.ChEng     = as.data.frame(tbl_ce)
rm(ce_query, tbl_ce)

ces_query = "select * from unified_channelengineering_summary"        # Channel Engineering data
tbl_ces   = tbl(con, sql(ces_query))
smc.ChEngSummary     = as.data.frame(tbl_ces)
rm(ces_query, tbl_ces)


#####-- Conform, combine --#####
smc.ChEng2 <- smc.ChEng %>%
  left_join(lustations.1[, c("stationid", "masterid", "huc", "county", "smcshed", "comid")], by=c("stationcode"="stationid")) %>%
  left_join(smc.ChEngSummary[, c("masterid", "channel_engineering_class")], by=c("masterid"="masterid"))

smc.asci2 <- smc.asci %>%
  left_join(lustations.1[, c("stationid", "masterid", "huc", "county", "smcshed", "comid")], by=c("stationcode"="stationid")) %>%
  filter(masterid %in% smc.ChEng2$masterid)

smc.csci2 <- smc.csci %>%
  left_join(lustations.1[, c("stationid", "masterid","latitude","longitude", "huc", "county", "smcshed", "comid")], by=c("stationcode"="stationid")) #%>%
  # filter(masterid %in% smc.ChEng2$masterid)

#write.csv(smc.ChEng2, "C:/Users/Jeffb/SCCWRP/SMC Stream Survey - SMC 2021-2025 redesign/Data/Working/DataRequests/KatieIrving/ChannelEngineering_071124.csv")
#write.csv(smc.asci2, "C:/Users/Jeffb/SCCWRP/SMC Stream Survey - SMC 2021-2025 redesign/Data/Working/DataRequests/KatieIrving/ASCI_at_ChannelEngineering_070924.csv")
#write.csv(smc.csci2, "C:/Users/Jeffb/SCCWRP/SMC Stream Survey - SMC 2021-2025 redesign/Data/Working/DataRequests/KatieIrving/CSCI_at_ChannelEngineering_070924.csv")


# Format Data -------------------------------------------------------------
names(smc.csci2)
csci <- smc.csci2 %>%
  select(masterid:comid, csci,sampledate, sampleyear, fieldreplicate)

write.csv(csci, "ignore/CSCI_CA_Aug2024.csv")

save(csci, file = "ignore/CSCI_CA_Aug2024.RData")
