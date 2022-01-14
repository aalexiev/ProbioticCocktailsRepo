## Graphing for OD600 data from CFS a/s assay on 2.7.21

library(ggplot2)
library(dplyr)
library(tidyr)
setwd("~/Documents/Projects/ProbioticCocktails/02_PreliminaryAnalysis/CFS as trial 2.7.21")


## load in concatenated data
CFS <- read.table("CFS_as_2.7.21_metadataLONG.txt", header = T, sep = "")

## Check if replicates were good replicates of each other
# will make box plots without collapsed replicates
levels(CFS$sampleID)
# no-Bd controls compared to live Bd
CFS_Bd <- CFS %>%
    dplyr::filter(sampleID %in% c("Live_Bd", "Heat_Bd", "SterileMedium"))
ggplot(CFS_Bd, aes(x = HoursElapsed, y = OD, color = sampleID)) +
    geom_point() +
    geom_line(aes(linetype = sampleID))


CFS_BTB_47 <- CFS %>%
    dplyr::filter(sampleID %in% c("live Bd", "heat-killed Bd", "BTB_47", "BTB_47+Bd"))
ggplot(CFS_BTB_47, aes(x = HoursSinceStart, y = OD600, color = sampleID)) +
    geom_point() +
    geom_line(aes(linetype = sampleID))


# Jliv treatments and controls
CFS_jliv <- CFS %>%
    dplyr::filter(sampleID %in% c("live Bd", "heat-killed Bd", "Jliv", "jliv+Bd"))
ggplot(CFS_jliv, aes(x = HoursSinceStart, y = OD600, color = sampleID)) +
    geom_point() +
    geom_line(aes(linetype = sampleID))


# combo treatments and controls
CFS_combo <- CFS %>%
    dplyr::filter(sampleID %in% c("live Bd", "heat-killed Bd", "BTB_47+Jliv", "BTB_47+Jliv+Bd"))
ggplot(CFS_combo, aes(x = HoursSinceStart, y = OD600, color = sampleID)) +
    geom_point() +
    geom_line(aes(linetype = sampleID))
# Take aways:
# again, combo treatment inhibited Bd only after time 4 (day 5)
# variation in the points increases as time goes on