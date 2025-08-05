### this is the first step script in identifying signals from the well-specific data
### this script will solely import the input data, pre-process it and assign score
### these scores are based mainly on 4 variables - PercentOccupied, SpotIntensity, SpotArea, BackgroundPenalty

#### package and data import point ####

## importing relevant packages -- same imports for all the steps
library(tidyverse)
library(ggcorrplot)
library(dplyr)
library(ggplot2)
library(viridis)
library(broom)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(magick)
library(readr)
library(tidyr)
library(corrplot)
library(scales)
library(ggpmisc)
library(ggpubr)

## sourcing functions
source("step2_makeCorrelationHeatmap.R")
source("step3_prepareCMVSpikeScores.R")
source("step4_makeBoxplots.R")

## importing the well-specific variable dataset
data = read_tsv("well_data.tsv")

#### pre-processing the import dataset and assigning scores ####

#selecting variables that will be used in the further analysis
feature_cols = c(
  "Count_SpotInWell",
  "AreaOccupied_AreaOccupied_SpotInWell",
  "Mean_SpotInWell_Intensity_MeanIntensity_Grey",
  "Mean_SpotInWell_AreaShape_Area",
  "Mean_Well_Intensity_IntegratedIntensity_Grey",
  "AreaOccupied_AreaOccupied_Well"
)

#designing new variables based on the already available variables
data = data %>%
  mutate(
    PercentOccupied = AreaOccupied_AreaOccupied_SpotInWell / AreaOccupied_AreaOccupied_Well,
    SpotIntensity = Mean_SpotInWell_Intensity_MeanIntensity_Grey,
    SpotArea = Mean_SpotInWell_AreaShape_Area,
    BackgroundPenalty = 1 / Mean_Well_Intensity_IntegratedIntensity_Grey,
    WellIntensity = Mean_Well_Intensity_IntegratedIntensity_Grey
  )

#extracting a subset of the well data and assigning NA or infinite values to zero
features_all = data %>%
  select(Count_SpotInWell, PercentOccupied, SpotIntensity, SpotArea, BackgroundPenalty, WellIntensity) %>%
  mutate(across(everything(), ~ ifelse(is.na(.) | !is.finite(.), 0, .)))

#calculating a proxy score composite of 4 variables - omitting counts of spot in well and well intensity
proxy_composite = features_all %>%
  select(PercentOccupied, SpotIntensity, SpotArea, BackgroundPenalty) %>%
  mutate(across(everything(), scale)) %>%
  rowwise() %>%
  mutate(proxy = mean(c_across(everything()))) %>%
  pull(proxy)

###!! weight definition for variables with a grid system

#grid definition for weighting the 4 variables - PercentOccupied, SpotIntensity, SpotArea, BackgroundPenalty
alpha_vals = seq(0.1, 1.5, by = 0.2)
beta_vals  = seq(0.5, 2.0, by = 0.3)
gamma_vals = seq(0.1, 1.5, by = 0.2)
delta_vals = seq(0.2, 2.0, by = 0.3)

#making a grid of the scores - all combinations of the four grids
grid_expand = expand.grid(alpha = alpha_vals, beta = beta_vals, gamma = gamma_vals, delta = delta_vals)

###!! identifying the combination of weights (grids) for the four variables that correlate best with proximity score

#going over all grid combinations, making score and checking for correlation
best_opt = grid_expand %>%
  rowwise() %>%
  mutate(score = list(
    log1p(features_all$Count_SpotInWell) *
      (features_all$PercentOccupied ^ alpha) *
      (features_all$SpotIntensity ^ beta) *
      (features_all$SpotArea ^ delta) /
      (features_all$BackgroundPenalty ^ gamma)
  ),
  correlation = cor(unlist(score), proxy_composite, method = "spearman")) %>%
  ungroup() %>%
  arrange(desc(correlation)) %>%
  slice(1)

##adding the five scores to the dataframe
data = data %>%
  mutate(
    Score_AllEmphasizeWellOccupancy = log1p(features_all$Count_SpotInWell) *
      (features_all$PercentOccupied ^ best_opt$alpha) *
      (features_all$SpotIntensity ^ best_opt$beta) *
      (features_all$SpotArea ^ best_opt$delta) /
      (features_all$BackgroundPenalty ^ best_opt$gamma),
    
    Score_AllEqualWeights = log1p(features_all$Count_SpotInWell) *
      (features_all$PercentOccupied ^ 1) *
      (features_all$SpotIntensity ^ 1) *
      (features_all$SpotArea ^ 1) /
      (features_all$BackgroundPenalty ^ 1),
    
    Score_AllOptimizedWeights = log1p(features_all$Count_SpotInWell) *
      (features_all$PercentOccupied ^ best_opt$alpha) *
      (features_all$SpotIntensity ^ best_opt$beta) *
      (features_all$SpotArea ^ best_opt$delta) /
      (features_all$BackgroundPenalty ^ best_opt$gamma),
    
    Score_OnlyOccupancyIntensity = features_all$PercentOccupied * features_all$WellIntensity,
    Score_OnlyOccupancy = features_all$PercentOccupied
    
    ## ~~ to check if to keep
    ## ~~ QualityFlag = if_else(rowSums(is.na(select(., all_of(feature_cols)))) > 0, "Low", "High")
  )

##choosing columns that are to be normalized
score_cols = c(
  "Score_AllEmphasizeWellOccupancy",
  "Score_AllEqualWeights",
  "Score_AllOptimizedWeights",
  "Score_OnlyOccupancyIntensity",
  "Score_OnlyOccupancy"
)

##creating a normalized score and a z-score -- we will use the normalized score in the end (0,1)
data = data %>%
  mutate(across(all_of(score_cols), list(
    
    #z-score
    
    Z = ~ {
      x = ifelse(is.na(.) | !is.finite(.), 0, .)
      scale(x)[, 1]
    },
    
    #normalized score
    
    Norm = ~ {
      x <- ifelse(is.na(.) | !is.finite(.), 0, .)
      rng <- range(x)
      if (diff(rng) == 0) rep(0, length(x)) else (x - rng[1]) / diff(rng)
    }
  )))

## writing the dataset along with the scores to a new file
write_tsv(data, "intermediate_files/well_data_withScores.tsv")

#### giving calls to functions for further processing ####

## Step 2 - for making the correlation and heatmaps (which includes image thumbnails)
step2_makeCorrelationHeatmap(data = data, image_base_path = "C:/Users/LindenHome/Desktop/ELISPOT") ## insert image path

## Step 3 - re-processing the data and making the CMV and Spike normalized scores with DMSO
elispot_data = step3_prepareCMVSpikeScores(well_data_scores = data)

## Step 4 - making the CMV and Spike plots with the normalized scores
step4_makeBoxplots(elispot_data = elispot_data)
