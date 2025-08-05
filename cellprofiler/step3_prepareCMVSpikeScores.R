### this is the third step script in identifying signals from the well-specific data
### this script will import the pre-processed dataset (which contains the scores)
### and prepare the files for making the spike and CMV specific plots

#writing a function that will take in the processed data from step 1 and 
#make the normalized scores for CMV and Spike
step3_prepareCMVSpikeScores = function(well_data_scores){
  
  #making a list of the score columns that will be used
  score_cols = c(
    "Score_AllEmphasizeWellOccupancy_Norm", "Score_AllEqualWeights_Norm",
    "Score_AllOptimizedWeights_Norm", "Score_OnlyOccupancyIntensity_Norm",
    "Score_OnlyOccupancy_Norm", 
    "Score_AllEmphasizeWellOccupancy_Z", "Score_AllEqualWeights_Z",
    "Score_AllOptimizedWeights_Z", "Score_OnlyOccupancyIntensity_Z",
    "Score_OnlyOccupancy_Z"
  )
  
  #listing the meta data columns
  meta_cols = c(
    "Age", "Gender", "Exc_Inc", "VoiceID", "Timepoint", "Antigen",
    "Cohort", "WHO", "Tumor stage", "Treatment intent", "Primary tumor localisation"
  )
  
  #combining spike 1 and spike 2 into a single columns
  well_data_scores = well_data_scores %>%
    mutate(Antigen = ifelse(Antigen %in% c("spike1", "spike2"), "spike", Antigen))
  
  #grouping individuals based on the meta data columns
  meta_data = well_data_scores %>%
    select(all_of(meta_cols)) %>%
    group_by(VoiceID, Timepoint, Antigen) %>%
    summarise(across(everything(), ~ first(.x)), .groups = "drop")  # Avoid duplication
  
  #taking the average of the replicates grouped on VOICE ID, timepoint and antigen (in that order)
  avg_replicates = well_data_scores %>%
    group_by(VoiceID, Timepoint, Antigen) %>%
    summarise(across(all_of(score_cols), mean, na.rm = TRUE, .names = "Avg_{.col}"), .groups = "drop")
  
  #extracting the dmso values which will be used for normalization
  dmso_avg = avg_replicates %>%
    filter(Antigen == "dmso") %>%
    select(VoiceID, Timepoint, starts_with("Avg_")) %>%
    rename_with(~ paste0("DMSO_", .), .cols = starts_with("Avg_"))
  
  #joining the dmso scores back to the average of the replicates data frame
  normalized_scores = avg_replicates %>%
    left_join(dmso_avg, by = c("VoiceID", "Timepoint"))
  
  #substracting the dmso values from the average replicates values
  for (score in score_cols) {
    norm_col = paste0("Norm_", score)
    avg_col = paste0("Avg_", score)
    dmso_col = paste0("DMSO_Avg_", score)
    normalized_scores[[norm_col]] = normalized_scores[[avg_col]] - normalized_scores[[dmso_col]]
  }
  
  #merging the normalized scores back with the metadata dataframe
  meta_data_withscores = normalized_scores %>%
    filter(Antigen != "dmso") %>%
    left_join(meta_data, by = c("VoiceID", "Timepoint", "Antigen")) %>%
    select(all_of(meta_cols), starts_with("Norm_"))
  
  #writing the scores along with metadata to a file
  write_tsv(meta_data_withscores, "intermediate_files/well_data_withScores_normDMSO_step3.tsv")
  
  #returning
  return(meta_data_withscores)

}
