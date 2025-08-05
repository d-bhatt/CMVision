### this is the fourth step script in identifying signals from the well-specific data
### this script will import the pre-processed normalized score dataset 
### and make CMV and Spike specific plots

#writing a function that will take in the processed data from step 3 and 
#make the CMV and Spike specific plots
step4_makeBoxplots = function(elispot_data){

  #filter data set to exclude individuals who were dropped
  elispot_data_filtered = elispot_data %>%
    filter(Exc_Inc == "I",
           Antigen %in% c("cmv", "spike"),
           !is.na(Cohort))
  
  #defining the score name and the label that will be used later for plotting
  score_cols = c("Norm_Score_AllOptimizedWeights_Norm")
  score_labels = c(
    Norm_Score_AllOptimizedWeights_Norm = "Optimized Weights")
  
  #### preparing the dataframes for plotting ####
  
  ## long pivoting to rename columns, check for NA scores, grouping on the three variables and 
  ## double checking for the two timepoints
  
  #long pivot to make the cmv dataframe that will be used later for plotting
  cmv_long = elispot_data_filtered %>%
    filter(Exc_Inc == "I") %>% 
    filter(Antigen == "cmv", Timepoint %in% c("T0", "T2")) %>%
    pivot_longer(cols = all_of(score_cols), names_to = "ScoreType", values_to = "Score") %>%
    filter(!is.na(Score)) %>%
    group_by(VoiceID, ScoreType, Cohort) %>%
    filter(n_distinct(Timepoint) == 2) %>%
    ungroup()
  
  #long pivot to make the spike dataframe that will be used later for plotting
  spike_long = elispot_data_filtered %>%
    filter(Exc_Inc == "I") %>% 
    filter(Antigen == "spike", Timepoint %in% c("T0", "T2")) %>%
    pivot_longer(cols = all_of(score_cols), names_to = "ScoreType", values_to = "Score") %>%
    filter(!is.na(Score)) %>%
    # Ensure both T0 and T2 are present for each VoiceID
    group_by(VoiceID, ScoreType, Cohort) %>%
    filter(n_distinct(Timepoint) == 2) %>%
    ungroup()
  
  #creating the dataframe for plotting CMV vs age at T2
  cmv_t2 = elispot_data_filtered %>%
    filter(Antigen == "cmv", Timepoint == "T2") %>%
    pivot_longer(cols = all_of(score_cols), names_to = "ScoreType", values_to = "Score")
  
  #creating the dataframe for plotting Spike vs age at T2
  spike_t2 = elispot_data_filtered %>%
    filter(Antigen == "spike", Timepoint == "T2") %>%
    pivot_longer(cols = all_of(score_cols), names_to = "ScoreType", values_to = "Score")
  
  #### plotting section ####
  
  #declaring a standard base theme to be used for plotting later on
  base_theme = theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  
  ## ! Plot 1 - CMV data at T2 - Response (Norm. score) vs Age
  cmv_age_plot =  ggplot(cmv_t2, aes(x = Age, y = Score, color = Cohort)) +
    geom_smooth(method = "lm", se = TRUE) +
    stat_fit_glance(
      method = "lm",
      aes(label = after_stat(paste0("R2: ", signif(r.squared, 3), "\nP: ", signif(p.value, 3)))),
      size = 3,
      color = "black"
    ) +
    facet_grid(. ~  Cohort , scales = "free", labeller = labeller(ScoreType = score_labels)) +
    scale_color_viridis_d(option = "C") +
    base_theme +
    labs(title = "CMV T2 Response vs Age", x = "Age", y = "Normalized Score")
  
  #saving the plot
  ggsave(filename = "plots/cmv_age_T2_plot.tiff", plot = cmv_age_plot, device = "tiff", width = 10, height = 3, dpi = 300)
  
  
  ## ! Plot 2 - Spike data at T2 - Response (Norm. score) vs Age
  spike_age_plot = ggplot(spike_t2, aes(x = Age, y = Score, color = Cohort)) +
    #geom_jitter(alpha = 1) +
    geom_smooth(method = "lm", se = TRUE) +
    stat_fit_glance(
      method = "lm",
      aes(label = after_stat(paste0("R2: ", signif(r.squared, 3), "\nP: ", signif(p.value, 3)))),
      size = 3,
      color = "black"
    ) +
    facet_grid(. ~ Cohort, scales = "free", labeller = labeller(ScoreType = score_labels)) +
    scale_color_viridis_d(option = "C") +
    base_theme +
    labs(title = "Spike T2 Response vs Age", x = "Age", y = "Normalized Score")
  
  #saving the plot
  ggsave(filename = "plots/spike_age_T2_plot.tiff", plot = spike_age_plot, device = "tiff", width = 10, height = 3, dpi = 300)
  
  
  ## ! Plot 3 - Spike Δ(T2−T0) vs CMV T2
  
  #preparing Spike Δ
  spike_delta = elispot_data_filtered %>%
    filter(Antigen == "spike", Timepoint %in% c("T0", "T2")) %>%
    pivot_longer(cols = all_of(score_cols), names_to = "ScoreType", values_to = "Score") %>%
    group_by(VoiceID, Timepoint, ScoreType, Cohort) %>%
    summarize(Score = mean(Score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Timepoint, values_from = Score) %>%
    filter(!is.na(T0), !is.na(T2)) %>%
    mutate(Spike_Delta = T2 - T0)
  
  #preparing CMV T2
  cmv_t2_avg = cmv_t2 %>%
    group_by(VoiceID, ScoreType, Cohort) %>%
    summarize(CMV_Score = mean(Score, na.rm = TRUE), .groups = "drop")
  
  #preparing dataframe for the correlation
  correlation_spike_cmv_t2 = inner_join(spike_delta, cmv_t2_avg, by = c("VoiceID", "ScoreType", "Cohort"))
  
  #plotting the correlation of Spike Δ and CMV T2
  spike_delta_cmv_plot = ggplot(correlation_spike_cmv_t2, aes(x = CMV_Score, y = Spike_Delta, color = Cohort)) +
    #geom_point(alpha = 0.6) +
    geom_smooth(method = "lm") +
    stat_fit_glance(
      method = "lm",
      aes(label = after_stat(paste0("R2: ", signif(r.squared, 3), "\nP: ", signif(p.value, 3)))),
      size = 3,
      color = "black"
    ) +
    facet_grid(. ~  Cohort, scales = "free", labeller = labeller(ScoreType = score_labels)) +
    scale_color_viridis_d(option = "C") +
    base_theme +
    labs(title = "Spike Δ(T2−T0) vs CMV T2 Score", x = "CMV Score (T2)", y = "Spike Score Change (T2 − T0)")
  
  #saving the plot
  ggsave(filename = "plots/spike_delta_CMV_T2_plot.tiff", plot = spike_delta_cmv_plot, device = "tiff", width = 10, height = 3, dpi = 300)
  
  
  ## ! Plot 4 - Spike paired/hybrid plot
  spike_hybrid_plot = ggplot(spike_long, aes(x = Timepoint, y = Score, color = Cohort)) +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.1, size = 1.5) +
    geom_line(aes(group = VoiceID), color = "gray60", alpha = 0.4) +
    geom_boxplot(aes(group = Timepoint), width = 0.3, alpha = 1, position = position_dodge(width = 0.6)) +
    facet_wrap(. ~ Cohort, ncol=4, nrow=5, labeller = labeller(ScoreType = score_labels)) +
    scale_color_viridis_d(option = "C") +
    scale_y_log10(
      limits = c(1e-10, 10),
      breaks = c(1e-10, 1e-5,1)
    )+
    base_theme +
    labs(title = "Spike Response: Paired and Distribution View", y = "Normalized Score", x = "Timepoint")
  
  #adding statistics to the plot
  spike_hybrid_plot = spike_hybrid_plot+
    # Statistical comparison: Wilcoxon test with p-value
    stat_compare_means(
      method = "wilcox.test",
      aes(group = Timepoint),
      #paired = TRUE,
      label = "p.format",
      label.x = 1.3,
      label.y = 0.95 * max(spike_long$Score, na.rm = TRUE),
      size = 3
    )
  
  #saving the plot
  ggsave(filename = "plots/spike_hybrid_plot.tiff", plot = spike_hybrid_plot, device = "tiff", width = 10, height = 3, dpi = 300)
  
  
  ## ! Plot 5 - CMV hybrid plot
  cmv_hybrid_plot = ggplot(cmv_long, aes(x = Timepoint, y = Score, color = Cohort)) +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.1, size = 1.5) +
    geom_line(aes(group = VoiceID), color = "gray60", alpha = 0.4) +
    geom_boxplot(aes(group = Timepoint), width = 0.3, alpha = 0.3, position = position_dodge(width = 0.6)) +
    facet_wrap(. ~ Cohort, scales = "free", ncol=4, nrow=5, labeller = labeller(ScoreType = score_labels)) +
    scale_color_viridis_d(option = "C") +
    scale_y_log10(
      limits = c(1e-10, 10),
      breaks = c(1e-10, 1e-5, 1)
    )+
    base_theme +
    labs(title = "CMV Response: Paired and Distribution View", y = "Normalized Score", x = "Timepoint")
  
  #adding statistics to the plot
  cmv_hybrid_plot = cmv_hybrid_plot+
    # Statistical comparison: Wilcoxon test with p-value
    stat_compare_means(
      method = "wilcox.test",
      aes(group = Timepoint),
      #paired = TRUE,
      label = "p.format",
      label.x = 1.3,
      label.y = 0.1 * max(cmv_long$Score, na.rm = TRUE),
      size = 3
    )
  
  #saving the plot
  ggsave(filename = "plots/cmv_hybrid_plot.tiff", plot = cmv_hybrid_plot, device = "tiff", width = 10, height = 3, dpi = 300)

}