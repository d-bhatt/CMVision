### this is the second step script in identifying signals from the well-specific data
### this script will import the pre-processed dataset (which contains the scores)
### and make heatmap and correlation plots

#writing a function that will take in the processed data from step 1 and make plots
step2_makeCorrelationHeatmap = function(data, image_base_path){

  # ## importing the well-specific variable dataset including the computed scores
  # data = read_tsv("well_data_withScores.tsv")
  
  #### Step 1 - making correlation plot ####
  
  #selecting columns and renaming them for better visualization
  cor_data = data %>%
    select(
      Score_AllEmphasizeWellOccupancy_Z,
      Score_AllEqualWeights_Z,
      Score_AllOptimizedWeights_Z,
      Score_OnlyOccupancyIntensity_Z,
      Score_OnlyOccupancy_Z,
      Count_SpotInWell,
      PercentOccupied,
      Mean_SpotInWell_AreaShape_Area,
      Mean_SpotInWell_Intensity_MeanIntensity_Grey,
      Mean_Well_Intensity_IntegratedIntensity_Grey
    ) %>%
    drop_na() %>%
    rename(
      `Emphasize Well Occupancy` = Score_AllEmphasizeWellOccupancy_Z,
      `Equal Weights` = Score_AllEqualWeights_Z,
      `Optimized Weights` = Score_AllOptimizedWeights_Z,
      `Well occupancy & intensity` = Score_OnlyOccupancyIntensity_Z,
      `Well occupancy` = Score_OnlyOccupancy_Z,
      `Spot per well` = Count_SpotInWell,
      `% Occupancy` = PercentOccupied,
      `Mean spot size` = Mean_SpotInWell_AreaShape_Area,
      `Mean spot intensity` = Mean_SpotInWell_Intensity_MeanIntensity_Grey,
      `Mean well intensity` = Mean_Well_Intensity_IntegratedIntensity_Grey
    )
  
  #computing the correlation matrix
  cor_matrix = cor(cor_data, method = "spearman")
  
  #computing the p-value matrix
  cor_matrix_pval = cor.mtest(cor_matrix)
  
  #plotting the correlation matrix
  cor_plot = ggcorrplot(cor_matrix[rownames(cor_matrix),], method = "circle", type = "lower", lab = TRUE, 
                        lab_size = 3, outline.color = "gray", hc.order = FALSE,
                        colors = c("blue", "white", "red"), legend.title = "ρ",
                        title = "Clustered Correlation Between Scores and ELISpot Features",
                        ggtheme = theme_minimal(), p.mat = cor_matrix_pval$p, insig = "blank")
  
  #adding colours to the y axis label
  cor_plot = cor_plot + theme(
    axis.text.y = element_text(color = c(rep("black", 5), rep("#00B0F0", 5)))
  )+ labs(subtitle = "Spearman correlation coefficient (ρ): Blue = negative, Red = positive. Blank = insignificant")
  
  #saving
  ggsave(plot = cor_plot, filename = "plots/correlation_plot.tiff", height = 8, width = 10, dpi = 300, device = "tiff")
  
  
  #### Step 2 - making heatmap -- based on the normalized score only (and not Z score) ####
  
  ### !! Part 1 - making two heatmaps - based on the scores and based on the original variables !! ###
  
  #defining score columns in the dataset
  score_columns = c(
    "Score_AllEmphasizeWellOccupancy_Norm",
    "Score_AllEqualWeights_Norm",
    "Score_AllOptimizedWeights_Norm",
    "Score_OnlyOccupancyIntensity_Norm",
    "Score_OnlyOccupancy_Norm",
    "Count_SpotInWell",
    "PercentOccupied",
    "Mean_SpotInWell_AreaShape_Area",
    "Mean_SpotInWell_Intensity_MeanIntensity_Grey",
    "Mean_Well_Intensity_IntegratedIntensity_Grey"
  )
  
  #subsetting the columns required for making the heatmap and centering the variables
  hm_data = data %>%
    select(all_of(score_columns)) %>%
    scale()
  
  #transposing the heatmap
  hm_matrix = t(hm_data)
  
  #splitting scores and features
  scores = hm_matrix[1:5, , drop = FALSE]
  features = hm_matrix[6:10, , drop = FALSE]
  
  #using Score_AllOptimizedWeights_Norm for ordering (row-wise)
  row_order = order(data$Score_AllOptimizedWeights_Norm)
  
  #plotting the heatmap - scores
  ht_scores = Heatmap(
    scores[, row_order],
    name = "Scores",
    col = colorRamp2(
      seq(0, 1, length.out = 5),
      c("#ECECEC", "#D3DCE4", "#AABBCF", "#7492B2", "#1E4F89")
    ),
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    cluster_rows = FALSE,
    show_heatmap_legend = FALSE,
    height = unit(4, "cm")
  )
  
  #plotting heatmap - features
  ht_features = Heatmap(
    features[, row_order],
    name = "Original Variables",
    col = colorRamp2(
      seq(0, 1, length.out = 5),
      c("#ECECEC", "#D3DCE4", "#AABBCF", "#7492B2", "#1E4F89")
    ),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_heatmap_legend = FALSE,
    height = unit(4, "cm")
  )
  
  ### !! Part 2 - making the image tiles that are representative of the scores and original variables !! ###

  #preparing image thumbnails - writing two interlinked lapply
  image_rows = lapply(score_columns, function(score_name) {
    score_vector = data[[score_name]]
    score_vector = ifelse(is.na(score_vector) | !is.finite(score_vector), 0, score_vector)
    quantile_positions = quantile(score_vector, probs = seq(0, 1, length.out = 40), na.rm = TRUE)
    indices = sapply(quantile_positions, function(q) which.min(abs(score_vector - q)))
    wells = data[indices, c("Image", "WellNum")]
    
    lapply(seq_len(nrow(wells)), function(i) {
      info = wells[i, ]
      subfolder = paste0("image", info$Image)
      filename = paste0("image", info$Image, "_Well_", info$WellNum, ".png")
      full_path = file.path(image_base_path, subfolder, filename)
      
      if (file.exists(full_path)) {
        img = image_read(full_path)
        img = image_scale(img, "x100")  # Resize height to 100px, preserving aspect ratio
        img = image_modulate(img, brightness = 120, saturation = 100, hue = 100)
        rasterGrob(as.raster(img), interpolate = TRUE)
      } else {
        nullGrob()
      }
      
    })
  })
  
  #preparing a fresh layout of the image thumbnails
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(length(image_rows) + 1, 1, heights = unit(c(rep(0.5, length(image_rows)), 5), "null"))))
  
  #going over all images and assigning positions
  for (row_idx in seq_along(image_rows)) {
    pushViewport(viewport(layout.pos.row = row_idx))
    row_images = image_rows[[row_idx]]
    for (i in seq_along(row_images)) {
      x_coord = unit((i - 0.5) / length(row_images), "npc")
      pushViewport(viewport(
        x = x_coord,
        y = unit(0.5, "npc"),
        width = unit(1 / length(row_images), "npc"),
        height = unit(1, "npc"),
        just = c("center", "center")
      ))
      grid.draw(row_images[[i]])
      popViewport()
    }
    popViewport()
  }
  
  #adding heatmaps to the layout 
  pushViewport(viewport(layout.pos.row = length(image_rows) + 1))
  draw(ht_scores %v% ht_features, newpage = FALSE)
  popViewport()

}
