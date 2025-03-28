# Purging per population  ------
# Load important packages 
pacman::p_load(ggplot2,openxlsx,RColorBrewer,dplyr,grid,gridExtra,ggrepel,
               reshape2,tidyr,paletteer,ggtree,stringr,lmerTest,lme4)

# Create a function for Lxy ----
Lxy <- function(freq_category_X, freq_category_Y, freq_intergenic_X, freq_intergenic_Y) {
  numerator = freq_category_X * (1 - freq_category_Y)
  denominator = freq_intergenic_X * (1 - freq_intergenic_Y)
  numerator / denominator
}


# Populations to analyze
populations <- c("BB", "DL", "NF", "WS", "VT", "SH")

# Create containers
all_rxy <- list()
all_plots <- list()

# Analysis loop
for(pop in populations) {
  cat("\n=== Processing population:", pop, "===\n")
  
  # Load data
  purging_data <- read.table(paste0("Data/Purging/ALL_", pop, ".tsv"), header = TRUE)
  
  # Calculate frequencies (display first 3 rows)
  purging_data <- purging_data %>%
    mutate(
      High_Freq = Alternate_HIGH/Total_HIGH,
      Moderate_Freq = Alternate_MODERATE/Total_MODERATE,
      Low_Freq = Alternate_LOW/Total_LOW,
      Modifier_Freq = Alternate_MODIFIER/Total_MODIFIER,
      Intergenic_Freq = Alternate_INTERGENIC/Total_INTERGENIC
    )
  
  cat("\nFrequency calculations:\n")
  print(head(purging_data, 3))
  
  # Reshape data
  purging_wide <- purging_data %>%
    select(Jackknife, Origin, 
           High_Freq, Moderate_Freq, Low_Freq, 
           Modifier_Freq, Intergenic_Freq) %>%
    pivot_wider(
      names_from = Origin,
      values_from = c(High_Freq, Moderate_Freq,
                      Low_Freq, Modifier_Freq, Intergenic_Freq),
      names_glue = "{.value}_{Origin}"
    )
  
  # Calculate Lxy and Rxy
  purging_wide <- purging_wide %>%
    mutate(
      L_High_Int_Nat = Lxy(High_Freq_introduced, High_Freq_native,
                           Intergenic_Freq_introduced, Intergenic_Freq_native),
      L_Moderate_Int_Nat = Lxy(Moderate_Freq_introduced, Moderate_Freq_native,
                               Intergenic_Freq_introduced, Intergenic_Freq_native),
      L_Low_Int_Nat = Lxy(Low_Freq_introduced, Low_Freq_native,
                          Intergenic_Freq_introduced, Intergenic_Freq_native),
      L_Modifier_Int_Nat = Lxy(Modifier_Freq_introduced, Modifier_Freq_native,
                               Intergenic_Freq_introduced, Intergenic_Freq_native),
      
      L_High_Nat_Int = Lxy(High_Freq_native, High_Freq_introduced,
                           Intergenic_Freq_native, Intergenic_Freq_introduced),
      L_Moderate_Nat_Int = Lxy(Moderate_Freq_native, Moderate_Freq_introduced,
                               Intergenic_Freq_native, Intergenic_Freq_introduced),
      L_Low_Nat_Int = Lxy(Low_Freq_native, Low_Freq_introduced,
                          Intergenic_Freq_native, Intergenic_Freq_introduced),
      L_Modifier_Nat_Int = Lxy(Modifier_Freq_native, Modifier_Freq_introduced,
                               Intergenic_Freq_native, Intergenic_Freq_introduced),
      
      Rxy_High = L_High_Int_Nat/L_High_Nat_Int,
      Rxy_Moderate = L_Moderate_Int_Nat/L_Moderate_Nat_Int,
      Rxy_Low = L_Low_Int_Nat/L_Low_Nat_Int,
      Rxy_Modifier = L_Modifier_Int_Nat/L_Modifier_Nat_Int
    )
  
  # Create and display Rxy dataframe
  rxy_df <- purging_wide %>%
    select(Rxy_High, Rxy_Moderate, Rxy_Low, Rxy_Modifier) %>%
    pivot_longer(
      cols = starts_with("Rxy_"),
      names_to = "Impact",
      values_to = "Rxy"
    ) %>%
    mutate(
      Impact = str_remove(Impact, "Rxy_"),
      Population = pop
    )
  
  cat("\nRxy summary for", pop, ":\n")
  print(summary(rxy_df$Rxy))
  
  # Store for combined analysis
  all_rxy[[pop]] <- rxy_df
  
  # Create and display plots
  box_plot <- ggplot(rxy_df, aes(x = Rxy, y = Impact, fill = Impact)) +
    geom_boxplot() +
    labs(title = paste("Rxy for", pop),
         x = "Rxy", y = "") +
    theme_minimal() +
    scale_x_continuous(limits = c(0.94, 1.06))
  
  density_plot <- ggplot(rxy_df, aes(x = Rxy, fill = Impact, color = Impact)) +
    geom_density(alpha = 0.5) +
    labs(title = paste("Rxy Distribution for", pop),
         x = "Rxy", y = "Density") +
    theme_minimal()
  
  # Display plots side by side
  grid.arrange(box_plot, density_plot, ncol = 2)
  
  # Pause between populations
  readline(prompt = "Press [Enter] to continue to next population...")
}

# Combined analysis
combined_rxy <- bind_rows(all_rxy)

# Display combined results
cat("\n=== Combined Results ===\n")
print(head(combined_rxy, 6))

# Combined visualization
combined_box <- ggplot(combined_rxy, aes(x = Rxy, y = Impact, fill = Population)) +
  geom_boxplot() +
  facet_wrap(~Population) +
  labs(title = "Comparative Rxy Across Populations",
       x = "Rxy", y = "") +
  theme_minimal()

combined_density <- ggplot(combined_rxy, aes(x = Rxy, fill = Population)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Impact, scales = "free_y") +
  labs(title = "Comparative Rxy Distributions",
       x = "Rxy", y = "Density") +
  theme_minimal()

# Display combined plots
print(combined_box)
print(combined_density)
