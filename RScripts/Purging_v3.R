# PURGING 2.0 ------
# Load required packages 
pacman::p_load(ggplot2, openxlsx, RColorBrewer, dplyr, grid, gridExtra, ggrepel,
               reshape2, tidyr, paletteer, ggtree, stringr, FSA, broom,purrr, 
               dunn.test)

# Load Italian Purging dataset
Ita_Purging <- read.table("Data/Purging/All_italian_freq.tsv", h=T)

# Calculate Lxy and Lyx ratios
Ita_Purging <- Ita_Purging %>%
  mutate(Lxy_HIGH = Sum_fxy_HIGH / Sum_fxy_INTERGENIC,
         Lyx_HIGH = Sum_fyx_HIGH / Sum_fyx_INTERGENIC,
         Lxy_MODERATE = Sum_fxy_MODERATE / Sum_fxy_INTERGENIC,
         Lyx_MODERATE = Sum_fyx_MODERATE / Sum_fyx_INTERGENIC,
         Lxy_LOW = Sum_fxy_LOW / Sum_fxy_INTERGENIC,
         Lyx_LOW = Sum_fyx_LOW / Sum_fyx_INTERGENIC,
         Lxy_MODIFIER = Sum_fxy_MODIFIER / Sum_fxy_INTERGENIC,
         Lyx_MODIFIER = Sum_fyx_MODIFIER / Sum_fyx_INTERGENIC,
         
         # Compute Rxy per mutation impact category
         Rxy_High = Lxy_HIGH / Lyx_HIGH,
         Rxy_Moderate = Lxy_MODERATE / Lyx_MODERATE,
         Rxy_Low = Lxy_LOW / Lyx_LOW,
         Rxy_Modifier = Lxy_MODIFIER / Lyx_MODIFIER)

# Select only relevant Rxy values
Rxy_Italian <- Ita_Purging %>%
  select(Rxy_High, Rxy_Moderate, Rxy_Low, Rxy_Modifier)

# Reshape data into long format for analysis
Rxy_Italian <- Rxy_Italian %>% 
  pivot_longer(cols = starts_with("Rxy_"),   # Select Rxy columns
               names_to = "Impact",          # Create column for categories
               values_to = "Rxy") %>%        # Assign values to Rxy
  mutate(Impact = str_remove(Impact, "Rxy_")) # Clean category names

# Ensure Impact is a factor before statistical tests
Rxy_Italian$Impact <- as.factor(Rxy_Italian$Impact)

# Kruskal-Wallis test
kruskal.test(Rxy ~ Impact, data = Rxy_Italian) # p-value < 2.2e-16
# Dunn’s test for pairwise comparisons
dunn.test(Rxy_Italian$Rxy, Rxy_Italian$Impact, method = "bonferroni")


# **Visualization - Boxplot**
ggplot(Rxy_Italian, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot() +
  labs(title = "Rxy by Impact Category - Italian",
       x = "Rxy",
       y = "Impact Category") +
  theme_minimal() +
  scale_x_continuous(limits = c(0.96, 1.04))

# **Visualization - Density Plot**
ggplot(Rxy_Italian, aes(x = Rxy, fill = Impact, color = Impact)) +
  geom_density(alpha = 0.5) +  # Overlay density curves
  labs(title = "Distribution of Rxy by Impact Category - Italian",
       x = "Rxy", y = "Density") +
  theme_minimal() +
  scale_x_continuous(limits = c(0.96, 1.04))

# All French ---- 
Fra_Purging <- read.table("Data/Purging/All_French_freq.tsv", h=T)

# Get Lxy and Lyx 
Fra_Purging <- Fra_Purging %>%
  mutate(Lxy_HIGH = Sum_fxy_HIGH/Sum_fxy_INTERGENIC,
         Lyx_HIGH = Sum_fyx_HIGH/Sum_fyx_INTERGENIC,
         Lxy_MODERATE = Sum_fxy_MODERATE/Sum_fxy_INTERGENIC,
         Lyx_MODERATE = Sum_fyx_MODERATE/Sum_fyx_INTERGENIC,
         Lxy_LOW = Sum_fxy_LOW/Sum_fxy_INTERGENIC,
         Lyx_LOW = Sum_fyx_LOW/Sum_fyx_INTERGENIC,
         Lxy_MODIFIER = Sum_fxy_MODIFIER/Sum_fxy_INTERGENIC,
         Lyx_MODIFIER = Sum_fyx_MODIFIER/Sum_fyx_INTERGENIC,
         
         # Get Rxy per impact 
         Rxy_High = Lxy_HIGH/Lyx_HIGH,
         Rxy_Moderate = Lxy_MODERATE/Lyx_MODERATE,
         Rxy_Low = Lxy_LOW/Lyx_LOW,
         Rxy_Modifier = Lxy_MODIFIER/Lyx_MODIFIER)

# Filter just Rxy values
Rxy_French<- Fra_Purging %>%
  select(Rxy_High,
         Rxy_Moderate,
         Rxy_Low,
         Rxy_Modifier)

# Transform data 
Rxy_French<- Rxy_French %>% 
  pivot_longer(cols = starts_with("Rxy_"),  # Select all R_ columns
               names_to = "Impact",       # New column for impact categories
               values_to = "Rxy") %>%      # New column for Rxy values
  mutate(Impact = str_remove(Impact, "Rxy_")) # Clean up impact names

# Plots 
ggplot(Rxy_French, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot() +
  labs(title = "Rxy by Impact Category - French",
       x = "Rxy",
       y = "Impact Category") +
  theme_minimal() +
  scale_x_continuous(limits = c(0.96, 1.04))

ggplot(Rxy_French, aes(x = Rxy, fill = Impact, color = Impact)) +
  geom_density(alpha = 0.5) +  # Overlay density curves
  labs( title = "Distribution of Rxy by Impact Category - French",
        x = "Rxy", y = "Density") +
  theme_minimal() +
  scale_x_continuous(limits = c(0.96, 1.04))

# Ensure Impact is a factor before statistical tests
Rxy_French$Impact <- as.factor(Rxy_French$Impact)

# Kruskal-Wallis test
kruskal.test(Rxy ~ Impact, data = Rxy_French) # p-value < 2.2e-16
# Dunn’s test for pairwise comparisons
dunn.test(Rxy_French$Rxy, Rxy_French$Impact, method = "bonferroni")


# All Native ---- 
Nat_Purging <- read.table("Data/Purging/All_natives_freq.tsv", h=T)

# Get Lxy and Lyx 
Nat_Purging <- Nat_Purging %>%
  mutate(Lxy_HIGH = Sum_fxy_HIGH/Sum_fxy_INTERGENIC,
         Lyx_HIGH = Sum_fyx_HIGH/Sum_fyx_INTERGENIC,
         Lxy_MODERATE = Sum_fxy_MODERATE/Sum_fxy_INTERGENIC,
         Lyx_MODERATE = Sum_fyx_MODERATE/Sum_fyx_INTERGENIC,
         Lxy_LOW = Sum_fxy_LOW/Sum_fxy_INTERGENIC,
         Lyx_LOW = Sum_fyx_LOW/Sum_fyx_INTERGENIC,
         Lxy_MODIFIER = Sum_fxy_MODIFIER/Sum_fxy_INTERGENIC,
         Lyx_MODIFIER = Sum_fyx_MODIFIER/Sum_fyx_INTERGENIC,
         
         # Get Rxy per impact 
         Rxy_High = Lxy_HIGH/Lyx_HIGH,
         Rxy_Moderate = Lxy_MODERATE/Lyx_MODERATE,
         Rxy_Low = Lxy_LOW/Lyx_LOW,
         Rxy_Modifier = Lxy_MODIFIER/Lyx_MODIFIER)

# Filter just Rxy values
Rxy_Natives<- Nat_Purging %>%
  select(Rxy_High,
         Rxy_Moderate,
         Rxy_Low,
         Rxy_Modifier)

# Transform data 
Rxy_Natives <- Rxy_Natives %>% 
  pivot_longer(cols = starts_with("Rxy_"),  # Select all R_ columns
               names_to = "Impact",       # New column for impact categories
               values_to = "Rxy") %>%      # New column for Rxy values
  mutate(Impact = str_remove(Impact, "Rxy_")) # Clean up impact names


# Plots 
ggplot(Rxy_Natives, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot() +
  labs(title = "Rxy by Impact Category - Natives",
       x = "Rxy",
       y = "Impact Category") +
  theme_minimal() +
  scale_x_continuous(limits = c(0.85, 1.15))
  
ggplot(Rxy_Natives, aes(x = Rxy, fill = Impact, color = Impact)) +
  geom_density(alpha = 0.5) +  # Overlay density curves
  labs( title = "Distribution of Rxy by Impact Category - Natives",
        x = "Rxy", y = "Density") +
  theme_minimal() +
  scale_x_continuous(limits = c(0.85, 1.15))


# Each population ----

# Load all required packages (excluding purrr for compatibility)
pacman::p_load( ggplot2, openxlsx, RColorBrewer, dplyr, grid, gridExtra, ggrepel,
  reshape2, tidyr, paletteer, ggtree, stringr, lmerTest, lme4, FSA,
  PMCMRplus, broom)

# ======================================================================
# Rxy ANALYSIS PIPELINE WITH WELCH'S ANOVA AND GAMES-HOWELL POST-HOC TESTS
# ======================================================================

# Define populations
italian_pops <- c("BB", "DL", "NF", "SH", "VT", "WS")
french_pops <- c("BU", "WB", "WE")

# List input files
freq_files <- list.files("Data/Purging/", pattern = "^All_.*_freq\\.tsv$", full.names = TRUE)

# Function to calculate Rxy per file
calculate_Rxy <- function(file_path) {
  pop <- stringr::str_extract(basename(file_path), "(?<=All_).*(?=_freq.tsv)")
  df <- read.table(file_path, header = TRUE)
  
  df <- df %>%
    mutate(Population = pop,
           Lxy_HIGH = Sum_fxy_HIGH / Sum_fxy_INTERGENIC,
           Lyx_HIGH = Sum_fyx_HIGH / Sum_fyx_INTERGENIC,
           Lxy_MODERATE = Sum_fxy_MODERATE / Sum_fxy_INTERGENIC,
           Lyx_MODERATE = Sum_fyx_MODERATE / Sum_fyx_INTERGENIC,
           Lxy_LOW = Sum_fxy_LOW / Sum_fxy_INTERGENIC,
           Lyx_LOW = Sum_fyx_LOW / Sum_fyx_INTERGENIC,
           Lxy_MODIFIER = Sum_fxy_MODIFIER / Sum_fxy_INTERGENIC,
           Lyx_MODIFIER = Sum_fyx_MODIFIER / Sum_fyx_INTERGENIC,
           Rxy_High = Lxy_HIGH / Lyx_HIGH,
           Rxy_Moderate = Lxy_MODERATE / Lyx_MODERATE,
           Rxy_Low = Lxy_LOW / Lyx_LOW,
           Rxy_Modifier = Lxy_MODIFIER / Lyx_MODIFIER) %>%
    select(Population, starts_with("Rxy_")) %>%
    pivot_longer(cols = starts_with("Rxy_"),
                 names_to = "Impact", values_to = "Rxy") %>%
    mutate(Impact = stringr::str_remove(Impact, "Rxy_"))
  
  return(df)
}

# Process all files
Rxy_all <- bind_rows(lapply(freq_files, calculate_Rxy)) %>%
  mutate(Group = case_when(
    Population %in% italian_pops ~ "Italian",
    Population %in% french_pops ~ "French",
    TRUE ~ "Native"
  ))

# ======================================================================
# STATISTICAL ANALYSIS FUNCTIONS
# ======================================================================

run_welch_anova <- function(data) {
  welch_result <- oneway.test(Rxy ~ Impact, data = data, var.equal = FALSE)
  
  if (welch_result$p.value < 0.05) {
    gh_result <- PMCMRplus::gamesHowellTest(Rxy ~ Impact, data = data)
    gh_tidy <- as_tibble(gh_result$p.value, rownames = "Comparison")
  } else {
    gh_tidy <- tibble(Comparison = "No significant differences")
  }
  
  list(
    ANOVA = broom::tidy(welch_result),
    PostHoc = gh_tidy
  )
}

analyze_population <- function(pop_name, full_data) {
  pop_data <- full_data %>% 
    filter(Population == pop_name) %>%
    mutate(Impact = as.factor(Impact))
  
  if (n_distinct(pop_data$Impact) < 2) {
    return(tibble(
      Population = pop_name,
      Test = "Insufficient groups",
      F.value = NA,
      p.value = NA,
      PostHoc = list(tibble(Note = "Need ≥2 impact categories"))
    ))
  }
  
  results <- run_welch_anova(pop_data)
  
  tibble(
    Population = pop_name,
    Test = "Welch's ANOVA",
    F.value = results$ANOVA$statistic,
    p.value = results$ANOVA$p.value,
    PostHoc = list(results$PostHoc)
  )
}

# ======================================================================
# RUN ANALYSIS FOR ALL POPULATIONS
# ======================================================================

italian_results <- do.call(rbind, lapply(italian_pops, analyze_population, full_data = Rxy_all))
french_results <- do.call(rbind, lapply(french_pops, analyze_population, full_data = Rxy_all))

all_results <- bind_rows(
  italian_results,
  french_results
) %>%
  mutate(Significant = ifelse(p.value < 0.05, "YES", "NO"))

# ======================================================================
# VISUALIZATION AND OUTPUT
# ======================================================================

anova_summary <- all_results %>%
  select(Population, F.value, p.value, Significant) %>%
  arrange(Significant, p.value)
print(anova_summary)

significant_comparisons <- all_results %>%
  filter(Significant == "YES") %>%
  unnest(PostHoc) %>%
  filter(!str_detect(Comparison, "No significant"))

if (nrow(significant_comparisons) > 0) {
  print(significant_comparisons)
} else {
  message("No significant pairwise differences found in any population")
}

plot_rxy_distribution <- function(data, title) {
  ggplot(data, aes(x = Rxy, y = Impact, fill = Population)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    facet_wrap(~Population, nrow = 2) +
    theme_minimal() +
    scale_x_continuous(limits = c(0.9, 1.1)) +
    labs(title = title, x = "Rxy", y = NULL)
}

print(plot_rxy_distribution(
  filter(Rxy_all, Group == "Italian"),
  "Comparative Rxy Across Italian Populations"
))

print(plot_rxy_distribution(
  filter(Rxy_all, Group == "French"),
  "Comparative Rxy Across French Populations"
))

