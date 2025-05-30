# PURGING 2.0 ------
# Load required packages 
pacman::p_load(ggplot2, openxlsx, RColorBrewer, dplyr, grid, gridExtra, ggrepel,
               reshape2, tidyr, paletteer, ggtree, stringr, lmerTest, lme4, FSA,
               PMCMRplus, broom)

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

# Get summary of the results after jackknifing 
calculate_minmax <- function(data) { 
  data %>%
    group_by(Impact) %>% 
    summarise(
      mean_Rxy = mean(Rxy, na.rm = TRUE), #Retrieve mean value 
      min_Rxy = min(Rxy, na.rm = TRUE), #Retrieve minimum value
      max_Rxy = max(Rxy, na.rm = TRUE), #Retrieve maximum value 
      .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), 
                  ~ formatC(., format = "f", digits = 3)))}

calculate_minmax(Rxy_Italian)


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

calculate_minmax(Rxy_French)


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
    TRUE ~ "Native" ))

# ======================================================================
# STATISTICAL ANALYSIS FUNCTIONS
# ======================================================================
population_results <- lapply(unique(Rxy_all$Population), function(pop) {
  Rxy_all %>%
    filter(Population == pop) %>%
    calculate_minmax() %>%
    mutate(Population = pop) %>%
    select(Population, everything())
})

# Name results properly
names(population_results) <- unique(Rxy_all$Population)

# Print formatted results
for (pop in names(population_results)) {
  cat("\n=== Results for", pop, "===\n")
  print(as_tibble(population_results[[pop]]))
}


# ======================================================================
# VISUALIZATION FOR EACH POPULATION
# ======================================================================

# Function to create a boxplot panel by population
plot_rxy_distribution <- function(data, title) {
  ggplot(data, aes(x = Rxy, y = Impact, fill = Impact)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    facet_wrap(~ Population, nrow = 2) +
    theme_minimal() +
    scale_x_continuous(limits = c(0.9, 1.1)) +
    labs(title = title, x = "Rxy", y = "Impact Category")
}

# Create panels for Italian and French populations
boxplot_Italian <- plot_rxy_distribution(filter(Rxy_all, Group == "Italian"),
                                         "Comparative Rxy Across Italian Populations")
boxplot_French  <- plot_rxy_distribution(filter(Rxy_all, Group == "French"),
                                         "Comparative Rxy Across French Populations")

# Print the boxplot panels
print(boxplot_Italian)
print(boxplot_French)

# Function to create a density plot panel by population
plot_rxy_density <- function(data, title) {
  ggplot(data, aes(x = Rxy, fill = Population, color = Population)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ Population, nrow = 2) +
    theme_minimal() +
    scale_x_continuous(limits = c(0.9, 1.1)) +
    labs(title = title, x = "Rxy", y = "Density")
}

# Create panels for Italian and French populations
density_Italian <- plot_rxy_density(filter(Rxy_all, Group == "Italian"),
                                    "Distribution of Rxy Across Italian Populations")
density_French  <- plot_rxy_density(filter(Rxy_all, Group == "French"),
                                    "Distribution of Rxy Across French Populations")

# Print the density plot panels
print(density_Italian)
print(density_French)

