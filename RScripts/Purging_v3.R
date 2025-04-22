# PURGING 2.0 ------
# Load important packages 
pacman::p_load(ggplot2,openxlsx,RColorBrewer,dplyr,grid,gridExtra,ggrepel,
               reshape2,tidyr,paletteer,ggtree,stringr,lmerTest,lme4,FSA)

# All Italian ---- 
Ita_Purging <- read.table("Data/Purging/All_italian_freq.tsv", h=T)

# Get Lxy and Lyx 
Ita_Purging <- Ita_Purging %>%
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
Rxy_Italian<- Ita_Purging %>%
  select(Rxy_High,
         Rxy_Moderate,
         Rxy_Low,
         Rxy_Modifier)

# Transform data 
Rxy_Italian <- Rxy_Italian %>% 
  pivot_longer(cols = starts_with("Rxy_"),  # Select all R_ columns
               names_to = "Impact",       # New column for impact categories
               values_to = "Rxy") %>%      # New column for Rxy values
  mutate(Impact = str_remove(Impact, "Rxy_")) # Clean up impact names

# Plots 
ggplot(Rxy_Italian, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot() +
  labs(title = "Rxy by Impact Category - Italian",
       x = "Rxy",
       y = "Impact Category") +
  theme_minimal() +
  scale_x_continuous(limits = c(0.96, 1.04))

ggplot(Rxy_Italian, aes(x = Rxy, fill = Impact, color = Impact)) +
  geom_density(alpha = 0.5) +  # Overlay density curves
  labs( title = "Distribution of Rxy by Impact Category - Italian",
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

# Define population groups
italian_pops <- c("BB", "DL", "NF", "SH", "VT", "WS")
french_pops <- c("BU", "WB", "WE")

# List input files
freq_files <- list.files("Data/Purging/", pattern = "^All_.*_freq\\.tsv$", full.names = TRUE)

# Function to calculate Rxy per file
calculate_Rxy <- function(file_path) {
  # Extract population name from file
  pop <- str_extract(basename(file_path), "(?<=All_).*(?=_freq.tsv)")
  
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
    mutate(Impact = str_remove(Impact, "Rxy_"))
  
  return(df)
}

# Apply function to all freq files
Rxy_all <- bind_rows(lapply(freq_files, calculate_Rxy))

# Assign population group
Rxy_all <- Rxy_all %>%
  mutate(Group = case_when(
    Population %in% italian_pops ~ "Italian",
    Population %in% french_pops ~ "French",
    TRUE ~ "Native"
  ))

# Optional: reorder population factor for better plot layout
Rxy_all$Population <- factor(Rxy_all$Population,
                             levels = c(italian_pops, french_pops))

# -----------------------
# PLOT FOR ITALIAN POPS
# -----------------------
Rxy_italian_only <- Rxy_all %>% filter(Group == "Italian")

ggplot(Rxy_italian_only, aes(x = Rxy, y = Impact, fill = Population)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  facet_wrap(~Population, nrow = 2) +
  theme_minimal() +
  scale_x_continuous(limits = c(0.9, 1.1)) +
  labs(title = "Comparative Rxy Across Italian Populations",
       x = "Rxy", y = NULL)

# -----------------------
#  PLOT FOR FRENCH POPS
# -----------------------
Rxy_french_only <- Rxy_all %>% filter(Group == "French")

ggplot(Rxy_french_only, aes(x = Rxy, y = Impact, fill = Population)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  facet_wrap(~Population, nrow = 1) +
  theme_minimal() +
  scale_x_continuous(limits = c(0.9, 1.1)) +
  labs(title = "Comparative Rxy Across French Populations",
       x = "Rxy", y = NULL)

