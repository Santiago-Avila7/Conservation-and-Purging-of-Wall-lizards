# Corrected Purging  ------
# Load important packages 
pacman::p_load(ggplot2,openxlsx,RColorBrewer,dplyr,grid,gridExtra,ggrepel,
               reshape2,tidyr,paletteer,ggtree,stringr,lmerTest,lme4)

# Create a function for Lxy ----
Lxy <- function(freq_category_X, freq_category_Y, freq_intergenic_X, freq_intergenic_Y) {
  numerator = freq_category_X * (1 - freq_category_Y)
  denominator = freq_intergenic_X * (1 - freq_intergenic_Y)
  numerator / denominator
}

# Italian data corrected by sample size - All origins annotation ----
Purging_data <- read.table("Data/Purging/All_runs_5th.tsv", h=T)

# Calculate allele frequency in each Jackknife 
Purging_data<- Purging_data %>%
  mutate(High_Freq = Alternate_HIGH/Total_HIGH,
         Moderate_Freq = Alternate_MODERATE/Total_MODERATE,
         Low_Freq=Alternate_LOW/Total_LOW,
         Modifier_Freq=Alternate_MODIFIER/Total_MODIFIER,
         Intergenic_Freq=Alternate_INTERGENIC/Total_INTERGENIC)

# Select and make long for Rxy calculation
Purging_wide <- Purging_data %>%
  select(Jackknife, Origin, 
         High_Freq, Moderate_Freq, Low_Freq, 
         Modifier_Freq, Intergenic_Freq) %>%
  pivot_wider(names_from = Origin,
              values_from = c(High_Freq, Moderate_Freq,
                              Low_Freq, Modifier_Freq, Intergenic_Freq),
              names_glue = "{.value}_{Origin}")

# Calculate Lxy and Lyx per category 
Purging_wide<- Purging_wide %>%
  mutate(
    # For introduced not native 
    L_High_Int_Nat = Lxy (High_Freq_introduced,High_Freq_native,
                          Intergenic_Freq_introduced,Intergenic_Freq_native),
    L_Moderate_Int_Nat = Lxy (Moderate_Freq_introduced,Moderate_Freq_native,
                              Intergenic_Freq_introduced,Intergenic_Freq_native),
    L_Low_Int_Nat = Lxy (Low_Freq_introduced,Low_Freq_native,
                         Intergenic_Freq_introduced,Intergenic_Freq_native),
    L_Modifier_Int_Nat = Lxy (Modifier_Freq_introduced,Modifier_Freq_native,
                              Intergenic_Freq_introduced,Intergenic_Freq_native),
    
    # For native not introduced 
    L_High_Nat_Int= Lxy(High_Freq_native,High_Freq_introduced,
                        Intergenic_Freq_native,Intergenic_Freq_introduced),
    L_Moderate_Nat_Int= Lxy(Moderate_Freq_native,Moderate_Freq_introduced,
                            Intergenic_Freq_native,Intergenic_Freq_introduced),
    L_Low_Nat_Int= Lxy(Low_Freq_native,Low_Freq_introduced,
                       Intergenic_Freq_native,Intergenic_Freq_introduced),
    L_Modifier_Nat_Int= Lxy(Modifier_Freq_native,Modifier_Freq_introduced,
                            Intergenic_Freq_native,Intergenic_Freq_introduced),
    #For Rxy
    Rxy_High = L_High_Int_Nat/L_High_Nat_Int,
    Rxy_Moderate = L_Moderate_Int_Nat/L_Moderate_Nat_Int,
    Rxy_Low = L_Low_Int_Nat/L_Low_Nat_Int,
    Rxy_Modifier = L_Modifier_Int_Nat/L_Modifier_Nat_Int)

# Select data 
Rxy_Italian<- Purging_wide %>%
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
  labs(title = "Rxy by Impact Category",
       x = "Rxy",
       y = "Impact Category") +
  theme_minimal()

ggplot(Rxy_Italian, aes(x = Rxy, fill = Impact, color = Impact)) +
  geom_density(alpha = 0.5) +  # Overlay density curves
  labs( title = "Distribution of Rxy by Impact Category",
        x = "Rxy", y = "Density") +
  theme_minimal()         


# Italian data all samples - All origins annotation ----

Purging_data <- read.table("Data/Purging/All_runs_6th.tsv", h=T)

# Calculate allele frequency in each Jackknife 
Purging_data<- Purging_data %>%
  mutate(High_Freq = Alternate_HIGH/Total_HIGH,
         Moderate_Freq = Alternate_MODERATE/Total_MODERATE,
         Low_Freq=Alternate_LOW/Total_LOW,
         Modifier_Freq=Alternate_MODIFIER/Total_MODIFIER,
         Intergenic_Freq=Alternate_INTERGENIC/Total_INTERGENIC)

# Select and make long for Rxy calculation
Purging_wide <- Purging_data %>%
  select(Jackknife, Origin, 
         High_Freq, Moderate_Freq, Low_Freq, 
         Modifier_Freq, Intergenic_Freq) %>%
  pivot_wider(names_from = Origin,
              values_from = c(High_Freq, Moderate_Freq,
                              Low_Freq, Modifier_Freq, Intergenic_Freq),
              names_glue = "{.value}_{Origin}")

# Calculate Lxy and Lyx per category 
Purging_wide<- Purging_wide %>%
  mutate(
    # For introduced not native 
    L_High_Int_Nat = Lxy (High_Freq_introduced,High_Freq_native,
                          Intergenic_Freq_introduced,Intergenic_Freq_native),
    L_Moderate_Int_Nat = Lxy (Moderate_Freq_introduced,Moderate_Freq_native,
                              Intergenic_Freq_introduced,Intergenic_Freq_native),
    L_Low_Int_Nat = Lxy (Low_Freq_introduced,Low_Freq_native,
                         Intergenic_Freq_introduced,Intergenic_Freq_native),
    L_Modifier_Int_Nat = Lxy (Modifier_Freq_introduced,Modifier_Freq_native,
                              Intergenic_Freq_introduced,Intergenic_Freq_native),
    
    # For native not introduced 
    L_High_Nat_Int= Lxy(High_Freq_native,High_Freq_introduced,
                        Intergenic_Freq_native,Intergenic_Freq_introduced),
    L_Moderate_Nat_Int= Lxy(Moderate_Freq_native,Moderate_Freq_introduced,
                            Intergenic_Freq_native,Intergenic_Freq_introduced),
    L_Low_Nat_Int= Lxy(Low_Freq_native,Low_Freq_introduced,
                       Intergenic_Freq_native,Intergenic_Freq_introduced),
    L_Modifier_Nat_Int= Lxy(Modifier_Freq_native,Modifier_Freq_introduced,
                            Intergenic_Freq_native,Intergenic_Freq_introduced),
    #For Rxy
    Rxy_High = L_High_Int_Nat/L_High_Nat_Int,
    Rxy_Moderate = L_Moderate_Int_Nat/L_Moderate_Nat_Int,
    Rxy_Low = L_Low_Int_Nat/L_Low_Nat_Int,
    Rxy_Modifier = L_Modifier_Int_Nat/L_Modifier_Nat_Int)

# Select data 
Rxy_Italian<- Purging_wide %>%
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
  labs(title = "Rxy by Impact Category",
       x = "Rxy",
       y = "Impact Category") +
  theme_minimal()

ggplot(Rxy_Italian, aes(x = Rxy, fill = Impact, color = Impact)) +
  geom_density(alpha = 0.5) +  # Overlay density curves
  labs( title = "Distribution of Rxy by Impact Category",
        x = "Rxy", y = "Density") +
  theme_minimal()         


# Italian data corrected by sample size - Only Italian annotation
Purging_data <- read.table("Data/Purging/All_runs_7th.tsv", h=T)

# Calculate allele frequency in each Jackknife 
Purging_data<- Purging_data %>%
  mutate(High_Freq = Alternate_HIGH/Total_HIGH,
         Moderate_Freq = Alternate_MODERATE/Total_MODERATE,
         Low_Freq=Alternate_LOW/Total_LOW,
         Modifier_Freq=Alternate_MODIFIER/Total_MODIFIER,
         Intergenic_Freq=Alternate_INTERGENIC/Total_INTERGENIC)

# Select and make long for Rxy calculation
Purging_wide <- Purging_data %>%
  select(Jackknife, Origin, 
         High_Freq, Moderate_Freq, Low_Freq, 
         Modifier_Freq, Intergenic_Freq) %>%
  pivot_wider(names_from = Origin,
              values_from = c(High_Freq, Moderate_Freq,
                              Low_Freq, Modifier_Freq, Intergenic_Freq),
              names_glue = "{.value}_{Origin}")

# Calculate Lxy and Lyx per category 
Purging_wide<- Purging_wide %>%
  mutate(
    # For introduced not native 
    L_High_Int_Nat = Lxy (High_Freq_introduced,High_Freq_native,
                          Intergenic_Freq_introduced,Intergenic_Freq_native),
    L_Moderate_Int_Nat = Lxy (Moderate_Freq_introduced,Moderate_Freq_native,
                              Intergenic_Freq_introduced,Intergenic_Freq_native),
    L_Low_Int_Nat = Lxy (Low_Freq_introduced,Low_Freq_native,
                         Intergenic_Freq_introduced,Intergenic_Freq_native),
    L_Modifier_Int_Nat = Lxy (Modifier_Freq_introduced,Modifier_Freq_native,
                              Intergenic_Freq_introduced,Intergenic_Freq_native),
    
    # For native not introduced 
    L_High_Nat_Int= Lxy(High_Freq_native,High_Freq_introduced,
                        Intergenic_Freq_native,Intergenic_Freq_introduced),
    L_Moderate_Nat_Int= Lxy(Moderate_Freq_native,Moderate_Freq_introduced,
                            Intergenic_Freq_native,Intergenic_Freq_introduced),
    L_Low_Nat_Int= Lxy(Low_Freq_native,Low_Freq_introduced,
                       Intergenic_Freq_native,Intergenic_Freq_introduced),
    L_Modifier_Nat_Int= Lxy(Modifier_Freq_native,Modifier_Freq_introduced,
                            Intergenic_Freq_native,Intergenic_Freq_introduced),
    #For Rxy
    Rxy_High = L_High_Int_Nat/L_High_Nat_Int,
    Rxy_Moderate = L_Moderate_Int_Nat/L_Moderate_Nat_Int,
    Rxy_Low = L_Low_Int_Nat/L_Low_Nat_Int,
    Rxy_Modifier = L_Modifier_Int_Nat/L_Modifier_Nat_Int)

# Select data 
Rxy_Italian<- Purging_wide %>%
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
  labs(title = "Rxy by Impact Category",
       x = "Rxy",
       y = "Impact Category") +
  theme_minimal()

ggplot(Rxy_Italian, aes(x = Rxy, fill = Impact, color = Impact)) +
  geom_density(alpha = 0.5) +  # Overlay density curves
  labs( title = "Distribution of Rxy by Impact Category",
        x = "Rxy", y = "Density") +
  theme_minimal()         
