# PURGING 2.0 ------
# Load important packages 
pacman::p_load(ggplot2,openxlsx,RColorBrewer,dplyr,grid,gridExtra,ggrepel,
               reshape2,tidyr,paletteer,ggtree,stringr,lmerTest,lme4,FSA,broom)
# Load data for Italian pops 2nd run ----
Purging_data <- read.table("Data/Purging/all_results.tsv", h=T)

# Calculate allele frequency in each Jackknife 
Purging_data<- Purging_data %>%
  mutate(High_Freq = HIGH/Total_Alleles,
         Moderate_Freq = MODERATE/Total_Alleles,
         Low_Freq=LOW/Total_Alleles,
         Modifier_Freq=MODIFIER/Total_Alleles,
         Intergenic_Freq=Intergenic_Alternate/Intergenic_Total)

# Select and make long for Rxy calculation
Purging_wide <- Purging_data %>%
  select(Jackknife, Origin, 
         High_Freq, Moderate_Freq, Low_Freq, 
         Modifier_Freq, Intergenic_Freq) %>%
  pivot_wider(names_from = Origin,
              values_from = c(High_Freq, Moderate_Freq,
                              Low_Freq, Modifier_Freq, Intergenic_Freq),
              names_glue = "{.value}_{Origin}")

# Create a function for Lxy ----
Lxy <- function(freq_category_X, freq_category_Y, freq_intergenic_X, freq_intergenic_Y) {
  numerator = freq_category_X * (1 - freq_category_Y)
  denominator = freq_intergenic_X * (1 - freq_intergenic_Y)
  numerator / denominator
}


# Calculate Lxy and Lyx per category 
Purging_wide<- Purging_wide %>%
  mutate(
    # For Introduced not Native 
    L_High_Int_Nat = Lxy (High_Freq_Introduced,High_Freq_Native,
                           Intergenic_Freq_Introduced,Intergenic_Freq_Native),
    L_Moderate_Int_Nat = Lxy (Moderate_Freq_Introduced,Moderate_Freq_Native,
                           Intergenic_Freq_Introduced,Intergenic_Freq_Native),
    L_Low_Int_Nat = Lxy (Low_Freq_Introduced,Low_Freq_Native,
                           Intergenic_Freq_Introduced,Intergenic_Freq_Native),
    L_Modifier_Int_Nat = Lxy (Modifier_Freq_Introduced,Modifier_Freq_Native,
                           Intergenic_Freq_Introduced,Intergenic_Freq_Native),
    
    # For Native not Introduced 
    L_High_Nat_Int= Lxy(High_Freq_Native,High_Freq_Introduced,
                         Intergenic_Freq_Native,Intergenic_Freq_Introduced),
    L_Moderate_Nat_Int= Lxy(Moderate_Freq_Native,Moderate_Freq_Introduced,
                         Intergenic_Freq_Native,Intergenic_Freq_Introduced),
    L_Low_Nat_Int= Lxy(Low_Freq_Native,Low_Freq_Introduced,
                         Intergenic_Freq_Native,Intergenic_Freq_Introduced),
    L_Modifier_Nat_Int= Lxy(Modifier_Freq_Native,Modifier_Freq_Introduced,
                         Intergenic_Freq_Native,Intergenic_Freq_Introduced),
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

# Stats 
# Shapiro-Wilk test for normality
shapiro.test(Rxy_Italian$Rxy)

# Visual check with Q-Q plot
qqnorm(Rxy_Italian$Rxy)
qqline(Rxy_Italian$Rxy)

# Applying non-parametric comparisons.
kruskal.test(Rxy ~ Impact, data = Rxy_Italian) # Significant 
dunnTest(Rxy ~ Impact, data = Rxy_Italian, method = "bonferroni")


#Test the difference against 1 
shapiro_results <- Rxy_Italian %>% 
  group_by(Impact) %>% 
  summarise(
    p_value = shapiro.test(Rxy)$p.value,
    .groups = "drop"
  )

print(shapiro_results)


# 3rd run ----
# Load data for Italian pops 
Purging_data <- read.table("Data/Purging/all_runs_3rd.tsv", h=T)

# Calculate allele frequency in each Jackknife 
Purging_data<- Purging_data %>%
  mutate(High_Freq = HIGH/Total_Alleles,
         Moderate_Freq = MODERATE/Total_Alleles,
         Low_Freq=LOW/Total_Alleles,
         Modifier_Freq=MODIFIER/Total_Alleles,
         Intergenic_Freq=Intergenic_Alternate/Intergenic_Total)

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

# Stats 

# Shapiro-Wilk test for normality
shapiro.test(Rxy_Italian$Rxy)

# Visual check with Q-Q plot
qqnorm(Rxy_Italian$Rxy)
qqline(Rxy_Italian$Rxy)

# Applying non-parametric comparisons.
kruskal.test(Rxy ~ Impact, data = Rxy_Italian) # Significant 
dunnTest(Rxy ~ Impact, data = Rxy_Italian, method = "bonferroni")


#Test the difference against 1 
shapiro_results <- Rxy_Italian %>% 
  group_by(Impact) %>% 
  summarise(
    p_value = shapiro.test(Rxy)$p.value,
    .groups = "drop"
  )

print(shapiro_results)

# 4th run ----
# Load data for Italian pops 
Purging_data <- read.table("Data/Purging/All_runs_4th.tsv", h=T)

# Calculate allele frequency in each Jackknife 
Purging_data<- Purging_data %>%
  mutate(High_Freq = HIGH/Total_Alleles,
         Moderate_Freq = MODERATE/Total_Alleles,
         Low_Freq=LOW/Total_Alleles,
         Modifier_Freq=MODIFIER/Total_Alleles,
         Intergenic_Freq=Intergenic_Alternate/Intergenic_Total)

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

# Stats 
# Shapiro-Wilk test for normality
shapiro.test(Rxy_Italian$Rxy)

# Visual check with Q-Q plot
qqnorm(Rxy_Italian$Rxy)
qqline(Rxy_Italian$Rxy)

# Applying non-parametric comparisons.
kruskal.test(Rxy ~ Impact, data = Rxy_Italian) # Significant 
dunnTest(Rxy ~ Impact, data = Rxy_Italian, method = "bonferroni")


#Test the difference against 1 
shapiro_results <- Rxy_Italian %>% 
  group_by(Impact) %>% 
  summarise(
    p_value = shapiro.test(Rxy)$p.value,
    .groups = "drop"
  )

print(shapiro_results)

# 5th run ----
# Load data for Italian pops 
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

# Stats ---- 
# Shapiro-Wilk test for normality
shapiro.test(Rxy_Italian$Rxy)

# Visual check with Q-Q plot
qqnorm(Rxy_Italian$Rxy)
qqline(Rxy_Italian$Rxy)

# Applying non-parametric comparisons.
kruskal.test(Rxy ~ Impact, data = Rxy_Italian) # Significant 
dunnTest(Rxy ~ Impact, data = Rxy_Italian, method = "bonferroni")


#Test the difference against 1 
shapiro_results <- Rxy_Italian %>% 
  group_by(Impact) %>% 
  summarise(
    p_value = shapiro.test(Rxy)$p.value,
    .groups = "drop"
  )

print(shapiro_results)

t_test_results <- Rxy_Italian %>%
  group_by(Impact) %>%
  do(tidy(t.test(.$Rxy, mu = 1))) %>%
  ungroup()

print(t_test_results)

#French inbreeding control (Removing DN)
# Load data 
Purging_data <- read.table("Data/Purging/All_French_IBC_jk.tsv", h=T)

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
Rxy_French<- Purging_wide %>%
  select(Rxy_High,
         Rxy_Moderate,
         Rxy_Low,
         Rxy_Modifier)

# Transform data 
Rxy_French <- Rxy_French %>% 
  pivot_longer(cols = starts_with("Rxy_"),  # Select all R_ columns
               names_to = "Impact",       # New column for impact categories
               values_to = "Rxy") %>%      # New column for Rxy values
  mutate(Impact = str_remove(Impact, "Rxy_")) # Clean up impact names

# Plots 
ggplot(Rxy_French, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot() +
  labs(title = "Rxy by Impact Category",
       x = "Rxy",
       y = "Impact Category") +
  theme_minimal() +
  scale_x_continuous(limits = c(0.94, 1.06))

ggplot(Rxy_Italian, aes(x = Rxy, fill = Impact, color = Impact)) +
  geom_density(alpha = 0.5) +  # Overlay density curves
  labs( title = "Distribution of Rxy by Impact Category",
        x = "Rxy", y = "Density") +
  theme_minimal()         

# Stats ---- 
# Shapiro-Wilk test for normality
shapiro.test(Rxy_Italian$Rxy)

# Visual check with Q-Q plot
qqnorm(Rxy_Italian$Rxy)
qqline(Rxy_Italian$Rxy)

# Applying non-parametric comparisons.
kruskal.test(Rxy ~ Impact, data = Rxy_Italian) # Significant 
dunnTest(Rxy ~ Impact, data = Rxy_Italian, method = "bonferroni")


#Test the difference against 1 
shapiro_results <- Rxy_Italian %>% 
  group_by(Impact) %>% 
  summarise(
    p_value = shapiro.test(Rxy)$p.value,
    .groups = "drop"
  )

print(shapiro_results)

t_test_results <- Rxy_Italian %>%
  group_by(Impact) %>%
  do(tidy(t.test(.$Rxy, mu = 1))) %>%
  ungroup()

print(t_test_results)
