# PURGING  ------
# Load important packages 
pacman::p_load(ggplot2,openxlsx,RColorBrewer,dplyr,grid,gridExtra,ggrepel,
               reshape2,tidyr,paletteer,ggtree,stringr,lmerTest,lme4,FSA)
# Load data 
Purging_data <- read.table("Data/Purging/Merged_jackknife.tsv", h=T)
Lizards <- readWorkbook("Data/Samples_Santiago.xlsx", sheet = 2)
Purging_data<- Purging_data %>% rename (ID=2)
Purging_full<-  Purging_data %>% 
  left_join(Lizards, by = "ID") %>% 
  relocate(Origin, .after = ID)

# Get all the intergenic sites data 
Intergenic_data<- read.table("Data/Purging/intergenic_allele_counts.tsv", h=T)
Intergenic_data<- Intergenic_data %>% rename (ID=1)
Intergenic_full<- merge (Intergenic_data, Lizards, by= "ID")

# Allele frequency of intergenic sites
Intergenic_AF <- Intergenic_full %>% 
  group_by(Origin) %>% 
  summarise(
    Total_Intergenic_Freq = sum(Total_Alternates)/sum(Total_Variants),
    .groups = "drop")

Intergenic_AF <- Intergenic_AF %>% 
  pivot_wider(
    names_from = Origin,              # Column to use for new column names
    values_from = Total_Intergenic_Freq)  # Column to use for values

# Modifier random runs (100k sites)
Modifier_data <- read.table("Data/Purging/merged_modifier.tsv", h=T)
Modifier_data<- Modifier_data %>% rename (ID=2)

# Merge 
Purging_full <- Purging_full %>%
  left_join(Modifier_data, by = c("Jackknife" = "Jackknife", "ID" = "ID"))

# Replace MODIFIER and Total_Alleles columns in the original dataset
Purging_full <- Purging_full %>%
  mutate( MODIFIER = ifelse(is.na(MODIFIER.y), MODIFIER.x, MODIFIER.y),
    Total_Alleles = Total_Alleles.x) %>%
  select(-MODIFIER.x, -MODIFIER.y, -Total_Alleles.x, -Total_Alleles.y)



#Summary of allele frequencies per origin and impact ----
Purging_sum <- Purging_full %>% 
  group_by(Jackknife, Origin) %>% 
  summarise(
    HIGH_freq = sum(HIGH) / sum(Total_Alleles),
    MODERATE_freq = sum(MODERATE) / sum(Total_Alleles),
    LOW_freq = sum(LOW) / sum(Total_Alleles),
    MODIFIER_freq = sum(MODIFIER) / sum(Total_Alleles),
    .groups = "drop"
  ) %>% 
  pivot_wider(names_from = Origin, 
              values_from = c(HIGH_freq, MODERATE_freq, LOW_freq, MODIFIER_freq),
              names_glue = "{.value}_{Origin}")

# Reorder columns for clarity
Purging_sum <- Purging_sum %>% 
  select(Jackknife, 
         starts_with("HIGH"), 
         starts_with("MODERATE"), 
         starts_with("LOW"), 
         starts_with("MODIFIER"))


# Functions

# Using the Jackknife data
Lxy <- function(freq_category_X, freq_category_Y, freq_intergenic_X, freq_intergenic_Y) {
  numerator = freq_category_X * (1 - freq_category_Y)
  denominator = freq_intergenic_X * (1 - freq_intergenic_Y)
  numerator / denominator
}

# Using all the intergenic sites  as denominator

# Get the frequency values of the frequency of all integenic sites 
fixed_intergenic <- c("Int-FRA" = 0.4168919,"Int-ITA" = 0.7652629,
                      "Nat-FRA" = 0.3275537,"Nat-ITA" = 0.7261919)

Lxy_fixed <- function(popX, popY, freq_category_X, freq_category_Y) {
  # Get fixed intergenic frequencies for the populations
  fixed_X <- fixed_intergenic[[popX]]
  fixed_Y <- fixed_intergenic[[popY]]
  
  # Calculate L(A,B)
  numerator = freq_category_X * (1 - freq_category_Y)
  denominator = fixed_X * (1 - fixed_Y)
  numerator / denominator
}


# Calculate Pairwise Lxy and Lyx and R Italian per Impact----

Purging_sum <- Purging_sum %>% 
  mutate(
    # For Int-IITA vs Nat-ITA
    L_High_ITA_IntNat = Lxy(`HIGH_freq_Int-ITA`,`HIGH_freq_Nat-ITA`,
                            `MODIFIER_freq_Int-ITA`,`MODIFIER_freq_Nat-ITA`),
    L_Moderate_ITA_IntNat = Lxy(`MODERATE_freq_Int-ITA`,`MODERATE_freq_Nat-ITA`,
                                `MODIFIER_freq_Int-ITA`,`MODIFIER_freq_Nat-ITA`),,
    L_Low_ITA_IntNat = Lxy(`LOW_freq_Int-ITA`,`LOW_freq_Nat-ITA`,
                           `MODIFIER_freq_Int-ITA`,`MODIFIER_freq_Nat-ITA`),,
    L_Modifier_ITA_IntNat = Lxy(`MODIFIER_freq_Int-ITA`,`MODIFIER_freq_Nat-ITA`,
                                `MODIFIER_freq_Int-ITA`, `MODIFIER_freq_Nat-ITA`),
    
    # For Int-ITA vs Nat-ITA
    L_High_ITA_NatInt = Lxy(`HIGH_freq_Nat-ITA`,`HIGH_freq_Int-ITA`,
                            `MODIFIER_freq_Nat-ITA`,`MODIFIER_freq_Int-ITA`),
    L_Moderate_ITA_NatInt = Lxy(`MODERATE_freq_Nat-ITA`,`MODERATE_freq_Int-ITA`,
                                `MODIFIER_freq_Nat-ITA`,`MODIFIER_freq_Int-ITA`),
    L_Low_ITA_NatInt = Lxy(`LOW_freq_Nat-ITA`,`LOW_freq_Int-ITA`,
                           `MODIFIER_freq_Nat-ITA`,`MODIFIER_freq_Int-ITA`),
    L_Modifier_ITA_NatInt = Lxy(`MODIFIER_freq_Nat-ITA`,`MODIFIER_freq_Int-ITA`,
                                `MODIFIER_freq_Nat-ITA`,`MODIFIER_freq_Int-ITA`),
    # For Rxy 
    Rxy_High_Ita = L_High_ITA_IntNat/L_High_ITA_NatInt,
    Rxy_Moderate_Ita = L_Moderate_ITA_IntNat/L_Moderate_ITA_NatInt,
    Rxy_Low_Ita = L_Low_ITA_IntNat/L_Low_ITA_NatInt,
    Rxy_Modifier_Ita = L_Modifier_ITA_IntNat/L_Modifier_ITA_NatInt)

Rxy_Italian<- Purging_sum %>%
  select(Rxy_High_Ita,
         Rxy_Moderate_Ita,
         Rxy_Low_Ita)

# Transform data 
Rxy_Italian <- Rxy_Italian %>% 
  pivot_longer(
    cols = starts_with("Rxy_"),  # Select all R_ columns
    names_to = "Impact",       # New column for impact categories
    values_to = "Rxy"          # New column for Rxy values
  ) %>% 
  mutate(
    Impact = str_remove(Impact, "Rxy_"),  # Clean up impact names
    Impact = str_remove(Impact, "_Ita")
  )

# Plot ---- 

ggplot(Rxy_Italian, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot() +
  labs(
    title = "Rxy by Impact Category",
    x = "Rxy",
    y = "Impact Category"
  ) +
  theme_minimal()

ggplot(Rxy_Italian, aes(x = Rxy, fill = Impact, color = Impact)) +
  geom_density(alpha = 0.5) +  # Overlay density curves
  labs(
    title = "Distribution of Rxy by Impact Category",
    x = "Rxy",
    y = "Density"
  ) +
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

# One-sample t-test
t_test_results <- Rxy_Italian %>% 
  group_by(Impact) %>% 
  summarise(
    t_statistic = t.test(Rxy, mu = 1)$statistic,
    p_value = t.test(Rxy, mu = 1)$p.value,
    .groups = "drop"
  )

print(t_test_results)


# Calculate Pairwise L(ab) and L(ba) and R Fra per Impact----

Purging_sum <- Purging_sum %>% 
  mutate(
    # For Int-FRA vs Nat-FRA
    L_High_FRA_IntNat = Lxy(`HIGH_freq_Int-FRA`,`HIGH_freq_Nat-FRA`,
                            `MODIFIER_freq_Int-FRA`,`MODIFIER_freq_Nat-FRA`),
    L_Moderate_FRA_IntNat = Lxy(`MODERATE_freq_Int-FRA`,`MODERATE_freq_Nat-FRA`,
                                `MODIFIER_freq_Int-FRA`,`MODIFIER_freq_Nat-FRA`),,
    L_Low_FRA_IntNat = Lxy(`LOW_freq_Int-FRA`,`LOW_freq_Nat-FRA`,
                           `MODIFIER_freq_Int-FRA`,`MODIFIER_freq_Nat-FRA`),,
    L_Modifier_FRA_IntNat = Lxy(`MODIFIER_freq_Int-FRA`,`MODIFIER_freq_Nat-FRA`,
                                `MODIFIER_freq_Int-FRA`, `MODIFIER_freq_Nat-FRA`),
    
    # For Int-FRA vs Nat-FRA
    L_High_FRA_NatInt = Lxy(`HIGH_freq_Nat-FRA`,`HIGH_freq_Int-FRA`,
                            `MODIFIER_freq_Nat-FRA`,`MODIFIER_freq_Int-FRA`),
    L_Moderate_FRA_NatInt = Lxy(`MODERATE_freq_Nat-FRA`,`MODERATE_freq_Int-FRA`,
                                `MODIFIER_freq_Nat-FRA`,`MODIFIER_freq_Int-FRA`),
    L_Low_FRA_NatInt = Lxy(`LOW_freq_Nat-FRA`,`LOW_freq_Int-FRA`,
                           `MODIFIER_freq_Nat-FRA`,`MODIFIER_freq_Int-FRA`),
    L_Modifier_FRA_NatInt = Lxy(`MODIFIER_freq_Nat-FRA`,`MODIFIER_freq_Int-FRA`,
                                `MODIFIER_freq_Nat-FRA`,`MODIFIER_freq_Int-FRA`),
    
    # For Rxy 
    Rxy_High_Fra = L_High_FRA_IntNat/L_High_FRA_NatInt,
    Rxy_Moderate_Fra = L_Moderate_FRA_IntNat/L_Moderate_FRA_NatInt,
    Rxy_Low_Fra = L_Low_FRA_IntNat/L_Low_FRA_NatInt,
    Rxy_Modifier_Fra = L_Modifier_FRA_IntNat/L_Modifier_FRA_NatInt)


Rxy_French<- Purging_sum %>%
  select(Rxy_High_Fra,
         Rxy_Moderate_Fra,
         Rxy_Low_Fra)

# Transform data 
Rxy_French <- Rxy_French %>% 
  pivot_longer(
    cols = starts_with("Rxy_"),  # Select all R_ columns
    names_to = "Impact",       # New column for impact categories
    values_to = "Rxy"          # New column for Rxy values
  ) %>% 
  mutate(
    Impact = str_remove(Impact, "Rxy_"),  # Clean up impact names
    Impact = str_remove(Impact, "_Fra")
  )

# Plot ---- 

ggplot(Rxy_French, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot() +
  labs(
    title = "Rxy by Impact Category",
    x = "Rxy",
    y = "Impact Category"
  ) +
  theme_minimal()

ggplot(Rxy_French, aes(x = Rxy, fill = Impact, color = Impact)) +
  geom_density(alpha = 0.5) +  # Overlay density curves
  labs(
    title = "Distribution of Rxy by Impact Category",
    x = "Rxy",
    y = "Density"
  ) +
  theme_minimal()

# Stats ---- 
# Shapiro-Wilk test for normality
shapiro.test(Rxy_French$Rxy)

# Visual check with Q-Q plot
qqnorm(Rxy_French$Rxy)
qqline(Rxy_French$Rxy)

# Non-parametric comparison 
kruskal.test(Rxy ~ Impact, data = Rxy_French)
dunn_results<-dunnTest(Rxy ~ Impact, data = Rxy_French, method = "bonferroni")
print(dunn_results)
# Test difference against 1 

shapiro_results <- Rxy_French %>% 
  group_by(Impact) %>% 
  summarise(
    p_value = shapiro.test(Rxy)$p.value,
    .groups = "drop"
  )

print(shapiro_results)

# One-sample t-test
t_test_results <- Rxy_French%>% 
  group_by(Impact) %>% 
  summarise(
    t_statistic = t.test(Rxy, mu = 1)$statistic,
    p_value = t.test(Rxy, mu = 1)$p.value,
    .groups = "drop"
  )

print(t_test_results)


# Size effect 
library(effsize)


High_values <- Rxy_French$Rxy[Rxy_French$Impact == "High"]
Low_values <- Rxy_French$Rxy[Rxy_French$Impact == "Low"]
Moderate_values <- Rxy_French$Rxy[Rxy_French$Impact == "Moderate"]


cliff.delta(High_values, Low_values)
cliff.delta(High_values, Moderate_values)
cliff.delta(Low_values, Moderate_values)




# Experiments ----

# Functions ----

Lxy <- function(freq_category_X, freq_category_Y, freq_intergenic_X, freq_intergenic_Y) {
  numerator = freq_category_X * (1 - freq_category_Y)
  denominator = freq_intergenic_X * (1 - freq_intergenic_Y)
  numerator / denominator
}

fixed_intergenic <- c("Int-FRA" = 0.4877882,"Int-ITA" = 0.4876248,
                      "Nat-FRA" = 0.4877862,"Nat-ITA" = 0.4876157)

Lxy_fixed <- function(popX, popY, freq_category_X, freq_category_Y) {
  # Get fixed intergenic frequencies for the populations
  fixed_X <- fixed_intergenic[[popX]]
  fixed_Y <- fixed_intergenic[[popY]]
  
  # Calculate L(A,B)
  numerator = freq_category_X * (1 - freq_category_Y)
  denominator = fixed_X * (1 - fixed_Y)
  numerator / denominator
}


Lxy(Purging_sum$`HIGH_freq_Int-ITA`,Purging_sum$`HIGH_freq_Nat-ITA`,
    Purging_sum$`MODIFIER_freq_Int-ITA`,Purging_sum$`MODERATE_freq_Nat-ITA`)

Lxy_fixed("Int-ITA","Nat-ITA",
          Purging_sum$`MODIFIER_freq_Int-ITA`,Purging_sum$`MODIFIER_freq_Nat-ITA`)


# Calculate Pairwise Lxy and Lyx and R Italian per Impact----

Purging_sum <- Purging_sum %>% 
  mutate(
    # For Int-ITA vs Nat-ITA
    L_High_ITA_IntNat = Lxy_fixed("Int-ITA", "Nat-ITA",
                                  `HIGH_freq_Int-ITA`, `HIGH_freq_Nat-ITA`),
    L_Moderate_ITA_IntNat = Lxy_fixed("Int-ITA", "Nat-ITA", 
                                      `MODERATE_freq_Int-ITA`, `MODERATE_freq_Nat-ITA`),
    L_Low_ITA_IntNat = Lxy_fixed("Int-ITA", "Nat-ITA", 
                                 `LOW_freq_Int-ITA`, `LOW_freq_Nat-ITA`),
    L_Modifier_ITA_IntNat = Lxy_fixed("Int-ITA", "Nat-ITA", 
                                      `MODIFIER_freq_Int-ITA`, `MODIFIER_freq_Nat-ITA`),
    
    # For Int-ITA vs Nat-ITA
    L_High_ITA_NatInt = Lxy_fixed("Nat-ITA", "Int-ITA",
                                  `HIGH_freq_Nat-ITA`, `HIGH_freq_Int-ITA`),
    L_Moderate_ITA_NatInt = Lxy_fixed("Nat-ITA", "Int-ITA",
                                      `MODERATE_freq_Nat-ITA`, `MODERATE_freq_Int-ITA`),
    L_Low_ITA_NatInt = Lxy_fixed("Nat-ITA", "Int-ITA",
                                 `LOW_freq_Nat-ITA`, `LOW_freq_Int-ITA`),
    L_Modifier_ITA_NatInt = Lxy_fixed("Nat-ITA", "Int-ITA",
                                      `MODIFIER_freq_Nat-ITA`, `MODIFIER_freq_Int-ITA`),
    
    # For Rxy 
    Rxy_High_Ita = L_High_ITA_IntNat/L_High_ITA_NatInt,
    Rxy_Moderate_Ita = L_Moderate_ITA_IntNat/L_Moderate_ITA_NatInt,
    Rxy_Low_Ita = L_Low_ITA_IntNat/L_Low_ITA_NatInt,
    Rxy_Modifier_Ita = L_Modifier_ITA_IntNat/L_Modifier_ITA_NatInt)

Rxy_Italian<- Purging_sum %>%
  select(Rxy_High_Ita,
         Rxy_Moderate_Ita,
         Rxy_Low_Ita,
         Rxy_Modifier_Ita)

# Transform data 
Rxy_Italian <- Rxy_Italian %>% 
  pivot_longer(
    cols = starts_with("Rxy_"),  # Select all R_ columns
    names_to = "Impact",       # New column for impact categories
    values_to = "Rxy"          # New column for Rxy values
  ) %>% 
  mutate(
    Impact = str_remove(Impact, "Rxy_"),  # Clean up impact names
    Impact = str_remove(Impact, "_Ita")
  )
  
ggplot(Rxy_Italian, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot() +
  labs(
    title = "Rxy by Impact Category",
    x = "Rxy",
    y = "Impact Category"
  ) +
  theme_minimal()

# French 
Purging_sum <- Purging_sum %>% 
  mutate(
    # For Int-FRA vs Nat-FRA
    L_High_FRA_IntNat = Lxy_fixed("Int-FRA", "Nat-FRA",
                                  `HIGH_freq_Int-FRA`, `HIGH_freq_Nat-FRA`),
    L_Moderate_FRA_IntNat = Lxy_fixed("Int-FRA", "Nat-FRA", 
                                      `MODERATE_freq_Int-FRA`, `MODERATE_freq_Nat-FRA`),
    L_Low_FRA_IntNat = Lxy_fixed("Int-FRA", "Nat-FRA", 
                                 `LOW_freq_Int-FRA`, `LOW_freq_Nat-FRA`),
    L_Modifier_FRA_IntNat = Lxy_fixed("Int-FRA", "Nat-FRA", 
                                      `MODIFIER_freq_Int-FRA`, `MODIFIER_freq_Nat-FRA`),
    
    # For Int-FRA vs Nat-FRA
    L_High_FRA_NatInt = Lxy_fixed("Nat-FRA", "Int-FRA",
                                  `HIGH_freq_Nat-FRA`, `HIGH_freq_Int-FRA`),
    L_Moderate_FRA_NatInt = Lxy_fixed("Nat-FRA", "Int-FRA",
                                      `MODERATE_freq_Nat-FRA`, `MODERATE_freq_Int-FRA`),
    L_Low_FRA_NatInt = Lxy_fixed("Nat-FRA", "Int-FRA",
                                 `LOW_freq_Nat-FRA`, `LOW_freq_Int-FRA`),
    L_Modifier_FRA_NatInt = Lxy_fixed("Nat-FRA", "Int-FRA",
                                      `MODIFIER_freq_Nat-FRA`, `MODIFIER_freq_Int-FRA`),
    
    # For Rxy 
    Rxy_High_FRA = L_High_FRA_IntNat/L_High_FRA_NatInt,
    Rxy_Moderate_FRA = L_Moderate_FRA_IntNat/L_Moderate_FRA_NatInt,
    Rxy_Low_FRA = L_Low_FRA_IntNat/L_Low_FRA_NatInt,
    Rxy_Modifier_FRA = L_Modifier_FRA_IntNat/L_Modifier_FRA_NatInt)

Rxy_French<- Purging_sum %>%
  select(Rxy_High_FRA,
         Rxy_Moderate_FRA,
         Rxy_Low_FRA,
         Rxy_Modifier_FRA)

# Transform data 
Rxy_French <- Rxy_French %>% 
  pivot_longer(
    cols = starts_with("Rxy_"),  # Select all R_ columns
    names_to = "Impact",       # New column for impact categories
    values_to = "Rxy"          # New column for Rxy values
  ) %>% 
  mutate(
    Impact = str_remove(Impact, "Rxy_"),  # Clean up impact names
    Impact = str_remove(Impact, "_Ita")
  )

ggplot(Rxy_French, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot() +
  labs(
    title = "Rxy by Impact Category",
    x = "Rxy",
    y = "Impact Category"
  ) +
  theme_minimal()


ggplot(Rxy_French, aes(x = Rxy, fill = Impact, color = Impact)) +
  geom_density(alpha = 0.5) +  # Overlay density curves
  labs(
    title = "Distribution of Rxy by Impact Category",
    x = "Rxy",
    y = "Density"
  ) +
  theme_minimal()

#Stats different from 0 -----
t_test_results <- Rxy_Italian %>%
  group_by(Impact) %>%
  summarise(
    t_test_p = t.test(Rxy, mu = 1)$p.value,  # One-sample t-test
    mean_Rxy = mean(Rxy),
    .groups = "drop"
  )

print(t_test_results)

pairwise.wilcox.test(Rxy_Italian$Rxy, Rxy_Italian$Impact, p.adjust.method = "bonferroni")

