# PURGING 2.0 ------
# Load important packages 
pacman::p_load(ggplot2,openxlsx,RColorBrewer,dplyr,grid,gridExtra,ggrepel,
               reshape2,tidyr,paletteer,ggtree,stringr,lmerTest,lme4,FSA,broom)
# Load freq calculations 
Purging <- read.table("Data/Purging/Italian_freq_results.tsv", h=T)

# Get Lxy and Lyx 
Purging <- Purging %>%
  mutate(Lxy_HIGH = Sum_fab_HIGH/Sum_fab_INTERGENIC,
         Lyx_HIGH = Sum_fba_HIGH/Sum_fba_INTERGENIC,
         Lxy_MODERATE = Sum_fab_MODERATE/Sum_fab_INTERGENIC,
         Lyx_MODERATE = Sum_fba_MODERATE/Sum_fba_INTERGENIC,
         Lxy_LOW = Sum_fab_LOW/Sum_fab_INTERGENIC,
         Lyx_LOW = Sum_fba_LOW/Sum_fba_INTERGENIC,
         Lxy_MODIFIER = Sum_fab_MODIFIER/Sum_fab_INTERGENIC,
         Lyx_MODIFIER = Sum_fba_MODIFIER/Sum_fba_INTERGENIC,
         
         # Get Rxy per impact 
         Rxy_High = Lxy_HIGH/Lyx_HIGH,
         Rxy_Moderate = Lxy_MODERATE/Lyx_MODERATE,
         Rxy_Low = Lxy_LOW/Lyx_LOW,
         Rxy_Modifier = Lxy_MODIFIER/Lyx_MODIFIER)

# Filter just Rxy values
Rxy_Italian<- Purging %>%
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
  theme_minimal() +
  scale_x_continuous(limits = c(0.96, 1.04))

ggplot(Rxy_Italian, aes(x = Rxy, fill = Impact, color = Impact)) +
  geom_density(alpha = 0.5) +  # Overlay density curves
  labs( title = "Distribution of Rxy by Impact Category",
        x = "Rxy", y = "Density") +
  theme_minimal() +
  scale_x_continuous(limits = c(0.96, 1.04))
