# Visualize the Gene diversity data ----
pacman::p_load(ggplot2,openxlsx,RColorBrewer,dplyr,gridExtra,ggrepel)

# 1. Heterozygocity ----

# Load the het data 
Het<- read.table("Data/PopGen/All_SNP_het.het",h=T)
Lizards <- readWorkbook("Data/Samples_Santiago.xlsx", sheet = 2)

# Use the select() function to drop the column
Het <- Het %>% select(-FID)
Het <- Het %>% rename(ID = 1)
Het$Exp_het<- 1-(Het$E.HOM./Het$N.NM.)
Het$Obs_het<- 1-(Het$O.HOM./Het$N.NM.)

Lizards_Het<- merge(Lizards,Het,by="ID")

# 1.1 Plots----
# Inbreeding coeficient
ggplot(Lizards_Het, aes(x = Origin, y = F , fill = Origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Boxplot to visualize distribution
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") + # Individual points
  scale_fill_manual(values= c("Nat-ITA" = "#006837",
                              "Nat-FRA" = "#A50026", 
                              "Int-ITA" = "#66BD63",
                              "Int-FRA" = "#F46D43")) +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 20) +
  labs(title = "Gene diversity", x = "Origin", y = "Inbreeding coeficient") +
  theme_minimal() +
  theme(legend.position = "none")

#Expected Heterozygocity
ggplot(Lizards_Het, aes(x = Origin, y = Exp_het , fill = Origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + 
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") + # Individual points
  scale_fill_manual(values= c("Nat-ITA" = "#006837",
                              "Nat-FRA" = "#A50026", 
                              "Int-ITA" = "#66BD63",
                              "Int-FRA" = "#F46D43")) +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 20) +
  labs(title = "Gene diversity", x = "Origin", y = "Expected heterozygocity") +
  theme_minimal() +
  theme(legend.position = "none")

# Observed Heterozygocity
ggplot(Lizards_Het, aes(x = Origin, y = Obs_het , fill = Origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + 
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") + # Individual points
  scale_fill_manual(values= c("Nat-ITA" = "#006837",
                              "Nat-FRA" = "#A50026", 
                              "Int-ITA" = "#66BD63",
                              "Int-FRA" = "#F46D43")) +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 20) +
  labs(title = "Gene diversity", x = "Origin", y = "Observed heterozygocity") +
  theme_minimal() +
  theme(legend.position = "none")


# 2. Nucleotide diversity ---- 
#Load data per population 

Int_Fra <- read.table("Data/PopGen/IntFra_nucleotide_diversity.windowed.pi", header = TRUE)
Int_Ita <- read.table("Data/PopGen/IntIta_nucleotide_diversity.windowed.pi", header = TRUE)
Nat_Fra <- read.table("Data/PopGen/NatFra_nucleotide_diversity.windowed.pi", header = TRUE)
Nat_Ita <- read.table("Data/PopGen/NatIta_nucleotide_diversity.windowed.pi", header = TRUE)

#Format info
Int_Fra$Origin <- "Int-FRA"
Int_Ita$Origin <- "Int-ITA"
Nat_Fra$Origin <- "Nat-FRA"
Nat_Ita$Origin <- "Nat-ITA"

combined_pi <- rbind(Int_Fra, Int_Ita, Nat_Fra, Nat_Ita)


ggplot(combined_pi, aes(x = Origin, y = PI, fill = Origin)) +
  geom_boxplot() +
  labs(title = "Nucleotide Diversity by Origin",
       x = "Origin",
       y = "Nucleotide Diversity (π)") +
  scale_fill_manual(values= c("Nat-ITA" = "#006837",
                              "Nat-FRA" = "#A50026", 
                              "Int-ITA" = "#66BD63",
                              "Int-FRA" = "#F46D43"))+
  theme_minimal()


# VCF tools confirmation  
Het_VCF<- read.table("Data/PopGen/All_SNP_VCFtools.het",h=T)
Het_VCF <- Het_VCF %>% rename(ID = 1)
Het_VCF$Exp_het<- 1-(Het_VCF$E.HOM./Het_VCF$N_SITES)
Lizards_H<- merge (Lizards, Het_VCF,by="ID")


ggplot(Lizards_H, aes(x = Origin, y = F , fill = Origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Boxplot to visualize distribution
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") + # Individual points
  scale_fill_manual(values= c("Nat-ITA" = "#006837",
                              "Nat-FRA" = "#A50026", 
                              "Int-ITA" = "#66BD63",
                              "Int-FRA" = "#F46D43")) +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 20) +
  labs(title = "Gene diversity", x = "Origin", y = "Inbreeding coeficient") +
  theme_minimal() +
  theme(legend.position = "none")
