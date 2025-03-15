#### RUNS OF HOMOZYGOSITY ----

## Packages load ###
pacman::p_load(ggplot2,openxlsx,RColorBrewer,dplyr,grid,gridExtra,ggrepel,
               reshape2,tidyr,paletteer,ggtree,stringr,gtools)

# 1. PLINK ----
# 1.1 Full data preparation ----
# Load origin data and ROH
Lizards <- readWorkbook("Data/Samples_Santiago.xlsx", sheet = 2)
Italian <- Lizards %>% filter(Origin %in% c("Int-ITA", "Nat-ITA"))
French <- Lizards %>% filter(Origin %in% c("Int-FRA", "Nat-FRA"))

#Load ROH data and edit table names
Ita_ROHs <-read.table("Data/PopGen/All_Italian_ROHs_1.hom",h=T)
Fra_ROHs <- read.table("Data/PopGen/All_French_ROHs_1.hom",h=T)
Ita_ROHs <- Ita_ROHs %>% select(-FID)
Ita_ROHs <- Ita_ROHs %>% rename(ID = 1)
Fra_ROHs <- Fra_ROHs %>% select(-FID)
Fra_ROHs <- Fra_ROHs %>% rename(ID = 1)

# Merge datasets 
Ita_ROHs<- merge(Italian,Ita_ROHs,by="ID")
Fra_ROHs<- merge(French,Fra_ROHs,by="ID")

# Remove the mixed origin sample
Fra_ROHs <- Fra_ROHs %>% filter(ID != "CW2214")

## Ploting 
# Histograms
# Italian 
ggplot(Ita_ROHs, aes(x = KB, fill = Origin)) +
  geom_histogram(bins = 30, color = "black") + 
  facet_wrap(~ Origin) + 
  scale_fill_manual(values = c("Nat-ITA" = "#006837", 
                               "Int-ITA" = "#66BD63")) + 
  labs(title = "Italian ROHs length - PLINK",
    x = "KB",
    y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
# French 
ggplot(Fra_ROHs, aes(x = KB, fill= Origin)) +
  geom_histogram(bins = 30, color = "black") + 
  facet_wrap(~ Origin) + 
  scale_fill_manual(values= c("Nat-FRA" = "#A50026", 
                              "Int-FRA" = "#F46D43"))+
  labs( title = "French ROHs length - PLINK",
        x = "KB",
        y = "Count")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

#  ***** Boxplots per ind ----
# Italian
ggplot(Ita_ROHs, aes(x = ID, y = KB, fill = Origin)) +
  geom_boxplot() +
  labs(title = "BCF length per sample",
       x = "ID",
       y = "Kb") +
  scale_fill_manual(values= c("Nat-ITA" = "#006837",
                              "Int-ITA" = "#66BD63"))+
  theme_minimal()

# French
ggplot(Fra_ROHs, aes(x = ID, y = KB, fill = Origin)) +
  geom_boxplot() +
  labs(title = "BCF length per sample",
       x = "ID",
       y = "Kb") +
  scale_fill_manual(values= c("Nat-FRA" = "#A50026", 
                              "Int-FRA" = "#F46D43"))+
  theme_minimal()

# ROH in the chromosomes per individual. 
#Transform chromosomes in factors to whole genome mapping. 

Ita_ROHs <- Ita_ROHs %>%
  mutate(CHR = factor(CHR, levels = sort(as.numeric(as.character(unique(CHR))))))
Fra_ROHs <- Fra_ROHs %>%
  mutate(CHR = factor(CHR, levels = sort(as.numeric(as.character(unique(CHR))))))

# Plot accordingly to genome position. 
ggplot(Ita_ROHs, aes(x=POS1, xend=POS2, y=ID, color=as.factor(Origin))) +
  geom_segment(aes(yend=ID), size=4) +  
  facet_wrap(~ CHR, scales="free_x") +   
  scale_color_manual(values= c("Nat-ITA" = "#006837", "Int-ITA" = "#66BD63")) +
  theme_minimal() +
  labs(x="Genomic Position", y="Sample", title="Italian ROH Regions per Chromosome and Individual - PLINK", color="Origin") +
  theme(strip.text = element_text(size=12), axis.text.y = element_text(size=8))+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(Fra_ROHs, aes(x=POS1, xend=POS2, y=ID, color=as.factor(Origin))) +
  geom_segment(aes(yend=ID), size=4) +  # Draw ROH regions
  facet_wrap(~ CHR, scales="free_x") +   # Separate chromosomes in panels, allow independent x-axis
  scale_color_manual(values= c("Nat-FRA" = "#A50026","Int-FRA" = "#F46D43"))+
  theme_minimal() +
  labs(x="Genomic Position", y="Sample", title="French ROH Regions per Chromosome and Individual - PLINK", color="Origin" ) +
  theme(strip.text = element_text(size=12), axis.text.y = element_text(size=8))+
  theme(plot.title = element_text(hjust = 0.5))

# 1.2 Summary data preparation ----
# Italian 
Ita_Sum_ROH <- read.table("Data/PopGen/All_Italian_ROHs_1.hom.indiv",h=T)
Ita_Sum_ROH <- Ita_Sum_ROH %>% select(-FID)
Ita_Sum_ROH <- Ita_Sum_ROH %>% rename(ID = 1)

Ita_Sum_ROH<- merge(Italian,Ita_Sum_ROH,by="ID")

# Calculation of FROH 
# Convert KB to BP
Ita_Sum_ROH$Total_ROH_BP <- Ita_Sum_ROH$KB * 1000
Ita_Sum_ROH$Total_Autosomal_BP <- 1423936391  # Total autosomal sites based on the sum of chrs (NCBI) 
# Calculate FROH per individual
Ita_Sum_ROH$FROH <- Ita_Sum_ROH$Total_ROH_BP / Ita_Sum_ROH$Total_Autosomal_BP

# Plot FROH 
ggplot(Ita_Sum_ROH, aes(x = Origin, y = FROH , fill = Origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Boxplot to visualize distribution
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black")+
  scale_fill_manual(values= c("Nat-ITA" = "#006837",
                             "Int-ITA" = "#66BD63")) +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 20) +
  labs(title = "Italian Inbreeding coefficient - PLINK", x = "Origin", y = "FROH") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# Plot ROH length 
ggplot(Ita_Sum_ROH, aes(x = Origin, y = KBAVG , fill = Origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black")+
  scale_fill_manual(values= c("Nat-ITA" = "#006837",
                              "Int-ITA" = "#66BD63")) +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 20) +
  labs(title = "Italian ROH length - PLINK", x = "Origin", y = "ROH length") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))
  
# French  
Fra_Sum_ROH <- read.table("Data/PopGen/All_French_ROHs_1.hom.indiv",h=T)
Fra_Sum_ROH <- Fra_Sum_ROH %>% select(-FID)
Fra_Sum_ROH <- Fra_Sum_ROH %>% rename(ID = 1)

Fra_Sum_ROH<- merge(French,Fra_Sum_ROH,by="ID")
Fra_Sum_ROH <- Fra_Sum_ROH %>% filter(ID != "CW2214")
# Calculation of FROH 
# Convert KB to BP
Fra_Sum_ROH$Total_ROH_BP <- Fra_Sum_ROH$KB * 1000
Fra_Sum_ROH$Total_Autosomal_BP <- 1423936391  # Total record in the VCF 
# Calculate FROH per individual
Fra_Sum_ROH$FROH <- Fra_Sum_ROH$Total_ROH_BP / Fra_Sum_ROH$Total_Autosomal_BP

# Plot FROH 
ggplot(Fra_Sum_ROH, aes(x = Origin, y = FROH , fill = Origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Boxplot to visualize distribution
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black")+
  scale_fill_manual(values= c("Nat-FRA" = "#A50026", 
                               "Int-FRA" = "#F46D43"))+
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 20) +
  labs(title = "French Inbreeding coefficient - PLINK", x = "Origin", y = "FROH") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))



# 2.BCFtools ----

# 2.1 Full Data preparation ----

# Italian
Ita_BCF <- read.table("Data/PopGen/All_Italian_ROH_BCFtools",h=F)

# Naming the columns and filtering for quality (Phred score) and minimum length (500k)
Ita_BCF <- Ita_BCF %>%
  select(-1) %>% # Remove the first column - a bunch of RG
  rename(ID = V2, CHR = V3, POS1 = V4, POS2 = V5, BP = V6, MARKERS = V7, QUALITY = V8) %>% 
  filter(QUALITY >= 90, BP >= 500000)  

# Merge the clean data with the origin
Ita_BCF<- merge(Italian,Ita_BCF,by="ID")

# Define names  of chromosome names to  numbers
Ita_BCF <- Ita_BCF %>%
  mutate(CHR = case_when(
    CHR == "CM014743.1" ~ "1",
    CHR == "CM014744.1" ~ "2",
    CHR == "CM014745.1" ~ "3",
    CHR == "CM014746.1" ~ "4",
    CHR == "CM014747.1" ~ "5",
    CHR == "CM014748.1" ~ "6",
    CHR == "CM014749.1" ~ "7",
    CHR == "CM014750.1" ~ "8",
    CHR == "CM014751.1" ~ "9",
    CHR == "CM014752.1" ~ "10",
    CHR == "CM014753.1" ~ "11",
    CHR == "CM014754.1" ~ "12",
    CHR == "CM014755.1" ~ "13",
    CHR == "CM014756.1" ~ "14",
    CHR == "CM014757.1" ~ "15",
    CHR == "CM014758.1" ~ "16",
    CHR == "CM014759.1" ~ "17",
    CHR == "CM014760.1" ~ "18",))

ggplot(Ita_BCF, aes(x = BP, fill = Origin)) +
  geom_histogram(bins = 30, color = "black") + 
  facet_wrap(~ Origin) + 
  scale_fill_manual(values = c("Nat-ITA" = "#006837", 
                               "Int-ITA" = "#66BD63")) + 
  labs(
    title = "Italian ROHs length - BCFtools",
    x = "PB",
    y = "Count") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# French 

Fra_BCF <- read.table("Data/PopGen/All_French_ROH_BCFtools",h=F)
Fra_BCF <- Fra_BCF %>%
  select(-1) %>% # Remove the first column - a bunch of RG
  rename(ID = V2, CHR = V3, POS1 = V4, POS2 = V5, BP = V6, MARKERS = V7, QUALITY = V8) %>% 
  filter(QUALITY >= 90, BP >= 500000)  

# Merge the clean data with the origin
Fra_BCF<- merge(French,Fra_BCF,by="ID")
Fra_BCF <- Fra_BCF %>% filter(ID != "CW2214")
# Change chromosome names 
Fra_BCF <- Fra_BCF %>%
  mutate(CHR = case_when(
    CHR == "CM014743.1" ~ "1",
    CHR == "CM014744.1" ~ "2",
    CHR == "CM014745.1" ~ "3",
    CHR == "CM014746.1" ~ "4",
    CHR == "CM014747.1" ~ "5",
    CHR == "CM014748.1" ~ "6",
    CHR == "CM014749.1" ~ "7",
    CHR == "CM014750.1" ~ "8",
    CHR == "CM014751.1" ~ "9",
    CHR == "CM014752.1" ~ "10",
    CHR == "CM014753.1" ~ "11",
    CHR == "CM014754.1" ~ "12",
    CHR == "CM014755.1" ~ "13",
    CHR == "CM014756.1" ~ "14",
    CHR == "CM014757.1" ~ "15",
    CHR == "CM014758.1" ~ "16",
    CHR == "CM014759.1" ~ "17",
    CHR == "CM014760.1" ~ "18",))

# Plot 
ggplot(Fra_BCF, aes(x = BP, fill= Origin)) +
  geom_histogram(bins = 30, color = "black") + 
  facet_wrap(~ Origin) + 
  scale_fill_manual(values= c("Nat-FRA" = "#A50026", 
                              "Int-FRA" = "#F46D43"))+
  labs(title = "French ROHs length - BCFtools",
    x = "BP",
    y = "Count") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# ROHs in the genome

Ita_BCF <- Ita_BCF %>%
  mutate(CHR = factor(CHR, levels = sort(as.numeric(as.character(unique(CHR))))))
Fra_BCF <- Fra_BCF %>%
  mutate(CHR = factor(CHR, levels = sort(as.numeric(as.character(unique(CHR))))))

ggplot(Ita_BCF, aes(x=POS1, xend=POS2, y=ID, color=as.factor(Origin))) +
  geom_segment(aes(yend=ID), size=2) +  # Draw ROH regions
  facet_wrap(~ CHR, scales="free_x") +   # Separate chromosomes in panels, allow independent x-axis
  scale_color_manual(values= c("Nat-ITA" = "#006837",
                               "Int-ITA" = "#66BD63")) +
  theme_minimal() +
  labs(x="Genomic Position", y="Sample", title="Italian ROH Regions per Chromosome and Individual- BCFtools", color="Origin") +
  theme(strip.text = element_text(size=12), axis.text.y = element_text(size=8))+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(Fra_BCF, aes(x=POS1, xend=POS2, y=ID, color=as.factor(Origin))) +
  geom_segment(aes(yend=ID), size=2) +  # Draw ROH regions
  facet_wrap(~ CHR, scales="free_x") +   # Separate chromosomes in panels, allow independent x-axis
  scale_color_manual(values= c("Nat-FRA" = "#A50026", 
                               "Int-FRA" = "#F46D43"))+
  theme_minimal() +
  labs(x="Genomic Position", y="Sample", title="French ROH Regions per Chromosome and Individual- BCFtools", color="Origin")  +
  theme(strip.text = element_text(size=12), axis.text.y = element_text(size=8))+
  theme(plot.title = element_text(hjust = 0.5))

# 2.2 Summary preparation----

# Italian Summary 
Ita_BCF_Summary <- Ita_BCF %>%
  group_by(ID,Origin) %>%
  summarise(
    Total_ROH_BP = sum(BP),  # Total length of ROHs
    Avg_ROH_BP = mean(BP),   # Average ROH length
    Total_Autosomal_BP =  1423936391,   # Constant column
    FROH = Total_ROH_BP / Total_Autosomal_BP) # FROH Calculation


# Plot 
ggplot(Ita_BCF_Summary, aes(x = Origin, y = FROH , fill = Origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Boxplot to visualize distribution
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black")+
  scale_fill_manual(values= c("Nat-ITA" = "#006837",
                              "Int-ITA" = "#66BD63")) +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 20) +
  labs(title = "Italian Inbreeding coefficient - BCFtools", x = "Origin", y = "FROH") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# French Summary 
Fra_BCF_Summary <- Fra_BCF %>%
  group_by(ID,Origin) %>%
  summarise(
    Total_ROH_BP = sum(BP),  # Total length of ROHs
    Avg_ROH_BP = mean(BP),   # Average ROH length
    Total_Autosomal_BP =  1423936391,   # Constant column
    FROH = Total_ROH_BP / Total_Autosomal_BP) 

# Plot 
ggplot(Fra_BCF_Summary, aes(x = Origin, y = FROH , fill = Origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Boxplot to visualize distribution
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black")+
  scale_fill_manual(values= c("Nat-FRA" = "#A50026", 
                              "Int-FRA" = "#F46D43"))+
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 20) +
  labs(title = "French Inbreeding coefficient - BCFtools", x = "Origin", y = "FROH") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

