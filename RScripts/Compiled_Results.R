######### Visualization of results ######
## Packages load ###
pacman::p_load(ggplot2,openxlsx,RColorBrewer,dplyr,grid,gridExtra,ggrepel,
               reshape2,tidyr,paletteer,ggtree,stringr,lmerTest,lme4,FSA)

#### 1. POPULATION STRUCTURE ####

## 1.1 PCA -----

# 1.1.1 Preparing the data set ----
#Load condensed sample data, sample + origin
Lizards <- readWorkbook("Data/Samples_Santiago.xlsx", sheet = 2)

# Using the PCA outputs from Plink
eval = read.table("Data/PopGen/All_Origins_PCA.eigenval", header=F)
evec = read.table("Data/PopGen/All_Origins_PCA.eigenvec", header=F)

# Merge the data with the origin. 
colnames(evec) <- c("ID", "RM", paste0("PC", 1:(ncol(evec) - 2)))
# Remove or ignore the redundant ID column
evec <- evec[, -2]  # removes the second column

# Calculate percentage variance explained for each PC
pve <- data.frame(PC = 1:20, pve = eval / sum(eval) * 100)
evec1.pc <- round(eval[1,1]/sum(eval)*100,digits=2)
evec2.pc <- round(eval[2,1]/sum(eval)*100,digits=2)
evec3.pc <- round(eval[3,1]/sum(eval)*100,digits=2)
evec4.pc <- round(eval[4,1]/sum(eval)*100,digits=2)
evec5.pc <- round(eval[5,1]/sum(eval)*100,digits=2)

# Merge PCA data and location
pca_data <- merge(evec, Lizards, by = "ID")

# 1.1.2 Plots including individual names and origin ----
p1<- ggplot(pca_data, aes(x = PC1, y = PC2, color = Origin)) +
  geom_point(size = 5) +
  labs(x = paste("PC1 (", evec1.pc, "%)", sep = ""),
       y = paste("PC2 (", evec2.pc, "%)", sep = "")) +
  scale_color_manual(values = c("Nat-ITA" = "#006837", 
                                "Nat-FRA" = "#A50026", 
                                "Int-ITA" = "#66BD63", 
                                "Int-FRA" = "#F46D43")) +
  theme_minimal()+
  theme(legend.position = "none", axis.title = element_text(size=20))


p2<- ggplot(pca_data, aes(x = PC1, y = PC3, color = Origin)) +
  geom_point(size = 5) +
  labs(x = paste("PC1 (", evec1.pc, "%)", sep = ""),
       y = paste("PC3 (", evec3.pc, "%)", sep = "")) +
  scale_color_manual(values = c("Nat-ITA" = "#006837", 
                                "Nat-FRA" = "#A50026", 
                                "Int-ITA" = "#66BD63", 
                                "Int-FRA" = "#F46D43")) +
 theme_minimal()+
  theme(legend.position = "none", axis.title = element_text(size=20))

p3<- ggplot(pca_data, aes(x = PC1, y = PC4, color = Origin)) +
  geom_point(size = 3) +
  labs(x = paste("PC1 (", evec1.pc, "%)", sep = ""),
       y = paste("PC4 (", evec4.pc, "%)", sep = "")) +
  scale_color_manual(values = c("Nat-ITA" = "#006837", 
                                "Nat-FRA" = "#A50026", 
                                "Int-ITA" = "#66BD63", 
                                "Int-FRA" = "#F46D43")) +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 20) +
  theme_minimal()+
  theme(legend.position = "bottom")

p4<- ggplot(pca_data, aes(x = PC1, y = PC5, color = Origin)) +
  geom_point(size = 3) +
  labs(x = paste("PC1 (", evec1.pc, "%)", sep = ""),
       y = paste("PC5 (", evec5.pc, "%)", sep = "")) +
  scale_color_manual(values = c("Nat-ITA" = "#006837", 
                                "Nat-FRA" = "#A50026", 
                                "Int-ITA" = "#66BD63", 
                                "Int-FRA" = "#F46D43")) +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 20) +
  theme_minimal()+
  theme(legend.position = "bottom")

grid.arrange(p1,p2, nrow = 1, ncol = 2,
             top = textGrob("PCA", gp = gpar(fontsize = 30)))


## 1.2. TREE ----

# 1.2.1 Data preparation ----
#Load the tree data
tree <- read.tree("Data/PopGen/All_Origins.treefile")

# Merge the tree data with the Lizard names and origin. 
tree_data <- full_join(as_tibble(tree), Lizards, by = c("label" = "ID"))



# 1.2.2 Plot the tree ----
p <- ggtree(tree, layout = "circular") %<+% tree_data +  # Plot the tree with the data. ( Circular or radial, creates a easy and not super ugly display)
  geom_tiplab(aes(color = Origin), size = 4, offset = 0.02, show.legend = F) + # Color tip labels by origin
  geom_tippoint(aes(color = Origin), size = 2, ) +       # Add colored points at tips
  scale_color_manual(values = c("Nat-ITA" = "#006837", 
                                "Nat-FRA" = "#A50026", 
                                "Int-ITA" = "#66BD63", 
                                "Int-FRA" = "#F46D43"),
                     labels = c("Introduced France", "Introduced Italy", 
                                "Native France", "Native Italy"))+ 
  guides(color = guide_legend (override.aes = list(size = 2),        # Make the color swatches larger
                               title = "Region of Origin")) + 
  theme(legend.position = "right", 
        legend.title = element_text(size=12),
        legend.text=element_text(size=10))
p


## 3. ADMIXTURE -----

# 1.3.1 Plot the CV error from the ADMIXTURE analysis----
cv_error <- read.table("Data/PopGen/CV_errors_summary.txt", h=F)
colnames (cv_error)<- c("rm","rm2","K","CV_error")
cv_error <- cv_error %>%
  mutate(K = str_extract(K, "\\d+"))
cv_error <- cv_error %>% select(-rm, -rm2)
plot(cv_error)

# 1.3.2 Organize the Data of K = 2,3 and 4 ----

# Load Admix data and merge with ORGANIZED population data 

# Load Admix data and merge with ORGANIZED population data 

Lizards_admix<- readWorkbook("Data/Samples_Santiago.xlsx", sheet = 3)
K2<- read.table("Data/PopGen/K2.Q", header=F)
K3<- read.table("Data/PopGen/K3.Q", header=F)
K4<- read.table("Data/PopGen/K4.Q", header=F)

# Create proper names for the columns in the Ks data sets and ID per sample.
K2<- cbind(Lizards_admix, K2)
colnames(K2) <- c("ID","Origin","Abbpop","Q1","Q2")
K3<- cbind(Lizards_admix, K3)
colnames(K3) <- c("ID","Origin","Abbpop","Q1","Q2","Q3")
K4<- cbind(Lizards_admix, K4)
colnames(K4) <- c("ID","Origin","Abbpop","Q1","Q2","Q3","Q4")

ordered_individuals <- K2 %>%
  arrange(-Q1) %>%
  pull(ID)  

# Apply ordering before pivoting
K2$ID <- factor(K2$ID, levels = ordered_individuals)
K3$ID <- factor(K3$ID, levels = ordered_individuals)
K4$ID <- factor(K4$ID, levels = ordered_individuals)

# Convert to long format
K2_long <- K2 %>% pivot_longer(cols = starts_with("Q"), names_to = "Ancestry", values_to = "Q_value")
K3_long <- K3 %>% pivot_longer(cols = starts_with("Q"), names_to = "Ancestry", values_to = "Q_value")
K4_long <- K4 %>% pivot_longer(cols = starts_with("Q"), names_to = "Ancestry", values_to = "Q_value")

# Add K identifier
K2_long$K <- "K2"
K3_long$K <- "K3"
K4_long$K <- "K4"

# Combine after ordering
admix_data <- bind_rows(K2_long, K3_long, K4_long)

# Re-apply factor ordering
admix_data$ID <- factor(admix_data$ID, levels = ordered_individuals)


# 1.3.3 Plot the ADMIXTURE ----
ggplot(admix_data, aes(x = ID, y = Q_value, fill = Ancestry)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_wrap(~K, ncol = 1) +
  scale_fill_manual(values = as.vector(paletteer_d("ggsci::default_jco"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Admixture Proportion", x = "Individual")


#### 2. GENETIC DIVERSITY ####

## 2.1. HETEROZYGOSITY ----

# Load the het data
Ita_Het<- read.table("Data/PopGen/All_Italian_het.tsv",h=T)
Fra_Het<- read.table("Data/PopGen/All_French_het.tsv",h=T)
Italian <- Lizards %>% filter(Origin %in% c("Int-ITA", "Nat-ITA"))
French <- Lizards %>% filter(Origin %in% c("Int-FRA", "Nat-FRA"))

# Merge data sets, remove sanity checks and the mixed origin sample
Ita_Het <- Ita_Het %>% rename(ID = 1)
Fra_Het <- Fra_Het %>% rename(ID = 1)
Ita_Het<- merge(Ita_Het,Italian,by="ID")
Fra_Het<- merge(Fra_Het,French,by="ID")
Ita_Het <- Ita_Het%>% select(-Total_Variant_Sites,-Missing_Genotypes)
Fra_Het <- Fra_Het%>% select(-Total_Variant_Sites,-Missing_Genotypes)
#Remove the sample CW2214 due to mixed origin 
Fra_Het <- Fra_Het[Fra_Het$ID != "CW2214", ]

# Calculate whole-genome observed heterozygosity based on the total amount of sites in the VCF per group. 
Ita_Het$Observed_Het<- Ita_Het$Heterozygous_Count/587596657 #Total number of Italian sites 
Fra_Het$Observed_Het<- Fra_Het$Heterozygous_Count/778409143 #Total number of French sites 

# Italian 
ggplot(Ita_Het, aes(x = Origin, y = Observed_Het , fill = Origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + 
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") + # Individual points
  scale_fill_manual(values= c("Nat-ITA" = "#006837",
                              "Int-ITA" = "#66BD63")) +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 20) +
  labs(title = "Whole-genome observed heterozygosity", x = "Origin", y = "Heterozygosity") +
  theme_minimal ()

# Significance test 
wilcox.test(Observed_Het ~ Origin, data = Ita_Het)

# French 
ggplot(Fra_Het, aes(x = Origin, y = Observed_Het , fill = Origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + 
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") + # Individual points
  scale_fill_manual(values= c("Nat-FRA" = "#A50026",
                              "Int-FRA" = "#F46D43")) +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 20) +
  labs(title = "Whole-genome observed heterozygosity", x = "Origin", y = "Heterozygosity") +
  theme_minimal() 

# Significance test 
wilcox.test(Observed_Het ~ Origin, data = Fra_Het)


# 2.2 ROHs ----
# No full up to date 

# Load ROH data and edit table names
Ita_ROHs <-read.table("Data/PopGen/All_Italian_ROHs_1.hom",h=T)
Fra_ROHs <- read.table("Data/PopGen/All_French_ROHs_1.hom",h=T)
Ita_ROHs <- Ita_ROHs %>% select(-FID)
Ita_ROHs <- Ita_ROHs %>% rename(ID = 1)
Fra_ROHs <- Fra_ROHs %>% select(-FID)
Fra_ROHs <- Fra_ROHs %>% rename(ID = 1)

Ita_ROHs<- merge(Italian,Ita_ROHs,by="ID")
Fra_ROHs<- merge(French,Fra_ROHs,by="ID")

# Remove the mixed origin sample
Fra_ROHs <- Fra_ROHs %>% filter(ID != "CW2214")


# Distribution 
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

# Distribution in chr 1 
Ita_ROHs <- Ita_ROHs %>%
  mutate(CHR = factor(CHR, levels = sort(as.numeric(as.character(unique(CHR))))))
Fra_ROHs <- Fra_ROHs %>%
  mutate(CHR = factor(CHR, levels = sort(as.numeric(as.character(unique(CHR))))))

# Filter for Chromosome 1
Ita_ROHs_chr1 <- Ita_ROHs %>% filter(CHR == 1)

# Plot only Chr 1
ggplot(Ita_ROHs_chr1, aes(x=POS1, xend=POS2, y=ID, color=as.factor(Origin))) +
  geom_segment(aes(yend=ID), linewidth =3) +  
  scale_color_manual(values= c("Nat-ITA" = "#006837", "Int-ITA" = "#66BD63")) +
  theme_minimal() +
  labs(x="Genomic Position", y="Sample", 
       title="Italian ROH Regions for Chromosome 1 - PLINK", 
       color="Origin") +
  theme(strip.text = element_text(size=12), 
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust = 0.5))

# Filter for Chromosome 1
Fra_ROHs_chr1 <- Fra_ROHs %>% filter(CHR == 1)

# Plot only Chr 1
ggplot(Fra_ROHs_chr1, aes(x=POS1, xend=POS2, y=ID, color=as.factor(Origin))) +
  geom_segment(aes(yend=ID), linewidth =3) +  
  scale_color_manual(values= c("Nat-FRA" = "#A50026","Int-FRA" = "#F46D43")) +
  theme_minimal() +
  labs(x="Genomic Position", y="Sample", 
       title="French ROH Regions for Chromosome 1 - PLINK", 
       color="Origin") +
  theme(strip.text = element_text(size=12), 
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust = 0.5))

# FROH 
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

# Test significance 
# Difference in the inbreeding coeficient
wilcox.test(FROH ~ Origin, data = Ita_Sum_ROH)
# Length of ROHs comparison 
Ita_model <- glmer(KB ~ Origin + (1 | ID), 
                   family = Gamma(link = "log"), 
                   data = Ita_ROHs)
summary(Ita_model)
anova(Ita_model)

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

# Test significance 
# Difference in the inbreeding coeficient
wilcox.test(FROH ~ Origin, data = Fra_Sum_ROH)
# Length of ROHs comparison 
Fra_model <- glmer(KB ~ Origin + (1 | ID), 
                   family = Gamma(link = "log"), 
                   data = Fra_ROHs)
summary(Fra_model)
anova(Fra_model)

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

# BCF tools.
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

# Distributions 

ggplot(Ita_BCF, aes(x = BP, fill = Origin)) +
  geom_histogram(bins = 30, color = "black") + 
  facet_wrap(~ Origin) + 
  scale_fill_manual(values = c("Nat-ITA" = "#006837", 
                               "Int-ITA" = "#66BD63")) + 
  labs(
    title = "Italian ROHs length - BCFtools",
    x = "BP",
    y = "Count") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

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
ggplot(Ita_BCF, aes(x = BP, fill = Origin)) +
  geom_histogram(bins = 30, color = "black") + 
  facet_wrap(~ Origin) + 
  scale_fill_manual(values = c("Nat-ITA" = "#006837", 
                               "Int-ITA" = "#66BD63")) + 
  labs(
    title = "Italian ROHs length - BCFtools",
    x = "BP",
    y = "Count") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

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

# Positions in chr 1 

Ita_BCF <- Ita_BCF %>%
  mutate(CHR = factor(CHR, levels = sort(as.numeric(as.character(unique(CHR))))))
Fra_BCF <- Fra_BCF %>%
  mutate(CHR = factor(CHR, levels = sort(as.numeric(as.character(unique(CHR))))))

# Filter for Chromosome 1
Ita_BCF_chr1 <- Ita_BCF %>% filter(CHR == 1)

# Plot only Chr 1
ggplot(Ita_BCF_chr1, aes(x=POS1, xend=POS2, y=ID, color=as.factor(Origin))) +
  geom_segment(aes(yend=ID), linewidth =3) +  
  scale_color_manual(values= c("Nat-ITA" = "#006837", "Int-ITA" = "#66BD63")) +
  theme_minimal() +
  labs(x="Genomic Position", y="Sample", 
       title="Italian ROH Regions for Chromosome 1 - BCFtools", 
       color="Origin") +
  theme(strip.text = element_text(size=12), 
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust = 0.5))

# Filter for Chromosome 1
Fra_BCF_chr1 <- Fra_BCF %>% filter(CHR == 1)

# Plot only Chr 1
ggplot(Fra_BCF_chr1, aes(x=POS1, xend=POS2, y=ID, color=as.factor(Origin))) +
  geom_segment(aes(yend=ID), linewidth =3) +  
  scale_color_manual(values= c("Nat-FRA" = "#A50026","Int-FRA" = "#F46D43")) +
  theme_minimal() +
  labs(x="Genomic Position", y="Sample", 
       title="French ROH Regions for Chromosome 1 - BCFtools", 
       color="Origin") +
  theme(strip.text = element_text(size=12), 
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust = 0.5))

# FROH 

# Italian Summary 
Ita_BCF_Summary <- Ita_BCF %>%
  group_by(ID,Origin) %>%
  summarise(
    Total_ROH_BP = sum(BP),  # Total length of ROHs
    Avg_ROH_BP = mean(BP),   # Average ROH length
    Total_Autosomal_BP =  1423936391,   # Constant column
    FROH = Total_ROH_BP / Total_Autosomal_BP) # FROH Calculation


# Test significance 
# Difference in the inbreeding coeficient
wilcox.test(FROH ~ Origin, data = Ita_BCF_Summary)
# Length of ROHs comparison 
Ita_BCF_model <- glmer(BP ~ Origin + (1 | ID), 
                       family = Gamma(link = "log"), 
                       data = Ita_BCF)
summary(Ita_BCF_model)
anova(Ita_BCF_model)

# Box Plot 
ggplot(Ita_BCF_Summary, aes(x = Origin, y = FROH , fill = Origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Boxplot to visualize distribution
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black")+
  scale_fill_manual(values= c("Nat-ITA" = "#006837",
                              "Int-ITA" = "#66BD63")) +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 20) +
  labs(title = "Italian Inbreeding coefficient - BCFtools", x = "Origin", y = "FROH") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# Reorder 
Ita_BCF_Summary$ID <- factor(Ita_BCF_Summary$ID, levels = Ita_BCF_Summary$ID[rev(order(Ita_BCF_Summary$FROH))])

# Create the barplot
ggplot(Ita_BCF_Summary, aes(y = ID, x = FROH, fill = Origin)) +
  geom_bar(stat = "identity") +
  labs(y = "Sample ID",
       x = "Percentage of ROH in Genome (FROH)",
       title = "Percentage of ROH in Genome per Sample") +
  scale_fill_manual(values = c("Nat-ITA" = "#006837", 
                               "Int-ITA" = "#66BD63")) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 15))

# French Summary 
Fra_BCF_Summary <- Fra_BCF %>%
  group_by(ID,Origin) %>%
  summarise(
    Total_ROH_BP = sum(BP),  # Total length of ROHs
    Avg_ROH_BP = mean(BP),   # Average ROH length
    Total_Autosomal_BP =  1423936391,   # Constant column
    FROH = Total_ROH_BP / Total_Autosomal_BP) 

# Test significance 
# Difference in the inbreeding coeficient
wilcox.test(FROH ~ Origin, data = Fra_BCF_Summary)
# Length of ROHs comparison 
Fra_BCF_model <- glmer(BP ~ Origin + (1 | ID), 
                       family = Gamma(link = "log"), 
                       data = Fra_BCF)
summary(Fra_BCF_model)
anova(Fra_BCF_model)

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

# Reorder 
Fra_BCF_Summary$ID <- factor(Fra_BCF_Summary$ID, levels = Fra_BCF_Summary$ID[rev(order(Fra_BCF_Summary$FROH))])

# Create the barplot
ggplot(Fra_BCF_Summary, aes(y = ID, x = FROH, fill = Origin)) +
  geom_bar(stat = "identity") +
  labs(y = "Sample ID",
       x = "Percentage of ROH in Genome (FROH)",
       title = "Percentage of ROH in Genome per Sample") +
  scale_fill_manual(values= c("Nat-FRA" = "#A50026", 
                              "Int-FRA" = "#F46D43"))+
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 15))
