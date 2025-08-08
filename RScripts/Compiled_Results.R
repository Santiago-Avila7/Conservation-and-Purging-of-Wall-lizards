######### Compilation of results and Visualization ######

## Packages ###
pacman::p_load(ggplot2,openxlsx,RColorBrewer,dplyr,grid,gridExtra,ggrepel,
               reshape2,tidyr,paletteer,ggtree,stringr,lmerTest,lme4,FSA,
               boot,patchwork,effsize,PMCMRplus, broom)

#### 1. POPULATION STRUCTURE ####

## 1.1 PCA -----

# 1.1.1 Preparing the data set ----
# Load condensed sample data, sample + origin
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
p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Origin)) +
  geom_point(size = 5) +
  labs(x = paste("PC1 (", evec1.pc, "%)", sep = ""),
       y = paste("PC2 (", evec2.pc, "%)", sep = "")) +
  scale_color_manual(
    values = c("Nat-ITA" = "#006837", 
               "Nat-FRA" = "#A50026", 
               "Int-ITA" = "#66BD63", 
               "Int-FRA" = "#F46D43"),
    labels = c("Nat-ITA" = "Native Italian", 
               "Nat-FRA" = "Native French", 
               "Int-ITA" = "Non-native Italian", 
               "Int-FRA" = "Non-native French")) +
  labs(subtitle = "A") +
  theme_minimal() +
  theme(legend.position = "bottom", axis.title = element_text(size = 12))

# Create p2 without legend but updated labels
p2 <- ggplot(pca_data, aes(x = PC1, y = PC3, color = Origin)) +
  geom_point(size = 5) +
  labs(x = paste("PC1 (", evec1.pc, "%)", sep = ""),
       y = paste("PC3 (", evec3.pc, "%)", sep = "")) +
  scale_color_manual(
    values = c("Nat-ITA" = "#006837", 
               "Nat-FRA" = "#A50026", 
               "Int-ITA" = "#66BD63", 
               "Int-FRA" = "#F46D43"),
    labels = c("Nat-ITA" = "Native Italian", 
               "Nat-FRA" = "Native French", 
               "Int-ITA" = "Non-native Italian", 
               "Int-FRA" = "Non-native French")) +
  labs(subtitle = "B") +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 12))

# Combine the plots with a shared label
PCA_plot <- p1 + p2 + 
  plot_layout(guides = "collect") + 
  plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5, size = 15))) &
  theme(legend.position = "bottom")

print(PCA_plot)


## 1.2. TREE ----

# 1.2.1 Data preparation ----
#Load the tree data
tree <- read.tree("Data/PopGen/All_Origins.treefile")

# Convert tree to tibble format
tree_data <- as_tibble(tree)

# Get correct node numbers for internal nodes
n_tips <- length(tree$tip.label)
internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)

# Create bootstrap data frame with proper node numbering
bootstrap_data <- data.frame(node = internal_nodes,
                             bootstrap = as.numeric(tree$node.label))# Convert to numeric

# Merge with main tree data
tree_data <- tree_data %>% 
  left_join(bootstrap_data, by = "node") %>% 
  # Merge with lizard metadata (tips only)
  left_join(Lizards, by = c("label" = "ID"))

# Create the tree plot
treeplot <- ggtree(tree, layout = "radial") %<+% tree_data +
  geom_tiplab2(aes(color = Origin), size = 3, offset = 0.03, align = TRUE, show.legend = FALSE) +
  geom_tippoint(aes(color = Origin), size = 2) +
  geom_nodelab(aes(label = round(bootstrap)), 
               color = "black", 
               hjust = -0.05,
               size = 2, 
               na.rm = TRUE,) +
  scale_color_manual(values = c("Nat-ITA" = "#006837", 
                                "Nat-FRA" = "#A50026", 
                                "Int-ITA" = "#66BD63", 
                                "Int-FRA" = "#F46D43"),
                     labels = c("Non-native France", "Non-native Italy", 
                                "Native France", "Native Italy")) +
  theme_minimal() +
  theme(axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)  )

print(treeplot)


## 1.3. ADMIXTURE -----

# 1.3.1 Plot the CV error from the ADMIXTURE analysis----
cv_error <- read.table("Data/PopGen/CV_errors_summary.txt", h=F)
colnames (cv_error)<- c("rm","rm2","K","CV_error")
cv_error <- cv_error %>%
  mutate(K = str_extract(K, "\\d+"))
cv_error <- cv_error %>% select(-rm, -rm2)
plot(cv_error) # Best supported K = 2 , 3 and 4


# 1.3.2 Organize the Data of K = 2,3 and 4 ----

# Load Admix data and merge with organized population data 

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

# Organize samples based on K2 (Best supported)
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


#### 2. GENETIC DIVERSITY AND INBREEDING ####

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

# Calculate whole-genome observed heterozygosity based on the total amount of sites in a full homozygous pseudogenome
Ita_Het$Observed_Het<- Ita_Het$Heterozygous_Count/1423936391 #Total number of autosomal sites 
Fra_Het$Observed_Het<- Fra_Het$Heterozygous_Count/1423936391 #Total number of autosomal sites 

# Italian 
ggplot(Ita_Het, aes(x = Origin, y = Observed_Het , fill = Origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + 
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") + # Individual points
  scale_fill_manual(values= c("Nat-ITA" = "#006837",
                              "Int-ITA" = "#66BD63")) +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 20) +
  labs(title = "Whole-genome observed heterozygosity", x = "Origin", y = "Heterozygosity") +
  theme_minimal ()

# Stats 
Ita_het_nat <- subset(Ita_Het, Origin == "Nat-ITA")$Observed_Het
Ita_het_int <- subset(Ita_Het, Origin == "Int-ITA")$Observed_Het

# Median difference
Ita_Het_diff <- median(Ita_het_nat) - median(Ita_het_int)

# Wilcoxon test
Ita_Het_wilcox <- wilcox.test(Observed_Het ~ Origin, data = Ita_Het)

# Cliff's delta
Ita_het_delta <- cliff.delta(Ita_het_nat, Ita_het_int)

# Output results
cat("Wilcoxon p:", Ita_Het_wilcox$p.value, "W:",Ita_Het_wilcox$statistic, "\n")
cat("Median difference (Nat-ITA - Int-ITA):", Ita_Het_diff, "\n")
print(Ita_het_delta)


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
Fra_het_nat <- subset(Fra_Het, Origin == "Nat-FRA")$Observed_Het
Fra_het_int <- subset(Fra_Het, Origin == "Int-FRA")$Observed_Het

# Median difference
Fra_Het_diff <- median(Fra_het_nat) - median(Fra_het_int)

# Wilcox test
Fra_Het_wilcox <- wilcox.test(Observed_Het ~ Origin, data = Fra_Het)

# Cliff's delta
Fra_het_delta <- cliff.delta(Fra_het_nat, Fra_het_int)

# Output results
cat("Wilcoxon p:", Fra_Het_wilcox$p.value, "W:",Fra_Het_wilcox$statistic, "\n")
cat("Median difference (Nat-Fra - Int-Fra):", Fra_Het_diff, "\n")
print(Fra_het_delta)

# Combined plot 
# Italian plot
italian_plot <- ggplot(Ita_Het, aes(x = Origin, y = Observed_Het, fill = Origin)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("Nat-ITA" = "#006837", 
                               "Int-ITA" = "#66BD63"),
                    labels = c("Nat-ITA" = "Native Italian", 
                               "Int-ITA" = "Non-native Italian")) +
  scale_x_discrete(
    labels = c("Nat-ITA" = "Native Italian", "Int-ITA" = "Non-native Italian")) +
  labs(subtitle = "A",
       x = "Origin",
       y = "Heterozygosity") +
  theme_minimal() +
  scale_y_continuous(limits = c(0.0007, 0.0019),
                     breaks = seq(0.0007, 0.0019, by = 0.0003)) +
  theme(plot.subtitle = element_text(size = 12),
        legend.position = "bottom")

# French plot
french_plot <- ggplot(Fra_Het, aes(x = Origin, y = Observed_Het, fill = Origin)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("Nat-FRA" = "#A50026", 
                               "Int-FRA" = "#F46D43"),
                    labels = c("Nat-FRA" = "Native French", 
                               "Int-FRA" = "Non-native French")) +
  scale_x_discrete(
    labels = c("Nat-FRA" = "Native French", "Int-FRA" = "Non-native French")) +
  labs(subtitle = "B",
       x = "Origin",
       y = "Heterozygosity") +
  theme_minimal() +
  scale_y_continuous(limits = c(0.0007, 0.0019),
                     breaks = seq(0.0007, 0.0019, by = 0.0003)) +
  theme(plot.subtitle = element_text(size = 12),
        legend.position = "bottom")

# Combine the plots
final_plot <- italian_plot + french_plot +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Display the combined plot
print(final_plot)


# 2.2 ROHs ----
# 2.2.1 BCF tools ----
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

# Count records per origin
Ita_BCF %>%
  count(Origin) # 1847 in Non-native pops, 73 in native pops 


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

# Count of records per origin 
Fra_BCF %>%
  count(Origin) # 1020 in Non-native pops and 561 in native pops  


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
  group_by(ID, Origin) %>%
  summarise(
    Total_ROH_BP = sum(BP),  # Total length of ROHs
    Avg_ROH_BP = mean(BP),   # Average ROH length
    ROH_Count = n(),         # Count of ROH records per sample
    Total_Autosomal_BP = 1423936391,  # Constant column - full hom Psudogenome 
    FROH = Total_ROH_BP / Total_Autosomal_BP, # FROH Calculation
    .groups = "drop")


# Test significance 
Ita_ROH_nat <- subset(Ita_BCF_Summary, Origin == "Nat-ITA")$FROH
Ita_ROH_int <- subset(Ita_BCF_Summary, Origin == "Int-ITA")$FROH

# Median differences
Ita_ROH_diff <- median(Ita_ROH_int) - median(Ita_ROH_nat)

# Wilcoxon test
Ita_ROH_wilcox <- wilcox.test(FROH ~ Origin, data = Ita_BCF_Summary)

# Cliff's delta
Ita_ROH_delta <- cliff.delta(Ita_ROH_int, Ita_ROH_nat)

# Output results
cat("Wilcoxon p:", Ita_ROH_wilcox$p.value, "W:",Ita_ROH_wilcox$statistic, "\n")
cat("Median difference (Nat-ITA - Int-ITA):", Ita_ROH_diff, "\n")
print(Ita_ROH_delta)

# Length of ROHs comparison 
Ita_BCF_model <- glmer(BP ~ Origin + (1 | ID), 
                       family = Gamma(link = "log"), 
                       data = Ita_BCF)
summary(Ita_BCF_model)

# Get the proportion of the effect 
exp(fixef(Ita_BCF_model)[-1]) #  ~ meaning 0.7X difference.

# Confidence intervals for the model 
# Function for bootstrapping 1000 replicates 
get_boot_ci <- function(model, term, nsim = 1000) {
  boot_fun <- function(.) {
    fixef(.)[term]
  }
  boot_ci <- bootMer(model, 
                     FUN = boot_fun, 
                     nsim = nsim,
                     type = "parametric")
  ci <- boot.ci(boot_ci, type = "perc", conf = 0.95)$percent[4:5]
  return(exp(ci))
}

# Italian confidence intervals
boot_ci_ita <- get_boot_ci(Ita_BCF_model, "OriginNat-ITA")
names(boot_ci_ita) <- c("Lower 95% CI", "Upper 95% CI")
boot_ci_ita # = 0.5056920 to 0.9194861  


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
Ita_BCF_Summary$ID <- factor(Ita_BCF_Summary$ID, levels = Ita_BCF_Summary$ID[rev(order(Ita_BCF_Summary$ROH_Count))])

# Create the barplot
ggplot(Ita_BCF_Summary, aes(y = ID, x = ROH_Count, fill = Origin)) +
  geom_bar(stat = "identity") +
  labs(y = "Sample ID",
       x = "Number of ROHs",
       title = "Total ROHs per Sample") +
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
    ROH_Count = n(), 
    Total_Autosomal_BP =  1423936391,   # Constant column full hom Psudogenome
    FROH = Total_ROH_BP / Total_Autosomal_BP) 

# Test significance 
# Difference in the inbreeding coeficient
Fra_ROH_nat <- subset(Fra_BCF_Summary, Origin == "Nat-FRA")$FROH
Fra_ROH_int <- subset(Fra_BCF_Summary, Origin == "Int-FRA")$FROH

# Median differences
Fra_ROH_diff <- median(Fra_ROH_int) - median(Fra_ROH_nat)

# Wilcoxon test
Fra_ROH_wilcox <- wilcox.test(FROH ~ Origin, data = Fra_BCF_Summary)

# Cliff's delta
Fra_ROH_delta <- cliff.delta(Fra_ROH_int, Fra_ROH_nat)

# Output results
cat("Wilcoxon p:", Fra_ROH_wilcox$p.value, "W:",Fra_ROH_wilcox$statistic, "\n")
cat("Median difference (Nat-Fra - Int-Fra):", Fra_ROH_diff, "\n")
print(Fra_ROH_delta)

# Length of ROHs comparison 
Fra_BCF_model <- glmer(BP ~ Origin + (1 | ID), 
                       family = Gamma(link = "log"), 
                       data = Fra_BCF)
summary(Fra_BCF_model)

# Proportion of difference  
exp(fixef(Fra_BCF_model)[-1])

# Confidence intervals 
# French model
boot_ci_fra <- get_boot_ci(Fra_BCF_model, "OriginNat-FRA")
names(boot_ci_fra) <- c("Lower 95% CI", "Upper 95% CI")
boot_ci_fra # = 0.5634890 to 0.8907801 

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
Fra_BCF_Summary$ID <- factor(Fra_BCF_Summary$ID, levels = Fra_BCF_Summary$ID[rev(order(Fra_BCF_Summary$ROH_Count))])

# Create the barplot
ggplot(Fra_BCF_Summary, aes(y = ID, x = ROH_Count, fill = Origin)) +
  geom_bar(stat = "identity") +
  labs(y = "Sample ID",
       x = "Number of ROHs",
       title = "Total ROHs per Sample") +
  scale_fill_manual(values= c("Nat-FRA" = "#A50026", 
                              "Int-FRA" = "#F46D43"))+
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 15))

# Combined plots
# Italian
# Boxplot
Ita_FROH_box <- ggplot(Ita_BCF_Summary, aes(x = Origin, y = FROH, fill = Origin)) +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c("Nat-ITA" = "#006837", "Int-ITA" = "#66BD63"),
    labels = c("Nat-ITA" = "Native Italian", "Int-ITA" = "Non-native Italian")) +
  scale_x_discrete(
    labels = c("Nat-ITA" = "Native Italian", "Int-ITA" = "Non-native Italian")) +
  labs(subtitle = "A",
       x = "Origin",
       y = "FROH") +
  theme_minimal() +
  theme(plot.subtitle = element_text(size = 12),
        legend.position = "none")

# Barplot
Ita_FROH_bar <- ggplot(Ita_BCF_Summary, aes(y = ID, x = ROH_Count, fill = Origin)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Nat-ITA" = "#006837", "Int-ITA" = "#66BD63"),
    labels = c("Nat-ITA" = "Native Italian", "Int-ITA" = "Non-native Italian")) +
  labs(subtitle = "B",
       y = "Sample",
       x = "Number of ROHs") +
  theme_minimal() +
  theme(plot.subtitle = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.position = "none")

# Combine Panels A and B into one top row 
Ita_FROH <- Ita_FROH_box + Ita_FROH_bar +
  plot_layout(guides = "collect") +
  theme(legend.position = "none")

# Order the samples in Panel C to match the order in Panel B 
ordered_ID_levels <- Ita_BCF_Summary %>% 
  arrange(desc(ROH_Count)) %>% 
  pull(ID) %>% 
  unique()

# Update your chromosome data so that the sample IDs become an ordered factor
Ita_BCF_chr1 <- Ita_BCF_chr1 %>%
  mutate(ID = factor(ID, levels = ordered_ID_levels))

# Chromosome 1 ROH segments 
Ita_ROH_chr1 <- ggplot(Ita_BCF_chr1, aes(x = POS1, xend = POS2, y = ID, color = as.factor(Origin))) +
  geom_segment(aes(yend = ID), linewidth = 3) +  
  scale_color_manual(values = c("Nat-ITA" = "#006837", "Int-ITA" = "#66BD63"),
    labels = c("Nat-ITA" = "Native Italian", "Int-ITA" = "Non-native Italian")) + 
  theme_minimal() +
  labs(subtitle = "C",
       x = "Genomic Position on Chromosome 1",
       y = "Sample",
       color = "Origin") +
  theme(strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 8),
    plot.subtitle = element_text(size = 12))

# Combine All Panels (A+B on top row, Panel C on bottom) ---
Ita_ALL <- Ita_FROH / Ita_ROH_chr1 +
  plot_layout(guides = "collect") 

# Display the final combined plot
print(Ita_ALL)

# French 
# Boxplot
Fra_FROH_box <- ggplot(Fra_BCF_Summary, aes(x = Origin, y = FROH, fill = Origin)) +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c("Nat-FRA" = "#A50026", "Int-FRA" = "#F46D43"),
    labels = c("Nat-FRA" = "Native French", "Int-FRA" = "Non-native French")) +
  scale_x_discrete(
    labels = c("Nat-FRA" = "Native French", "Int-FRA" = "Non-native French")) +
  labs(subtitle = "A",
       x = "Origin",
       y = "FROH") +
  theme_minimal() +
  theme(plot.subtitle = element_text(size = 12),
    legend.position = "none")

# Barplot for French samples
Fra_FROH_bar <- ggplot(Fra_BCF_Summary, aes(y = ID, x = ROH_Count, fill = Origin)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Nat-FRA" = "#A50026", "Int-FRA" = "#F46D43"),
    labels = c("Nat-FRA" = "Native French", "Int-FRA" = "Non-native French")) +
  labs(subtitle = "B",
       y = "Sample",
       x = "Number of ROHs") +
  theme_minimal() +
  theme(plot.subtitle = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    legend.position = "none")

# Combine Panels A and B into one top row
Fra_FROH <- Fra_FROH_box + Fra_FROH_bar +
  plot_layout(guides = "collect") +
  theme(legend.position = "none")

# Order the samples in Panel C to match the order in Panel B
ordered_ID_levels <- Fra_BCF_Summary %>% 
  arrange(desc(ROH_Count)) %>% 
  pull(ID) %>% 
  unique()

# Update the chromosome data so that the sample IDs become an ordered factor
Fra_BCF_chr1 <- Fra_BCF_chr1 %>%
  mutate(ID = factor(ID, levels = ordered_ID_levels))

# Chromosome 1 ROH segments for French samples
Fra_ROH_chr1 <- ggplot(Fra_BCF_chr1, aes(x = POS1, xend = POS2, y = ID, color = as.factor(Origin))) +
  geom_segment(aes(yend = ID), linewidth = 3) +  
  scale_color_manual(values = c("Nat-FRA" = "#A50026", "Int-FRA" = "#F46D43"),
    labels = c("Nat-FRA" = "Native French", "Int-FRA" = "Non-native French")) + 
  theme_minimal() +
  labs(subtitle = "C",
       x = "Genomic Position on Chromosome 1",
       y = "Sample",
       color = "Origin") +
  theme(strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 8),
    plot.subtitle = element_text(size = 12))

# Combine All Panels (Panels A+B on top, Panel C on bottom)
Fra_ALL <- Fra_FROH / Fra_ROH_chr1 +
  plot_layout(guides = "collect")

# Display the final combined plot
print(Fra_ALL)


# 2.2.2 PLINK ---- 
# PLINK was used to contrast with the results form BCFTools

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

# Length of ROHs comparison 
Ita_model <- glmer(KB ~ Origin + (1 | ID), 
                   family = Gamma(link = "log"), 
                   data = Ita_ROHs)
summary(Ita_model)

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

# Length of ROHs comparison 
Fra_model <- glmer(KB ~ Origin + (1 | ID), 
                   family = Gamma(link = "log"), 
                   data = Fra_ROHs)
summary(Fra_model)

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

# 3. PURGING  -----

# 3.1 Calculate Rxy for Italian-origin samples ----

# Load Italian Purging dataset - Relative frequencies and Jackknifing information
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

# Function to retrieve summary of the results after jackknifing 
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

# Visualize Italian-origin 
# Boxplot
Ita_plot<- ggplot(Rxy_Italian, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot(fill = "#66BD63",  
               outlier.shape = NA, 
               alpha = 0.8) +
  labs(subtitle =  "A - Italian",
       x = "Rxy",
       y = "Impact Category") +
  theme_minimal() +
  scale_x_continuous(limits = c(0.96, 1.04))

print (Ita_plot)

# Density Plot
ggplot(Rxy_Italian, aes(x = Rxy, fill = Impact, color = Impact)) +
  geom_density(alpha = 0.5) + 
  labs(title = "Distribution of Rxy by Impact Category - Italian",
       x = "Rxy", y = "Density") +
  theme_minimal() +
  scale_x_continuous(limits = c(0.96, 1.04))

# 3.2 Calculate Rxy for French-origin samples ----
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

# Visualize French-origin 
# BoxPlot 
Fra_plot <- ggplot(Rxy_French, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot(fill = "#F46D43",  
               outlier.shape = NA, 
               alpha = 0.8) +
  labs(subtitle = "B - French",
       x = "Rxy",
       y = "Impact Category") +
  theme_minimal() +
  scale_x_continuous(limits = c(0.96, 1.04))

print(Fra_plot)  

# Density plot 
ggplot(Rxy_French, aes(x = Rxy, fill = Impact, color = Impact)) +
  geom_density(alpha = 0.5) +  
  labs( title = "Distribution of Rxy by Impact Category - French",
        x = "Rxy", y = "Density") +
  theme_minimal() +
  scale_x_continuous(limits = c(0.96, 1.04))

# Confirm data as factors ans check min max values after jackknifing
Rxy_French$Impact <- as.factor(Rxy_French$Impact)
calculate_minmax(Rxy_French)

# Combined plot 
Purging_plot<- Ita_plot + Fra_plot +
  plot_layout(guides = "collect")

print (Purging_plot)


# 3.3 Calculate Rxy for each population ----

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

# Gather all the results 
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

# Generate plots per population
# Create Italian plot with green color
boxplot_Italian <- Rxy_all %>%
  filter(Group == "Italian") %>%
  ggplot(aes(x = Rxy, y = Impact)) +
  geom_boxplot(fill = "#66BD63",  # Green for Italian
               outlier.shape = NA, 
               alpha = 0.6) +
  facet_wrap(~ Population, nrow = 2) +
  theme_minimal() +
  scale_x_continuous(limits = c(0.9, 1.1)) +
  labs(x = "Rxy", y = "Impact Category") +
  theme(legend.position = "none")

# Create French plot with orange color
boxplot_French <- Rxy_all %>%
  filter(Group == "French") %>%
  ggplot(aes(x = Rxy, y = Impact)) +
  geom_boxplot(fill = "#F46D43",  # Orange for French
               outlier.shape = NA, 
               alpha = 0.6) +
  facet_wrap(~ Population, nrow = 2) +
  theme_minimal() +
  scale_x_continuous(limits = c(0.9, 1.1)) +
  labs(x = "Rxy", y = "Impact Category") +
  theme(legend.position = "none")

# Display plots
boxplot_Italian
boxplot_French


######## END #########