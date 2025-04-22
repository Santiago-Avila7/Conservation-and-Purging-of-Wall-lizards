#### Visualization of the PCA ---- #####
pacman::p_load(ggplot2,openxlsx,RColorBrewer,dplyr,grid,gridExtra,ggrepel,
               reshape2,tidyr,paletteer,ggtree,stringr,lmerTest,lme4,patchwork)

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
               "Int-ITA" = "Introduced Italian", 
               "Int-FRA" = "Introduced French")
  ) +
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
               "Int-ITA" = "Introduced Italian", 
               "Int-FRA" = "Introduced French")
  ) +
  theme_minimal() +
  theme(legend.position = "none", axis.title = element_text(size = 12))

# Combine the plots with a shared legend and centered title
final_plot <- p1 + p2 + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = "PCA", 
    theme = theme(plot.title = element_text(hjust = 0.5, size = 15))
  ) & theme(legend.position = "bottom")

# Display the final plot
print(final_plot)
