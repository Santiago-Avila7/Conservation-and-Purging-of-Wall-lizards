#### Visualization of the PCA ---- #####
pacman::p_load(ggplot2,openxlsx,RColorBrewer,dplyr,gridExtra,ggrepel)

#Load sample data 
Lizards <- readWorkbook("Data/Samples_Santiago.xlsx", sheet = 2)

#Using the PCA outputs from the plink code
eval = read.table("Data/PopGen/Plink_PCA.eigenval", header=F)
evec = read.table("Data/PopGen/Plink_PCA.eigenvec", header=F)

#Merge the data with the origin. 
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

#Plot using location information 
p1<- ggplot(pca_data, aes(x = PC1, y = PC2, color = Origin)) +
  geom_point(size = 3) +
  labs(x = paste("PC1 (", evec1.pc, "%)", sep = ""),
       y = paste("PC2 (", evec2.pc, "%)", sep = "")) +
  scale_color_manual(values = c("Nat-ITA" = "#006837", 
                                "Nat-FRA" = "#A50026", 
                                "Int-ITA" = "#66BD63", 
                                "Int-FRA" = "#F46D43")) +
  theme_minimal()+
  theme(legend.position = "bottom")

p2<- ggplot(pca_data, aes(x = PC1, y = PC3, color = Origin)) +
  geom_point(size = 3) +
  labs(x = paste("PC1 (", evec1.pc, "%)", sep = ""),
       y = paste("PC3 (", evec3.pc, "%)", sep = "")) +
  scale_color_manual(values = c("Nat-ITA" = "#006837", 
                                "Nat-FRA" = "#A50026", 
                                "Int-ITA" = "#66BD63", 
                                "Int-FRA" = "#F46D43")) +
  theme_minimal()+
  theme(legend.position = "bottom")

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


grid.arrange(p1,p2, nrow = 1, ncol = 2,top="PCA")





