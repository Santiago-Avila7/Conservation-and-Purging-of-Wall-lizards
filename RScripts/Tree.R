# See the tree
pacman::p_load(ggplot2,openxlsx,RColorBrewer,dplyr,gridExtra,reshape2,tidyr,paletteer,ggtree)

#Load the tree data 
tree_file <- "Data/PopGen/All_SNP_final.min4.phy.varsites.phy.treefile"  # Replace with your tree file
tree <- read.tree(tree_file)

# Load Lizard data 
Lizards <- readWorkbook("Data/Samples_Santiago.xlsx", sheet = 2)

# Merge metadata with tree tip labels
# Ensure the 'Sample' column in metadata matches tip labels in your tree
tree_data <- full_join(as_tibble(tree), Lizards, by = c("label" = "ID"))

# Visualize the tree
p <- ggtree(tree, layout = "radial") %<+% tree_data +  # Attach metadata to the tree
  geom_tiplab(aes(color = Origin), size = 3, offset = 0.01) +         # Color tip labels by origin
  geom_tippoint(aes(color = Origin), size = 2) +       # Add colored points at tips
  scale_color_manual(values = c("Nat-ITA" = "#006837", 
                                "Nat-FRA" = "#A50026", 
                                "Int-ITA" = "#66BD63", 
                                "Int-FRA" = "#F46D43")) +  # Adjust colors as needed
  theme(legend.position = "right")                    
p
