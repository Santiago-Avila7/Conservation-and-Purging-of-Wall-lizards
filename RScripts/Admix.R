## See ADMIXTURE data 
pacman::p_load(ggplot2,openxlsx,RColorBrewer,dplyr,gridExtra,reshape2,tidyr,paletteer,stringr)

#CV error plots
cv_error <- read.table("Data/PopGen/CV_errors_summary.txt", h=F)
colnames (cv_error)<- c("rm","rm2","K","CV_error")
cv_error <- cv_error %>%
  mutate(K = str_extract(K, "\\d+"))
cv_error <- cv_error %>% select(-rm, -rm2)
plot(cv_error)

#Load Admix data and merge with population data 

Lizards_admix<- readWorkbook("Data/Samples_Santiago.xlsx", sheet = 3)
K2<- read.table("Data/PopGen/K2.Q", header=F)
K3<- read.table("Data/PopGen/K3.Q", header=F)
K4<- read.table("Data/PopGen/K4.Q", header=F)

K2<- cbind(Lizards_admix, K2)
colnames(K2) <- c("ID","Origin","Abbpop","Q1","Q2")
K3<- cbind(Lizards_admix, K3)
colnames(K3) <- c("ID","Origin","Abbpop","Q1","Q2","Q3")
K4<- cbind(Lizards_admix, K4)
colnames(K4) <- c("ID","Origin","Abbpop","Q1","Q2","Q3","Q4")

K2_long <- K2 %>% pivot_longer(cols = starts_with("Q"), names_to = "Ancestry", values_to = "Q_value")
K3_long <- K3 %>% pivot_longer(cols = starts_with("Q"), names_to = "Ancestry", values_to = "Q_value")
K4_long <- K4 %>% pivot_longer(cols = starts_with("Q"), names_to = "Ancestry", values_to = "Q_value")

# Add a K identifier column
K2_long$K <- "K2"
K3_long$K <- "K3"
K4_long$K <- "K4"

# Combine all datasets into one long dataframe
all_data <- bind_rows(K2_long, K3_long, K4_long)

# Order individuals based on a chosen Q value (e.g., average Q1)
ordered_individuals <- all_data %>%
  filter(Ancestry == "Q1") %>%
  group_by(ID) %>%
  summarize(mean_Q1 = mean(Q_value)) %>%
  arrange(-mean_Q1) %>%
  pull(ID)

# Convert ID to a factor to control the order in the plot
all_data$ID <- factor(all_data$ID, levels = ordered_individuals)

# Plot
ggplot(all_data, aes(x = ID, y = Q_value, fill = Ancestry)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_wrap(~K, ncol = 1) +
  scale_fill_manual(values = as.vector(paletteer_d("ggsci::default_jco"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Admixture Proportion", x = "Individual")









