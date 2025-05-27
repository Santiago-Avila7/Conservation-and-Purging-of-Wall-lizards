##### MAPPING THE SAMPLES ######
pacman::p_load(openxlsx,rnaturalearth,rnaturalearthdata,rnaturalearthhires,
               ggrepel,ggplot2,RColorBrewer,dplyr,patchwork)


# Read the Excel file - Get one sample per location ----
LizardsData <- readWorkbook("Data/Samples_Santiago.xlsx", sheet = 1)
UKLizards <- subset(LizardsData, NativeVSIntro=="Intro")
UKLizards <- UKLizards %>%
  group_by(Longitude, Latitude) %>%
  slice(1)
NativeLizards <- subset(LizardsData, NativeVSIntro=="Native")

# Create the maps ----
#UK
UK<- ne_countries(country = "united kingdom", returnclass = "sf", scale="large")

map_UK <- ggplot(data = UK) +
  geom_sf(fill = "white", color = "black", linewidth = 0.8) +
  geom_point(data = UKLizards, aes(x = Longitude, y = Latitude, colour = Origin),
             size = 5, shape = 20) +
  coord_sf(xlim = c(-6.0, 2.0), ylim = c(50.0, 53.0), expand = FALSE) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_colour_manual(values = c("FRA" = "#F46D43", "ITA" = "#66BD63")) +
  theme(panel.background = element_rect(fill = "#BFD4FF"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "bottom") +
  geom_text_repel(data = UKLizards, aes(label = Abbpop, x = Longitude, y = Latitude),
                  box.padding = 0.7,
                  point.padding = 0,
                  size = 5,
                  segment.color = 'grey50') +
  labs(subtitle = "A")

#France + Italy
world <- ne_countries(scale = "medium", returnclass = "sf")
france_italy <- world %>% filter(admin %in% c("France", "Italy","Spain","Germany","Austria","Switzerland","Croatia"))

map_native <- ggplot(data = france_italy) +
  geom_sf(fill = "white", color = "black", linewidth = 0.8) +
  geom_point(data = NativeLizards, aes(x = Longitude, y = Latitude, colour = Origin),
             size = 5, shape = 20) +
  coord_sf(xlim = c(-5, 14), ylim = c(41, 49), expand = FALSE) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_colour_manual(values = c("FRA" = "#A50026", "ITA" = "#006837")) +
  theme(panel.background = element_rect(fill = "#BFD4FF"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "bottom") +
  geom_text_repel(data = NativeLizards, aes(label = Abbpop, x = Longitude, y = Latitude),
                  box.padding = 0.7,
                  point.padding = 0,
                  size = 5,
                  segment.color = 'grey50') +
  labs(subtitle = "B")

# Combine the maps and display them together
maps <- map_UK + map_native + plot_layout(ncol = 2)
print(maps)


