library(dplyr)
library(ggplot2)

K.annual <- read.csv("min.known.alive.forPVA.2Apr2025.csv") %>% 
  select(-X)
K.annual_sum <- K.annual 
K.annual_sum$Min_known <-round(rowSums(K.annual[,2:7], na.rm = T))

ggplot(K.annual_sum, aes(x=Year, y=Min_known)) +
  geom_line(color= "grey") +
  geom_point(shape=21, color="black", fill="#69b3a2", size=4) +
  labs(
    x = "Year",
    y = "Minimum number known alive across AUs",
    #title = "Minimum number known alive for the snail kite population used to calibrate carrying
#capacity."
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5),
    strip.text = element_text(face = "bold")
  )
