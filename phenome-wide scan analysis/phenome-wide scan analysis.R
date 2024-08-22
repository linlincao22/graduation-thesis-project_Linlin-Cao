library(ggplot2)
library(dplyr)

# Read the data
miR_130a_3p_atlas <- read.table("130a_3p.csv", sep = ",", header = TRUE)

data <- miR_130a_3p_atlas %>%
  filter(!is.na(Domain) & !is.na(Trait) & !is.na(P.value) & !is.na(N)) %>%
  mutate(P = as.numeric(P.value), N = as.numeric(N)) %>%
  mutate(logP = -log10(P)) %>%
  group_by(Domain) %>%
  mutate(trait_position = row_number()) %>%
  ungroup() %>%
  mutate(Domain_factor = factor(Domain, levels = unique(Domain)))


# Plot with adjusted positioning of points to avoid overlapping
miR_130a_3p_plot <- ggplot(data, aes(x = Domain_factor, y = logP, color = Domain, size = N)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5,  jitter.height = 0.02, dodge.width = 0.7), alpha = 0.9) +
  scale_size(range = c(1, 10), labels = scales::scientific) +
  labs(title = "PheWAS Plot for rs685149",
       y = "-log10(P)",
       x = "Domain",
       size = "Sample Size") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 20),  # Adjust font size for x-axis text
    axis.text.y = element_text(size = 14),  # Adjust font size for y-axis text
    axis.title.y = element_text(size = 20),  # Adjust font size for y-axis title
    plot.title = element_text(size = 20),  # Adjust font size for title
    legend.title = element_text(size = 20),  # Adjust font size for legend title
    legend.text = element_text(size = 20),  # Adjust font size for legend text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    panel.grid = element_blank()    # Remove gridlines
  ) +
  theme(text=element_text(size=16,  family="Times New Roman")) +
  scale_y_continuous(limits = c(0,15), breaks = c(2.5, 5, 7.5, 10, 12.5, 15) ) +
  geom_hline(yintercept = 7.3, linetype = "dashed", color = "black") +
  guides(
    size = guide_legend(
      override.aes = list(shape = 16, color = "black", fill = "white", stroke = 0),  # Black dots, white background, no border
      title = "Sample Size"
    ),
    color = guide_legend(
      override.aes = list(size = 5)  # Increase size of dots in Domain legend
    )
  )

miR_130a_3p_plot
ggsave("phewas plot revised y axis/rs685149_version1.jpeg", miR_130a_3p_plot, width = 12, height = 5, dpi = 2000)



# Read the data
miR_132_3p_atlas <- read.table("132_3p.csv", sep = ",", header = TRUE)

data <- miR_132_3p_atlas %>%
  filter(!is.na(Domain) & !is.na(Trait) & !is.na(P.value) & !is.na(N)) %>%
  mutate(P = as.numeric(P.value), N = as.numeric(N)) %>%
  mutate(logP = -log10(P)) %>%
  group_by(Domain) %>%
  mutate(trait_position = row_number()) %>%
  ungroup() %>%
  mutate(Domain_factor = factor(Domain, levels = unique(Domain)))


# Plot with adjusted positioning of points to avoid overlapping
miR_132_3p_plot <- ggplot(data, aes(x = Domain_factor, y = logP, color = Domain, size = N)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.6,  jitter.height = 0.01, dodge.width = 0.7), alpha = 0.9) +
  scale_size(range = c(1, 10), labels = scales::scientific) +
  labs(title = "PheWAS Plot for rs3803809",
       y = "-log10(P)",
       x = "Domain",
       size = "Sample Size") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 20),  # Adjust font size for x-axis text
    axis.text.y = element_text(size = 14),  # Adjust font size for y-axis text
    axis.title.y = element_text(size = 20),  # Adjust font size for y-axis title
    plot.title = element_text(size = 20),  # Adjust font size for title
    legend.title = element_text(size = 20),  # Adjust font size for legend title
    legend.text = element_text(size = 20),  # Adjust font size for legend text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    panel.grid = element_blank() 
  ) +
  theme(text=element_text(size=16,  family="Times New Roman")) +
  scale_y_continuous(limits = c(0,60), breaks = c(5,10,15,20,25,30,35,40,45,50,55,60) ) +
  geom_hline(yintercept = 7.3, linetype = "dashed", color = "black") +
  guides(
    size = guide_legend(
      override.aes = list(shape = 16, color = "black", fill = "white", stroke = 0),  # Black dots, white background, no border
      title = "Sample Size"
    ),
    color = guide_legend(
      override.aes = list(size = 3)  # Increase size of dots in Domain legend
    )
  )

miR_132_3p_plot
ggsave("phewas plot revised y axis/rs3803809_version1.jpeg", miR_132_3p_plot, width = 12, height = 5, dpi = 2000)


View(data)






# Read the data
miR_133a_3p_atlas <- read.table("133a_3p.csv", sep = ",", header = TRUE)

data <- miR_133a_3p_atlas %>%
  filter(!is.na(Domain) & !is.na(Trait) & !is.na(P.value) & !is.na(N)) %>%
  mutate(P = as.numeric(P.value), N = as.numeric(N)) %>%
  mutate(logP = -log10(P)) %>%
  group_by(Domain) %>%
  mutate(trait_position = row_number()) %>%
  ungroup() %>%
  mutate(Domain_factor = factor(Domain, levels = unique(Domain)))


# Plot with adjusted positioning of points to avoid overlapping
miR_133a_3p_plot <- ggplot(data, aes(x = Domain_factor, y = logP, color = Domain, size = N)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5,  jitter.height = 0.01, dodge.width = 0.7), alpha = 0.9) +
  scale_size(range = c(1, 10), labels = scales::scientific) +
  labs(title = "PheWAS Plot for rs12479469",
       y = "-log10(P)",
       x = "Domain",
       size = "Sample Size") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 20),  # Adjust font size for x-axis text
    axis.text.y = element_text(size = 14),  # Adjust font size for y-axis text
    axis.title.y = element_text(size = 20),  # Adjust font size for y-axis title
    plot.title = element_text(size = 20),  # Adjust font size for title
    legend.title = element_text(size = 20),  # Adjust font size for legend title
    legend.text = element_text(size = 20),  # Adjust font size for legend text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    panel.grid = element_blank()    # Remove gridlines 
  ) +
  theme(text=element_text(size=16,  family="Times New Roman")) +
  scale_y_continuous(limits = c(0,15), breaks = c(2.5,5,7.5,10,12.5,15) ) +
  geom_hline(yintercept = 7.3, linetype = "dashed", color = "black") +
  guides(
    size = guide_legend(
      override.aes = list(shape = 16, color = "black", fill = "white", stroke = 0),  # Black dots, white background, no border
      title = "Sample Size"
    ),
    color = guide_legend(
      override.aes = list(size = 3)  # Increase size of dots in Domain legend
    )
  ) 

miR_133a_3p_plot
ggsave("phewas plot revised y axis/rs12479469_version1.jpeg", miR_133a_3p_plot, width = 12, height = 5, dpi = 2000)

View(data)





