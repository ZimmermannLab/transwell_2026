### Script for viability & TEER data visualization ###


# libraries 
library(tidyverse)
library(readxl)

# plotting theme 
cell_plot_theme <- theme(
  text = element_text(family = "Arial", size = 10), # Set Arial font and text size
  axis.text = element_text(size = 6), # Adjust axis text size
  axis.title = element_text(size = 8), # Adjust axis title size
  legend.text = element_text(size = 6), # Adjust legend text size
  legend.title = element_text(size = 6), # Adjust legend title size
  plot.title = element_text(size = 12, hjust = 0.5), # Center title
  panel.background = element_blank(), # Remove background color
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Add border
  panel.grid.major = element_blank(), # Subtle grid
  panel.grid.minor = element_blank(), # No minor grid
  plot.margin = margin(5, 5, 5, 5) # Adjust margins
)

#### Viability data ####
df_viab = read.table(file.path("input_folder", "AB012I_cellviab_24well.txt"), skip = 39, sep = ",", header = T)
df_viab = na.omit(df_viab)

# calculate the relative viability to GMM 
df_viab = df_viab %>% 
  mutate(gmm_value = mean(`X560.590`[group == "GMM"])) %>% 
  mutate(rel_viab = `X560.590` / gmm_value)

p = ggplot(df_viab, aes(x=group, y=rel_viab, fill = group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0), size = 0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Fluoresence 560/590") + 
  theme(text = element_text(size = 20)) + 
  theme_minimal() +
  ylab("Relative signal 560/590") +
  theme(axis.title.x = element_blank()) +
  #scale_x_discrete(labels = strain_lab) + 
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#DDA0DD", "#5F9EA0", "#E75480","#F08080", "#FFD1DC", "#FFA07A")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, 
                                   face = "italic")) + 
  cell_plot_theme
print(p)

ggsave("./Plots/AB012I_viability.tiff", plot = p, dpi = 600, width = 3.5, height = 2.5, units = "in")
ggsave("./Plots/AB012I_viability.png", plot = p, dpi = 600, width = 3.5, height = 2.5, units = "in")

### TEER data ####
library(xlsx)
df_teer = read_excel(file.path("input_folder", "AB012I_TEER_24well.xlsx"))

ggplot(df_teer, aes(x=measurement, y=TEER)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  facet_wrap(~group)
  
df_teer$measurement <- factor(df_teer$measurement, 
                              levels = c("before", "after_4H", "after_8H"))


df_teer %>% dplyr::filter(measurement != "after_8H") %>% 
  ggplot(., aes(x=measurement, y=TEER, fill=measurement)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  geom_jitter(width = 0.2) + 
  ylab("Electrical resistance (Ohm)") +
  facet_wrap(~group, nrow = 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(text = element_text(size = 20)) + 
  scale_fill_manual(labels = c("before" = "0 hours", 
                               "after_4H" = "4 hours"), 
                    values = c("before" = "#FFC0CB", "after_4H" = "#5F9EA0")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(
    labels = c("before" = "0 hours", "after_4H" = "4 hours")) +
  cell_plot_theme

ggsave("./Plots/AB012I_TEER.tiff", dpi = 600, width = 4.5, height = 2, units = "in")
ggsave("./Plots/AB012I_TEER.png", dpi = 600, width = 4.5, height = 2.5, units = "in")

df_teer %>% 
  ggplot(., aes(x=measurement, y=TEER, fill=measurement)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  geom_jitter(width = 0.2) + 
  ylab("Electrical resistance (Ohm)") +
  facet_wrap(~group, nrow = 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(text = element_text(size = 20)) + 
  scale_fill_manual(labels = c("before" = "0 hours", 
                               "after_4H" = "4 hours", 
                               "after_8H" = "8 hours (4 hours recovery)"), 
                    values = c("before" = "#FFC0CB", "after_4H" = "#5F9EA0", "after_8H" = "orange")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(
    labels = c("before" = "0 hours", "after_4H" = "4 hours", "after_8H" = "after recovery")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1)) +
  cell_plot_theme

ggsave("./Plots/AB012I_TEER_w8hours.tiff", dpi = 600, width = 4.5, height = 2.5, units = "in")



