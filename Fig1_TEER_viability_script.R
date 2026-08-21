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

strain_lab = c("An" = "A. naeslundi", 
               "Ao" = "A. omnicolens", 
               "Bt" = "B. thetaiotaomicron", 
               "Bu" = "B. uniformis", 
               "Cc" = "C. comes", 
               "Cr" = "C. ramosum", 
               "Cs" = "C. scindens", 
               "Df" = "D. formicigenerans",
               "Ef" = "E. faecalis", 
               "Gh" = "G. haemolysans", 
               "GMM" = "media control")

#### Viability data ####
df_viab = read.csv(file.path("input_folder", "Caco2_viability_data.csv"))

p = ggplot(df_viab, aes(x=condition, y=rel_viab, fill = bacterial_time)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = bacterial_time), position = position_jitterdodge(jitter.width = 0), size = 0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Relative fluo 560/590") + 
  theme(text = element_text(size = 20)) + 
  theme_minimal() +
  labs(fill = "Bacterial incubation time") + 
  theme(axis.title.x = element_blank()) +
  scale_fill_manual(values = c("16_hour" = "#FF7F7F", 
                               "24_hour" = "#E75480")) + 
  scale_x_discrete(labels = strain_lab) + 
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, 
                                   face = "italic")) + 
  cell_plot_theme
print(p)
### TEER data ####
df_teer = read.csv(file.path("input_folder", "Caco2_TEER.csv"))

# Create custom labels for the x-axis
custom_labels <- unlist(lapply(
  c("GMM", "An", "Ao", "Bt", "Bu", "Cc", "Cr", "Cs", "Df", "Ef", "Gh"),
  function(cond) {
    c(
      paste0(cond, "\n16H"),
      paste0(cond, "\n24H"),
      paste0(cond, "\n16H"),
      paste0(cond, "\n24H")
    )
  }
))

# plot all together 
p = df_teer %>% 
  filter(!(condition == "Cr" & name == "T0" & TEER < 1000)) %>% 
  mutate(
    # Create the interaction factor with explicit ordering
    interaction_factor = factor(
      paste(condition, name, exp_group, sep = "_"),
      levels = unlist(lapply(
        c("GMM", "An", "Ao", "Bt", "Bu", "Cc", "Cr", "Cs", "Df", "Ef", "Gh"),
        function(cond) {
          c(
            paste(cond, "T0", "group_16H", sep = "_"),
            paste(cond, "T0", "group_24H", sep = "_"),
            paste(cond, "T4_HBSS", "group_16H", sep = "_"),
            paste(cond, "T4_HBSS", "group_24H", sep = "_")
          )
        }
      ))
    ),
    # Create a grouping variable for conditions
    condition_factor = factor(condition, levels = c("GMM", "An", "Ao", "Bt", "Bu", "Cc", "Cr", "Cs", "Df", "Ef", "Gh"))
  ) %>% 
  ggplot(aes(x = interaction_factor, y = TEER, fill = name)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = name), position = position_jitterdodge(jitter.width = 0), size=0.7) +
  scale_fill_manual(labels = c("T0" = "0 hours", 
                               "T4_HBSS" = "4 hours"), 
                    values = c("T0" = "#FFC0CB", "T4_HBSS" = "#5F9EA0")) + 
  facet_grid(~condition_factor, scales = "free_x", space = "free_x") +
  scale_x_discrete(
    labels = function(x) {
      # Extract only the time part from the interaction factor
      sub(".*_(group_16H|group_24H)$", "\\1", x) %>%
        sub("group_", "", .) # Remove "group_" prefix
    }
  ) +
  labs(fill = "Exposure time of \nCaco2 cells to \nbacterial supernatant", 
       x = "Bacterial species and cultivation time until supernatant collection") + 
  theme_minimal() +
  theme(legend.key.size = unit(0.3, "cm")) +
  cell_plot_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.8))

print(p)
