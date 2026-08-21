### Script for plotting of annotated hits ### 

library(tidyverse)
library(xcms)
library(MSnbase)

# get data on all time points 
AB007_files = list.files(path = "input_folder", full.names = TRUE, recursive = TRUE,
                         pattern = "Pool29|Pool30|Pool31|Pool32|HBSS")
source("./EIC_all_annot_function.R")

# select the molecular features of interest
df_eic = read.csv(file.path("input_folder", "AB007_baso_ttest_candidates.csv"))

library(BiocParallel)
bp <- SerialParam()
register(bp)

#### fit EICs ####
# plots theme 
cell_plot_theme <- theme(
  text = element_text(family = "Arial", size = 10), # Set Arial font and text size
  axis.text = element_text(size = 8), # Adjust axis text size
  axis.title = element_text(size = 10), # Adjust axis title size
  legend.text = element_text(size = 8), # Adjust legend text size
  legend.title = element_text(size = 8), # Adjust legend title size
  plot.title = element_text(size = 12, hjust = 0.5), # Center title
  panel.background = element_blank(), # Remove background color
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Add border
  panel.grid.major = element_blank(), # Subtle grid
  panel.grid.minor = element_blank(), # No minor grid
  plot.margin = margin(5, 5, 5, 5) # Adjust margins
)

## FT7682 Vanillic acid ####

df_temp = df_eic %>% filter(ion == "FT7682")

# first third 
candidates_df = EIC_all_mzmine(df_drugs = df_temp[1,], 
                               file_list = AB007_files, 
                               output_folder = "./EIC_output/")

write.csv(candidates_df, "./FT7682_ion.csv")
candidates_df = read.csv(file.path("input_folder", "FT7682_ion.csv"))

candidates_df$sample_name[28]
# add metadata 
candidates_df$side = sapply(strsplit(as.character(candidates_df$sample_name), "_"), '[', 7)
candidates_df$time = sapply(strsplit(as.character(candidates_df$sample_name), "_"), '[', 4)
candidates_df$pool = sapply(strsplit(as.character(candidates_df$sample_name), "_"), '[', 6)
candidates_df$strain = sapply(strsplit(as.character(candidates_df$sample_name), "_"), '[', 2)
candidates_df$bac_time = sapply(strsplit(as.character(candidates_df$sample_name), "_"), '[', 3)
candidates_df$time = as.numeric(gsub("T", "", candidates_df$time))

unique(candidates_df$ion)

# clean-up df 
# get rid of double peaks 
temp = candidates_df %>% 
  dplyr::group_by(sample_name, ion) %>% 
  dplyr::count(sample_name, sort=T) %>% 
  filter(n >1)
length(unique(temp$sample_name))
hist(temp$n)
hist(candidates_df$sn)

# use selection based on higher signal to noise ratio
fil_df = candidates_df %>% 
  dplyr::group_by(sample_name, ion, column) %>% 
  filter(sn == max(sn))

strain_labs = c("Cramosum" = "C. ramosum", 
                "Dformi" = "D. formicigenerans", 
                "Ghaemo" = "G. haemolysans")

p = candidates_df %>% 
  filter(strain %in% c("Cramosum", "Dformi", "Ghaemo")) %>% 
  filter(bac_time == "TP8") %>% 
  ggplot(., aes(x=time, y=intb, color=side)) + 
  geom_smooth(se = F, method = "loess", linewidth = 0.9) +
  geom_point() + 
  labs(x = "Time (minutes)", 
       y = "Peak area Intensity") + 
  theme_minimal() + 
  #facet_wrap(~condition+drug+strain) + 
  facet_grid(side~strain, scales = "free_y", 
             labeller = labeller(strain = strain_labs)) + 
  scale_color_manual(values = c("api" = "lightblue", "baso" = "lightpink")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE), n.breaks = 4) +
  theme(strip.text = element_text(face = "italic")) +
  cell_plot_theme +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm")
  )

print(p)
ggsave("./Plots/FT7682_vanillic_acid.tiff", plot = p, dpi = 300, width = 3.5, height = 2.5, units = "in")



## FT5225 5−hydroxynorvaline ####

df_temp = df_eic %>% filter(ion == "FT5225")

# first third 
candidates_df = EIC_all_mzmine(df_drugs = df_temp[1,], 
                               file_list = AB007_files, 
                               output_folder = "./EIC_output/")

write.csv(candidates_df, "./FT5225_ion.csv")
candidates_df = read.csv(file.path("input_folder", "FT5225_ion.csv"))

candidates_df$sample_name[28]
# add metadata 
candidates_df$side = sapply(strsplit(as.character(candidates_df$sample_name), "_"), '[', 7)
candidates_df$time = sapply(strsplit(as.character(candidates_df$sample_name), "_"), '[', 4)
candidates_df$pool = sapply(strsplit(as.character(candidates_df$sample_name), "_"), '[', 6)
candidates_df$strain = sapply(strsplit(as.character(candidates_df$sample_name), "_"), '[', 2)
candidates_df$bac_time = sapply(strsplit(as.character(candidates_df$sample_name), "_"), '[', 3)
candidates_df$time = as.numeric(gsub("T", "", candidates_df$time))

unique(candidates_df$ion)

# clean-up df 
# get rid of double peaks 
temp = candidates_df %>% 
  dplyr::group_by(sample_name, ion) %>% 
  dplyr::count(sample_name, sort=T) %>% 
  filter(n >1)
length(unique(temp$sample_name))
hist(temp$n)
hist(candidates_df$sn)

# use selection based on higher signal to noise ratio
fil_df = candidates_df %>% 
  dplyr::group_by(sample_name, ion, column) %>% 
  filter(sn == max(sn))

p = candidates_df %>% 
  #filter(strain %in% c("Cramosum", "Dformi", "Ghaemo")) %>% 
  filter(bac_time == "TP8") %>% 
  ggplot(., aes(x=time, y=intb, color=side)) + 
  geom_smooth(se = F, method = "loess", linewidth = 0.9) +
  geom_point() + 
  labs(x = "Time (minutes)", 
       y = "Peak area Intensity") + 
  theme_minimal() + 
  #facet_wrap(~condition+drug+strain) + 
  facet_grid(side~strain, scales = "free_y", 
             labeller = labeller(strain = strain_labs)) + 
  scale_color_manual(values = c("api" = "lightblue", "baso" = "lightpink")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(strip.text = element_text(face = "italic")) +
  cell_plot_theme

print(p)


## api hits####

## hits: FT11512 Rimboxo, FT15725 Dethiobiotin, FT2029 Cyclo(D−Arg−L−Pro)
# FT14576 Decarestrictine D, FT2029 Cyclo(D−Arg−L−Pro), FT3430 6−hydroxytryprostatin

df_eic = read.csv(file.path("input_folder", "AB007_api_ttest_candidates.csv"))
df_temp = df_eic %>% 
  filter(ion %in% c("FT11512", "FT15725", "FT2029", "FT14576", "FT2029", "FT3430")) %>% 
  dplyr::select(ion, row.ID, row.m.z, row.retention.time) %>% 
  distinct()

# fit ions
candidates_df = EIC_all_mzmine(df_drugs = df_temp, 
                               file_list = AB007_files, 
                               output_folder = "./EIC_output/")

write.csv(candidates_df, "./api_ions.csv")
#candidates_df = read.csv("FT5225_ion.csv")

candidates_df$sample_name[28]
# add metadata 
candidates_df$side = sapply(strsplit(as.character(candidates_df$sample_name), "_"), '[', 7)
candidates_df$time = sapply(strsplit(as.character(candidates_df$sample_name), "_"), '[', 4)
candidates_df$pool = sapply(strsplit(as.character(candidates_df$sample_name), "_"), '[', 6)
candidates_df$strain = sapply(strsplit(as.character(candidates_df$sample_name), "_"), '[', 2)
candidates_df$bac_time = sapply(strsplit(as.character(candidates_df$sample_name), "_"), '[', 3)
candidates_df$time = as.numeric(gsub("T", "", candidates_df$time))

unique(candidates_df$ion)

# clean-up df 
# get rid of double peaks 
temp = candidates_df %>% 
  dplyr::group_by(sample_name, ion) %>% 
  dplyr::count(sample_name, sort=T) %>% 
  filter(n >1)
length(unique(temp$sample_name))
hist(temp$n)
hist(candidates_df$sn)

# use selection based on higher signal to noise ratio
fil_df = candidates_df %>% 
  dplyr::group_by(sample_name, ion, column) %>% 
  filter(sn == max(sn))
unique(fil_df$ion)

strain_labs = c("Anaes" = "A. naeslundi", 
                "Ef2" = "E. faecalis", 
                "Ghaemo" = "G. haemolysans")

p = fil_df %>% 
  filter(ion == "FT2029") %>% 
  filter(strain %in% c("Ghaemo","Anaes","Ef2")) %>% 
  filter(bac_time == "TP8") %>% 
  ggplot(., aes(x=time, y=intb, color=side)) + 
  geom_smooth(se = F, method = "loess", linewidth = 0.9) +
  geom_point() + 
  labs(x = "Time (minutes)", 
       y = "Peak area Intensity") + 
  theme_minimal() + 
  #facet_wrap(~condition+drug+strain) + 
  facet_grid(side~strain, scales = "free_y", 
             labeller = labeller(strain = strain_labs)) + 
  scale_color_manual(values = c("api" = "lightblue", "baso" = "lightpink")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE), n.breaks = 4) +
  theme(strip.text = element_text(face = "italic")) +
  cell_plot_theme +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm")
  )

print(p)
ggsave("./Plots/FT2029_cycloArgPro.tiff", plot = p, dpi = 300, width = 3.5, height = 2.5, units = "in")
