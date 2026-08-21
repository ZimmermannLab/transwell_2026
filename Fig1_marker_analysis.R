#### Script for AB012G marker analysis ####

## Set up 
####Load Libraries#####
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(ggforce)

#### Data import and cleaning ####
data_list = file.path(
  "input_folder",
  c("markers_AB012G1_AB012G5_AB012G6.csv", "markers_AB012G2_AB012G3.csv")
)

#read them
screen <- lapply(data_list, read.csv)


# run the following data cleaning functions for the screen list 
# extract the targeted drugs names (+ remove spaces)
targeted_drugs = function(df){
  gsub(" ", "_", gsub(" Results","",grep("Results",colnames(df),fixed=T,value=T),fixed = T), fixed = T)
}

target_drugs = lapply(screen, targeted_drugs)

# remove unnecessary columns 
remove_cols = function(df){
  
  if (any(df[1] == "!") == TRUE){ df = df[, -c(1,2)] }
  else {df = df}
}

screen = lapply(screen, remove_cols)



# update column names
colnms <- list()
for (i in seq_along(target_drugs)) { 
  colnms[[i]] <- as.character(c("Name", "Data_file", "Type", "level", "timestamp", 
                                as.vector(sapply(target_drugs[[i]], function(x) paste(x,c("_RT","_Area"),sep="")))))
}


for (i in seq_along(screen)){
  colnames(screen[[i]]) <- colnms[[i]]
}

# remove first row 
# set empty to NA
cleaning_funct = function(df){
  df = df %>% 
    slice(-1) %>% 
    mutate_all(list(~na_if(.,"")))
}

screen = lapply(screen, cleaning_funct)

# merge the dataframes in the screen list together
merged_df = bind_rows(screen)

# Add columns with timepoints & conditions 
merged_df$strain = sapply(strsplit(as.character(merged_df$Name), "_"), '[', 2)
merged_df$TP = sapply(strsplit(as.character(merged_df$Name), "_"), '[', 3)
merged_df$rep = sapply(strsplit(as.character(merged_df$Name), "_"), '[', 4)
merged_df$side = sapply(strsplit(as.character(merged_df$Name), "_"), '[', 5)
merged_df$exp = sapply(strsplit(as.character(merged_df$Name), "_"), '[', 1)

# Change area variables to numeric
cols = grep("_Area", colnames(merged_df))
merged_df[, cols] = apply(merged_df[,cols], 2, function(x) as.numeric(as.character(x)))

# Change TPs to numeric 
merged_df$time = as.numeric(gsub("T", "", merged_df$TP))


# delete waterwashes and IS samples 
merged_df = merged_df[!grepl("waterwash|IS|Water", merged_df$Name),]

groups_exp = c("group_24H" = c("AB012G1", "AB012G5"), 
               "group_16H" = c("AB012G2", "AB012G3", "AB012G6"))

merged_df = merged_df %>% 
  mutate(exp_group = ifelse(exp %in% c("AB012G1", "AB012G5"), "group_24H", "group_16H"))


cols = grep("marker.*_Area", colnames(merged_df))

# plot raw data 
pdf("AB012Gall_raw_markers.pdf")
for (i in cols){
  
  y_var = sym(colnames(merged_df)[i])
  
  # Default line plot
  p <-  ggplot(merged_df, aes(x=time, y=!!y_var, color = side, shape = rep, 
                              group = interaction(side, rep, strain))) +
    geom_point() +
    geom_line(aes(group=interaction(rep, side, strain, exp))) +
    labs(title= colnames(merged_df)[i], 
         subtitle="Timecourse", 
         x = "Time (minutes)", 
         y = "Peak area Intensity (raw)") + 
    theme_minimal() + 
    facet_wrap(~strain+exp_group)
  
  print(p)
  
}
dev.off()



# define colors for same plotting in different groups 
tar_colors = c("An" = "deeppink",
               "Ao" = "darkblue",
               "Bt" = "skyblue", 
               "Bu" = "brown", 
               "Cc" = "darkolivegreen2", 
               "Cr" = "darkorange", 
               "Cs" = "gold2", 
               "Df" = "darkorchid", 
               "Ef" = "darkgreen", 
               "Gh" = "chartreuse", 
               "GMM" = "red"
)

cb_palette <- c(
  "#FFB6C1", # Light Pink
  "#F08080",# Light Coral
  "#ADD8E6", # Light Blue
  "#E75480", # Dark Pink
  "#FFD1DC", # Pastel Pink
  "#A2C2E0", # Pale Blue
  "#FF69B4", # Pink
  "#B0E0E6", # Powder Blue
  "#FF7F7F", # Salmon Pink
  "#87CEEB", # Sky Blue
  "#DDA0DD", # Plum
  "#AFEEEE", # Pale Turquoise
  "#FFC0CB", # Classic Pink
  "#5F9EA0"  # Cadet Blue
)


### plot baso side #### 
pdf("AB012Gall_baso_markers.pdf")
for (i in cols){
  
  y_var = sym(colnames(merged_df)[i])
  
  # Default line plot
  tmp = merged_df %>% filter(side == "baso")
  p <-  ggplot(tmp, aes(x=time, y=!!y_var, color = strain, shape = rep, 
                              group = interaction(side, rep, strain, exp))) +
    geom_point() +
    geom_line(aes(group=interaction(rep, side, strain, exp))) +
    labs(title= colnames(merged_df)[i], 
         subtitle="Timecourse", 
         x = "Time (minutes)", 
         y = "Peak area Intensity (raw)") + 
    theme_minimal() + 
    scale_color_manual(values = tar_colors) +
    facet_wrap(~exp_group)
  
  print(p)
  
}
dev.off()


#### Scale to TP0 api #### 
scaled = merged_df %>% 
  dplyr::select(Name, strain, exp, exp_group, time, rep, side, all_of(grep("marker.*_Area", colnames(merged_df))))

scaled = pivot_longer(scaled, cols = 8:14)

scaled$permeable = ifelse(scaled$name %in% non_per, "non_permeable", "permeable")

to_divide_by = scaled %>% 
  filter(side == "api") %>% 
  group_by(name, strain, exp_group) %>% 
  slice(which.min(time))

# calculate mean of api TP0 for all 4 reps 
to_divide_by = rename(to_divide_by, api_t00 = value)
to_divide_by = to_divide_by %>% 
  group_by(name, strain, exp_group) %>% 
  mutate(mean_api_t00 = mean(api_t00, na.rm = T)) %>% 
  select(name, strain, exp, exp_group, mean_api_t00) %>% 
  distinct()

scaled = left_join(scaled, to_divide_by)

scaled = scaled %>% 
  mutate(Scaled_Area = value/mean_api_t00)

### plot markers baso only for non-permeable ####
non_per = c("marker_nadolol.Results_Area", "marker_terbutaline.Results_Area", "marker_etoposide.Results_Area")

tmp = merged_df %>% filter(side == "baso") %>% 
  dplyr::select(Name, strain, exp, exp_group, time, rep, all_of(non_per))

tmp = pivot_longer(tmp, cols = 7:9)

# filter out crazy peaks 
tmp = tmp %>% 
  filter(value < 3.5e5)

tmp = tmp %>% 
  filter(!(name == "marker_etoposide.Results_Area" & value > 3.5e4))

# change names 
tmp$name = as.factor(tmp$name)
levels(tmp$name)
levels(tmp$name) = c("Etoposide", "Nadolol", "Terbutaline")

p <-  ggplot(tmp, aes(x=time, y=value, color = strain, 
                      group = interaction(rep, strain, exp))) +
  geom_point() +
  geom_line(aes(group=interaction(rep, strain, exp))) +
  labs(title= "Minimally permeable marker compounds", 
       subtitle="Timecourse", 
       x = "Time (minutes)", 
       y = "Peak area Intensity (raw)") + 
  theme_minimal() + 
  scale_color_manual(values = tar_colors) +
  facet_wrap(~exp_group+name)

print(p)


### Figure for GMM ####
cell_plot_theme <- theme(
  text = element_text(family = "Arial", size = 10), # Set Arial font and text size
  axis.text = element_text(size = 10), # Adjust axis text size
  axis.title = element_text(size = 12), # Adjust axis title size
  legend.text = element_text(size = 10), # Adjust legend text size
  legend.title = element_text(size = 10), # Adjust legend title size
  plot.title = element_text(size = 12, hjust = 0.5), # Center title
  panel.background = element_blank(), # Remove background color
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Add border
  panel.grid.major = element_blank(), # Subtle grid
  panel.grid.minor = element_blank(), # No minor grid
  plot.margin = margin(5, 5, 5, 5) # Adjust margins
)

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

tmp = merged_df %>% 
  filter(strain == "GMM") %>% 
  dplyr::select(Name, strain, exp, exp_group, time, rep, side, all_of(grep("marker.*_Area", colnames(merged_df))))

tmp = pivot_longer(tmp, cols = 8:14)

#tmp$permeable = ifelse(tmp$name %in% non_per, "non_permeable", "permeable")

fig1 = tmp %>% 
 filter(name %in% c("marker_antipyrine.Results_Area", "marker_nadolol.Results_Area"))

drug_labs = c("marker_antipyrine.Results_Area" = "Antipyrine", 
              "marker_nadolol.Results_Area" = "Nadolol")

p <-  ggplot(fig1, aes(x=time, y=value, color = side)) +
  geom_point(aes(shape = name), show.legend=F, size=0.7) +
  #geom_line(aes(group=interaction(rep, side, exp))) +
  geom_smooth(se = F, method = "auto", linewidth = 0.8, aes(group = interaction(side, name))) +
  labs(x = "Time (minutes)", 
       y = "Peak Area Intensity") + 
  theme_minimal() + 
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("api" = "lightblue", "baso" = "lightpink"),
    labels = c("api" = "Apical", "baso" = "Basolateral")
  ) +
  #scale_color_manual(values = tar_colors) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  facet_grid(name~exp_group, scales = "free_y", labeller = labeller(name = drug_labs)) +
  scale_y_continuous(labels = scales::scientific) +
  cell_plot_theme

print(p)
ggsave("AB012_gmm_markers_1d.tiff", plot = p, dpi = 300, width = 3, height = 3, units = "in")

### all for suppl. figure ###
drug_labs = c("marker_antipyrine.Results_Area" = "Antipyrine", 
              "marker_nadolol.Results_Area" = "Nadolol", 
              "marker_terbutaline.Results_Area" = "Terbutaline", 
              "marker_propanolol.Results_Area" = "Propanolol", 
              "marker_haloperidol.Results_Area" = "Haloperidol", 
              "marker_etoposide.Results_Area" = "Etoposide", 
              "marker_warfarin.Results_Area" = "Warfarin")


p <-  ggplot(tmp, aes(x=time, y=value, color = side)) +
  geom_point(show.legend=F, size=0.7) +
  #geom_line(aes(group=interaction(rep, side, exp))) +
  geom_smooth(se = F, method = "auto", linewidth = 0.8, aes(group = interaction(side, name))) +
  labs(x = "Time (minutes)", 
       y = "Peak Area Intensity") + 
  theme_minimal() + 
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.position = "top") +
  scale_color_manual(values = c("api" = "lightblue", "baso" = "lightpink")) +
  #scale_color_manual(values = tar_colors) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  facet_grid(name~exp_group, scales = "free_y", labeller = labeller(name = drug_labs)) +
  cell_plot_theme
print(p)
ggsave("AB012_gmm_markers_1d.tiff", plot = p, dpi = 300, width = 3, height = 2.5, units = "in")


### Nadolol figure ### 
tmp = merged_df %>% 
  dplyr::select(Name, strain, exp, exp_group, time, rep, side, all_of(grep("marker.*_Area", colnames(merged_df))))

tmp = pivot_longer(tmp, cols = 8:14)

fig1 = tmp %>% 
  filter(name %in% c("marker_nadolol.Results_Area"))

unique(fig1$exp)
fig1 = fig1 %>% 
  filter(exp %in% c('AB012G1', "AB012G2", "AB012G3"))

# filter out crazy peaks 
fig1 = fig1 %>% 
  filter(value < 3.5e5)

drug_labs = c("marker_terbutaline.Results_Area" = "Terbutaline", 
              "marker_nadolol.Results_Area" = "Nadolol")


p <-  fig1 %>% filter(side == "baso") %>%
  filter(strain %in% c("Cr", "Gh", "Cs", "Bt", "Cc", "GMM")) %>% 
  ggplot(., aes(x=time, y=value, color = strain)) +
  #geom_point(size=0.7) +
  #geom_line(aes(group=interaction(rep, strain, exp))) +
  geom_smooth(se = F, method = "loess", linewidth = 0.9) +
  #stat_summary(fun.data ='mean_se', geom = "smooth", se = TRUE, alpha=0.1) +
  labs(x = "Time (minutes)", 
       y = "Peak Area Intensity", 
       color = "Species") + 
  theme_minimal() + 
  #scale_color_manual(values = c("api" = "lightblue", "baso" = "lightpink")) +
  scale_color_manual(values = cb_palette, 
                     labels = c(
                       "Cr" = expression(italic("C. ramosum")),
                       "Gh" = expression(italic("G. haemolysans")),
                       "Cs" = expression(italic("C. scindens")),
                       "Bt" = expression(italic("B. thetaiotaomicron")),
                       "Cc" = expression(italic("C. comes")),
                       "GMM" = expression("media control")
                     )) + # Italicize legend labels) +
  theme(
    legend.text.align = 0) + # Align text to the left
  scale_x_continuous(breaks = c(0, 120, 240)) +
  theme(
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm")
  ) +
  facet_grid(name~exp_group, labeller = labeller(name = drug_labs, strain = strain_labs)) +
  cell_plot_theme

print(p)
ggsave("AB012_baso_markers_1e.tiff", plot = p, dpi = 300, width = 4, height = 2, units = "in")

### all for suppl fig ###
tmp = merged_df %>% 
  dplyr::select(Name, strain, exp, exp_group, time, rep, side, all_of(grep("marker.*_Area", colnames(merged_df))))

tmp = pivot_longer(tmp, cols = 8:14)

fig1 = tmp %>% 
  filter(name %in% c("marker_nadolol.Results_Area", 
                     "marker_terbutaline.Results_Area", 
                     "marker_etoposide.Results_Area"))

unique(fig1$exp)
fig1 = fig1 %>% 
  filter(exp %in% c('AB012G1', "AB012G2", "AB012G3"))

# filter out crazy peaks 
fig1 = fig1 %>% 
  filter(value < 3.5e5) %>% 
  filter(!(name == "marker_etoposide.Results_Area" & value > 2e4))

drug_labs = c("marker_terbutaline.Results_Area" = "Terbutaline", 
              "marker_nadolol.Results_Area" = "Nadolol", 
              "marker_etoposide.Results_Area" = "Etoposide")


p <-  fig1 %>% filter(side == "baso") %>%
  #filter(strain %in% c("Cr", "Gh", "Cs", "Bt", "Cc", "GMM")) %>% 
  ggplot(., aes(x=time, y=value, color = strain)) +
  #geom_point(size=0.7) +
  #geom_line(aes(group=interaction(rep, strain, exp))) +
  geom_smooth(se = F, method = "loess", linewidth = 0.9) +
  #stat_summary(fun.data ='mean_se', geom = "smooth", se = TRUE, alpha=0.1) +
  labs(x = "Time (minutes)", 
       y = "Peak Area Intensity", 
       color = "Species") + 
  theme_minimal() + 
  #scale_color_manual(values = c("api" = "lightblue", "baso" = "lightpink")) +
  scale_color_manual(values = cb_palette, 
                     labels = c(
                       "Cr" = expression(italic("C. ramosum")),
                       "Gh" = expression(italic("G. haemolysans")),
                       "Cs" = expression(italic("C. scindens")),
                       "Bt" = expression(italic("B. thetaiotaomicron")),
                       "Cc" = expression(italic("C. comes")),
                       "GMM" = expression("media control"), 
                       "An" = expression(italic("A. naeslundi")), 
                       "Ao" = expression(italic("A. omnicolens")), 
                       "Bu" = expression(italic("B. uniformis")), 
                       "Df" = expression(italic("D. formicigenerans")), 
                       "Ef" = expression(italic("E. faecalis"))
                     )) + # Italicize legend labels) +
  theme(
    legend.text.align = 0) + # Align text to the left
  scale_x_continuous(breaks = c(0, 120, 240)) +
  theme(
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm")
  ) +
  facet_grid(name~exp_group, labeller = labeller(name = drug_labs, strain = strain_labs), 
             scales = "free") +
  cell_plot_theme

print(p)
ggsave("AB012_baso_markers_suppl2.tiff", plot = p, dpi = 300, width = 4, height = 4, units = "in")
