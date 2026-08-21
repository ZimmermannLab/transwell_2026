#### Script to analyse transwell + drugs experiment ### 

## Set up 
####Load Libraries#####
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(ggforce)

cell_plot_theme <- theme(
  text = element_text(family = "Arial", size = 8), # Set Arial font and text size
  axis.text = element_text(size = 6), # Adjust axis text size
  axis.title = element_text(size = 6), # Adjust axis title size
  legend.text = element_text(size = 6), # Adjust legend text size
  legend.title = element_text(size = 6), # Adjust legend title size
  plot.title = element_text(size = 12, hjust = 0.5), # Center title
  panel.background = element_blank(), # Remove background color
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Add border
  panel.grid.major = element_blank(), # Subtle grid
  panel.grid.minor = element_blank(), # No minor grid
  plot.margin = margin(5, 5, 5, 5) # Adjust margins
)

#### Data import and cleaning ####
data_list = file.path("input_folder", "Data_AB012I_transwell_drugs.csv")

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
  colnms[[i]] <- as.character(c("Sample_name",
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
merged_df$group = sapply(strsplit(as.character(merged_df$Sample_name), "_"), '[', 2)
merged_df$TP = sapply(strsplit(as.character(merged_df$Sample_name), "_"), '[', 3)
merged_df$rep = sapply(strsplit(as.character(merged_df$Sample_name), "_"), '[', 6)
merged_df$side = sapply(strsplit(as.character(merged_df$Sample_name), "_"), '[', 7)
merged_df$time = sapply(strsplit(as.character(merged_df$Sample_name), "_"), '[', 4)
merged_df$condition = sapply(strsplit(as.character(merged_df$Sample_name), "_"), '[', 5)


# Change area variables to numeric
cols = grep("_Area", colnames(merged_df))
merged_df[, cols] = apply(merged_df[,cols], 2, function(x) as.numeric(as.character(x)))

# Change TPs to numeric 
merged_df$time = as.numeric(gsub("T", "", merged_df$time))

cols = grep("marker.*_Area", colnames(merged_df))

#### Plot markers ####
# plot raw data 
pdf("./Plots/AB012_transwell2_raw_markers.pdf")
for (i in cols){
  
  y_var = sym(colnames(merged_df)[i])
  
  # Default line plot
  p <-  ggplot(merged_df, aes(x=time, y=!!y_var, color = side, shape = rep, 
                              group = interaction(side, rep, TP))) +
    geom_point() +
    geom_line(aes(group=interaction(rep, side, TP))) +
    labs(title= colnames(merged_df)[i], 
         subtitle="Timecourse", 
         x = "Time (minutes)", 
         y = "Peak area Intensity (raw)") + 
    theme_minimal() + 
    facet_wrap(~group)
  
  print(p)
  
}
dev.off()


#### Plot IS ####
IS_cols = grep("IS.*_Area", colnames(merged_df))

# plot raw data 
pdf("./Plots/AB012_transwell2_raw_ISsignal.pdf")
for (i in IS_cols){
  
  y_var = sym(colnames(merged_df)[i])
  
  # Default line plot
  p <-  ggplot(merged_df, aes(x=time, y=!!y_var, color = side, shape = rep, 
                              group = interaction(side, rep, TP, group))) +
    geom_point() +
    geom_line(aes(group=interaction(rep, side, TP, group))) +
    labs(title= colnames(merged_df)[i], 
         subtitle="Timecourse", 
         x = "Time (minutes)", 
         y = "Peak area Intensity (raw)") + 
    theme_minimal() 
  
  print(p)
  
}
dev.off()


#### Plot drugs & metabolites ####
##### Restructure df ####
# pivot merged df 
area_cols = grep("_Area$", colnames(merged_df), value = T)
my_prep_df = merged_df[, c("Sample_name", "condition", "TP", "rep", "time", "group", "side", area_cols)]
my_prep_df = pivot_longer(my_prep_df, 
                          cols = grep("_Area$", colnames(my_prep_df), value = T), 
                          values_to = "Area", 
                          names_to = "drug")

my_prep_df = my_prep_df %>% 
  filter(!grepl("water|IS", Sample_name))



##### Roxatidine acetate ####
colnames(merged_df)
d = "drug_roxatidine.Results_Area"
m = "met_roxatidine.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m) 


drug_labs = c("drug_roxatidine.Results_Area" = "Roxatidine acetate", 
              "met_roxatidine.Results_Area" = "Roxatidine")
strain_labs = c("GMM" = "Sterile control", 
                "MB002" = "donor 1", 
                "MB003" = "donor 2", 
                "MB005" = "donor 3", 
                "mouse" = "mouse")


temp_df %>% 
  #filter(TP == "TP16") %>% 
  #filter(condition == "drugs") %>% 
  ggplot(., aes(x=time, y=Area, color = side, shape = rep, 
                            group = interaction(side, rep, condition))) +
  geom_point() +
  geom_line(aes(group=interaction(rep, side, condition))) +
  labs(
       x = "Time (minutes)", 
       y = "Peak area Intensity (raw)") + 
  theme_minimal() + 
  facet_grid(drug+condition~group+TP, 
             labeller = labeller(drug = drug_labs, group = strain_labs),
             scales = "free_y") +
  scale_color_manual(
    values = c("api" = "lightblue", "baso" = "lightpink"),
    labels = c("api" = "Apical", "baso" = "Basolateral")
  ) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific) +
  cell_plot_theme

summ_df <- temp_df %>%
  filter(TP == "TP16",
         condition == "drugs",
         group != "mouse") %>%
  group_by(side, drug, group, time) %>%
  summarise(
    mean_area = mean(Area, na.rm = TRUE),
    se_area = sd(Area, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# --- Plot for API ---
p_api <- summ_df %>%
  filter(side == "api") %>%
  ggplot(aes(x = time, y = mean_area, 
             group = group, color = side)) +
  
  # SE ribbon (not smoothed)
  geom_ribbon(
    aes(ymin = mean_area - se_area,
        ymax = mean_area + se_area,
        fill = side),
    alpha = 0.2,
    color = NA
  ) +
  
  # Smoothed mean line
  geom_smooth(
    method = "loess",
    se = FALSE,
    aes(linetype = group),
    linewidth = 1.2
  ) +
  
  facet_wrap(~ drug,
             labeller = labeller(drug = drug_labs),
             scales = "free_y") +
  scale_color_manual(values = c("api" = "lightblue")) +
  scale_fill_manual(values = c("api" = "lightblue")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific, n.breaks = 4) +
  labs(x = "", y = "Peak area (API)") +
  theme_minimal() +
  cell_plot_theme +
  theme(legend.position = "none")

# --- Plot for BASO ---
p_baso <- summ_df %>%
  filter(side == "baso") %>%
  ggplot(aes(x = time, y = mean_area, 
             group = group, color = side)) +
  
  geom_ribbon(
    aes(ymin = mean_area - se_area,
        ymax = mean_area + se_area,
        fill = side),
    alpha = 0.2,
    color = NA
  ) +
  
  geom_smooth(
    method = "loess",
    se = FALSE,
    aes(linetype = group),
    linewidth = 1.2
  ) +
  
  facet_wrap(~ drug,
             labeller = labeller(drug = drug_labs),
             scales = "fixed") +
  scale_color_manual(values = c("baso" = "lightpink")) +
  scale_fill_manual(values = c("baso" = "lightpink")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific) +
  labs(x = "Time (minutes)", y = "Peak area (BASO)") +
  theme_minimal() +
  cell_plot_theme +
  theme(legend.position = "none")


# --- Combine with patchwork ---
p_api / p_baso

ggsave("./Plots/roxace_transwell.tiff", dpi = 600, width = 3.5, height = 3, units = "in")


##### Nicergoline ####
d = "drug_nicergoline.Results_Area"
m = "met_nicergoline.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m)


drug_labs = c("drug_nicergoline.Results_Area" = "Nicergoline", 
              "met_nicergoline.Results_Area" = "Nicergoline metabolite")

temp_df %>% 
  #filter(TP == "TP16") %>% 
  filter(condition == "drugs") %>% 
  ggplot(., aes(x=time, y=Area, color = side, shape = rep, 
                group = interaction(side, rep, condition))) +
  geom_point() +
  geom_line(aes(group=interaction(rep, side, condition))) +
  labs(
    x = "Time (minutes)", 
    y = "Peak area Intensity (raw)") + 
  theme_minimal() + 
  facet_grid(drug+condition~group+TP, 
             labeller = labeller(drug = drug_labs, group = strain_labs),
             scales = "free_y") +
  scale_color_manual(
    values = c("api" = "lightblue", "baso" = "lightpink"),
    labels = c("api" = "Apical", "baso" = "Basolateral")
  ) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific) +
  cell_plot_theme


##### Bisacodyl ####
d = "drug_bisacodyl.Results_Area"
m = "met_bisacodyl278.Results_Area"
m1 = "met_bisacodyl320.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m | drug == m1)


drug_labs = c("drug_bisacodyl.Results_Area" = "Bisacodyl", 
              "met_bisacodyl278.Results_Area" = "Metabolite (278 m/z)", 
              "met_bisacodyl320.Results_Area" = "Metabolite (320 m/z)")


temp_df %>% 
  #filter(TP == "TP16") %>% 
  #filter(condition == "drugs") %>% 
  ggplot(., aes(x=time, y=Area, color = side, shape = rep, 
                group = interaction(side, rep, condition))) +
  geom_point() +
  geom_line(aes(group=interaction(rep, side, condition))) +
  labs(
    x = "Time (minutes)", 
    y = "Peak area Intensity (raw)") + 
  theme_minimal() + 
  facet_grid(drug+condition~group+TP, 
             labeller = labeller(drug = drug_labs, group = strain_labs),
             scales = "free_y") +
  scale_color_manual(
    values = c("api" = "lightblue", "baso" = "lightpink"),
    labels = c("api" = "Apical", "baso" = "Basolateral")
  ) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific) +
  cell_plot_theme



##### MMF ####
d = "drug_mycophenolatemofetil.Results_Area"
m = "met_mmf.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m) 

drug_labs = c("drug_mycophenolatemofetil.Results_Area" = "Mycophenolate Mofetil", 
              "met_mmf.Results_Area" = "Mycophenolic acid")

temp_df %>% 
  filter(TP == "TP16") %>% 
  filter(condition == "drugs") %>% 
  filter(!(drug == "drug_mycophenolatemofetil.Results_Area" & Area < 6000)) %>% 
  ggplot(., aes(x=time, y=Area, color = side, 
                group = interaction(side, rep, condition))) +
  #geom_line(aes(group=interaction(rep, side, condition, group), linetype = group)) +
  geom_smooth(aes(linetype = group, group = interaction(side, condition, group)), method = "loess", linewidth = 1.5, se = F) +
  labs(
    x = "Time (minutes)", 
    y = "Peak area Intensity (raw)") + 
  theme_minimal() + 
  facet_grid(side~condition+drug, 
              labeller = labeller(drug = drug_labs, group = strain_labs),
              scales = "free") +
  scale_color_manual(
  values = c("api" = "lightblue", "baso" = "lightpink"),
  labels = c("api" = "Apical", "baso" = "Basolateral")
  ) +
  #scale_color_manual(values = c("#DDA0DD", "#5F9EA0", "#E75480","#F08080", "#FFD1DC", "#FFA07A")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific) +
  cell_plot_theme

ggsave("./Plots/MMF_transwell.tiff", dpi = 600, width = 3.5, height = 3, units = "in")

library(patchwork)

# --- Plot for API ---
p_api <- temp_df %>% 
  filter(TP == "TP16",
         condition == "drugs",
         side == "api",
         !(drug == "drug_mycophenolatemofetil.Results_Area" & Area < 6000)) %>% 
  filter(group != "mouse") %>% 
  ggplot(aes(x = time, y = Area, color = side,
             group = interaction(rep, condition))) +
  geom_smooth(aes(linetype = group, 
                  group = interaction(side, condition, group)), method = "loess", linewidth = 1.5, se = F) +
  facet_wrap(~ drug,
             labeller = labeller(drug = drug_labs),
             scales = "free_y") +
  labs(x = "", y = "Peak area (API)") +
  scale_color_manual(values = c("api" = "lightblue")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific) +
  theme_minimal() +
  cell_plot_theme +
  theme(legend.position = "none")

# --- Plot for BASO ---
p_baso <- temp_df %>% 
  filter(TP == "TP16",
         condition == "drugs",
         side == "baso",
         !(drug == "drug_mycophenolatemofetil.Results_Area" & Area < 6000)) %>% 
  filter(group != "mouse") %>% 
  ggplot(aes(x = time, y = Area, color = side,
             group = interaction(rep, condition))) +
  geom_smooth(aes(linetype = group, 
                  group = interaction(side, condition, group)), method = "loess", linewidth = 1.5, se = F) +
  facet_wrap(~ drug,
             labeller = labeller(drug = drug_labs),
             scales = "free_y") +
  labs(x = "Time (minutes)", y = "Peak area (BASO)") +
  scale_color_manual(values = c("baso" = "lightpink")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific) +
  theme_minimal() +
  cell_plot_theme + 
  theme(legend.position = "none")

# --- Combine with patchwork ---
p_api / p_baso

ggsave("./Plots/MMF_transwell2.tiff", dpi = 600, width = 3.5, height = 3, units = "in")


##### Carfilzomib ####
d = "drug_carfilzomib.Results_Area"
m = "met_carfilzomib420.Results_Area"
m1 = "met_carfilzomib442.Results_Area"
m2 = "met_carfilzomib737.Results_Area" 
m3 = "met_carfilzomib319.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug %in% c(m, m1, m2, m3, d)) 


drug_labs = c("drug_carfilzomib.Results_Area" = "Carfilzomib", 
              "met_carfilzomib420.Results_Area" = "Metabolite (420 m/z)", 
              "met_carfilzomib442.Results_Area" = "Metabolite (442 m/z)", 
              "met_carfilzomib737.Results_Area" = "Metabolite (737 m/z)", 
              "met_carfilzomib319.Results_Area" = "Metabolite (319 m/z)")

temp_df %>% 
  filter(TP == "TP16") %>% 
  filter(condition == "drugs") %>% 
  ggplot(., aes(x=time, y=Area, color = side, shape = rep, 
                group = interaction(side, rep, condition))) +
  geom_point() +
  geom_line(aes(group=interaction(rep, side, condition))) +
  labs(
    x = "Time (minutes)", 
    y = "Peak area Intensity (raw)") + 
  theme_minimal() + 
  facet_wrap(drug+condition~group+TP, 
             labeller = labeller(drug = drug_labs, group = strain_labs), 
             nrow = 5, 
             scales = "free_y") +
  scale_color_manual(
    values = c("api" = "lightblue", "baso" = "lightpink"),
    labels = c("api" = "Apical", "baso" = "Basolateral")
  ) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific) +
  cell_plot_theme

# --- Plot for API ---
p_api <- temp_df %>% 
  filter(TP == "TP16",
         condition == "drugs",
         side == "api",
         drug %in% c(d, m3)) %>% 
  filter(group != "mouse") %>% 
  ggplot(aes(x = time, y = Area, color = side,
             group = interaction(rep, condition))) +
  geom_smooth(aes(linetype = group, 
                  group = interaction(side, condition, group)), method = "loess", linewidth = 1.5, se = F) +
  facet_wrap(~ drug,
             labeller = labeller(drug = drug_labs),
             scales = "free_y") +
  labs(x = "", y = "Peak area (API)") +
  scale_color_manual(values = c("api" = "lightblue")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific) +
  theme_minimal() +
  cell_plot_theme +
  theme(legend.position = "none")

# --- Plot for BASO ---
p_baso <- temp_df %>% 
  filter(TP == "TP16",
         condition == "drugs",
         side == "baso",
         drug %in% c(d, m3)) %>% 
  filter(group != "mouse") %>% 
  ggplot(aes(x = time, y = Area, color = side,
             group = interaction(rep, condition))) +
  geom_smooth(aes(linetype = group, 
                  group = interaction(side, condition, group)), method = "loess", linewidth = 1.5, se = F) +
  facet_wrap(~ drug,
             labeller = labeller(drug = drug_labs),
             scales = "free_y") +
  labs(x = "Time (minutes)", y = "Peak area (BASO)") +
  scale_color_manual(values = c("baso" = "lightpink")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific) +
  theme_minimal() +
  cell_plot_theme + 
  theme(legend.position = "none")

# --- Combine with patchwork ---
p_api / p_baso

##### Deflazacort ####
d = "drug_deflazacort.Results_Area"
m = "met_deflazacort.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m) 

drug_labs = c("drug_deflazacort.Results_Area" = "Deflazacort", 
              "met_deflazacort.Results_Area" = "Deflazacort metabolite")


temp_df %>% 
  #filter(TP == "TP16") %>% 
  #filter(condition == "drugs") %>% 
  ggplot(., aes(x=time, y=Area, color = side, shape = rep, 
                group = interaction(side, rep, condition))) +
  geom_point() +
  geom_line(aes(group=interaction(rep, side, condition))) +
  labs(
    x = "Time (minutes)", 
    y = "Peak area Intensity (raw)") + 
  theme_minimal() + 
  facet_grid(drug+condition~group+TP, 
             labeller = labeller(drug = drug_labs, group = strain_labs),
             scales = "free_y") +
  scale_color_manual(
    values = c("api" = "lightblue", "baso" = "lightpink"),
    labels = c("api" = "Apical", "baso" = "Basolateral")
  ) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific) +
  cell_plot_theme



##### Racecadotril ####
d = "drug_racecadotril.Results_Area"
m = "met_racecadotril.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m)


drug_labs = c("drug_racecadotril.Results_Area" = "Racecadotril", 
              "met_racecadotril.Results_Area" = "Racecadotril metabolite")

temp_df %>% 
  #filter(TP == "TP16") %>% 
  filter(condition == "drugs") %>% 
  ggplot(., aes(x=time, y=Area, color = side, shape = rep, 
                group = interaction(side, rep, condition))) +
  geom_point() +
  geom_line(aes(group=interaction(rep, side, condition))) +
  labs(
    x = "Time (minutes)", 
    y = "Peak area Intensity (raw)") + 
  theme_minimal() + 
  facet_grid(drug+condition~group+TP, 
             labeller = labeller(drug = drug_labs, group = strain_labs),
             scales = "free_y") +
  scale_color_manual(
    values = c("api" = "lightblue", "baso" = "lightpink"),
    labels = c("api" = "Apical", "baso" = "Basolateral")
  ) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific) +
  cell_plot_theme


##### Azilsartan ####
d = "drug_azilsartan.Results_Area"
m = "met_azilsartan.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m)


drug_labs = c("drug_azilsartan.Results_Area" = "Azilsartan", 
              "met_azilsartan.Results_Area" = "Azilsartan metabolite")

temp_df %>% 
  #filter(TP == "TP16") %>% 
  #filter(condition == "drugs") %>% 
  ggplot(., aes(x=time, y=Area, color = side, shape = rep, 
                group = interaction(side, rep, condition))) +
  geom_point() +
  geom_line(aes(group=interaction(rep, side, condition))) +
  labs(
    x = "Time (minutes)", 
    y = "Peak area Intensity (raw)") + 
  theme_minimal() + 
  facet_grid(drug+condition~group+TP, 
             labeller = labeller(drug = drug_labs, group = strain_labs),
             scales = "free_y") +
  scale_color_manual(
    values = c("api" = "lightblue", "baso" = "lightpink"),
    labels = c("api" = "Apical", "baso" = "Basolateral")
  ) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific) +
  cell_plot_theme


##### Dexamethasone ####
data_list = file.path("input_folder", "Data_AB012I_transwell_drugs.csv")

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
  colnms[[i]] <- as.character(c("Sample_name", "data_file", "type", "level", "timestamp",
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
merged_df$group = sapply(strsplit(as.character(merged_df$Sample_name), "_"), '[', 2)
merged_df$TP = sapply(strsplit(as.character(merged_df$Sample_name), "_"), '[', 3)
merged_df$rep = sapply(strsplit(as.character(merged_df$Sample_name), "_"), '[', 6)
merged_df$side = sapply(strsplit(as.character(merged_df$Sample_name), "_"), '[', 7)
merged_df$time = sapply(strsplit(as.character(merged_df$Sample_name), "_"), '[', 4)
merged_df$condition = sapply(strsplit(as.character(merged_df$Sample_name), "_"), '[', 5)


# Change area variables to numeric
cols = grep("_Area", colnames(merged_df))
merged_df[, cols] = apply(merged_df[,cols], 2, function(x) as.numeric(as.character(x)))

# Change TPs to numeric 
merged_df$time = as.numeric(gsub("T", "", merged_df$time))

# pivot merged df 
area_cols = grep("_Area$", colnames(merged_df), value = T)
my_prep_df = merged_df[, c("Sample_name", "condition", "TP", "rep", "time", "group", "side", area_cols)]
my_prep_df = pivot_longer(my_prep_df, 
                          cols = grep("_Area$", colnames(my_prep_df), value = T), 
                          values_to = "Area", 
                          names_to = "drug")

my_prep_df = my_prep_df %>% 
  filter(!grepl("water|IS", Sample_name))


unique(my_prep_df$drug)
d = "drug_dexamethasone.Results_Area"
m = "met_dexamethasone.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m)


drug_labs = c("drug_dexamethasone.Results_Area" = "Dexamethasone", 
              "met_dexamethasone.Results_Area" = "Dexamethasone metabolite")

temp_df %>% 
  #filter(TP == "TP16") %>% 
  #filter(condition == "drugs") %>% 
  ggplot(., aes(x=time, y=Area, color = side, shape = rep, 
                group = interaction(side, rep, condition))) +
  geom_point() +
  geom_line(aes(group=interaction(rep, side, condition))) +
  labs(
    x = "Time (minutes)", 
    y = "Peak area Intensity (raw)") + 
  theme_minimal() + 
  facet_grid(drug+condition~group+TP, 
             labeller = labeller(drug = drug_labs, group = strain_labs),
             scales = "free_y") +
  scale_color_manual(
    values = c("api" = "lightblue", "baso" = "lightpink"),
    labels = c("api" = "Apical", "baso" = "Basolateral")
  ) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific) +
  cell_plot_theme


temp_df %>% 
  filter(TP == "TP16") %>% 
  filter(condition == "drugs") %>% 
  ggplot(., aes(x=time, y=Area, color = side, 
                group = interaction(side, rep, condition))) +
  #geom_line(aes(group=interaction(rep, side, condition, group), linetype = group)) +
  geom_smooth(aes(linetype = group, group = interaction(side, condition, group)), method = "loess", linewidth = 1.5, se = F) +
  labs(
    x = "Time (minutes)", 
    y = "Peak area Intensity (raw)") + 
  theme_minimal() + 
  facet_grid(side~condition+drug, 
             labeller = labeller(drug = drug_labs, group = strain_labs),
             scales = "free") +
  scale_color_manual(
    values = c("api" = "lightblue", "baso" = "lightpink"),
    labels = c("api" = "Apical", "baso" = "Basolateral")
  ) +
  #scale_color_manual(values = c("#DDA0DD", "#5F9EA0", "#E75480","#F08080", "#FFD1DC", "#FFA07A")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific) +
  cell_plot_theme

ggsave("./Plots/dexamethasone_transwell.tiff", dpi = 600, width = 3.5, height = 3, units = "in")
ggsave("./Plots/dexamethasone_transwell.png", dpi = 600, width = 3.5, height = 3, units = "in")





