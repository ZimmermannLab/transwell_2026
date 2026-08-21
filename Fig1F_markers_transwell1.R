#### Script to analyse marker compouns in transwell experiment 1 ####

## Set up 
####Load Libraries#####
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(ggforce)

#### Data import and cleaning ####
data_list = file.path("input_folder", "Data_AB012I_markers_transwell1.csv")

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
merged_df$rep = sapply(strsplit(as.character(merged_df$Sample_name), "_"), '[', 4)
merged_df$side = sapply(strsplit(as.character(merged_df$Sample_name), "_"), '[', 5)
merged_df$exp = sapply(strsplit(as.character(merged_df$Sample_name), "_"), '[', 1)

# Change area variables to numeric
cols = grep("_Area", colnames(merged_df))
merged_df[, cols] = apply(merged_df[,cols], 2, function(x) as.numeric(as.character(x)))

# Change TPs to numeric 
merged_df$time = as.numeric(gsub("T", "", merged_df$TP))

cols = grep("marker.*_Area", colnames(merged_df))

# plot raw data 
pdf("./Plots/AB012I_transwell1_raw_markers.pdf")
for (i in cols){
  
  y_var = sym(colnames(merged_df)[i])
  
  # Default line plot
  p <-  ggplot(merged_df, aes(x=time, y=!!y_var, color = side, shape = rep, 
                              group = interaction(side, rep))) +
    geom_point() +
    geom_line(aes(group=interaction(rep, side))) +
    labs(title= colnames(merged_df)[i], 
         subtitle="Timecourse", 
         x = "Time (minutes)", 
         y = "Peak area Intensity (raw)") + 
    theme_minimal() + 
    facet_wrap(~group)
  
  print(p)
  
}
dev.off()



# delete waterwashes and IS samples 
merged_df = merged_df[!grepl("waterwash|IS|Water", merged_df$Sample_name),]


### plotting for manuscript ###

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
  dplyr::select(Sample_name, group, exp, time, rep, side, all_of(grep("marker.*_Area", colnames(merged_df))))

tmp = pivot_longer(tmp, cols = grep("marker.*_Area", colnames(tmp)))

non_per = c("marker_nadolol.Results_Area", "marker_terbutaline.Results_Area", "marker_etoposide.Results_Area")
tmp$permeable = ifelse(tmp$name %in% non_per, "non_permeable", "permeable")

fig1 = tmp %>% 
  filter(name %in% c("marker_nadolol.Results_Area"))

drug_labs = c("marker_antipyrine.Results_Area" = "Antipyrine", 
              "marker_nadolol.Results_Area" = "Nadolol")

p <-  fig1 %>% 
  #filter(group %in% c("GMM", "MB002")) %>% 
  ggplot(., aes(x=time, y=value, color = side, shape = group)) +
  geom_point(aes(shape = group), size=0.7) +
  #geom_line(aes(group=interaction(rep, side, exp))) +
  geom_smooth(se = F, method = "auto", linewidth = 0.8, 
              aes(group = interaction(side, name, group))) +
  labs(x = "Time (minutes)", 
       y = "Peak Area Intensity") + 
  theme_minimal() + 
  theme(legend.position = "bottom") +
  scale_color_manual(
    values = c("api" = "lightblue", "baso" = "lightpink")) +
  #scale_color_manual(values = tar_colors) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  facet_grid(~name, scales = "free_y", labeller = labeller(name = drug_labs)) +
  scale_y_continuous(labels = scales::scientific, breaks = c(0e0, 1e6, 2e6)) +
  cell_plot_theme +
  guides(color = "none")

print(p)
ggsave("./Plots/AB012I_markers_transwell1_nadolol_allgroups.tiff", plot = p, dpi = 600, width = 3, height = 2.5, units = "in")

### all for suppl. figure ###
drug_labs = c("marker_antipyrine.Results_Area" = "Antipyrine", 
              "marker_nadolol.Results_Area" = "Nadolol", 
              "marker_terbutaline.Results_Area" = "Terbutaline", 
              "marker_propanolol.Results_Area" = "Propranolol", 
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
  facet_grid(name~group, scales = "free_y", labeller = labeller(name = drug_labs)) +
  scale_y_continuous(labels = scales::scientific) +
  theme(legend.position = "none") +
  cell_plot_theme
print(p)
ggsave("AB012I_allmarkers_transwell1.tiff", plot = p, dpi = 600, width = 6, height = 6.5, units = "in")


tmp %>% 
  filter(group == "GMM") %>% 
ggplot(., aes(x=time, y=value, color = side)) +
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
  facet_grid(name~group, scales = "free_y", labeller = labeller(name = drug_labs)) +
  scale_y_continuous(labels = scales::scientific, n.breaks = 3) +
  cell_plot_theme


#### Calculate Peff for markers ####
calculate_peff <- function(tmp,
                           marker_name,
                           group = "GMM",
                           side_baso = "baso",
                           side_donor = "api",
                           time_start = 30,
                           time_end = 150,
                           volume = 100,
                           surface = 0.11) {
  
  # basolateral data for the chosen marker
  baso_data <- tmp %>%
    dplyr::filter(
      group == group,
      side == side_baso,
      name == marker_name
    )
  
  # donor concentration at time 0
  donor_conc <- mean(
    tmp$value[
      tmp$side == side_donor &
        tmp$name == marker_name &
        tmp$time == 0
    ],
    na.rm = TRUE
  )
  
  # dC/dt
  dC_dt <- mean(baso_data$value[baso_data$time == time_end], na.rm = TRUE) -
    mean(baso_data$value[baso_data$time == time_start], na.rm = TRUE)
  
  # Peff
  Peff <- (volume / (surface * donor_conc)) * dC_dt
  
  return(Peff)
}

peff_table <- tmp %>%
  distinct(name) %>%
  mutate(
    Peff = map_dbl(
      name,
      ~ calculate_peff(
        tmp = tmp,
        marker_name = .x
      )
    )
  )

peff_table$FA = c(73, 35, NA, 90, NA, 50, 98)

ggplot(peff_table, aes(x=log(Peff), y=FA)) +
  geom_point()

