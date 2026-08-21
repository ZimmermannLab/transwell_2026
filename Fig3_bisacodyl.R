# Set up 
######Load Libraries########
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(ggforce)

#### Data import and cleaning ####

#read them
screen <- list(read.csv(file.path("input_folder", "AB012_drugs_mets_table.csv")))


# run the following data cleaning functions for the screen list 
# extract the targeted drugs names (+ remove spaces)
targeted_drugs = function(df){
  gsub(" ", "_", gsub(" Results","",grep("Results",colnames(df),fixed=T,value=T),fixed = T), fixed = T)
}


target_drugs = lapply(screen, targeted_drugs)

# update column names
colnms <- list()
for (i in seq_along(target_drugs)) { 
  colnms[[i]] <- as.character(c("sample_name", "data_file", "type", "level", "time_stamp", 
                                as.vector(sapply(target_drugs[[i]], function(x) paste(x,c("_RT","_Area"),sep="")))))
}
# for table3
for (i in seq_along(target_drugs)) { 
  colnms[[i]] <- as.character(c("sample_name",
                                as.vector(sapply(target_drugs[[i]], function(x) paste(x,c("_RT","_Area"),sep="")))))
}
#
# for table carfil
colnms <- list()
for (i in seq_along(target_drugs)) { 
  colnms[[i]] <- as.character(c("nothing1", "nothing2", "sample_name", "data_file", "type", "level", "time_stamp", 
                                as.vector(sapply(target_drugs[[i]], function(x) paste(x,c("_RT","_Area"),sep="")))))
}
#

for (i in seq_along(screen)){
  colnames(screen[[i]]) <- colnms[[i]]
}

# remove first row 
# set empty to NA
cleaning_funct = function(df){
  df = df %>% 
    dplyr::slice(-1) %>% 
    mutate_all(list(~na_if(.,"")))
}

screen = lapply(screen, cleaning_funct)

lapply(screen, colnames)

# get cleaned df 
merged_df <- bind_rows(screen, .id = "column_label")
merged_df$sample_name[6]

# Add columns with timepoints & conditions 
merged_df$strain = sapply(strsplit(as.character(merged_df$sample_name), "_"), '[', 2)
merged_df$condition = sapply(strsplit(as.character(merged_df$sample_name), "_"), '[', 3)
merged_df$TPs = sapply(strsplit(as.character(merged_df$sample_name), "_"), '[', 4)
merged_df$Rep = sapply(strsplit(as.character(merged_df$sample_name), "_"), '[', 6)
merged_df$side = sapply(strsplit(as.character(merged_df$sample_name), "_"), '[', 7)
merged_df$pool = sapply(strsplit(as.character(merged_df$sample_name), "_"), '[', 5)

# Change area variables to numeric
cols = grep("_Area", colnames(merged_df))
merged_df[, cols] = apply(merged_df[,cols], 2, function(x) as.numeric(as.character(x)))

# Change the timepoints to numeric 
merged_df$Time = as.numeric(gsub( "T", "", as.character(merged_df$TPs)))


#### Define plotting specs ####
library(showtext)

# Enable Arial font
font_add("Arial", regular = "arial.ttf", bold = "arialbd.ttf") # Ensure Arial font files are accessible
showtext_auto()
# Function to standardize ggplot theme
cell_plot_theme <- theme(
  text = element_text(family = "Arial", size = 6), # Set Arial font and text size
  axis.text = element_text(size = 8), # Adjust axis text size
  axis.title = element_text(size = 8), # Adjust axis title size
  legend.text = element_text(size = 8), # Adjust legend text size
  legend.title = element_text(size = 8), # Adjust legend title size
  plot.title = element_text(size = 12, hjust = 0.5), # Center title
  panel.background = element_blank(), # Remove background color
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Add border
  panel.grid.major = element_blank(), # Subtle grid
  panel.grid.minor = element_blank(), # No minor grid
  plot.margin = margin(5, 5, 5, 5) # Adjust margins
)

#### Bisacodyl ####
unique(my_prep_df$drug)

d = "drug_bisacodyl.Results_Area"
m = "met_bisacodyl_1.Results_Area"


temp_df = my_prep_df %>% 
  filter(drug == d | drug == m) %>%
  filter(drug == d & strain == "nobug" | drug == m & strain == "Ef") %>% 
  #filter(Area > 1000) %>% 
  na.omit()

t_df = temp_df %>% 
  filter(condition == "TP7")

# Smooth the data by averaging the value at each time point with its two flanking points
t_df <- t_df %>%
  group_by(side, drug, strain) %>%
  arrange(Time) %>%
  mutate(smooth_value = (Area + lag(Area, 1, order_by = Time) + lead(Area, 1, order_by = Time)) / 3) %>%
  ungroup()


drug_labs = c("drug_bisacodyl.Results_Area" = "Bisacodyl", 
              "met_bisacodyl_1.Results_Area" = "Active drug metabolite")
strain_labs = c("nobug" = "Sterile control", 
                "Ef" = "E. faecalis")

# Define empty panels with proper data types
empty_panels <- data.frame(
  Time = c(0, 0), # Use valid numeric values
  Area = c(NA_real_, NA_real_), # Ensure Area is numeric
  side = c("aa", "aa")) # Character column


# Default line plot
p <- temp_df %>% 
  #filter(side == "baso") %>% 
  filter(strain %in% c("Ef", "nobug") & condition == "TP7") %>% 
  ggplot(., aes(x=Time, y=Area, color=side, group=side))+
  #stat_summary(fun.data ='mean_se', geom = "smooth", se = TRUE, alpha=0.1) +
  geom_point() +
  geom_smooth(method = "loess", linewidth = 1.5, se = F) +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  labs(#title= "Bisacodyl and metabolite", 
    #subtitle="Timecourse", 
    x = "Time (minutes)", 
    y = "Peak area Intensity") + 
  theme_minimal() + 
  #facet_wrap(~condition+drug+strain) + 
  facet_grid(side~drug+strain, scales = "free_y", 
             labeller = labeller(drug = drug_labs, strain = strain_labs, 
                                 side = function(labels) rep("", length(labels))), 
  ) + 
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific, n.breaks = 3) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("api" = "lightblue", "baso" = "lightpink")) +
  theme(
    strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1), lineheight = 0.8), # Adjust facet label size and margin
    panel.spacing = unit(0.5, "lines")) +
  cell_plot_theme +
  geom_blank(data = empty_panels)
print(p)
ggsave("bisacodyl.tiff", plot = p, dpi = 300, width = 3, height = 3.5, units = "in") # Cell guidelines recommend these dimensions
