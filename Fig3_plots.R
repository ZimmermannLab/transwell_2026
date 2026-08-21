# Set up 
######Load Libraries########
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(ggforce)

#### Data import and cleaning ####

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


#### Plot raw values ####

# pivot merged df 
area_cols = grep("_Area$", colnames(merged_df), value = T)
my_prep_df = merged_df[, c("sample_name", "condition", "side", "Time", "Rep", "strain", "pool", area_cols)]
my_prep_df = pivot_longer(my_prep_df, 
                          cols = grep("_Area$", colnames(my_prep_df), value = T), 
                          values_to = "Area", 
                          names_to = "drug")

# check drugs only 
drugs = my_prep_df$drug[grep("drug_|met_", unique(my_prep_df$drug))]

# Save plots in pdf
pdf("./plots/AB012E3_raw_values_api&baso.pdf")
for (i in drugs){
  
  temp_df = my_prep_df %>% 
    filter(drug == i) 
  temp_df = temp_df[!is.na(temp_df$Rep), ]
  
  # Default line plot
  
  p <-  ggplot(temp_df, aes(x=Time, y=Area, 
                            group=interaction(condition, side, strain), 
                            color=side)) +
    geom_point() +
    geom_line(aes(group=interaction(Rep, side, strain))) +
    labs(title= i, 
         subtitle="Timecourse", 
         x = "Time (minutes)", 
         y = "Peak area Intensity") + 
    theme_minimal() + 
    facet_wrap(~strain+condition) + 
    geom_blank()
  
  print(p)
  
}

dev.off()


# plot single combination of drug-bug 


#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

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


#### Roxatine ace #####
unique(my_prep_df$drug)
d = "drug_roxatidine.Results_Area"
m = "met_roxatidine.Results_Area"
s = "Bu"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m) %>% 
  filter(strain == s | strain == "nobug" & pool == "pool2")


t_df = temp_df %>% 
  filter(drug == d & strain == "nobug" | drug == m & condition == "TP8" & strain == "Bu")

drug_labs = c("drug_roxatidine.Results_Area" = "Roxatine Acetate", 
              "met_roxatidine.Results_Area" = "Bacterial metabolite:\n deacetyl-Roxatine \nacetate")
bug_labs = c("nobug" = "Sterile control", 
             "Bu"  = "B. uniformis")

# Smooth the data by averaging the value at each time point with its two flanking points
t_df <- t_df %>%
  group_by(side, drug, strain) %>%
  arrange(Time) %>%
  mutate(smooth_value = (Area + lag(Area, 1, order_by = Time) + lead(Area, 1, order_by = Time)) / 3) %>%
  ungroup()

# Default line plot
p1 <- ggplot(t_df, aes(x=Time, y=Area, color=side, group=side))+
  geom_point() +
  #stat_summary(aes(y = smooth_value, group = side), fun=mean, geom="line", linewidth = 1.5) +
  #stat_summary(fun.data ='mean_se', geom = "smooth", se = TRUE, alpha=0.1) +
  geom_smooth(se = F, method = "loess", linewidth = 1.5) +
  labs(#title= "Roxatidine Acetate and metabolite", 
       #subtitle="Timecourse", 
       x = "Time (minutes)", 
       y = "Peak area Intensity") + 
  theme_minimal() + 
  theme(legend.position = "none") +
  facet_grid(side~drug+strain, labeller = labeller(drug = drug_labs, 
                                                   strain = bug_labs, 
                                                   side = function(labels) rep("", length(labels)))) + 
  scale_color_manual(values = c("api" = "lightblue", "baso" = "lightpink")) + 
  theme(
    strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1), lineheight = 0.8), # Adjust facet label size and margin
    panel.spacing = unit(0.5, "lines")) + # Reduce space between panels+
  cell_plot_theme
print(p1)
# Save the plot
ggsave("roxatidine.tiff", plot = p, dpi = 300, width = 4.5, height = 3.5, units = "in") # Cell guidelines recommend these dimensions


#### Azilsartan ####
unique(my_prep_df$drug)
d = "drug_azilsartan.Results_Area"
m = "met_azilsartan.Results_Area"
#library(ggbreak)

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m) %>% 
  na.omit()

t_df = temp_df %>% 
  filter(drug == d & condition == "TP0" | drug == m & condition == "TP7") %>% 
  #filter(Area > 5000) %>% 
  filter(strain == "Ef")

drug_labs = c("drug_azilsartan.Results_Area" = "Azilsartan", 
              "met_azilsartan.Results_Area" = "Bacterial metabolite")
condition_labs = c("TP0" = "0 hour bacterial incubation", 
             "TP7"  = "7 hour bacterial incubation")

# Define empty panels with proper data types
empty_panels <- data.frame(
  Time = c(0, 0), # Use valid numeric values
  Area = c(NA_real_, NA_real_), # Ensure Area is numeric
  side = c("aa", "aa")) # Character column

# make figure of legend only 
ggplot(t_df, aes(x=Time, y=Area, color=side, group=side))+
  geom_point() +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  geom_smooth(se = F, method = "loess", linewidth = 1.5) +
  scale_color_manual(
    values = c("api" = "lightblue", "baso" = "lightpink"),
    labels = c("api" = "Apical", "baso" = "Basolateral")
  ) +
  cell_plot_theme +
  labs(color = "Side") # Ensure legend title is set
ggsave("legend.tiff", dpi = 300, width = 3.5, height = 3, units = "in") # Cell guidelines recommend these dimensions

  

# Default line plot
p2 <-ggplot(t_df, aes(x=Time, y=Area, color=side, group=side))+
  geom_point() +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  geom_smooth(se = F, method = "loess", linewidth = 1.5) +
  #stat_summary(fun.data ='mean_se', geom = "smooth", se = TRUE, alpha=0.1) +
  labs(#title= "Azilsartan and metabolite", 
       #subtitle="Timecourse", 
       x = "Time (minutes)", 
       y = "Peak area Intensity") + 
  theme_minimal() + 
  theme(legend.position = "none", 
        strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1), lineheight = 0.8), # Adjust facet label size and margin
        panel.spacing = unit(0.5, "lines") # Reduce space between panels
  ) +
  scale_color_manual(values = c("api" = "lightblue", "baso" = "lightpink")) +
  scale_y_continuous(labels = scales::scientific, breaks = c(0, 1.5e4, 3e4, 6e5, 8e5)) + # Scientific notation for Y-axis
  scale_x_continuous(breaks = c(0, 120, 240)) +    # Specific ticks for X-axis
  facet_grid(
    side ~ drug + condition, 
    labeller = labeller(
      drug = drug_labs, 
      condition = condition_labs, 
      side = function(labels) rep("", length(labels))
    ), 
    scales = "free_y", 
    switch = "y"
  ) +
  cell_plot_theme +
  geom_blank(data = empty_panels)

print(p2)
# Save the plot
ggsave("azilsartan.tiff", plot = p2, dpi = 300, width = 3.25, height = 3, units = "in") # Cell guidelines recommend these dimensions

library(cowplot)
p1 = p1 + 
  theme(plot.margin = margin(0, 0, 0, 0))
p2 = p2 + 
  theme(plot.margin = margin(0, 0, 0, 0))
p_1 = plot_grid(p1, p2)
ggsave("fig3A-B.tiff", plot = p_1, dpi = 300, width = 6.5, height = 2.5, units = "in") # Cell guidelines recommend these dimensions


#### Carfilzomib ####
unique(my_prep_df$drug)
#library(ggbreak)

d = "drug_carfilzomib.Results_Area"
m = "met_carfilzomib1.Results_Area"
m1 = "met_carfilzomib2.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m | drug == m1) %>% 
  filter(Area > 5000) %>% 
  na.omit()

t_df = temp_df %>% 
  filter(drug == d & condition == "TP0" | 
           drug == m1 & condition %in% c("TP7", "TP8"))

# Smooth the data by averaging the value at each time point with its two flanking points
t_df <- t_df %>%
  group_by(side, drug, strain) %>%
  arrange(Time) %>%
  mutate(smooth_value = (Area + lag(Area, 1, order_by = Time) + lead(Area, 1, order_by = Time)) / 3) %>%
  ungroup()

drug_labs = c("drug_carfilzomib.Results_Area" = "Carfilzomib", 
              "met_carfilzomib2.Results_Area" = "Bacterial metabolite")
condition_labs = c("TP0" = "0 hour bacterial incubation", 
                   "TP8"  = "8 hour bacterial incubation")

empty_panels <- data.frame(
  Time = c(0, 0), # Use valid numeric values
  Area = c(NA_real_, NA_real_), # Ensure Area is numeric
  side = c("aa", "aa")) # Character column

# Default line plot
p <- t_df %>% 
  filter(strain %in% c("Bu")) %>% 
  filter(!(Time == 120 & drug == d)) %>% 
  add_row(
    side = "baso", 
    condition = "TP0", 
    drug = "drug_carfilzomib.Results_Area", 
    Time = NA, 
    Area = NA, 
    strain = "Bu"
  ) %>% 
  ggplot(., aes(x=Time, y=Area, color=side, group=side))+
  geom_smooth(method = "loess", linewidth = 1.5, se = F) +
  geom_point() +
  #stat_summary(fun.data ='mean_se', geom = "smooth", se = TRUE, alpha=0.1) +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  labs(#title= "Carfilzomib and metabolite", 
       #subtitle="Timecourse", 
       x = "Time (minutes)", 
       y = "Peak area Intensity") + 
  theme_minimal() + 
  facet_grid(side~drug+condition, 
             labeller = labeller(drug = drug_labs, condition = condition_labs, 
                                 side = function(labels) rep("", length(labels))), 
             scales = "free_y") +
  theme(legend.position = "none",
    strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1), lineheight = 0.8), # Adjust facet label size and margin
    panel.spacing = unit(0.5, "lines")) +
  scale_color_manual(values = c("api" = "lightblue", "baso" = "lightpink")) +
  cell_plot_theme + 
  scale_y_continuous(labels = scales::scientific, breaks = c(0, 1e4, 2e4, 3e5, 6e5)) + # Scientific notation for Y-axis
  scale_x_continuous(breaks = c(0, 120, 240)) +
  geom_blank(data = empty_panels)

print(p)
ggsave("carfilzomib.tiff", plot = p, dpi = 300, width = 3.25, height = 3, units = "in") # Cell guidelines recommend these dimensions



#### Nicergoline ####
unique(my_prep_df$drug)

d = "drug_nicergoline.Results_Area"
m = "met_nicergoline.Results_Area"


temp_df = my_prep_df %>% 
  filter(drug == d | drug == m) %>% 
  filter(Area > 5500) %>% 
  na.omit()

t_df = temp_df %>% 
  filter(condition == "TP8" & strain %in% c("nobug", "Bu"))

# Smooth the data by averaging the value at each time point with its two flanking points
t_df <- t_df %>%
  group_by(side, drug, strain) %>%
  arrange(Time) %>%
  mutate(smooth_value = (Area + lag(Area, 1, order_by = Time) + lead(Area, 1, order_by = Time)) / 3) %>%
  ungroup()

drug_labs = c("drug_nicergoline.Results_Area" = "Nicergoline", 
              "met_nicergoline.Results_Area" = "Bacterial metabolite")
condition_labs = c("TP0" = "0 hour bacterial incubation", 
                   "TP8"  = "8 hour bacterial incubation")
strain_labs = c("Bu" = "B. uniformis", "nobug" = "Sterile control")

# Define empty panels with proper data types
empty_panels <- data.frame(
  Time = c(0, 0), # Use valid numeric values
  Area = c(NA_real_, NA_real_), # Ensure Area is numeric
  side = c("api", "baso"))# Character column

# Default line plot
p_n <- t_df %>% 
  filter(drug == "drug_nicergoline.Results_Area") %>% 
  #filter(!(strain == "nobug" &) %>% 
  ggplot(., aes(x=Time, y=Area, color=side, group=side))+
  geom_point() +
  #stat_summary(fun.data ='mean_se', geom = "smooth", se = TRUE, alpha=0.1) +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  geom_smooth(method = "loess", linewidth = 1.5, se = F) +
  labs(#title= "Nicergoline and metabolite", 
       #subtitle="Timecourse", 
       x = "Time (minutes)", 
       y = "Peak area Intensity") + 
  theme_minimal() + 
  #facet_wrap(~condition+drug+strain) + 
  facet_grid(side~drug+strain,
             labeller = labeller(drug = drug_labs, condition = condition_labs, 
                                 strain = strain_labs, side = function(labels) rep("", length(labels))), 
            scales = "free") + 
  #scale_y_break(breaks = c(1e5, 2e5)) +
  #geom_smooth(se = F, method = "auto", linewidth = 1.5) +
  scale_color_manual(values = c("api" = "lightblue", "baso" = "lightpink")) +
  theme(legend.position = "none",
    strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1), lineheight = 0.8), # Adjust facet label size and margin
    panel.spacing = unit(0.5, "lines")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific, n.breaks = 3) +
  cell_plot_theme + 
  geom_blank(data = empty_panels)
print(p_n)
ggsave("nicergoline.tiff", plot = p_n, dpi = 300, width = 3, height = 2, units = "in") # Cell guidelines recommend these dimensions


p_m <- t_df %>% 
  filter(drug == "met_nicergoline.Results_Area") %>% 
  #filter(!(strain == "nobug" &) %>% 
  ggplot(., aes(x=Time, y=Area, color=side, group=side))+
  geom_point() +
  #stat_summary(fun.data ='mean_se', geom = "smooth", se = TRUE, alpha=0.1) +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  geom_smooth(method = "loess", linewidth = 1.5, se = F) +
  labs(#title= "Nicergoline and metabolite", 
    #subtitle="Timecourse", 
    x = "Time (minutes)", 
    y = "Peak area Intensity") + 
  theme_minimal() + 
  #facet_wrap(~condition+drug+strain) + 
  facet_grid(side~drug+strain,
             labeller = labeller(drug = drug_labs, condition = condition_labs, 
                                 strain = strain_labs, side = function(labels) rep("", length(labels))), 
             scales = "free") + 
  #scale_y_break(breaks = c(1e5, 2e5)) +
  #geom_smooth(se = F, method = "auto", linewidth = 1.5) +
  scale_color_manual(values = c("api" = "lightblue", "baso" = "lightpink")) +
  theme(legend.position = "none",
        strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1), lineheight = 0.8), # Adjust facet label size and margin
        panel.spacing = unit(0.5, "lines")) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  scale_y_continuous(labels = scales::scientific, n.breaks = 3) +
  cell_plot_theme + 
  geom_blank(data = empty_panels)
print(p_m)
ggsave("nicergoline_met.tiff", plot = p_m, dpi = 300, width = 3, height = 2, units = "in") # Cell guidelines recommend these dimensions



