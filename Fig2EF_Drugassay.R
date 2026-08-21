## Set up 
######Load Libraries########
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(ggforce)

#### Data import and cleaning ####
data_list = file.path("input_folder", "Data_AB012I_drugassay3.csv")

#read them
screen <- lapply(data_list, read.csv)

# run the following data cleaning functions for the screen list 
# extract the targeted drugs names (+ remove spaces)
targeted_drugs = function(df){
  gsub(" ", "_", gsub(" Results","",grep("Results",colnames(df),fixed=T,value=T),fixed = T), fixed = T)
}

target_drugs = lapply(screen, targeted_drugs)

# update column names
colnms <- list()
for (i in seq_along(target_drugs)) { 
  colnms[[i]] <- as.character(c("Name", 
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

# get cleaned df 
merged_df <- bind_rows(screen, .id = "column_label")

# Add columns with timepoints & conditions 
merged_df$group = sapply(strsplit(as.character(merged_df$Name), "_"), '[', 2)
merged_df$condition = sapply(strsplit(as.character(merged_df$Name), "_"), '[', 3)
merged_df$TPs = sapply(strsplit(as.character(merged_df$Name), "_"), '[', 4)
merged_df$well = sapply(strsplit(as.character(merged_df$Name), "_"), '[', 5)
merged_df$rep = sapply(strsplit(as.character(merged_df$Name), "_"), '[', 6)

# plot all for overview 
# Change area variables to numeric
cols = grep("_Area", colnames(merged_df))
merged_df[, cols] = apply(merged_df[,cols], 2, function(x) as.numeric(as.character(x)))

drug_cols = grep("(drug|met).*Area", colnames(merged_df))

merged_df = merged_df %>% 
  dplyr::filter(TPs != "all")

pdf("./Plots/drugassay_raw.pdf")
for (i in drug_cols){
  
  # Default line plot
  p <-  ggplot(merged_df, aes_string(x="TPs", y=colnames(merged_df)[i], shape="rep", color = "condition")) +
    geom_point() +
    labs(title= colnames(merged_df)[i], 
         subtitle="Timecourse", 
         x = "Time (hours)", 
         y = "Peak area Intensity") + 
    theme_minimal() + 
    facet_wrap(~group, scales = "free") + 
    geom_blank()
  
  print(p)
  
}

dev.off()


#### Restructure df ####
# pivot merged df 
area_cols = grep("_Area$", colnames(merged_df), value = T)
my_prep_df = merged_df[, c("Name", "condition", "TPs", "rep", "well", "group", area_cols)]
my_prep_df = pivot_longer(my_prep_df, 
                          cols = grep("_Area$", colnames(my_prep_df), value = T), 
                          values_to = "Area", 
                          names_to = "drug")

my_prep_df$TPs = as.character(gsub("TP", "", my_prep_df$TPs))


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


#### Calculate thresholds ####
# define threshold: calculate mean of each ion of the DMSO samples and the drug-GMM samples 
# threshold = mean + 2 stdev 

# calculate thresholds
my_prep_df <- my_prep_df %>%
  group_by(drug) %>%
  mutate(
    mean_DMSO = median(Area[condition == "DMSO"], na.rm = TRUE), 
    sd2_DMSO   = 2 * sd(Area[condition == "DMSO"], na.rm = TRUE),
  )

my_prep_df$threshold = my_prep_df$mean_DMSO + my_prep_df$sd2_DMSO

# add area with thesholds 
my_prep_df$corr_area = ifelse(my_prep_df$Area < my_prep_df$threshold, NA, my_prep_df$Area)


#### make boxplot-pairs ####


### roxatidine ####
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

# Default line plot
p <- temp_df %>% 
  filter(condition != "DMSO") %>% 
  ggplot(., aes(x=TPs, y=Area, color = group))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter()) +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  labs(x = "Time (hours)", 
    y = "Concentration (a.u.)") + 
  theme_minimal() + 
  #facet_wrap(~condition+drug+strain) + 
  facet_grid(drug~group, scales = "free_y", 
             labeller = labeller(drug = drug_labs, group = strain_labs)) +
  scale_y_continuous(labels = scales::scientific) +
  #theme(legend.position = "none") +
  scale_color_manual(values = c("#DDA0DD", "#5F9EA0", "#E75480","#F08080", "#FFD1DC", "#FFA07A"), 
                     labels = strain_labs) +
  theme(
    legend.position = "top",            # place legend above
    legend.direction = "horizontal",    # make it horizontal
  ) +
  cell_plot_theme 
print(p)

ggsave("./Plots/drugassays/roxace.tiff", plot = p, dpi = 600, width = 3.5, height = 2.5, units = "in")


### Nicergoline ####
colnames(merged_df)
d = "drug_nicergoline.Results_Area"
m = "met_nicergoline.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m)


drug_labs = c("drug_nicergoline.Results_Area" = "Nicergoline", 
              "met_nicergoline.Results_Area" = "Nicergoline metabolite")
strain_labs = c("GMM" = "Sterile control", 
                "MB002" = "donor 1", 
                "MB003" = "donor 2", 
                "MB005" = "donor 3", 
                "mouse" = "mouse")

# Default line plot
p <- temp_df %>% 
  filter(condition != "DMSO") %>% 
  ggplot(., aes(x=TPs, y=Area, color = group))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter()) +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  labs(x = "Time (hours)", 
       y = "Concentration (a.u.)") + 
  theme_minimal() + 
  #facet_wrap(~condition+drug+strain) + 
  facet_grid(drug~group, scales = "free_y", 
             labeller = labeller(drug = drug_labs, group = strain_labs)) +
  scale_y_continuous(labels = scales::scientific, n.breaks = 4) +
  #theme(legend.position = "none") +
  scale_color_manual(values = c("#DDA0DD", "#5F9EA0", "#E75480","#F08080", "#FFD1DC", "#FFA07A"), 
                     labels = strain_labs) +
  # theme(
  #   legend.position = "top",            # place legend above
  #   legend.direction = "horizontal",    # make it horizontal
  # ) +
  theme(legend.position = "none") +
  cell_plot_theme 
print(p)

ggsave("./Plots/drugassays/nicergoline_nolegend.tiff", plot = p, dpi = 600, width = 3.5, height = 3, units = "in")


### bisacodyl ####
colnames(merged_df)
d = "drug_bisacodyl.Results_Area"
m = "met_bisacodyl278.Results_Area"
m1 = "met_bisacodyl320.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m | drug == m1)


drug_labs = c("drug_bisacodyl.Results_Area" = "Bisacodyl", 
              "met_bisacodyl278.Results_Area" = "Metabolite (278 m/z)", 
              "met_bisacodyl320.Results_Area" = "Metabolite (320 m/z)")
strain_labs = c("GMM" = "Sterile control", 
                "MB002" = "donor 1", 
                "MB003" = "donor 2", 
                "MB005" = "donor 3", 
                "mouse" = "mouse")

# Default line plot
p <- temp_df %>% 
  filter(condition != "DMSO") %>% 
  ggplot(., aes(x=TPs, y=Area, color = group))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter()) +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  labs(x = "Time (hours)", 
       y = "Concentration (a.u.)") + 
  theme_minimal() + 
  facet_grid(drug~group, scales = "free", 
             labeller = labeller(drug = drug_labs, group = strain_labs)) +
  scale_y_continuous(labels = scales::scientific) +
  #theme(legend.position = "none") +
  scale_color_manual(values = c("#DDA0DD", "#5F9EA0", "#E75480","#F08080", "#FFD1DC", "#FFA07A"), 
                     labels = strain_labs) +
  theme(
    legend.position = "top",            # place legend above
    legend.direction = "horizontal",    # make it horizontal
  ) +
  cell_plot_theme 
print(p)

ggsave("./Plots/drugassays/bisacodyl.tiff", plot = p, dpi = 600, width = 3.5, height = 2.5, units = "in")

### mmf ####
colnames(merged_df)
d = "drug_mycophenolatemofetil.Results_Area"
m = "met_mmf.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m) 


drug_labs = c("drug_mycophenolatemofetil.Results_Area" = "Mycophenolate Mofetil", 
              "met_mmf.Results_Area" = "Mycophenolic acid")
strain_labs = c("GMM" = "Sterile control", 
                "MB002" = "donor 1", 
                "MB003" = "donor 2", 
                "MB005" = "donor 3", 
                "mouse" = "mouse")

# Default line plot
p <- temp_df %>% 
  filter(condition != "DMSO") %>% 
  ggplot(., aes(x=TPs, y=Area, color = group))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter()) +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  labs(x = "Time (hours)", 
       y = "Concentration (a.u.)") + 
  theme_minimal() + 
  facet_grid(drug~group, scales = "free_y", 
             labeller = labeller(drug = drug_labs, group = strain_labs)) +
  scale_y_continuous(labels = scales::scientific, n.breaks = 3) +
  #theme(legend.position = "none") +
  scale_color_manual(values = c("#DDA0DD", "#5F9EA0", "#E75480","#F08080", "#FFD1DC", "#FFA07A"), 
                     labels = strain_labs) +
  theme(legend.position = "none") +
  cell_plot_theme 
print(p)

ggsave("./Plots/drugassays/mmf.tiff", plot = p, dpi = 600, width = 3.5, height = 3, units = "in")


### Carfilzomib ####
colnames(merged_df)
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
strain_labs = c("GMM" = "Sterile control", 
                "MB002" = "donor 1", 
                "MB003" = "donor 2", 
                "MB005" = "donor 3", 
                "mouse" = "mouse")

# Default line plot
p <- temp_df %>% 
  filter(condition != "DMSO") %>% 
  filter(drug %in% c(d, m3)) %>% 
  ggplot(., aes(x=TPs, y=Area, color = group))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter()) +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  labs(x = "Time (hours)", 
       y = "Concentration (a.u.)") + 
  theme_minimal() + 
  facet_grid(drug~group+condition, scales = "free", 
             labeller = labeller(drug = drug_labs, group = strain_labs)) +
  scale_y_continuous(labels = scales::scientific) +
  #theme(legend.position = "none") +
  scale_color_manual(values = c("#DDA0DD", "#5F9EA0", "#E75480","#F08080", "#FFD1DC", "#FFA07A"), 
                     labels = strain_labs) +
  theme(
    legend.position = "top",            # place legend above
    legend.direction = "horizontal",    # make it horizontal
  ) +
  cell_plot_theme 
print(p)

ggsave("./Plots/drugassays/carfilzomib_main.tiff", plot = p, dpi = 600, width = 3.5, height = 2.5, units = "in")
#ggsave("./Plots/drugassays/carfilzomib_all.tiff", plot = p, dpi = 600, width = 3.5, height = 4.5, units = "in")


### Methylprednisolone ####
colnames(merged_df)
d = "drug_methylprednisolone.Results_Area"
m = "met_methylprednisolone1.Results_Area"
m1 = "met_methylprednisolone2.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug %in% c(m, m1, d)) 


drug_labs = c("drug_methylprednisolone.Results_Area" = "Methylprednisolone", 
              "met_methylprednisolone1.Results_Area" = "Methylprednisolone M1", 
              "met_methylprednisolone2.Results_Area" = "Methylprednisolone M2")
strain_labs = c("GMM" = "Sterile control", 
                "MB002" = "donor 1", 
                "MB003" = "donor 2", 
                "MB005" = "donor 3", 
                "mouse" = "mouse")

# Default line plot
p <- temp_df %>% 
  filter(condition != "DMSO") %>% 
  ggplot(., aes(x=TPs, y=corr_area, color = group))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter()) +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  labs(x = "Time (hours)", 
       y = "Concentration (a.u.)") + 
  theme_minimal() + 
  facet_grid(drug~group+condition, scales = "free", 
             labeller = labeller(drug = drug_labs, group = strain_labs)) +
  scale_y_continuous(labels = scales::scientific) +
  #theme(legend.position = "none") +
  scale_color_manual(values = c("#DDA0DD", "#5F9EA0", "#E75480","#F08080", "#FFD1DC", "#FFA07A")) +
  cell_plot_theme 
print(p)



### Deflazacort ####
colnames(merged_df)
d = "drug_deflazacort.Results_Area"
m = "met_deflazacort.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m) 


drug_labs = c("drug_deflazacort.Results_Area" = "Deflazacort", 
              "met_deflazacort.Results_Area" = "Deflazacort metabolite")
strain_labs = c("GMM" = "Sterile control", 
                "MB002" = "donor 1", 
                "MB003" = "donor 2", 
                "MB005" = "donor 3", 
                "mouse" = "mouse")

# Default line plot
p <- temp_df %>% 
  filter(condition != "DMSO") %>% 
  ggplot(., aes(x=TPs, y=Area, color = group))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter()) +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  labs(x = "Time (hours)", 
       y = "Concentration (a.u.)") + 
  theme_minimal() + 
  facet_grid(drug~group+condition, scales = "free", 
             labeller = labeller(drug = drug_labs, group = strain_labs)) +
  scale_y_continuous(labels = scales::scientific) +
  #theme(legend.position = "none") +
  scale_color_manual(values = c("#DDA0DD", "#5F9EA0", "#E75480","#F08080", "#FFD1DC", "#FFA07A")) +
  theme(
    legend.position = "top",            # place legend above
    legend.direction = "horizontal",    # make it horizontal
  ) +
  cell_plot_theme 
print(p)

ggsave("./Plots/drugassays/deflazacort.tiff", plot = p, dpi = 600, width = 3.5, height = 4.5, units = "in")



### Racecadotril ####
colnames(merged_df)
d = "drug_racecadotril.Results_Area"
m = "met_racecadotril.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m)


drug_labs = c("drug_racecadotril.Results_Area" = "Racecadotril", 
              "met_racecadotril.Results_Area" = "Racecadotril metabolite")
strain_labs = c("GMM" = "Sterile control", 
                "MB002" = "donor 1", 
                "MB003" = "donor 2", 
                "MB005" = "donor 3", 
                "mouse" = "mouse")

# Default line plot
p <- temp_df %>% 
  filter(condition != "DMSO") %>% 
  ggplot(., aes(x=TPs, y=corr_area, color = group))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter()) +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  labs(x = "Time (hours)", 
       y = "Concentration (a.u.)") + 
  theme_minimal() + 
  facet_grid(drug~group+condition, scales = "free", 
             labeller = labeller(drug = drug_labs, group = strain_labs)) +
  scale_y_continuous(labels = scales::scientific) +
  #theme(legend.position = "none") +
  scale_color_manual(values = c("#DDA0DD", "#5F9EA0", "#E75480","#F08080", "#FFD1DC", "#FFA07A")) +
  cell_plot_theme 
print(p)

### Ataluren ####
colnames(merged_df)
d = "drug_ataluren.Results_Area"
m = "met_ataluren.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m)


drug_labs = c("drug_ataluren.Results_Area" = "Ataluren", 
              "met_ataluren.Results_Area" = "Ataluren metabolite")
strain_labs = c("GMM" = "Sterile control", 
                "MB002" = "donor 1", 
                "MB003" = "donor 2", 
                "MB005" = "donor 3", 
                "mouse" = "mouse")

# Default line plot
p <- temp_df %>% 
  filter(condition != "DMSO") %>% 
  ggplot(., aes(x=TPs, y=corr_area, color = group))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter()) +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  labs(x = "Time (hours)", 
       y = "Concentration (a.u.)") + 
  theme_minimal() + 
  facet_grid(drug~group+condition, scales = "free", 
             labeller = labeller(drug = drug_labs, group = strain_labs)) +
  scale_y_continuous(labels = scales::scientific) +
  #theme(legend.position = "none") +
  scale_color_manual(values = c("#DDA0DD", "#5F9EA0", "#E75480","#F08080", "#FFD1DC", "#FFA07A")) +
  cell_plot_theme 
print(p)


### Dexamethasone ####
colnames(merged_df)
d = "drug_dexamethasone.Results_Area"
m = "met_dexamethasone.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m)


drug_labs = c("drug_dexamethasone.Results_Area" = "Dexamethasone", 
              "met_dexamethasone.Results_Area" = "Dexamethasone metabolite")
strain_labs = c("GMM" = "Sterile control", 
                "MB002" = "donor 1", 
                "MB003" = "donor 2", 
                "MB005" = "donor 3", 
                "mouse" = "mouse")

# Default line plot
p <- temp_df %>% 
  filter(condition != "DMSO") %>% 
  ggplot(., aes(x=TPs, y=Area, color = group))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter()) +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  labs(x = "Time (hours)", 
       y = "Concentration (a.u.)") + 
  theme_minimal() + 
  facet_grid(drug~group+condition, scales = "free", 
             labeller = labeller(drug = drug_labs, group = strain_labs)) +
  scale_y_continuous(labels = scales::scientific) +
  #theme(legend.position = "none") +
  scale_color_manual(values = c("#DDA0DD", "#5F9EA0", "#E75480","#F08080", "#FFD1DC", "#FFA07A"), 
                     labels = strain_labs) +
  theme(
    legend.position = "top",            # place legend above
    legend.direction = "horizontal",    # make it horizontal
  ) +
  cell_plot_theme 
print(p)

ggsave("./Plots/drugassays/dexamethasone.tiff", plot = p, dpi = 600, width = 3.5, height = 4.5, units = "in")



### Lenvatinib ####



### Azilsartan ####
colnames(merged_df)
d = "drug_azilsartan.Results_Area"
m = "met_azilsartan.Results_Area"

temp_df = my_prep_df %>% 
  filter(drug == d | drug == m)


drug_labs = c("drug_azilsartan.Results_Area" = "Azilsartan", 
              "met_azilsartan.Results_Area" = "Azilsartan metabolite")
strain_labs = c("GMM" = "Sterile control", 
                "MB002" = "donor 1", 
                "MB003" = "donor 2", 
                "MB005" = "donor 3", 
                "mouse" = "mouse")

# Default line plot
p <- temp_df %>% 
  filter(condition != "DMSO") %>% 
  ggplot(., aes(x=TPs, y=Area, color = group))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter()) +
  #stat_summary(aes(y = Area, group = side), fun=mean, geom="line") +
  labs(x = "Time (hours)", 
       y = "Concentration (a.u.)") + 
  theme_minimal() + 
  facet_grid(drug~group+condition, scales = "free", 
             labeller = labeller(drug = drug_labs, group = strain_labs)) +
  scale_y_continuous(labels = scales::scientific) +
  #theme(legend.position = "none") +
  scale_color_manual(values = c("#DDA0DD", "#5F9EA0", "#E75480","#F08080", "#FFD1DC", "#FFA07A"), 
                     labels = strain_labs) +
  theme(
    legend.position = "top",            # place legend above
    legend.direction = "horizontal",    # make it horizontal
  ) +
  cell_plot_theme 
print(p)

ggsave("./Plots/drugassays/azilsartan.tiff", plot = p, dpi = 600, width = 3.5, height = 2.5, units = "in")


#### Make heatmaps ####
# calculate the fold change 
#spliting the df based on time
my_prep_df = my_prep_df %>% 
  filter(grepl("met_|drug_", drug))

df1 =my_prep_df%>%filter(TPs == "TP00")
df2 =my_prep_df%>%filter(TPs == "TP16")

df1 = df1 %>% select(Name, drug, condition, group, Area)
df2 = df2 %>% select(Name, drug, condition, group, Area)

#renaming the column
names(df1)[names(df1) == "Area"] <- "00hr"
names(df2)[names(df2) == "Area"] <- "16hr"

#Joining two df
df3 = full_join(df1, df2, by = c("drug","condition","group"))
# exclude combos where I'm missing one of the two timepoints
df3 = df3[!is.na("Name.x") & !is.na("Name.y"),]


#computing the FC based on t=0hr
df_FC = df3 %>% group_by(drug, group, condition) %>% mutate(FC=log2(mean(`16hr`, na.rm = T)/mean(`00hr`, na.rm = T)))

# compute the pvalue associated to FC on t=0hr
df_FC = df_FC %>% group_by(drug, group, condition) %>% mutate(FC_pval = tryCatch({
  t.test(`2hr`,`0hr`, alternative = "less")$p.value
},
error=function(cond) {
  return(1)
})
)

# transform to percentage 
df_FC = df_FC %>% 
  mutate(percentage = (100-((2^(FC))*100))) 


### heatmap drugs ####
library(ComplexHeatmap)

results_final = df_FC %>% 
  dplyr::select(c(drug, group, FC, percentage, condition)) %>% 
  filter(grepl("drug_", drug)) %>% 
  filter(condition != "DMSO") %>% 
  mutate(percentage = if_else(percentage < 0, 0, percentage)) %>% 
  distinct() %>% 
  ungroup()

results_final$percentage [is.na(results_final$percentage) | results_final$percentage < 0] <- 0

results_mat = results_final %>% 
  select(drug, group, percentage) %>% 
  pivot_wider(names_from = group, 
              values_from = percentage, 
              values_fill = 0)



drugs = results_mat$drug

results_mat = as.matrix(results_mat[2:ncol(results_mat)])
rownames(results_mat) = drugs

colnames(results_mat)

Heatmap(results_mat, name = "% of depletion", 
        col = colorRampPalette(c("black", "white"))(100), 
        row_title = "Drugs tested", 
        column_title = "Microbiome sample", 
        row_names_side = "left",            
        row_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(10, "cm"),
        heatmap_legend_param = 
          list(title_gp = gpar(fontsize = 12, lineheight = 1.5)), 
        column_title_gp = 
          gpar(fontsize = 12, lineheight = 1.2), 
        cluster_columns = F, 
        cluster_rows = T, 
        
)


### heatmap metabolites ####
library(ComplexHeatmap)

results_final = df_FC %>% 
  dplyr::select(c(drug, group, FC, percentage, condition)) %>% 
  filter(grepl("met", drug)) %>% 
  filter(condition != "DMSO") %>% 
  mutate(percentage = if_else(percentage > 0, 0, percentage)) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(percentage = log10(abs(percentage)))

results_final$percentage [is.na(results_final$percentage) | results_final$percentage < 0] <- 0

results_mat = results_final %>% 
  select(drug, group, percentage) %>% 
  pivot_wider(names_from = group, 
              values_from = percentage, 
              values_fill = 0)



drugs = results_mat$drug

results_mat = as.matrix(results_mat[2:ncol(results_mat)])
rownames(results_mat) = drugs

colnames(results_mat)

Heatmap(results_mat, name = "log10(% of production)", 
        col = colorRampPalette(c("black", "white"))(100), 
        row_title = "Drugs tested", 
        column_title = "Microbiome sample", 
        row_names_side = "left",            
        row_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(10, "cm"),
        heatmap_legend_param = 
          list(title_gp = gpar(fontsize = 12, lineheight = 1.5)), 
        column_title_gp = 
          gpar(fontsize = 12, lineheight = 1.2), 
        cluster_columns = F, 
        cluster_rows = T, 
        
)
