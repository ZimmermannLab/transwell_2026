### Script for DMSO pools api baso ###

library(tidyverse)

cell_plot_theme <- theme(
  text = element_text(family = "Arial", size = 10), # Set Arial font and text size
  axis.text = element_text(size = 8), # Adjust axis text size
  axis.title = element_text(size = 10), # Adjust axis title size
  legend.text = element_text(size = 10), # Adjust legend text size
  legend.title = element_text(size = 10), # Adjust legend title size
  plot.title = element_text(size = 12, hjust = 0.5), # Center title
  panel.background = element_blank(), # Remove background color
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Add border
  panel.grid.major = element_blank(), # Subtle grid
  panel.grid.minor = element_blank(), # No minor grid
  plot.margin = margin(5, 5, 5, 5) # Adjust margins
)


### Read df ####
#load("Afterpeakdetection_dmsopools.RData")
df = read.csv(file.path("input_folder", "mzmine_dmso_api_baso__quant.csv"))

# tidy data 
df$feature_ID = paste0("FT", df$row.ID)
# separate the peak_info from the areas 
mzmine_peak = df[, c(1:13)]
mzmine_peak$feature_ID = df$feature_ID
mzmine_area = df %>% 
  dplyr::select(contains("pos"))
mzmine_area = cbind(df$feature_ID, mzmine_area)

### Annotation ####
# do annotation with MiMeDB 
hydrogen = 1.007825

# calculate only for which best ion is empty 
neutral_mass = function(df){ 
  
  df$neutral.M.mass[df$best.ion == ""] = df$row.m.z + hydrogen
  return(df)
}
mzmine_peak = neutral_mass(mzmine_peak)

# read in MiMeDB database
library(data.table)
MiMeDB = fread(file.path("input_folder", "mimedb_metabolites_v20240319.csv"))
MiMeDB$moldb_mono_mass = as.numeric(MiMeDB$moldb_mono_mass)
# change NA in exact mass to 0 
MiMeDB$moldb_mono_mass[is.na(MiMeDB$moldb_mono_mass)] <- 0

# define function to annotate peaks based on mass 
add_annot_MiMeDB = function(peaks){ 
  
  peaks$annot_mime = NA
  
  for(j in 1:nrow(peaks)){ 
    
    cur_mass = peaks$neutral.M.mass[j]
    
    mass_shift = abs(cur_mass - MiMeDB$moldb_mono_mass)
    
    # set filter on MASS; absolute (2mDA) or relative (20ppm)  
    pos = (mass_shift <= 0.002 | mass_shift <= (20*10^-6*cur_mass))
    
    if(sum(pos)==1){
      
      # we actually have an overlap!
      peaks$annot_mime[j] = MiMeDB$mime_id[pos]
      
    }else if(sum(pos, na.rm = T)>1){
      
      peaks$annot_mime[j] = paste0(unique(MiMeDB$mime_id[pos]), collapse = "; ")
      
    }
  } 
  
  return(peaks)
  
} 

mzmine_peak = add_annot_MiMeDB(mzmine_peak)
# how many ions have an annot? 
sum(!is.na(mzmine_peak$annot_mime))

### add annotations to ion df ###
peak_min = mzmine_peak %>% 
  dplyr::select(feature_ID, row.ID, row.m.z, row.retention.time, neutral.M.mass, annot_mime)
peak_min = na.omit(peak_min)
mzmine_area_annot = left_join(mzmine_area, peak_min, by=c("df$feature_ID" = "feature_ID"))
sum(!is.na(mzmine_area_annot$annot_mime))

# keep only the ones with an annotation 
colnames(mzmine_area)[1] <- "feature_ID"
mzmine_area_annot = mzmine_area %>% 
  dplyr::filter(feature_ID %in% peak_min$feature_ID)



### Stat test: strain against nobug ####
#### api only ####
source("compare_tonnobug_mzmine.R")
cols_to_delete = c("sample_name", "time", "pool", "strain")

api_df = mzmine_area_annot %>% 
  dplyr::select(c(feature_ID, grep("api", colnames(mzmine_area_annot))))
results = compare_tonobug_mzmine(complete_df = api_df, cols_to_delete = cols_to_delete)
api_results = results

# bind results together
results = lapply(results, function(x) as.data.frame(x))
results =bind_rows(results, .id = "strain")


# look at p-val distribution 
ggplot(results, aes(x=padj)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~strain)

range(results$pval, na.rm = T)

# plot volcano plots 
ggplot(results, aes(x=FC, y=-log10(padj))) + 
  geom_point() + 
  facet_wrap(~strain)

# select positive FC features
# threshold of >5 (produced metabolites)
selection = results %>% 
  filter(FC > 5 & padj < 0.05)
unique(selection$strain)
table(selection$strain)

# how many unique ions? 
length(unique(selection$ion))

# get peak info on ions 
ions_info_nobug = left_join(selection, 
                            peak_min, 
                            by = c("ion" = "feature_ID"))
table(ions_info_nobug$strain) 

# change RT from minutes to seconds 
ions_info_nobug$row.retention.time = ions_info_nobug$row.retention.time * 60

all_ions = ions_info_nobug
#all_ions = rbind(ions_info, ions_info_nobug)
length(unique(all_ions$ion))
write.csv(all_ions, "AB007_api_ttest_candidates.csv")

### Fit EICs in AB007 ####
library(mzR)
library(xcms)
library(MSnbase)

AB007_files = list.files(path = "input_folder", full.names = TRUE, recursive = TRUE, pattern = "\\.mzML$")
source("EIC_all_mzmine_function.R")

# get only unique ions 
all_ions = all_ions %>% 
  dplyr::select(c(ion, row.ID, row.m.z, row.retention.time, neutral.M.mass, annot_mime)) %>% 
  distinct()

ion_df = EIC_all_mzmine(df_drugs = all_ions, 
                        file_list = AB007_files, 
                        output_folder = "./EIC_output_mzmine/")
ion_df_api= ion_df
length(unique(ion_df$ion))
write.csv(ion_df, "AB007_untar_api_iondf.csv")
ion_df = read.csv("AB007_untar_api_iondf.csv", row.names = 1)

### Prep df #####
# add metadata 
ion_df$side = sapply(strsplit(as.character(ion_df$sample_name), "_"), '[', 7)
ion_df$time = sapply(strsplit(as.character(ion_df$sample_name), "_"), '[', 4)
ion_df$pool = sapply(strsplit(as.character(ion_df$sample_name), "_"), '[', 6)
ion_df$strain = sapply(strsplit(as.character(ion_df$sample_name), "_"), '[', 2)
ion_df$time = as.numeric(gsub("T", "", ion_df$time))

# delete the non-integrated peaks 
ion_df = na.omit(ion_df)
length(unique(ion_df$ion))

# clean-up df 
# get rid of double peaks 
temp = ion_df %>% 
  group_by(sample_name, ion) %>% 
  dplyr::count(sample_name, sort=T) %>% 
  filter(n >1)
length(unique(temp$sample_name))
hist(temp$n)

# use selection based on higher signal to noise ratio
ion_df = ion_df %>% 
  group_by(sample_name, ion) %>% 
  filter(sn == max(sn))

selection = ion_df

# check for more than 3 replicates in more than 2 time points 
selection = selection %>% 
  group_by(ion, time, strain, side) %>%
  mutate(reps = n_distinct(pool)) %>% 
  filter(reps >=3) %>% 
  group_by(ion, strain) %>% 
  mutate(tps = n_distinct(time)) %>% 
  filter(tps >= 2) %>% 
  ungroup()

# how many samples and ions?
length(unique(selection$ion))
length(unique(ion_df$ion))
unique(selection$strain)

### Filter out ions that are also in nobug ####
# how many ions also have a signal in nobug? 
nobug_features = selection %>% 
  dplyr::filter(strain == "nobug")

length(unique(nobug_features$ion))

hist(nobug_features$intb[nobug_features$intb < 1e6])

# delete ion features that also have a signal in nobug (as they would not be bacterially produced)
selection = selection %>% 
  dplyr::filter(!(ion %in% unique(nobug_features$ion)))

# how many features remain? 
length(unique(selection$ion))
# 261 ions are microbe-specific (+ meet mean intensity threshold)

#### save api ions 
api_ions = all_ions[all_ions$ion %in% unique(selection$ion), ]
write.csv(api_ions, "./output/api_ions_info.csv")

# add annot name 
peak_min$first_annot = sapply(strsplit(as.character(peak_min$annot_mime), ";"), '[', 1)
selection = left_join(selection, peak_min[,c(1,7)], by = c("ion" = "feature_ID"))
selection = left_join(selection, MiMeDB[, c(2:3)], by = c("first_annot" = "mime_id"))


### Calculate FCs of each ion per strain #### 
# calculate FC over api side 
api_t0 = selection %>% 
  ungroup() %>% 
  group_by(strain) %>% 
  filter(side == "api") %>% 
  filter(time == min(time)) %>% 
  dplyr::rename(t0 = intb) %>% 
  dplyr::select(strain, side, ion, pool, t0)

api_t2 = selection %>% 
  ungroup() %>% 
  group_by(strain) %>% 
  filter(side == "api") %>% 
  filter(time == max(time)) %>% 
  dplyr::rename(t2 = intb) %>% 
  dplyr::select(strain, side, ion, pool, t2)
# join the two time points
api_FC = inner_join(api_t0, api_t2)
api_FC = api_FC %>% 
  group_by(ion, strain) %>% 
  mutate(api_FC = mean(t2) / mean(t0))
length(unique(api_FC$ion))

# log2 for api_fc 
api_FC$log2FC = log2(api_FC$api_FC)
unique(api_FC$side)

# how many have a negative api FC? 
neg_api = api_FC %>% 
  dplyr::select(strain, ion, log2FC) %>% 
  filter(log2FC < 0) %>% 
  distinct()
length(unique(neg_api$ion))

# which of the negative FC ions can also be measured on baso side? 
baso_ions = selection %>% 
  filter(side == "baso") %>% 
  select(strain, reps, tps, name, ion) %>% 
  distinct()
passing_through = inner_join(neg_api, baso_ions, by = c("strain", "ion"))
length(unique(passing_through$ion))

## plot passing through ions ##
pdf("./Plots/Ab007_untar_passing_ions.pdf")
for(i in unique(passing_through$ion)){ 
  
  t =  selection %>% 
    filter(ion == i)
  
  s = all_ions$strain[all_ions$ion == i]
  
  p = ggplot(t, aes(x=time, y=intb, color=side)) + 
    geom_point() + 
    geom_line(aes(group=interaction(side, pool))) + 
    facet_wrap(~strain) + 
    ggtitle(paste0("feature: ", i, " ", t$name[1])) + 
    labs(subtitle = paste("Hit for:", paste(s, collapse = ",")))
  
  print(p)
  
}
dev.off()


### Which negative FC and do not pass? ###
neg_api = neg_api %>% 
  filter(!(ion %in% unique(passing_through$ion)))
length(unique(neg_api$ion))

pdf("./Plots/Ab007_untar_nonpass_ions.pdf")
for(i in unique(neg_api$ion)){ 
  
  t =  selection %>% 
    filter(ion == i)
  
  s = ions_info_nobug$strain[ions_info_nobug$ion == i]
  
  p = ggplot(t, aes(x=time, y=intb, color=side)) + 
    geom_point() + 
    geom_line(aes(group=interaction(side, pool))) + 
    facet_wrap(~strain) + 
    ggtitle(paste0("feature: ", i, " ", t$name[1])) + 
    labs(subtitle = paste("Hit for:", paste(s, collapse = ",")))
  
  print(p)
  
}
dev.off()

## hits: FT11512 Rimboxo, FT15725 Dethiobiotin, FT2029 Cyclo(D−Arg−L−Pro)
# FT14576 Decarestrictine D, FT2029 Cyclo(D−Arg−L−Pro), FT3430 6−hydroxytryprostatin

### Make heatmap api non pass ####
library(ComplexHeatmap)

length(unique(api_FC$ion))
mat = neg_api %>% 
  filter(strain != "nobug") %>% 
  left_join(., selection[, c(11,15,19)], by = c("ion", "strain")) %>% 
  group_by(ion, strain) %>% 
  mutate(mean_FC = mean(log2FC)) %>% 
  ungroup() %>% 
  dplyr::select(strain, name, mean_FC) %>% 
  distinct()
mat = setDF(mat)

# select only ions that are in at least 3 species 
tmp = mat %>% 
  group_by(name) %>% 
  mutate(no_strain = n_distinct(strain)) %>% 
  filter(no_strain > 2) %>% 
  dplyr::select(-no_strain)

mat = pivot_wider(tmp, values_from = mean_FC, values_fn = {mean})
rn = mat$strain
cn = colnames(mat)


# Calculate the number of non-NA rows for each column
non_na_counts <- colSums(!is.na(mat))
# Reorder columns based on the count of non-NA rows
mat <- mat[, order(non_na_counts, decreasing = TRUE)]


mat = as.matrix(mat[-1])
rownames(mat) = rn

rownames(mat)
rownames(mat) <- c("D. formicigenerans", "A. omnicolens", "E. faecalis",
                    "G. haemolysans", "C. scindens",
                     "A. naeslundi", "C. ramosum",
                   "B. thetaiotaomicron", "C. comes", "B. uniformis")


ht <- Heatmap(
  mat, 
  name = "Log2 FC", 
  #col = colorRampPalette(c("black", "white"))(100),
  col = colorRampPalette(c("#FFB6C1", "#5F9EA0"))(100), 
  row_title = "Bacterial species", 
  column_title = "Metabolites", 
  row_names_gp = gpar(fontsize = 6, fontface = "italic"), # Italic row labels
  column_names_gp = gpar(fontsize = 4), 
  column_names_rot = 45,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6, lineheight = 1.2),
    labels_gp = gpar(fontsize = 6)
  ), 
  column_title_gp = gpar(fontsize = 8, lineheight = 1.2),
  row_title_gp = gpar(fontsize = 8, lineheight = 1.2),
  cluster_columns = FALSE, 
  cluster_rows = FALSE
)

tiff("./Plots/AB007_api_heatmap_fig4_fil.tiff", width = 6.5, height = 3, units = "in", res = 300)
draw(ht)
dev.off()


### T-test for baso samples #####
source("compare_tonnobug_mzmine.R")
cols_to_delete = c("sample_name", "time", "pool", "strain")

baso_df = mzmine_area_annot %>% 
  dplyr::select(c(feature_ID, grep("baso", colnames(mzmine_area_annot))))
baso_results = compare_tonobug_mzmine(complete_df = baso_df, cols_to_delete = cols_to_delete)

results = baso_results
# bind results together
results = lapply(results, function(x) as.data.frame(x))
results =bind_rows(results, .id = "strain")

# look at p-val distribution 
ggplot(results, aes(x=padj)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~strain)

range(results$pval, na.rm = T)

# plot volcano plots 
ggplot(results, aes(x=FC, y=-log10(padj))) + 
  geom_point() + 
  facet_wrap(~strain)

# select positive FC features
# threshold of >5 (produced metabolites)
selection = results %>% 
  filter(FC > 3 & padj < 0.05)
unique(selection$strain)
table(selection$strain)

# how many unique ions? 
length(unique(selection$ion))

# get peak info on ions 
ions_info_nobug = left_join(selection, 
                            peak_min, 
                            by = c("ion" = "feature_ID"))
table(ions_info_nobug$strain) 

# change RT from minutes to seconds 
ions_info_nobug$row.retention.time = ions_info_nobug$row.retention.time * 60

baso_ions = ions_info_nobug
#all_ions = rbind(ions_info, ions_info_nobug)
length(unique(baso_ions$ion))
write.csv(baso_ions, "AB007_baso_ttest_candidates.csv")

integrate_baso = baso_ions %>% 
  dplyr::select(ion, row.ID, row.m.z, row.retention.time, annot_mime) %>% 
  distinct()

### Fit EICs in AB007 ####
library(mzR)
library(xcms)
library(MSnbase)

AB007_files = list.files(path = "input_folder", full.names = TRUE, recursive = TRUE, pattern = "\\.mzML$")
source("EIC_all_mzmine_function.R")

baso_ion_df = EIC_all_mzmine(df_drugs = integrate_baso, 
                        file_list = AB007_files, 
                        output_folder = "./EIC_output_mzmine/")
baso_ion_save = baso_ion_df
write.csv(baso_ion_df, "./AB007_untar_baso_iondf.csv")

baso_ion_df = read.csv("AB007_untar_baso_iondf.csv")

### Prep df #####
# add metadata 
baso_ion_df$side = sapply(strsplit(as.character(baso_ion_df$sample_name), "_"), '[', 7)
baso_ion_df$time = sapply(strsplit(as.character(baso_ion_df$sample_name), "_"), '[', 4)
baso_ion_df$pool = sapply(strsplit(as.character(baso_ion_df$sample_name), "_"), '[', 6)
baso_ion_df$strain = sapply(strsplit(as.character(baso_ion_df$sample_name), "_"), '[', 2)
baso_ion_df$time = as.numeric(gsub("T", "", baso_ion_df$time))

# delete the non-integrated peaks 
baso_ion_df = na.omit(baso_ion_df)
length(unique(baso_ion_df$ion))

# clean-up df 
# get rid of double peaks 
temp = baso_ion_df %>% 
  group_by(sample_name, ion) %>% 
  dplyr::count(sample_name, sort=T) %>% 
  filter(n >1)
length(unique(temp$sample_name))
hist(temp$n)

# use selection based on higher signal to noise ratio
baso_ion_df = baso_ion_df %>% 
  group_by(sample_name, ion) %>% 
  filter(sn == max(sn))

fil_baso = baso_ion_df

# check for more than 3 replicates in more than 2 time points 
fil_baso = fil_baso %>% 
  group_by(ion, time, strain, side) %>%
  mutate(reps = n_distinct(pool)) %>% 
  filter(reps >=3) %>% 
  group_by(ion, strain) %>% 
  mutate(tps = n_distinct(time)) %>% 
  filter(tps >= 2) %>% 
  ungroup()

# how many samples and ions?
length(unique(fil_baso$ion))
length(unique(baso_ion_df$ion))
unique(fil_baso$strain)

# add annot name 
fil_baso = left_join(fil_baso, peak_min[,c(1,7)], by = c("ion" = "feature_ID"))
fil_baso = left_join(fil_baso, MiMeDB[, c(2:3)], by = c("first_annot" = "mime_id"))


### Find increasing features @baso side ####
# calculate FC over baso side 
baso_t0 = fil_baso %>% 
  ungroup() %>% 
  group_by(strain) %>% 
  filter(side == "baso") %>% 
  filter(time == min(time)) %>% 
  dplyr::rename(t0 = intb) %>% 
  dplyr::select(strain, side, ion, pool, t0)

baso_t2 = fil_baso %>% 
  ungroup() %>% 
  group_by(strain) %>% 
  filter(side == "baso") %>% 
  filter(time == max(time)) %>% 
  dplyr::rename(t2 = intb) %>% 
  dplyr::select(strain, side, ion, pool, t2)
# join the two time points
baso_FC = inner_join(baso_t0, baso_t2)

# Calculate the log2FC
baso_FC = baso_FC %>% 
  group_by(ion, strain) %>% 
  mutate(baso_FC = log2(mean(t2)) / log2(mean(t0)))

ggplot(baso_FC, aes(x=baso_FC)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~strain)

## Keep only features with log2FC > x ##
baso_FC = baso_FC %>% 
  dplyr::filter(baso_FC > 1) %>% 
  dplyr::select(-c(pool, t0, t2)) %>% 
  distinct()

length(unique(baso_FC$ion))

# visuzalize results 
ggplot(baso_FC, aes(x=baso_FC)) + 
  geom_histogram(bins = 20) +
  facet_wrap(~strain)

# add annot to baso_FC 
baso_FC = left_join(baso_FC, peak_min[,c(1,7)], by = c("ion" = "feature_ID"))
baso_FC = left_join(baso_FC, MiMeDB[, c(2:3)], by = c("first_annot" = "mime_id"))

length(unique(baso_FC$ion))

# get baso ions info 
baso_ions_info = integrate_baso[integrate_baso$ion %in% unique(baso_FC$ion), ]
write.csv(baso_ions_info, "./output/baso_ions_info.csv")

# how many ions are shared? 
table(baso_FC$strain, baso_FC$ion)
baso_FC %>% 
  group_by(ion) %>% 
  count() %>% 
  filter(n>1)

baso_multiple = baso_FC %>% 
  group_by(ion) %>% 
  count() %>% 
  filter(n>2) %>% 
  left_join(., baso_FC[, c(3,6)]) %>% 
  ungroup %>% 
  group_by(name) %>% 
  slice(1)

library(forcats)  
# Reorder the `name` variable based on the values in `n`
baso_multiple$name <- fct_reorder(baso_multiple$name, baso_multiple$n)


p = ggplot(baso_multiple, aes(x=name, y=n)) + 
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  labs(y="Number of species") + 
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_y_continuous(breaks = c(0, 5, 10)) + 
  cell_plot_theme
print(p)
#ggsave("AB007_baso_shared_4d.tiff", plot = p, dpi = 300, width = 4, height = 4, units = "in")


