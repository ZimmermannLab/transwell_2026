## AB005 drug screen analysis ##

# Set-up 
library("data.table")
library("dplyr")
library("ggplot2")
library("tidyverse")
library("RColorBrewer")
library("stringr")
library("hrbrthemes")
library( "gplots" )


#######Import data###########
load(file.path("input_folder", "Fig2A_drugs_all_raw_with_metadata.rds"))

#### make sure to only select one peak per drug ####
temp = my_df %>% 
  group_by(sample_name, drug) %>% 
  count(well_location, sort=T) %>% 
  filter(n >1)
length(unique(temp$drug))
length(unique(temp$sample_name))
hist(temp$n)

### Select one peak per drug based on closest RT ###
# read-in pooling scheme with RT for stds 
# list with mass and rt for the molecules
pooling_scheme <- read.csv(file.path("input_folder", "NEW_pooling_scheme_info_with_IS.csv"))
pooling_scheme <- na.omit(pooling_scheme) # new version of doc has NAs

# merge pooling scheme with my_df 
tmp = left_join(my_df, pooling_scheme[, c(1, 30:35)], by = c("drug" = "drug"))
tmp$RT.NR = tmp$RT.NR * 60

# calculate delta RT between peak detected and RT.NR 
# pick only peak with smallest delta RT when drug has >1 peak in sample
test = tmp %>% 
  mutate(deltaRT = rt - RT.NR) %>% 
  mutate(deltaRT = abs(deltaRT)) %>% 
  group_by(sample_name, drug) %>% 
  filter(deltaRT == min(deltaRT))

length(unique(test$drug))
length(unique(tmp$drug))
table(my_df$strain)
table(test$strain)

# plot raw data 
pdf("./Plot/rawdata_single_drugs_scalesfree.pdf")
for (d in unique(test$Product.Name)){ 
  
  temp_df = test %>% 
    filter(Product.Name == d)
  
  p <-  ggplot(temp_df, aes(x=TPs, y=into, color = pool, 
                            group=sample_name)) +
    geom_point() +
    geom_line(aes(group=pool)) +
    labs(title= paste0("Timecourse ", d), 
         x = "Time (hours)", 
         y = "Raw peak area") + 
    theme_minimal() + 
    facet_wrap(~strain, scales = "free") + 
    geom_blank()
  
  print(p)
}

dev.off()

#### Check integration: how many drugs have 4 reps at TP0? + distribution of areas ####
integration_check = test %>% group_by(Product.Name, strain, TPs) %>% 
  summarise(counts = n())

integration_check %>% filter(!grepl("IS_", Product.Name)) %>%
  filter(TPs == "TP0") %>% 
  ggplot(., aes(x=counts)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~strain)

# are always the same drugs not integrated? 
integration_check %>% filter(!grepl("IS_", Product.Name)) %>%
  filter(TPs == "TP0") %>% 
  filter(counts < 4) %>% 
  ggplot(., aes(x=strain, y=Product.Name)) + 
  geom_tile(aes(fill=counts)) + 
  scale_fill_continuous(type = "gradient") + 
  theme(axis.text.y = element_text(size = 5)) # y-axis text size

# how does this relate to total dataset? 
integration_check %>% filter(!grepl("IS_", Product.Name)) %>%
  filter(TPs == "TP0") %>% 
  ggplot(., aes(x=strain, y=Product.Name)) + 
  geom_tile(aes(fill=counts))


# overlap between drugs per strain 
integration_check %>% 
  filter(TPs == "TP0") %>% 
  mutate(integration = ifelse(counts >= 4, "4reps", "<4reps")) %>% 
  ggplot(., aes(x=integration)) + 
  geom_bar(stat = "count") + 
  facet_wrap(~strain)

# how many drugs are not detected at TP0? (at least 4 reps)
integration_check %>% 
  filter(TPs == "TP0") %>% 
  filter(counts < 4) %>% 
  ungroup() %>% 
  select(Product.Name) %>% 
  summarise(n_distinct(Product.Name))
# 241 drugs are included 
  
# for how many are we mising "only" 1 rep? 
integration_check %>% 
  filter(TPs == "TP0") %>% 
  filter(counts < 3) %>% 
  ungroup() %>% 
  select(Product.Name) %>% 
  summarise(n_distinct(Product.Name))
# 172


# what if I kick out Cscindens & Dformi? 
integration_check %>% 
  filter(TPs == "TP0") %>% 
  filter(!(strain == "Cscindens" | strain == "Dformi")) %>% 
  filter(counts < 4) %>% 
  ungroup() %>% 
  select(Product.Name) %>% 
  summarise(n_distinct(Product.Name))

integration_check %>% 
  filter(TPs == "TP0") %>% 
  filter(!(strain == "Cscindens" | strain == "Dformi")) %>%
  filter(counts < 4) %>% 
  ggplot(., aes(x=strain, y=Product.Name)) + 
  geom_tile(aes(fill=counts)) + 
  scale_fill_continuous(type = "gradient") + 
  theme(axis.text.y = element_text(size = 5)) # y-axis text size

# check Ef signals 
integration_check %>% 
  filter(TPs == "TP0") %>% 
  filter(!grepl("IS_", Product.Name)) %>%
  filter((strain == "Efaecalis" | strain == "Ef2")) %>%
  ggplot(., aes(x=strain, y=Product.Name)) + 
  geom_tile(aes(fill=counts)) + 
  scale_fill_continuous(type = "gradient") + 
  theme(axis.text.y = element_text(size = 5)) # y-axis text size


######Normalize by IS###########################

my_prep_df = test
# Divide samples into Ctrl VS samples VS internal standards
setDT(my_prep_df)
my_prep_df = my_prep_df[grep("AB005", my_prep_df$sample_name, fixed = T), Status := "SAMPLE"]
my_prep_df = my_prep_df[grep("IS_",drug,fixed=T), Status := "INT_STD"] 


save(my_prep_df,file = "Env/pre_normalization.rds")


# extract Drugs labeled as "Internal Standard"
IS_drugs<-grep("IS_", unique(my_prep_df$drug), fixed = T,value = T)
# select the IS you want
IS_only_df = my_prep_df[drug%in%IS_drugs,]

# visualize IS Area
library(tidyverse)
library(ggplot2)
setDF(IS_only_df)
IS_only_df %>% ggplot(aes(x=strain, y=intb, fill=TPs)) + 
  geom_boxplot() + facet_wrap(~drug, scale="free") + theme(axis.text.x = element_text(angle = 90)) + 
  theme(axis.text.x = element_text(size = 6))
#ggsave(filename = 'Plot/IS_Area_by_strain_in_time.png', dpi=300, width = 25, height = 20)

# kick out sulfametoxazole as IS 
IS_only_df = IS_only_df %>% 
  filter(drug != "IS_sulfamethoxazole")


# read function 
# Historical normalization-function source calls were removed; the required
# normalization code is included below so the analysis does not depend on lab paths.

# median values over all IS samples per experiment  
median_IS  = IS_only_df %>%
  group_by(interaction(strain, drug)) %>%
  mutate(across(
    .cols = intb,
    .fns = list(Median =  median, FC = ~.x / median(.x, na.rm = T)), na.rm = TRUE, 
    .names = "{col}_{fn}"))


# Let's have a look at a boxplot of logged fold changes
boxplot(median_IS$intb_FC[median_IS$drug == "IS_caffeine"])
boxplot(median_IS$intb_FC[median_IS$drug == "IS_tolfenamic_acid"])
boxplot(median_IS$intb_FC[median_IS$drug == "IS_warfarin"])
boxplot(median_IS$intb_FC[median_IS$drug == "IS_ipriflavone"])



#Define samples in which FC is >2 or <1/2 of median
# toBeRemoved <- median_IS %>% 
#   filter(
#     if_any(ends_with("FC"), ~ . > 2 | . < 0.5)
#   )

#Remove those time points from further analysis
#my_df <- anti_join(median_IS, toBeRemoved)
#median_IS <- anti_join(median_IS, toBeRemoved)

# Let's have a look at the IS Area again
median_IS %>% ggplot(aes(x=strain, y=intb_FC, fill=TPs)) + 
  geom_boxplot() + facet_wrap(~drug, scale="free") + theme(axis.text.x = element_text(angle = 90)) + 
  theme(axis.text.x = element_text(size = 6))


# Compute correction factor to normalize measurements
# calculate the mean of the median_FC per well  
median_IS <- median_IS %>%
  group_by(interaction(strain, pool)) %>%
  mutate(across(
    .cols = intb_FC,
    .fns = list(FC_ISmedian = ~mean(.x, na.rm = T)), 
    .names = "{col}_{fn}")) %>% 
  ungroup()

# Let's have a look at the IS Area again
median_IS %>% ggplot(aes(x=strain, y=intb_FC_FC_ISmedian, fill=TPs)) + 
  geom_boxplot() + facet_wrap(~drug, scale="free") + theme(axis.text.x = element_text(angle = 90)) + 
  theme(axis.text.x = element_text(size = 6))
  
# make minimal median_IS
median_IS <- median_IS %>% 
    select(sample_name, intb_FC_FC_ISmedian) %>% 
  distinct()

# Summarize FCs from different IS compounds per sample
# median_IS <- median_IS %>% 
#   group_by(sample_name) %>%
#   summarise(ISMean_per_well = mean(intb_FC_FC_ISmean, na.rm=T))


#For each strain, pool and timepoint, divide the raw area of the drug (in that pool) by the FC to generate the corrected area value
my_norm_df <- left_join(my_prep_df, median_IS, by = "sample_name")
my_norm_df$correctedArea <- my_norm_df$intb / my_norm_df$intb_FC_FC_ISmedian


# Check corrected area for IS compounds 
# select the IS you want
IS_only_df = my_norm_df[drug %in% IS_drugs,]

# kick out Ipriflavone as IS 
IS_only_df = IS_only_df %>% 
  filter(drug != "IS_sulfamethoxazole")

# visualize IS Area
setDF(IS_only_df)
IS_only_df %>% ggplot(aes(x=strain, y=correctedArea, fill=TPs)) + 
  geom_boxplot() + facet_wrap(~drug, scale="free") + theme(axis.text.x = element_text(angle = 90)) + 
  theme(axis.text.x = element_text(size = 6))
ggsave(filename = 'Plot/IS_Area_by_strain_in_time_normalized.png', dpi=300, width = 25, height = 20)


# save the file
save(my_norm_df, file = "Env/post_normalization.rds")
save(my_norm_df, file = "Env/post_normalization.csv")
write.csv(my_norm_df, "AB005_tar_data_norm_240308.csv")

length(unique(my_norm_df$drug))

# kick out bad runs
#my_norm_df = my_norm_df %>% 
 # filter(!(strain == "Dformi" | strain == "Cscindens" | strain == "Ef2"))

##### Plot drug/time plots ####
pdf("./Plot/AB005_ISnorm_singlestrains_240321.pdf")
for (d in unique(my_norm_df$Product.Name)){ 

  temp_df = my_norm_df %>% 
    filter(Product.Name == d)
  
  p <-  ggplot(temp_df, aes(x=TPs, y=correctedArea, color = pool, 
                          group=interaction(pool, sample_name))) +
    geom_point() +
    geom_line(aes(group=interaction(pool))) +
    labs(title= paste0("Timecourse ", d), 
        x = "Time (hours)", 
        y = "IS_normalized peak area") + 
    theme_minimal() + 
    facet_wrap(~strain) +
    geom_blank()

print(p)

} 
dev.off()


#### Check data: assess noise in bwteen pools ####
# plot boxplot of all pools per TP per strain 


pdf("./Plot/pool_intensity_boxplot.pdf")
for (i in levels(as.factor(my_norm_df$strain))){ 
  temp_df = my_norm_df %>% 
    filter(strain == i)
  
  p = ggplot(temp_df, aes(x=well_location, y=correctedArea))+ 
    geom_boxplot() + 
    labs(title= paste0("Boxplot of intensities per pool ", i), 
         x = "Pool location", 
         y = "IS_normalized peak area") + 
    theme_minimal() + 
    facet_wrap(~TPs) + 
    geom_blank()
  
  print(p)
}
dev.off()



##### Identify degraded drugs for TP8 ##### 
# Optional checkpoint: load("Env/post_normalization.rds")

#selecting the desired column for further processing
my_cols=c("sample_name","Product.Name","rt","intb","pool","strain","TPs","well_location", "correctedArea")
df = my_norm_df[,..my_cols]

# kick out bad run
df = df %>% 
 filter(!(strain == "Efaecalis"))

unique(df$strain)

# keep only drugs with 3 reps
df = na.omit(df)
setDT(df)
to_keep=df[,.N, by=.(strain, Product.Name, TPs)][N>=3,]

# how many drugs per strain? 
to_keep %>% 
  group_by(strain) %>%
  summarise(count = n_distinct(Product.Name))

# Semi-join 
# (return only rows of big that have a 
# match in small)
setDT(df)
df=df[na.omit(df[to_keep, on =.(strain, Product.Name, TPs), which=TRUE])]


setDF(df)
# Check: plot random subsets of drug for each strain to check data quality 
##### Plot drug/time plots ####
pdf("./Plot/AB005_timecourses_selected_perstrain.pdf")
for (s in unique(df$strain)){ 
  
  temp_df = df %>% 
    filter(strain == s)
  
  length(unique(temp_df$Product.Name))
  
  # select random subset of drugs 
  drugs = sample(temp_df$Product.Name, 10)
  tp = temp_df[temp_df$Product.Name %in% drugs, ]
  
  p <-  ggplot(tp, aes(x=TPs, y=correctedArea, color = pool, 
                            group=interaction(pool, sample_name))) +
    geom_point() +
    geom_line(aes(group=interaction(pool))) +
    labs(title= paste0("Timecourses random subset: ", s), 
         x = "Time (hours)", 
         y = "IS_normalized peak area") + 
    facet_wrap(~Product.Name, scales = "free") +
    theme_minimal() + 
    geom_blank()
  
  print(p)
  
} 

dev.off()

#### Calcuate FC for TP8 ####
df = df %>% 
  filter(!grepl("IS_", Product.Name))

setDT(df)

#spliting the df based on time
df1 =df%>%filter(TPs == "TP0")
df2 =df%>%filter(TPs == "TP8")


df1 = df1 %>% select(sample_name,Product.Name, pool, strain, correctedArea)
df2 = df2 %>% select(sample_name,Product.Name, pool, strain, correctedArea)

#renaming the column
names(df1)[names(df1) == "correctedArea"] <- "0hr"
names(df2)[names(df2) == "correctedArea"] <- "8hr"

#Joining two df
library(tidyverse)
df3 = full_join(df1, df2, by = c("Product.Name","pool","strain"))
# exclude combos where I'm missing one of the two timepoints
df3 = df3[!is.na(sample_name.x) & !is.na(sample_name.y),]

#computing the FC based on t=0hr
df_FC = df3 %>% group_by(Product.Name, strain) %>% mutate(FC=log2(mean(`8hr`, na.rm = T)/mean(`0hr`, na.rm = T)))

# compute the pvalue associated to FC on t=0hr
df_FC = df_FC %>% group_by(Product.Name, strain) %>% mutate(FC_pval = tryCatch({
  t.test(`8hr`,`0hr`, alternative = "less")$p.value
},
error=function(cond) {
  return(1)
})
)

# adjust p-values with fdr
df_temp = df_FC %>% select(Product.Name, strain, FC, FC_pval) %>% distinct()

df_temp = df_temp %>% 
  group_by(Product.Name) %>% 
  mutate(padj_strain = p.adjust(FC_pval, method = "BH"))

df_FC_corrected = left_join(df_FC, df_temp, by = c("Product.Name","strain","FC","FC_pval"))
hist(df_FC_corrected$padj_strain)

#pval_thresh = 0.05

# let's take a look at the FC with at least negative drug degradation
result1  = df_FC_corrected %>% filter(FC < 0 & padj_strain < 0.2) %>% 
  select(Product.Name, strain, FC, FC_pval, padj_strain) %>%
  distinct()

result1 %>% 
  group_by(strain) %>% 
  count(strain)

# add TP name 
result1$TP = "TP8"



##### Identify degraded drugs for TP2 ##### 
#spliting the df based on time
df1 =df%>%filter(TPs == "TP0")
df2 =df%>%filter(TPs == "TP2")

df1 = df1 %>% select(sample_name,Product.Name, pool, strain, correctedArea)
df2 = df2 %>% select(sample_name,Product.Name, pool, strain, correctedArea)

#renaming the column
names(df1)[names(df1) == "correctedArea"] <- "0hr"
names(df2)[names(df2) == "correctedArea"] <- "2hr"

#Joining two df
library(tidyverse)
df3 = full_join(df1, df2, by = c("Product.Name","pool","strain"))
# exclude combos where I'm missing one of the two timepoints
df3 = df3[!is.na(sample_name.x) & !is.na(sample_name.y),]


#computing the FC based on t=0hr
df_FC = df3 %>% group_by(Product.Name, strain) %>% mutate(FC=log2(mean(`2hr`, na.rm = T)/mean(`0hr`, na.rm = T)))

# compute the pvalue associated to FC on t=0hr
df_FC = df_FC %>% group_by(Product.Name, strain) %>% mutate(FC_pval = tryCatch({
  t.test(`2hr`,`0hr`, alternative = "less")$p.value
},
error=function(cond) {
  return(1)
})
)

# adjust p-values with fdr
df_temp = df_FC %>% select(Product.Name, strain, FC, FC_pval) %>% distinct()

df_temp = df_temp %>% 
  group_by(Product.Name) %>% 
  mutate(padj_strain = p.adjust(FC_pval, method = "BH"))

df_FC_corrected_TP2 = left_join(df_FC, df_temp, by = c("Product.Name","strain","FC","FC_pval"))


# let's take a look at the FC with least x% drug degradation
result2  = df_FC_corrected_TP2 %>% filter(FC < 0 & padj_strain < 0.2) %>% 
  select(Product.Name, strain, FC, FC_pval, padj_strain) %>%
  distinct()

# add TP info  
result2$TP = "TP2"

# how many drugs identified as degraded?
result2 %>% 
  group_by(strain) %>% 
  count(strain)

# rbind result1 (TP8) and result2 (TP2)
results = rbind(result1, result2)


## filter to remove Efaecalis ## 
results = results %>% 
  filter(strain != "Efaecalis")

# how many drugs are degraded? 
length(unique(results$Product.Name))

#### plot barplot per strain & per drug ####

# merge results for TP2 and TP8 
df_FC_corrected$TP = "TP8"
df_FC_corrected_TP2$TP = "TP2"
df_FCs = rbind(df_FC_corrected, df_FC_corrected_TP2)

# per drug 
p_list = list()

for(d in unique(df_FCs$Product.Name)){
  
  tmp = df_FCs %>% 
    filter(Product.Name == d)
  
  tmp <- tmp %>% mutate(pval = ifelse(padj_strain<=0.2,"signif","non_signif"))
  
  
  p = ggplot(tmp, aes(x=strain, y=FC, fill = pval)) +
    geom_bar(stat='identity') +
    # scale_y_log10(limits = c(10^-2, 10^1)) +
    coord_flip() +
    scale_fill_manual(name = "p-value", labels = c(">0.1", "<=0.1"),
                      values = c("gray70",
                                 "steelblue1","tomato1",
                                 "gray60",
                                 "steelblue4","tomato4"))+#"#BAB0AC", "#E15759")) +
    ggtitle(d) + 
    facet_wrap(~TP)
  
  #plot(p)
  
  
  p_list[[d]] = p 
  
}


ggsave("./Plot/AB005_all_TPs.pdf", 
       gridExtra::marrangeGrob(grobs = p_list, nrow = 2, ncol = 2), 
       height = 6, width = 10)


#### plot normalized data for all normal-time hits #### 
## plot hits per strain ## 
p_list = list()
for (s in unique(results$strain)){ 
  
  temp_df = df %>% 
    filter(strain == s)
  
  # select drugs that are hits 
  drugs = results[results$strain == s, ]
  drugs = drugs$Product.Name
  tp = temp_df[temp_df$Product.Name %in% drugs, ]
  
  
  for (d in drugs){
    
    tmp = tp %>% 
      filter(Product.Name == d)
    
    p <-  ggplot(tmp, aes(x=TPs, y=correctedArea, color = pool, 
                          group=interaction(pool, sample_name))) +
      geom_point() +
      geom_line(aes(group=interaction(pool))) +
      labs(title= paste0("Timecourses: ", s, "_", d), 
           x = "Time (hours)", 
           y = "IS_normalized peak area") + 
      theme_minimal() + 
      geom_blank()
    
    #print(p)
    
    p_list[[paste0(s, "_", d)]] = p
    
  }
  
} 

print(p_list[[1]])

ggsave("./Plot/AB005_hits.pdf", 
       gridExtra::marrangeGrob(grobs = p_list, nrow = 2, ncol = 2), 
       height = 6, width = 8)


#### check for fast degraded drugs ####

# we will identify them by comparing t0 with drug to t0 in controls
# compute the pvalue associated to FC on t=0hr 
df_FC_quick = c()

for(d in unique(df_FC_corrected$Product.Name)){
  
  temp = df_FC_corrected %>% filter(Product.Name == d & strain != "nobug") %>% select(-`8hr`)
  ctrl = df_FC_corrected %>% filter(Product.Name == d & strain == "nobug") %>% select(`0hr`)
  
  temp = temp %>% group_by(strain) %>% mutate(FC_quick = log2(mean(`0hr`, na.rm = T)/mean(ctrl$`0hr`, na.rm = T)),
                                                        FC_quick_pval = tryCatch({
                                                          t.test(`0hr`, ctrl$`0hr`, alternative = "two.sided")$p.value
                                                        },
                                                        error=function(cond) {
                                                          return(1)
                                                        })
  )
  
  df_FC_quick = rbind(df_FC_quick, temp)
  
}


df_FC_quick$FC_quick_padj = p.adjust(df_FC_quick$FC_quick_pval, method = "BH")

# let's consider a drug-bug combo a quick metabolizer if it respects several conditions:
result3  = df_FC_quick %>% filter(FC_quick <= -2 & FC_quick_padj < 0.05) %>% 
  select(Product.Name, strain, FC_quick, FC_quick_pval, FC_quick_padj) %>%
  distinct()

# add TP info  
result3$TP = "quick"

# how many drugs identified as degraded?
result3 %>% 
  group_by(strain) %>% 
  count(strain)

# make colnames equal 
colnames(results)
colnames(result3) <- c("Product.Name", "strain", "FC", "FC_pval", "padj_strain", "TP")

results = rbind(results, result3)
length(unique(results$Product.Name))

#### plot timecourses of results ####
# disregard quick hits 
results = rbind(result1, result2)

# merge df with results to only contain the "hits"
final_selection = left_join(results, df, by = c("strain", "Product.Name"))

# make a p-val selection 
hist(final_selection$padj_strain)

pval_thresh = 0.2

final_selection = final_selection %>% 
 filter(padj_strain <= pval_thresh)

final_selection %>% 
  select(Product.Name, strain, TP) %>% 
  distinct() %>% 
  group_by(strain, TP) %>% 
  count(strain) %>% 
  print(n=30)

pdf("./Plot/AB005_timecourses_hitsonly.pdf")
for (d in unique(final_selection$Product.Name)){ 
  
  temp_df = final_selection %>% 
    filter(Product.Name == d) 
  
  p <-  ggplot(temp_df, aes(x=TPs, y=correctedArea, color = pool, 
                            group=interaction(pool, sample_name))) +
    geom_point() +
    geom_line(aes(group=interaction(pool))) +
    labs(title= paste0("Timecourse ", d),
         subtitle = paste("Ajd. p-vals:", toString(unique(temp_df$padj_strain))),
         x = "Time (hours)", 
         y = "IS_normalized peak area") + 
    theme_minimal() + 
    facet_wrap(~strain + TP, scales = "free") +
    geom_blank()
  
  print(p)
  
} 

dev.off()

save.image("AB005_tar_results.RData")

write.csv(final_selection, "AB005_degraded_drugs_areas.csv")

# final numbers 
final_selection %>% 
  select(Product.Name, strain, TP) %>% 
  distinct() %>% 
  group_by(strain, TP) %>% 
  count(strain) %>% 
  print(n=30)

####


#### Transform FC to percentage #### 
# calculate percentage 
#results = rbind(result1, result2)

#results = rbind(result1, result2)
#results = rbind(results, result3)

results = results %>% 
  mutate(percentage = (100-((2^(FC))*100))) 

saved_results = results
range(results$percentage)

#### P-val distributions ####
hist(results$FC_pval, breaks = 30)
hist(results$padj_strain, breaks = 30)


ggplot(results, aes(x=padj_strain)) + 
  geom_histogram() + 
  facet_wrap(~strain)

length(unique(final_selection$Product.Name))
length(unique(results$Product.Name))

#### Filter results  ####
# filter for at least 10% degradation 
hist(results$percentage)
results = results %>% 
  filter(percentage > 10)

results %>% select(Product.Name, strain, TP) %>% 
  distinct() %>% 
  group_by(strain, TP) %>% 
  count(strain) %>% 
  print(n=30)

# filter for spontaneous degradation per time point 
results_nobug = results %>% 
  filter(strain == "nobug") %>% 
  dplyr::select(Product.Name, control = percentage, TP)

# Add control degradation to the original dataframe
fil_results <- results %>%
  left_join(results_nobug, by = c("Product.Name", "TP"))

# Annotate rows that pass the filter (degradation > nobug degradation)
fil_results <- fil_results %>%
  mutate(PassesFilter = ifelse(
    strain.x != "nobug" & (!is.na(control) & percentage > control),
    "yes",
    "no"
  ))

# now delete the rows that have strain.y = "nobug" & PassesFilter == F
to_delete = fil_results %>% 
  dplyr::filter((strain.y == "nobug" & PassesFilter == "no"))

fil_results = anti_join(fil_results, to_delete)

fil_results %>% 
  select(Product.Name, strain.x, TP) %>% 
  distinct() %>% 
  group_by(strain.x, TP) %>% 
  count(strain.x) %>% 
  print(n=30)
length(unique(fil_results$Product.Name))

# filter single hits for Cscindens & Dformi 
# results_fil1 = results %>% 
#   group_by(Product.Name) %>% 
#   add_count(Product.Name) %>% 
#   filter((strain == "Cscindens" | strain == "Dformi") & TP == "quick" & n == 1) %>% 
#   select(strain, Product.Name)
# 
# results_fil2 = results %>% 
#   filter(strain == "Cscindens" | strain == "Dformi") %>% 
#   group_by(Product.Name) %>% 
#   add_count(Product.Name) %>% 
#   filter(TP == "quick" & n == 2) %>% 
#   select(-n)
#   
# # anti join the results 
# results_final = anti_join(results, results_fil1, 
#                           by = c("Product.Name", "strain"))
# results_final = anti_join(results_final, results_fil2)

hist(fil_results$padj_strain)

#### plot heatmap 2024-03-21 ####
library(ComplexHeatmap)

results_final = fil_results %>% 
  dplyr::select(c(Product.Name, strain = strain.x, FC, FC_pval, padj_strain, TP, percentage))
write.csv(results_final, "./AB005_results_p02_10percent_241126.csv")
results_mat = results_final %>% 
  select(Product.Name, strain, percentage, TP) %>% 
  pivot_wider(names_from = strain, 
              values_from = percentage, 
              values_fill = 0)

# Calculate the number of non-zero rows for each column
non_zero_counts <- colSums(results_mat != 0)
# Reorder columns based on the non-zero counts in descending order
results_mat <- results_mat[, order(non_zero_counts, decreasing = TRUE)]

# Calculate row sums (or use another metric, like row maximums)
row_sums <- rowSums(results_mat[, 3:ncol(results_mat)])
# Reorder rows based on the row sums
results_mat <- results_mat[order(row_sums, decreasing = TRUE), ]


length(unique(results_final$Product.Name))
drugs = results_mat$Product.Name
TPs = results_mat$TP

results_mat = as.matrix(results_mat[3:ncol(results_mat)])
#rownames(results_mat) = drugs

colnames(results_mat)
colnames(results_mat) <- c("C. ramosum", "C. comes", "B. uniformis", "G. haemolysans",
                           "B. thetaiotaomicron", "E. faecalis", "A. naeslundi", 
                           "A. omnicolens", "D. formicigenerans",  "C. scindens")




png(file = './Plot/results_Heatmap_allTP_pval02_241126.png',
    width = 22, height = 14, units = "cm", res=800)
Heatmap(results_mat, name = "% of degradation", 
        col = colorRampPalette(c("black", "white"))(100), 
        right_annotation = rowAnnotation(TP = TPs, 
                                         col = list(TP = c("TP8" = "#C2185B", "TP2" = "#F06292"))), 
        row_title = "Drugs tested", 
        column_title = "Bacterial strains", 
        row_names_gp = gpar(fontsize = 4), 
        heatmap_legend_param = 
          list(title_gp = gpar(fontsize = 12, lineheight = 1.5)), 
        column_title_gp = 
          gpar(fontsize = 12, lineheight = 1.2), 
        cluster_columns = F
        
)
dev.off()

#### heatmap for figure 2A ####
# Heatmap creation with italic column labels
ht <- Heatmap(
  results_mat, 
  name = "% of degradation", 
  col = colorRampPalette(c("black", "white"))(100), 
  right_annotation = rowAnnotation(
    TP = TPs, 
    col = list(TP = c(
      "TP8" = "#C2185B", 
      "TP2" = "#F06292"
    ))
  ), 
  row_title = "Drugs tested", 
  column_title = "Bacterial strains", 
  row_names_gp = gpar(fontsize = 4), 
  column_names_gp = gpar(fontsize = 8, fontface = "italic"), # Italic column labels
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

tiff("./Plot/AB005_heatmap_fig2.tiff", width = 6.5, height = 3, units = "in", res = 300)
draw(ht)
dev.off()

#### Compare hits to mapping paper #### 
mapping_results = read.csv(file.path("input_folder", "Mapping_BugDrug_all_binary.csv"))
mapping_results = mapping_results %>% 
  mutate(deg_mapping = rowSums(select(., where(is.numeric)))) %>% 
  dplyr::select(DrugName, deg_mapping) %>% 
  mutate(deg_mapping = ifelse(deg_mapping > 1, 1, 0))



# make similar binary table for my screen 
screen_binary = pooling_scheme %>% 
  dplyr::select(Product.Name, M.SD)

screen_binary$deg_screen = ifelse(screen_binary$Product.Name %in% unique(results_final$Product.Name), 
                                  1, 0)
sum(screen_binary$deg_screen)

# fuzzy join based on the drug names from both screens 
library(fuzzyjoin)
mapping_join = fuzzyjoin::stringdist_left_join(screen_binary, mapping_results, 
                                          by=c("Product.Name" = "DrugName"), 
                                          ignore_case = TRUE, 
                                          method = "jw", 
                                          max_dist = 99, 
                                          distance_col = "dist") %>%  
  group_by(Product.Name) %>% 
  slice_min(order_by = dist, n=1)

hist(mapping_join$dist)
mapping_join = mapping_join %>% 
  filter(dist < 0.1)

mapping_join$overlap = ifelse(mapping_join$deg_screen == mapping_join$deg_mapping, 
                              T, F)

#### plot single examples #### 
d= "Azilsartan"

temp_df = df %>% 
  filter(Product.Name == d) %>% 
  mutate(time = as.numeric(sub("TP", "", TPs)))

temp_df %>% 
  filter(strain %in% c("Btheta", "Buni", "Ccomes", "Cramosum", "nobug")) %>% 
  ggplot(., aes(x= time, y= correctedArea)) + 
  geom_point() + 
  stat_summary(fun.data ='mean_se', geom = "smooth", se = TRUE, alpha=0.1) + 
  theme_minimal() + 
  facet_wrap(~strain) +
  labs(x="Incubation time", 
       y="Peak AUC") +
  ggtitle(paste0("Time course "), d)


d= "Roxatidine Acetate HCl"
temp_df = df %>% 
  filter(Product.Name == d) %>% 
  mutate(time = as.numeric(sub("TP", "", TPs)))

temp_df %>% 
  filter(strain %in% c("Buni", "nobug")) %>% 
  ggplot(., aes(x= time, y= correctedArea)) + 
  geom_point() + 
  stat_summary(fun.data ='mean_se', geom = "smooth", se = TRUE, alpha=0.1) + 
  theme_minimal() + 
  facet_wrap(~strain) +
  labs(x="Incubation time", 
       y="Peak AUC") +
  ggtitle("Roxatidine Acetate & bacterial incubation")

d = "Carfilzomib (PR-171)"
temp_df = df %>% 
  filter(Product.Name == d) %>% 
  mutate(time = as.numeric(sub("TP", "", TPs)))

temp_df %>% 
  filter(strain %in% c("Btheta", "Ccomes", "Aomnicolens", "nobug")) %>% 
  ggplot(., aes(x= time, y= correctedArea)) + 
  geom_point() + 
  stat_summary(fun.data ='mean_se', geom = "smooth", se = TRUE, alpha=0.1) + 
  theme_minimal() + 
  facet_wrap(~strain) +
  labs(x="Incubation time", 
       y="Peak AUC") + 
  ggtitle(paste0("Time course "), d)
















