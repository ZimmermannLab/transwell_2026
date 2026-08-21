# FUNCTION to do wilcoxon test per feature per strain


library(data.table)

compare_tonobug_mzmine = function(complete_df, cols_to_delete){ 
  
  my_df = complete_df
  
  ##Define the class factor ##
  sample_code = unlist(lapply(colnames(my_df)[2:dim(my_df)[2]],
                              function(x) paste0(unlist(strsplit(x, "_"))[2:6], collapse = "_")))
  
  # transpose df 
  my_dataset_t = setDT(my_df)
  my_dataset_t = data.table::transpose(my_df, make.names = 1, keep.names = "rn")
  
  # set abundances columns as numeric 
  my_dataset_t[, 2:dim(my_dataset_t)[2]] <- 
    lapply(my_dataset_t[,2:dim(my_dataset_t)[2]], as.numeric)
  
  
  # log-transform data 
  prelog_data = my_dataset_t[,2:dim(my_dataset_t)[2]]
  # impute NA to 500
  prelog_data[prelog_data == 0] <- 1
  log_data <- log(prelog_data, 2)
  
  
  # include pool & TP metadata 
  sample_name = my_dataset_t$rn
  log_data = cbind(sample_name, log_data)
  
  #plot values per sample
  # df_plot = log_data %>%
  #   pivot_longer(cols = names(log_data)[names(log_data) %like% "FT"],
  #                values_to = "area",
  #                names_to = "ion")
  # 
  # ggplot(df_plot, aes(x=sample_name, y=area)) +
  #   geom_boxplot(alpha= 0.8) +
  #   #geom_point(aes(shape=organ, group=expID), size=2, position=position_jitterdodge()) +
  #   labs(title = paste("all values")) +
  #   xlab("sample") +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  log_data$strain = sapply(strsplit(as.character(log_data$sample_name), "_"), '[', 2)
  log_data$pool = sapply(strsplit(as.character(log_data$sample_name), "_"), '[', 6)
  
  
  
  results = list()
  # run T-test per strain & calculate fold change 
  for(s in unique(log_data$strain)){ 
    
    # select data per drug 
    df1 = log_data %>%
      ungroup() %>% 
      filter(strain == s) %>% 
      dplyr::select(-any_of(cols_to_delete)) %>% 
      as.matrix() 
    
    df2 = log_data %>%  
      ungroup() %>% 
      filter(strain == "nobug") %>% 
      dplyr::select(-any_of(cols_to_delete)) %>% 
      as.matrix()
    
    # Create an empty matrix
    pvalue = matrix(NA, ncol(df1), ncol=4)
    
    for(i in 1:length(colnames(df1))){
      FC = mean(df1[,i], na.rm = T) - mean(df2[,i], na.rm=T) # FC is calculated by subtraction, because values are already log2-transformed
      test_result = t.test(df1[,i], df2[,i], alternative = "two.sided", paired=F) 
      pvalue[i, 1] <- test_result$p.value
      pvalue[i, 2] <- test_result$statistic
      pvalue[i, 4] <- FC
      row.names(pvalue) <- colnames(df1)
    }
    
    # do p-value adjustment
    pvalue[, 3] <- stats::p.adjust(pvalue[, 1], method="BH",
                                   n=length(na.omit(pvalue[, 1])))
    
    colnames(pvalue) = c("pval", "statistic", "padj", "FC")
    pvalue = as.data.frame(pvalue)
    
    pvalue = rownames_to_column(pvalue, var="ion")
    
    results[[s]] = pvalue
    
    
  }
  
  return(results)
  
} 


