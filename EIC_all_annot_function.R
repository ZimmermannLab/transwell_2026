# function to get targeted EICs for all strains from mzmine preprocessed data
library(xcms)

# function to determine ppm range 
ppm_fun = function(exact_mass){
  ppm = ((100/10^6) * exact_mass)
  mzr = c(exact_mass - ppm, exact_mass + ppm)
  return(mzr)
}


## Defining parameters ##
# centwave algorithm
cwp <- CentWaveParam(
  peakwidth = c(1, 8),
  prefilter = c(1, 5000),
  integrate = 1,
  mzdiff = 0.01, 
  extendLengthMSW = TRUE
)


# merging neighboring peaks 
mpp <- MergeNeighboringPeaksParam(
  expandRt = 4,
  minProp = 0.3
)

EIC_all_mzmine <- function(df_drugs, file_list, output_folder){ 
  
  # Setting the folder to save the plots
  my_output_folder = output_folder
  if (!file.exists(my_output_folder)){
    dir.create(file.path(output_folder))
  }
  
  # create dataframe to store results 
  res_df <- data.frame()
  
  # extract samples for each strain
  samples = file_list
  
  # see which ions belong to each strain 
  df_selection = df_drugs
  
  # Creating the annotated  metadata df
  pd <- data.frame(sample_name = sub(basename(samples), pattern = ".mzML",
                                     replacement = "", fixed = TRUE),
                   stringsAsFactors = FALSE)
  
  raw_data1 <- readMSData(files = samples, 
                          msLevel. = 1, mode = "onDisk", pdata = new("NAnnotatedDataFrame", pd))
  
  for (i in unique(df_selection$ion)){
    
    tryCatch( { 
    mdrug <- df_selection['row.m.z'][df_selection['ion'] == i]
    mzr <- ppm_fun(mdrug)
    rtr <- c(df_selection['row.retention.time'][df_selection['ion'] == i] - 20, df_selection['row.retention.time'][df_selection['ion'] == i] + 20)
    
    # Getting raw chromatograms
    chr_raw <- xcms::chromatogram(raw_data1, rt = rtr, mz = mzr, aggregationFun = "max")
    plot(chr_raw)
    
    # Peak detection on the XIC
    xchr <- findChromPeaks(chr_raw, param = cwp)
    
    # Getting the results and refining peaks
    merged_peaks <- refineChromPeaks(xchr, mpp)
    
    temp <- as.data.frame(chromPeaks(merged_peaks))
    sample_info = as.data.frame(merged_peaks@phenoData@data)
    sample_info = rownames_to_column(sample_info, "column")
    res <- merge(temp, sample_info)
    
    if (nrow(res) == 0) {
      res[1,] = 0
      res$ion <- i
    }
    
    res$ion <- i 
    res_df <- base::rbind(res_df, res)
    }, error = function(e){
      message(paste("Error in ion", i, ":", e$message))})
    
    # plot_fun = function(){ 
    #   png(file = paste0(output_folder, "/", i, "_", "_EIC_int", ".png"))
    #   plot_d = plot(merged_peaks)
    #   title(sub = i)
    #   dev.off()}
    # 
    # ifelse(res[1,] == 0, 
    #        print(paste0("Peak could not be integrated, check manually: ", i)), 
    #        plot_fun())
    
    
    
  } # close value (ion) for loop 
  
  
  return(res_df)
  
} # close function