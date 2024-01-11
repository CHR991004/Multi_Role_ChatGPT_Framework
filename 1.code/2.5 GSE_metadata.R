# Using the R to batch read metadata into a list

# Load required library
library(openxlsx)
library(GEOquery)

# Read the Excel file containing GEO IDs
results_df <- openxlsx::read.xlsx("alldrug_GSE.xlsx")

# Extract GSE IDs from the dataframe
gse_ids <- results_df$accession
unique_gse_ids <- unique(gse_ids)

# Initialize variables
gse_list <- list()
error_char <- c()

# Directory for storing temporary files
if (!dir.exists("temp")) {
  dir.create("temp")
}

# Iterate over each unique GSE ID
for (i in 1:length(unique_gse_ids)) {
  print(i)
  gse_id <- unique_gse_ids[i]
  
  # Initialize retry parameters
  max_retries <- 5
  download_successful <- FALSE
  
  # Try to download data, with retries
  for (attempt in 1:max_retries) {
    tryCatch({
      gse_data <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE, destdir = "./temp")
      
      # Check if matrix format data is returned
      if (length(gse_data) > 0) {
        gse_list[[gse_id]] <- gse_data[[1]]
        download_successful <- TRUE
        break # Break the loop if download is successful
      }
    }, error = function(e) {
      message("Attempt ", attempt, " failed for ", gse_id, ": ", e$message)
    })
    
    if (download_successful) {
      break # Continue to next GSE ID if current download is successful
    }
  }
  
  # Record failed GSE ID after maximum retries
  if (!download_successful) {
    message("No Series Matrix File found after multiple attempts for ", gse_id)
    error_char <- c(error_char, gse_id)
  }
}

# Filter the original dataframe for successfully downloaded GSE IDs
results_df_ingselist <- results_df[results_df$accession %in% names(gse_list),]

# Save the filtered dataframe and the list of GSE data
# save(results_df_ingselist, gse_list, file = "ATCdrug_metadata.Rdata")

# If you need to write the error GSE IDs to a supplementary file
write.table(error_char, file = "GSE_id_error.xls", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
