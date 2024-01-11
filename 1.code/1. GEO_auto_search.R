# install.packages("rentrez")
# install.packages("XML")
# install.packages("rentrez")
# install.packages("GEOquery")
# BiocManager::install('GEOquery')


library(dplyr)
library(rentrez)
library(XML)
library(GEOquery)
library(httr)
library(jsonlite)
library(stringi)
library(stringr)

# Function to concatenate multiple values with a comma
concatenate_values <- function(value) {
  if (is.null(value)) {
    return(NA)
  } else if (is.character(value) && length(value) > 1) {
    return(paste(value, collapse=","))
  } else {
    return(value)
  }
}

ATC_result <- read.table("ATC_result.txt",sep = "\t",header = T)

all_df<-matrix(ncol = 6,dimnames = list(NA,c('accession','gpl','title','summary','gdstype','drug')))%>%as.data.frame()

no_results <- c()

# Loop through each drug in ATC_result

for (x in 1:nrow(ATC_result)) {
  print(x)
  query_drug <- paste0('"',ATC_result$drug_name[x],'"')
  
  # Construct the query
  query <- paste0(query_drug, '[Title] AND ("Homo sapiens"[porgn] OR "Mus musculus"[porgn] OR "Rattus norvegicus"[porgn]) AND ("gse"[Filter] AND ("Expression profiling by array"[Filter] OR "Expression profiling by high throughput sequencing"[Filter]))')
  
  # Retry logic
  attempt <- 1
  max_attempts <- 5
  successful <- FALSE
  
  while (attempt <= max_attempts && !successful) {
    tryCatch({
      # Perform the search
      search_results <- entrez_search(db = "gds", term = query, retmax = 1)
      total_results <- search_results$count
      search_results <- entrez_search(db = "gds", term = query, retmax = total_results, use_history = TRUE)
      
      # If search results are empty, add to no_results and break the loop
      if (search_results$count == 0) {
        no_results <- c(no_results, query_drug)
        break
      }
      
      # Continue with processing logic...
      # Initialize an empty data frame for results
      results_df <- data.frame(accession = character(), gpl = character(), title = character(), summary = character(), gdstype = character(), stringsAsFactors = FALSE)
      
      # Adjust batch size as needed
      batch_size <- 10
      
      for (i in seq(1, total_results, by = batch_size)) {
        batch_ids <- search_results$ids[i:min(i + batch_size - 1, total_results)]
        # Fetch summaries for the batch of IDs
        if (length(batch_ids) > 0) {
          if (length(batch_ids) == 1) {
            fetched_summaries <- list()
            fetched_summaries[[batch_ids]] <- entrez_summary(db = "gds", id = batch_ids)
          } else {
            fetched_summaries <- entrez_summary(db = "gds", id = batch_ids)
          }
          # Process each summary and add to the data frame
          for (id in names(fetched_summaries)) {
            summary_item <- fetched_summaries[[id]]
            
            GPL <- str_replace_all(summary_item$gpl, pattern = ";", replacement = ",GPL")
            GPL <- paste("GPL", GPL, sep = "")
            results_df <- rbind(results_df, data.frame(
              accession = summary_item$accession,
              gpl = concatenate_values(GPL),
              title = summary_item$title,
              summary = summary_item$summary,
              gdstype = concatenate_values(summary_item$gdstype),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
      results_df <- cbind(results_df, drug = rep(query_drug, nrow(results_df)))
      all_df <- rbind(all_df, results_df)
      successful <- TRUE
    }, error = function(e) {
      cat("Attempt", attempt, "failed:", e$message, "\n")
      attempt <- attempt + 1
      Sys.sleep(60)
    })
  }
  
  # If all attempts fail, log the drug and go to the next iteration
  if (!successful) {
    no_results <- c(no_results, query_drug)
  }
}
all_df_2 <- all_df[!duplicated(all_df$accession),]
all_df <- all_df[!duplicated(all_df$accession),]
openxlsx::write.xlsx(all_df,file = "alldrug_GSE.xlsx")
save(all_df,ATC_result,concatenate_values,fetch_with_retries,file = "all_df_data.Rdata")
