# Note: After converting to URLs,
# They can be used with IDM or other download software to achieve multi-threaded batch downloading,
# thus realizing higher efficiency.

# Load required libraries
require(doParallel)
library(stringr)
library(GEOquery)
library(xml2)
library(parallel)
library(openxlsx)
source("function.R") # Source additional functions from a file named "function.R"

options(timeout=60) # Set the timeout to 60 seconds

# Detect the number of cores on the machine and use two less than the maximum
n.cores <- detectCores() - 2

# Set the input file name
input_file <- "GSE_id_supple.xls"

# Read the input file. Adjust 'header' based on whether the file has a header row
all_GSE <- read.table(input_file, sep = "\t", header = FALSE)

# Extract unique GEO IDs using a regular expression
GEO <- unique(unlist(str_match_all(all_GSE[,1], pattern = "GSE[0-9]*")))

# Initialize an empty vector to store results
merge <- c()

# Create a cluster for parallel processing
clust <- makeCluster(n.cores)

# Parallel processing: Apply a function over GEO IDs
a <- parLapply(clust, GEO, fun = url, merge, getDirListing)

# Stop the cluster after processing
stopCluster(clust)

# Unlist and combine the results
a_char <- unlist(a)

# Write the results to an Excel file
write.table(a_char, file = "GSE_id_url.xls", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Register the parallel backend to use multicore, in case of foreach loops
registerDoParallel(n.cores)

# Parallel downloading of data
foreach(i = 1:length(a)) %dopar% {
  try(download_fun(a = a, i = i))
}

# Stop the implicit cluster created by foreach
stopImplicitCluster()
