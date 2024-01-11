library(openxlsx)
require(doParallel)
library(stringr)
library(GEOquery)
library(xml2)
library(parallel)
root_dir <- "H:\\chat_project\\drugGPT\\6. auto_array"
code_dir <- paste0(root_dir,"\\code")
file_dir <- paste0(root_dir,"\\file")

setwd(file_dir)
all_GSE <- read.xlsx("single_drug_perturbations-v1.0.xlsx",sheet = 1)
all_GSE <- read.table("gene.txt",sep = "\t",header = F)

setwd(code_dir)
source("function.R")

GEO<- unique(unlist(str_match_all(all_GSE[,1],pattern = "GSE[0-9]*")))
merge<-c()
url<- function(GEO,merge,getDirListing){
	require(doParallel)
	library(stringr)
	library(GEOquery)
	library(xml2)
	library(parallel)
	stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
	gdsurl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/"
	b = getDirListing(sprintf(gdsurl, stub, GEO))
	ret <- list()
	for (x in 1:length(b))
	{
		ret[[x]]<-sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s", 
																				stub, GEO, b[x])
	}
	merge<-c(merge,unlist(ret))
	return(merge)
}


# 16 核
clust <- makeCluster(6)
a <- parLapply(clust, GEO, fun = url,merge,getDirListing)
stopCluster(clust)

output <- unlist(a)
GSE_name <- str_match(output,"GSE[0-9]*_series_matrix.txt.gz|GSE[0-9]*\\-GPL[0-9]*_series_matrix.txt.gz")

# write.table(output,"all_url.txt",sep = "\t",col.names = F,row.names = F,quote = F)

registerDoParallel(6)
foreach(i=1:length(output)) %dopar% try(download.file(output[i],destfile = paste0("E:\\0.HR\\1_tasks\\past\\0.auto_array\\6-20\\series\\",GSE_name[i],"_series_matrix.txt.gz")))
stopImplicitCluster()

#批量下载GPL
library(GEOquery)
library(stringr)

rt <- openxlsx::read.xlsx(xlsxFile = "file/final_results_qced_ann.xlsx",sheet = 1)
rt <- rt[rt$Exp_type=="Expression profiling by array",]
#rt=read.xlsx("GSEid.xlsx",1,header=F)
rt1=rt$GPL
rt1<-unique(rt1)
setwd("./file/GPL/")
finished<-list.files(pattern = "GPL[0-9]*",full.names = F)
finished_name<-str_extract(finished,"GPL[0-9]*")
last_file<-setdiff(rt1,finished_name)
for (i in last_file) {
  success <- FALSE
  while (!success) {
    tryCatch({
      a <- getGEO(i, destdir = '.', getGPL = TRUE, AnnotGPL = TRUE)
      output <- a@dataTable@table
      file_path <- paste(".\\", i, ".txt", sep="")
      write.table(output, file_path, sep = "\t", row.names = FALSE, col.names = TRUE)
      success <- TRUE
    }, error = function(e) {
      cat("Error in downloading:", e$message, "\n")
      cat("Retrying in 10 seconds...\n")
      Sys.sleep(10)  # 等待10秒
    })
  }
}

