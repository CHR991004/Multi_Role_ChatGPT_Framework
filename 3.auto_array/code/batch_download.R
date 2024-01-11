require(doParallel)
library(stringr)
library(GEOquery)
library(xml2)
library(parallel)
library(openxlsx)

input_file <- "AML_array_group.xlsx"#需修改
setwd("E:\\0.HR\\1_tasks\\1.blood\\blood_fig\\array")#需修改

n.cores <- 8#获得最大核数，或者自行设置

all_GSE<-read.xlsx(input_file,sheet = 1)
GEO<- unique(unlist(str_match_all(all_GSE$GSE,pattern = "GSE[0-9]*")))
getDirListing <- function(url) {
	# Takes a URL and returns a character vector of filenames
	a <- xml2::read_html(url)
	fnames = grep('^G',xml_text(xml_find_all(a,'//a/@href')),value=TRUE)
	return(fnames)
}

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


clust <- makeCluster(n.cores)
a <- parLapply(clust, GEO, fun = url,merge,getDirListing)
stopCluster(clust)
output_url <- unlist(a)
num<-str_split(input_file,".xlsx")[[1]][1]
write.table(output_url,paste0(num,"_outputURL.txt"),sep = "\t",row.names = F,col.names = F,quote = F)

registerDoParallel(n.cores)
foreach(i=1:length(a)) %dopar% try(download.file(a[[i]],destfile = paste0("E:\\0.HR\\1_tasks\\1.blood\\blood_fig\\array\\series\\",GEO[i],"_series_matrix.txt.gz")))
stopImplicitCluster()
