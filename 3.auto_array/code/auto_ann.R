library(dplyr)
library(stringr)
library(limma)
library(rvest)
library(GEOquery)
library(data.table)
rm(list = ls())
root <- "H:\\chat_project\\drugGPT\\6. auto_array\\file"
array_series_dir <- "H:\\chat_project\\drugGPT\\6. auto_array\\file\\series"
array_plat_dir <- "H:\\chat_project\\drugGPT\\6. auto_array\\file\\GPL"
output_dir <- "H:\\chat_project\\drugGPT\\6. auto_array\\file\\genematrix"
other_ID<-c()

setwd(root)
array <- openxlsx::read.xlsx("final_1215.xlsx",sheet = 1)
no_matrix <- c()
for (i in 1:nrow(array)) {
  print(i)
  drug_ID <- array$drug_id[i]
  # drug_ID <- str_replace_all(drug_ID,pattern = "_NA",replacement = "")
  # drug_ID <- str_replace_all(drug_ID,pattern = " ",replacement = "_")
  GSE_ID <- array$gse_id[i]
  GSM_one <- str_match(string = array$ctrl_ids[i],pattern = "(GSM[0-9]*).*")[,2]
  GPL_ID <- array$GPL[i]
  setwd(array_series_dir)
  all_series_file <- list.files(array_series_dir,".gz")
  series_file_name <- grep(pattern = paste0("^",GSE_ID),x = all_series_file,value = T)
  if (length(grep("GPL",series_file_name)) != 0) {
    file_location <- as.numeric(grep(GPL_ID,series_file_name))
  }else{
    file_location <- 1
  }
  if (length(file_location)==0) {
    series_file_name <- grep(pattern = paste0("^", GSE_ID, "_series_matrix\\.txt\\.gz$"), x = all_series_file, value = TRUE)
  }
  
  GSE_array<-try(GEOquery:::parseGSEMatrix(series_file_name[file_location], destdir = array_series_dir, 
                                 AnnotGPL = F, getGPL = F)[["eset"]]@assayData[["exprs"]],silent = T)
  if ("try-error" %in% class(GSE_array)) {
    series_matrix<-fread(series_file_name,sep = "\t",header = T,
                          skip = "!series_matrix_table_begin",quote = "")
  }else{
    series_matrix <- as.data.frame(cbind(ID_REF=rownames(GSE_array),GSE_array))
  }
  if(nrow(series_matrix)==0){
    no_matrix <- c(no_matrix, GSE_ID);next
  }
  colnames(series_matrix)<-gsub('"',"",colnames(series_matrix))
  series_matrix$ID_REF<-gsub('"',"",series_matrix$ID_REF)
  # series_matrix$ID_REF<-str_match(series_matrix$ID_REF,pattern = "[0-9]+")
  setwd(array_plat_dir)
  all_gpl_file <- list.files(array_plat_dir,".txt")
  gpl_file_name <- grep(pattern = GPL_ID,x = all_gpl_file,value = T)
  if(length(gpl_file_name) == 0){
    other_ID<c(other_ID,GSE_ID);next
  }
  gpl_matrix <- read.table(gpl_file_name,sep = "\t",header = T,check.names = F,comment.char = "#",fill = T,quote = "")
  symbol_col <- try(grep(pattern = "symbol",
                         colnames(gpl_matrix),
                         ignore.case = T,value = T))
  if(length(symbol_col) == 0){
    other_ID<c(other_ID,GSE_ID);next
  }else if(length(symbol_col) != 1){
    symbol_col <- symbol_col[-grep(pattern = "UniGene symbol",
         symbol_col,
         ignore.case = T,value = F)]
  }
  ann <- gpl_matrix[,c("ID",symbol_col)]

  geneMatrix<-merge(ann,series_matrix,by.x = "ID",by.y = "ID_REF")[,-1]
  temMatrix<-as.matrix(geneMatrix[geneMatrix[,symbol_col]!="",])
  rownames(temMatrix)<-temMatrix[,1]
  geneMatrix_rep_out<-temMatrix[,-1]%>%avereps()%>%as.data.frame()
  geneMatrix_output<-cbind(genes=rownames(geneMatrix_rep_out),geneMatrix_rep_out)
  setwd(output_dir)
  write.table(geneMatrix_output,paste(drug_ID,"_gene_matrix.txt",sep = ""),sep = "\t",row.names = F,col.names = T,quote = F)

}

