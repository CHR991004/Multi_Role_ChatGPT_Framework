group_data <- openxlsx::read.xlsx("./7. auto_RNA-seq/final_results_qced_ann.xlsx")
group_data <- group_data[group_data$Exp_type=="Expression profiling by high throughput sequencing",]
gse_list <- group_data$gse_id
library("dplyr")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("data.table")
library(RCurl)
# 初始化错误和缺失文件的列表
errors_list <- c()
missing_files_list <- c()

# 循环处理每个GSE
for (GSE_number in gse_list) {
  tryCatch({
    # 构建下载路径
    urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
    path <- paste(urld, paste("acc=", GSE_number, sep=""), paste("file=", GSE_number, "_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep=""), sep="&")
    
    # 检查文件是否存在
    if (!url.exists(path)) {
      missing_files_list <- c(missing_files_list, GSE_number)
      next  # 跳过当前循环
    }
    
    # 下载文件
    download.file(path, destfile = paste("./7. auto_RNA-seq/matrix/",GSE_number, "_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep=""))
    
    # 读取和处理数据
    raw_counts <- as.matrix(data.table::fread(paste("./7. auto_RNA-seq/matrix/",GSE_number, "_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep=""), header=T, colClasses="integer"), rownames=1)
    rownames(raw_counts) <- mapIds(org.Hs.eg.db, keys = rownames(raw_counts), column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
    
    raw_counts <- raw_counts %>%
      .[!is.na(rownames(.)), ] %>%
      as.data.frame()
    
    gene_counts <- cbind(genes = rownames(raw_counts), raw_counts)
    write.table(gene_counts, file = paste(GSE_number, "_count.xls", sep=""), sep = "\t", row.names = F, col.names = T)
    
  }, error = function(e) {
    # 记录错误
    errors_list <- c(errors_list, GSE_number)
  })
}
