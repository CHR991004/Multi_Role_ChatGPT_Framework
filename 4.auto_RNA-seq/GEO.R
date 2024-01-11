library("dplyr")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("data.table")

# 1. 自定义输入GSE号
GSE_number <- "GSE71014"

# 2. 分析得到对应的gene_counts
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, paste("acc=", GSE_number, sep=""), paste("file=", GSE_number, "_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep=""), sep="&")
download.file(path, destfile = paste(GSE_number, "_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep=""))

raw_counts <- as.matrix(data.table::fread(paste(GSE_number, "_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep=""), header=T, colClasses="integer"), rownames=1)
rownames(raw_counts) <- mapIds(org.Hs.eg.db, keys = rownames(raw_counts), column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

raw_counts <- raw_counts %>%
  .[!is.na(rownames(.)), ] %>%
  as.data.frame()

gene_counts <- cbind(genes = rownames(raw_counts), raw_counts)
write.table(gene_counts,file = "GSE89774_count.xls",sep = "\t",row.names = F,col.names = T)
