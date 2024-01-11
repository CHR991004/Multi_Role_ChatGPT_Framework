library(limma)
library(impute)
library(stringr)
library(biomaRt)
library(openxlsx)
library(tidyr)
rm(list = ls())
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")
mus = useMart(biomart="ensembl",dataset="mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org")
rat = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host="https://dec2021.archive.ensembl.org")
# saveRDS(object = ensembl,file = "mus_ensembl.RDS")
# saveRDS(object = human,file = "human_ensembl.RDS")
# ensembl <- readRDS("E:\\0.HR\\1_tasks\\2021\\herb\\mus_ensembl.RDS")
# human <- readRDS("E:\\0.HR\\1_tasks\\2021\\herb\\human_ensembl.RDS")
logFoldChange=0
adjustP=0.05

infor_matrix <- read.xlsx("./file/final_results_qced_final_ann.xlsx")
matrix_now <- list.files("E:\\0.HR\\1_tasks\\past\\0.auto_array\\D2D\\D2D_matrix")
matrix_ID <- str_split(matrix_now,"_",simplify = T)[,1]
output_name <- str_split(matrix_now,".txt",simplify = T)[,1]

setwd("E:\\0.HR\\1_tasks\\past\\0.auto_array\\D2D\\D2D_matrix")#需修改
for (i in 1316:length(matrix_ID)) {
		HBEXP_ID <- matrix_ID[i]
		infor_location <- as.numeric(grep(pattern = paste0("^",HBEXP_ID,"$"),x = rownames(infor_matrix)))
		Herb_name <- infor_matrix$perturbation[infor_location]
		con_group <- str_split(string = infor_matrix$control_GSM[infor_location] ,pattern = "[;]",simplify = T)[1,]
		con_group <-  con_group[con_group!=""]
		drug_group <- str_split(string = infor_matrix$test_GSM[infor_location] ,pattern = "[;]",simplify = T)[1,] 
		drug_group <-  drug_group[drug_group!=""]
		if ((length(con_group)==1)|(length(drug_group)==1)) {
		  next
		}
		setwd("E:\\0.HR\\1_tasks\\past\\0.auto_array\\D2D\\D2D_matrix")
		if (infor_matrix$organism[infor_location]=="Homo sapiens") {
		  rt=read.table(matrix_now[i],sep="\t",header=T,check.names = F)
		  rt <- na.omit(rt)
		  rt<-rt[,colnames(rt)!="UniGene symbol"]
		}else if(infor_matrix$organism[infor_location]=="Mus musculus"){
			rt=read.table(matrix_now[i],sep="\t",header=T,check.names = F)
			rt <- na.omit(rt)
			rt<-rt[,colnames(rt)!="UniGene symbol"]
			genes <- rt$genes
			geneENS_mouse<- getLDS(attributes=c("mgi_symbol","hgnc_symbol"),
																										attributesL = "hgnc_symbol",
																										values=genes,
																										mart=mus,martL = human,uniqueRows=T)
			rt = merge(geneENS_mouse,rt,by.x ="MGI.symbol",by.y="genes")
			rt <- rt[,3:ncol(rt)]
		
		}else{
			rt=read.table(matrix_now[i],sep="\t",header=T,check.names = F)
			rt <- na.omit(rt)
			rt<-rt[,colnames(rt)!="UniGene symbol"]
			genes <- rt$genes
			geneENS_rat<- getLDS(attributes = c("rgd_symbol"), 
																										filters = "rgd_symbol", 
																										values = genes , 
																										mart = rat, 
																										attributesL = c("hgnc_symbol"), 
																										martL = human, 
																										uniqueRows=T)
			rt = merge(geneENS_rat,rt,by.x ="RGD.symbol",by.y="genes")
			rt <- rt[,2:ncol(rt)]
		}
		rt = as.matrix(rt)
		rownames(rt) = rt[,1]
		exp =rt[,2:ncol(rt)]
		dimnames = list(rownames(exp),colnames(exp))
		exp = matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
		exp <- exp[,c(drug_group,con_group)]
		exp <- na.omit(exp)
		exp <- subset(exp, rowMeans(exp) > 0)
		# exp <- exp[apply(exp,1,min)>0,]
		# exp <- exp[apply(exp,1,sd)>0,]

		#impute missing expression data
		mat=impute.knn(exp)
		rt=mat$data
		rt=avereps(rt)
		#normalize
		setwd("E:\\0.HR\\1_tasks\\past\\0.auto_array\\D2D\\diff")
		rt=normalizeBetweenArrays(as.matrix(rt))
		
		qx <- as.numeric(quantile(rt, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
		LogC <- (qx[5] > 100) ||
			(qx[6]-qx[1] > 50 && qx[2] > 0)
		if (LogC) { rt[which(rt <= 0)] <- NaN
		rt <- log2(rt) }
		#differential
		#class <- c("con","con","treat","con","treat","treat")
		class <- c(rep("ITP",length(drug_group)),rep("nom",length(con_group)))    #??Ҫ?޸?
		design <- model.matrix(~0+factor(class))
		colnames(design) <- c("ITP","nom")
		fit <- lmFit(rt,design)
		cont.matrix<-makeContrasts(ITP-nom,levels=design)
		fit2 <- contrasts.fit(fit, cont.matrix)
		fit2 <- eBayes(fit2)
		
		allDiff=topTable(fit2,adjust='fdr',number=200000)
		#write table
		logFoldChange = mean(abs(allDiff$logFC))+2*sd(abs(allDiff$logFC))
		diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & P.Value < adjustP )), ]
		diffSig <- cbind(rownames(diffSig),diffSig)
		colnames(diffSig)[1]<-"genes"
		diffSig<-separate_rows(data = diffSig,genes,sep = "///")
		write.table(diffSig,file=paste(output_name[i],".xls",sep = ""),sep="\t",quote=F,row.names = F)
		
}

