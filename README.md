
本仓库记录了从原始数据获取到Multi Role ChatGPT Framework（MRCF）搭建的全流程代码和用法介绍，专业的开发人员可用于以下用途：

- GEO数据库自动挖掘
- 针对单一领域的深入的转录组学数据挖掘
- 搭建特定领域的MRCF，用于增加ChatGPT回复的准确性和一致性

**note：** 为了代码的简洁性，`x.5 function.R`是`x. *.R`的补充函数或操作。对于GEO数据库分析来说，可以直接开始运行`4. Optimal_classifier.R`。

代码目录为:
```
├─1.code
│      1. GEO_auto_search.R
│      1.5 function.R
│      2. GEOdownload_url.R
│      2.5. GSE_metadata.R
│      3. modify_prompt_code.R
│      3. Prompt_train.R
│      4. Optimal_classifier.R
│      4.5 function.R
│
├─2.supple_data
│      ATC_result.txt
│      qc_prompt.txt
│      requirement.txt
│      verify_prompt_v1.0.txt
│
├─3.auto_array
│  │  auto_array.Rproj
│  │  GPT_batch_limma.R
│  │
│  └─code
│          auto_ann.R
│          batch_download.R
│          error_check.R
│          function.R
│          GEO_auto_download.R
│          limma_Pvalue_batch.R
│          main_code.R
│
└─4.auto_RNA-seq
        DESeq2auto.R
        GEO.R
        RNAseq_download.R
```

# 1.code

## 1. GEO_auto_search.R

此代码用于批量检索关键词对应的数据集，获取GSE号等信息。本研究以所有的ATC drug name为例。

### input

`ATC_result`是一个带有`drug_name`列的data.frame，示例如下：
```
ATC_result <- read.table("ATC_result.txt",sep = "\t",header = T)

> head(ATC_result)
      drug_name ATC_3
1 fluprednidene  D07C
2   antibiotics  D07C
3     fentonium  A03B
4 obiltoxaximab  J06B
5   trifarotene  D10A
6     metformin  A10B
```

在实际运用中，drug_name列可替换为其他检索词

### process

- 数据处理
如果筛选标准和下面一致，则无需修改其他代码：
1. 物种：人/小鼠/大鼠
2. 测序类型：Expression profiling by array or Expression profiling by high throughput sequencing

### output
得到的结果中，主要记录了数据集的`accession`，测序平台`gpl`，测序类型`gdstype`，以及检索词`drug`。当然还有该数据集的`title`和`summary`（因为长度过长而省略）
```
> head(all_df[,c(1,2,5,6)])
  accession      gpl
1 GSE231451 GPL25947
2 GSE160369 GPL20301
3 GSE205530 GPL21103
4 GSE199117 GPL25947
5 GSE218370 GPL24247
6 GSE210906 GPL21626
                                             gdstype          drug
1 Expression profiling by high throughput sequencing "antibiotics"
2 Expression profiling by high throughput sequencing "antibiotics"
3 Expression profiling by high throughput sequencing "antibiotics"
4 Expression profiling by high throughput sequencing "antibiotics"
5 Expression profiling by high throughput sequencing "antibiotics"
6 Expression profiling by high throughput sequencing "antibiotics"
```

## 2. GEOdownload_url.R

根据`GSE accession`批量获取`*_series_matrix.txt.gz`的下载url，使用这些url能够借助IDM等工具实现多线程下载，能够比`GEOquery`包更快。`*_series_matrix.txt.gz`记录了芯片数据的表达谱矩阵以及所有数据集的`sample metadata`

### input

示例数据如下：
```
all_GSE <- read.table(input_file, sep = "\t", header = FALSE)

> head(all_GSE)
         V1
1 GSE105417
2 GSE128649
3  GSE58819
4  GSE37464
5  GSE24278
6  GSE18009
```

### process

基于`GEOquery`包批量获取`*_series_matrix.txt.gz`的下载url

### output

示例输出如下：
```
> head(a_char)
                                                                                                 V1
1 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE105nnn/GSE105417/matrix/GSE105417_series_matrix.txt.gz
2 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE128nnn/GSE128649/matrix/GSE128649_series_matrix.txt.gz
3    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58819/matrix/GSE58819_series_matrix.txt.gz
4    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE37nnn/GSE37464/matrix/GSE37464_series_matrix.txt.gz
5    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE24nnn/GSE24278/matrix/GSE24278_series_matrix.txt.gz
6    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18009/matrix/GSE18009_series_matrix.txt.gz
```

## 2.5 GSE_metadata.R

根据`*_series_matrix.txt.gz`批量获取每个数据集的`sample metadata`。

### input

示例`results_df`输入与`1. GEO_auto_search.R`一致，仅需提供`accession`列即可

### procession

基于`GEOquery`的`getGEO` function来获取`sample metadata`

### output

最终得到`gse_list`

```
> summary(gse_list[1:5])
          Length Class         Mode
GSE172061 1      ExpressionSet S4  
GSE122934 1      ExpressionSet S4  
GSE139547 1      ExpressionSet S4  
GSE76003  1      ExpressionSet S4  
GSE60567  1      ExpressionSet S4 
```
详细信息可在R语言中查看
![](https://pic-wx.oss-cn-beijing.aliyuncs.com/202402051134446.png)

## 3. modify_prompt_code.R

### input

- 1. 金标准文件

```
> train_data <- openxlsx::read.xlsx("GEO_train.xlsx")

> str(train_data)
'data.frame':	3879 obs. of  6 variables:
 $ geo_id   : chr  "GSE466" "GSE11208" "GSE11348" "GSE3586" ...
 $ cell_type: chr  "Muscle - Striated (Skeletal) (MMHCC)" "Ganglioneuroblastoma" "Nose" "Myocardial tissue" ...
 $ ctrl_ids : chr  "GSM4372|GSM4373|GSM4374|GSM4375|GSM4376" "GSM282582|GSM282583|GSM282587|GSM282588|GSM282589" "GSM286649|GSM286650|GSM286651|GSM286655|GSM286656|GSM286657|GSM286658|GSM286659|GSM286660|GSM286661|GSM286662|G"| __truncated__ "GSM82393|GSM82394|GSM82395|GSM82396|GSM82397|GSM82398|GSM82399|GSM82400|GSM82401|GSM82402|GSM82403|GSM82404|GSM"| __truncated__ ...
 $ pert_ids : chr  "GSM4377|GSM4378|GSM4379|GSM4380|GSM4381" "GSM282584|GSM282585|GSM282586|GSM282590|GSM282591|GSM282592" "GSM286646|GSM286647|GSM286648|GSM286652|GSM286653|GSM286654|GSM286667|GSM286668|GSM286669|GSM286679|GSM286680|G"| __truncated__ "GSM82408|GSM82409|GSM82410|GSM82411|GSM82412|GSM82413|GSM82414|GSM82415|GSM82416|GSM82417|GSM82418|GSM82419|GSM82420" ...
 $ type     : chr  "disease" "disease" "disease" "disease" ...
 $ pert_name: chr  "Duchenne muscular dystrophy" "Nicotine addiction" "Rhinovirus infection" "dilated cardiomyopathy" ...
```

- 参数设置

`api_key`: 例如`sk-yj0SfCYBRHGETByP9e1a3b136bA9*********************`
`max_attempts`：最大ChatGPT尝试次数（避免过度消费）
`expected_column_count`：期待的列数（以金标准为参考）
`original_prompt`：用户根据自己经验编写的`prompt`

### process

使用代码设定的流程，将`Original Annotator`生成的结果与金标准比较，用以修改`original prompt`。

### output

当输出结果符合条件后则停止迭代，示例输出如下：

```
> str(optimized_prompts[1:5])
List of 5
 $ GSE89899 : chr "Additional Requirements:1.If the patient is on medication, the type is the drug\n2.single ctrl_ids and pert_ids"| __truncated__
 $ GSE111118: chr "Additional Requirements:1.If the patient is on medication, the type is the drug\n2.single ctrl_ids and pert_ids"| __truncated__
 $ GSE193541: chr "Additional Requirements:1.If the patient is on medication, the type is the drug\n2.single ctrl_ids and pert_ids"| __truncated__
 $ GSE110104: chr "Additional Requirements:1.If the patient is on medication, the type is the drug\n2.single ctrl_ids and pert_ids"| __truncated__
 $ GSE33329 : chr "Additional Requirements:1.If the patient is on medication, the type is the drug\n2.single ctrl_ids and pert_ids"| __truncated__
```

## 4. Optimal_classifier.R

应用`optimal prompt`进行用户特定领域的分析。对于GEO数据库分析来说，可以直接进行这一步。

### input

1. `gse_list`(与`2.5 GSE_metadata.R`一致)
2. `GSE accession`（例如`GSE89899`）
3. 其余参数设置与`3. modify_prompt_code.R`一致

### process

无需修改其他的代码，即可应用`Optimal Annotator`进行实际分析

### output

- 未质控的分类结果

```
> head(final_results_df)
     gse_id   cell_type              ctrl_ids
1 GSE235910 Oral Cavity GSM7511788|GSM7511789
2 GSE235910 Oral Cavity            GSM7511794
3 GSE235910 Oral Cavity GSM7511796|GSM7511797
4 GSE235910 Oral Cavity            GSM7511799
5 GSE235910 Oral Cavity            GSM7511802
6 GSE235910 Oral Cavity GSM7511804|GSM7511807
                          pert_ids type               pert_name
1 GSM7511790|GSM7511792|GSM7511793 drug Afatinib(40mg, 14 days)
2                       GSM7511795 drug Afatinib(40mg, 14 days)
3            GSM7511798|GSM7511801 drug Afatinib(40mg, 14 days)
4                       GSM7511800 drug Afatinib(40mg, 14 days)
5                       GSM7511803 drug Afatinib(40mg, 14 days)
6            GSM7511808|GSM7511813 drug Afatinib(40mg, 14 days)
```

- 质控后的分类结果

```
> head(final_results_df)
     gse_id                     cell_type
1  GSE12343               cultured pineal
2  GSE12343               cultured pineal
3  GSE12343               cultured pineal
4 GSE124457                         Kelly
5 GSE172188 Rheumatoid arthritis synovium
6 GSE172188 Rheumatoid arthritis synovium
                                                ctrl_ids
1                          GSM310055|GSM310056|GSM310057
2                          GSM310055|GSM310056|GSM310057
3                          GSM310055|GSM310056|GSM310057
4                       GSM3533555|GSM3533556|GSM3533557
5 GSM5243597|GSM5243609|GSM5243611|GSM5243613|GSM5243615
6 GSM5243599|GSM5243601|GSM5243603|GSM5243605|GSM5243607
                                                pert_ids type
1                          GSM310058|GSM310059|GSM310060 drug
2                          GSM310061|GSM310062|GSM310063 drug
3                          GSM310064|GSM310065|GSM310066 drug
4                       GSM3533558|GSM3533559|GSM3533560 drug
5 GSM5243598|GSM5243610|GSM5243612|GSM5243614|GSM5243616 drug
6 GSM5243600|GSM5243602|GSM5243604|GSM5243606|GSM5243608 drug
                            pert_name
1                      norepinephrine
2                      dibutyryl-cAMP
3                           forskolin
4                                  P3
5 Abatacept(125mg per week, 16 weeks)
6 Abatacept(125mg per week, 16 weeks)
                                      qc_results    reasons
1 condition 1(+), condition 2(+), condition 3(+) NA, NA, NA
2 condition 1(+), condition 2(+), condition 3(+) NA, NA, NA
3 condition 1(+), condition 2(+), condition 3(+) NA, NA, NA
4 condition 1(+), condition 2(+), condition 3(+) NA, NA, NA
5 condition 1(+), condition 2(+), condition 3(+) NA, NA, NA
6 condition 1(+), condition 2(+), condition 3(+) NA, NA, NA
```

# 2.supple_data

储存了`1. code`中所需的补充文件

# 3.auto_array

microarray自动化的差异表达分析流程

## input

- 1. final_results_df(由`Optimal_classifier.R`生成)
- 2. `*_series_matrix.txt.gz`(由`2. GEOdownload_url.R`生成)

## process

1. 使用`auto_ann.R`将`probe_id`注释为`gene_symbol`
2. 使用`limma_Pvalue_batch.R`，基于limma包对microarray数据集进行批量的差异表达分析

## output

示例输出如下：

```
genes	logFC	AveExpr	t	P.Value	adj.P.Val	B
ADAMTS5	-3.909486083	10.28677279	-17.46052198	9.28E-09	0.000129264	10.34793892
EDIL3	3.069385556	8.399433556	17.07907903	1.15E-08	0.000129264	10.17773547
C10orf96	-2.539028139	8.523583069	-16.01644675	2.12E-08	0.000129264	9.673360857
```

# 4.auto_array

## input
- 1. final_results_df(由`Optimal_classifier.R`生成)
- 2. `*_count.xls`(由`RNAseq_download.R`(GREIN database)或`GEO.R`(GEO database)生成)

## process

使用`DESeq2auto.R`对RNA-seq数据集按照特定分组进行批量的差异表达分析

## output

差异表达分析结果示例如下：

```
genes	logFC	baseMean	lfcSE	stat	P.Value	adj.P.Val
A2M-AS1	0.190352991	7.852281849	0.778205726	0.244604974	0.806762312	0.899349499
A2ML1	-0.403786097	6.578866775	0.791732104	-0.510003441	0.610049051	0.771574343
```

# 运行环境

```
> sessionInfo()
R version 4.3.2 (2023-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=Chinese (Simplified Han)_Hong Kong SAR.utf8 
[2] LC_CTYPE=Chinese (Simplified Han)_Hong Kong SAR.utf8   
[3] LC_MONETARY=Chinese (Simplified Han)_Hong Kong SAR.utf8
[4] LC_NUMERIC=C                                           
[5] LC_TIME=Chinese (Simplified Han)_Hong Kong SAR.utf8    

time zone: Asia/Shanghai
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods  
[7] base     

other attached packages:
[1] Biobase_2.62.0      BiocGenerics_0.48.1

loaded via a namespace (and not attached):
[1] compiler_4.3.2    tools_4.3.2       rstudioapi_0.15.0
[4] Rcpp_1.0.11       stringi_1.8.3     zip_2.3.0        
[7] openxlsx_4.2.5.2
```

----

> 如果在使用中有其他问题，可联系`202101421521@b.sxmu.edu.cn`
