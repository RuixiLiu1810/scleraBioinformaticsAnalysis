#代码持续更新
#代码更新百度网盘：通过网盘分享的文件：GEO代码更新
# 链接: https://pan.baidu.com/s/1w9C8-_766G2FPd8CO5m4jg?pwd=GEOO 提取码: GEOO 
# --来自百度网盘超级会员v6的分享

#一、R包下载方法####
options(timeout = 88888888)
#1.R官方
install.packages('glue') #安装 必须加引号（示例R包，可以不装）
library(glue) #加载 可以不加引号
remove.packages('glue') #删除

# https://cran.r-project.org/ 
#手动安装:Software > Packages > Table of available packages, sorted by name > ctrl+F 搜索
#https://cran.r-project.org/src/contrib/Archive

#2.BiocManager
install.packages('BiocManager') 
BiocManager::install('GEOquery')
library(GEOquery)

#手动：直接搜索GEOquery关键词进BiocManager官网下载tar.gz
#https://bioconductor.org/packages/release/BiocViews.html#___Software


#3.github
# https://cran.r-project.org/bin/windows/Rtools/   Rtools对C++包编译
install.packages('devtools') 
library(devtools)
#（示例R包，不用装）
devtools::install_github("schymane/ReSOLUTION") 

#手动：直接搜索R包名，进入github，Releases/Packages 下载tar.gz文件
#手动：code > zip下载 > new project > install
devtools::install_local("E:/R_do/临时文件/ReSOLUTION-master.zip") 


#开始前，需要安装的R包
install.packages('tidyverse')
install.packages('BiocManager')
BiocManager::install("GEOquery") #先安装BiocManager
install.packages("limma")


#二、GEO芯片数据整理####
#GEO官网：https://www.ncbi.nlm.nih.gov/
#下载Series Matrix File(s)文件 及Platforms：GPL文件
#清空环境
rm(list=ls())
# BiocManager::install("GEOquery")
library(GEOquery)
library(tidyverse)
setwd('E:/R_do/阿伦生信教学代码/GEO/GEO数据下载与整理/') #设置工作路径
gset <- getGEO("GSE299988", destdir = ".", AnnotGPL = F, getGPL = F)#手动下载，自动读取
class(gset) 
gset[[1]] #查看信息
#提取基因表达数据  
exp <- exprs(gset[[1]])
#1.标准化####
boxplot(exp, outline=F,las=2)
library(limma)
exp <- normalizeBetweenArrays(exp)
boxplot(exp, outline=F,las=2)
#2.归一化####
range(exp)#查看range范围，如果数值过大，需要归一化，log2处理，如果<100,不做log处理
exp <- log2(exp+1)
range(exp)
range(exp,na.rm = T) # 有NA值，运行此代码
#3.基因注释，GPL导入####
gpl<-read.table('GPL21185-21174.txt',header=TRUE,fill=T,sep='\t',comment.char='#',stringsAsFactors=FALSE,
                quote='')  #去掉前几列没用的行
View(gpl)
colnames(gpl) #查看列名
ids <- gpl[,c('ID','GENE_SYMBOL')]  #需要查看表格填写正确的列名
colnames(ids) <- c('probe_id','symbol') #重命名
library(stringr)
ids$symbol <- trimws(str_split(ids$symbol,'//',simplify = T)[,1]) #确认是//还是///，需要提取第几个元素，#trimws 去除字符前后的空格
ids <- ids[ids$symbol!='',]#删除空白
ids <- ids[ids$symbol!='---',]#删除空白
ids <- ids[ids$symbol!='--- ',]#删除空白
ids <- ids[ids$probe_id %in% rownames(exp),]
ids$probe_id <- as.character(ids$probe_id) #确保probe_id列是字符型
exp2 <- exp[ids$probe_id,]
table(rownames(exp2)==ids$probe_id)
exp3 <- cbind(ids,exp2)
#4.基因去重####
exp3 <- aggregate( . ~ symbol,exp3[,-1],max) #max最大，或者mean平均 随你
rownames(exp3) <- exp3$symbol
exp3 <- exp3[,-c(1)]

#5.获取临床数据####
rt1 <- pData(gset[[1]])
#保存总临床数据表
write.csv(rt1, file='GSE299988_cli.csv',row.names = T)
#对样本进行分组
#1.疾病分组
colnames(rt1)
group_list <- ifelse(rt1$`tissue types:ch1`== 'Normal' , "Normal","Tumor")
# group_list <- ifelse(str_detect(rt1$`tissue types:ch1`, "Nor"), "Normal","Tumor") #方法二：字符选择法
# group_list <- ifelse(str_detect(rt1$`title`, "normal_adjace"), "Normal",ifelse(str_detect(rt1$`title`, "positiveLN_RAIrefractiv"), "P_Tumor","N_Tumor"))
group_list
rt1$group <- group_list
rt2 <- rt1 %>% dplyr::select("group")
#2.分组排列
rt3 <- rt2 %>% arrange(rt2) 
#如果疾病在前面就用desc，确保对照组在前，疾病在后
# rt3 <- rt2 %>% arrange(desc(group)) #desc 倒置

#6.数据保存####
#exp3和rt3 列行一一对应
exp3 <- exp3[,rownames(rt3)]
identical(rownames(rt3),colnames(exp3))
# [1] TRUE
#数据保存
save(rt3,exp3,file = "GSE299988.rda")

#确认
rm(list=ls())
load( "GSE299988.rda")


#三、GEO芯片数据（非symbol基因名）处理####
# 举例：GSE200539
rm(list=ls())
library(GEOquery)
library(tidyverse)
setwd('E:/R_do/阿伦生信教学代码/GEO/GEO数据下载与整理/') #设置工作路径
gset <- getGEO("GSE200539", destdir = ".", AnnotGPL = F, getGPL = F)#手动下载，自动读取
class(gset) 
gset[[1]] #查看信息
#提取基因表达数据  
exp <- exprs(gset[[1]])
#1.标准化####
boxplot(exp, outline=F,las=2)
# library(limma)
# exp <- normalizeBetweenArrays(exp)
# boxplot(exp, outline=F,las=2)
#2.归一化####
range(exp)#查看range范围，如果数值过大，需要归一化，log2处理
# exp <- log2(exp+1)
# range(exp)

#3.基因注释，soft导入####
GPL <- getGEO(filename = 'GPL23159_family.soft.gz',destdir = '.')
#获取注释信息
gpl <- GPL@dataTable@table
View(gpl)
colnames(gpl) #查看列名

#3.基因注释，GPL导入####
gpl<-read.table('GPL23159-184565.txt',header=TRUE,fill=T,sep='\t',comment.char='#',stringsAsFactors=FALSE,
                quote='')  #去掉前几列没用的行
View(gpl)
colnames(gpl) #查看列名
ids <- gpl[,c('ID','SPOT_ID.1')]  #需要查看表格填写正确的列名
colnames(ids) <- c('probe_id','symbol') #重命名
library(stringr)
ids$symbol <- trimws(str_split(ids$symbol,'//',simplify = T)[,1]) #确认是//还是///，需要提取第几个元素，#trimws 去除字符前后的空格
#转换
library(org.Hs.eg.db)
library(AnnotationDbi)
symbol_ids <- mapIds(
  x = org.Hs.eg.db,
  keys = ids$symbol,       # 输入需要转化的列
  keytype = "REFSEQ",    # 明确指定输入类型
  column = "SYMBOL",  # 转换为目标类型：ENTREZID、ALIAS、REFSEQ、ENSEMBL、SYMBOL、
  multiVals = "first"   # 多个匹配时取第一个
)
symbol_ids
ids$symbol <- symbol_ids
ids <- na.omit(ids) #删除NA

ids$probe_id <- as.character(ids$probe_id)#保证探针列为字符型

exp2 <- exp[ids$probe_id,]
table(rownames(exp2)==ids$probe_id)
exp3 <- cbind(ids,exp2)
#4.基因去重####
exp3 <- aggregate( . ~ symbol,exp3[,-1],max) #max最大，或者mean平均 随你
rownames(exp3) <- exp3$symbol
exp3 <- exp3[,-c(1)]

#5.获取临床数据####
rt1 <- pData(gset[[1]])
#保存总临床数据表
write.csv(rt1, file='GSE200539_cli.csv',row.names = T)
#对样本进行分组
#1.疾病分组
colnames(rt1)
group_list <- ifelse(str_detect(rt1$`title`, "CTL"),"CTL",
                     ifelse(str_detect(rt1$`title`, "IDMF"),"IDMF",
                            ifelse(str_detect(rt1$`title`, "DRF"),"DRF",
                            'MMF')))
group_list
rt1$group <- group_list
rt2 <- rt1 %>% dplyr::select("group")
#2.分组排列
rt3 <- rt2 %>% arrange(rt2) 
#如果疾病在前面就用desc，确保对照组在前，疾病在后
# rt3 <- rt2 %>% arrange(desc(group)) #desc 倒置

#6.数据保存####
#exp3和rt3 列行一一对应
exp3 <- exp3[,rownames(rt3)]
identical(rownames(rt3),colnames(exp3))
# [1] TRUE
#数据保存
save(rt3,exp3,file = "GSE200539.rda")

#确认
rm(list=ls())
load( "GSE200539.rda")


#四、Affymetrix原始芯片数据处理####
setwd('E:/R_do/阿伦生信教学代码/GEO/Affymetrix原始数据处理/') #设置工作路径
rm(list=ls())
library(tidyverse)
# BiocManager::install("oligo")
library(oligo)
#1.数据导入####
FileName <- list.files(path = './GSE200539_RAW')
FileName
sample_id <- str_sub(FileName,1,10)
sample_id
rt <- data.frame(FileName=FileName,row.names =sample_id,
                        group=rep(c('CTL','IDMF','DRF','MMF'),each=3))

cel_files <- list.celfiles( './GSE200539_RAW', full.names = TRUE,listGzipped = T)  # 获取文件列表
cel_files
cel_data <- read.celfiles(cel_files)  # 读取数据
#2.质控处理####
#BiocManager::install('arrayQualityMetrics')
library(arrayQualityMetrics)
library(GEOquery) 
expr <- exprs(cel_data)
colnames(expr) <- sample_id
eset <- ExpressionSet(assayData = expr,phenoData = AnnotatedDataFrame(rt))
arrayQualityMetrics(eset,
                    outdir='GSE200539_QulityMetrics',
                    force=T,do.logtransform=T,intgroup = 'group')
#rma
eset_rma <- oligo::rma(cel_data)
class(eset_rma) 
#[1] "ExpressionSet"
exp <- exprs(eset_rma)
head(exp)
colnames(exp) <- str_sub(colnames(exp),1,10)
dim(exp)
boxplot(exp, outline=FALSE, notch=T,las=2)
range(exp)
#3.基因注释，GPL导入####
gpl<-read.table('GPL23159-184565.txt',header=TRUE,fill=T,sep='\t',comment.char='#',stringsAsFactors=FALSE,
                quote='')  #去掉前几列没用的行
View(gpl)
colnames(gpl) #查看列名
ids <- gpl[,c('ID','SPOT_ID.1')]  #需要查看表格填写正确的列名
colnames(ids) <- c('probe_id','symbol') #重命名
library(stringr)
ids$symbol <- trimws(str_split(ids$symbol,'//',simplify = T)[,1]) #确认是//还是///，需要提取第几个元素，#trimws 去除字符前后的空格
#转换
library(org.Hs.eg.db)
library(AnnotationDbi)
symbol_ids <- mapIds(
  x = org.Hs.eg.db,
  keys = ids$symbol,       # 输入需要转化的列
  keytype = "REFSEQ",    # 明确指定输入类型
  column = "SYMBOL",  # 转换为目标类型：ENTREZID、ALIAS、REFSEQ、ENSEMBL、SYMBOL、
  multiVals = "first"   # 多个匹配时取第一个
)
symbol_ids
ids$symbol <- symbol_ids
ids <- na.omit(ids) #删除NA

exp2 <- exp[ids$probe_id,]
table(rownames(exp2)==ids$probe_id)
exp3 <- cbind(ids,exp2)
#4.基因去重####
exp3 <- aggregate( . ~ symbol,exp3[,-1],max) #max最大，或者mean平均 随你
rownames(exp3) <- exp3$symbol
exp3 <- exp3[,-c(1)]

#5.数据保存####
save(exp3,file = "GSE200539_raw.rda")

exp33 <- exp3
load( "GSE200539.rda")

#五、Illumina原始数据处理####
#清空环境
rm(list=ls())
# BiocManager::install("GEOquery")
library(GEOquery)
library(tidyverse)
library(limma)
setwd('E:/R_do/阿伦生信教学代码/GEO/ilumina原始数据处理/') #设置工作路径
gset <- getGEO("GSE68004", destdir = ".", AnnotGPL = F, getGPL = F)#手动下载，自动读取
class(gset) 
gset[[1]] #查看信息
#提取基因表达数据  
exp <- exprs(gset[[1]])
#1.标准化####
boxplot(exp, outline=F,las=2)
library(limma)
exp <- normalizeBetweenArrays(exp)
boxplot(exp, outline=F,las=2)
#2.归一化####
range(exp)#查看range范围，如果数值过大，需要归一化，log2处理
range(exp,na.rm = T)
# [1]   -22.82553 17174.43372


#1.观察数据结构####
dat1 <- data.table::fread('GSE68004_non-normalized_ProbeRowSignalDetectionDataset_5383_143523_1.txt')
dat2 <- data.table::fread('GSE68004_non-normalized_ProbeRowSignalDetectionDataset_5383_143529_1.txt')
#2.读取数据####
x1 <- read.ilmn(files = 'GSE68004_non-normalized_ProbeRowSignalDetectionDataset_5383_143523_1.txt',
                expr = '750',#数据列名 相同的内容
                other.columns='Detection', #p值
                probeid='ID_REF') #芯片列名


x2 <- read.ilmn(files = 'GSE68004_non-normalized_ProbeRowSignalDetectionDataset_5383_143529_1.txt',
                expr = '586',#数据列名 相同的开头
                other.columns='Detection', #p值
                probeid='ID_REF') #芯片列名
#合并
pro_id <- intersect(rownames(x1),rownames(x2))
x1 <- x1[pro_id,]
x2 <- x2[pro_id,]
x <- cbind(x1, x2)

#3.校准化-log处理####
y <- neqc(x) 
#4.提取表达数据####
exp <- y$E
boxplot(exp,las=3,outline=F)
range(exp)
#[1]  4.161602 14.037369
#去批次
library(sva)
batch <- c(rep(1, ncol(x1)), rep(2, ncol(x2))) #输入批次信息
exp_combat <- ComBat(dat = exp, batch = batch, mod = NULL, par.prior = TRUE)
boxplot(exp_combat,las=3,outline=F)
range(exp_combat)
#[1]  3.807521 14.370650
#5.整理临床数据
rt <- pData(gset[[1]])
colnames(rt)
rt <- dplyr::select(rt,'title','source_name_ch1')
sample_id <- substring(rt$title,nchar(rt$title)-8)#提取字符-8就是最后8个
sample_id
rownames(rt) <- sample_id

rt$group <- trimws(str_split(rt$source_name_ch1,',',simplify = T)[,2])
rt <- dplyr::select(rt,group)

#保持行名列名一致
exp_combat <- exp_combat[,rownames(rt)]
identical(rownames(rt),colnames(exp_combat))
#[1] TRUE
save(exp_combat,rt,file = 'exp_combat.rad')

#质控检查
rt1 <- rt[colnames(x$E), , drop = FALSE]
library(arrayQualityMetrics)
colnames(x$E) <- rownames(rt1)
eset <- ExpressionSet(assayData = x$E,phenoData = AnnotatedDataFrame(data = rt1))
arrayQualityMetrics(eset,
                    outdir='GSE68004_quality',
                    force=T,
                    intgroup='group',
                    do.logtransform=T)


rm(list=ls())
load('exp_combat.rad')

#六、GEO高通量测序数据整理####
rm(list=ls())
setwd('E:/R_do/阿伦生信教学代码/GEO/GEO高通量测序/') #设置工作路径
library(tidyverse)
exp <- data.table::fread("GSE280609_sclerarna_count.csv", header = T)
table(is.na(exp)) #查看是否有缺失值na
exp[is.na(exp)] <- 0
table(is.na(exp)) #查看是否有缺失值na
exp <- column_to_rownames(exp,var = 'Gene')

#Count 去重
exp <- exp[exp$gene_id != "" & !is.na(exp$gene_id), ]       #去除空基因名
exp <- rowsum(as_data_frame(exp[,-1]), group = exp$gene_id) #以gene_id分组求和，实现去重

#质量控制
colSums(exp)       #基本规模
colMeans(exp == 0) #零表达比例


#1.设置分组信息####
sample_id <- colnames(exp)
sample_id
group <- rep(c('control','myopia'),each=5)
group
rt <- data.frame(row.names =sample_id,group=group)
identical(rownames(rt),colnames(exp))
# [1] TRUE

#2.Dseq2分析####
# BiocManager::install('DESeq2')
library(DESeq2)
#创建DESeq2分组对象
rt$group <- factor(rt$group)
rt$group
# 构建DESeq2对象
dds <- DESeqDataSetFromMatrix(
  countData = exp,
  colData = rt,
  design = ~ group
)
#进行低表达基因过滤(自定义)
filter <- rowSums(counts(dds)) > 1  # 筛选每行加起来大于1的基因
dds <- dds[filter, ]
#差异表达分析
dds <- DESeq(dds)
#指定具体的对比组，格式为c("因子名", "实验组", "对照组")
res <- results(dds, contrast = c("group", "myopia", "control"))
#提取差异表达结果
res1 <- as.data.frame(res)
write.table(res1, file='DEG1.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
#结果筛选padj<0.05 & abs(log2FoldChange)>=1
degs <- filter(res1,padj<0.05 & abs(log2FoldChange)>=0.5)#筛选padj<0.05及abs(log2FoldChange)>=1的基因
view(degs)

#绘制离散图
plotDispEsts(dds)
#p值直方图,直观看下p值显著基因的多少
hist(res$padj)
#2.1火山图####
library(ggplot2)
library(cowplot)
library(ggrepel)
#自定义log2FoldChange     
log2FoldChange <- 1
type1 = (res1$padj < 0.05)&(res1$log2FoldChange < -log2FoldChange)
type2 = (res1$padj < 0.05)&(res1$log2FoldChange > log2FoldChange)
res1$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
res1 <- res1[!is.na(res1$change),]
write.table(res1, file='DEG2.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
table(res1$change)
# DOWN   NOT    UP 
# 29 42737    58
DEG <- mutate(res1,Gene_symbol=row.names(res1))
UPDEG <- subset(DEG,change=='UP')
UPDEG_5 <- top_n(x = UPDEG,n = -5,wt = padj)
DOWNDEG <- subset(DEG,change=='DOWN')
DOWNDEG_5 <- top_n(x = DOWNDEG,n = -5,wt = padj)

p <- ggplot(data = DEG,
            aes(x = log2FoldChange,
                y = -log10(padj))) +
  geom_point(alpha = 0.5, size = 4.5,
             aes(color = change)) +
  ylab("-log10(padj)") +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), lty = 5, col = "black", lwd = 0.5) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5) +
  theme_half_open() +
  # 使用 geom_text_repel 避免文本重叠，并设置文本在点上方
  geom_text_repel(data = UPDEG_5, aes(label = Gene_symbol), vjust = 1.5, size = 3,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.35, "lines"),
                  segment.color = 'grey50',
                  nudge_y = 0.1,  # 垂直向上移动文本
                  direction = "y") +  # 限制文本只在垂直方向排列
  geom_text_repel(data = DOWNDEG_5, aes(label = Gene_symbol), vjust = 1.5, size = 3,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.35, "lines"),
                  segment.color = 'grey50',
                  nudge_y = 0.1,  # 垂直向上移动文本
                  direction = "y")  # 限制文本只在垂直方向排列
p
pdf(file = 'volcano.pdf',width = 6,height = 6)
p
dev.off()

#3.count数据vst转化####
vst <- vst(dds,blind = F)
exp1 <- assay(vst)
range(exp1)
exp1 <- as.data.frame(exp1)
save(exp,exp1,rt,file = 'GSE280609.rda')
rm(list=ls())
load('GSE280609.rda')



#七、limma差异分析####
rm(list = ls())
setwd('E:/R_do/阿伦生信教学代码/GEO/差异分析/') #设置工作路径
library(limma)
library(pheatmap)
library(tidyverse)
library(ggpubr)
#数据读取
load('GSE299988.rda')
rt <- rt3
exp <- exp3
str(exp) 
exp <- exp %>% mutate(across(everything(), as.numeric))
str(exp)
#确认行名列名一致#确认行名列名一致exp33
identical(rownames(rt),colnames(exp))
#1.火山图####
group_list <- factor(rt$group,levels = c("Normal","Tumor"))
group_list
design <- model.matrix(~group_list)
#比较矩阵命名
design
#线性模型拟合
fit <- lmFit(exp,design)
#贝叶斯检验
fit2 <- eBayes(fit)
deg <- topTable(fit2, coef = 2, number = Inf)
DEG = na.omit(deg) #differently expressed genes
write.table(DEG, file='DEG1.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
#进行注释 fold change
#通常logFC设置1,adj.P.Val < 0.05
logFC_cutoff <- 1
type1 = (DEG$adj.P.Val < 0.05)&(DEG$logFC < -logFC_cutoff)
type2 = (DEG$adj.P.Val < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
write.table(DEG, file='DEG2.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
table(DEG$change)

library(ggplot2)
library(cowplot)
library(ggrepel)
DEG <- mutate(DEG,Gene_symbol=row.names(DEG))
UPDEG <- subset(DEG,change=='UP')
UPDEG_5 <- top_n(x = UPDEG,n = -5,wt = P.Value)
DOWNDEG <- subset(DEG,change=='DOWN')
DOWNDEG_5 <- top_n(x = DOWNDEG,n = -5,wt = P.Value)
p <- ggplot(data = DEG,
            aes(x = logFC,
                y = -log10(adj.P.Val))) +
  geom_point(alpha = 0.5, size = 4.5,
             aes(color = change)) +
  ylab("-log10(adj.P.Val)") +
  scale_color_manual(values = c("#1F77B4", "grey", "#FF7F0E")) +
  geom_vline(xintercept = c(-1, 1), lty = 5, col = "black", lwd = 0.5) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.5) +
  theme_half_open() +
  # 使用 geom_text_repel 避免文本重叠，并设置文本在点上方
  geom_text_repel(data = UPDEG_5, aes(label = Gene_symbol), vjust = 1.5, size = 3,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.35, "lines"),
                  segment.color = 'grey50',
                  nudge_y = 0.1,  # 垂直向上移动文本
                  direction = "y") +  # 限制文本只在垂直方向排列
  geom_text_repel(data = DOWNDEG_5, aes(label = Gene_symbol), vjust = 1.5, size = 3,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.35, "lines"),
                  segment.color = 'grey50',
                  nudge_y = 0.1,  # 垂直向上移动文本
                  direction = "y")  # 限制文本只在垂直方向排列
p

pdf(file = 'volcano.pdf',width = 10,height = 6)
p
dev.off()

#2.热图####
library(pheatmap)
identical(rownames(rt),colnames(exp))
DEG <- read.table("DEG2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
table(DEG$change)

diff <- rownames(DEG)[DEG$change !="NOT"]
exp_diff <- exp[diff,]
range(exp_diff)
#作图，如果文件损坏，dev.off() 运行至null后再次运行
pdf(file = 'heatmap.pdf',width=6,height=20)
pheatmap(exp_diff,
         annotation_col=rt,
         color = colorRampPalette((c("blue","white","red")))(100),
         cluster_cols = F,
         cluster_rows = T,
         show_colnames = F,
         show_rownames = T ,
         border_color = NA,
         scale='row') #scale='row' 进行标准化
dev.off()
#特定展示
library(tidyverse)
diff =  DEG[DEG$change !="NOT",]
up <- diff %>% top_n(25,logFC)
dw <- diff %>% top_n(-25,logFC)
all <-c(rownames(up),rownames(dw))
exp_diff2 <- exp[all,]
pdf(file = 'heatmap2.pdf',width=6,height=8)
pheatmap(exp_diff2,
         annotation_col=rt,
         color = colorRampPalette((c("blue","white","red")))(100),
         cluster_cols = F,
         cluster_rows = F,
         show_colnames = F,
         show_rownames = T,
         border_color = NA,
         scale='row') #scale='row' 进行标准化
dev.off()


#八、富集分析####
rm(list = ls())
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
setwd('E:/R_do/阿伦生信教学代码/GEO/差异分析/')#设置工作路径
DEG <- read.table("DEG2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
table(DEG$change)

diff <- rownames(DEG)[DEG$change!='NOT']

#将基因ID从Symbol转换为EntrezID
gene_entrez <- bitr(diff, fromType ="SYMBOL", toType ="ENTREZID", OrgDb = org.Hs.eg.db)
#可以查看一下转换比
cat("所有基因ID转换", nrow(gene_entrez),"/", length(diff))
untrans <- setdiff(diff, gene_entrez$SYMBOL)
untrans
isValid <- "C9orf92" %in% keys(org.Hs.eg.db, keytype = "SYMBOL")
isValid #查看T还是F

# 使用mapIds函数转换
library(AnnotationDbi)
entrez_ids <- mapIds(
  x = org.Hs.eg.db,
  keys = untrans,       # 输入的基因别名列表
  keytype = "ALIAS",    # 明确指定输入类型为别名
  column = "ENTREZID",  # 转换为ENTREZID作为中间步骤
  multiVals = "first"   # 多个匹配时取第一个
)
entrez_df <- data.frame(ENTREZID = entrez_ids, check.names = FALSE)
entrez_df <- na.omit(entrez_df)
entrez_df <- entrez_df %>% rownames_to_column(var = 'SYMBOL')

#合并
gene_entrez2 <- rbind(gene_entrez,entrez_df)

#1.GO富集####
GO <- enrichGO(gene = gene_entrez2$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType ="ENTREZID",
               ont ="ALL", # "BP", "MF", "CC"或"ALL"
               pAdjustMethod ="BH",
               pvalueCutoff =0.4,
               qvalueCutoff =0.4,readable = T)

GO_res <- GO@result
write.table(GO_res,file="GO_res.txt",sep="\t",quote=F,row.names = F)
#简单作图
barplot(GO, showCategory = 20)
dotplot(GO, showCategory = 20)
#分类展示
barplot(GO, showCategory = 5,,split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
#展示特定通路
barplot(GO, showCategory = GO_select,,split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')
#2.KEGG分析####
KEGG <- enrichKEGG(gene = gene_entrez2$ENTREZID,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05, 
                   qvalueCutoff =0.05) 
KEGG_res <- KEGG@result
write.table(KEGG_res,file="KEGG_res.txt",sep="\t",quote=F,row.names = F)
#简单作图
barplot(KEGG, showCategory = 20)
dotplot(KEGG, showCategory = 20)
#3.弦图####

library(GOplot)
library(ggplot2)
library(tidyr)

DEG <- read.table("DEG2.txt",sep = "\t",check.names = F, stringsAsFactors = F,header = T,row.names = 1)
DEG <- mutate(DEG,Gene_symbol=rownames(DEG))

KEGG_res <- read.table("KEGG_res.txt", sep = "\t",header = TRUE, stringsAsFactors = FALSE)
KEGG <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") 
KEGG_res <- KEGG@result
kegg_select <- KEGG_res[KEGG_res$Description %in% c("AMPK signaling pathway","Oxidative phosphorylation"),]
kegg_list <- strsplit(kegg_select$geneID, "/")

## 构建 Term-Gene
n <- nrow(kegg_select)

a11 <- do.call(
  rbind,
  lapply(seq_len(n), function(i) {
    data.frame(
      Term  = kegg_select$Description[i],
      gene  = kegg_list[[i]],
      stringsAsFactors = FALSE
    )
  })
)

a11$value <- 1
a12 <- spread(a11, key = "Term", value = "value", fill = 0)
rownames(a12) <- a12$gene
a12 <- a12[, -1, drop = FALSE]

## 取交集并按 logFC 排序
id <- intersect(rownames(a12), DEG$Gene_symbol)
DEG2 <- DEG[id, , drop = FALSE]
DEG2 <- DEG2[order(DEG2$log2FoldChange, decreasing = TRUE), , drop = FALSE]

a12 <- a12[rownames(DEG2), , drop = FALSE]
a12$logFC <- DEG2$log2FoldChange

##绘制弦图
GOChord(
  a12,
  space        = 0.001,     # 基因之间的间距
  gene.order   = "logFC", # 按 logFC 对基因排序
  gene.space   = 0.25,    # 基因名跟圆圈的相对距离
  gene.size    = 3,       # 基因名字体大小
  border.size  = 0.1,     # 线条粗细
  process.label= 8        # term 字体大小
)

#通路图可视化
#BiocManager::install('pathview')
library(pathview)
DEG <- rownames_to_column(DEG,var = 'SYMBOL')
pgenes <- inner_join(DEG,gene_entrez2,by = 'SYMBOL')
pgenes <- dplyr::select(pgenes,ENTREZID,P.Value)

pvalue <- -log10(pgenes$P.Value)
gene_names <- as.character(pgenes$ENTREZID)

gene_data <- pvalue
names(gene_data) <- gene_names
#选择兴趣的通路
hsa04978 <- pathview(gene.data  = gene_data,
                     pathway.id = "hsa04978", #通路名字
                     species    = "hsa", #代表人类
                     out.suffix = "pathview",  
                     limit      = list(gene= 10, cpd= 1)) #对基因数量限制，cpd化合物数量


#九、GSEA富集分析####
# https://www.gsea-msigdb.org/gsea/msigdb/index.jsp
#1.常规GSEA####
rm(list = ls())
library(tidyverse)
library(limma)
setwd('E:/R_do/阿伦生信教学代码/GEO/差异分析/')#设置工作路径
load('GSE299988.rda')
rt <- rt3
exp <- exp3
str(exp) 
exp <- exp %>% mutate(across(everything(), as.numeric))
str(exp)
#确认行名列名一致#确认行名列名一致
identical(rownames(rt),colnames(exp))
#组间差异分析
group_list <- factor(rt$group,levels = c("Normal","Tumor"))
group_list
design <- model.matrix(~group_list)
#比较矩阵命名
design
#线性模型拟合
fit <- lmFit(exp,design)
#贝叶斯检验
fit2 <- eBayes(fit)
deg <- topTable(fit2, coef = 2, number = Inf)
DEG = na.omit(deg) 
DEG <- rownames_to_column(DEG,var = 'Gene')

#整理输入数据
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
#转换ENTREZID
# genelist <- bitr(DEG$Gene, fromType="SYMBOL",
#                  toType="ENTREZID", OrgDb='org.Hs.eg.db')
# DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))
geneList <- DEG[,2]
names(geneList) <- as.character(DEG[,'Gene_symbol'])
#排序
geneList <- sort(geneList, decreasing = TRUE)
geneList

#读取gmt文件
gmt <- read.gmt('h.all.v2026.1.Hs.symbols.gmt')
GSEA <- GSEA(geneList,
             TERM2GENE = gmt, 
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             minGSSize = 20,
             maxGSSize = 500) 

res <- GSEA@result
write.csv(res,file = 'H_GSEA.csv')

p <- gseaplot2(GSEA,geneSetID = 'HALLMARK_ALLOGRAFT_REJECTION',
               title = 'HALLMARK_ALLOGRAFT_REJECTION',color = 'red')
p
p1 <- gseaplot2(GSEA,geneSetID = 1:3,pvalue_table = T, subplots = 1:2)
p1  


#2.单基因GSEA####
rm(list = ls())
library(tidyverse)
setwd('E:/scleraBioinformaticsAnalysis/Datasets/3_Myopia_Tissues/High/GSE280609')#设置工作路径
load('GSE280609.rda')
rt <- rt3
exp <- exp3
str(exp) 
exp <- exp %>% mutate(across(everything(), as.numeric))
str(exp)
#确认行名列名一致#确认行名列名一致
identical(rownames(rt),colnames(exp))

#提取 自定义基因,整理高低表达组
rt1 <- as.data.frame(t(exp['PRKAA1',]))
rt1$group <- ifelse(rt1$PRKAA1 > median(rt1$PRKAA1),'up','down')
rt1 <- arrange(rt1,'PRKAA1',group)
identical(rownames(rt1),colnames(exp))
exp <- exp[,rownames(rt1)]
identical(rownames(rt1),colnames(exp))
# [1] TRUE
#组间差异分析
group_list <- factor(rt1$group,levels = c("down","up"))
group_list
design <- model.matrix(~group_list)
#比较矩阵命名
design
#线性模型拟合
fit <- lmFit(exp,design)
#贝叶斯检验
fit2 <- eBayes(fit)
deg <- topTable(fit2, coef = 2, number = Inf)
DEG = na.omit(deg) 
DEG <- rownames_to_column(DEG,var = 'Gene')
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
#转换ENTREZID
# genelist <- bitr(DEG$Gene, fromType="SYMBOL",
#                  toType="ENTREZID", OrgDb='org.Hs.eg.db')
# DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))
geneList <- DEG[,2]
names(geneList) <- as.character(DEG[,'Gene_symbol'])
#排序
geneList <- sort(geneList, decreasing = TRUE)
geneList
#读取gmt文件
gmt <- read.gmt('c2.cp.kegg_medicus.v2026.1.Hs.symbols.gmt')
GSEA <- GSEA(geneList,
             TERM2GENE = gmt,
             pvalueCutoff = 0.05, #改成1获得全部结果
             pAdjustMethod = "BH",
             minGSSize = 20,
             maxGSSize = 500) 
res <- GSEA@result
write.csv(res,file = 'KEGG_GSEA.csv')
p <- gseaplot2(GSEA,geneSetID = 'KEGG_MEDICUS_REFERENCE_COPI_VESICLE_FORMATION',
               title = 'KEGG_MEDICUS_REFERENCE_COPI_VESICLE_FORMATION',color = 'red')
p
p1 <- gseaplot2(GSEA,geneSetID = 1:3)
p1  

p2 <- gseaplot2(GSEA,
                geneSetID = c('KEGG_MEDICUS_REFERENCE_WNT_SIGNALING_PATHWAY',
                              'KEGG_MEDICUS_REFERENCE_ITGA_B_FAK_CDC42_SIGNALING_PATHWAY'))
p2

#3.是不是可以自创基因集？
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
gmt <- read.gmt('GOBP_MITOPHAGY.v2025.1.Hs.gmt')



#十、韦恩图VN_plot####
rm(list = ls())
library(tidyverse)
setwd('E:/R_do/阿伦生信教学代码/GEO/VN图/')#设置工作路径
#1.数据读取####
DEG <- read.table("DEG2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
table(DEG$change)
diff <- rownames(DEG)[DEG$change!='NOT']

dat1 <- read.csv('ferroptosis.csv')
dat1 <- dat1$Gene.Symbol

dat2 <- read.csv('mitophage.csv')
dat2 <- dat2$Gene.Symbol

a <- list(DEGs=diff,
          ferroptosis=dat1,
          mitophage=dat2)

#演示数据，sample随机抽样
library(RColorBrewer)
display.brewer.pal(9,'Set1')
brewer.pal(9, "Set1")[3]
"#4DAF4A"
display.brewer.pal(8,'Accent')

'#E41A1C'（红）  '#377EB8'（蓝）  '#4DAF4A'（绿）   '#984EA3'（紫）
'#FF7F00'（橙）  '#FFFF33'（黄）  '#A65628'（棕）   '#F781BF'（粉）
'#1B9E77'（深绿）'#D95F02'（深橙）'#7570B3'（深紫） '#E7298A'（深粉）
'#66A61E'（亮绿）'#E6AB02'（亮黄）'#A6761D'（棕黄） '#666666'（深灰）
'#A6CEE3'（浅蓝）'#1F78B4'（中蓝）'#B2DF8A'（浅绿） '#33A02C'（中绿）
'#FB9A99'（浅红）'#E31A1C'（中红）'#FDBF6F'（浅橙） '#FF7F00'（中橙）

set.seed(123) #随机种子
a <- list(a=sample(1:100, 35),
          b=sample(1:100, 45),
          c=sample(1:100, 55),
          d=sample(1:100, 80),
          e=sample(1:100, 90),
          f=sample(1:100, 100))
#2.VN作图####
#2.1 ggvenn包 最多4组
library(ggvenn)
ggvenn(a)


ggvenn(
  a,
  columns = NULL,
  show_elements = FALSE,
  show_percentage = TRUE,
  digits = 1,
  fill_color = c("blue", "yellow", "green", "red"),
  fill_alpha = 0.3,
  stroke_color = "black",
  stroke_alpha = 1,
  stroke_size = F,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 6,
  text_color = "black",
  text_size = 3,
  label_sep = ",",
  count_column = NULL,
  show_outside = c("auto", "none", "always"),
  auto_scale = FALSE
)
#2.2 venn包 最多7组  
# BiocManager::install('venn')
library(venn) 
venn(a)

{
venn(
  a,
  ilabels = "counts",  # 显示数字
  zcolor = c("#FF9999", "#66B2FF", "#99FF99"),  # 自定义柔和花瓣色（粉、蓝、绿）
  opacity = 0.6,       # 透明度
  box = FALSE,         # 是否除外边框
  ilcs = 1,            # 圈内数字大小
  sncs = 1.2,          # 组名大小
  border = "white",    # 边缘颜色
  lwd = 1.1            # 边框线条粗细
)
  }

#3.输出重叠结果####
res <- Reduce(intersect, a)
res
write.csv(res,file = 'venn.csv',row.names = F)

#十一、箱线图展示差异####
rm(list=ls())
setwd('E:/R_do/阿伦生信教学代码/GEO/VN图') #设置工作路径
library(tidyverse)
library(ggplot2)
library(ggpubr)
load('GSE299988.rda')

#数据准备
exp <- exp3
str(exp)
exp <- exp %>% mutate(across(everything(), as.numeric))
str(exp)
#分组信息
rt <- rt3

#差异表达可视化验证（箱线图）
exps <- exp %>% t() %>% as.data.frame() %>% dplyr::select("PPARGC1A",'SNCA')
identical(rownames(exps),rownames(rt))
# [1] TRUE 

rts <- cbind(exps,rt)

# 转换数据为长格式，方便绘图
rts_long <- rts %>% pivot_longer(cols = c("PPARGC1A",'SNCA'),
                                 names_to = "genes",values_to = "Expression")

#箱线图 gglot2
library(ggplot2)
p <- ggplot(data = rts_long,aes(x = genes,y = Expression,
                                fill = group))+
  geom_boxplot() +stat_compare_means(method ='wilcox.test') +theme_classic()+
  scale_fill_manual(values = c("#e67e22", "#3498db"))
p


#小提琴图 ggpubr
library(ggplot2)
library(ggpubr)
p <- ggviolin(rts_long,x = "genes",y = "Expression",
              fill = "group",
              palette = c("#e67e22", "#3498db"),
              add.params = list(fill = "white")) +
  stat_compare_means(aes(group = group),
                     method = "wilcox.test", #method='wilcox.test' ;method = "t.test"
                     label = "p.signif",symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                                           symbols = c("***","**","*","ns")))
p


'#E41A1C'（红）  '#377EB8'（蓝）  '#4DAF4A'（绿）   '#984EA3'（紫）
'#FF7F00'（橙）  '#FFFF33'（黄）  '#A65628'（棕）   '#F781BF'（粉）
'#1B9E77'（深绿）'#D95F02'（深橙）'#7570B3'（深紫） '#E7298A'（深粉）
'#66A61E'（亮绿）'#E6AB02'（亮黄）'#A6761D'（棕黄） '#666666'（深灰）
'#A6CEE3'（浅蓝）'#1F78B4'（中蓝）'#B2DF8A'（浅绿） '#33A02C'（中绿）
'#FB9A99'（浅红）'#E31A1C'（中红）'#FDBF6F'（浅橙） '#FF7F00'（中橙）

mycol <- c('#984EA3','#B2DF8A')

mycol <- c("#e67e22", "#3498db")

{
p1 <- ggviolin(rts_long, x = "genes", y = "Expression", fill = "group", 
               palette = mycol, legend = "top", alpha = 0.5,
               add = "boxplot", 
               add.params = list(width = 0.3, #箱图的宽度
                                 alpha = 0.8, outlier.shape = NA),  # 调整箱线图宽度和透明度
               outlier.shape = NA, color = alpha("black", 0.5)) +             
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***","**","*","ns"))) + 
  geom_jitter(aes(color = group),
              position = position_jitterdodge(
                jitter.width = 0.25,
                dodge.width = 0.8),   #散点的位置
              size = 1.2,
              alpha = 0.3) +
  scale_color_manual(values = mycol) +
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color = NA),
        panel.border = element_rect(fill = NA, color = NA),
        axis.text.x = element_text(size = 10, face = "plain", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, face = "plain", hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5, size = 10, face = "bold"))
}

p1


#十二、WGCNA####
rm(list = ls())
options(stringsAsFactors = FALSE)
library(tidyverse)
#install.packages("WGCNA")
library(WGCNA)
setwd('E:/R_do/阿伦生信教学代码/GEO/WGCNA/')#设置工作路径
load('GSE84844.rda')

#1.数据准备####

exp <- as.data.frame(t(exp))
#输入的基因:可以是全部基因、或者高变基因、也可以是差异基因，建议至少至少2000以上
#样本：15对15以上

#1.1 选取高变基因，var方差前百分之
vars_res <- apply(exp,2,var)
# 计算百分位数截止值
per_res <- quantile(vars_res, probs = seq(0, 1, 0.25)) # quantile生成分位数
# 选取方差位于前25%的基因
per_gene <- exp[, which(vars_res > per_res[4])]
datExpr0 <- data.matrix(per_gene)

#1.2 选取高变基因，绝对偏差中位数排名前5000的基因
# mad_res <- apply(exp, 2, mad, na.rm = TRUE)
# mad_gene <- exp[, order(mad_res, decreasing = TRUE)[1:5000], drop = FALSE]
# datExpr0 <- data.matrix(mad_gene)

#2.数据检查####
#检查样本和基因是否ok
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

#如果没有达标就需要筛选，[1] TRUE 也能运行，不影响结果
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#进行样本聚类
sampleTree = hclust(dist(datExpr0), method = "average")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
#剪切线
abline(h = 70, col = "red")

#切除未聚类的样本
clust = cutreeStatic(sampleTree, cutHeight = 70, minSize = 10)#注意：cutHeight = 这里才是真正切割
table(clust)
keepSamples = (clust==1)
#注意查看切割数量0有多少 1有多少
datExpr0 = datExpr0[keepSamples, ]


#表型数据分组
table(rt)
traitData <- data.frame(control=c(rep(1,30),rep(0,30)),
                        pSS=c(rep(0,30),rep(1,30)))
row.names(traitData) <- rownames(exp)

sameSample <- intersect(rownames(datExpr0), rownames(traitData))
datExpr <- datExpr0[sameSample,]
datTraits <- traitData[sameSample,]

#去除后再聚类
sampleTree2 = hclust(dist(datExpr), method = "average")
plot(sampleTree2)
#标记颜色
traitColors = numbers2colors(datTraits, signed = FALSE)
#3.第一部分样本可视化作图####
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
save(datExpr, datTraits, file = "WGCNA_data1.rda")

#4.网络构建####
rm(list = ls())
load(file = "WGCNA_data1.rda")
#power值散点图
allowWGCNAThreads()   #多线程工作
powers = c(c(1:10), seq(from = 12, to=20, by=2)) #幂指数范围1到20
powers
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft
#作图
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
#拟合指数与power值散点图，确保散点在0.8以上，最好0.9以上
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") #根据自己选择修改，如0.9则不变
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#自动计算softPower
softPower <- sft$powerEstimate
softPower
#查看最佳softPower值
#选择power
softPower = 6
adjacency = adjacency(datExpr, power = softPower)

#TOM矩阵
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
#基因聚类
geneTree = hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
##动态剪切模块识别
minModuleSize = 30  #最小模块基因数
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
#标记颜色
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
#5.第二部分基因聚类作图####
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#6.相似模块聚类####
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#剪切高度，0.25标识相异性（代表75%以上相似性模块合并）及MEDissThres越小代表要求相似性越高，最终模块越多
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
#相似模块合并（也可不做，根据结果调整）
merge = mergeCloseModules(datExpr, dynamicColors,cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
plotDendroAndColors(geneTree, mergedColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
save(MEs, moduleLabels, moduleColors, geneTree, file = "WGCNA_data2.rda")

#7.第三部分模块与特征热图####
#模块与表型数据热图
rm(list=ls())
load('WGCNA_data1.rda')
load('WGCNA_data2.rda')
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
#热图展示
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
sizeGrWindow(10,6)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50), #greenWhiteRed
               textMatrix = textMatrix,
               cex.lab.y = 0.5,
               cex.lab.x = 0.5,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#模块选择相关性0.4以上，P<0.05
#结果不好，调节高变基因，如25%至50%

#8.模块提取####
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
#完成基因模热图
#进行模块基因输出
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

#9.批量输出散点图####
for (trait in traitNames){
  traitColumn=match(trait,traitNames)  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      outPdf=paste(trait, "_", module,".pdf",sep="")
      pdf(file=outPdf,width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      abline(v=0.8,h=0.5,col="red")
      dev.off()
    }
  }
}

#10.批量输出模块基因####
for (mod in 1:nrow(table(moduleColors)))
{  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0(modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}

#11.其他补充图####
rm(list=ls())
load('WGCNA_data1.rda')
load('WGCNA_data2.rda')
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#11.1 TOM聚类图####
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);#power与之前保持一致
#增加强度
plotTOM = dissTOM^7;
#对角线设置为NA
diag(plotTOM) = NA;
#绘图（十分缓慢，建议进行nSelect挑选）
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

#挑选基因，快速作图
nSelect = 400
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
sizeGrWindow(9,9)
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

#11.2 特征模块热图####
rm(list=ls())
load('WGCNA_data1.rda')
load('WGCNA_data2.rda')
#重新计算模块特征基因
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
#提取特征
Traits = as.data.frame(datTraits$pSS); #选择特征
names(Traits) = "Traits"
#将特征添加到现有模块特征基因中
MET = orderMEs(cbind(MEs, Traits))
#绘制特征基因与性状之间的关系图
par(cex = 1)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

#绘制树状图
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
#绘制热图矩阵
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
#WGCNA结束

#十三、Lasso降维####
rm(list=ls())
library(tidyverse)
setwd('E:/R_do/阿伦生信教学代码/GEO/Lasso/')#设置工作路径
load('GSE84844.rda')

dif <- read.table('DEG2.txt',header = TRUE, sep = "\t",row.names = 1)
difs <- rownames(dif)[dif$change!='NOT']
#自己的数据重新命名dat和rt，方便后续代码一致
dat <- t(exp[difs,])
rt <- rt

rt$group <- factor(rt$group, levels = c("control", "pSS"))
identical(rownames(rt),rownames(dat))

# install.packages("glmnet")
library(glmnet)
#设置随机种子
set.seed(123)
#构建Lasso模型x，y
x <- as.matrix(dat)
y <- rt$group
fit <- glmnet(x,y,family='binomial', alpha=1)

pdf(file="lambda.pdf", height= 4.5, width= 4.5)
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit = cv.glmnet(x, y, alpha = 1, family="binomial",nfolds = 10)

pdf(file="deviance.pdf", height= 4.5, width= 4.5)
plot(cvfit)
dev.off()
#输出lasso结果
coef <- coef(fit,s=cvfit$lambda.min)
index <- which(coef!=0)
actCoef <- coef[index]
lassogenes <- row.names(coef)[index]
lassogenes
geneCoef <- cbind(Gene=lassogenes,Coef=actCoef)
geneCoef <- as.data.frame(geneCoef[-1,])

#coef柱状图
order <- arrange(geneCoef,Coef)
order$Coef <- as.numeric(order$Coef)
par(mar = c(5, 8, 4, 2)) #下左上右
barplot(order$Coef,names.arg = order$Gene,
        horiz=TRUE,las=1,
        col=ifelse(order$Coef > 0,"tomato","steelblue"),
        xlab="Coefficient Value",
        main="Biomarker Impact Direction",
        cex.names = 0.9,cex.lab = 1,cex.main = 1)

#十四、随机森林####
rm(list=ls())
library(tidyverse)
setwd('E:/R_do/阿伦生信教学代码/GEO/RF/')#设置工作路径
load('GSE84844.rda')

dif <- read.table('DEG2.txt',header = TRUE, sep = "\t",row.names = 1)
difs <- rownames(dif)[dif$change!='NOT']
#自己的数据重新命名dat和rt，方便后续代码一致
dat <- t(exp[difs,])%>% as.data.frame()
rt <- rt

rt$group <- factor(rt$group, levels = c("control", "pSS"))
# rt$group <- c(rep(0, 30), rep(1, 30)) #也可以改成0和1
identical(rownames(rt),rownames(dat))

#1.RF建模####
#install.packages('randomForest')
library(randomForest)
#设置随机种子
set.seed(123)
fit <- randomForest(x = dat, #表达数据
                    y = as.factor(rt$group),#分组数据
                    importance = TRUE,ntree = 500,proximity=TRUE)

plot(fit) 

#2.基因气泡图可视化####
pdf(file="RF_top10.pdf", height= 4.5, width= 5)

par(mar=c(5,2,2,2)) 
varImpPlot(fit, n.var=10, scale=FALSE,type = 2,main="Importance of Variables for top 10 gens",cex =0.9)

dev.off()
#3.数据提取####
res_rf <- as.data.frame(importance(fit, scale=FALSE))
top10_genes <- res_rf %>% slice_max(MeanDecreaseGini,n=10) %>% rownames()
top10_genes

top10_genes <- data.frame(gene=top10_genes)
write.csv(top10_genes, file= "RF.csv")

#十五、GSVA分析####
rm(list=ls())
library(tidyverse)
library(clusterProfiler)
setwd('E:/R_do/阿伦生信教学代码/GEO/GSVA/')#设置工作路径
load('GSE84844.rda')
gmt <- read.gmt('h.all.v2025.1.Hs.symbols.gmt') #读取gmt数据
geneset <- split(gmt$gene,gmt$term) #切割为list
#1.GSVA分析####
# BiocManager::install("GSVA")
library(GSVA)
gsvaParam <- gsvaParam(as.matrix(exp),geneset,kcdf = "Gaussian",#"Gaussian", "Poisson"
                  minSize = 3) 
gsva <- gsva(gsvaParam)
#2.1 GSVA热图####
library(pheatmap)
# rownames(gsva) <- substr(rownames(gsva),10,nchar(rownames(gsva))) #提取geneset第10个字符以后的字符
identical(rownames(rt),colnames(gsva))

pheatmap(gsva,
         annotation_col=rt,
         color = colorRampPalette((c("blue","white","red")))(100),
         cluster_cols = F,
         cluster_rows = T,
         show_colnames = F,
         show_rownames = T ,
         border_color = 'black',#不加边框 NA
         scale='row', #scale='row' 进行标准化
         fontsize = 8,
         fontsize_row = 6,
         fontsize_col = 8) 


#2.2 GSVA的通路差异分析####
library(limma)
group_list <- factor(rt$group,levels = c("control","pSS"))#分组
group_list
design <- model.matrix(~group_list)
#比较矩阵命名
design
#线性模型拟合
fit <- lmFit(gsva,design)
#贝叶斯检验
fit2 <- eBayes(fit)
deg <- topTable(fit2, coef = 2, number = Inf)
DEG <- na.omit(deg) #differently expressed genes
write.table(DEG, file='GSVA_DEG.txt',sep = "\t",row.names = T,col.names = NA,quote = F)

#筛选adj.P.Val < 0.05的通路
# DEG <- filter(DEG,adj.P.Val < 0.05)
#作图
DEG$geneset <- rownames(DEG)
DEG <- arrange(DEG,desc(t))
cutoff <- 1.96 #自定义t 值分界
DEG$change <- ifelse(DEG$t > cutoff, "UP", ifelse(DEG$t < -cutoff, "DOWN", "NOT")) 

ggplot(DEG, aes(x = t, y = reorder(geneset, t), fill = change)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("UP" = "blue", "DOWN" = "green", "NOT" = "gray")) +
  xlab("t value of GSVA score") +
  ylab("") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 5))

#美化版
# DEG$geneset <- substr(DEG$geneset,10,nchar(DEG$geneset)) #提取geneset第10个字符以后的字符
p <- {
ggplot(DEG, aes(x = t, y = reorder(geneset, t), fill = change)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = c(-cutoff,cutoff), #添加两侧竖线
             color = 'white',
             linetype = 2, #虚线
             size = 0.5)+ #大小
  geom_text(
    data = DEG,  # 添加通路名
    aes(
      x = ifelse(t > 0, -0.1, 0.1),  
      y = reorder(geneset, t),       
      label = geneset,               # 标签内容：基因集名称
      hjust = ifelse(t > 0, 1, 0)),    
    size = 2.5,  # 标签文字大小
    color = "black",  # 标签颜色
    vjust = 0.5)+       # 垂直居中对齐
  scale_fill_manual(values = c("UP" = "blue", "DOWN" = "green", "NOT" = "gray")) + #定义颜色
  xlab("t value of GSVA score") + #x轴名
  ylab("") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),  # 隐藏原y轴标签
    axis.ticks.y = element_blank(), # 隐藏原y轴刻度
    axis.line.y = element_blank(),  # 隐藏原y轴线
    axis.line.x = element_line(size = 0.4)) +  # x轴线粗细
  xlim(min(DEG$t) - 0.5, max(DEG$t) + 0.5) # 定义x轴边距
}
p

#十六、ssGSEA免疫浸润####
rm(list=ls())
library(tidyverse)
library(clusterProfiler)
setwd('E:/R_do/阿伦生信教学代码/GEO/GSVA/')#设置工作路径
# install.packages("readxl")
data <- readxl::read_xlsx('imm.xlsx',col_names = T)
cellMarker <- split(data$Metagene,data$`Cell type`)
save(cellMarker,file = '28cellMarker.Rdata')

rm(list=ls())
load('GSE84844.rda')
load('28cellMarker.Rdata')
#1.ssGSEA分析####
# BiocManager::install("GSVA")
library(GSVA)
ssgsea_res <- gsva(
  expr = as.matrix(exp),    # 表达矩阵（行=基因，列=样本）
  cellMarker,               # 基因集列表（cellMarker）
  method = "ssgsea",        # 指定 ssGSEA 算法
  kcdf = "Gaussian")        # 连续数据使用Gaussian；count数据使用Poisson
ssgsea_res <- as.data.frame(ssgsea_res)
#新版GSVA运行下面
ssgseaParam <- ssgseaParam(as.matrix(exp),cellMarker) 
ssgsea_res <- gsva(ssgseaParam)
ssgsea_res <- as.data.frame(ssgsea_res)
#2.1免疫浸润箱线图####
ssgsea_t <- t(ssgsea_res)
identical(rownames(ssgsea_t),rownames(rt))
rts <- cbind(rt,ssgsea_t)

rts_long <- rts %>% pivot_longer(cols = colnames(rts)[-1],
                                 names_to = "cellMarker",values_to = "Score")
# 绘制小提琴图并添加 Wilcoxon 检验的 p 值（用星号表示）
library(ggplot2)
library(ggpubr)
p <- ggboxplot(data = rts_long, x="cellMarker", y="Score", fill = "group",
            ylab="ssGSEA Score",
            xlab="",
            legend.title = 'group',
            palette = c("#E74C3C", "#3498DB"),
            width=0.6, add = "none",
            outlier.shape = NA #不显示异常值
            )+
            rotate_x_text(45)+stat_compare_means(aes(group = group),method="wilcox",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
                               label = "p.signif")+ theme(
            axis.text.x = element_text(size = 10, # x轴文字大小
            color = "black",    # 文字颜色
            vjust = 1           # x轴文字垂直对齐
            ))
p
#2.1免疫浸润相关性展示####

#剔除对照,保留疾病:根据自己的分组进行修改
table(rt$group)  
ssgsea_res <- ssgsea_res[,rownames(rt)[rt$group=='pSS']]
exp <- exp[,rownames(rt)[rt$group=='pSS']]
#行列转置
ssgsea_t <- t(ssgsea_res)
exp_t <- t(exp)
identical(rownames(ssgsea_t),rownames(exp_t))
#加载R包
library(psych)
library(corrplot)
#2.1.1免疫细胞间相关性####
immu_cor <- corr.test(ssgsea_t,method ="spearman")
p <- immu_cor$p
r <- immu_cor$r

corrplot(r,tl.cex = 0.5)
#美化版
{
corrplot(r, method ="pie",
         type = "upper",#上半部分显示
         tl.cex =0.5,#标签字号
         tl.col ="black",#标签颜色
         tl.pos ="lt",#文字标签在左边
         tl.srt =90,#调整上方标签角度
         p.mat = p,
         sig.level = c(0.001,0.01,0.05),
         insig ="label_sig",
         pch.cex =1) #叉号大小
corrplot(r, method ="number",
         type = "lower",#下半部分显示
         tl.pos ="n", #取消文字标签显示
         cl.pos ="n",#取消颜色图例显示
         number.cex = 0.4,#数字大小
         add = TRUE #增加
) 
}
#2.1.2基因与免疫细胞的相关性####
#选取相关基因
genes <- c('PTPRC','IL7','CD3E','A1BG','A2M','RPS5','CD40')#自定义
exp_t <- exp_t[,genes]
#相关性分析
cor <- corr.test(exp_t,ssgsea_t,method ="spearman")
r<- cor$r #提取相关性系数
p <- cor$p #提取相关性系数的p值

#定义p值
sig <- function(p_val) {
  dplyr::case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01 ~ "**",
    p_val < 0.05 ~ "*",
    TRUE ~ "ns"
  )}

sig_p <- apply(p, c(1, 2), sig)

display_matrix <- matrix(
  paste0(round(r, 2), "\n(", sig_p, ")"),  # 内容
  nrow = nrow(r),                          # 行数
  ncol = ncol(r)                          # 列数
  )
#热图展示
library(pheatmap)
library(RColorBrewer)
pheatmap(
  r,
  scale = "none",  # 标准化
  color = colorRampPalette((c("blue","white","red")))(100),  # 红蓝渐变（红=正相关，蓝=负相关）
  breaks = seq(-1, 1, length.out = 101),
  show_rownames = TRUE,  # 显示基因名
  show_colnames = TRUE,  # 显示免疫细胞名
  display_numbers = display_matrix,
  treeheight_row = 0,    # 隐藏基因聚类树（4个基因无需聚类）
  treeheight_col = 10,   # 保留免疫细胞聚类树
  cellwidth = 20,        # 调整单元格宽度
  cellheight = 25,       # 调整单元格高度
  main = "Correlation between Genes and Immune Cell ssGSEA Scores"# 标题
  )

#3.单独相关性散点图####
exps <- cbind(exp_t,ssgsea_t)
colnames(exps)
#输入感兴趣的基因与通路
target_gene <- 'IL7'
target_cell <- "Activated CD8 T cell"

p <- ggplot(data = exps, aes(x = .data[[target_gene]], y = .data[[target_cell]])) +
              geom_point(alpha = 0.7, color = "black", fill = "#999999", size = 2) +
  geom_smooth(method = "lm",      # 线性拟合（若想非线性用method="loess"）
             se = TRUE,          # 显示95%置信区间（更严谨，可改为FALSE隐藏）
             color = "#2E86AB",  # 拟合线颜色
             linewidth = 1,      # 线宽
             fill = "#2E86AB",   # 置信区间填充色
             alpha = 0.2         # 置信区间透明度（不遮挡散点）
             ) +
  stat_cor(method = "spearman",        # 与corr.test方法一致
           label.x.npc = 0.55,         # 标注相关性x位置（0=左，1=右）
           label.y.npc = 0.05,         # 标注相关性Y位置（0=下，1=上）
           size = 4,                   # 文字大小（比之前更大，更易读）
           color = "#A23B72",          # 文字颜色（醒目）
           fontface = "bold",          # 文字加粗
           r.accuracy = 0.01,          # 相关系数保留2位小数
           p.accuracy = 0.001,         # P值保留3位小数（<0.001会显示p<0.001）
           label.sep = ","        # 自定义分隔符（更简洁）
            ) +
           labs( x = paste0(target_gene, " Expression"),  # X轴：基因表达量（英文更通用）
                y = paste0(target_cell, " SSGSEA Score"),    # Y轴：免疫细胞评分（明确指标类型）
                title = paste0(target_gene, " vs. ", target_cell),  # 标题简洁明了
                subtitle = "Spearman Correlation Analysis"     # 副标题：说明分析类型
           ) +
  theme_bw() +  # 白色背景+黑色边框（比theme_minimal更适合科研论文）
  theme(
    plot.title = element_text(
      hjust = 0.5, size = 15, face = "bold", color = "#2E86AB"
    ),  # 标题居中、加粗、配色呼应
    plot.subtitle = element_text(
      hjust = 0.5, size = 12, color = "#666666"
    ),  # 副标题居中、灰色
    axis.title = element_text(
      size = 13, face = "bold", color = "#333333"
    ),  # 坐标轴标题加粗
    axis.text = element_text(
      size = 11, color = "#555555"
    ),  # 坐标轴刻度文字清晰
    panel.grid = element_blank(),  # 隐藏网格线（减少干扰）
    legend.position = "none",       # 无图例（无需图例，直接标注）
    plot.margin = margin(t = 5, r = 5,b = 5, l = 5, unit = "mm")
    )
# 显示图片
print(p)


#十七、多数据集合并去批次####
rm(list=ls())
setwd('E:/R_do/阿伦生信教学代码/GEO/合并去批次/') #设置工作路径
library(tidyverse)
#1.数据读取####
load( "GSE55235.rda")
load( "GSE55457.rda")
load( "GSE77298.rda")
load( "GSE82107.rda")
#定义批次
GSE55235_rt$batch <- rep('GSE55235',30)
GSE55457_rt$batch <- rep('GSE55457',33)
GSE77298_rt$batch <- rep('GSE77298',23)
GSE82107_rt$batch <- rep('GSE82107',17)
#合并
rt <- rbind(GSE55235_rt,GSE55457_rt,GSE77298_rt,GSE82107_rt)
table(rt)
# batch
# group     GSE55235 GSE55457 GSE77298 GSE82107
# control       10       10        7        7
# OA            20       23       16       10
library(limma)
#如果已经normalizeBetweenArrays不需要运行
GSE55235 <- normalizeBetweenArrays(GSE55235)
boxplot(GSE55235, outline=F,las=2)

GSE55457 <- normalizeBetweenArrays(GSE55457)
boxplot(GSE55457, outline=F,las=2)

GSE77298 <- normalizeBetweenArrays(GSE77298)
boxplot(GSE77298, outline=F,las=2)

GSE82107 <- normalizeBetweenArrays(GSE82107)
boxplot(GSE82107, outline=F,las=2)
#2.取基因集共有的基因####
exp1 <- GSE55235 %>% as.data.frame %>% rownames_to_column("gene") 
exp2 <- GSE55457 %>% as.data.frame %>%  rownames_to_column("gene") 
exp3 <- GSE77298 %>% as.data.frame %>%  rownames_to_column("gene") 
exp4 <- GSE82107 %>% as.data.frame %>%  rownames_to_column("gene")
co_genes <- Reduce(intersect,list(exp1$gene,
                                  exp2$gene,
                                  exp3$gene,
                                  exp4$gene))
#3.merge合并####
exp <- merge(exp1,exp2)
exp <- merge(exp,exp3)
exp <- merge(exp,exp4)
exp <- column_to_rownames(exp,var = 'gene')
identical(colnames(exp),rownames(rt))
#合并后未去批次的箱线图
boxplot(exp, outline=F, las=2, 
        col = factor(rt$batch))
legend("topright", legend=unique(rt$batch), 
       fill=unique(factor(rt$batch)), cex=0.5)

#4.combat合并去批次####
identical(colnames(exp),rownames(rt))
# BiocManager::install('sva')
library(sva)
#定义生物学分组和批次组别
group <- rt$group
batch <- rt$batch
mod <- model.matrix(~group,data = rt)
combine <- ComBat(exp,batch=batch,mod=mod,par.prior=T)
boxplot(combine, outline=F, las=2, 
        col = factor(rt$batch))
legend("topright", legend=unique(rt$batch), 
       fill=unique(factor(rt$batch)), cex=0.5)
#合并后去批次的箱线图
combine <- normalizeBetweenArrays(combine)
boxplot(combine, outline=F, las=2, 
        col = factor(rt$batch))
legend("topright", legend=unique(rt$batch), 
       fill=unique(factor(rt$batch)), cex=0.5)

#5.PCA分析####
#PCA前
library(FactoMineR)
library(factoextra)
exp_t <- t(exp)#转置 行和列
before_pca <- PCA(exp_t,graph = F)
p1 <- fviz_pca_ind(before_pca,geom.ind = 'point',
                   col.ind = batch,palette = 'jco',
                   addEllipses = T, #添加圆的置信区间
                   legend.title='Batchs',title='Before_PCA')
p1

#PCA后
combine_t <- t(combine)#转置 行和列
after_pca <- PCA(combine_t,graph = F)
p2<- fviz_pca_ind(after_pca,geom.ind = 'point', col.ind = batch,palette = 'jco',#col.ind = group看看
                  addEllipses = T, #添加圆的置信区间
                  legend.title='Batchs',title='After_PCA')
p2

#差异分析
library(limma)
rts <- arrange(rt,desc(group))
exps <- combine[,rownames(rts)]
identical(colnames(exps),rownames(rts))
group_list <- factor(rts$group,levels = c("control","OA"))
group_list
design <- model.matrix(~group_list)
#比较矩阵命名
design
#线性模型拟合
fit <- lmFit(exps,design)
#贝叶斯检验
fit2 <- eBayes(fit)
deg <- topTable(fit2, coef = 2, number = Inf)
DEG <- na.omit(deg) 
#进行注释 fold change
#通常logFC设置1,adj.P.Val < 0.05
logFC_cutoff <- 1
type1 = (DEG$adj.P.Val < 0.05)&(DEG$logFC < -logFC_cutoff)
type2 = (DEG$adj.P.Val < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)

# DOWN   NOT    UP 
# 85 12972   180

DEG$genes <- rownames(DEG)

#十八、ROC曲线####

#单个变量的ROC####
rm(list=ls())
library(tidyverse)
setwd('E:/R_do/阿伦生信教学代码/GEO/ROC/')#设置工作路径
load('GSE84844.rda')
library(pROC) 

# gene <- c('PTPRC')
# expr <- t(exp[gene,])
# rts <- cbind(rt,expr)

#组别定义为0和1
gene <- c('PTPRC')
expr <- t(exp[gene,])

table(rt$group)
rt$status <- c(rep(0,30),rep(1,30))
rts <- cbind(rt,expr)
rts <- rts[,-1]

colnames(rts)

roc_obj <- roc(status~PTPRC, data = rts)

plot.roc(roc_obj,print.auc=TRUE)

plot.roc(roc_obj,
     print.auc=TRUE,   #输出AUC值
     auc.polygon = TRUE, #绘制AUC阴影区
     auc.polygon.col = adjustcolor("#2E86AB", alpha.f = 0.1), # 绘制AUC阴影区颜色
     main = "ROC CURVE",   #图形标题
     col= "#2E86AB",lwd = 2,    #曲线颜色
     print.thres = F,print.thres.col="red",print.thres.cex=0.5,#是否显示阈值：cut-off值字体的颜色、大小
     identity.col="#2E86AB",identity.lty = 2,identity.lwd = 1, # 对角线颜色、线及宽度
     legacy.axes=TRUE,asp = NA)

auc_value <- auc(roc_obj)
auc_value
ci_auc <- ci.auc(roc_obj)
ci_auc
text(x = 0.2, y = 0.1,labels = paste0("95% CI: ", round(ci_auc[1], 3), " - ", round(ci_auc[3], 3)),
     col = "red", cex = 0.8)

#多变量的ROC####
rm(list=ls())
library(tidyverse)
setwd('E:/R_do/阿伦生信教学代码/GEO/ROC/')#设置工作路径
load('GSE84844.rda')
library(pROC) 
genes <- c('PTPRC','A1CF','OSTC','OSTCP1')
expr <- t(exp[genes,])

table(rt$group)
rt$status <- c(rep(0,30),rep(1,30))
rts <- cbind(rt,expr)
rts <- rts[,-1]
colnames(rts)

formula <- as.formula(paste("status~", paste(genes, collapse = "+")))
formula
roc_list <- roc(formula, data = rts)
ggroc(roc_list)

#美化+AUC
auc_values <- sapply(roc_list, function(x) round(auc(x), 3)) 
auc_values
auc_df <- data.frame(genes = names(auc_values),  # 基因名
          auc = auc_values, # AUC值
          x = 0.85,   # 所有AUC的X坐标（可调整，如0.6/0.8）
          y = seq(0.4, 0.1, by = -0.1)  # Y坐标错开（0.4→0.3→0.2→0.1），避免标注重叠
)

ggroc(roc_list, legacy.axes = TRUE, size = 0.8) + #线条粗细
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") + #对角线
  scale_color_brewer(palette = "Set1",
                     name = "Genes",          # 图例标题
                     labels = genes           # 图例标签（纯基因名，无多余字符）
  ) +
  labs(title = "ROC Curves for Multiple Genes",
       x = "1 - Specificity",
       y = "Sensitivity") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), #主题大小位置
        legend.position = "right")+ #图例位置
  geom_text(data = auc_df,
    aes(x = x, y = y, label = paste0(genes,"：AUC = ", auc), color = genes ), #设置位置，字符，颜色
    size = 3,        # 字体大小
    fontface = "bold",# 加粗
    inherit.aes = FALSE,show.legend = FALSE #不投射图例以及轴
  ) 

#十九、SVM-RFE算法####
rm(list=ls())
setwd('E:/R_do/阿伦生信教学代码/GEO/RF/')#设置工作路径
load('GSE84844.rda')
#加载R包
library(tidyverse)
library(e1071)
#1.数据整理####
exp <- exp  #重新命名 行是基因，列是样本的基因表达数据
rt <- rt    #重新命名分组数据
#提取差异基因
dif <- read.table('DEG2.txt',header = TRUE, sep = "\t",row.names = 1)
difs <- rownames(dif)[dif$change!='NOT']
exp <- exp[difs,]

#合并基因表达和分组数据
input <- cbind(t(exp),rt) %>% dplyr::select(group,everything())
str(input) #查看分组是否为因子
input$group <- as.factor(input$group)
input$group #确认

#加载自定义SVM-RFE函数
source("msvmRFE.R")

#2.SVM-RFE特征选择####
#执行SVM-RFE,设置k值，halve.above值，用于切割分配数据
svmRFE(input, k = 10, halve.above = 50) 
#设置折交叉验证，可以是5，也可以10，都行，折越大，运行越慢
n_folds <- 10
num_samples <- nrow(input) 
#交叉验证折叠
set.seed(123) #设置种子
fold_index <- rep(1:n_folds, length.out = num_samples)  
fold_index <- sample(fold_index)  
fold_list <- lapply(1:n_folds, function(x) which(fold_index == x))  
# 执行交叉验证
set.seed(123) #设置种子 
cv_results <- lapply(fold_list,svmRFE.wrap,input,k = 10,halve.above = 50)

#特征排序结果输出保存
feature_ranking <- WriteFeatures(cv_results, input, save = FALSE)
# 保存特征排序结果
write.table(feature_ranking, file = "SVM_feature_ranking.txt", sep = "\t",
            quote = FALSE, row.names = T)

#3.交叉验证性能评估####
#评估特征数量范围
max_features <- ifelse(ncol(input) - 1 > 30, 30, ncol(input) - 1)
max_features

# max_features <- ncol(input) - 1
# max_features
#每完成10%输出一次进度
set.seed(123) #设置种子
performance_list <- lapply(1:max_features, function(features) {
  if(features %% round(max_features/10) == 0) {
    cat(sprintf(">> 进度: %.0f%% [%s]\n", 
                100*features/max_features, 
                Sys.time()))
  }
  FeatSweep.wrap(features, cv_results, input)
})

#4.性能可视化####
#提取误差率
error_rates <- sapply(performance_list, function(x) {
  if (is.null(x)) return(NA) else x$error
})

baseline_error <- min(prop.table(table(input$group)))
PlotErrors(error_rates, no.info = baseline_error)

pdf("SVM-RFE错误率.pdf", width = 7, height = 5)
PlotErrors(error_rates, no.info = baseline_error)
dev.off()

which.min(error_rates)
#计算准确率
baseline_accuracy <- 1 - baseline_error
accuracy_rates <- 1 - error_rates
Plotaccuracy(accuracy_rates, no.info = baseline_accuracy)

pdf("SVM-RFE准确率.pdf", width = 7, height = 5)
Plotaccuracy(accuracy_rates, no.info = baseline_accuracy)
dev.off()

which.max(accuracy_rates)

#5.提取特征基因####
SVM_Genes <- feature_ranking[1:which.min(error_rates),'FeatureName']
SVM_Genes

#保存
write.table(SVM_Genes, file = "SVM-RFE特征基因.txt", sep = "\t", 
            quote = FALSE, row.names = T, col.names = FALSE)



#其他图
#特征数量与性能关系图
pdf("SVM-RFE特征数量与性能关系图.pdf", width = 8, height = 6)
plot(1:max_features, accuracy_rates, type = "b", col = "skyblue", lwd = 2,
     xlab = "Number of Features", ylab = "Performance Metric", 
     main = "SVM-RFE Feature Optimization",
     ylim = c(0, 1))
lines(1:max_features, error_rates, type = "b", col = "firebrick", lwd = 2)
abline(h = baseline_accuracy, col = "gray50", lty = 2)
abline(h = baseline_error, col = "gray50", lty = 2)
legend("right", legend = c("Accuracy", "Error Rate", "Random Baseline"),
       col = c("skyblue", "firebrick", "gray50"), lwd = 2, lty = c(1, 1, 2),
       cex = 0.8,bty = "n")
dev.off()
#条形图
features <- head(feature_ranking,  which.min(error_rates))
pdf("SVM-RFE特征条形图.pdf", width = 8, height = 6)
ggplot(features, aes(x = reorder(FeatureName, AvgRank), 
                         y = AvgRank, fill = AvgRank)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "gray50", high = "firebrick") +
  coord_flip() +
  labs(title = "Top Features Selected by SVM-RFE", 
       x = "Gene Name", y = "Average Rank (Lower = More Important)") +
  theme_minimal(base_size = 12) +theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()


#二十、IOBR免疫微环境####
rm(list=ls())
library(tidyverse)
#用devtools安装IOBR
# https://github.com/IOBR/IOBR 点击Releases手动安装tar.gz或代码安装
# https://iobr.github.io/book/ 说明书
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

depens<-c('tibble', 'survival', 'survminer', 'limma', "DESeq2","devtools", 'limSolve', 'GSVA', 'e1071', 'preprocessCore', 
          "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor",  "timeROC", "pracma", "factoextra", 
          "FactoMineR", "WGCNA", "patchwork", 'ggplot2', "biomaRt", 'ggpubr', 'ComplexHeatmap')
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
}
if (!requireNamespace("IOBR", quietly = TRUE))
  devtools::install_github("IOBR/IOBR") #包太大 容易报错 手动安装

library(IOBR) #加载

#1.基因集评分####
#特征评分选择.
signature_score_calculation_methods
names(signature_tme) #免疫微环境相关基因集
names(signature_metabolism)#代谢相关基因集
names(signature_tumor)#生物医学基础相关基因集

names(signature_collection) #全部基因集
signature <- signature_collection #查看全部基因集
signature$APM #基因提取

signature_collection_citation #基因集来源文献
citation <- signature_collection_citation #基因集来源文献
view(citation) #查看基因集来源文献

names(go_bp) #GO
names(kegg) #KEGG
names(hallmark) #hallmark
#数据举例
setwd('E:/R_do/阿伦生信教学代码/GEO/GSVA/')#设置工作路径
load('GSE84844.rda')

sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = exp,
                             signature       = signature_collection,
                             method          = "pca",
                             mini_gene_count = 2)

sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = exp,
                             signature       = go_bp, 
                             method          = "ssgsea",
                             mini_gene_count = 2)



#2.免疫微环境—反卷积####
#tme_deconvolution_methods 查看算法
tme_deconvolution_methods

#算法介绍具体查阅文献

help(deconvo_tme)

#cibersort
cibersort<-deconvo_tme(eset = exp, method = "cibersort", #算法名
                       arrays = TRUE,  #是否是芯片数据,Currently affects 'CIBERSORT', 'svr' and 'xCell'
                       tumor = F, #是否是肿瘤, Currently affects 'EPIC'
                       perm = 10 )  #perm理论越大越精准但是耗时长
#保存
write.table(cibersort,file = 'cibersort.txt',sep = '\t',row.names = T,col.names = NA,quote = F)

#绘制比例图
res<-cell_bar_plot(input = cibersort, features = colnames(cibersort)[2:23], title = "CIBERSORT Cell Fraction",
                   coord_filp =T  #coord_filp =T横向,F竖向
)

help(iobr_cor_plot)
rt <- rownames_to_column(rt,var = 'ID')
identical(rt$ID,cibersort$ID)

results <- iobr_cor_plot(pdata_group = rt, id1 = "ID",
                         feature_data = cibersort, id2 = "ID",
                         group = "group",is_target_continuous = F,
                         show_plot = TRUE,ProjectID = "IOBR")

sig_group <- sig_group

#自定义分组箱线图展示
cibersort <- cibersort[,1:23]
cibersort <- column_to_rownames(cibersort,var = 'ID')
rt <- column_to_rownames(rt,var = 'ID')
rts <- merge(rt,cibersort)

rts_long <- rts %>% pivot_longer(cols = colnames(rts)[-1],
                                 names_to = "cellMarker",values_to = "Score")
# 绘制小提琴图并添加 Wilcoxon 检验的 p 值（用星号表示）
library(ggplot2)
library(ggpubr)
p <- ggboxplot(data = rts_long, x="cellMarker", y="Score", fill = "group",
               ylab="Infiltration levels",
               xlab="",
               legend.title = 'group',
               palette = c("#E74C3C", "#3498DB"),
               width=0.6, add = "none",
               outlier.shape = NA #不显示异常值
)+
  rotate_x_text(45)+stat_compare_means(aes(group = group),method="wilcox",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
                                       label = "p.signif")+ theme(
                                         axis.text.x = element_text(size = 10, # x轴文字大小
                                                                    color = "black",    # 文字颜色
                                                                    vjust = 1           # x轴文字垂直对齐
                                         ))
p

#xcell
xcell<-deconvo_tme(eset = exp, method = "xcell", #算法名
                   arrays = TRUE,  #是否是芯片数据,Currently affects 'CIBERSORT', 'svr' and 'xCell'
                   tumor = F, #是否是肿瘤, Currently affects 'EPIC'
                   perm = 10 )  #perm理论越大越精准但是耗时长

#estimate
estimate<-deconvo_tme(eset = exp, method = "estimate", #算法名
                  arrays = TRUE,  #是否是芯片数据,Currently affects 'CIBERSORT', 'svr' and 'xCell'
                  tumor = F, #是否是肿瘤, Currently affects 'EPIC'
                  perm = 10 )  #perm理论越大越精准但是耗时长

#EPIC
epic<-deconvo_tme(eset = exp, method = "epic", arrays = TRUE)
#MCPcounter
mcp<-deconvo_tme(eset = exp, method = "mcpcounter")
#xcell
xcell<-deconvo_tme(eset = exp, method = "xcell", arrays = TRUE)
#TIMER
timer<-deconvo_tme(eset = exp, method = "timer", group_list = rep("stad",dim(eset_acrg)[2]))
#quanTIseq
quantiseq<-deconvo_tme(eset = exp, tumor = TRUE, arrays = TRUE, scale_mrna = TRUE, method = "quantiseq")
#IPS
ips<-deconvo_tme(eset = exp, method = "ips", plot= FALSE)

#3.亚组cluster####
cell<-deconvo_tme(eset = exp, method = "cibersort", #算法名
                       arrays = TRUE,  #是否是芯片数据,Currently affects 'CIBERSORT', 'svr' and 'xCell'
                       tumor = F, #是否是肿瘤, Currently affects 'EPIC'
                       perm = 10 )  #perm理论越大越精准但是耗时长
tme <- tme_cluster(input = cell, features = colnames(cell)[2:23], id = "ID", scale = TRUE, method = "kmeans", max.nc = 5)
colnames(tme) <- gsub(colnames(tme), pattern = "_CIBERSORT", replacement = "")
#热图
res <- sig_heatmap(input = tme, features = colnames(tme)[3:ncol(tme)], group = "cluster", path = "result", palette = 6)
#箱线图
cols <- c("#E74C3C", "#3498DB", "#2ECC71", "#F39C12")
p <- sig_box(tme, variable = "cluster", signature = "Mast_cells_activated", 
              jitter = TRUE, cols =  cols, show_pvalue = TRUE, size_of_pvalue = 4)



#二十一、列线图可视化及校准曲线、DCA、ROC曲线评估####
rm(list=ls())
setwd('E:/R_do/阿伦生信教学代码/GEO/GSVA/')#设置工作路径
load('GSE84844.rda')

library(tidyverse)
library(rms)    #构建逻辑回归模型R包
library(pROC)   #绘制ROC曲线R包
library(rmda)   #绘制DCA曲线R包

#1.数据整理####
genes <- c("COX5B","CRKL","ATP5H","PARP14") 
exp <- exp[genes,] #提取目标基因
exp <- t(exp) #转置
dat <- cbind(rt,exp) #表达与分组合并
# dat$group <- ifelse(dat$group == "control", 0, 1) 

ddist <- datadist(dat) #计算并存储数据集中变量的统计量
options(datadist = "ddist")   #使用函数datadist()将数据打包
ddist #查看参考组是否正确

#2.模型拟合####
colnames(dat)
dat$group <- as.factor(dat$group)
paste(genes, collapse = "+")
fit <- lrm(group~COX5B+CRKL+ATP5H+PARP14, data=dat, x=T, y=T)
fit   #查看模型拟合结果

#3.列线图可视化####
nom <- nomogram(fit, fun=plogis,
                fun.at=c(0.001,0.01,0.1,0.2,0.4,0.7,0.95,0.999), #预测概率的刻度
                lp=F, funlabel="Probability of SS")
plot(nom,xfrac = 0.25,      #刻度与标度的距离
     cex.axis = 0.8,        #刻度上数值的大小
     col.grid = c("gray90") #刻度网格线
     ) 

#4.ROC####
roc_obj <- roc(group~predict(fit), data = dat)

plot.roc(roc_obj,
         print.auc=TRUE,   #输出AUC值
         auc.polygon = TRUE,auc.polygon.col = adjustcolor("#2E86AB", alpha.f = 0.1), # 绘制AUC阴影区及颜色
         main = "ROC CURVE",   #图形标题
         col= "#2E86AB",lwd = 2,    #曲线颜色
         print.thres = F,print.thres.col="red",print.thres.cex=0.5,#是否显示阈值：cut-off值字体的颜色、大小
         identity.col="#2E86AB",identity.lty = 2,identity.lwd = 1,legacy.axes=TRUE,asp = NA)# 对角线颜色、线及宽度

auc_value <- auc(roc_obj)
auc_value
ci_auc <- ci.auc(roc_obj)
ci_auc
text(x = 0.2, y = 0.1,labels = paste0("95% CI: ", round(ci_auc[1], 3), " - ", round(ci_auc[3], 3)),
     col = "red", cex = 0.8)

#5.校准曲线####
set.seed(123) 
cal <- calibrate(fit, method="boot", B=1000)
pdf('calibration.pdf',width = 6,height = 6)
plot(cal,xlab='Predicted Probability',ylab='Actual Probability',sub=F)
dev.off()

#AI结果
# 指标	优秀	   良好	       一般	      较差
# MAE	  <0.05	  0.05-0.10	  0.10-0.15	  >0.15
# MSE	  <0.005	0.005-0.01	0.01-0.02	  >0.02

#6.DCA曲线####
dat$group <- ifelse(dat$group == "control", 0, 1) 
dca_obj<- decision_curve(data= dat,  #输入数据
                       group~COX5B+CRKL+ATP5H+PARP14,# 定义自变量与因变量
                       family = binomial(link ='logit'),
                       thresholds= seq(0,1, by = 0.01),
                       confidence.intervals = 0.95)
pdf('DCA.pdf',width = 6,height = 6)
plot_decision_curve(dca_obj,
                    curve.names="Nonadherence prediction nomogram",   #曲线名称
                    xlab="Threshold probability",   #x轴名称
                    cost.benefit.axis =FALSE, col= "skyblue",
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()

#7.测试集验证####
#数据整理，统一格式
load('GSE299988.rda')
genes <- c("COX5B","CRKL","ATP5H","PARP14") 
new_dat <- exp3[genes,] #提取目标基因
new_dat <- t(new_dat) #转置
new_dat <- cbind(rt3,new_dat) #表达与分组合并

new_dat$pred_prob <- predict(fit, newdata = new_dat, type = "fitted")
range(new_dat$pred_prob) # 正常输出：0~1

#ROC曲线验证
new_roc_obj <- roc(group ~ pred_prob, data = new_dat)
plot.roc(new_roc_obj,
         print.auc=TRUE,   #输出AUC值
         auc.polygon = TRUE,auc.polygon.col = adjustcolor("#2E86AB", alpha.f = 0.1), # 绘制AUC阴影区及颜色
         main = "ROC CURVE(New Dataset)",   #图形标题
         col= "#2E86AB",lwd = 2,    #曲线颜色
         print.thres = F,print.thres.col="red",print.thres.cex=0.5,#是否显示阈值：cut-off值字体的颜色、大小
         identity.col="#2E86AB",identity.lty = 2,identity.lwd = 1,legacy.axes=TRUE,asp = NA)# 对角线颜色、线及宽度

#校准曲线验证
new_dat$group <- ifelse(new_dat$group == "control", 0, 1)
pdf('new_calibration.pdf',width = 6,height = 6)
val.prob(
  p = new_dat$pred_prob,       # 测试集的预测概率（fit模型输出的fitted值）
  y = as.numeric(new_dat$group), # 验证集的真实结局（0=control，1=pSS）
  xlab = "Predicted Probability",ylab = "Actual Probability")
dev.off()

#DCA曲线验证
new_dca_obj<- decision_curve(data= new_dat,  #输入数据
                         formula = group ~ pred_prob,# 定义自变量与因变量
                         family = binomial(link ='logit'),
                         thresholds= seq(0,1, by = 0.01),
                         confidence.intervals = 0.95)

pdf('new_DCA.pdf',width = 6,height = 6)
plot_decision_curve(new_dca_obj,
                    curve.names="Nonadherence prediction nomogram",   #曲线名称
                    xlab="Threshold probability",   #x轴名称
                    cost.benefit.axis =FALSE, col= "skyblue",
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()



#二十二、相关性/共表达分析####
rm(list=ls())
setwd('E:/R_do/阿伦生信教学代码/GEO/GEO数据下载与整理')
library(tidyverse)

load('GSE299988.rda')
str(exp3) 
exp <- exp3
rt <- rt3

#相关性分析

#保留疾病组
rt <- filter(rt,group=='Tumor')
exp <- exp[,rownames(rt)]
exp <- t(exp)

library(psych)
library(corrplot)
library(ggpubr)
#1.常规相关性####
PTPRC <- exp[,'PTPRC',drop=F] 
A1BG <- exp[,'A1BG',drop=F] 
cor <- corr.test(PTPRC,A1BG,method = 'spearman') 
r <- cor$r #相关性值
p <- cor$p #p值

dat <- cbind(PTPRC,A1BG)
p <- ggplot(data = dat, aes(x = PTPRC, y = A1BG)) +
  geom_point(alpha = 0.7, color = "black", fill = "#999999", size = 2) +
  geom_smooth(method = "lm",      # 线性拟合（若想非线性用method="loess"）
              se = TRUE,          # 显示95%置信区间
              color = "#2E86AB",  # 拟合线颜色
              linewidth = 1,      # 线宽
              fill = "#2E86AB",   # 置信区间填充色
              alpha = 0.2         # 置信区间透明度（不遮挡散点）
  ) +
  stat_cor(method = "spearman",        # 与corr.test方法一致
           label.x.npc = 0.55,         # 标注相关性x位置（0=左，1=右）
           label.y.npc = 0.05,         # 标注相关性Y位置（0=下，1=上）
           size = 4,                   # 文字大小
           color = "#A23B72",          # 文字颜色
           fontface = "bold",          # 文字加粗
           r.accuracy = 0.01,          # 相关系数保留2位小数
           p.accuracy = 0.001,         # P值保留3位小数
           label.sep = ","        # 自定义分隔符
  ) +theme_classic() 
p

#2.批量相关性####
#选择你感兴趣的基因 举例：PTPRC
Gene <- 'PTPRC'
gene <- exp[,Gene,drop=F]
EXP <- exp[,colnames(exp) != Gene,drop=F]
#循环
length(colnames(EXP))
res <- data.frame()
for (i in 1:length(colnames(EXP))) {
  genes <- EXP[,i,drop=F]
  cor <- corr.test(gene,genes,method = 'spearman')
  r <- cor$r
  p <- cor$p
  out <- data.frame(gene=colnames(gene),
                    genes=colnames(genes),
                    r值=cor$r[1],
                    p值=cor$p[1])
  res <- rbind(res,out)
}
View(res)

#筛选
#自定义 r和p
high_r <- res %>% filter(abs(r值)>0.8,p值<0.05)%>% arrange(desc(r值))

#3.批量出图####
Gene <- 'PTPRC'
gene <- exp[,Gene,drop=F]
EXP <- exp[,colnames(exp) != Gene,drop=F]
length(colnames(exp))
res <- data.frame()
for (i in 1:length(colnames(EXP))) {
  genes <- EXP[,i,drop=F]
  cor <- corr.test(gene,genes,method = 'spearman')
  r <- cor$r
  p <- cor$p
  if((abs(r)>0.8 & p<0.05)) {
    dat <- as.data.frame(cbind(gene,genes))
    p <- ggplot(data = dat, aes(x = gene, y = genes)) +
      geom_point(alpha = 0.7, color = "black", fill = "#999999", size = 2) +
      geom_smooth(method = "lm",      # 线性拟合（若想非线性用method="loess"）
                  se = TRUE,          # 显示95%置信区间
                  color = "#2E86AB",  # 拟合线颜色
                  linewidth = 1,      # 线宽
                  fill = "#2E86AB",   # 置信区间填充色
                  alpha = 0.2         # 置信区间透明度
      ) +
      stat_cor(method = "spearman",        # 与corr.test方法一致
               label.x.npc = 0.55,         # 标注相关性x位置（0=左，1=右）
               label.y.npc = 0.05,         # 标注相关性Y位置（0=下，1=上）
               size = 4,                   # 文字大小
               color = "#A23B72",          # 文字颜色
               fontface = "bold",          # 文字加粗
               r.accuracy = 0.01,          # 相关系数保留2位小数
               p.accuracy = 0.001,         # P值保留3位小数
               label.sep = ","        # 自定义分隔符
      ) +theme_classic()+labs( x = colnames(gene),y = colnames(genes)) 
    name <- paste0(colnames(gene),'_',colnames(genes))
    pdf(file = paste0(name, ".pdf"),width=5,height=5)
    print(p)
    dev.off()
    out <- data.frame(gene=colnames(gene),
                      genes=colnames(genes),
                      r值=cor$r[1],
                      p值=cor$p[1])
    res <- rbind(res,out)
  }}
View(res)

#二十三、安捷伦原始数据处理####
#示例数据：GSE289533
rm(list=ls())
setwd('E:/R_do/阿伦生信教学代码/GEO/安捷伦/') #设置工作路径
library(tidyverse)
library(limma)
#1.创建分组数据####
Filename <- list.files(path = './GSE289533_RAW/')
Filename
sample_id <- substr(Filename,1,10)
sample_id
group <- c(rep('control',3),
           rep('Cyclobutrifluram',3))
group

rt <- data.frame(row.names =sample_id,
                 group=group)

#2.原始数据读取####
data <- read.maimages(files = Filename,
                      source = 'agilent',
                      path = './GSE289533_RAW/', #进一步指定工作路径
                      names = sample_id,
                      other.columns = 'gIsWellAboveBG',
                      green.only = T) #cy3单色选择T,双色芯片cy3+cy5选择F
data$rt <- rt #把分组数据放入data文件

#校准前箱线图
boxplot(data$E,las=2,outline=F,col=factor(rt$group))

#3.背景校准####
bg <- backgroundCorrect(RG=data,
                        method = 'normexp',
                        offset = 50,
                        normexp.method = 'mle')
#标准化
#单色cy3使用normalizeBetweenArrays 双色cy3+cy5使用normalizeWithinArrays
exp <- normalizeBetweenArrays(bg,method = 'quantile')
View(exp$E) 
range(exp$E)
#校准后箱线图
boxplot(exp$E,las=2,outline=F,col=factor(rt$group))

table(exp$genes$ControlType)
# 0是检测数据，1是positive对照，-1是negative对照

#增加探针
# 筛选仅检测探针（排除对照探针）
keep_probe <- exp$genes$ControlType == 0
exp_E <- exp$E[keep_probe, ] # 只保留检测探针的表达矩阵
# 给表达矩阵添加探针名（仅保留有效探针）
rownames(exp_E) <- exp$genes$ProbeName[keep_probe]
# 转为数据框
exp <- as.data.frame(exp_E)
range(exp) # 此时exp仅包含有效检测探针，无对照

#4.基因注释，GPL导入####
gpl<-read.table('GPL21163-3202.txt',header=TRUE,fill=T,sep='\t',comment.char='#',stringsAsFactors=FALSE,
                quote='')  #去掉前几列没用的行
View(gpl)
colnames(gpl) #查看列名
ids <- gpl[,c('ID','GENE_SYMBOL')]  #需要查看表格填写正确的列名
colnames(ids) <- c('probe_id','symbol') #重命名
ids <- ids[ids$symbol!='',]#删除空白
ids <- ids[ids$probe_id %in% rownames(exp),]
exp2 <- exp[ids$probe_id,]
table(rownames(exp2)==ids$probe_id)
exp3 <- cbind(ids,exp2)
#5.基因去重####
exp3 <- aggregate( . ~ symbol,exp3[,-1],max) #max最大，或者mean平均 随你
rownames(exp3) <- exp3$symbol
exp3 <- exp3[,-c(1)]

#6.数据保存####
save(exp3,rt,file = 'GSE289533.rda')
rm(list = ls())
load('GSE289533.rda')
