#一、环境配置####
library(GEOquery)
library(tidyverse)
library(DESeq2)
library(limma)
library(pheatmap)
library(ggpubr)
library(cowplot)

#色板
morandi <- c(
  "#6C8E8A", # 灰青
  "#8FA6A0", # 雾青绿
  "#A3B18A", # 灰橄榄
  "#B7B7A4", # 亚麻灰
  "#C2B8A3", # 沙灰
  "#BFAE9C", # 米棕
  "#C7A9A3", # 灰陶粉
  "#B7A3B5", # 灰藕紫
  "#9FA8C3", # 灰蓝紫
  "#8DA3B3", # 雾钢蓝
  "#AEB8B1", # 冷雾灰
  "#C7C7C7", # 中性灰
  "#9C9C8A", # 岩灰橄榄
  "#B5C0B7", # 鼠尾草灰
  "#C8BFB6", # 灰陶土
  "#A6A29B"  # 暖岩灰
)
palette(morandi)

rm(list=ls())
library(GEOquery)
library(tidyverse)
setwd('E:/scleraBioinformaticsAnalysis/Datasets/3_Myopia_Tissues/High/GSE200053') #设置工作
gset <- getGEO("GSE200053", destdir = ".", AnnotGPL = F, getGPL = F)
class(gset)
gset[[1]]#查看信息
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
GPL <- getGEO(filename = 'GPL17021_family.soft.gz',destdir ='.')
#获取注释信息
gpl <- GPL@dataTable@table
View(gpl)
colnames(gpl)#查看列名

#二、高通量数据(Count)分析####

rm(list = ls())
setwd('E:/scleraBioinformaticsAnalysis/Datasets/3_Myopia_Tissues/High/GSE280609')
library(tidyverse)
exp <- data.table::fread("GSE280609_sclerarna_count.csv",header = T)
table(is.na(exp))
exp[is.na(exp)] <- 0
table(is.na(exp))

#如果 RAW_COUNTS 有重复
library(data.table)
dt <- as.data.table(exp)
dt <- dt[nzchar(gene_id)]
exp <- dt[, lapply(.SD, sum), by=gene_id]

exp <- column_to_rownames(exp,var ='gene_id')

#1.设置分组信息####
sample_id <- colnames(exp)
sample_id
group <- rep(c('control','myopia'),each=5)
group
rt <- data.frame(row.names =sample_id,group=group)
identical(rownames(rt),colnames (exp))
# [1] TRUE

#2.DEseg2分析####
# BiocManager::install('DESeq2')
library(DESeq2)
#创建DESeq2分组对象
rt$group <- factor(rt$group)
rt$group
#构建DESeg2对象
dds <- DESeqDataSetFromMatrix(
  countData = exp,
  colData = rt,
  design = ~ group
)
#进行低表达基因过滤(自定义)
filter <- rowSums(counts(dds))>1 #筛选每行加起来大于1的基因
dds <- dds[filter, ]
#差异表达分析
dds <- DESeq(dds)#指定具体的对比组，格式为c("因子名"，"实验组"，"对照组")
res <- results(dds,contrast = c("group","myopia","control"))
#提取差异表达结果
res1 <- as.data.frame(res)
#结果筛选padj<0.05 & abs(log2Foldchange)>=1
degs <- filter(res1,padj<0.05 & abs(log2FoldChange)>=1)
view(degs)
#绘制离散图
plotDispEsts(dds)#p值直方图，直观看下p值显著基因的多少
hist(res$padj)

#3.count数据vst转化####
vst <- vst(dds,blind = F)
exp1 <- assay(vst)
range(exp1)
exp1 <- as.data.frame(exp1)
save(exp,exp1,rt,degs,file = 'GSE280609.rda')
rm(list = ls())
load('GSE280609.rda')

#4.火山图####
DEG = na.omit(res1) #differently expressed genes
write.table(DEG, file='DEG.txt',sep = "\t",row.names = T,col.names = NA, quote = F)
#进行注释 foldchange

#通常1ogFC设置1,adj.P.Val < 0.05
logFC_cutoff <- 1

type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)

DEG$change = ifelse (type1, "DOWN", ifelse (type2, "UP", "NOT"))
write.table(DEG, file='DEG2.txt',sep = "\t",row.names = T,col.names = NA, quote = F)
table(DEG$change)

library(ggplot2)
library(cowplot)
library(ggrepel)
DEG <- mutate(DEG,Gene_symbol=row.names (DEG))
UPDEG <- subset (DEG,change=='UP')
UPDEG_5 <- top_n(x = UPDEG,n = -5,wt = padj)
DOWNDEG <- subset (DEG, change=='DOWN')
DOWNDEG_5 <- top_n(x = DOWNDEG,n = -5,wt = padj)

p <- ggplot(
  data = DEG,
  aes(x = log2FoldChange, y = -log10(padj))
) +
  geom_point(alpha = 0.5, size = 4.5,
             aes(color = change)) +
  ylab("-log10(padj)") +
  scale_color_manual(values = c("#1F77B4","grey","#FF7F0E")) +   
  
  geom_vline(xintercept = c(-2, 2),
             lty = 5, col ="black", linewidth = 0.5) +           
  
  geom_hline(yintercept = -log10(0.05),
             lty = 4, col = "black", linewidth = 0.5) +          
  
  theme_half_open() +                                             
  
  geom_text_repel(
    data = UPDEG_5,
    aes(label = Gene_symbol),
    vjust = 1.5, size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.35, "lines"),
    segment.color = 'grey50',
    nudge_y = 0.1,
    direction = "y"
  ) +                                                             
  
  geom_text_repel(
    data = DOWNDEG_5,
    aes(label = Gene_symbol),
    vjust = 1.5, size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.35, "lines"),
    segment.color = 'grey50',
    nudge_y = 0.1,
    direction = "y"
  )

pdf('volcano.pdf', width = 10, height = 6)
p
dev.off()

#5.热图####
library(pheatmap)
identical(rownames(rt),colnames(exp))
DEG <- read.table("DEG2.txt",sep = "\t",check.names = F, stringsAsFactors = F,header = T,row.names = 1)
table(DEG$change)
#作图
diff <- rownames(DEG)[DEG$change !="NOT"]
exp_diff <- exp[diff, ]
range(exp_diff)

#展示特定基因
diff = DEG[DEG$change !="NOT",]
up <- diff %>% top_n(25,logFC)
dw <- diff %>% top_n(-25,logFC)
all <-c(rownames(up),rownames(dw))
exp_diff2<- exp[all,]

pdf(file = 'heatmap.pdf',width=6,height=20)
pheatmap(exp_diff,
         annotation_col=rt,
         color = colorRampPalette((c("blue","white","red")))(100),
         cluster_cols = F,
         cluster_rows = T,
         show_colnames = F,
         show_rownames = T,
         border_color = NA,
         scale = 'row')
dev.off() #scale='row'进行标准化

#三、limma差异分析####
#FPKM -> TPM
exp <- log2(exp + 1)

library(limma)
library(pheatmap)
library(tidyverse)
library(ggpubr)
#数据读取
load('GSE299988.rda')
str(exp)
exp <- exp %>% mutate(across(everything(), as.numeric))
str(exp)
#确认行名列名一致
identical(rownames(rt),colnames (exp))
#1.火山图####
group_list <-factor(rt$group,levels = c("control","myopia"))
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
write.table(DEG, file='DEG1.txt',sep = "\t",row.names = T,col.names = NA, quote = F)
#进行注释 foldchange
#通常1ogFC设置1,adj.P.Val < 0.05
logFC_cutoff <- 1
type1 = (DEG$adj.P.Val < 0.05)&(DEG$logFC < -logFC_cutoff)
type2 = (DEG$adj.P.Val < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse (type1, "DOWN", ifelse (type2, "UP", "NOT"))
write.table(DEG, file='DEG2.txt',sep = "\t",row.names = T,col.names = NA, quote = F)
table(DEG$change)
library(ggplot2)
library(cowplot)
library(ggrepel)
DEG <- mutate(DEG,Gene_symbol=row.names (DEG))
UPDEG <- subset (DEG,change=='UP')
UPDEG_5 <- top_n(x = UPDEG,n = -5,wt = P.Value)
DOWNDEG <- subset (DEG, change=='DOWN')
DOWNDEG_5 <- top_n(x = DOWNDEG,n = -5,wt = P.Value)

p <- ggplot(
  data = DEG,
  aes(x = logFC, y = -log10(adj.P.Val))
) +
  geom_point(alpha = 0.5, size = 4.5,
             aes(color = change)) +
  ylab("-log10(adj.P.Val)") +
  scale_color_manual(values = c("#1F77B4","grey","#FF7F0E")) +   
  
  geom_vline(xintercept = c(-1, 1),
             lty = 5, col ="black", linewidth = 0.5) +           
  
  geom_hline(yintercept = -log10(0.05),
             lty = 4, col = "black", linewidth = 0.5) +          
  
  theme_half_open() +                                             
  
  geom_text_repel(
    data = UPDEG_5,
    aes(label = Gene_symbol),
    vjust = 1.5, size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.35, "lines"),
    segment.color = 'grey50',
    nudge_y = 0.1,
    direction = "y"
  ) +                                                             
  
  geom_text_repel(
    data = DOWNDEG_5,
    aes(label = Gene_symbol),
    vjust = 1.5, size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.35, "lines"),
    segment.color = 'grey50',
    nudge_y = 0.1,
    direction = "y"
  )

pdf('volcano.pdf', width = 10, height = 6)
p
dev.off()

#2.热图####
library(pheatmap)
identical(rownames(rt),colnames(exp))
DEG <- read.table("DEG2.txt",sep = "\t",check.names = F, stringsAsFactors = F,header = T,row.names = 1)
table(DEG$change)
#作图
diff <- rownames(DEG)[DEG$change !="NOT"]
exp_diff <- exp[diff, ]
range(exp_diff)

#展示特定基因
diff = DEG[DEG$change !="NOT",]
up <- diff %>% top_n(25,log2FoldChange)
dw <- diff %>% top_n(-25,log2FoldChange)
all <-c(rownames(up),rownames(dw))
exp_diff2<- exp[all,]

pdf(file = 'heatmap.pdf',width=6,height=20)
pheatmap(exp_diff,
         annotation_col=rt,
         color = colorRampPalette((c("blue","white","red")))(100),
         cluster_cols = F,
         cluster_rows = T,
         show_colnames = F,
         show_rownames = T,
         border_color = NA,
         scale = 'row')
dev.off() #scale='row'进行标准化

#四、富集分析####

# BiocManager::install('clusterProfiler')
rm(list = ls())
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db) #Homo Sapien
library(org.Mm.eg.db) #Mus Musmuculus
library(enrichplot)
library(ggplot2)

setwd('E:/scleraBioinformaticsAnalysis/Datasets/3_Myopia_Tissues/High/GSE280609')
DEG <- read.table("DEG2.txt",sep = "\t",check.names = F,stringsAsFactors = F, header = T,row.names=1)
table(DEG$change)
diff <- rownames (DEG) [DEG$change!='NOT']
#将基因ID从Symbol转换为EntrezID
gene_entrez <- bitr(diff, fromType ="SYMBOL",toType="ENTREZID",OrgDb = org.Hs.eg.db)
cat("所有基因ID转换",nrow(gene_entrez),"/",length(diff))
untrans <- setdiff(diff,gene_entrez$SYMBOL)
untrans
#查看基因名是否为官方名
isValid <- "c9orf92" %in% keys(org.Hs.eg.db, keytype="SYMBOL")
isValid#查看T还是F

#使用mapIds函数转换
library(AnnotationDbi)
entrez_ids <- mapIds(
  x = org.Hs.eg.db,
  keys = untrans,       #输入的基因别名列表
  keytype = "ALIAS",    #明确指定输入类型为别名
  column = "ENTREZID",  #转换为ENTREZID作为中间步骤
  multivals = "first"   #多个匹配时取第一个 
)


entrez_df <- data.frame(ENTREZID= entrez_ids,check.names = FALSE)
entrez_df <- na.omit(entrez_df)
entrez_df <- entrez_df %>% rownames_to_column(var = 'SYMBOL')
#合并
gene_entrez2 <- rbind(gene_entrez,entrez_df)

#对于不存在OrgDb的物种
library(biomaRt)
diff <- rownames (DEG) [DEG$change!='NOT']

mart_cpo <- useEnsembl(biomart = "genes",dataset = "cporcellus_gene_ensembl")
mart_hsa <- useEnsembl(biomart = "genes",dataset = "hsapiens_gene_ensembl")

cpo_map <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters    = "external_gene_name",
  values     = diff,
  mart       = mart_cpo
)

attrs_cpo <- listAttributes(mart_cpo)$name
grep("hsapiens|ortholog", attrs_cpo, value=TRUE)[1:30]

orth1 <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "hsapiens_homolog_ensembl_gene",              # 关键：human ortholog Ensembl gene
    "hsapiens_homolog_associated_gene_name"      # human symbol（可选））
  ),
  filters = "ensembl_gene_id",
  values  = unique(cpo_map$ensembl_gene_id),
  mart    = mart_cpo
)

hsa_entrez <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = unique(orth1$hsapiens_homolog_ensembl_gene),
  mart       = mart_hsa
)

# 合并
orth <- merge(orth1, cpo_map, by="ensembl_gene_id", all.x=TRUE)
orth <- merge(orth, hsa_entrez, by.x="hsapiens_homolog_ensembl_gene", by.y="ensembl_gene_id", all.x=TRUE)

# 结果：cpo_symbol -> human_symbol -> human_entrez
orth <- orth[, c("external_gene_name", "hsapiens_homolog_associated_gene_name", "hgnc_symbol", "entrezgene_id")]
colnames(orth)[1] <- "cpo_symbol"
head(orth)

gene_entrez2 <- na.omit(orth)
colnames(gene_entrez2) <- c("cpo_symbol", "hsa_name","hsa_symbol", "ENTREZID")


#1.GO富集####
GO <- enrichGO(gene = gene_entrez2$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType ="ENTREZID",
               ont ="ALL",            # "BP","MF"，"CC"或"ALL"
               pAdjustMethod ="BH",
               pvalueCutoff =0.1,
               qvalueCutoff =0.2,readable = T)
GO_res <- GO@result
write.table(GO_res, file="GO_res.txt",sep="\t",quote=F,row.names = F)
#简单作图
barplot(GO, showCategory = 20)
dotplot(GO, showCategory = 20)
#分类展示
barplot(GO, showCategory = 5,split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')

#展示指定条目
terms_keep <- GO_res$Description[ grep("mito", GO_res$Description, ignore.case = TRUE) ]

barplot(GO, showCategory = terms_keep, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free")

#2.KEGG分析####
KEGG <- enrichKEGG(gene = gene_entrez2$ENTREZID,
                   organism='hsa', #小鼠为'mmu'
                   pvalueCutoff = 0.2,
                   qvalueCutoff = 0.2)
KEGG <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType = "ENTREZID") 
KEGG_res <- KEGG@result
write.table(KEGG_res, file="KEGG_res.txt",sep="\t",quote=F,row.names = F)
#简单作图
kdf <- as.data.frame(KEGG)
ids <- kdf$ID[kdf$category != "Human Diseases"]

dotplot(KEGG, showCategory = ids)

barplot(KEGG,showCategory = 20)
dotplot(KEGG, showCategory = 20)
#通路图可视
化#BiocManager::install('pathview')
library(pathview)
DEG <- rownames_to_column(DEG,Var = 'SYMBOL')
pgenes <- inner_join(DEG,gene_entrez2,by = 'SYMBOL')
pgenes <- dplyr::select(pgenes,ENTREZID,P.Value)

pvalue <- -log10(pgenes$P.Value)
gene_names <- as.character(pgenes$ENTREzID)

gene_data <- pvalue
names(gene_data) <- gene_names
#选择兴趣的通路
hsa04978 <- pathview(gene.data=gene_names,
                     pathway.id="hsa04978",#通路名字
                     species ="hsa",      #代表人类
                     out.suffix = "pathview",
                     limit= list(gene= 10, cpd= 1)) #基因数量限制，cpd化合物数量限制

#3.弦图####

library(GOplot)
library(ggplot2)
library(tidyr)

DEG <- read.table("DEG2.txt",sep = "\t",check.names = F, stringsAsFactors = F,header = T,row.names = 1)
KEGG_res <- read.table("KEGG_res.txt", sep = "\t",header = TRUE, stringsAsFactors = FALSE)
KEGG <- setReadable(KEGG, OrgDb = org.Mm.eg.db, keyType = "ENTREZID") 
KEGG_res <- KEGG@result
kegg_list <- strsplit(KEGG_res$geneID, "/")

## 构建 Term-Gene
top_n <- 10
top_n <- min(top_n, nrow(KEGG))

a11 <- do.call(
  rbind,
  lapply(seq_len(top_n), function(i) {
    data.frame(
      Term  = KEGG$Description[i],
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
DEG2 <- DEG2[order(DEG2$logFC, decreasing = TRUE), , drop = FALSE]

a12 <- a12[rownames(DEG2), , drop = FALSE]
a12$logFC <- DEG2$logFC

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

#五、GSEA富集分析####
# https://www.gsea-msigdb.org/gsea/msigdb/index.jsp

#1.常规GSEA####
rm(list = ls())
library(tidyverse)
library(limma)
setwd('E:/scleraBioinformaticsAnalysis/Datasets/3_Myopia_Tissues/High/GSE200053')#设置工作路径
load('GSE299988.rda')
str(exp)
exp <- exp %>% mutate(across(everything(), as.numeric))
str(exp)
#确认行名列名一致
identical(rownames(rt),colnames(exp))

#组间差异分析
group_list <-factor(rt$group,levels = c("Normal","Tumor"))
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
DEG <- rownames_to_column(DEG,var= 'Gene')
#整理输入数据
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
#转换ENTREZID
# genelist <- bitr(DEG$Gene, fromType="SYMBOL",
#                  toType="ENTREZID", OrgDb='org.Hs.eg.db')
# DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))
geneList <- DEG[,'logFC']
names(geneList) <- as.character(DEG[,'Gene'])
#排序:
geneList <- sort(geneList, decreasing= TRUE)
geneList

#读取gmt文件
gmt <- read.gmt('mh.all.v2026.1.Mm.symbols.gmt')
GSEA <- GSEA(geneList,
             TERM2GENE = gmt,
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             minGSSize = 5,
             maxGSSize = 500)
res <- GSEA@result
write.csv(res,file = 'H_GSEA.csv')
p <- gseaplot2(GSEA,geneSetID = 'HALLMARK_ALLOGRAFT_REJECTION',
               title = 'HALLMARK_ALLOGRAFT_REJECTION',color ="red")
p
p1<- gseaplot2(GSEA,geneSetID = 1:3) #geneSetID可以选择多个
p1

#2.单基因GSEA####
rm(list = ls())
library(tidyverse)
library(limma)
setwd('E:/scleraBioinformaticsAnalysis/Datasets/3_Myopia_Tissues/High/GSE200053')#设置工作路径
load('GSE299988.rda')
str(exp)
exp <- exp %>% mutate(across(everything(), as.numeric))
str(exp)
#确认行名列名一致
identical(rownames(rt),colnames(exp))

#提取 自定义基因，整理高低表达组
rtl <- as.data.frame(t(exp['PTPRc',]))
rt1$group <- ifelse(rt1$PTPRC > median(rt1$PTPRc),'up','down')
rtl <- arrange(rtl,'PTPRc',group)
identical (rownames (rt1),colnames(exp))
#if FALSE
exp <- exp[,rownames(rt1)]
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
DEG <- rownames_to_column(DEG,var= 'Gene')
#整理输入数据
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
#转换ENTREZID
# genelist <- bitr(DEG$Gene, fromType="SYMBOL",
#                  toType="ENTREZID", OrgDb='org.Hs.eg.db')
# DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))
#排序:
geneList <- sort(geneList, decreasing= TRUE)
geneList

#读取gmt文件
gmt <- read.gmt('h.al1.v2025.1.Hs.symbols.gmt')
GSEA <- GSEA(geneList,
             TERM2GENE = gmt,
             pvaluecutoff = 0.05,
             pAdjustMethod = "BH",
             minGSSize = 20,
             maxGSSize = 500)
res <- GSEA@result
write.csv(res,file = 'H_GSEA.csv')
p <- gseaplot2(GSEA,geneSetID = 'HALLMARK_ALLOGRAFT_REJECTION',
               title = 'HALLMARK_ALLOGRAFT_REJECTION',color ="red")
p
p1<- gseaplot2(GSEA,geneSetID = 1:3) #
p1




#四、故障排查####
#pdf()未正常生成图像文件
graphics.off()
