#一、环境配置####
library(GEOquery)
library(tidyverse)
library(DESeq2)
library(limma)
library(pheatmap)
library(ggpubr)
library(cowplot)

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
setwd('E:/scleraBioinformaticsAnalysis/Datasets/3_Myopia_Tissues/High/GSE200053')
library(tidyverse)
exp <- data.table::fread("GSE200053_mRNA_sclera_normalized_counts.txt",header = T)
table(is.na(exp))
exp[is.na(exp)] <- 0
table(is.na(exp))
exp <- column_to_rownames(exp,var ='Gene Symbol')

#1.设置分组信息####
sample_id <- colnames(exp)
sample_id
group <- rep(c('L','HLL','HL'),each=6)
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
res <- results(dds,contrast = c("group","H","L"))
#提取差异表达结果
res1 <- as.data.frame(res)
#结果筛选padj<0.05 & abs(log2Foldchange)>=1
degs <- filter(res1,padj<0.05 & abs(log2Foldchange)>=1)
view(degs)
#绘制离散图
plotDispEsts(dds)#p值直方图，直观看下p值显著基因的多少
hist(res$padj)

#3.count数据vst转化####
vst <- vst(dds,blind = F)
expl <- assay(vst)
range(exp1)
exp1 <- as.data.frame(exp1)
save(exp,exp1,rt,file = 'GSE309139.rda')
rm(list = ls())
load('GSE309139.rda')

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

#四、富集分析####

# BiocManager::install('clusterProfiler')
rm(list = ls())
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

setwd('E:/scleraBioinformaticsAnalysis/Datasets/3_Myopia_Tissues/High/GSE200053')
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
  Column = "ENTREZID",  #转换为ENTREZID作为中间步骤
  multivals = "first"   #多个匹配时取第一个 
)

entrez_df<- data.frame(ENTREZID= entrez_ids,check.names = FAentrez_df<- na.omit(entrez_df)<- entrez_df %>% rownames_to_column(var = 'SYMBoL')entrez_df
                       #合并
                       gene_entrez2 <- rbind(gene_entrez,entrez_df)

#四、故障排查####
#pdf()未正常生成图像文件
graphics.off()
