rm(list = ls())
setwd("E:\scleraBioinformaticsAnalysis\Datasets\3_Myopia_Tissues\High\GSE200053")
library(tidyverse)
exp <- data.table::fread("GSE309139_raw_counts.txt", header = table(is.na(exp))
exp[is.na(exp)] <- 0
table(is.na(exp))
exp <- column_to_rownames(exp,var ='Gene')