###############################################
# This script takes the transNOAH cancer data set 
# and applies the test by calculating the 
# proposed statistic, and then proceeds
# to draw conclusions. Additionally, it applies
# comparison tests from Zhong and Chen 2011
###############################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GEOquery")
# libraries
library(GEOquery)
library(data.table)
library(glmnet)
library(SIHR)
#Set the working directory of this script in advance.
source("Real_functions.R")
# data loading for transNOAH to be done once and saving
#Sys.setenv("VROOM_CONNECTION_SIZE" = 500000)
#readr::local_edition(1)
#transNOAH = GEOquery::getGEO("GSE50948")
#save(transNOAH, file = "transNOAH.RData")

# loading semi-processed transNOAH data
load("transNOAH.rData")

# processing data
tn = transNOAH$GSE50948_series_matrix.txt.gz
data = assayData(tn)$exprs
#View(Biobase::fData(tn))
####################Process the data if you do not directly load the HDGO.csv file.
#GOstep1 = Biobase::fData(tn)$`Gene Ontology Molecular Function`
#uniqueGO1 = unique(unlist(strsplit(unlist(GOstep1), "[^0-9]+")))
#uniqueGO = sort(uniqueGO1[nchar(uniqueGO1) == 7])

# get only the subset of info where gene title contains cancer
#xsubset = Biobase::fData(tn)[which(Biobase::fData(tn)$`Gene Title` %like% "*cancer*"|Biobase::fData(tn)$`Target Description` %like% "*cancer*"),]$ID
#GOterms = rep(0,length(uniqueGO))
#for(i in 1:length(uniqueGO)){
#  xsubset = Biobase::fData(tn)[which(Biobase::fData(tn)$`Gene Ontology Molecular Function` %like% uniqueGO[i]),]$ID
#  GOterms[i] = length(xsubset)
#}
#HDGO = data.frame("GO_term" = uniqueGO[which(GOterms > 158 & GOterms < 2092)], "terms" = GOterms[which(GOterms > 158 & GOterms < 2092)])
#HDGO = HDGO[order(HDGO$terms),]
#write.csv(HDGO,"HDGO.csv")
####################
HDGO = read.csv("HDGO.csv")
# loading the data
yrows = which(rownames(data) %in% c("204531_s_at", "211851_x_at"))

y1 = data[yrows[1],]
ybar1 = mean(y1)

x1 = data[-(which(rownames(data) %in% c("204531_s_at", "211851_x_at"))),]

Record_y1_stat =  NULL
Record_probe = matrix(rep(0,2),1,2)
for(j in 1:125){
  time0 = Sys.time()
  xpick = Biobase::fData(tn)[which(Biobase::fData(tn)$`Gene Ontology Molecular Function` %like% HDGO$GO_term[j]),]$ID
  x = x1[which(rownames(x1) %in% xpick),]
  probe = rownames(x)
  n = ncol(x)
  p = nrow(x)
  r = p
  C = diag(1,p)
  Omega = diag(1,p)
  xbar = rowSums(x)/n
  ybar1 = mean(y1)
  y11 = y1 - ybar1
  xnew = x
  for(i in 1:n){
    xnew[,i] = (x[,i] - xbar)
  }
  X = t(xnew)
  gamma = rep(0,p)
  Result0 = Inference2(X,y11,"Full",1)
  stat_full = Result0$stat_npe
  #######################The Proposed method
  Result3 = Infer_decor_full_global(X,y11)
  Power_stat = Result3$stat^2
  asd1 = ifelse(Power_stat>=2*log(p) + 2*sqrt(log(n))*log(log(p)),Power_stat,0)
  probe_active = probe[which(asd1!=0)]
  length_id = 1
  if(length(probe_active)!=0)
  {
    probe_id = c(HDGO$GO_term[j],probe_active)
    maxlength = max(dim(Record_probe)[2],length(probe_id))
    length(probe_id) = maxlength
    length_id = max(length_id,length(probe_id))
    Record_probe = rbind(Record_probe,probe_id)
  }
  PE =  sqrt(p)*sum(asd1)
  #######################Comparison Method
  # Zhong and Chen 2011
  TstatsZC1 = zhongchen2011(X = x, y = y1, beta = 0, delta = 0, Sigma = NULL, small.sig = NULL, T0 = FALSE)
#####################
  Record_y1_stat = rbind(Record_y1_stat, c(HDGO$GO_term[j],abs(stat_full), abs(stat_full)+PE,abs(TstatsZC1[2])))
  print(c(j,HDGO$GO_term[j],stat_full, stat_full+PE,TstatsZC1[2],probe_active))
  write.csv(Record_probe,"FULL_Y1_probe.csv")
  write.csv(Record_y1_stat,paste("FULL_Y1_stat.csv",sep=""))
}
#Record_y1_stat = read.csv("FULL_Y1_stat.csv",row.names = 1)
###Without Bonferroni correction
apply(Record_y1_stat[,-c(1:2)], 2, function(x) HDGO$GO_term[(which(x>=qnorm(1-0.05/2)))])
apply(Record_y1_stat[,-c(1:2)], 2, function(x) length(which(x>=qnorm(1-0.05/2))))
###With Bonferroni correction
apply(Record_y1_stat[,-c(1:2)], 2, function(x) HDGO$GO_term[(which(x>=qnorm(1-0.05/250)))])

apply(Record_y1_stat[,-c(1:2)], 2, function(x) length(which(x>=qnorm(1-0.05/250))))
#the id 16779 51082 46983  8134  3824 16301