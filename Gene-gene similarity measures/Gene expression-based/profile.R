## 20200624bywyq
setwd("F:/R workspace/Gene-gene similarity measures/Gene expression-based/")

drugall = read.csv("trt_sh.cgs_L7epdata_extract.csv", header = T, row.names = 1)
genesymbol = read.csv("12328genesymbol.csv", header = F)
drugall = apply(drugall, 2, rank, ties.method = "first")

ncol(geneall)
ES_matrix_A = matrix(data = NA, ncol = 1766, nrow = 1766)
for(i in 1:(ncol(drugall)-1)){
  up_i = drugall[,i][1:250]
  down_i = drugall[,i][12079:12328]
  for(j in (i+1):ncol(drugall)){
    up_es=GSEA.EnrichmentScore(gene.list = drugall[,j],gene.set = up_i,weighted.score.type = 0,correl.vector = genesymbol)[[1]]#上调的signature与PRL之间的ES
    down_es=GSEA.EnrichmentScore(gene.list = drugall[,j],gene.set = down_i,weighted.score.type = 0,correl.vector = cor_vect)[[1]]#下调的signature与PRL之间的ES
    ES_matrix_A[i,j] = ((up_es-down_es)/2)
  }
  print(i)
}
write.csv(ES_matrix_A, file = "ES_matrix_A.csv")

ES_matrix_B = matrix(data = NA, ncol = 1766, nrow = 1766)
for(i in 1:(ncol(drugall)-1)){
  for(j in (i+1):ncol(drugall)){
    up_j = drugall[,j][1:250]
    down_j = drugall[,j][12079:12328]
    up_es=GSEA.EnrichmentScore(gene.list = drugall[,i],gene.set = up_j,weighted.score.type = 0,correl.vector = genesymbol)[[1]]#上调的signature与PRL之间的ES
    down_es=GSEA.EnrichmentScore(gene.list = drugall[,i],gene.set = down_j,weighted.score.type = 0,correl.vector = cor_vect)[[1]]#下调的signature与PRL之间的ES
    ES_matrix_B[i,j] = ((up_es-down_es)/2)
  }
  print(i)
}
write.csv(ES_matrix_B, file = "ES_matrix_B.csv")

diag(ES_matrix_A) = 1
diag(ES_matrix_B) = 1
save(ES_matrix_A, file = "ES_matrix_A.Rdata")
save(ES_matrix_B, file = "ES_matrix_B.Rdata")

for(i in 1:1765){
  for(j in (i+1):1766){
    ES_matrix_A[j,i] = ES_matrix_A[i,j]
    ES_matrix_B[j,i] = ES_matrix_B[i,j]
  }
}

ES_matrix = (ES_matrix_A + ES_matrix_B)/2
write.csv(ES_matrix, file = "ES_matrix.csv")


GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. 
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
  #   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
  #   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
  no.tag.indicator <- 1 - tag.indicator 
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  }
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector**alpha)
  sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > - min.ES) {
    #      ES <- max.ES
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    #      ES <- min.ES
    ES <- signif(min.ES, digits=5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}




