#!/usr/bin/Rscript

setwd("F:/R workspace/Gene-gene similarity measures/PPI-based/")
####Calculate the distance of genes in the PPI network####
library(dplyr)
library(doParallel)
library(igraph)
ppi=read.table("ppi.tsv",header = T,sep="\t",check.names = F,quote = "")
ppi=as.matrix(ppi[,1:2])
# ppi_gene is our target gene (genes in our dataset)
ppi_gene=as.matrix(read.table("target_Entrez_ID.txt",header = F,sep="\t",check.names = F,quote = ""))
all_gene=unique(c(ppi[,1],ppi[,2]))
rep_gene=intersect(all_gene,ppi_gene)
ppi_gene=rep_gene
rm(rep_gene)
# ppi_gene=c(as.numeric(ppi[,1]),as.numeric(ppi[,2]))
# ppi_gene=as.matrix(unique(ppi_gene))
# ppi_gene=sample(ppi_gene,size=100)
# Create an undirected graph
g <- make_graph(t(ppi[,1:2]),directed = F)
a=rep(1,nrow(ppi))
g <- set_graph_attr(g,'weight',a) 
ppi_distance=matrix(0,length(ppi_gene),length(ppi_gene))
colnames(ppi_distance)=rownames(ppi_distance)=ppi_gene
print("start loop")
cl <- makeCluster(4)
registerDoParallel(cl)
ppi_distance= foreach (i = 1:length(ppi_gene),.combine=rbind) %dopar%
  {
    ppi_row=NULL
    for(j in 1:length(ppi_gene))
    {
      gene_1=ppi_gene[i]
      gene_2=ppi_gene[j]
      ppi_path=igraph::shortest_paths(g,gene_1,gene_2,weights = a)
      path=length(unlist(ppi_path[[1]]))-1
      ppi_row=c(ppi_row,path)
    }
    ppi_row
  }
stopCluster(cl)
print("stop loop")
write.table(ppi_distance,file = "ppi_distance.txt",quote = F,sep = "\t",col.names = T,row.names = T)
write.csv(ppi_distance, file = "ppi_distance.csv")
save(ppi_distance,file = "ppi_distance.Rdata")
print("wirte distance")
####Calculate gene-gene similarity based on ppi distance####
A=0.9*exp(1)
b=1
ppi_sim=matrix(0,length(ppi_gene),length(ppi_gene))
colnames(ppi_sim)=rownames(ppi_sim)=ppi_gene
timestart<-Sys.time()
for(i in 1:nrow(ppi_sim))
{
  print(i)
  for(j in i:nrow(ppi_sim))
  {
    ppi_sim[i,j]=A*exp(-b*ppi_distance[i,j])
    ppi_sim[j,i]=ppi_sim[i,j]
  }
}
diag(ppi_sim)=1
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)
save(ppi_sim,file="ppi_sim.Rdata")
write.table(ppi_sim,file="ppi_sim.txt",sep ="\t",quote = F,row.names = T,col.names = T)
write.csv(ppi_sim,file = "ppi_sim.csv")
save.image(file = "workspace.RData")
