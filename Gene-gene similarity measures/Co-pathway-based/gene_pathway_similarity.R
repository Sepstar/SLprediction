setwd("F:/R workspace/Gene-gene similarity measures/Co-pathway-based/")

#calculate gene pathway similarity
gene_pathway_matrix = read.csv(file = "gene_pathway_matrix.csv", row.names = 1)
gene_id = unique(rownames(gene_pathway_matrix))
gene_gene_s=matrix(0,length(gene_id),length(gene_id))
rownames(gene_gene_s)=colnames(gene_gene_s)=gene_id
for(i in 1:nrow(gene_gene_s))
{
  print(i)    
  gene_1=gene_pathway_matrix[i,]
  gene_1[which(gene_1==0)]=NA
  for(j in i:nrow(gene_gene_s))
  {
    print(paste("j =",j))
    gene_2=gene_pathway_matrix[j,]
    gene_2[which(gene_2==0)]=NA
    a=length(which(gene_1==gene_2))
    b=length(which(gene_1==1))+length(which(gene_2==1))-a
    gene_gene_s[i,j]=a/b
    gene_gene_s[j,i]=gene_gene_s[i,j]
  }
}
write.csv(gene_gene_s, file = "gene_gene_s.csv")
