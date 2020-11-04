setwd("/home/SL/knn")
.libPaths("/home/SL/Rpackage/")
library(caret) 
set.seed(777)

sim = read.csv(file = "affini_fusion.csv", row.names = 1)
colnames(sim) = rownames(sim)
diag(sim) = 0
postiveset = read.csv(file = "new_postiveset.csv", header = F)
# edge_topn = edge_fusion_all[order(edge_fusion_all[,5], decreasing = T),][1:n,]  
# postiveset = postiveset[order(postiveset[,3], decreasing = T),]
negtiveset = read.csv(file = "new_negtiveset.csv", header = F)
randomnum = sample(1:nrow(negtiveset), nrow(postiveset))
negtiveset = negtiveset[randomnum,]
postiveset$label = rep(1, nrow(postiveset))
negtiveset$label = rep(0, nrow(negtiveset))
alldataset = rbind(postiveset, negtiveset)
write.csv(alldataset, file = "alldataset.csv")

pair_sim <- function(drug_d1, drug_d2, drugp_d1, drugp_d2)
{
  f1 = sqrt(sim[drug_d1,drugp_d1]*sim[drug_d2,drugp_d2])
  f2 = sqrt(sim[drug_d1,drugp_d2]*sim[drug_d2,drugp_d1])
  f = max(f1,f2)
  return(f)
}

folds = createFolds(y = as.factor(alldataset$label), k = 10)
save(folds, file = "folds_affini_fusion.Rdata")
for (x in 1:10) {
  traindata = alldataset[-folds[[x]],]
  testdata = alldataset[folds[[x]],]
  testdata$score = rep(0, nrow(testdata))
  scoreall_all = data.frame()
  for (i in 1:nrow(testdata)) {
    scoreall = c()
    for (j in 1:nrow(traindata)) {
      score = pair_sim(as.character(testdata[i,1]),as.character(testdata[i,2]),as.character(traindata[j,1]),as.character(traindata[j,2]))
      scoreall = c(scoreall, score)
      print(paste("x=",x,",testi=",i,",testj=",j,sep = ""))
    }
    scoreall_all = rbind(scoreall_all, scoreall)
    index = order(scoreall, decreasing=TRUE)[1:10]
    testdata[i,4] = sum(traindata[index,3])
  }
  write.csv(scoreall_all, file = paste("scoreall", x, ".csv", sep = ""))
  write.csv(testdata, file = paste("testdata", x, ".csv", sep = ""))
  write.csv(traindata, file = paste("traindata", x, ".csv", sep = ""))
}






