setwd("F:/R workspace/ROCandPR/")
library(ROCR)

#### roc ####
# kindex = 19
allres = data.frame()
for(k in 1:10){
  postiveset = read.csv(file = paste("F:/R workspace/KNN/results/fusion/fusion_testdata",k,".csv",sep=""), header = T, row.names = 1)
  postiveset = postiveset[order(postiveset[,3], decreasing = T),]

  # traindata = read.csv(file = paste("F:/R workspace/KNN/results/fusion/fusion_traindata",k,".csv",sep=""), header = T, row.names = 1)
  # scoreall = read.csv(file = paste("F:/R workspace/KNN/results/fusion/fusion_scoreall",k,".csv",sep=""), header = T, row.names = 1)
  # 
  # for (i in 1:nrow(scoreall)) {
  #   scoretemp = scoreall[i,]
  #   index = order(scoretemp, decreasing=TRUE)[1:kindex]
  #   postiveset[i,5] = sum(traindata[index,3])
  #   print(paste("k",k,",i=",i,sep = ""))
  # }
  allres = rbind(allres, postiveset)
}
pos_adjall = allres
pos_adjall$predicted=1*(pos_adjall$score>9)
x = as.numeric(pos_adjall$score)
center <- x - min(as.numeric(pos_adjall$score))
R <- max(as.numeric(pos_adjall$score)) - min(as.numeric(pos_adjall$score))
newscore <- center/R
pos_adjall$newscore = newscore
pred <-prediction(pos_adjall$newscore,pos_adjall$label)
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='#377EB8', 
     # main=paste("ROC curve (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)


#### pr ####
allres = data.frame()
for(k in 1:10){
  postiveset = read.csv(file = paste("F:/R workspace/KNN/results/fusion/fusion_testdata",k,".csv",sep=""), header = T, row.names = 1)
  postiveset = postiveset[order(postiveset[,3], decreasing = T),]
  
  # traindata = read.csv(file = paste("F:/R workspace/KNN/results/fusion/fusion_traindata",k,".csv",sep=""), header = T, row.names = 1)
  # scoreall = read.csv(file = paste("F:/R workspace/KNN/results/fusion/fusion_scoreall",k,".csv",sep=""), header = T, row.names = 1)
  # 
  # for (i in 1:nrow(scoreall)) {
  #   scoretemp = scoreall[i,]
  #   index = order(scoretemp, decreasing=TRUE)[1:kindex]
  #   postiveset[i,5] = sum(traindata[index,3])
  #   print(paste("k",k,",i=",i,sep = ""))
  # }
  allres = rbind(allres, postiveset)
}
pos_adjall = allres
pos_adjall$predicted=1*(pos_adjall$score>9)
x = as.numeric(pos_adjall$score)
center <- x - min(as.numeric(pos_adjall$score))
R <- max(as.numeric(pos_adjall$score)) - min(as.numeric(pos_adjall$score)) 
newscore <- center/R
pos_adjall$newscore = newscore
pred2 <-prediction(pos_adjall$newscore,pos_adjall$label)
perf2 <- performance(pred2, "prec", "rec")
plot(perf2,
     xlim=c(0,1), ylim=c(0,1),col='#377EB8',
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)

