setwd("F:/R workspace/SNF/")
library(SNFtool)

BPsim = as.matrix(read.csv(file = "new_BPsim.csv", header = T, row.names = 1))
CCsim = as.matrix(read.csv(file = "new_CCsim.csv", header = T, row.names = 1))
MFsim = as.matrix(read.csv(file = "new_MFsim.csv", header = T, row.names = 1))
profilesim = as.matrix(read.csv(file = "new_profilesim.csv", header = T, row.names = 1))
pathwaysim = as.matrix(read.csv(file = "new_pathwaysim.csv", header = T, row.names = 1))
ppisim = as.matrix(read.csv(file = "new_ppisim.csv", header = T, row.names = 1))
seqsim = as.matrix(read.csv(file = "new_seqsim.csv", header = T, row.names = 1))

colnames(BPsim) = rownames(BPsim)
colnames(CCsim) = rownames(CCsim)
colnames(MFsim) = rownames(MFsim)
colnames(profilesim) = rownames(profilesim)
colnames(pathwaysim) = rownames(pathwaysim)
colnames(ppisim) = rownames(ppisim)
colnames(seqsim) = rownames(seqsim)
diag(BPsim) = 1
diag(CCsim) = 1
diag(MFsim) = 1
diag(profilesim) = 1
diag(pathwaysim) = 1
diag(ppisim) = 1
diag(seqsim) = 1

K = 20; # number of neighbors, usually (10~30)
alpha = 0.5; # hyperparameter, usually (0.3~0.8)
T = 20; # Number of Iterations, usually (10~20)
affni_BP = affinityMatrix(1-BPsim, K, alpha)
affni_CC = affinityMatrix(1-CCsim, K, alpha)
affni_MF = affinityMatrix(1-MFsim, K, alpha)
affni_profile = affinityMatrix(1-profilesim, K, alpha)
affni_pathway = affinityMatrix(1-pathwaysim, K, alpha)
affni_ppi = affinityMatrix(1-ppisim, K, alpha)
affni_seq = affinityMatrix(1-seqsim, K, alpha )

affni_fusion = SNF(list(affni_BP,affni_CC,affni_MF, affni_profile, affni_pathway, affni_ppi, affni_seq), K, T)
rownames(affni_fusion) = rownames(BPsim)
colnames(affni_fusion) = rownames(BPsim)
diag(affni_fusion) = 1

write.csv(affni_fusion, file = "affini_fusion.csv")

