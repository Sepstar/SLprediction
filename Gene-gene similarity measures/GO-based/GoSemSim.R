setwd("F:/R workspace/Gene-gene similarity measures/GO-based/")
library(GOSemSim)

entrezid = read.csv(file = "Entrez_ID.csv",header = F)
entrezid = unlist(entrezid)
d.CC = godata('org.Hs.eg.db', ont="CC", computeIC=FALSE)
d.BP = godata('org.Hs.eg.db', ont="BP", computeIC=FALSE)
d.MF = godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)

geneBPsim = mgeneSim(entrezid, semData=d.BP, measure="Wang")
write.csv(geneBPsim, file = "BPsim.csv")
geneCCsim = mgeneSim(entrezid, semData=d.CC, measure="Wang")
write.csv(geneCCsim, file = "CCsim.csv")
geneMFsim = mgeneSim(entrezid, semData=d.MF, measure="Wang")
write.csv(geneMFsim, file = "MFsim.csv")

