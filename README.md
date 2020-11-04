# SLprediction

## “Synthetic Lethal Interactions Prediction Based on Multiple Similarity Measures Fusion”

Abstract：The synthetic lethality (SL) relationship arises when a combination of deficiencies in two genes leads to cell death, whereas a deficiency in either one of the two genes does not. The survival of the mutant tumor cells depends on the SL partner genes of the mutant gene, so the cancer cells could be selectively killed by inhibiting the SL partners of the oncogenic genes but normal cells not. Therefore, developing SL pairs identification methods is increasingly needed for cancer targeted therapy. In this paper, we proposed a new approach based on similarity fusion to predict SL pairs. Multiple types of gene similarity measures are integrated and k-NN algorithm are applied to achieve the similarity-based classification task between gene pairs. As a similarity-based method, our method demonstrated excellent performance in multiple experiments. Besides the effectiveness of our method, the ease of use and expansibility can also make our method more widely used in practice.

# Get Started

## Run Example

Gene-gene similarity measures: calculate seven gene–gene similarity measures including the similarities based on gene expression profile, gene encoded protein sequence, PPI network, co-pathway, Gene Ontology Biological Process (GOBP), Gene Ontology Cellular Component (GOCC) and Gene Ontology Molecular Function (GOMF).

SNF.R：integrate the similarity measures based on the gene expression profile, protein sequence, protein–protein interaction (PPI) network, co-pathway and Gene Ontology (GO).

knn.R：the main procedure of SL prediction. The similarity-based classification task between gene pairs.

ROCandPR：Draw ROC and PR curves.
