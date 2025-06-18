# Tool application

In this section, we will use the human infectious IEDB dataset as an example showing how we apply each tool on the data. 

## GIANA and ClusTCR
These two methods share the same input files. GIANA can be installed from <https://github.com/s175573/GIANA>. 
iSMART can be installed from <https://github.com/s175573/iSMART>. 
```
# read the processed IEDB file
tmp <- readRDS("path_to_IEDB_data")
d_giana <- data.frame(aminoAcid=tmp$Beta_cdr3, vGeneName=tmp$Beta_Vgene, `count..templates.reads.`=1) 
# save the data for GIANA
write.table(d_giana, paste0("../data/GIANA/", "human_infectious", ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
```
Then go to the installed GIANA folder and run the following code in the terminal
```
python /GIANA4.py -f ../data/GIANA/chuman_infectious.tsv -o ../out/GIANA/
```
For iSMART, go to the installed iSMART folder and run the following code in the terminal.
```
python iSMARTv3.py -f ../data/GIANA/human_infectious.tsv -o ../out/iSMART/
```

## GLIPH2
```
# read the processed IEDB file
tmp <- readRDS("path_to_IEDB_data")
d_gliph <- data.frame(CDR3b=tmp$Beta_cdr3, TRBV=tmp$Beta_Vgene, TRBJ=tmp$Beta_Jgene, CDR3a=tmp$Alpha_cdr3,
                        `subject:condition`=NA, count=1) %>% 
    distinct(CDR3b, .keep_all = TRUE)
write.table(d_gliph, paste0("../data/GLIPH2/", "human_infectious", ".tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
```
We run GLIPH2 through its website at <http://50.255.35.37:8080/>. 

## TCRex
```
# read the processed IEDB file
tmp <- readRDS("path_to_IEDB_data")
d_tcrex <- data.frame(CDR3_beta=tmp$Beta_cdr3,	TRBJ_gene=tmp$Beta_Jgene,	TRBV_gene=tmp$Beta_Vgene) 
write.table(d_tcrex, paste0("../data/TCRex/", "human_infectious", ".tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
```
We run TCRex through its website at <https://tcrex.biodatamining.be/>.

## TCRbase
```
# read the processed IEDB file
tmp <- readRDS("path_to_IEDB_data")
d_tcrbase <- data.frame(cdr1a=tmp$Alpha_cdr1, cdr2a=tmp$Alpha_cdr2, cdr3a=tmp$Alpha_cdr3,
                          cdr1b=tmp$Beta_cdr1, cdr2b=tmp$Beta_cdr2, cdr3b=tmp$Beta_cdr3) %>% 
    distinct(cdr3b, .keep_all = TRUE)
d_tcrbase$cdr1a = "XXXXX"
d_tcrbase$cdr2a = "XXXXX"
d_tcrbase$cdr1b = "XXXXX"
d_tcrbase$cdr2b = "XXXXX"
write.table(d_tcrbase, paste0("../data/TCRbase/", "human_infectious", ".tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
```
We run TCRbase through its website at <https://services.healthtech.dtu.dk/services/TCRbase-1.0/>. 

## TCR-BERT
TCR-BERT can be installed from <https://github.com/wukevin/tcr-bert>.
```
# read the processed IEDB file
tmp <- readRDS("path_to_IEDB_data")
d_tcrbert <- data.frame(cdr3b=tmp$Beta_cdr3) %>% 
    distinct(cdr3b, .keep_all = TRUE)
d_tcrbert$group <- 1
write.table(d_tcrbert, paste0("../data/TCRbert/", "human_infectious", ".tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
```
Then go to the installed TCR-BERT folder, and run the following code in the terminal.
```
python bin/embed_and_cluster.py ../data/TCRbert/human_infectious.tsv ../out/TCRbert/human_infectious.tsv -r 32 -g 0

```
## TCRMatch
TCRMatch can be installed from <https://github.com/IEDB/TCRMatch>.
```
# read the processed IEDB file
tmp <- readRDS("path_to_IEDB_data")
d_tcrmatch <- unique(tmp$Beta_cdr3)
write.table(d_tcrmatch, paste0("../data/TCRMatch/", "human_infectious", ".txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
``` 
Then go to the installed TCRMatch folder, and run the following code in the terminal. 
```
```
## pMTnet
```
# read the processed IEDB file
tmp <- readRDS("path_to_IEDB_data")
d_pmtnet_v1 <- data.frame(CDR3=tmp$Beta_cdr3, Antigen=tmp$Epitope, HLA=tmp$MHC) %>% 
    distinct(CDR3, .keep_all = TRUE)

d_pmtnet_v1 <- d_pmtnet_v1[!d_pmtnet_v1$HLA %like% "human", ]
d_pmtnet_v1$HLA <- gsub("\\,.*", "", d_pmtnet_v1$HLA)
d_pmtnet_v1$HLA <- gsub("HLA-", "", d_pmtnet_v1$HLA)
d_pmtnet_v1 <- d_pmtnet_v1[!d_pmtnet_v1$HLA %like% " ", ]
d_pmtnet_v1 <- d_pmtnet_v1[!d_pmtnet_v1$HLA %like% "\\/", ]
d_pmtnet_v1 <- d_pmtnet_v1 %>% arrange(desc(HLA))
d_pmtnet_v1 <- na.omit(d_pmtnet_v1)
  
table <- as.data.frame(table(d_pmtnet_v1$Antigen)) %>% arrange(desc(Freq))
d_pmtnet_v1 <- d_pmtnet_v1[d_pmtnet_v1$Antigen %in% table$Var1[1:10], ]
out <- d_pmtnet_v1[rep(seq_len(nrow(d_pmtnet_v1)), each = 10), ]
out$Antigen <- rep(table$Var1[1:10], nrow(d_pmtnet_v1))
write.csv(out, paste0("../data/pMTnet/", "human_infectious", ".csv"), 
            row.names = FALSE, quote = FALSE)
```
We run the pMTnet through its website at <https://dbai.biohpc.swmed.edu/pmtnet/analysis-base.php>.

## SETE
SETE can be installed from <https://github.com/wonanut/SETE>. 
```
# read the processed IEDB file
tmp <- readRDS("path_to_IEDB_data")
d_sete <- na.omit(data.frame(epitope=tmp$Epitope, vb_gene=tmp$Beta_Vgene, cdr3b=tmp$Beta_cdr3)) 
```
Then go to the installed SETE folder, we can run the following Python code. 
```
import numpy as np
import pandas as pd
import collections
import matplotlib.pyplot as plt

from SETE import *
from itertools import cycle
from sklearn.metrics import roc_curve, auc, accuracy_score, precision_score, recall_score
from sklearn.model_selection import train_test_split, KFold, StratifiedKFold
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier, OneVsOneClassifier
from scipy import interp
from sklearn.base import clone
from sklearn.decomposition import PCA, KernelPCA
from sklearn.cross_decomposition import CCA

from sklearn.feature_selection import SelectKBest, chi2
from sklearn.ensemble import GradientBoostingClassifier

import warnings

warnings.filterwarnings('ignore')

f = os.listdir("../data/SETE")

for file in f:
    X, y, epiname_list = data_preprocess("../data/SETE/" + file, 3)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.35, random_state=666)
    classifier = GradientBoostingClassifier(learning_rate=0.1,
                                        min_samples_leaf=20, max_features='sqrt', subsample=0.8,
                                        random_state=666, n_estimators=70, max_depth=11,
                                        min_samples_split=60, loss="log_loss"
                                        )

    X_train, X_test = pca_analyse(X_train, X_test, 0.9)
    classifier.fit(X_train, y_train)
    classifier.score(X_test, y_test)
    pred = classifier.predict(X_test)
    result = pd.DataFrame(np.array([pred, y_test]).T, columns=["pred", "true"])
    result.to_csv("../out/SETE/" + file)
```
## NetTCR2.0
The NetTCR can be installed from <https://github.com/mnielLab/NetTCR-2.0>. 
```
# read the processed IEDB file
tmp <- readRDS("path_to_IEDB_data")
d_tmp <- data.frame(CDR3a=sub(".", "", tmp$Alpha_cdr3),
                      CDR3b=sub(".", "", tmp$Beta_cdr3), 
                      peptide=tmp$Epitope)
d_tmp <- d_tmp[!grepl(" ", d_tmp$peptide), ]
d_tmp <- d_tmp[!grepl("-", d_tmp$peptide), ]
d_tmp <- d_tmp[!grepl("#", d_tmp$CDR3a), ]
d_tmp <- d_tmp[d_tmp$peptide != "carbamazepine", ]
d_tmp <- d_tmp[!grepl("X", d_tmp$CDR3a), ]
d_tmp <- d_tmp[!grepl("\\*", d_tmp$CDR3a), ]
d_tmp <- na.omit(d_tmp)
n_train <- nrow(d_tmp) * 0.65
idx <- sample(seq(nrow(d_tmp)), round(n_train))
d_train <- d_tmp[idx, ]
d_train$binder <- 0
d_test <- d_tmp[-idx, ]
d_test <- d_test[nchar(d_test$peptide) <= 9,]
d_train <- d_train[d_train$peptide %in% d_test$peptide, ]
d_test <- d_test[d_test$peptide %in% d_train$peptide, ]

write.table(d_train, paste0("../data/NetTCR/train_", "human_infectious", ".csv"), 
              sep = ",", row.names = FALSE, quote = FALSE)
write.table(d_test, paste0("../data/NetTCR/test_", "human_infectious", ".csv"), 
              sep = ",", row.names = FALSE, quote = FALSE)
```
Then we can go to the installed NetTCR folder and run the following code in the terminal. 
```
python NetTCR-2.0/nettcr.py --trainfile ../data/NetTCR/train_human_infectious.csv --test ../data_new/NetTCR/test_human_infectious.csv
 -o ../out/NetTCR/human_infectious.csv
```
## tcrdist3
Please refer to this website <https://tcrdist3.readthedocs.io/en/latest/> for tcrdist3. 
```
# read the processed IEDB file
d_f <- readRDS("path_to_IEDB_data")
tmp <- data.frame(subject=d_f$IEDB_id, epitope=d_f$Epitope, count=1, 
                                v_a_gene=d_f$Alpha_Vgene, j_a_gene=d_f$Alpha_Jgene,
                                cdr3_a_aa=d_f$Alpha_cdr3, v_b_gene=d_f$Beta_Vgene, 
                                j_b_gene=d_f$Beta_Jgene, cdr3_b_aa=d_f$Beta_cdr3) %>%
    distinct(cdr3_b_aa, .keep_all = TRUE) %>%
    distinct(cdr3_a_aa, .keep_all = TRUE)
tmp <- tmp[nchar(tmp$epitope) <= 9,]
d_tcrdist3 <- tmp
d_tcrdist3$clone_id <- seq(nrow(d_tcrdist3))
write.table(d_tcrdist3, paste0("../data/tcrdist3/", "human_infectious", ".tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
```
Then we can run the following python code. 
```
import pandas as pd
import numpy as np
import os
from tcrdist.repertoire import TCRrep
from tcrdist.tree import TCRtree
import pwseqdist as pw
from tcrdist.rep_diff import neighborhood_diff
from sklearn.cluster import AgglomerativeClustering

path = "../../data/tcrdist3/"
f = os.listdir(path)

for file in f:
    df = pd.read_csv(path + file, sep="\t")
    tr = TCRrep(cell_df = df, 
            organism = 'mouse', 
            chains = ['beta'], 
            db_file = 'alphabeta_gammadelta_db.tsv')
    distance = tr.pw_beta
    model = AgglomerativeClustering(affinity='precomputed', n_clusters=10, linkage='complete').fit(distance)
    tr.clone_df["labels"] = model.labels_
    tr.clone_df.to_csv("../out/tcrdist3/" + file, sep="\t")
```

