# Database
## IEDB datasets
All the IEDB datasets were downloaded from https://www.iedb.org/. The parameters for each datasets are as following. 
 **Human infectious dataset** 
 ```
 Include Positive Assays
 No B cell assays
 No MHC assays
 Host: Homo sapiens (human)
 Disease Data: Infectious Disease
 Export Type: Full, all data columns
 ```
 **Human allergic dataset** 
```
Include Positive Assays
No B cell assays
No MHC assays
Host: Homo sapiens (human)
Disease Data: Allergic Disease
Export Type: Full, all data columns
```
**Human autoimmune dataset**
```
Include Positive Assays
No B cell assays
No MHC assays
Host: Homo sapiens (human)
Disease Data: Autoimmune Disease
Export Type: Full, all data columns
```
**Mouse dataset**
```
Include Positive Assays
No B cell assays
No MHC assays
MHC Restriction Type: Class I
Export Type: Full, all data columns
```
**IEDB data processing**
```
f <- list.files("../iedb_datasets/")
iedb_data <- NULL
for(i in seq(f)){
  d <- openxlsx::read.xlsx(paste0("../datasets/", f[i]), colNames=FALSE)
  # column 24: cdr3 alpha curated 
  d <- d[!is.na(d[, 24]), ]
  # column 53: cdr3 beta curated 
  d <- d[!is.na(d[, 53]), ]
  
  d_f <- d[, c(2, 7:9, 12, 16, 20, 24, 30, 36, 45, 49, 53, 59, 65)]
  names(d_f) <- c("IEDB_id", "Epitope", "Source_molecule", "Source_organism", "MHC",
                  "Alpha_Vgene", "Alpha_Jgene", "Alpha_cdr3", "Alpha_cdr1", "Alpha_cdr2",
                  "Beta_Vgene", "Beta_Jgene", "Beta_cdr3", "Beta_cdr1", "Beta_cdr2")
  d_f <- d_f[-c(1:2), ]
  d_f <- d_f[!is.na(d_f$Beta_Vgene), ]
  d_f <- d_f[!is.na(d_f$Beta_Jgene), ]
  d_f <- d_f[!is.na(d_f$Alpha_Vgene), ]
  d_f <- d_f[!is.na(d_f$Alpha_Jgene), ]
  d_f <- d_f[!duplicated(d_f$Beta_cdr3), ]
  iedb_data[[i]] <- d_f
}
```
## McPAS dataset (cancer VS. non-cancer)
The dataset was downloaded from https://friedmanlab.weizmann.ac.il/McPAS-TCR/. 

**McPAS data processing**
```
d <- read.csv("../datasets/McPAS-TCR.csv")
d <- d[d$Species == "Human", ]
d <- d[!is.na(d$CDR3.alpha.aa), ]
d <- d[!is.na(d$CDR3.beta.aa), ]
d <- d[!is.na(d$TRAJ), ]
d <- d[!is.na(d$TRAV), ]
d <- d[!is.na(d$TRBJ), ]
d <- d[!is.na(d$TRBV), ]
d <- d[!duplicated(d$CDR3.beta.aa), ]
d_cancer <- d[d$Category == "Cancer", ]
d_noncancer <- d[d$Category != "Cancer", ]
```
## VDJdb dataset (Human influenza, SARS-cov-2, CMV, and mouse)
The data set was downloaded from https://vdjdb.cdr3.net/.

**VDJdb data processing**
```
# vdjdb
d <- read.csv("../datasets/vdjdb_full.txt", sep="\t")
d <- d[d$cdr3.alpha != "", ]
d <- d[d$cdr3.beta != "", ]
d <- d[d$v.alpha != "", ]
d <- d[d$j.alpha != "", ]
d <- d[d$v.beta != "", ]
d <- d[d$j.beta != "", ]
d_f <- d[d$species == "HomoSapiens", ]
vdj_human_data <- list()
vdj_human_data[["Influenza"]] <- d_f[d_f$antigen.species %like% "Influenza", ]
vdj_human_data[["SARS-cov-2"]] <- d_f[d_f$antigen.species %like% "SARS-CoV-2", ]
vdj_human_data[["CMV"]] <- d_f[d_f$antigen.species == "CMV", ]
vdj_mouse_data <- d[d$species == "MusMusculus", ]
```

# Real data 
## 10x nsclc
```
d <- read.csv("../datasets/vdj_v1_hs_nsclc_t_all_contig_annotations.csv")
d <- d[d$cdr3 != "None", ]
d_alpha <- d[d$chain == "TRA", ]
d_beta <- d[d$chain == "TRB", ]
d_new <- d_beta %>% left_join(d_alpha, join_by(barcode))
d_new <- d_new[!duplicated(d_new$cdr3.x), ]
human_nsclc <- na.omit(d_new)
```
## Su covid-19 
```
d <- read.csv("../datasets/su_TCR_cd8_AllPatients.csv")
d <- d[d$TRB_cdr3 != "None", ]
d <- d[d$TRA_cdr3 != "None", ]
d <- d[!duplicated(d$TRB_cdr3), ]
d <- d[!is.na(d$patient_subID), ]
set.seed(2333)
d <- d[sample(nrow(d), 10000), ]
su_covid19 <- d
```
## Zheng pancancer
```
d <- read.csv("../datasets/zheng_pancancer.txt", sep="\t")
d <- d[d$cdr3 != "None", ]
d <- d[d$v_gene != "None", ]
d_alpha <- d[d$chain == "TRA",]
d_beta <- d[d$chain == "TRB", ]
d_new <- na.omit(d_beta %>% left_join(d_alpha, join_by(barcode)))
d_new <- d_new[!duplicated(d_new$cdr3.x), ]
set.seed(2333)
d_new <- d_new[sample(nrow(d_new), 10000), ]
zheng_pancancer <- d_new
```
## ctcl
```
d <- read.csv("../datasets/CTCL-TCRab.csv")
d <- d[d$cdr3 != "None", ]
d_alpha <- d[d$chain == "TRA", ]
d_beta <- d[d$chain == "TRB", ]
d_new <- na.omit(d_beta %>% left_join(d_alpha, join_by(barcode)))
d_new <- d_new[!duplicated(d_new$cdr3.x), ]
mimitou_ctcl <- d_new
```
## Breast
```
d <- read.csv("../datasets/azizi_breasttumor.csv")
d <- d[d$cdr3 != "None", ]
d_alpha <- d[d$chain == "TRA", ]
d_beta <- d[d$chain == "TRB", ]
d_new <- na.omit(d_beta %>% left_join(d_alpha, join_by(barcode)))
d_new <- d_new[!duplicated(d_new$cdr3.x), ]
azizi_breast <- d_new
```
## Longitudianl melanoma
```

```

