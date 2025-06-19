In this section, we will still use the human infectious dataset as an example 
to show how we process the results getting from different tools. 

## Set up the code before we process the result from each tool
```
# load libraries 
library(dplyr)
library(readxl)
library(funtimes)
library(aricode)
library(ggplot2)
f <- "human_infectious.csv"
```
## GIANA
```
# read the results
d_re <- read.table(paste0("../out/GIANA/", f), sep="\t")
# read the input file 
d_input <- read.table(paste0("../data/GIANA/", f), sep="\t")

d <- d_re %>% left_join(d_inout, by = c("V1" = "Beta_cdr3"))
pct <- nrow(d) / nrow(d_input)
prty <- purity(d$Epitope, d$V2)$pur
ari <- ARI(d$Epitope, d$V2)
nmi <- NMI(d$Epitope, d$V2)

giana_measure <- list(pct, prty, ari, nmi)
```
## iSMART
```
# read the results
d_re <- read.table(paste0("../out/iSMART/", f), sep="\t", header=1)
# read the input file 
d_input <- read.table(paste0("../data/GIANA/", f), sep="\t")

d <- d_re %>% left_join(d_input, by = c("aminoAcid" = "Beta_cdr3"))
pct <- nrow(d) / nrow(d_input)
d <- d[!is.na(d$Epitope), ]
prty <- purity(d$Epitope, d$Group)$pur
ari <- ARI(d$Epitope, d$Group)
nmi <- NMI(d$Epitope, d$Group)

ismart_measure <- list(pct, prty, ari, nmi)
```
## clustTCR
```
# read the results
d_re <- read.table(paste0("../out/clusTCR/", f), sep="\t", header=1)
# read the input file 
d_input <- read.table(paste0("../data/GIANA/", f), sep="\t")

d <- d_re %>% left_join(d_input, by = c("aminoAcid" = "Beta_cdr3"))
pct <- nrow(d) / nrow(d_input)
d <- d[!is.na(d$Epitope), ]
prty <- purity(d$Epitope, d$Group)$pur
ari <- ARI(d$Epitope, d$Group)
nmi <- NMI(d$Epitope, d$Group)

clustcr_measure <- list(pct, prty, ari, nmi)
```
## TCRex
```
# read the results
d_re <- read.table(paste0("../out/cTCRex/", f), sep="\t", header=1)
# read the input file 
d_input <- read.table(paste0("../data/TCRex/", f), sep="\t")

d <- d_re %>% group_by(CDR3_beta) %>% arrange(desc(score)) %>% slice(1)
d <- d %>% left_join(d_input, by = c("CDR3_beta" = "Beta_cdr3"))
  
pct <- nrow(d) / nrow(d_input)
d <- d[!is.na(d$Epitope), ]
prty <- purity(d$Epitope, d$epitope)$pur
ari <- ARI(d$Epitope, d$epitope)
nmi <- NMI(d$Epitope, d$epitope)

tcrex_measure <- list(pct, prty, ari, nmi)
```
## GLIPH2
```
# read the results
d_re <- read.csv(paste0("../out/GLIPH2/", f))
# read the input file 
d_input <- read.table(paste0("../data/GLIPH2/", f), sep="\t")

d <- d_re %>% left_join(d_input, by = c("TcRb" = "Beta_cdr3"))
pct <- nrow(d) / nrow(d_input)
d <- d[!is.na(d$Epitope), ]
prty <- purity(d$Epitope, d$pattern)$pur
ari <- ARI(d$Epitope, d$pattern)
nmi <- NMI(d$Epitope, d$pattern)

gliph2_measure <- list(pct, prty, ari, nmi)
```
## TCR-BERT
```
# read the results
d_re <- read.csv(paste0("../out/TCRbert/", x), header=FALSE)
# read the input file 
d_input <- read.table(paste0("../data/TCRbert/", f), sep="\t")

splitList <- split(d_re, 1:nrow(d_re))
splitList <- lapply(seq(splitList), function(i) data.frame(cdr3=unlist(splitList[[i]]), group=i))
re_bert <- do.call(rbind, splitList)
re_bert <- re_bert[re_bert$cdr3 != "", ]
d <- re_bert %>% left_join(d_input, by = c("cdr3" = "Beta_cdr3"))
pct <- nrow(d) / nrow(input_f[[i]])
d <- d[!is.na(d$Epitope), ]
prty <- purity(d$Epitope, d$group)$pur
ari <- ARI(d$Epitope, d$group)
nmi <- NMI(d$Epitope, d$group)

tcrbert_measure <- list(pct, prty, ari, nmi)
```
## pMTnet
```
# read the results
d_re <- read.csv(paste0("../out/pMTnet/", x))
# read the input file 
d_input <- read.csv(paste0("../data/pMTnet/", x))

d <- d_re %>% group_by(CDR3) %>% arrange(Rank) %>% slice(1)
d <- d %>% left_join(d_input, by = c("CDR3" = "Beta_cdr3"))

pct <- nrow(d) / nrow(d_input)
d <- d[!is.na(d$Epitope), ]
prty <- purity(d$Epitope, d$Antigen)$pur
ari <- ARI(d$Epitope, d$Antigen)
nmi <- NMI(d$Epitope, d$Antigen)

pmtnet_measure <- list(pct, prty, ari, nmi)
```
## SETE
```
# read the results
d_re <- read.csv(paste0("../out/SETE/", x))
# read the input file 
d_input <- read.csv(paste0("../data/SETE/", x))

d <- d_re
pct <- nrow(d) / nrow(d_input)
prty <- purity(d$true, d$pred)$pur
ari <- ARI(d$true, d$pred)
nmi <- NMI(d$true, d$pred)
 
sete_measure <- list(pct, prty, ari, nmi)
```
## NetTCR2.0
```
# read the results
d_re <- read.csv(paste0("../out/NetTCR/", x))
# read the input file 
d_input <- read.csv(paste0("../data/NetTCR/", x))

d_re$group <- paste0(d_re$CDR3a, "+", re_nettcr$CDR3b)
d <- re_nettcr %>% group_by(group) %>% arrange(desc(prediction)) %>% slice(1)
  
d_input <- d_input %>% mutate(CDR3a=sub(".", "", data_f[[i]]$Alpha_cdr3)) %>% 
                                          mutate(CDR3b=sub(".", "", data_f[[i]]$Beta_cdr3))
d <- d %>% left_join(d_input)
  
pct <- nrow(d) / nrow(d_input)
d <- d[!is.na(d$Epitope), ]
prty <- purity(d$Epitope, d$peptide)$pur
ari <- ARI(d$Epitope, d$peptide)
nmi <- NMI(d$Epitope, d$peptide)

nettcr_measure <- list(pct, prty, ari, nmi)
```
## tcrdist3
```
# read the results
d_re <- read.table(paste0("../out/tcrdist3/", x), sep = "\t", header=1)
# read the input file 
d_input <- read.csv(paste0("../data/tcrdist3/", x))

d <- d_re
pct <- nrow(d) / nrow(d_input)
prty <- purity(d$true, d$pred)$pur
ari <- ARI(d$true, d$pred)
nmi <- NMI(d$true, d$pred)
 
tcrdist3_measure <- list(pct, prty, ari, nmi)
```
## TCRmatch
```
# read the results
d_re <- read.table(paste0("../out/cTCRex/", f), sep="\t", header=1)
# read the input file 
d_input <- read.table(paste0("../data/TCRex/", f), sep="\t")

pct <- nrow(d) / nrow(d_input)
cdr3 <-lapply(
    d$trimmed_input_sequence,
    function(x) d_input$Beta_cdr3[grepl(x, d_input$Beta_cdr3, fixed = TRUE)]
    )
cdr3 <- lapply(cdr3, function(x) x[1])
d$cdr3 <- unlist(cdr3)
d <- d %>% left_join(d_input, by = c("cdr3" = "Beta_cdr3"))
prty <- purity(d$Epitope, d$epitope)$pur
ari <- ARI(d$Epitope, d$epitope)
nmi <- NMI(d$Epitope, d$epitope)

tcrmatch_measure <- list(pct, prty, ari, nmi)
```
## Create the plot 
```
# combine the results 
measure_l <- list(
  tcrex_measure, ismart_measure, gliph2_measure, sete_measure, clustcr_measure, 
  giana_measure, nettcr_measure, pmtnet_measure, tcrbert_measure, tcrmatch_measure, 
  tcrdist3_measure
)
names(measure_l) <- c("TCRex", "iSMART", "GLIPH2", "SETE", "clusTCR", "GIANA", 
                      "NetTCR", "pMTnet", "TCR-BER", "TCRMatch", "tcrdist3")
measure_l <- lapply(measure_l, function(x) {names(x) <- c("Percentage", "Purity", "ARI", "NMI"); x})
measure_l <- lapply(measure_l, function(x) as.data.frame(do.call(cbind, x)))
measure_l <- lapply(seq(measure_l), function(x) {measure_l[[x]]$Method <- names(measure_l)[x]; measure_l[[x]]})

measure_d <- do.call(rbind, measure_l)
measure_d <- measure_d[measure_d$NMI != "NaN", ]
measure_d$`Percentage` <- as.numeric(measure_d$`Percentage`)
measure_d$Purity <- as.numeric(measure_d$Purity)
measure_d$ARI <- as.numeric(measure_d$ARI)
measure_d$NMI <- as.numeric(measure_d$NMI)

# set the colors 
colscale <- scale_fill_manual(
  values = c( "clusTCR" = "#FED789FF", "GIANA" = "#E64B35FF", "GLIPH2" = "#4DBBD5FF", "iSMART" = "#00A087FF", 
             "NetTCR" = "#3C5488FF","pMTnet" = "#F39B7FFF", "SETE" = "#8491B4FF",
             "TCR-BERT" = "#91D1C2FF", "tcrdist3" = "#DC0000FF", "TCRex" = "#7E6148FF",
             "TCRMatch" = "#B09C85FF")
  )
# bar plot 
measure_d_long <- reshape2::melt(measure_d)
ggplot(measure_d_long) + 
	geom_col(aes(x=variable, y=value, fill=Method), position="dodge") + 
         colscale + 
         theme_classic() + 
      	 theme(axis.text.x = element_text(size=12),
	 	axis.text.y = element_text(size=12),
            	axis.title.x = element_blank(), 
            	axis.title.y = element_blank(), legend.position = "none")
	
```