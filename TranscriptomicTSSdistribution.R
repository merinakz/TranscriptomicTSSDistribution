#!/usr/bin/env Rscript
rm(list = ls())
library(dplyr)
library(tibble)
####################################
#Read your TSSpreditor file output.
####################################
Ksg <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5'seq_061922/Quinnoutput/KsgData_Cleaned.csv", header = T)
top <- Ksg %>% filter(SuperStrand == "+")
comp <- Ksg %>% filter(SuperStrand == "-")

####################################
#Read your annotation file
#################################### 

genes <- read.delim("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/3seq/Nikhil-T4pipeline-master/NC_003028.v3.17.ncrna.genes", header = F)
genes <- genes %>% add_column(UTR="", .before = "V2") %>% 
  add_column(termUTR = "", .after = "V3")
genes <- genes[complete.cases(genes),]
for(i in 1:nrow(genes)){
  genes$UTR[i] <- genes$V2[i] - 150
  genes$termUTR[i] <- genes$V3[i] + 150
}
genes.top <- genes %>% filter(V4 == "+")
colnames(genes.top) <- c("genome", "UTR", "from", "to","termUTR", "strand", "name", "old.name", "new.name", "WP", "rfam", "other" )
genes.comp <- genes %>% filter(V4 == "-")
colnames(genes.comp) <- c("genome", "termUTR", "to", "from","UTR", "strand", "name", "old.name", "new.name", "WP", "rfam", "other" )

genes.top$UTR <- as.numeric(genes.top$UTR)
genes.top$from <- as.numeric(genes.top$from)
genes.top$to <- as.numeric(genes.top$to)
genes.top$termUTR <- as.numeric(genes.top$termUTR)
genes.top <- genes.top[complete.cases(genes.top),]
genes.comp$UTR <- as.numeric(genes.comp$UTR)
genes.comp$from <- as.numeric(genes.comp$from)
genes.comp$to <- as.numeric(genes.comp$to)
genes.comp$termUTR <- as.numeric(genes.comp$termUTR)
genes.comp <- genes.comp[complete.cases(genes.comp),]
comp$SuperPos <- as.numeric(comp$SuperPos)
top$SuperPos <- as.numeric(top$SuperPos)

for(i in 2:nrow(genes.top)){
  if(genes.top$UTR[i] < genes.top$to[i-1]){
    genes.top$UTR[i] <- genes.top$to[i-1] + 50
  }
  if(genes.top$termUTR[i-1] > genes.top$from[i]){
    genes.top$termUTR[i-1] <- genes.top$from[i] - 50
  }
}

for(i in 2:nrow(genes.comp)){
  if(genes.comp$UTR[i-1] > genes.comp$to[i]){
    genes.comp$UTR[i-1] <- genes.comp$to[i] - 50
  }
  if(genes.comp$termUTR[i] < genes.comp$from[i-1]){
    genes.comp$termUTR[i] <- genes.comp$from[i-1] + 50
  }
}

rownames(genes.top) <- 1:nrow(genes.top)
rownames(genes.comp) <- 1:nrow(genes.comp)

comp <- comp %>% add_column(TSSUTR ="", .after = "SuperPos") %>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")
top <- top %>% add_column(TSSUTR ="", .after = "SuperPos")%>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")

####################################
#Annotate your file
####################################

for(i in 1:nrow(top)){
  for(j in 1:nrow(genes.top)){
    if((top$SuperPos[i] >= genes.top$from[j]) & (top$SuperPos[i] <= genes.top$to[j])){
      top$Locus[i] <- genes.top$old.name[j]
    }
    if((top$SuperPos[i] >= genes.top$to[j]) & (top$SuperPos[i] <= genes.top$termUTR[j])){
      top$termUTR[i] <- genes.top$old.name[j]
    }
    if((top$SuperPos[i] >= genes.top$UTR[j]) & (top$SuperPos[i] <= genes.top$from[j])){
      top$TSSUTR[i] <- genes.top$old.name[j]
    }
  }
}
for(i in 1:nrow(comp)){
  for(j in 1:nrow(genes.comp)){
    if((comp$SuperPos[i] <= genes.comp$from[j]) & (comp$SuperPos[i] >= genes.comp$to[j])){
      comp$Locus[i] <- genes.comp$old.name[j]
    }
    if((comp$SuperPos[i] <= genes.comp$to[j]) & (comp$SuperPos[i] >= genes.comp$termUTR[j])){
      comp$termUTR[i] <- genes.comp$old.name[j]
    }
    if((comp$SuperPos[i] <= genes.comp$UTR[j]) & (comp$SuperPos[i] >= genes.comp$from[j])){
      comp$TSSUTR[i] <- genes.comp$old.name[j]
    }
  }
}

Ksg <- rbind(comp, top)
rm(comp,top)

Ksg <- Ksg %>% arrange(SuperPos) %>% 
  filter(!Genome == "NDCt0")
Ksggenes <- Ksg
Ksg <- Ksg[,-c(15:33)]

###################### NEW CASE 1 ################################ 
######################      No Change      ###################### 
##################################################################  
KsgNDCt60 <- Ksg %>% filter(Genome == "NDCt60" & enriched == 0)
Ksg1q <-  Ksg %>% filter(Genome == "Ksg1q_t60" & enriched == 0)
Ksg3q <- Ksg %>% filter(Genome == "Ksg3q_t60" & enriched == 0)
case1 <- rbind(KsgNDCt60,Ksg1q,Ksg3q) 
case1 <- arrange(case1, SuperPos)
rm(KsgNDCt60,Ksg1q, Ksg3q)
case <- case1
temp <- table(case$SuperPos)
case1 <- case[case$SuperPos %in% names(temp[temp==3]),]
# &
KsgNDCt60 <- Ksg %>% filter(Genome == "NDCt60" & enriched == 1)
Ksg1q <-  Ksg %>% filter(Genome == "Ksg1q_t60" & enriched == 1)
Ksg3q <- Ksg %>% filter(Genome == "Ksg3q_t60" & enriched == 1)
case2 <- rbind(KsgNDCt60,Ksg1q,Ksg3q) 
case2 <- arrange(case2,SuperPos)
rm(KsgNDCt60, Ksg1q,Ksg3q)
case <- case2
temp <- table(case$SuperPos)
case2 <- case[case$SuperPos %in% names(temp[temp==3]),]
case1 <- rbind(case1, case2)

###################### NEW CASE 2 ###########################
######################      Only NDC Enriched      ############## 
##################################################################  
KsgNDCt60 <- Ksg %>% filter(Genome == "NDCt60" & enriched == 1)
Ksg1q <-  Ksg %>% filter(Genome == "Ksg1q_t60" & enriched == 0)
Ksg3q <- Ksg %>% filter(Genome == "Ksg3q_t60" & enriched == 0)
case2 <- rbind(KsgNDCt60,Ksg1q,Ksg3q) 
case2 <- arrange(case2, SuperPos)
rm(KsgNDCt60, Ksg1q, Ksg3q)
case <- case2
temp <- table(case$SuperPos)
case2 <- case[case$SuperPos %in% names(temp[temp==3]),]

###################### NEW CASE 3 ###########################
######################      Only 1Q Enriched      ############## 
##################################################################  
KsgNDCt60 <- Ksg %>% filter(Genome == "NDCt60" & enriched == 0)
Ksg1q <-  Ksg %>% filter(Genome == "Ksg1q_t60" & enriched == 1)
Ksg3q <- Ksg %>% filter(Genome == "Ksg3q_t60" & enriched == 0)
case3 <- rbind(KsgNDCt60, Ksg1q,Ksg3q)
case3 <- arrange(case3,SuperPos)
rm(KsgNDCt60, Ksg1q, Ksg3q)
case <- case3
temp <- table(case$SuperPos)
case3 <- case[case$SuperPos %in% names(temp[temp==3]),]

###################### NEW CASE 4 ###########################
######################      Only 3Q Enriched      ############## 
##################################################################  
KsgNDCt60 <- Ksg %>% filter(Genome == "NDCt60" & enriched == 0)
Ksg1q <-  Ksg %>% filter(Genome == "Ksg1q_t60" & enriched == 0)
Ksg3q <- Ksg %>% filter(Genome == "Ksg3q_t60" & enriched == 1)
case4 <- rbind(KsgNDCt60,Ksg1q,Ksg3q) 
case4 <- arrange(case4, SuperPos)
rm(KsgNDCt60, Ksg1q, Ksg3q)
case <- case4
temp <- table(case$SuperPos)
case4 <- case[case$SuperPos %in% names(temp[temp==3]),]

###################### NEW CASE 5 ###########################
######################      Both 1Q & 3Q Enriched      ############## 
##################################################################  
KsgNDCt60 <- Ksg %>% filter(Genome == "NDCt60" & enriched == 0)
Ksg1q <-  Ksg %>% filter(Genome == "Ksg1q_t60" & enriched == 1)
Ksg3q <- Ksg %>% filter(Genome == "Ksg3q_t60" & enriched == 1)
case5 <- rbind(KsgNDCt60, Ksg1q,Ksg3q)
case5 <- arrange(case5, SuperPos)
rm(KsgNDCt60, Ksg1q,Ksg3q)
case <- case5
temp <- table(case$SuperPos)
case5 <- case[case$SuperPos %in% names(temp[temp==3]),]

################################

rownames(case1) <- 1:nrow(case1)
rownames(case2) <- 1:nrow(case2)
rownames(case3) <- 1:nrow(case3)
rownames(case4) <- 1:nrow(case4)
rownames(case5) <- 1:nrow(case5)

#################################
# Case 1 Transcriptomic Distribution
#################################
genic.case1 <- case1[((!case1$Locus == "") & (case1$termUTR == "") &(case1$TSSUTR == "") ),]
UTR.case1 <- case1[((!case1$TSSUTR == "") & (case1$Locus == "") & (case1$termUTR == "")),]
intergenic.case1 <- case1[((case1$Locus == "") & (case1$TSSUTR == "") & (case1$termUTR == "")),]
termUTR.case1 <- case1[((case1$Locus == "") & (case1$TSSUTR == "") & (!case1$termUTR == "")),]
other1 <- case1[((!case1$Locus == "") & (!case1$TSSUTR == "") & (!case1$termUTR == "")),]
other2 <- case1[((case1$Locus == "") & (!case1$TSSUTR == "") & (!case1$termUTR == "")),]
other3 <- case1[((!case1$Locus == "") & (case1$TSSUTR == "") & (!case1$termUTR == "")),]
other4 <- case1[((!case1$Locus == "") & (!case1$TSSUTR == "") & (case1$termUTR == "")),]
other.case1 <- rbind(other1, other2, other3, other4)

#################################
# Case 2 Transcriptomic Distribution
#################################
genic.case2 <- case2[((!case2$Locus == "") & (case2$termUTR == "") &(case2$TSSUTR == "") ),]
UTR.case2 <- case2[((!case2$TSSUTR == "") & (case2$Locus == "") & (case2$termUTR == "")),]
intergenic.case2 <- case2[((case2$Locus == "") & (case2$TSSUTR == "") & (case2$termUTR == "")),]
termUTR.case2 <- case2[((case2$Locus == "") & (case2$TSSUTR == "") & (!case2$termUTR == "")),]
other1 <- case2[((!case2$Locus == "") & (!case2$TSSUTR == "") & (!case2$termUTR == "")),]
other2 <- case2[((case2$Locus == "") & (!case2$TSSUTR == "") & (!case2$termUTR == "")),]
other3 <- case2[((!case2$Locus == "") & (case2$TSSUTR == "") & (!case2$termUTR == "")),]
other4 <- case2[((!case2$Locus == "") & (!case2$TSSUTR == "") & (case2$termUTR == "")),]
other.case2 <- rbind(other1, other2, other3, other4)

#################################
# Case 3 Transcriptomic Distribution
#################################
genic.case3 <- case3[((!case3$Locus == "") & (case3$termUTR == "") &(case3$TSSUTR == "") ),]
UTR.case3 <- case3[((!case3$TSSUTR == "") & (case3$Locus == "") & (case3$termUTR == "")),]
intergenic.case3 <- case3[((case3$Locus == "") & (case3$TSSUTR == "") & (case3$termUTR == "")),]
termUTR.case3 <- case3[((case3$Locus == "") & (case3$TSSUTR == "") & (!case3$termUTR == "")),]
other1 <- case3[((!case3$Locus == "") & (!case3$TSSUTR == "") & (!case3$termUTR == "")),]
other2 <- case3[((case3$Locus == "") & (!case3$TSSUTR == "") & (!case3$termUTR == "")),]
other3 <- case3[((!case3$Locus == "") & (case3$TSSUTR == "") & (!case3$termUTR == "")),]
other4 <- case3[((!case3$Locus == "") & (!case3$TSSUTR == "") & (case3$termUTR == "")),]
other.case3 <- rbind(other1, other2, other3, other4)

#################################
# Case 4 Transcriptomic Distribution
#################################
genic.case4 <- case4[((!case4$Locus == "") & (case4$termUTR == "") &(case4$TSSUTR == "") ),]
UTR.case4 <- case4[((!case4$TSSUTR == "") & (case4$Locus == "") & (case4$termUTR == "")),]
intergenic.case4 <- case4[((case4$Locus == "") & (case4$TSSUTR == "") & (case4$termUTR == "")),]
termUTR.case4 <- case4[((case4$Locus == "") & (case4$TSSUTR == "") & (!case4$termUTR == "")),]
other1 <- case4[((!case4$Locus == "") & (!case4$TSSUTR == "") & (!case4$termUTR == "")),]
other2 <- case4[((case4$Locus == "") & (!case4$TSSUTR == "") & (!case4$termUTR == "")),]
other3 <- case4[((!case4$Locus == "") & (case4$TSSUTR == "") & (!case4$termUTR == "")),]
other4 <- case4[((!case4$Locus == "") & (!case4$TSSUTR == "") & (case4$termUTR == "")),]
other.case4 <- rbind(other1, other2, other3, other4)

#################################
# Case 5 Transcriptomic Distribution
#################################
genic.case5 <- case5[((!case5$Locus == "") & (case5$termUTR == "") &(case5$TSSUTR == "") ),]
UTR.case5 <- case5[((!case5$TSSUTR == "") & (case5$Locus == "") & (case5$termUTR == "")),]
intergenic.case5 <- case5[((case5$Locus == "") & (case5$TSSUTR == "") & (case5$termUTR == "")),]
termUTR.case5 <- case5[((case5$Locus == "") & (case5$TSSUTR == "") & (!case5$termUTR == "")),]
other1 <- case5[((!case5$Locus == "") & (!case5$TSSUTR == "") & (!case5$termUTR == "")),]
other2 <- case5[((case5$Locus == "") & (!case5$TSSUTR == "") & (!case5$termUTR == "")),]
other3 <- case5[((!case5$Locus == "") & (case5$TSSUTR == "") & (!case5$termUTR == "")),]
other4 <- case5[((!case5$Locus == "") & (!case5$TSSUTR == "") & (case5$termUTR == "")),]
other.case5 <- rbind(other1, other2, other3, other4)

EnrichedOutput <- as.data.frame(matrix(nrow = 5,ncol = 13))
colnames(EnrichedOutput) <- c("cases", "totalnumber", "genic.number", "genicpercentage", "intergenic.number", "intergenic.percentage", "UTR.number", "UTR.percentage", "termUTR.number", "term.UTR.percentage","other.number","other.percentage", "total.percentage")

EnrichedOutput$cases <- c("case1", "case2", "case3", "case4", "case5")

EnrichedOutput$totalnumber <- c(nrow(case1), nrow(case2), nrow(case3), nrow(case4), nrow(case5))

EnrichedOutput$genic.number <- c(nrow(genic.case1), nrow(genic.case2), nrow(genic.case3), nrow(genic.case4), nrow(genic.case5))

EnrichedOutput$intergenic.number <- c(nrow(intergenic.case1), nrow(intergenic.case2), nrow(intergenic.case3), nrow(intergenic.case4), nrow(intergenic.case5))

EnrichedOutput$UTR.number <- c(nrow(UTR.case1), nrow(UTR.case2), nrow(UTR.case3), nrow(UTR.case4), nrow(UTR.case5))

EnrichedOutput$termUTR.number <- c(nrow(termUTR.case1), nrow(termUTR.case2), nrow(termUTR.case3), nrow(termUTR.case4), nrow(termUTR.case5))

EnrichedOutput$other.number <- c(nrow(other.case1), nrow(other.case2), nrow(other.case3), nrow(other.case4), nrow(other.case5))

for(i in 1:nrow(EnrichedOutput)){
  EnrichedOutput$genicpercentage[i] <- 100*EnrichedOutput$genic.number[i]/EnrichedOutput$totalnumber[i]
}

for(i in 1:nrow(EnrichedOutput)){
  EnrichedOutput$intergenic.percentage[i] <- 100*EnrichedOutput$intergenic.number[i]/EnrichedOutput$totalnumber[i]
}

for(i in 1:nrow(EnrichedOutput)){
  EnrichedOutput$UTR.percentage[i] <- 100*EnrichedOutput$UTR.number[i]/EnrichedOutput$totalnumber[i]
}

for(i in 1:nrow(EnrichedOutput)){
  EnrichedOutput$term.UTR.percentage[i] <- 100*EnrichedOutput$termUTR.number[i]/EnrichedOutput$totalnumber[i]
}

for(i in 1:nrow(EnrichedOutput)){
  EnrichedOutput$other.percentage[i] <- 100*EnrichedOutput$other.number[i]/EnrichedOutput$totalnumber[i]
}

for(i in 1:nrow(EnrichedOutput)){
  EnrichedOutput$total.percentage[i] <- EnrichedOutput$genicpercentage[i] + EnrichedOutput$intergenic.percentage[i] + EnrichedOutput$UTR.percentage[i] + EnrichedOutput$term.UTR.percentage[i] + EnrichedOutput$other.percentage[i]
}


write.csv(EnrichedOutput, file = "EnrichedOutput-Ksg150UTR.csv", sep = "", col.names = F, row.names = F)
