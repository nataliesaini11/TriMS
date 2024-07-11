#CUSTOMTYPE <- "gCn_to_A"
#setwd("/Volumes/sambashare/Saini_Lab/Natalie/test/flanks/")

library(stringr)
filename <- list.files(pattern = "_seq.txt")
filename
filename_1 <- str_remove(filename, "_seq.txt")

file <- read.table(filename, sep="\t", header=T)
head(file)

library(stringr)

file$"G_to_T"=0
for (i in 1:nrow(file)){
if (file[i , "Reference_Allele"]=="G" && file[i, "Tumor_Seq_Allele2"]=="T") {
  file[i,"G_to_T"]<-1} else { file[i, "G_to_T"]<-0}
}

file$"C_to_A"=0
for (i in 1:nrow(file)){
  if (file[i , "Reference_Allele"]=="C" && file[i, "Tumor_Seq_Allele2"]=="A") {
    file[i,"C_to_A"]<-1} else { file[i, "C_to_A"]<-0}
}

file$"gCn_to_A"=0
for (i in 1:nrow(file)){
  if (str_detect(file[i , "Trinucleotide_Motif"],"GC[A,T,G,C]") && file[i, "Tumor_Seq_Allele2"]=="A") {
    file[i,"gCn_to_A"]<-1} else { file[i, "gCn_to_A"]<-0}
}

file$"nGc_to_T"=0
for (i in 1:nrow(file)){
  if (str_detect(file[i , "Trinucleotide_Motif"],"[A,T,G,C]GC") && file[i, "Tumor_Seq_Allele2"]=="T") {
    file[i,"nGc_to_T"]<-1} else { file[i, "nGc_to_T"]<-0}
}

file$"g"=0
for (i in 1:nrow(file)){
    count <- (str_count(file[i , "Flanks"],"G"))
              if(file[i , "Reference_Allele"]=="G") {count <- count-1}
    file[i, "g"] <- count
}

file$"c"=0
for (i in 1:nrow(file)){
  count <- (str_count(file[i , "Flanks"],"C"))
  if(file[i , "Reference_Allele"]=="C") {count <- count-1}
  file[i, "c"] <- count
}

file$"gcn"=0
for (i in 1:nrow(file)){
  count <- (str_count(file[i , "Flanks"],"GC[A,T,C,G]"))
  if(str_detect(file[i , "Trinucleotide_Motif"],"GC[A,T,G,C]")) {count <- count-1}
  file[i, "gcn"] <- count
}

file$"ngc"=0
for (i in 1:nrow(file)){
  count <- (str_count(file[i , "Flanks"],"[A,T,C,G]GC"))
  if(str_detect(file[i , "Trinucleotide_Motif"],"[A,T,G,C]GC")) {count <- count-1}
  file[i, "ngc"] <- count
}

write.table(file, file = paste(filename_1, "_seq1.txt", sep = ""), sep="\t", col.names = T, row.names = F)

#### Now combine data for calculating sample-specific enrichment, mutloads, fisher, fisher-pcorr.

file <- read.table(file = paste(filename_1, "_seq1.txt", sep = ""), sep="\t", header=T)

head(file)

calc <- data.frame(Tumor_Sample_Barcode = unique(file$Tumor_Sample_Barcode))
calc$Tumor_Sample_Barcode = as.character(calc$Tumor_Sample_Barcode)

calc$"G_to_T"=0
calc$"C_to_A"=0
calc$"gCn_to_A"=0
calc$"nGc_to_T"=0
calc$"g"=0
calc$"c"=0
calc$"gcn"=0
calc$"ngc"=0
for(i in 1:nrow(calc)){
  sub <- subset(file, Tumor_Sample_Barcode==calc[i, "Tumor_Sample_Barcode"])
  calc[i,"G_to_T"]=sum(sub$"G_to_T")
  calc[i,"C_to_A"]=sum(sub$"C_to_A")
  calc[i,"gCn_to_A"]=sum(sub$"gCn_to_A")
  calc[i,"nGc_to_T"]=sum(sub$"nGc_to_T")
  calc[i,"g"]=sum(sub$"g")
  calc[i,"c"]=sum(sub$"c")
  calc[i,"gcn"]=sum(sub$"gcn")
  calc[i,"ngc"]=sum(sub$"ngc")
}
fin_row = nrow(calc)+1
calc[fin_row, "Tumor_Sample_Barcode"] = "Total"
calc[fin_row,"G_to_T"]=sum(file$"G_to_T")
calc[fin_row,"C_to_A"]=sum(file$"C_to_A")
calc[fin_row,"gCn_to_A"]=sum(file$"gCn_to_A")
calc[fin_row,"nGc_to_T"]=sum(file$"nGc_to_T")
calc[fin_row,"g"]=sum(file$"g")
calc[fin_row,"c"]=sum(file$"c")
calc[fin_row,"gcn"]=sum(file$"gcn")
calc[fin_row,"ngc"]=sum(file$"ngc")

#### Now only keep the mutatios in pyrimidines in the quotes

calc$"C_to_A_revcomp" <- calc$C_to_A + calc$G_to_T
calc$"gCn_to_A_revcomp" <- calc$gCn_to_A + calc$nGc_to_T
calc$"c_revcomp" <- calc$g + calc$c
calc$"gcn_revcomp" <- calc$gcn + calc$ngc

calc$"gCn_to_A_Enrich" <- ((calc$gCn_to_A_revcomp)*(calc$c_revcomp))/((calc$C_to_A_revcomp)*(calc$gcn_revcomp))

for (i in 1:nrow(calc)) {
  x=c(calc[i, "gCn_to_A_revcomp"], (calc[i, "C_to_A_revcomp"] - calc[i, "gCn_to_A_revcomp"]))
  y=c((calc[i, "gcn_revcomp"]), (calc[i, "c_revcomp"] - calc[i, "gcn_revcomp"]))
  df=as.matrix(rbind(x,y))
  calc[i, "gCn_to_A_Enrich_Fisher_Pval"] = (fisher.test(df,alternative = "greater", conf.level = 0.95))$p.value
  }                                                 

BH_Corr_Pval <- p.adjust(calc$gCn_to_A_Enrich_Fisher_Pval, method="BH", n=nrow(calc))
calc$"BH_Corr_Pval" = BH_Corr_Pval

calc$"gCn_to_A_MutLoad" = 0

for(i in 1:nrow(calc)){
  if(calc[i, "BH_Corr_Pval"]<=0.05){
    calc[i, "gCn_to_A_MutLoad"] = ((calc[i, "gCn_to_A_revcomp"] * (calc[i, "gCn_to_A_Enrich"] - 1))/calc[i,"gCn_to_A_Enrich"])
    }
}

write.table(calc, file = paste(filename_1, "_fisherPCorr.txt", sep = ""), sep="\t", col.names = T, row.names = F)



