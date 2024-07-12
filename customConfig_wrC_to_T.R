#CUSTOMTYPE <- "wrC_to_T"
#setwd("/Volumes/sambashare/Saini_Lab/Natalie/test/flanks/")

library(stringr)
filename <- list.files(pattern = "_seq.txt")
filename
filename_1 <- str_remove(filename, "_seq.txt")

file <- read.table(filename, sep="\t", header=T)
head(file)

library(stringr)

file$"C_to_T"=0
for (i in 1:nrow(file)){
if (file[i , "Reference_Allele"]=="C" && file[i, "Tumor_Seq_Allele2"]=="T") {
  file[i,"C_to_T"]<-1} else { file[i, "C_to_T"]<-0}
}

file$"G_to_A"=0
for (i in 1:nrow(file)){
  if (file[i , "Reference_Allele"]=="G" && file[i, "Tumor_Seq_Allele2"]=="A") {
    file[i,"G_to_A"]<-1} else { file[i, "G_to_A"]<-0}
}

file$"wrC_to_T"=0
for (i in 1:nrow(file)){
  if (str_detect(file[i , "Trinucleotide_Motif_for_AID_UPSTREAM"],"[A,T][A,G]C") && file[i, "Tumor_Seq_Allele2"]=="T") {
    file[i,"wrC_to_T"]<-1} else { file[i, "wrC_to_T"]<-0}}

file$"Gys_to_A"=0
for (i in 1:nrow(file)){
  if (str_detect(file[i , "Trinucleotide_Motif_for_AID_DOWNSTREAM"],"G[T,C][C,G]") && file[i, "Tumor_Seq_Allele2"]=="A") {
    file[i,"Gys_to_A"]<-1} else { file[i, "Gys_to_A"]<-0}
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

file$"wrc"=0
for (i in 1:nrow(file)){
  count <- (str_count(file[i , "Flanks"],"[A,T][A,G]C"))
  if(str_detect(file[i , "Trinucleotide_Motif_for_AID_UPSTREAM"],"[A,T][A,G]C")) {count <- count-1}
  file[i, "wrc"] <- count
}

file$"gys"=0
for (i in 1:nrow(file)){
  count <- (str_count(file[i , "Flanks"],"G[T,C][C,G]"))
  if(str_detect(file[i , "Trinucleotide_Motif_for_AID_DOWNSTREAM"],"G[T,C][C,G]")) {count <- count-1}
  file[i, "gys"] <- count
}

write.table(file, file = paste(filename_1, "_seq1.txt", sep = ""), sep="\t", col.names = T, row.names = F)

#### Now combine data for calculating sample-specific enrichment, mutloads, fisher, fisher-pcorr.

file <- read.table(file = paste(filename_1, "_seq1.txt", sep = ""), sep="\t", header=T)

head(file)

calc <- data.frame(Tumor_Sample_Barcode = unique(file$Tumor_Sample_Barcode))
calc$Tumor_Sample_Barcode = as.character(calc$Tumor_Sample_Barcode)

calc$"G_to_A"=0
calc$"C_to_T"=0
calc$"wrC_to_T"=0
calc$"Gys_to_A"=0
calc$"g"=0
calc$"c"=0
calc$"wrc"=0
calc$"gys"=0
for(i in 1:nrow(calc)){
  sub <- subset(file, Tumor_Sample_Barcode==calc[i, "Tumor_Sample_Barcode"])
  calc[i,"G_to_A"]=sum(sub$"G_to_A")
  calc[i,"C_to_T"]=sum(sub$"C_to_T")
  calc[i,"wrC_to_T"]=sum(sub$"wrC_to_T")
  calc[i,"Gys_to_A"]=sum(sub$"Gys_to_A")
  calc[i,"g"]=sum(sub$"g")
  calc[i,"c"]=sum(sub$"c")
  calc[i,"wrc"]=sum(sub$"wrc")
  calc[i,"gys"]=sum(sub$"gys")
}
fin_row = nrow(calc)+1
calc[fin_row, "Tumor_Sample_Barcode"] = "Total"
calc[fin_row,"G_to_A"]=sum(file$"G_to_A")
calc[fin_row,"C_to_T"]=sum(file$"C_to_T")
calc[fin_row,"wrC_to_T"]=sum(file$"wrC_to_T")
calc[fin_row,"Gys_to_A"]=sum(file$"Gys_to_A")
calc[fin_row,"g"]=sum(file$"g")
calc[fin_row,"c"]=sum(file$"c")
calc[fin_row,"wrc"]=sum(file$"wrc")
calc[fin_row,"gys"]=sum(file$"gys")

#### Now only keep the mutatios in pyrimidines in the quotes

calc$"C_to_T_revcomp" <- calc$C_to_T + calc$G_to_A
calc$"wrC_to_T_revcomp" <- calc$wrC_to_T + calc$Gys_to_A
calc$"c_revcomp" <- calc$g + calc$c
calc$"wrc_revcomp" <- calc$wrc + calc$gys

calc$"wrC_to_T_Enrich" <- ((calc$wrC_to_T_revcomp)*(calc$c_revcomp))/((calc$C_to_T_revcomp)*(calc$wrc_revcomp))

calc$"wrC_to_T_Enrich_Fisher_Pval"=0


for (i in 1:nrow(calc)) {
  x=c(calc[i, "wrC_to_T_revcomp"], (calc[i, "C_to_T_revcomp"]-calc[i,"wrC_to_T_revcomp"]))
  y=c((calc[i, "wrc_revcomp"]), (calc[i, "c_revcomp"] - calc[i, "wrc_revcomp"]))
  df=as.matrix(rbind(x,y))
  calc[i, "wrC_to_T_Enrich_Fisher_Pval"] = (fisher.test(df,alternative = "greater", conf.level = 0.95))$p.value
  }                                                 

BH_Corr_Pval <- p.adjust(calc$wrC_to_T_Enrich_Fisher_Pval, method="BH", n=nrow(calc))
calc$"BH_Corr_Pval" = BH_Corr_Pval

calc$"wrC_to_T_MutLoad" = 0

for(i in 1:nrow(calc)){
  if(calc[i, "BH_Corr_Pval"]<=0.05){
    calc[i, "wrC_to_T_MutLoad"] = ((calc[i, "wrC_to_T_revcomp"] * (calc[i, "wrC_to_T_Enrich"] - 1))/calc[i,"wrC_to_T_Enrich"])
    }
}

write.table(calc, file = paste(filename_1, "_fisherPCorr.txt", sep = ""), sep="\t", col.names = T, row.names = F)



