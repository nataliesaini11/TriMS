#CUSTOMTYPE <- "nTw_to_N"
#setwd("/Volumes/SDB/Saini_Lab/Natalie/scleroderma/Mutation_signature/res_nTw_to_N/")

library(stringr)
filename <- list.files(pattern = "_seq.txt")
filename
filename_1 <- str_remove(filename, "_seq.txt")

file <- read.table(filename, sep="\t", header=T)
head(file)

library(stringr)

file$"T_to_N"=0
for (i in 1:nrow(file)){
if (file[i , "Reference_Allele"]=="T") {
  file[i,"T_to_N"]<-1} else { file[i, "T_to_N"]<-0}
}

file$"A_to_N"=0
for (i in 1:nrow(file)){
  if (file[i , "Reference_Allele"]=="A" ) {
    file[i,"A_to_N"]<-1} else { file[i, "A_to_N"]<-0}
}

file$"nTw_to_N"=0
for (i in 1:nrow(file)){
  if (str_detect(file[i , "Trinucleotide_Motif"],"[A,T,G,C]T[A,T]")) {
    file[i,"nTw_to_N"]<-1} else { file[i, "nTw_to_N"]<-0}}

file$"wAn_to_N"=0
for (i in 1:nrow(file)){
  if (str_detect(file[i , "Trinucleotide_Motif"],"[T,A]A[A,T,C,G]")) {
    file[i,"wAn_to_N"]<-1} else { file[i, "wAn_to_N"]<-0}
}

file$"t"=0
for (i in 1:nrow(file)){
    count <- (str_count(file[i , "Flanks"],"T"))
              if(file[i , "Reference_Allele"]=="T") {count <- count-1}
    file[i, "t"] <- count
}

file$"a"=0
for (i in 1:nrow(file)){
  count <- (str_count(file[i , "Flanks"],"A"))
  if(file[i , "Reference_Allele"]=="A") {count <- count-1}
  file[i, "a"] <- count
}

file$"ntw"=0
for (i in 1:nrow(file)){
  count <- (str_count(file[i , "Flanks"],"[A,T,G,C]T[A,T]"))
  if(str_detect(file[i , "Trinucleotide_Motif"],"[A,T,G,C]T[A,T]")) {count <- count-1}
  file[i, "ntw"] <- count
}

file$"wan"=0
for (i in 1:nrow(file)){
  count <- (str_count(file[i , "Flanks"],"[T,A]A[A,T,C,G]"))
  if(str_detect(file[i , "Trinucleotide_Motif"],"[T,A]A[A,T,C,G]")) {count <- count-1}
  file[i, "wan"] <- count
}

write.table(file, file = paste(filename_1, "_seq1.txt", sep = ""), sep="\t", col.names = T, row.names = F)

#### Now combine data for calculating sample-specific enrichment, mutloads, fisher, fisher-pcorr.

file <- read.table(file = paste(filename_1, "_seq1.txt", sep = ""), sep="\t", header=T)

head(file)

calc <- data.frame(Tumor_Sample_Barcode = unique(file$Tumor_Sample_Barcode))
calc$Tumor_Sample_Barcode = as.character(calc$Tumor_Sample_Barcode)

calc$"A_to_N"=0
calc$"T_to_N"=0
calc$"nTw_to_N"=0
calc$"wAn_to_N"=0
calc$"a"=0
calc$"t"=0
calc$"ntw"=0
calc$"wan"=0
for(i in 1:nrow(calc)){
  sub <- subset(file, Tumor_Sample_Barcode==calc[i, "Tumor_Sample_Barcode"])
  calc[i,"A_to_N"]=sum(sub$"A_to_N")
  calc[i,"T_to_N"]=sum(sub$"T_to_N")
  calc[i,"nTw_to_N"]=sum(sub$"nTw_to_N")
  calc[i,"wAn_to_N"]=sum(sub$"wAn_to_N")
  calc[i,"a"]=sum(sub$"a")
  calc[i,"t"]=sum(sub$"t")
  calc[i,"ntw"]=sum(sub$"ntw")
  calc[i,"wan"]=sum(sub$"wan")
}
fin_row = nrow(calc)+1
calc[fin_row, "Tumor_Sample_Barcode"] = "Total"
calc[fin_row,"A_to_N"]=sum(file$"A_to_N")
calc[fin_row,"T_to_N"]=sum(file$"T_to_N")
calc[fin_row,"nTw_to_N"]=sum(file$"nTw_to_N")
calc[fin_row,"wAn_to_N"]=sum(file$"wAn_to_N")
calc[fin_row,"a"]=sum(file$"a")
calc[fin_row,"t"]=sum(file$"t")
calc[fin_row,"ntw"]=sum(file$"ntw")
calc[fin_row,"wan"]=sum(file$"wan")

#### Now only keep the mutatios in pyrimidines in the quotes

calc$"T_to_N_revcomp" <- calc$T_to_N + calc$A_to_N
calc$"nTw_to_N_revcomp" <- calc$nTw_to_N + calc$wAn_to_N
calc$"t_revcomp" <- calc$a + calc$t
calc$"ntw_revcomp" <- calc$ntw + calc$wan

calc$"nTw_to_N_Enrich" <- ((calc$nTw_to_N_revcomp)*(calc$t_revcomp))/((calc$T_to_N_revcomp)*(calc$ntw_revcomp))

calc$"nTw_to_N_Enrich_Fisher_Pval"=0


for (i in 1:nrow(calc)) {
  x=c(calc[i, "nTw_to_N_revcomp"], (calc[i, "T_to_N_revcomp"] - calc[i, "nTw_to_N_revcomp"]))
  y=c((calc[i, "ntw_revcomp"]), (calc[i, "t_revcomp"] - calc[i, "ntw_revcomp"]))
  df=as.matrix(rbind(x,y))
  calc[i, "nTw_to_N_Enrich_Fisher_Pval"] = (fisher.test(df,alternative = "greater", conf.level = 0.95))$p.value
  }                                                 

BH_Corr_Pval <- p.adjust(calc$nTw_to_N_Enrich_Fisher_Pval, method="BH", n=nrow(calc))
calc$"BH_Corr_Pval" = BH_Corr_Pval

calc$"nTw_to_N_MutLoad" = 0

for(i in 1:nrow(calc)){
  if(calc[i, "BH_Corr_Pval"]<=0.05){
    calc[i, "nTw_to_N_MutLoad"] = ((calc[i, "nTw_to_N_revcomp"] * (calc[i, "nTw_to_N_Enrich"] - 1))/calc[i,"nTw_to_N_Enrich"])
    }
}

write.table(calc, file = paste(filename_1, "_fisherPCorr.txt", sep = ""), sep="\t", col.names = T, row.names = F)



