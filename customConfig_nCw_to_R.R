#CUSTOMTYPE <- "nCw_to_R"
#setwd("/Volumes/SDB/Saini_Lab/Natalie/scleroderma/Mutation_signature/res_nCw_to_N/")

library(stringr)
filename <- list.files(pattern = "_seq.txt")
filename
filename_1 <- str_remove(filename, "_seq.txt")

file <- read.table(filename, sep="\t", header=T)
head(file)

library(stringr)

file$"C_to_R"=0
for (i in 1:nrow(file)){
if (file[i , "Reference_Allele"]=="C") {
  file[i,"C_to_R"]<-1} else { file[i, "C_to_R"]<-0}
}

file$"G_to_Y"=0
for (i in 1:nrow(file)){
  if (file[i , "Reference_Allele"]=="G" ) {
    file[i,"G_to_Y"]<-1} else { file[i, "G_to_Y"]<-0}
}

file$"nCw_to_R"=0
for (i in 1:nrow(file)){
  if (str_detect(file[i , "Trinucleotide_Motif"],"[A,T,G,C]C[A,T]")) {
    file[i,"nCw_to_R"]<-1} else { file[i, "nCw_to_R"]<-0}}

file$"wGn_to_Y"=0
for (i in 1:nrow(file)){
  if (str_detect(file[i , "Trinucleotide_Motif"],"[T,A]G[A,T,C,G]")) {
    file[i,"wGn_to_Y"]<-1} else { file[i, "wGn_to_Y"]<-0}
}

file$"c"=0
for (i in 1:nrow(file)){
    count <- (str_count(file[i , "Flanks"],"C"))
              if(file[i , "Reference_Allele"]=="C") {count <- count-1}
    file[i, "c"] <- count
}

file$"g"=0
for (i in 1:nrow(file)){
  count <- (str_count(file[i , "Flanks"],"G"))
  if(file[i , "Reference_Allele"]=="G") {count <- count-1}
  file[i, "g"] <- count
}

file$"ncw"=0
for (i in 1:nrow(file)){
  count <- (str_count(file[i , "Flanks"],"[A,T,G,C]C[A,T]"))
  if(str_detect(file[i , "Trinucleotide_Motif"],"[A,T,G,C]C[A,T]")) {count <- count-1}
  file[i, "ncw"] <- count
}

file$"wgn"=0
for (i in 1:nrow(file)){
  count <- (str_count(file[i , "Flanks"],"[T,A]G[A,T,C,G]"))
  if(str_detect(file[i , "Trinucleotide_Motif"],"[T,A]G[A,T,C,G]")) {count <- count-1}
  file[i, "wgn"] <- count
}

write.table(file, file = paste(filename_1, "_seq1.txt", sep = ""), sep="\t", col.names = T, row.names = F)

#### Now combine data for calculating sample-specific enrichment, mutloads, fisher, fisher-pcorr.

file <- read.table(file = paste(filename_1, "_seq1.txt", sep = ""), sep="\t", header=T)

head(file)

calc <- data.frame(Tumor_Sample_Barcode = unique(file$Tumor_Sample_Barcode))
calc$Tumor_Sample_Barcode = as.character(calc$Tumor_Sample_Barcode)

calc$"G_to_Y"=0
calc$"C_to_R"=0
calc$"nCw_to_R"=0
calc$"wGn_to_Y"=0
calc$"g"=0
calc$"c"=0
calc$"ncw"=0
calc$"wgn"=0
for(i in 1:nrow(calc)){
  sub <- subset(file, Tumor_Sample_Barcode==calc[i, "Tumor_Sample_Barcode"])
  calc[i,"G_to_Y"]=sum(sub$"G_to_Y")
  calc[i,"C_to_R"]=sum(sub$"C_to_R")
  calc[i,"nCw_to_R"]=sum(sub$"nCw_to_R")
  calc[i,"wGn_to_Y"]=sum(sub$"wGn_to_Y")
  calc[i,"g"]=sum(sub$"g")
  calc[i,"c"]=sum(sub$"c")
  calc[i,"ncw"]=sum(sub$"ncw")
  calc[i,"wgn"]=sum(sub$"wgn")
}
fin_row = nrow(calc)+1
calc[fin_row, "Tumor_Sample_Barcode"] = "Total"
calc[fin_row,"G_to_Y"]=sum(file$"G_to_Y")
calc[fin_row,"C_to_R"]=sum(file$"C_to_R")
calc[fin_row,"nCw_to_R"]=sum(file$"nCw_to_R")
calc[fin_row,"wGn_to_Y"]=sum(file$"wGn_to_Y")
calc[fin_row,"g"]=sum(file$"g")
calc[fin_row,"c"]=sum(file$"c")
calc[fin_row,"ncw"]=sum(file$"ncw")
calc[fin_row,"wgn"]=sum(file$"wgn")

#### Now only keep the mutatios in pyrimidines in the quotes

calc$"C_to_R_revcomp" <- calc$C_to_R + calc$G_to_Y
calc$"nCw_to_R_revcomp" <- calc$nCw_to_R + calc$wGn_to_Y
calc$"c_revcomp" <- calc$g + calc$c
calc$"ncw_revcomp" <- calc$ncw + calc$wgn

calc$"nCw_to_R_Enrich" <- ((calc$nCw_to_R_revcomp)*(calc$c_revcomp))/((calc$C_to_R_revcomp)*(calc$ncw_revcomp))

calc$"nCw_to_R_Enrich_Fisher_Pval"=0


for (i in 1:nrow(calc)) {
  x=c(calc[i, "nCw_to_R_revcomp"], (calc[i, "C_to_R_revcomp"] - calc[i, "nCw_to_R_revcomp"]))
  y=c((calc[i, "ncw_revcomp"]), (calc[i, "c_revcomp"] - calc[i, "ncw_revcomp"]))
  df=as.matrix(rbind(x,y))
  calc[i, "nCw_to_R_Enrich_Fisher_Pval"] = (fisher.test(df,alternative = "greater", conf.level = 0.95))$p.value
  }                                                 

BH_Corr_Pval <- p.adjust(calc$nCw_to_R_Enrich_Fisher_Pval, method="BH", n=nrow(calc))
calc$"BH_Corr_Pval" = BH_Corr_Pval

calc$"nCw_to_R_MutLoad" = 0

for(i in 1:nrow(calc)){
  if(calc[i, "BH_Corr_Pval"]<=0.05){
    calc[i, "nCw_to_R_MutLoad"] = ((calc[i, "nCw_to_R_revcomp"] * (calc[i, "nCw_to_R_Enrich"] - 1))/calc[i,"nCw_to_R_Enrich"])
    }
}

write.table(calc, file = paste(filename_1, "_fisherPCorr.txt", sep = ""), sep="\t", col.names = T, row.names = F)



