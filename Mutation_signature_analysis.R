#Signature of Mutation Analysis -  SigMA
#Natalie's edited files for running signatures
# Enter following commands"
# DATADIR="/home/sainilab/sambashare/Saini_Lab/Natalie/test/flanks/"
# FILESDIR = "where you want the output files to be"
#FILENAME="Name of the file.txt"
#REFPATH="Path to reference" Note - for yeast use the edited files with Chromosome names written as >1 , >2...
#PARENTDIR="Where the scripts are"
# Edit mutation file so that Chromosome names are written as 1, 2, .... mit.
#CUSTOMTYPE="tCn_to_A"
#source(paste0(PARENTDIR,"Mutation_signature_analysis.R"))
#main()

CUSTOMMOTIFS <- TRUE

# execute this function to launch the pipeline

main <- function(
 # scriptDir = paste0(PARENTDIR, "code/"),
  dataDir = DATADIR,
  # use this option for the default APOBEC pipeline
  filesDir = FILESDIR,
  # use this option for the pipeline using a custom motif
  #filesDir = paste0(PARENTDIR, "work/res_yCn/"),
  refFasta = REFPATH,
  filename = FILENAME,
  scriptDir = PARENTDIR,
  customtype = CUSTOMTYPE
)
{
  library(stringr)
  
  #FILESDIR <<- filesDir
  DATADIR <<- dataDir
  # Reference Genome
  REF_FA <<- refFasta
 
  cat(FILESDIR, "\n")
  setwd(FILESDIR)
  
  FILENAME <<- filename
 
  SCRIPTDIR <<- scriptDir
  #### Getting flanks

#setwd("/Volumes/SDB/Saini_Lab/Natalie/scleroderma/Mutation_signature/")
file_big <- read.table(paste(dataDir, filename, sep=""), sep="\t", header=T)
#file_big <- read.table("scleroderma_correct_combined_snps.txt", sep="\t", header=T)

file <- subset(file_big, Variant_Type == "SNP")

#bedfile <- cbind((file$Chromosome), (file$Start_position)-21, (file$End_position)+20)

bedfile <- cbind(file$Chromosome, 
                 as.character(file$Start_position - 21), 
                 format(file$End_position + 20, scientific = FALSE, trim = TRUE, digits = 12))
nrow(bedfile)

write.table(bedfile, file="bedfile.bed", col.names = F, row.names = F, sep="\t", quote = F)

cmdString <- paste("bedtools getfasta -fi ", REF_FA, " -bed bedfile.bed -tab > bedfile_seq.bed", sep="")
cat("executing:", cmdString, "\n")
system(cmdString)

file_big <- read.table(paste(dataDir, filename, sep=""), sep="\t", header=T)
file <- subset(file_big, Variant_Type == "SNP")
#bedfile <- cbind((file$Chromosome), (file$Start_position)-2, (file$End_position)+1)
bedfile <- cbind(file$Chromosome, 
                 as.character(file$Start_position - 2), 
                 format(file$End_position + 1, scientific = FALSE, trim = TRUE, digits = 12))
write.table(bedfile, file="bedfile_trinuc.bed", col.names = F, row.names = F, sep="\t", quote = F)

cmdString <- paste("bedtools getfasta -fi ", REF_FA, " -bed bedfile_trinuc.bed -tab > bedfile_trinuc_seq.bed", sep="")
cat("executing:", cmdString, "\n")
system(cmdString)

### For AID signature
file_big <- read.table(paste(dataDir, filename, sep=""), sep="\t", header=T)
file <- subset(file_big, Variant_Type == "SNP")
#bedfile <- cbind((file$Chromosome), (file$Start_position)-3, (file$End_position))
bedfile <- cbind(file$Chromosome, 
                 as.character(file$Start_position - 3), 
                 format(file$End_position, scientific = FALSE, trim = TRUE, digits = 12))
write.table(bedfile, file="upstream_bedfile_trinuc_for_AID.bed", col.names = F, row.names = F, sep="\t", quote = F)

cmdString <- paste("bedtools getfasta -fi ", REF_FA, " -bed upstream_bedfile_trinuc_for_AID.bed -tab > upstream_bedfile_trinuc_seq_for_AID.bed", sep="")
cat("executing:", cmdString, "\n")
system(cmdString)

file_big <- read.table(paste(dataDir, filename, sep=""), sep="\t", header=T)
file <- subset(file_big, Variant_Type == "SNP")
#bedfile <- cbind((file$Chromosome), (file$Start_position-1), (file$End_position+2))
bedfile <- cbind(file$Chromosome, 
                 as.character(file$Start_position - 1), 
                 format(file$End_position + 2, scientific = FALSE, trim = TRUE, digits = 12))
write.table(bedfile, file="downstream_bedfile_trinuc_for_AID.bed", col.names = F, row.names = F, sep="\t", quote = F)

cmdString <- paste("bedtools getfasta -fi ", REF_FA, " -bed downstream_bedfile_trinuc_for_AID.bed -tab > downstream_bedfile_trinuc_seq_for_AID.bed", sep="")
cat("executing:", cmdString, "\n")
system(cmdString)

bedfileseq <- read.table("bedfile_seq.bed", sep = "\t", header=F)
bedfiletrinucseq <- read.table("bedfile_trinuc_seq.bed", sep = "\t", header=F)
upbedfiletrinucseqforAID <- read.table("upstream_bedfile_trinuc_seq_for_AID.bed", sep = "\t", header=F)
downbedfiletrinucseqforAID <- read.table("downstream_bedfile_trinuc_seq_for_AID.bed", sep = "\t", header=F)

Flanks <- bedfileseq$V2
Trinucleotide_Motif <- bedfiletrinucseq$V2
Trinucleotide_Motif_for_AID_UPSTREAM <- upbedfiletrinucseqforAID$V2
Trinucleotide_Motif_for_AID_DOWNSTREAM <- downbedfiletrinucseqforAID$V2

file_seq <- cbind(file, Flanks, Trinucleotide_Motif, Trinucleotide_Motif_for_AID_UPSTREAM, Trinucleotide_Motif_for_AID_DOWNSTREAM)
head(file_seq)

file_seq$Flanks <- toupper(file_seq$Flanks)
file_seq$Trinucleotide_Motif <- toupper(file_seq$Trinucleotide_Motif)
file_seq$Trinucleotide_Motif_for_AID_UPSTREAM <- toupper(file_seq$Trinucleotide_Motif_for_AID_UPSTREAM)
file_seq$Trinucleotide_Motif_for_AID_DOWNSTREAM <- toupper(file_seq$Trinucleotide_Motif_for_AID_DOWNSTREAM)

filename_1 <- str_remove(filename, ".txt")

write.table(file_seq, file = paste(filename_1, "_seq.txt", sep = ""), sep="\t", col.names=T, row.names=F, quote=F)

if (CUSTOMMOTIFS) {
  customConfigScript <- paste("customConfig_", customtype, ".R", sep="")
  source(paste(SCRIPTDIR, "/", customConfigScript, sep=""))
}

}


