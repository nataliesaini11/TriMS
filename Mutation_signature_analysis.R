#Trinucleotide Mutation Signature - TriMS
## Enter following commands"
# DATADIR="PATH_TO_DATA_DIR"
# FILESDIR = "where you want the output files to be"
#FILENAME="Name of the file.txt"
#REFPATH="Path to reference"
#PARENTDIR="Where the scripts are"
#CUSTOMTYPE="gCn_to_A"
#source(paste0(PARENTDIR,"Mutation_signature_analysis.R"))
#main()

CUSTOMMOTIFS <- TRUE

# execute this function to launch the pipeline

main <- function(
  dataDir = DATADIR,
  filesDir = FILESDIR,
  refFasta = REFPATH,
  filename = FILENAME,
  scriptDir = PARENTDIR,
  customtype = CUSTOMTYPE
)
{
  library(stringr)
  
    DATADIR <<- dataDir
    REF_FA <<- refFasta
 
  cat(FILESDIR, "\n")
  setwd(FILESDIR)
  
  FILENAME <<- filename
 
  SCRIPTDIR <<- scriptDir
  #### Getting flanks
file_big <- read.table(paste(dataDir, filename, sep=""), sep="\t", header=T)

file <- subset(file_big, Variant_Type == "SNP")

bedfile <- cbind((file$Chromosome), (file$Start_position)-21, (file$End_position)+20)
nrow(bedfile)

write.table(bedfile, file="bedfile.bed", col.names = F, row.names = F, sep="\t", quote = F)

cmdString <- paste("bedtools getfasta -fi ", REF_FA, " -bed bedfile.bed -tab > bedfile_seq.bed", sep="")
cat("executing:", cmdString, "\n")
system(cmdString)

file_big <- read.table(paste(dataDir, filename, sep=""), sep="\t", header=T)
file <- subset(file_big, Variant_Type == "SNP")
bedfile <- cbind((file$Chromosome), (file$Start_position)-2, (file$End_position)+1)
write.table(bedfile, file="bedfile_trinuc.bed", col.names = F, row.names = F, sep="\t", quote = F)

cmdString <- paste("bedtools getfasta -fi ", REF_FA, " -bed bedfile_trinuc.bed -tab > bedfile_trinuc_seq.bed", sep="")
cat("executing:", cmdString, "\n")
system(cmdString)

bedfileseq <- read.table("bedfile_seq.bed", sep = "\t", header=F)
bedfiletrinucseq <- read.table("bedfile_trinuc_seq.bed", sep = "\t", header=F)

Flanks <- bedfileseq$V2
Trinucleotide_Motif <- bedfiletrinucseq$V2
file_seq <- cbind(file, Flanks, Trinucleotide_Motif)
head(file_seq)

file_seq$Flanks <- toupper(file_seq$Flanks)
file_seq$Trinucleotide_Motif <- toupper(file_seq$Trinucleotide_Motif)

filename_1 <- str_remove(filename, ".txt")

write.table(file_seq, file = paste(filename_1, "_seq.txt", sep = ""), sep="\t", col.names=T, row.names=F, quote=F)

if (CUSTOMMOTIFS) {
  customConfigScript <- paste("customConfig_", customtype, ".R", sep="")
  source(paste(SCRIPTDIR, "/", customConfigScript, sep=""))
}

}


