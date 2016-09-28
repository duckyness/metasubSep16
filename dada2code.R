#load dada2
library(dada2)
#load shortread
library(ShortRead)
#load ggplot2
library(ggplot2)
#define path as folder containing data
path <- "rachel_train_Sep16_Sep28/demultiplexed"
#define fns to list files in path
fns <- list.files(path)
fns
#define fastqs
fastqs <- fns[grepl(".fastq.gz$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnRs <- fastqs[grepl(".r1", fastqs)] # Just the reverse read files
fnFs <- fastqs[grepl(".r2", fastqs)] # Just the forward read files
# Get sample names from the first part of the forward read filenames
sample.names <- sapply(strsplit(fnFs, ".r2"), `[`, 1)
# Fully specify the path for the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
#examine quality profiles of forward and reverse reads
plotQualityProfile(fnFs[[1]])
plotQualityProfile(fnRs[[1]])
# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
# Filter
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                    trimLeft=c(10, 10), truncLen=c(200,270), 
                    maxN=0, maxEE=5, truncQ=2, 
                    compress=TRUE, verbose=TRUE)}
#dereplicate filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
#inspect
derepFs
derepRs
#perform joint sample inference and error rate estimation
dadaFs <- dada(derepFs, err=NULL, selfConsist = TRUE)
dadaRs <- dada(derepRs, err=NULL, selfConsist = TRUE, MAX_CONSIST = 20)
#inspect
dadaFs
dadaRs
#Visualize estimated error rates, nominalQ=true gives red line
plotErrors(dadaFs, nominalQ=TRUE)
plotErrors(dadaRs, nominalQ=TRUE)
#merge forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
head(mergers[[2]])
head(mergers[[3]])
head(mergers[[4]])
head(mergers[[5]])
#construct sequence table
seqtab <- makeSequenceTable(mergers[names(mergers) != "HMP"])
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(colnames(seqtab)))
#remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "rachel_train_Aug16_demux_5sept/rdp_train_set_14.fa.gz")
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
unname(head(taxa))
#DADA2 accuracy on mock community
unqs.mock <- getUniques(removeBimeraDenovo(mergers[["filtered/HMP"]], verbose=TRUE))
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
mockRef <- readFasta(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, as.character(sread(mockRef))))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
