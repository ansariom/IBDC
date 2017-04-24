#!/usr/bin/Rscript
source('/nfs0/BPP/Megraw_Lab/mitra/Projects/3PEAT_model/r_scripts/read_annotations.R')
#file <- "peat.full_original_dataset.TSS.strandREannotated_NAMED.csv" 
#file <- "peat.full_original_dataset.TSS.strandREannotated_NAMED.csv"
args = commandArgs(trailingOnly=TRUE)
file <- args[1]
outfile <- args[2]
min_reads <- as.numeric(args[3])

##### Pull NP Tags #####
read_np <- function() {
    reads_per_tag <- 100
    outfile <<- "TSS100Narrow_PEAT_STRAND.txt"

    np_annos <- read_np_gene_annotations(reads_per_tag, file)

    # we only want tags in the following locations:
    # - TSS, 5'utr, <250, <500
    locations <- c("tss", "5'utr", "<250", "<500")
    filtered_annos <- subset(np_annos, TranscriptLocation %in% locations)

    return(filtered_annos )
}

read_br <- function() {
    reads_per_tag <- 100
    outfile <<- "TSS100Broad_PEAT_STRAND.txt"

    np_annos <- read_br_gene_annotations(reads_per_tag, file)

    # we only want tags in the following locations:
    # - TSS, 5'utr, <250, <500
    locations <- c("tss", "5'utr", "<250", "<500")
    filtered_annos <- subset(np_annos, TranscriptLocation %in% locations)

    return(filtered_annos )
}

read_wp <- function() {
    reads_per_tag <- 100
    outfile <<- "TSS100Weak_PEAT_STRAND.txt"

    np_annos <- read_wp_gene_annotations(reads_per_tag, file)

    # we only want tags in the following locations:
    # - TSS, 5'utr, <250, <500
    locations <- c("tss", "5'utr", "<250", "<500")
    filtered_annos <- subset(np_annos, TranscriptLocation %in% locations)

    return(filtered_annos )
}

read_all <- function() {
    reads_per_tag <- min_reads
    if (is.null(outfile)){
    	outfile <<- paste("TSS", min_reads, "All.txt")
    }

    np_annos <- read_gene_annotations(reads_per_tag, file)
    # we only want tags in the following locations:
    # - TSS, 5'utr, <250, <500
    locations <- c("tss", "5'utr", "<250", "<500")
    filtered_annos <- subset(np_annos, TranscriptLocation %in% locations)

    return(filtered_annos )
}

read_mirna <- function() {
    reads_per_tag <- 1
    outfile <<- "mirnaTSS_PEAT_STRAND.txt"

    # just pull out every miRNA annotation
    annos <- read_mirna_annotations(reads_per_tag, file)
    return(annos)
}

#### Pull NPs ####
#filtered_annos <- read_np()

#### Pull BRs ####
#filtered_annos <- read_br()

#### Pull WPs ####
#filtered_annos <- read_wp()

#### Pull miRNAs ####
#filtered_annos <- read_mirna()

#### Pull All ####

filtered_annos <- read_all()

selected_cols <- subset(filtered_annos, select=c(GeneName, Chromosome, Strand, ModeLocation))

##### Write Tags #####
#csv_outfile <- paste(outfile, "_ALL",min_reads , ".csv", sep="")
write.table(filtered_annos, file=outfile, row.names=FALSE,  quote=FALSE, sep=",")
#write.table(selected_cols, file=outfile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

