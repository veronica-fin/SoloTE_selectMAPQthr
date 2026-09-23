args <- commandArgs(trailingOnly=TRUE)

counts_filename <- args[1]
outname <- args[2]

# Load data - Expecting: Feature, Barcode, UMI
counts <- read.delim(counts_filename, header=TRUE, sep="\t")
colnames(counts) <- c("Gene_id", "barcode", "umi")

# Create unique lists for features and barcodes
# This ensures that 'Gene;SoloTE' is treated as one unique feature
genes <- data.frame(id=unique(counts$Gene_id), name=unique(counts$Gene_id))
barcodes <- data.frame(id=unique(counts$barcode))

# Map the names to integer indices for the Matrix Market format
counts$matchgene <- match(counts$Gene_id, genes$id)
counts$matchbarcode <- match(counts$barcode, barcodes$id)

# Create Output Directory
dir.create(outname, showWarnings = FALSE)

# Write 10x compatible files
write.table(genes, file=paste0(outname, "/features.tsv"), 
            quote=FALSE, row.names=FALSE, sep="\t", col.names=FALSE)
write.table(barcodes, file=paste0(outname, "/barcodes.tsv"), 
            quote=FALSE, row.names=FALSE, sep="\t", col.names=FALSE)

# Write Matrix Market Header
cat("%%MatrixMarket matrix coordinate integer general\n%\n", file=paste0(outname, "/matrix.mtx"))
cat(nrow(genes), " ", nrow(barcodes), " ", nrow(counts), "\n", 
    file=paste0(outname, "/matrix.mtx"), append=TRUE, sep="")

# Write Data: [Feature Index] [Barcode Index] [Count]
write.table(counts[, c("matchgene", "matchbarcode", "umi")], 
            file=paste0(outname, "/matrix.mtx"), 
            sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)