#Load table from genome annotation
gff_raw <- read.table("../LocusZoom-Plasmodium/Pfalciparum3D7.gff",
                      header    = FALSE,
                      sep       = "\t",
                      comment   = "#",       
                      quote     = "",
                      stringsAsFactors = FALSE)


names(gff_raw) <- c("seqnames", "source", "type", "start", "end",
                    "score", "strand", "phase", "attributes")

parse_gff_attr <- function(attributes, key) {
  pattern <- paste0("(?<=", key, "=)[^;]+") #Create function to extract from the table using characters to match the pattern
  match <- regmatches(attributes, regexpr(pattern, attributes, perl = TRUE))
  ifelse(length(match) == 0, NA, match)
}

gff_raw$ID          <- sapply(gff_raw$attributes, parse_gff_attr, key = "ID")
gff_raw$Name        <- sapply(gff_raw$attributes, parse_gff_attr, key = "Name")
gff_raw$description <- sapply(gff_raw$attributes, parse_gff_attr, key = "description")
gff_raw$Parent      <- sapply(gff_raw$attributes, parse_gff_attr, key = "Parent")
gff_raw$attributes  <- NULL
gff_df <- gff_raw

valid_chroms <- c(paste0("Pf3D7_", sprintf("%02d", 1:14), "_v3"),
                  "Pf3D7_API_v3",
                  "Pf_M76611")


genes <- gff_df[gff_df$type == "protein_coding_gene" & gff_df$seqnames %in% valid_chroms, ]

genes <- genes[, c("seqnames", "start", "end", "strand", "ID", "Name", "description")]
names(genes) <- c("Chrom", "Start", "End", "Strand", "GeneID", "Gene", "Description")

cds <- gff_df[gff_df$type == "CDS", ]
cds$Parent <- sapply(cds$Parent, function(x) x[1])

mrna <- gff_df[gff_df$type == "mRNA", ]
mrna$Parent <- sapply(mrna$Parent, function(x) x[1])
mrna_gene_map <- mrna[, c("ID", "Parent")]
names(mrna_gene_map) <- c("mRNA_ID", "GeneID")

cds <- merge(cds, mrna_gene_map, by.x = "Parent", by.y = "mRNA_ID", all.x = TRUE)

cds_summary <- data.frame(
  GeneID   = tapply(cds$GeneID, cds$GeneID, function(x) x[1]),
  cdsStart = tapply(cds$start,  cds$GeneID, min),
  cdsEnd   = tapply(cds$end,    cds$GeneID, max),
  row.names = NULL
)

genes <- merge(genes, cds_summary, by = "GeneID", all.x = TRUE)
genes$cdsStart[is.na(genes$cdsStart)] <- 0
genes$cdsEnd[is.na(genes$cdsEnd)]     <- 0

genes$GeneLength <- abs(genes$End - genes$Start)
genes$cdsLength  <- abs(genes$cdsEnd - genes$cdsStart)
genes$Coding     <- ifelse(genes$cdsLength > 0, "Coding", "Non-Coding")

genes <- genes[order(genes$GeneLength, decreasing = TRUE), ]
genes <- genes[!duplicated(genes$GeneID), ]

genes$Gene[is.na(genes$Gene)] <- genes$GeneID[is.na(genes$Gene)]

chrom_map <- c(setNames(paste0("chr", 1:14),
                        paste0("Pf3D7_", sprintf("%02d", 1:14), "_v3")),
               "Pf3D7_API_v3" = "chrAPI",
               "Pf_M76611"    = "chrMT")
genes$Chrom <- chrom_map[genes$Chrom]

chrom_order <- c(paste0("chr", 1:14), "chrAPI", "chrMT")
genes$Chrom <- factor(genes$Chrom, levels = chrom_order)
genes <- genes[order(genes$Chrom, genes$Start), ]
genes$Chrom <- as.character(genes$Chrom)

genes <- genes[, c("Gene", "Chrom", "Start", "End", "Coding")]

write.table(genes,
            file      = "PlasmoDB_Pfalciparum3D7_UniqueGeneList.txt",
            quote     = FALSE,
            row.names = FALSE,
            sep       = "\t",
            na        = "")

message("Done")