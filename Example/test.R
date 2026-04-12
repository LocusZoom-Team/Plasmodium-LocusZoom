#!/usr/bin/env Rscript
# ============================================================================
# LocusZoom Plot Generator for Plasmodium falciparum
# ============================================================================
# 
# USAGE:
# Rscript test.R <assoc_file> <ld_file> <snp> <genes_file> [secondary_snp] [offset] [title] [output_base]

# Load required packages
suppressPackageStartupMessages({
  library(vroom)
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
assoc_file <- args[1]
ld_file <- args[2]
snp <- args[3]
genes_file <- args[4]
secondary_snps <- ifelse(length(args) >= 5, args[5], NA)
offset <- ifelse(length(args) >= 6, as.numeric(args[6]), 50000)
title <- ifelse(length(args) >= 7, args[7], "LocusZoom Plot")
output_file <- ifelse(length(args) >= 8, args[8], "LocusZoom_Plot.jpg")

# Source the function
source("../functions/locus_zoom.R")

# Check files
if(!file.exists(assoc_file)) stop("Association file not found")
if(!file.exists(ld_file)) stop("LD file not found")
if(!file.exists(genes_file)) stop("Genes file not found")

# Read files
cat("\n========================================\n")
cat("Reading files...\n")

data <- vroom(assoc_file, delim = "\t", show_col_types = FALSE, comment = "")
data <- as.data.frame(data)
cat("Association:", nrow(data), "rows\n")

ld <- read.table(ld_file, header = TRUE, stringsAsFactors = FALSE)
cat("LD file:", nrow(ld), "rows\n")

# Display LD SNP range
ld$BP <- as.numeric(gsub(".*:", "", ld$SNP_B))
cat("\nLD SNP position range:", min(ld$BP), "to", max(ld$BP), "\n")
cat("LD SNP span:", max(ld$BP) - min(ld$BP), "bp\n")

# Find lead SNP position
lead_snp_data <- data[data$SNP == snp, ]
if(nrow(lead_snp_data) > 0) {
  lead_pos <- lead_snp_data$BP[1]
  cat("\nLead SNP position:", lead_pos, "\n")
  needed_left <- lead_pos - min(ld$BP)
  needed_right <- max(ld$BP) - lead_pos
  cat("Recommended minimum offset:", max(needed_left, needed_right), "bp\n")
}

genes <- read.delim(genes_file, stringsAsFactors = FALSE, header = TRUE)
cat("Genes:", nrow(genes), "rows\n")

# Generate plot
cat("\n========================================\n")
cat("Generating LocusZoom plot...\n")
cat("  - Auto-offset will adjust to capture all LD SNPs\n")
cat("  - LD colors based on unique R2 values\n")
cat("========================================\n")

locus.zoom(data = data, snp = snp, ld.file = ld, offset_bp = offset,
           genes.data = genes, plot.title = title, nominal = 3, significant = 4,
           file.name = output_file, secondary.snp = secondary_snps,
           plot.type = "jpg", sig.type = "P", nonhuman = TRUE)

