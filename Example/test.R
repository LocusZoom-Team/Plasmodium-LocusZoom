library(dplyr)
library(vroom)

example.assoc.log <- read.delim("example.assoc.log", stringsAsFactors = FALSE, header = TRUE)
Example.ld <- read.table("Example.ld", stringsAsFactors = FALSE, header = TRUE) ###########################only thing we need to change
Example.genes <- read.delim("Example.genes", stringsAsFactors = FALSE, header = TRUE)

# load the locuszoom function into R
source("functions/locus_zoom.R")

single = "rs2540781"
secondary = c("rs2540781", "rs8053279")
#############################################we will need to change this so that user will specify different SNP names that isnt RS, also change gene list
# create a LocusZoom-like plot ##################################################################################
source("functions/locus_zoom.R")
locus.zoom(data = example.assoc.log, snp = "rs1008400", ld.file = Example.ld, offset = 200000, genes.data = Example.genes, noncoding = FALSE, plot.title = "Association of FTO with BMI in Europeans", nominal = 6, significant = 7.3, file.name = "Example.jpg", secondary.snp = NA, population = "EUR", sig.type = "P")
#################################################################################################################


Stuff below this can likely be deleted

# Testing with one secondary variant
source("functions/locus_zoom.R")
locus.zoom(data = example.assoc.log, snp = "rs1008400", ld.file = Example.ld, offset = 200000, genes.data = Example.genes, noncoding = FALSE, plot.title = "Association of FTO with BMI in Europeans", nominal = 6, significant = 7.3, file.name = "Example.jpg", secondary.snp = single, population = "EUR", sig.type = "P")

# Testing with two secondary variant
source("functions/locus_zoom.R")
locus.zoom(data = example.assoc.log, snp = "rs1008400", ld.file = Example.ld, offset = 200000, genes.data = Example.genes, noncoding = FALSE, plot.title = "Association of FTO with BMI in Europeans", nominal = 6, significant = 7.3, file.name = "Example.jpg", secondary.snp = secondary, population = "EUR", sig.type = "P")

# Testing with a larger data set #####DELETE
Example.assoc.linear <- read.delim("eur_chr4.txt", stringsAsFactors = FALSE, header = TRUE)
Example.ld <- read.table("rs145179124.ld", stringsAsFactors = FALSE, header = TRUE)

example.assoc.log$P = example.assoc.log$P + .Machine$double.xmin

source("functions/locus_zoom.R")
locus.zoom(data = example.assoc.log, snp = "rs145179124", ld.file = Example.ld, offset = 500000, genes.data = Example.genes, noncoding = FALSE, plot.title = "EUR gout", nominal = 6, significant = 7.3, file.name = "rs145179124.jpg", secondary.snp = NA, population = "EUR", sig.type = "P")

# Try locating other loci as secondary SNP:
loci = read.delim('EUR_meta_full1_clean_rsid_0.01.clumped.clean', stringsAsFactors = F, header = T) %>% select(1:5)

example.assoc.log$P[which(example.assoc.log$SNP != 'rs145179124')] = example.assoc.log$P[which(example.assoc.log$SNP != 'rs145179124')] + .Machine$double.xmin
secondary = loci$SNP[which(loci$SNP != 'rs145179124')]

source("functions/locus_zoom.R")
locus.zoom(data = example.assoc.log, snp = "rs145179124", ld.file = Example.ld, offset = 500000, genes.data = Example.genes, noncoding = FALSE, plot.title = "EUR gout", nominal = 6, significant = 7.3, file.name = "rs145179124.jpg", secondary.snp = secondary, population = "EUR", sig.type = "P", secondary.label = T)

# Try with SLC2A9: ######DELETE
Example.assoc.linear <- read.delim("eur_chr4.txt", stringsAsFactors = FALSE, header = TRUE)
Example.ld <- read.table("rs76242518.ld", stringsAsFactors = FALSE, header = TRUE)

Example.assoc.linear$P = Example.assoc.linear$P + .Machine$double.xmin

Example.assoc.linear$P[which(Example.assoc.linear$SNP != 'rs76242518')] = Example.assoc.linear$P[which(Example.assoc.linear$SNP != 'rs76242518')] + .Machine$double.xmin
secondary = loci$SNP[which(loci$SNP != 'rs76242518')]

source("functions/locus_zoom.R")
locus.zoom(data = Example.assoc.linear, snp = "rs76242518", ld.file = Example.ld, offset = 500000, genes.data = Example.genes, noncoding = FALSE, plot.title = "EUR gout", nominal = 6, significant = 7.3, file.name = "rs76242518.jpg", secondary.snp = secondary, population = "EUR", sig.type = "P", secondary.label = T, nplots = 2)

