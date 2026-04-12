# Locus Zoom for Plasmodium falciparum
# Adapted from Tanya Major & Riku Takei's original code
# Modified for Plasmodium falciparum by Gabe, Belinda and Kashif
# April 2026 - FINAL VERSION

#### Function to make LocusZoom like plots for P. falciparum ####
locus.zoom <- function(data = NULL, snp = NA, gene = NA, region = NA, 
                       ld.file = NULL, offset_bp = 50000, 
                       genes.data = NULL, psuedogenes = FALSE, 
                       RNAs = FALSE, plot.title = NULL, 
                       plot.type = "jpg", nominal = 3, 
                       significant = 4, 
                       file.name = NULL, secondary.snp = NA, 
                       secondary.label = FALSE, secondary.circle = TRUE, 
                       genes.pvalue = NULL, colour.genes = FALSE, 
                       sig.type = "P", 
                       nplots = FALSE, ignore.lead = FALSE, 
                       rsid.check = FALSE, nonhuman = TRUE) {
  
  # Load scales library
  if("scales" %in% data.frame(installed.packages())[, "Package"]){
    library("scales")
  } else{
    stop("This function requires the package 'scales' to run.")
  }
  
  # Helper function to convert P. falciparum chromosome names
  convert_pf_chr <- function(chr) {
    chr <- as.character(chr)
    if(grepl("Pf3D7_", chr)) {
      num <- gsub("Pf3D7_0*(\\d+)_.*", "\\1", chr)
      return(as.character(as.numeric(num)))
    } else if(grepl("^0", chr)) {
      return(as.character(as.numeric(chr)))
    }
    return(chr)
  }
  
  # Helper function to extract position from SNP name
  extract_position <- function(snp_name) {
    as.numeric(gsub(".*:", "", snp_name))
  }
  
  # Load Data
  if (is.data.frame(data)) {
    lead.data = data
  } else {
    lead.data = data[[1]]
  }
  
  # Clean column names
  names(lead.data) <- gsub("^#", "", names(lead.data))
  names(lead.data) <- gsub("\\.", "", names(lead.data))
  
  # Map column names
  if("CHROM" %in% names(lead.data)) names(lead.data)[names(lead.data) == "CHROM"] <- "CHR"
  if("POS" %in% names(lead.data)) names(lead.data)[names(lead.data) == "POS"] <- "BP"
  if("ID" %in% names(lead.data)) names(lead.data)[names(lead.data) == "ID"] <- "SNP"
  
  # Error check
  if(!all(c("CHR", "BP", "SNP", "P") %in% names(lead.data))){
    stop(paste("Data must contain CHR, BP, SNP, and P columns.\n",
               "Available:", paste(names(lead.data), collapse=", ")))
  }
  
  lead.data$SNP = as.character(lead.data$SNP)
  lead.data$CHR <- sapply(lead.data$CHR, convert_pf_chr)
  lead.data$BP <- as.numeric(lead.data$BP)
  lead.data$P <- as.numeric(lead.data$P)
  lead.data$P[lead.data$P == 0 | is.na(lead.data$P)] <- .Machine$double.xmin
  lead.data$logP <- -log10(lead.data$P)
  
  # Check gene data
  if(!is.null(genes.data)) {
    if(!all(c("Gene", "Chrom", "Start", "End") %in% names(genes.data))){
      stop("Genes data must contain Gene, Chrom, Start, and End columns")
    }
    genes.data$Chrom <- sapply(genes.data$Chrom, convert_pf_chr)
    genes.data$Start <- as.numeric(genes.data$Start)
    genes.data$End <- as.numeric(genes.data$End)
    if("Coding" %in% names(genes.data)) {
      genes.data$Coding <- as.character(genes.data$Coding)
    }
  }
  
  # Get initial region from SNP or gene
  if (all(is.na(region))) {
    region = get.region(lead.data, snp, genes.data, gene) 
  }
  
  # Extract lead SNP info given by the user
  if (!is.na(snp)) {
    if(grepl("Pf3D7_", snp)) {
      snp_parts <- strsplit(snp, ":")[[1]]
      snp_chr_num <- gsub("Pf3D7_0*(\\d+)_.*", "\\1", snp_parts[1])
      snp_pos <- as.numeric(snp_parts[2])
      user_lead_snp <- snp
      user_lead_pos <- snp_pos
      user_lead_chr <- as.character(snp_chr_num)
    } else if(grepl(":", snp)) {
      snp_parts <- strsplit(snp, ":")[[1]]
      user_lead_snp <- snp
      user_lead_pos <- as.numeric(snp_parts[2])
      user_lead_chr <- as.character(snp_parts[1])
    } else {
      user_lead_snp <- snp
      snp_ind <- which(lead.data$SNP == snp)
      if(length(snp_ind) > 0) {
        user_lead_pos <- lead.data$BP[snp_ind[1]]
        user_lead_chr <- as.character(lead.data$CHR[snp_ind[1]])
      } else {
        user_lead_pos <- as.numeric(region[2])
        user_lead_chr <- as.character(region[1])
      }
    }
  } else {
    user_lead_snp <- NA
    user_lead_pos <- as.numeric(region[2])
    user_lead_chr <- as.character(region[1])
  }
  
  lead_chr_temp = user_lead_chr
  lead_pos_temp = user_lead_pos
  
  cat("\n=== INITIAL REGION ===\n")
  cat("Lead SNP:", user_lead_snp, "\n")
  cat("Position:", lead_pos_temp, "\n")
  cat("User-provided offset:", offset_bp, "bp\n")
  
  # AUTO-ADJUST OFFSET BASED ON LD FILE
  adjusted_offset <- offset_bp
  
  if (!is.null(ld.file) && nrow(ld.file) > 0) {
    cat("\n=== ANALYZING LD FILE ===\n")
    
    snp_b_col <- grep("SNP_B$", names(ld.file), ignore.case = TRUE, value = TRUE)[1]
    if(is.na(snp_b_col)) {
      snp_b_col <- grep("SNP_B|SNP2", names(ld.file), ignore.case = TRUE, value = TRUE)[1]
    }
    
    if(!is.na(snp_b_col)) {
      ld_positions <- extract_position(ld.file[[snp_b_col]])
      min_ld_pos <- min(ld_positions, na.rm = TRUE)
      max_ld_pos <- max(ld_positions, na.rm = TRUE)
      
      needed_left_offset <- lead_pos_temp - min_ld_pos
      needed_right_offset <- max_ld_pos - lead_pos_temp
      needed_offset <- max(needed_left_offset, needed_right_offset)
      
      if(needed_offset > offset_bp) {
        adjusted_offset <- ceiling(needed_offset * 1.1)
        cat("Offset increased to:", adjusted_offset, "bp to capture all LD SNPs\n")
      } else {
        cat("Current offset sufficient to capture all LD SNPs\n")
      }
    }
  }
  
  # Define region with adjusted offset
  region_start = lead_pos_temp - adjusted_offset
  region_end = lead_pos_temp + adjusted_offset
  region_chr = lead_chr_temp
  
  # Filter genes to region
  if(!is.null(genes.data)) {
    genes.data = genes.data[genes.data$Chrom == region_chr, ]
    genes.data = genes.data[genes.data$End > region_start, ]
    genes.data = genes.data[genes.data$Start < region_end, ]
    
    if(!psuedogenes && "Coding" %in% names(genes.data)) {
      genes.data = genes.data[!(genes.data$Coding %in% c("pseudogene", "Non-Coding")), ]
    }
    if(!RNAs && "Coding" %in% names(genes.data)) {
      genes.data = genes.data[!(genes.data$Coding %in% c("lncRNA", "ncRNA")), ]
    }
  }
  
  # Subset data to region
  data = subset.data(lead.data, region_chr, region_start, region_end)
  data$logP <- -log10(data$P)
  
  if(nrow(data) == 0) {
    stop(paste("No SNPs found in region", region_chr, ":", region_start, "-", region_end))
  }
  
  cat("\n=== GWAS DATA IN REGION ===\n")
  cat("Total SNPs in region:", nrow(data), "\n")
  
  # Get lead variant - USE USER-SPECIFIED SNP
  if (!is.na(snp)) {
    lead.ind = which(data$SNP == user_lead_snp)
    if(length(lead.ind) == 0) {
      lead.ind = which.min(abs(data$BP - user_lead_pos))
      cat("Using closest SNP at position", data$BP[lead.ind], "as lead\n")
    }
  } else {
    lead.ind = which(data$logP %in% max(data$logP, na.rm = TRUE))[1]
  }
  
  if(length(lead.ind) == 0) {
    lead.ind = 1
  }
  
  lead.snp = data$SNP[lead.ind]
  lead.chr = as.character(data$CHR[lead.ind])
  lead.pos = data$BP[lead.ind]
  lead.logp = data$logP[lead.ind]
  
  cat("\n=== LEAD SNP FOR PLOTTING ===\n")
  cat("Lead SNP:", lead.snp, "\n")
  cat("Position:", lead.pos, "\n")
  cat("-log10(P):", round(lead.logp, 2), "\n")
  
  # PROCESS LD FILE
  if (!is.null(ld.file) && nrow(ld.file) > 0) {
    
    cat("\n=== LD FILE PROCESSING ===\n")
    
    r2_col <- grep("^R2$", names(ld.file), ignore.case = TRUE, value = TRUE)[1]
    if(is.na(r2_col)) {
      r2_col <- grep("R2|R_SQ", names(ld.file), ignore.case = TRUE, value = TRUE)[1]
    }
    
    snp_b_col <- grep("SNP_B$", names(ld.file), ignore.case = TRUE, value = TRUE)[1]
    if(is.na(snp_b_col)) {
      snp_b_col <- grep("SNP_B|SNP2", names(ld.file), ignore.case = TRUE, value = TRUE)[1]
    }
    
    if(!is.na(r2_col) && !is.na(snp_b_col)) {
      
      ld_data <- ld.file[, c(snp_b_col, r2_col)]
      names(ld_data) <- c("SNP", "R2")
      ld_data$R2 <- as.numeric(as.character(ld_data$R2))
      ld_data$R2[is.na(ld_data$R2)] <- 0
      
      ld_data <- ld_data[order(ld_data$R2, decreasing = TRUE), ]
      ld_data <- ld_data[!duplicated(ld_data$SNP), ]
      
      if(!lead.snp %in% ld_data$SNP) {
        ld_data <- rbind(ld_data, data.frame(SNP = lead.snp, R2 = 1))
      }
      
      cat("Total SNPs in LD file:", nrow(ld_data), "\n")
      
      # PRINT SUMMARY OF LD SNPs
      cat("\n=== LD SNPs SUMMARY ===\n")
      cat(sprintf("%-30s %10s %12s\n", "SNP", "Position", "R2"))
      cat(rep("-", 55), "\n", sep="")
      
      ld_data_sorted <- ld_data[order(extract_position(ld_data$SNP)), ]
      for(i in 1:nrow(ld_data_sorted)) {
        pos <- extract_position(ld_data_sorted$SNP[i])
        r2 <- ld_data_sorted$R2[i]
        cat(sprintf("%-30s %10d %12.6f\n", ld_data_sorted$SNP[i], pos, r2))
      }
      cat(rep("-", 55), "\n", sep="")
      
      # Get unique R2 values
      unique_r2 <- sort(unique(ld_data$R2))
      n_unique <- length(unique_r2)
      
      # Create color palette
      if(n_unique == 1) {
        exact_colors <- "#B2182B"
      } else {
        exact_colors <- colorRampPalette(c(
          "#2166AC", "#d56bb7ae", "#67A9CF", "#d6bf25", "#D1E5F0",
          "#adcf16", "#F4A582", "#D6604D", "#09a40e", "#0e0e0e"
        ))(n_unique)
      }
      
      color_map <- setNames(exact_colors, as.character(sort(unique_r2)))
      
      LD.colours <- data.frame(
        R2_value = sort(unique_r2),
        Colour = exact_colors,
        stringsAsFactors = FALSE
      )
      
      data.plot <- merge(data, ld_data, by = "SNP", all.x = TRUE)
      
      data.plot$Colour <- sapply(1:nrow(data.plot), function(i) {
        r2 <- data.plot$R2[i]
        if(is.na(r2)) {
          return("#D3D3D3")
        } else {
          r2_char <- as.character(r2)
          if(r2_char %in% names(color_map)) {
            return(color_map[r2_char])
          } else {
            closest_idx <- which.min(abs(unique_r2 - r2))
            return(exact_colors[closest_idx])
          }
        }
      })
      
      data.plot$has_ld <- !is.na(data.plot$R2)
      
      cat("\n=== COLOR ASSIGNMENT ===\n")
      cat("SNPs with LD colors:", sum(data.plot$has_ld), "\n")
      cat("SNPs without LD (gray):", sum(!data.plot$has_ld), "\n")
      
    } else {
      data.plot = data
      data.plot$Colour = "#D3D3D3"
      data.plot$has_ld = FALSE
      LD.colours <- data.frame(R2_value = 0, Colour = "#D3D3D3")
    }
  } else {
    data.plot = data
    data.plot$Colour = "#D3D3D3"
    data.plot$has_ld = FALSE
    LD.colours <- data.frame(R2_value = 0, Colour = "#D3D3D3")
  }
  
  # Mark lead SNP
  data.plot$is_lead <- (data.plot$SNP == lead.snp)
  
  # Reorder for plotting (gray first, then colored, then lead on top)
  gray_snps <- data.plot[!data.plot$has_ld & !data.plot$is_lead, ]
  colored_snps <- data.plot[data.plot$has_ld & !data.plot$is_lead, ]
  lead_snp_row <- data.plot[data.plot$is_lead, ]
  
  gray_snps <- gray_snps[order(gray_snps$BP), ]
  colored_snps <- colored_snps[order(colored_snps$BP), ]
  
  data.plot <- rbind(gray_snps, colored_snps, lead_snp_row)
  
  ### Make Plot ###
  npanel = ifelse(nplots, length(data), 1)
  plot.height = (npanel * 80) + 50
  
  if(plot.type == "jpg"){
    jpeg(width = 200, height = plot.height, units = "mm", res = 300, filename = file.name)
  } else if(plot.type == "png") {
    png(width = 200, height = plot.height, units = "mm", res = 300, filename = file.name)
  } else if(plot.type == "svg") {
    svg(width = (200 / 25.4), height = (plot.height / 25.4), filename = file.name)
  } else if(plot.type == "pdf") {
    pdf(file.name, width = 200/25.4, height = plot.height/25.4)
  } else if(plot.type != "view_only"){
    stop("Plot type must be jpg, png, svg, pdf, or view_only")
  }
  
  mat.row = (2 * npanel) + 1
  locus.par = c(4, 15)
  layout(matrix(c(1:mat.row), byrow = TRUE), heights = c(rep(locus.par, npanel), 10))
  
  x.min = region_start
  x.max = region_end
  
  # Create locus plot
  if (npanel > 1) {
    for (i in 1:npanel) {
      tmp.dat = data.plot[[i]]
      y.max = max(tmp.dat$logP, 8, na.rm = TRUE)
      plot.var = c(y.max, x.min, x.max, lead.snp, nominal, significant)
      plot.locus(data.plot = tmp.dat, plot.title = names(data.plot)[i], 
                 secondary.snp = secondary.snp, secondary.label = secondary.label, 
                 secondary.circle = secondary.circle, sig.type = sig.type, plot.var = plot.var,
                 ld.colours = LD.colours)
      rm(tmp.dat)
    }
  } else {
    y.max = max(data.plot$logP, 8, na.rm = TRUE)
    plot.var = c(y.max, x.min, x.max, lead.snp, nominal, significant)
    plot.locus(data.plot = data.plot, plot.title = plot.title, 
               secondary.snp = secondary.snp, secondary.label = secondary.label, 
               secondary.circle = secondary.circle, sig.type = sig.type, plot.var = plot.var,
               ld.colours = LD.colours)
  }
  
  # Plot Gene tracks
  par(mar = c(4, 4, 0.5, 8), mgp = c(2, 1, 0), xpd = FALSE)
  
  n_genes <- ifelse(is.null(genes.data), 0, nrow(genes.data))
  if(n_genes > 20) {
    track.max = 5
    font.size = 0.55
  } else if(n_genes > 10) {
    track.max = 4
    font.size = 0.65
  } else {
    track.max = 3
    font.size = 0.75
  }
  
  plot(1, type = "n", yaxt = "n", 
       xlab = paste("Position on Chromosome", lead.chr, "(bp)"), 
       ylab = "", xlim = c(x.min, x.max), ylim = c(0, track.max), 
       xaxt = "n", cex.lab = 1)
  
  x_marks = axTicks(side = 1)
  axis(side = 1, at = x_marks, 
       labels = format(x_marks, scientific = FALSE, big.mark = ",", trim = TRUE),
       cex.axis = 0.8)
  
  if(!is.null(genes.data) && nrow(genes.data) != 0) {
    
    if(colour.genes && !is.null(genes.pvalue)) {
      if(!all(c("Gene", "P") %in% names(genes.pvalue))){
        stop("genes.pvalue must contain Gene and P columns")
      }
      genes.pvalue = genes.pvalue[genes.pvalue$Gene %in% genes.data$Gene, ]
      genes.data = merge.gene.colour(genes.data, genes.pvalue, GENE.colours)
    } else {
      genes.data$Colour = "#7F7F7F"
    }
    
    if(n_genes > 20) {
      y_positions = rep(c(4.5, 3.5, 2.5, 1.5, 0.5), length.out = n_genes)
      genes.data$Y = y_positions
    } else if(n_genes > 10) {
      y_positions = rep(c(3.5, 2.5, 1.5, 0.5), length.out = n_genes)
      genes.data$Y = y_positions
    } else {
      y_positions = rep(c(2.5, 1.5, 0.5), length.out = n_genes)
      genes.data$Y = y_positions
    }
    
    for (i in 1:nrow(genes.data)) {
      rect(xleft = max(genes.data$Start[i], x.min), 
           xright = min(genes.data$End[i], x.max),
           ybottom = genes.data$Y[i] - 0.1, 
           ytop = genes.data$Y[i] + 0.1,
           col = alpha(genes.data$Colour[i], 0.6), 
           border = "black", lwd = 0.8)
      
      mid_x = (max(genes.data$Start[i], x.min) + min(genes.data$End[i], x.max)) / 2
      text(x = mid_x, 
           y = genes.data$Y[i] + 0.25, 
           labels = genes.data$Gene[i], 
           font = 3, 
           cex = font.size,
           srt = 0,
           adj = 0.5)
    }
  }
  
  if( plot.type != "view_only" ) {
    dev.off()
  }
  
  cat("\n=== PLOT COMPLETE ===\n")
  cat("Plot saved to:", file.name, "\n")
}

#### Function to create LocusZoom style plot  ####
plot.locus <- function(data.plot = NULL, plot.title = NULL, secondary.snp = NA, 
                       secondary.label = FALSE, secondary.circle = TRUE, 
                       sig.type = "P", plot.var = NULL, ld.colours = NULL) {
  
  y.max = as.numeric(plot.var[1])
  x.min = as.numeric(plot.var[2])
  x.max = as.numeric(plot.var[3])
  lead.snp = plot.var[4]
  nominal = as.numeric(plot.var[5])
  significant = as.numeric(plot.var[6])
  
  # Add 10% padding to y-axis
  y_padding <- max(1, y.max * 0.1)
  y_limit <- c(0, y.max + y_padding)
  
  # Plot SNP presence
  par(mar = c(0, 4, 2, 8), mgp = c(2, 1, 0), xpd = FALSE)
  plot(x = data.plot$BP, y = rep(1, times = nrow(data.plot)), axes = FALSE, 
       pch = "|", xlab = "", ylab = "Plotted\nSNPs", las = 2, 
       xlim = c(x.min, x.max), cex.lab = 0.7, cex.axis = 0.6, 
       col = alpha(colour = "black", alpha = 0.15))
  title(plot.title, line = 0, cex.main = 1)
  
  # Plot Manhattan/LocusZoom (original style)
  par(mar = c(0, 4, 0, 8), mgp = c(2, 1, 0), xpd = FALSE)
  ylab = expression(-log[10](italic(P)))
  
  # Plot all SNPs
  plot(x = data.plot$BP, y = data.plot$logP, ylim = y_limit, 
       pch = 19, col = alpha(as.character(data.plot$Colour), 0.8), 
       xlab = "", ylab = ylab, cex = 1.2, xaxt = "n", 
       xlim = c(x.min, x.max))
  
  # Draw significance lines 
  abline(h = nominal, col = "blue", lty = "dashed", lwd = 1.5)
  abline(h = significant, col = "red", lty = "dashed", lwd = 1.5)
  
  
  
  # Plot the lead SNP (purple diamond, larger and bold)
  if (lead.snp %in% data.plot$SNP) {
    ind = which(data.plot$SNP == lead.snp)
    lead.pos = data.plot$BP[ind]
    lead.logp = data.plot$logP[ind]
    points(x = lead.pos, y = lead.logp, pch = 18, cex = 2.5, col = "#9400D3", lwd = 2)
    text(x = lead.pos, y = lead.logp + (y.max * 0.02), 
         labels = lead.snp, pos = 3, cex = 0.7, font = 2)
  }
  
  # Plot secondary SNPs if provided
  if(any(!is.na(secondary.snp))){
    secondary.snp <- secondary.snp[secondary.snp != lead.snp]
    check = which(data.plot$SNP %in% secondary.snp)
    if (length(check) != 0) {
      secondary.data = data.plot[data.plot$SNP %in% c(secondary.snp, lead.snp), ]
      plot.secondary.point(data = secondary.data, snps = secondary.data$SNP, 
                          lead.snp = lead.snp, plot.data = data.plot, 
                          plot.var = plot.var, label = secondary.label, 
                          circle = secondary.circle)
    }
  }
  
  # Add CLEAN LEGEND 
  par(xpd = TRUE)
  
  has_ld_colors <- !is.null(ld.colours) && nrow(ld.colours) > 1 && 
                   any(data.plot$has_ld == TRUE, na.rm = TRUE)
  
  if(has_ld_colors) {
    # Get R2 values for legend
    r2_values <- ld.colours$R2_value[ld.colours$R2_value > 0]
    r2_colors <- ld.colours$Colour[ld.colours$R2_value > 0]
    n_r2 <- length(r2_values)
    
    # Format R2 labels nicely
    format_r2 <- function(x) {
      if(x < 0.0001) return(sprintf("%.2e", x))
      if(x < 0.01) return(sprintf("%.4f", x))
      if(x < 0.1) return(sprintf("%.3f", x))
      if(x < 1) return(sprintf("%.3f", x))
      return(sprintf("%.1f", x))
    }
    
    # Select legend entries (max 10 for clean display)
    if(n_r2 <= 10) {
      legend_entries <- 1:n_r2
      legend_cex <- 0.6
      legend_ncol <- 1
    } else {
      legend_entries <- round(seq(1, n_r2, length.out = 10))
      legend_cex <- 0.55
      legend_ncol <- 1
    }
    
    legend_colors <- r2_colors[legend_entries]
    legend_labels <- sapply(r2_values[legend_entries], format_r2)
    
    # LD legend (top)
    legend(x = "topright", 
           legend = legend_labels, 
           col = legend_colors, 
           fill = legend_colors, 
           border = legend_colors, 
           pt.cex = 1, 
           cex = legend_cex, 
           bg = "white", 
           box.lwd = 0.5, 
           title = expression(bold(r^2)), 
           ncol = legend_ncol,
           inset = c(-0.12, 0))
    
    # Calculate position for next legend
    legend_height <- length(legend_entries) * 0.12
    legend_y_offset <- legend_height + 0.05
    
    # No LD data legend (middle)
    legend(x = "topright",
           legend = c("No LD data"),
           col = c("#D3D3D3"),
           fill = c("#D3D3D3"),
           border = c("#D3D3D3"),
           pt.cex = 1,
           cex = 0.6,
           bg = "white",
           box.lwd = 0.5,
           title = "",
           inset = c(-0.12, legend_y_offset))
    
    # Lead SNP legend (bottom)
    legend(x = "topright",
           legend = c("Lead SNP"),
           col = c("#9400D3"),
           pch = c(18),
           pt.cex = 1.5,
           cex = 0.6,
           bg = "white",
           box.lwd = 0.5,
           title = "",
           inset = c(-0.12, legend_y_offset * 2))
    
  } else {
    # Simple legend when no LD data
    legend(x = "topright",
           legend = c("SNPs", "Lead SNP"),
           col = c("#D3D3D3", "#9400D3"),
           pch = c(19, 18),
           pt.cex = c(1.2, 1.5),
           cex = 0.65,
           bg = "white",
           box.lwd = 0.5,
           title = "Legend",
           inset = c(-0.12, 0))
  }
}

#### Helper functions ####

get.region <- function(snp.dat, snp, gene.dat, gene) {
  if (!is.na(snp)) {
    if(grepl("Pf3D7_", snp)) {
      snp_parts <- strsplit(snp, ":")[[1]]
      snp_chr_num <- gsub("Pf3D7_0*(\\d+)_.*", "\\1", snp_parts[1])
      snp_pos <- as.numeric(snp_parts[2])
      region = c(snp_chr_num, snp_pos, snp_pos)
    } else if(grepl(":", snp)) {
      snp_parts <- strsplit(snp, ":")[[1]]
      region = c(snp_parts[1], as.numeric(snp_parts[2]), as.numeric(snp_parts[2]))
    } else {
      snp.ind = which(snp.dat$SNP == snp)
      if(length(snp.ind) == 0){
        stop(paste0("SNP ", snp, " not found"))
      }
      region = c(snp.dat$CHR[snp.ind], snp.dat$BP[snp.ind], snp.dat$BP[snp.ind])
    }
  } else if (!is.na(gene) && !is.null(gene.dat)) {
    gene.ind = which(gene.dat$Gene == gene)
    if(length(gene.ind) == 0){
      stop(paste0("Gene ", gene, " not found"))
    }
    region = c(gene.dat$Chrom[gene.ind], gene.dat$Start[gene.ind], gene.dat$End[gene.ind])
  } else {
    stop("Must specify SNP or Gene")
  }
  return(region)
}

subset.data <- function(data, chr, start, end) {
  res = data
  res = res[as.character(res$CHR) == as.character(chr), ]
  res = res[res$BP >= start, ]
  res = res[res$BP <= end, ]
  return(res)
}

plot.secondary.point <- function(data, snps, lead.snp, plot.data, plot.var, 
                                  label = FALSE, circle = TRUE) {
  
  lead.ind = which(data$SNP == lead.snp)
  if(length(lead.ind) > 0) {
    lead.pos = data$BP[lead.ind]
  } else {
    lead.pos = mean(as.numeric(plot.var[2:3]))
  }
  
  data = data[which(data$SNP != lead.snp), ]
  data = data[data$BP >= as.numeric(plot.var[2]) & data$BP <= as.numeric(plot.var[3]), ]
  
  if (circle && nrow(data) > 0) {
    points(x = data$BP, y = data$logP, pch = 1, cex = 1.8, col = "#FF0000", lwd = 1.5)
  }
  
  if (label && nrow(data) > 0) {
    x.min = as.numeric(plot.var[2])
    x.max = as.numeric(plot.var[3])
    x.offset = abs(x.max - x.min) / 150 * 15
    
    data$label.x.offset[data$BP < lead.pos] = data$BP[data$BP < lead.pos] - x.offset / 3
    data$side[data$BP < lead.pos] = 2
    data$label.x.offset[data$BP >= lead.pos] = data$BP[data$BP >= lead.pos] + x.offset / 3
    data$side[data$BP >= lead.pos] = 4
    data$label.y.offset = NA
    
    for(s in data$SNP){
      ind = which(data$SNP == s)
      pos = data$BP[ind]
      logp = data$logP[ind]
      surrounding.data = plot.data[plot.data$BP > (pos - (abs(x.max - x.min) / 6)) & 
                                   plot.data$BP < (pos + (abs(x.max - x.min) / 6)), ]
      surrounding.data = surrounding.data[surrounding.data$logP > logp & 
                                          surrounding.data$logP < logp + 5, ]
      if(length(surrounding.data[, 1]) == 0){
        data$label.y.offset[ind] = (logp + 2) * 1.03
      } else {
        data$label.y.offset[ind] = max(surrounding.data$logP) * 1.03
      }
    }
    
    for(s in data$SNP){
      ind = which(data$SNP == s)
      pos = data$BP[ind]
      logp = data$logP[ind]
      label.x = data$label.x.offset[ind]
      label.y = data$label.y.offset[ind]
      side = data$side[ind]
      
      text(x = label.x, y = (label.y * 1.01), labels = s, cex = 0.65, pos = side, offset = 0.2)
      segments(x0 = pos, x1 = label.x, y0 = logp, y1 = label.y * 1.01, col = "gray50", lwd = 0.8)
    }
  }
}

merge.gene.colour <- function(data, pvalues, GENE.colours) {
  res = merge(data, pvalues, by = "Gene", all.x = TRUE)
  res$gene.col = cut(res$P, breaks = c(0, 1e-20, 1e-15, 1e-10, 2.6e-6, 1), 
                     labels = c("<1e-20", ">1e-20", ">1e-15", ">1e-10", ">2.6e-6"), 
                     include.lowest = TRUE, right = TRUE)
  res = merge(res, GENE.colours, by.x = "gene.col", by.y = "Threshold", all.x = TRUE)
  res$Colour[is.na(res$Colour)] = "#7F7F7F"
  res = res[order(res$Start), ]
  return(res)
}
