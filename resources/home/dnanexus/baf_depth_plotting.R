# SCRIPT TO PLOT BAF AND DEPTH WITH KARYOPLOTER
#
# @input - list of files in tsv format with "CHROM", "POS", "DP" and "AD"
# @output - png plot of BAF and binned Depth values for each file
#
# Constraints
##################
# Depth values that are higher than the max value in Y axis will be plotted in the BAF plot
# Only .tsv files can be provided
# Chromosome names need to be provided as it defaults to include "chr"

# Load required modules
##################################
library(stringr, quietly = TRUE)
library(karyoploteR, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(polars, quietly = TRUE)
library(argparse, quietly = TRUE)
library(data.table, quietly = TRUE)

# Get configurable inputs via command line
####################################################

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("--vcf", type="character", default="*.vcf.tsv",
    help="Input VCF-derived file in tsv format [default %(default)s]")
parser$add_argument("--gvcf", type="character", default="*.gvcf.tsv",
    help="Input GVCF-derived file in tsv format [default %(default)s]")
parser$add_argument("--min_baf", type="double", default=0.04,
    help="Minimum BAF threshold displayed [default %(default)s]")
parser$add_argument("--max_baf", type="double", default=0.96,
    help="Maximum BAF threshold displayed [default %(default)s]")
parser$add_argument("--max_depth", type="double", default=0.9, 
    help = "Max depth to be shown on plot [default %(default)s]")
parser$add_argument("--min_depth", type="integer", default=5,
    help="Minimum depth allowed [default %(default)s]")
parser$add_argument("--bin_size", type="integer",
    help="Bin size for reducing noise in depth plot")
parser$add_argument("--chr_names", type="character", default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y",
    help="Chromosome names [default %(default)s]")
parser$add_argument("--genome", type="character", default="hg19",
    help="Genome build for plotKaryotype function [default %(default)s]")
parser$add_argument("--symmetry", type="logical", default=TRUE,
    help="Whether to plot BAF symmetrically [default %(default)s]")
parser$add_argument("--output_tsv", type="logical", default=FALSE,
    help="Whether to write TSV outputs [default %(default)s]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

# Validate percentile parameter
if (args$max_depth < 0 || args$max_depth > 1) {
  stop("max_depth must be a percentile value between 0 and 1")
}

# get sample name & check they're the same for the two inputs
SAMPLE_NAME <- str_split_1(sub("\\.vcf\\.tsv$", "", basename(args$vcf)), "_")[1]
if (SAMPLE_NAME != str_split_1(sub("\\.gvcf\\.tsv$", "", basename(args$gvcf)), "_")[1]) {
  stop("Sample names in VCF and GVCF files do not match")
}

# Configurables
##################

# Input file names
VCF_FILE <- args$vcf
GVCF_FILE <- args$gvcf

# Minimum and maximum thresholds of BAFs to be plotted
MIN_BAF <- args$min_baf
MAX_BAF <- args$max_baf

# Adjusts the Y-axis plot for mean_depth values
MAX_DEPTH_PCT <- args$max_depth

# Chromosome labels to feature in the plot X-axis
CHR_NAMES <- strsplit(args$chr_names, ",")[[1]]

# Minimum depth
MIN_DEPTH <- args$min_depth

# Bin size for depth plot
BIN_SIZE <- args$bin_size

# Genome build for plotKaryotype function
GENOME <- args$genome

# Option to make BAF plot symmetrical
SYMMETRY <- isTRUE(args$symmetry)


# List of functions
##################################

# Function to read files into dfs
# @parameter file
# returns df
read_to_df <- function(file) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  # check file is not empty
  if (file.info(file)$size == 0) {
    # make an empty dataframe
    print(paste("File", file, "is empty. Returning an empty dataframe"))
    empty_df <- data.frame(Chr = character(0), Position = numeric(0), Depth = numeric(0), BAF = numeric(0))
    return(empty_df)
  }
  else {
    df <- data.table::fread(file, header = FALSE, sep = "\t", na.strings = c("NA", ".", ""))
    if (ncol(df) == 4) {
    colnames(df) <- c("Chr", "Position", "Depth", "Allele_Depth")
    # Remove rows with NA in Depth or Allele_Depth
    df <- df[!is.na(df$Allele_Depth) & !is.na(df$Depth), ]
    } else if (ncol(df) == 3) {
      colnames(df) <- c("Chr", "Position", "Depth")
      # Remove rows with NA in Depth
      df <- df[!is.na(df$Depth), ]
    }
    else {
      stop("Invalid TSV format: Expected 3 or 4 columns (CHROM, POS, DP[, AD])")
    }
    # Check Position & Depth are numeric
    if (!is.numeric(df$Position) || !is.numeric(df$Depth)) {
        stop("Invalid TSV: Position and Depth must be numeric")
    }
    # Return Dataframe
    return(as.data.frame(df))
}
}

# Function to calculate baf 
# @parameter df
# @parameter sym - logical: whether to add symmetrical BAF values
# Returns df 
calculate_baf <- function(df, sym) {
  df[c("Ref_AD", "Alt_AD")] <- str_split_fixed(df$Allele_Depth, ",", 3)[, 1:2]
  df$RAF <- as.numeric(df$Ref_AD)
  df$BAF <- as.numeric(df$Alt_AD)
  # Avoid division by zero
  df$RAF <- ifelse(df$Depth > 0, as.numeric(df$Ref_AD) / df$Depth, NA)
  df$BAF <- ifelse((as.numeric(df$Ref_AD) + as.numeric(df$Alt_AD)) > 0, 
                   as.numeric(df$Alt_AD) / (as.numeric(df$Ref_AD) + as.numeric(df$Alt_AD)), NA) #calculate BAF strictly relative to the read counts of the two main alleles
      # Add symmetrical values if required
  if (isTRUE(sym)) {
      symmetric_df <- df
      symmetric_df$BAF <- 1 - df$BAF
      df <- rbind(df, symmetric_df)
    }
  df <- df[!is.na(df$BAF), ]
  
  return(as.data.frame(df))
} 

# Function to return dfs with binned depth
# @parameter df
# @parameter bin_size - integer: length of window for variant aggregation
# returns df_binned

bin_df <- function(df, bin_size) {
  # Convert to data.table by reference (no memory copy)
  setDT(df)
  # Sort in place
  setorder(df, Chr, Position)
  # Create bin groups and calculate mean depth simultaneously
  df[, bin := floor(Position / bin_size) * bin_size]
  df_binned <- df[, .(Position = bin[1], mean_depth = mean(Depth, na.rm = TRUE)), by = .(Chr, bin)]
  # Remove the temporary bin column
  df_binned[, bin := NULL]
  return(as.data.frame(df_binned))
}


# Function to create snp.data.baf format as input to karyoploteR plot BAF
# @parameter df
# returns snp.data for baf plot

get_snp_data_BAF <- function(df) {
  snp.data.baf <- toGRanges(df[, c("Chr", "Position", "Position", "Depth", "BAF")])
  names(mcols(snp.data.baf)) <- c("Depth", "BAF")
  seqlevelsStyle(snp.data.baf) <- "UCSC" # convert to "chr" to comply with USCS rules

  return(snp.data.baf)
}

# Function to create snp.data.depth format as input to karyoploteR plot binned depth
# @parameter df, use a binned df
# returns snp.data for depth plot

get_snp_data_Depth <- function(df) {
  snp.data.depth <- toGRanges(df[, c("Chr", "Position", "Position", "mean_depth")])
  names(mcols(snp.data.depth)) <- c("mean_depth")
  seqlevelsStyle(snp.data.depth) <- "UCSC" # convert to "chr" to comply with USCS rules

  return(snp.data.depth)
}

# Function to create a karyoploteR plot of BAF vs binned Depth
# @parameters snp.data.baf, snp.data.depth, file_name
# @parameter max_depth - integer : maximum depth value for y-axis in plots
# @parameter chr_names - vector with chromosomes to be listed
# @parameters min_baf and max_baf - integers : include only variants in range min_baf < BAF < max_baf
# returns plot

get_plot <- function(snp.data.baf, snp.data.depth, file_name, max_depth, chr_names, min_baf, max_baf, genome_build, median_depths) {
  file_name_png <- paste0(SAMPLE_NAME, ".baf.png")
  png(file_name_png, width = 15, height = 5, units = "in", res = 600)
  plot_parameters <- getDefaultPlotParams(plot.type = 4)
  plot_parameters$data1inmargin <- 2
  baf_depth_plot <- plotKaryotype(genome = genome_build, plot.type = 4, ideogram.plotter = NULL, plot.params = plot_parameters, labels.plotter = NULL)
  kpAddChromosomeNames(baf_depth_plot, chr.names = chr_names)
  kpAddCytobandsAsLine(baf_depth_plot) # Add centromeres

  # Alpha transparency added here to prevent blobs
  col_baf <- adjustcolor("darkorange2", alpha.f = 0.3)
  col_depth_norm <- adjustcolor("darkblue", alpha.f = 0.3)
  col_depth_high <- adjustcolor("magenta", alpha.f = 0.3)

  # top graph
  baf_threshold <- which(snp.data.baf$BAF > min_baf & snp.data.baf$BAF <= max_baf)
  modified_high_depth <- snp.data.depth$mean_depth > max_depth 
  snp.data.depth$mean_depth <- pmin(snp.data.depth$mean_depth, max_depth) 

  modified_depth <- ifelse(modified_high_depth, col_depth_high, col_depth_norm) 

  kpAxis(baf_depth_plot, r0 = 0.55, r1 = 1, tick.pos = c(0, 0.25, 0.5, 0.75, 1))
  kpAbline(baf_depth_plot, h=c(0.25, 0.5, 0.75), lty = 0.5, r0 =0.55, r1=1)

  if (length(snp.data.baf) > 0) {
    # Reduced point size (cex)
    kpPoints(baf_depth_plot,
      data = snp.data.baf[baf_threshold], y = snp.data.baf[baf_threshold]$BAF,
      cex = 0.2, r0 = 0.55, r1 = 1, col = col_baf
    )
    title <- paste0("BAF vs Depth for ", SAMPLE_NAME)
  } else {
    title <- paste0("Provided VCFs contain insufficient data for plotting - ", SAMPLE_NAME)
  }
  # bottom graph
  kpAxis(baf_depth_plot, r0 = 0, r1 = 0.45, ymax = max_depth, ymin = 0)
  kpPoints(baf_depth_plot,
    data = snp.data.depth, y = snp.data.depth$mean_depth,
    cex = 0.2, r0 = 0, r1 = 0.45, ymax = max_depth, ymin = 0, col = modified_depth
  )

  # add horizontal line at median depth for each chromosome
  df <- data.frame(x = snp.data.depth)
  for (chr in unique(df$x.seqnames)) {
    # Subset the vector using the chromosome name
    median_depth <- median_depths[chr]
    if (!is.na(median_depth)) {
      prop <- median_depth / max_depth # scale to the bottom plot
      kpAbline(baf_depth_plot, h=prop, chr=chr, col = "darkred", lwd = 3, r0 = 0, r1 = 0.45)
    }
  }
  kpAddMainTitle(baf_depth_plot, main = title)
  kpAddChromosomeSeparators(baf_depth_plot, col = "darkgray", lty = 3, data.panel = "all")

  dev.off()
}

# call functions
####################

# read tsv file into df for BAF plot
df_vcf <- read_to_df(VCF_FILE)

# read tsv file into df for DEPTH plot
df_gvcf <- read_to_df(GVCF_FILE)

# Filter out low depth rows for BAF plotting
df_filtered <- df_vcf[df_vcf$Depth >= MIN_DEPTH, ]
if (nrow(df_filtered) == 0) {
  print("No rows with Depth >= MIN_DEPTH found. Returning original df_vcf.")
  df_filtered <- df_vcf
}

# DOWNSAMPLE BAF for readability and memory reduction
# 50,000 points is more than enough to clearly see zygosity without blacking out the plot

set.seed(42)

MAX_PLOT_POINTS <- 50000
if (nrow(df_filtered) > MAX_PLOT_POINTS) {
  print(paste("Downsampling BAF data from", nrow(df_filtered), "to", MAX_PLOT_POINTS, "points for plot readability..."))
  df_filtered <- df_filtered[sample(nrow(df_filtered), MAX_PLOT_POINTS), ]
}

# calculate BAF values and add symmetrical values if required
if (! (nrow(df_filtered) == 0)) {
  df_filtered <- calculate_baf(df_filtered, SYMMETRY)
}

if (args$output_tsv) {
  print("User requested TSVs. Writing calculated data before plotting to disk...")
  write.table(df_filtered, file=paste0(SAMPLE_NAME, ".vcf.baf.tsv"), quote=FALSE, sep='\t', col.names = NA)
  write.table(df_gvcf, file=paste0(SAMPLE_NAME, ".gvcf.baf.tsv"), quote=FALSE, sep='\t', col.names = NA)
}

# get quantiles for plotting limits
MAX_DEPTH <- round(quantile(df_gvcf$Depth, probs = MAX_DEPTH_PCT, names = FALSE), digits = 0)

# dynamic bin size choice
if (is.null(BIN_SIZE)) {
  if ((length(df_gvcf$Depth) / 2000) >= 1) {
    BIN_SIZE <- round(length(df_gvcf$Depth) / 2000)
  } else {
    BIN_SIZE <- 1
  }
}

# aggregate gvcf df into binned df for depth plot
df_binned <- bin_df(df_gvcf, BIN_SIZE)
snp.data.depth_full <- get_snp_data_Depth(df_binned)
df_full <- data.frame(x=snp.data.depth_full)
median_depths <- tapply(df_full$x.mean_depth, df_full$x.seqnames, median, na.rm = TRUE)

# DOWNSAMPLE BINNED DATA for readability and memory reduction
# 50,000 points is more than enough to clearly see depth distribution without blacking out the plot
MAX_BINNED_POINTS <- 50000
if (nrow(df_binned) > MAX_BINNED_POINTS) {
  print(paste("Downsampling binned depth data from", nrow(df_binned), "to", MAX_BINNED_POINTS, "points for plot readability..."))
  # Use systematic sampling to maintain genomic distribution (taking every Nth row)
  step_size <- ceiling(nrow(df_binned) / MAX_BINNED_POINTS)
  df_binned <- df_binned[seq(1, nrow(df_binned), by = step_size), ]
}

# MEMORY CLEAR: Removed massive intermediate write.tables here.
rm(df_vcf, df_gvcf)
gc()

print("Converting filtered dataframe to snp.data.baf format for karyoploter...")
snp.data.baf <- get_snp_data_BAF(df_filtered)

# convert dfs for depth plotting into snp.data for karyoploter
print("Converting binned dataframe to snp.data.depth format for karyoploter...")
snp.data.depth <- get_snp_data_Depth(df_binned)

# generate plots and save them
print("Generating BAF vs Depth plot...")
get_plot(snp.data.baf, snp.data.depth, SAMPLE_NAME, MAX_DEPTH, CHR_NAMES, MIN_BAF, MAX_BAF, GENOME, median_depths)
