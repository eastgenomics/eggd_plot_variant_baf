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

# Get configurable inputs via command line
####################################################

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("--vcf", type="character", default="*.vcf.tsv",
    help="Input VCF-derived file in tsv format [default %(default)s]")
parser$add_argument("--min_baf", type="double", default=0.04,
    help="Minimum BAF threshold displayed [default %(default)s]")
parser$add_argument("--max_baf", type="double", default=0.96,
    help="Maximum BAF threshold displayed [default %(default)s]")
parser$add_argument("--max_depth", type="double", default=0.9, 
    help = "Max depth to be shown on plot [default %(default)s]")
parser$add_argument("--min_depth", type="integer", default=5,
    help="Minimum depth allowed [default %(default)s]")
parser$add_argument("--chr_names", type="character", default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y",
    help="Chromosome names [default %(default)s]")
parser$add_argument("--genome", type="character", default="hg19",
    help="Genome build for plotKaryotype function [default %(default)s]")
parser$add_argument("--symmetry", type="logical", default=TRUE,
    help="Whether to plot BAF symmetrically [default %(default)s]")
                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

# Configurables
##################

# Input file name
VCF_FILE <- args$vcf

# Minimum and maximum thresholds of BAFs to be plotted
MIN_BAF <- args$min_baf
MAX_BAF <- args$max_baf

# Adjusts the Y-axis plot for mean_depth values
MAX_DEPTH <- args$max_depth

# Chromosome labels to feature in the plot X-axis
CHR_NAMES <- strsplit(args$chr_names, ",")[[1]]

# Minimum depth
MIN_DEPTH <- args$min_depth

# Genome build for plotKaryotype function
GENOME <- args$genome

# Option to make BAF plot symmetrical
SYMMETRY <- toupper(args$symmetry)


# List of functions
##################################

# Function to read files into dfs
# @parameter file
# returns df

read_to_df <- function(file) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  df <- read.table(file = file, header = FALSE)
  if (ncol(df) != 4) {
    stop("Invalid TSV format: Expected 4 columns (CHROM, POS, DP, AD)")
  }
  colnames(df) <- c("Chr", "Position", "Depth", "Allele_Depth")
  df <- df[!is.na(df$Allele_Depth) & !is.na(df$Depth), ]
  if (!all(sapply(df[c("Position", "Depth")], is.numeric))) {
    stop("Invalid TSV format: Position and Depth must be numeric")
  }

  # calculate BAF
  df[c("Ref_AD", "Alt_AD")] <- str_split_fixed(df$Allele_Depth, ",", 2)
  df$RAF <- as.numeric(df$Ref_AD)
  df$BAF <- as.numeric(df$Alt_AD)
  # Avoid division by zero
  df$RAF <- ifelse(df$Depth > 0, as.numeric(df$Ref_AD) / df$Depth, NA)
  df$BAF <- ifelse(df$Depth > 0, as.numeric(df$Alt_AD) / df$Depth, NA)

  # Add symmetrical values if required
  if (SYMMETRY) {
    symmetric_df <- df
    symmetric_df$BAF <- 1 - df$BAF
    df <- rbind(df, symmetric_df)
  }

  return(df)
}

# Function to return dfs with binned depth
# @parameter df
# @parameter bin_size - integer: length of window for variant aggregation
# returns df_binned

bin_df <- function(df, bin_size) {
  polars_df <- as_polars_df(df[order(df$Chr, df$Position), ])
  rolling_df <- polars_df$rolling(
    index_column = "Position",
    group_by = "Chr",
    period = paste(as.character(bin_size), "i", sep = "")
  )$agg(mean_depth = pl$mean("Depth"))$gather_every(bin_size)
  df_binned <- rolling_df$to_data_frame()

  return(df_binned)
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

get_plot <- function(snp.data.baf, snp.data.depth, file_name, max_depth, chr_names, min_baf, max_baf, genome_build) {
  file_name_png <- paste0(sub(".tsv", "", file_name), ".png")
  png(file_name_png, width = 15, height = 5, units = "in", res = 600)
  plot_parameters <- getDefaultPlotParams(plot.type = 4)
  plot_parameters$data1inmargin <- 2
  baf_depth_plot <- plotKaryotype(genome = genome_build, plot.type = 4, ideogram.plotter = NULL, plot.params = plot_parameters, labels.plotter = NULL)
  kpAddChromosomeNames(baf_depth_plot, chr.names = chr_names)
  kpAddCytobandsAsLine(baf_depth_plot) # Add centromeres
  # top graph
  baf_threshold <- which(snp.data.baf$BAF > min_baf & snp.data.baf$BAF <= max_baf)
  modified_high_depth <- snp.data.depth$mean_depth > max_depth # get values above max_depth
  snp.data.depth$mean_depth <- pmin(snp.data.depth$mean_depth, max_depth) # assign the max to max_depth
  modified_depth <- ifelse(
    modified_high_depth, 'darkgreen', 'darkblue'
  ) # Assign colors based on the mean_depth
  kpAxis(baf_depth_plot, r0 = 0.55, r1 = 1, tick.pos = c(0, 0.25, 0.5, 0.75, 1))
  kpAbline(baf_depth_plot, h=c(0.25, 0.5, 0.75), lty = 0.5, r0 =0.55, r1=1)
  kpPoints(baf_depth_plot,
    data = snp.data.baf[baf_threshold], y = snp.data.baf[baf_threshold]$BAF,
    cex = 0.5, r0 = 0.55, r1 = 1, col = "darkorange2"
  )
  # bottom graph
  kpAxis(baf_depth_plot, r0 = 0, r1 = 0.45, ymax = max_depth, ymin = 0)
  kpPoints(baf_depth_plot,
    data = snp.data.depth, y = snp.data.depth$mean_depth,
    cex = 0.5, r0 = 0, r1 = 0.45, ymax = max_depth, ymin = 0, col = modified_depth
  )
  kpAddMainTitle(baf_depth_plot, main = "BAF vs Depth")
  kpAddChromosomeSeparators(baf_depth_plot, col = "darkgray", lty = 3, data.panel = "all")

  dev.off()
}

# call functions
####################

# read tsv file into df for BAF plot
df_trimmed <- read_to_df(VCF_FILE)

# get quantiles for plotting limits
MAX_DEPTH <- quantile(df_trimmed$Depth, probs = MAX_DEPTH, names = FALSE)

# Filter out low depth rows
df_filtered <- df_trimmed[df_trimmed$Depth >= MIN_DEPTH, ]

# make tsvs for testing if needed
base <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(VCF_FILE)))
write.table(df_filtered, file=paste0(base, ".baf.tsv"), quote=FALSE, sep='\t', col.names = NA)

# dynamic bin size choice
BIN_SIZE <- ifelse(length(df_filtered$Depth)/3000 >= 1, round(length(df_filtered$Depth)/3000), 1)

# read bed file into binned df for depth plot
df_binned <- bin_df(df_trimmed, BIN_SIZE)

# convert dfs for baf plotting into snp.data for karyoploter
snp.data.baf <- get_snp_data_BAF(df_filtered)

# convert dfs for depth plotting into snp.data for karyoploter
snp.data.depth <- get_snp_data_Depth(df_binned)

# generate plots and save them
get_plot(snp.data.baf, snp.data.depth, VCF_FILE, MAX_DEPTH, CHR_NAMES, MIN_BAF, MAX_BAF, GENOME)

