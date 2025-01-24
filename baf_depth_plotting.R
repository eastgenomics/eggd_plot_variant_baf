# SCRIPT TO PLOT BAF AND DEPTH WITH KARYOPLOTER
#
# @input - list of files in tsv format with "CHROM", "POS", "DP" and "AD"
# @output - png plot of BAF and binned Depth values for each file
#
# Constraints
##################
# The genome parameter in the plotkaryotype function defaults to hg19
# Depth values that are higher than the max value in Y axis will be plotted in the BAF plot
# Only .tsv files can be provided
# Chromosome names need to be provided as it defaults to include "chr"

# Configurables
##################

# Minimum and maximum thresholds of BAFs to be plotted
MIN_BAF <- 0.04
MAX_BAF <- 0.96

# Variants with depth below this threshold are excluded
MIN_DEPTH <- 50

# Used to aggregate depth values for smoother visualization
BIN_SIZE <- 1000

# Adjusts the Y-axis plot for mean_depth values
MAX_DEPTH_PLOT <- 750

# Chromosome labels to feature in the plot X-axis
CHR_NAMES <- c(paste0(1:22), "X", "Y")

# Load required modules
##################################
library(stringr, quietly = TRUE)
library(karyoploteR, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(polars, quietly = TRUE)

# List of functions
##################################

# Function to read files into dfs
# @parameter file
# @parameter min_depth, minimum depth for variant filtering
# returns trimmed_df

read_to_df <- function(file, min_depth = MIN_DEPTH) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  df <- read.table(file = file, header = FALSE)
  if (ncol(df) != 4) {
    stop("Invalid TSV format: Expected 4 columns (CHROM, POS, DP, AD)")
  }
  colnames(df) <- c("Chr", "Position", "Depth", "Allele_Depth")
  if (!all(sapply(df[c("Position", "Depth")], is.numeric))) {
    stop("Invalid TSV format: Position and Depth must be numeric")
  }

  # calculate BAF
  df[c("Ref_AD", "Alt_DP")] <- str_split_fixed(df$Allele_Depth, ",", 2)
  df$RAF <- as.numeric(df$Ref_AD)
  df$BAF <- as.numeric(df$Alt_DP)
  # Avoid division by zero
  df$RAF <- ifelse(df$Depth > 0, as.numeric(df$Ref_AD) / df$Depth, NA)
  df$BAF <- ifelse(df$Depth > 0, as.numeric(df$Alt_DP) / df$Depth, NA)
  df_trimmed <- df[df$Depth > min_depth, ]

  return(df_trimmed)
}

# Function to return dfs with binned depth
# @parameter df
# @parameter bin_size - integer: length of window for variant aggregation
# returns df_binned

bin_df <- function(df, bin_size = BIN_SIZE) {
  polars_df <- as_polars_df(df)
  rolling_df <- polars_df$rolling(
    index_column = "Position",
    group_by = "Chr",
    period = paste(as.character(bin_size), "i", sep = "")
  )$agg(mean_depth = pl$mean("Depth"))$gather_every(bin_size)
  df_binned <- rolling_df$to_data_frame()

  return(df_binned)
}

# Function to create snp.data.baf format as input to karyoploteR plot BAF
# @parameter df, use a trimmed df
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
# @parameters snp.data.bf, snp.data.depth, file_name
# @parameter max_depth_plot - integer : maximum depth value for y-axis in plots
# @parameter chr_names - vector with chromosomes to be listed
# @parameters min_baf and max_baf - integers : include only variants in range 0.4 < BAF < 0.96
# returns plot

get_plot <- function(snp.data.baf, snp.data.depth, file_name, max_depth_plot = MAX_DEPTH_PLOT, chr_names = CHR_NAMES, min_baf = MIN_BAF, max_baf = MAX_BAF) {
  file_name_png <- paste0(sub(".tsv", "", file_name), ".png")
  png(file_name_png, width = 15, height = 5, units = "in", res = 600)
  plot_parameters <- getDefaultPlotParams(plot.type = 4)
  plot_parameters$data1inmargin <- 2
  baf_depth_plot <- plotKaryotype(plot.type = 4, ideogram.plotter = NULL, plot.params = plot_parameters, labels.plotter = NULL)
  kpAddChromosomeNames(baf_depth_plot, chr.names = chr_names)
  kpAddCytobandsAsLine(baf_depth_plot) # Add centromers
  # top graph
  baf_threshold <- which(snp.data.baf$BAF > min_baf & snp.data.baf$BAF < max_baf)
  kpAxis(baf_depth_plot, r0 = 0.55, r1 = 1, tick.pos = c(0, 0.25, 0.5, 0.75, 1))
  kpPoints(baf_depth_plot,
    data = snp.data.baf[baf_threshold], y = snp.data.baf[baf_threshold]$BAF,
    cex = 0.5, r0 = 0.55, r1 = 1, col = "darkorange2"
  )
  # bottom graph
  kpAxis(baf_depth_plot, r0 = 0, r1 = 0.45, ymax = max_depth_plot, ymin = 0)
  kpPoints(baf_depth_plot,
    data = snp.data.depth, y = snp.data.depth$mean_depth,
    cex = 0.5, r0 = 0, r1 = 0.45, ymax = max_depth_plot, ymin = 0, col = "darkblue"
  )
  kpAddMainTitle(baf_depth_plot, main = "BAF vs Depth")
  kpAddChromosomeSeparators(baf_depth_plot, col = "darkgray", lty = 3, data.panel = "all")

  dev.off()
}

# call functions
####################

# list bed files for plotting
gvcf_files <- list.files(path = ".", pattern = ".tsv")

# read bed files into dfs for BAF plot
df_list <- list()

for (file in gvcf_files) {
  df_trimmed <- read_to_df(file)
  df_list[[length(df_list) + 1]] <- df_trimmed
}

# read bed files into binned dfs for depth plot
df_binned_list <- list()

for (df in df_list) {
  df_binned <- bin_df(df)
  df_binned_list[[length(df_binned_list) + 1]] <- df_binned
}

# convert dfs for baf plotting into snp.data for karyoploter
snp.data.baf_list <- list()
for (df in df_list) {
  snp.data.baf <- get_snp_data_BAF(df)
  snp.data.baf_list[[length(snp.data.baf_list) + 1]] <- snp.data.baf
}

# convert dfs for depth plotting into snp.data for karyoploter
snp.data.depth_list <- list()
for (df in df_binned_list) {
  snp.data.depth <- get_snp_data_Depth(df)
  snp.data.depth_list[[length(snp.data.depth_list) + 1]] <- snp.data.depth
}

# generate plots and save them
mapply(function(snp.data.baf, snp.data.depth, file_name) {
  get_plot(snp.data.baf, snp.data.depth, file_name)
}, snp.data.baf_list, snp.data.depth_list, gvcf_files)