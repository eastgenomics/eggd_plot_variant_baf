#.libPaths('~/R/library')
library(stringr)
library(karyoploteR)
library(biomaRt)
library(dplyr)
library("polars")
polars_code_completion_activate()

#read_files
gvcf_files <- list.files(path=".", pattern="genome.vcf.bed")


print(gvcf_files)


read_to_df <- function(file) {
  
  df <- read.table(file = file, header=FALSE)
  colnames(df) <- c("Chr", "Position", "Depth", "Allele_Depth")
  df$Chr<-sub("chr","",df$Chr)
  df[c('Ref_AD', 'Alt_DP', 'Alt_DP2')] <- str_split_fixed(df$Allele_Depth, ',', 3)
  df$Ref_AD <- as.numeric(df$Ref_AD)
  df$Alt_DP <- as.numeric(df$Alt_DP)
  df$Alt_DP2 <- as.numeric(df$Alt_DP2)
  df$RAF <- df$Ref_AD/df$Depth
  df$BAF <- df$Alt_DP/df$Depth
  df$BAF2 <- df$Alt_DP2/df$Depth
  df$log2 <- log2(df$Depth)
  df_trimmed <- df[df$Depth > 50, ]
  
  return(df_trimmed)
  
}

df_list <- list()

for (file in gvcf_files) {
  df_trimmed <- read_to_df(file)
  df_list[[length(df_list) + 1]] <- df_trimmed 
}

print(df_list)

#

bin_df <- function(df) {
  
  polars_df <- as_polars_df(df)
  rolling_df <-  polars_df$rolling(
    index_column='Position', 
    group_by='Chr',
    period="1000i")$agg(mean_depth = pl$mean('Depth'))$gather_every(1000)
  
  df_binned <- rolling_df$to_data_frame()
  
  return(df_binned)
  
}



df_binned_list <- list()

for (df in df_list) {
  df_binned <- bin_df(df)
  df_binned_list[[length(df_binned_list) + 1]] <- df_binned 
}


print(df_binned_list)
  
#snpdata object for each df

get_snp_data_BAF <- function(df) {

  snp.data.baf <- toGRanges(df[,c("Chr", "Position", "Position", "Depth", "BAF")])
  names(mcols(snp.data.baf)) <- c("Depth", "BAF")
  seqlevelsStyle(snp.data.baf) <- "UCSC"

  
  return(snp.data.baf)
  
}

get_snp_data_Depth <- function(df) {
  snp.data.depth <- toGRanges(df[,c("Chr", "Position", "Position", "mean_depth")])
  names(mcols(snp.data.depth)) <- c("mean_depth")
  seqlevelsStyle(snp.data.depth) <- "UCSC"
  
  return(snp.data.depth)
  
}

snp.data.baf_list <- list()
for (df in df_list) {
  snp.data.baf <- get_snp_data_BAF(df)
  snp.data.baf_list[[length(snp.data.baf_list) + 1]] <- snp.data.baf 
}



snp.data.depth_list <- list()
for (df in df_binned_list) {
  snp.data.depth <- get_snp_data_Depth(df)
  snp.data.depth_list[[length(snp.data.depth_list) + 1]] <- snp.data.depth 
}

get_plot <- function(snp.data.baf, snp.data.depth, file_name) {
  
  file_name_png <- paste0(sub(".bed", "", file_name), ".png")
  png(file_name_png, width = 2500, height = 750)
  
  pp <- getDefaultPlotParams(plot.type = 4)
  pp$data1inmargin <- 2
  
  kp <- plotKaryotype(plot.type = 4, ideogram.plotter = NULL, plot.params = pp)
  
  kpAddCytobandsAsLine(kp)
  kpAxis(kp, r0=0.55, r1=1)
  kpPoints(kp, data=snp.data.baf, y=snp.data.baf$BAF, cex=0.3, r0=0.55, r1=1, col='darkgrey')
  kpAxis(kp, r0=0, r1=0.45, ymax=750, ymin = 5)
  kpPoints(kp, data=snp.data.depth, y=snp.data.depth$mean_depth, cex=0.3, r0=0, r1=0.45, ymax=750, ymin=0, col="darkblue")
  kpAddMainTitle(kp, main='BAF vs Depth')
  kpAddChromosomeSeparators(kp, col="gray", lty=3, data.panel="all")
  
  dev.off()
  
}

mapply(function(snp.data.baf, snp.data.depth, file_name) {
  get_plot(snp.data.baf, snp.data.depth, file_name)
}, snp.data.baf_list, snp.data.depth_list, gvcf_files)
