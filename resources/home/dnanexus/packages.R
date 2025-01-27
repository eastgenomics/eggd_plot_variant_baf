libPath <- "~/R/library"

if (!dir.exists(libPath)) {
  dir.create(libPath, recursive = TRUE)
}

.libPaths(libPath)
install.packages(c("stringr", "dplyr"), lib = libPath)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = libPath)
}
BiocManager::install("karyoploteR", lib = libPath)

install.packages(
  "https://github.com/pola-rs/r-polars/releases/latest/download/polars__x86_64-pc-linux-gnu.gz",
  repos = NULL,
  lib = libPath
)

