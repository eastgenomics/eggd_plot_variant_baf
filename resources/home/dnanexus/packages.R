libPath <- "~/R/library"

if (!dir.exists(libPath)) {
  dir.create(libPath, recursive = TRUE)
}

install.packages("stringr", version = "1.5.1", lib = libPath, repos = "https://cloud.r-project.org")
install.packages("dplyr", version = "1.1.4", lib = libPath, repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = libPath, repos = "https://cloud.r-project.org")
}

BiocManager::install("karyoploteR", version = "1.28.0", lib = libPath)

install.packages(
  "https://github.com/pola-rs/r-polars/releases/download/0.22.0/polars__x86_64-pc-linux-gnu.gz",
  repos = NULL,
  lib = libPath
)