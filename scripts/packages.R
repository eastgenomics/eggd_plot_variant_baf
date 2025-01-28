libPath <- "~/R/library"

if (!dir.exists(libPath)) {
  dir.create(libPath, recursive = TRUE)
}

.libPaths(libPath)

# Install CRAN packages with explicit versions
install.packages("stringr", version = "1.5.1", lib = libPath, repos = "https://cloud.r-project.org")
install.packages("dplyr", version = "1.1.4", lib = libPath, repos = "https://cloud.r-project.org")

# Install BiocManager if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = libPath, repos = "https://cloud.r-project.org")
}


# Install karyoploteR using Bioconductor v3.18 (latest in R v4.3)
# https://bioconductor.org/news/bioc_3_18_release/
# BiocManager installs the version of packages that was released in the specified bioconductor release
BiocManager::install(
  "karyoploteR",
  version = "3.18",  # Bioconductor release version (not package version); installs karyoploterR v1.28.0
  lib = libPath,
)

# Install polars from GitHub release
install.packages(
  "https://github.com/pola-rs/r-polars/releases/download/v0.22.0/polars__x86_64-pc-linux-gnu.gz",
  repos = NULL,
  lib = libPath
)