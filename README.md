# eggd_plot_variant_baf

## What does this app do?
Plots BAF and depth of variants from given VCF


## What data are required for this app to run?
**Required input files:**
1. A VCF file (`.vcf`) - containing the merged small variants we want to plot.
2. R packages (`.tar.gz`) - compressed tarball of R packages needed to generate plot.
<br>
**R Packages and Versions:**

- `stringr` (v1.5.1)
- `dplyr` (v1.1.4)
- `karyoploteR` (v1.28.0)
- `polars` (v0.22.0)

**How to build the package**
The package was built on Ubuntu 24.04 and R v4.3. Below are the steps:
1. Update and install the required dependencies:
```
sudo apt-get update
sudo apt-get install -y libssl-dev libxml2-dev gcc pkg-config

```
2. Run the R script to install the required R packages:
`Rscript scripts/packages.R`

3. Compress the R library folder:
`tar -czvf R_packages.tar.gz R/library`


## What does this app output?
This app outputs:
- `{prefix}.png` : Image of the generated plot in PNG format.


## How to run this app from command line?
```
dx run eggd_plot_variant_baf \
-ivcf=file-xxxx \
-ipackages=file-xxxx  \
--destination="output/eggd_plot_variant_baf"
```

## Notes
The current version of the plotting has the following constraints:
- plotting of variant depth is currently hard limited at DP<50, has a upper Y axis limit of 750 and plots the mean depth across 1000 consecutive variants
- plotting of the BAF is hard limited at < 0.04 and > 0.96