# eggd_plot_variant_baf

## What does this app do?
Plots BAF and depth of variants from given VCF/GVCF pair.


## What data are required for this app to run?
**Required input files:**
1. A VCF file (`.vcf`) - containing the VEP filtered variants for the BAF plot.
2. A GVCF file (`.gvcf`) - containing the merged small variants for the depth plot.
3. R packages (`.tar.gz`) - compressed tarball of R packages needed to generate plot.
<br>

**R Packages and Versions:**

- `stringr` (v1.5.1)
- `dplyr` (v1.1.4)
- `karyoploteR` (v1.28.0)
- `polars` (v0.22.0)
- `argparse` (v2.2.5)

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
-igvcf=file-xxxx \
-imax_depth=0.98 \
-imin_depth=20 \
-imin_baf=0 \
-imax_baf=1
-ibin_size=500
-ipackages=file-xxxx  \
--destination="output/eggd_plot_variant_baf"
```

## Notes
The current version must provide a value to the `bin_size` option as automatic scaling is not yet functional.