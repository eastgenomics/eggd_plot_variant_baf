# eggd_plot_variant_baf

## What does this app do?
Plots BAF and depth of variants from given VCF/GVCF pair. The BAF plot is the top panel and the Depth plot is the bottom panel.


## What data are required for this app to run?
**Required input files:**
1. A VCF file (`.vcf`) - Provides allele counts (AD) and depth (DP) for the BAF plot.
2. A GVCF file (`.gvcf`) - containing the depth across all sites for the depth plot.
3. R packages (`.tar.gz`) - compressed tarball of R packages needed to generate plot.
<br>

### Optional input Parameters:

#### Parameters for BAF plot:
- `min_baf` : The minimum BAF for the baf plot. Default is 0.04
- `max_baf` : The maximum BAF for the baf plot. Default is 0.96
- `min_depth` : The minimum depth for the BAF plot. Default is 5
- `symmetry` : Input to plot symmetrical points on the BAF plot. Default is true.

#### Parameters for Depth plot:
- `bin_size` : The bin size to use for plotting depth.
- `max_depth` : The maximum depth for the plots used to set the y-axis. Points above this cut are drawn in magenta. Default is 0.9.


#### Common Parameters for both plots:
- `bed_filter` : The BED file used to filter the input VCFs.
- `genome` : Genome build for plotting. Default is hg19, alternative value is hg38.
- `min_qual` : The minimum QUAL of variants for plotting. Default is 0
- `chr_names` : Chromosome names used for axis labels and to filter VCF contigs.
- `output_tsv` : Whether to output a TSV file of the BAF dataframe for testing purposes. Default is False.

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

### If output_tsv is set to True:
- `prefix`.vcf.baf.tsv : tsv that contains the fields `Chr, Position, Depth, Ref_AD, Alt_AD, RAF, BAF`
- `prefix`.filtered.baf.tsv tsv that contains the fields `Chr, Position, Depth, Ref_AD, Alt_AD, RAF, BAF` with variants with Depth >= min_depth and between `min_baf` and `max_baf`.
- `prefix`.gcf.baf.tsv : Contains the Chr, Position and Depth from the GVCF input.
- `prefix`.binned.baf.tsv : Contains the binned version of GVCF depth data. Contains `Chr, Position, mean_depth`





## How to run this app from command line?
```
dx run eggd_plot_variant_baf \
-ivcf=file-xxxx \
-igvcf=file-xxxx \
-imax_depth=0.98 \
-imin_depth=20 \
-imin_baf=0 \
-imax_baf=1
-ipackages=file-xxxx  \
--destination="output/eggd_plot_variant_baf"
```
