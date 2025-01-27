# eggd_plot_variant_baf

## What does this app do?
Plots BAF and depth of variants from given VCF


## What data are required for this app to run?
Required input files:
1. A VCF file (`.vcf`) - containing the merged small variants we want to plot.
2. R packages (`.tar.gz`) - compressed tarball of R packages needed to generate plot.


## What does this app output?
This app outputs:
- `{prefix}.png` : Image of the generated plot in PNG format.


## How to run this app from command line?
```
dx run eggd_mosdepth \
-ivcf=file-xxxx \
-ipackages=file-xxxx  \
--destination="output/eggd_plot_variant_baf"
```

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

#### This app was made by EMEE GLH
