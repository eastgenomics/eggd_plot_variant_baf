#!/bin/bash
set -exo pipefail

main() {

    sudo apt-get install -y bcftools

    #Rscript packages.R
    tar -xzvf R_packages.tar.gz
   

    dx-download-all-inputs

    bcftools query -f '%CHROM\t%POS\t%INFO/DP\t[ %AD]\n' $vcf_path -o "$vcf_prefix.vcf.bed"

    Rscript bafs_plot.R
    
    mkdir -p out/baf_plot
    mv *.png out/baf_plot

    dx-upload-all-outputs

    echo "Done."
}