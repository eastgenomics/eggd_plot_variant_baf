#!/bin/bash

# prefixes all lines of commands written to stdout with datetime
PS4='\000[$(date)]\011'
export TZ=Europe/London
set -exo pipefail

# set frequency of instance usage in logs to 10 seconds
kill $(ps aux | grep pcp-dstat | head -n1 | awk '{print $2}')
/usr/bin/dx-dstat 10

main() {

    dx-download-all-inputs

    sudo dpkg -i libtinfo5_6.2-0ubuntu2_amd64.deb
    sudo dpkg -i libncurses5_6.2-0ubuntu2_amd64.deb
    sudo dpkg -i libssl1.1_1.1.1f-1ubuntu2_amd64.deb

    tar -xzf $packages_path
    echo "R_LIBS_USER=~/R/library" >> ~/.Renviron

    bcftools query -f '%CHROM\t%POS\t%INFO/DP\t[ %AD]\n' $vcf_path -o "$vcf_prefix.vcf.tsv"

    if ! Rscript baf_depth_plotting.R; then
    echo "Error: BAF plotting failed with exit code $?" >&2
    exit 1
    fi
    
    mkdir -p out/baf_plot
    mv *.png out/baf_plot

    dx-upload-all-outputs

    echo "Done."
}