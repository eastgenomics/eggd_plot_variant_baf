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
    bcftools query -f '%CHROM\t%POS\t%INFO/DP\t[ %AD]\n' $gvcf_path -o "$gvcf_prefix.gvcf.tsv"

    # construct optional argument string
    options=""
    options+="--min_baf $min_baf "
    options+="--max_baf $max_baf "
    options+="--max_depth $max_depth "
    options+="--min_depth $min_depth "
    options+="--bin_size $bin_size "
    options+="--chr_names $chr_names "
    options+="--genome $genome "
    options+="--symmetry $symmetry "


    # Run R script with error handling
    if ! Rscript baf_depth_plotting.R --vcf "$vcf_prefix.vcf.tsv" --gvcf "$gvcf_prefix.gvcf.tsv" $options; then
    echo "Error: BAF plotting failed with exit code $?" >&2
    exit 1
    fi
    
    # Deal with output
    mkdir -p out/baf_plot
    mv *.png out/baf_plot

    if [ "$output_tsv" = true ]; then
        mkdir -p out/tsv
        mv *.baf.tsv out/tsv
    fi

    dx-upload-all-outputs

    echo "Done."
}