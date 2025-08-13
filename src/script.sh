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

    # Derive prefixes for output filenames 
    vcf_prefix=$(basename "$vcf_path"  | sed -E 's/\.vcf(\.gz)?$//')
    gvcf_prefix=$(basename "$gvcf_path"| sed -E 's/\.vcf(\.gz)?$//')

    # Adapt region strings based on contig style
    vcf_regions="$chr_names"
    if bcftools view -h "$vcf_path" | grep -q '##contig=<ID=chr'; then
        vcf_regions=$(awk -v RS=, -v ORS=, '{print "chr"$0}' <<< "$chr_names")
        vcf_regions="${vcf_regions%,}"
    fi

    gvcf_regions="$chr_names"
    if bcftools view -h "$gvcf_path" | grep -q '##contig=<ID=chr'; then
        gvcf_regions=$(awk -v RS=, -v ORS=, '{print "chr"$0}' <<< "$chr_names")
        gvcf_regions="${gvcf_regions%,}"
    fi

    echo "gvcf_regions = $gvcf_regions"

    # Compress and index if needed
    for var in vcf_path gvcf_path; do
        path="${!var}"
        if [[ "$path" != *.vcf.gz ]]; then
            bgzip -c "$path" > "${path}.gz"
            printf -v "$var" '%s' "${path}.gz"
        fi
        bcftools index -t "${!var}"
    done

    # Query VCF for CHROM, POS, depth and allele counts
    bcftools query -r "$vcf_regions" \
        -f '%CHROM\t%POS\t%INFO/DP\t[ %AD]\n' \
        "$vcf_path" -o "${vcf_prefix}.vcf.tsv"

    # Query gVCF with fallback logic FORMAT/MIN_DP -> FORMAT/DP ->  INFO/DP
    bcftools query -u -r "$gvcf_regions" \
        -f '%CHROM\t%POS\t[%MIN_DP]\t[%DP]\t%INFO/DP\n' "$gvcf_path" | \
    awk -F'\t' 'BEGIN{OFS="\t"}{
        dp=$3; 
        if(dp=="."||dp=="") dp=$4; 
        if(dp=="."||dp=="") dp=$5; 
        print $1,$2,dp
    }' > "${gvcf_prefix}.gvcf.tsv"

    # Fail early if gVCF TSV is empty
    if [ ! -s "${gvcf_prefix}.gvcf.tsv" ]; then
        echo "Error: gVCF TSV is empty — check region string or file contents" >&2
        exit 1
    fi


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