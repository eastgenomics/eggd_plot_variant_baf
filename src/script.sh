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

    # Ensure input chromosome names do not have 'chr' prefix
    chr_names_normalized=$(echo "$chr_names" | awk -v RS=, -v ORS=, '{sub(/^chr/, ""); print}' | sed 's/,$//')

    # Adapt region strings based on contig style
    vcf_regions="$chr_names_normalized"
    if bcftools view -h "$vcf_path" | grep -q '##contig=<ID=chr'; then
        vcf_regions=$(awk -v RS=, -v ORS=, '{print "chr"$0}' <<< "$chr_names_normalized")
        vcf_regions="${vcf_regions%,}"
    fi

    gvcf_regions="$chr_names_normalized"
    if bcftools view -h "$gvcf_path" | grep -q '##contig=<ID=chr'; then
        gvcf_regions=$(awk -v RS=, -v ORS=, '{print "chr"$0}' <<< "$chr_names_normalized")
        gvcf_regions="${gvcf_regions%,}"
    fi

    # Remove all trailing newlines from region strings
    gvcf_regions=$(echo "$gvcf_regions" | tr -d '\n')

    # Compress and index if needed
    for var in vcf_path gvcf_path; do
        path="${!var}"
        if [[ "$path" != *.gz ]]; then
            bgzip -c "$path" > "${path}.gz"
            printf -v "$var" '%s' "${path}.gz"
        fi
        bcftools index -t "${!var}"
    done

    if [[ "$min_qual" -gt 0 ]]; then
        echo "Filtering VCF: keeping QUAL >= ${min_qual}"
        bcftools view -i "QUAL>=${min_qual}" "$vcf_path" -Oz -o tmp.vcf.gz && mv tmp.vcf.gz "$vcf_path"
        bcftools index -f -t "$vcf_path"

        echo "Filtering gVCF: keeping QUAL >= ${min_qual} for variants, all ref blocks"
        # This keeps the REF blocks intact while filtering variants
        bcftools view -i "(ALT != \".\" && QUAL>=${min_qual}) || ALT = \".\"" "$gvcf_path" -Oz -o tmp.gvcf.gz && mv tmp.gvcf.gz "$gvcf_path"
        bcftools index -f -t "$gvcf_path"
    fi

    # Query VCF for CHROM, POS, depth and allele counts
    bcftools query -r "$vcf_regions" \
        -f '%CHROM\t%POS\t%INFO/DP\t[ %AD]\n' \
        "$vcf_path" -o "${vcf_prefix}.vcf.tsv"

    # Query gVCF with fallback logic FORMAT/DP ->  INFO/DP
    bcftools query -u -r "$gvcf_regions" \
        -f '%CHROM\t%POS\t[%DP]\t%INFO/DP\n' "$gvcf_path" | \
    awk -F'\t' 'BEGIN{OFS="\t"}{
        dp=$3; 
        if(dp=="."||dp=="") dp=$4; 
        print $1,$2,dp
    }' > "${gvcf_prefix}.gvcf.tsv"

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