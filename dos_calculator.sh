#!/bin/bash
#
# Calculate per-gene Direction of Selection (DoS) statistic for each population
#
# DoS = (Dn / (Dn + Ds)) - (Pn / (Pn + Ps))
#
# DoS > 0 suggests positive (Darwinian) selection
# DoS < 0 suggests negative (purifying) selection

set -euo pipefail

# Default values
MIN_SITES=4
OUTPUT="dos_results.tsv"

# Usage function
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required:
  -n, --nonsyn-vcf FILE      VCF file with nonsynonymous variants
  -s, --syn-vcf FILE         VCF file with synonymous variants
  -p, --populations FILE     Tab-delimited: population_name sample1,sample2,...
  
Optional:
  -o, --outgroup SAMPLES     Comma-separated list of outgroup sample names
  -O, --output FILE          Output file name (default: dos_results.tsv)
  -m, --min-sites INT        Minimum Pn + Dn for inclusion (default: 4)
  -g, --gene-field NAME      INFO field name for gene (default: GENE)
  -h, --help                 Show this help message

Example:
  $0 -n nonsyn.vcf.gz -s syn.vcf.gz -p populations.txt -o outgroup1,outgroup2

Population file format (tab-delimited):
  POP1	sample1,sample2,sample3
  POP2	sample4,sample5,sample6
EOF
    exit 1
}

# Parse arguments
NONSYN_VCF=""
SYN_VCF=""
POP_FILE=""
OUTGROUP=""
GENE_FIELD="GENE"

while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--nonsyn-vcf)
            NONSYN_VCF="$2"
            shift 2
            ;;
        -s|--syn-vcf)
            SYN_VCF="$2"
            shift 2
            ;;
        -p|--populations)
            POP_FILE="$2"
            shift 2
            ;;
        -o|--outgroup)
            OUTGROUP="$2"
            shift 2
            ;;
        -O|--output)
            OUTPUT="$2"
            shift 2
            ;;
        -m|--min-sites)
            MIN_SITES="$2"
            shift 2
            ;;
        -g|--gene-field)
            GENE_FIELD="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check required arguments
if [[ -z "$NONSYN_VCF" ]] || [[ -z "$SYN_VCF" ]] || [[ -z "$POP_FILE" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check if required tools are available
for tool in bcftools awk; do
    if ! command -v $tool &> /dev/null; then
        echo "Error: $tool is required but not installed"
        exit 1
    fi
done

echo "=== DoS Calculator ===" >&2
echo "Nonsynonymous VCF: $NONSYN_VCF" >&2
echo "Synonymous VCF: $SYN_VCF" >&2
echo "Populations file: $POP_FILE" >&2
echo "Outgroup: ${OUTGROUP:-none}" >&2
echo "Min sites: $MIN_SITES" >&2
echo "" >&2

# Create temporary directory
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Read population assignments
declare -A POP_SAMPLES
while IFS=$'\t' read -r pop samples; do
    POP_SAMPLES[$pop]=$samples
done < "$POP_FILE"

echo "Found ${#POP_SAMPLES[@]} populations" >&2

# Function to extract gene from VCF INFO field
extract_gene() {
    local info=$1
    echo "$info" | grep -oP "${GENE_FIELD}=\K[^;]+" || echo "NA"
}

# Function to count variants for a population
count_variants() {
    local vcf=$1
    local pop=$2
    local samples=$3
    local outgroup_samples=$4
    local output_file=$5
    
    # Get sample column indices
    local header=$(bcftools view -h "$vcf" | tail -1)
    local sample_cols=""
    
    IFS=',' read -ra SAMPLE_ARRAY <<< "$samples"
    for sample in "${SAMPLE_ARRAY[@]}"; do
        local col=$(echo "$header" | tr '\t' '\n' | grep -n "^${sample}$" | cut -d: -f1)
        if [[ -n "$col" ]]; then
            sample_cols="$sample_cols $col"
        fi
    done
    
    # Get outgroup column indices
    local outgroup_cols=""
    if [[ -n "$outgroup_samples" ]]; then
        IFS=',' read -ra OUT_ARRAY <<< "$outgroup_samples"
        for sample in "${OUT_ARRAY[@]}"; do
            local col=$(echo "$header" | tr '\t' '\n' | grep -n "^${sample}$" | cut -d: -f1)
            if [[ -n "$col" ]]; then
                outgroup_cols="$outgroup_cols $col"
            fi
        done
    fi
    
    # Debug output
    echo "    Sample columns: $sample_cols" >&2
    echo "    Outgroup columns: $outgroup_cols" >&2
    
    # Count total variants in VCF
    local total_vars=$(bcftools view -H "$vcf" | wc -l)
    echo "    Total variants in VCF: $total_vars" >&2
    
    # Process VCF
    bcftools view -H "$vcf" | awk -v pop="$pop" -v gene_field="$GENE_FIELD" \
        -v sample_cols="$sample_cols" -v outgroup_cols="$outgroup_cols" '
    BEGIN {
        split(sample_cols, scols, " ")
        split(outgroup_cols, ocols, " ")
        
        # Remove empty elements
        j = 0
        for (i in scols) {
            if (scols[i] != "") {
                j++
                sample_idx[j] = scols[i]
            }
        }
        n_samples = j
        
        j = 0
        for (i in ocols) {
            if (ocols[i] != "") {
                j++
                outgroup_idx[j] = ocols[i]
            }
        }
        n_outgroup = j
        
        variants_processed = 0
        variants_with_gene = 0
    }
    {
        variants_processed++
        
        # Extract gene name from INFO field
        gene = "NA"
        
        # Try SnpEff ANN field (newer format) - most common
        # ANN=A|missense_variant|MODERATE|gene_name|gene_id|...
        if ($8 ~ /ANN=/) {
            match($8, /ANN=[^;]*/)
            ann_str = substr($8, RSTART+4, RLENGTH-4)
            # Handle multiple annotations separated by comma
            split(ann_str, ann_entries, ",")
            # Take first annotation
            split(ann_entries[1], ann_parts, "|")
            if (length(ann_parts) >= 4 && ann_parts[4] != "") {
                gene = ann_parts[4]
            }
        }
        
        # Try standard GENE field
        if (gene == "NA" && $8 ~ gene_field "=") {
            n = split($8, info_fields, ";")
            for (i = 1; i <= n; i++) {
                if (info_fields[i] ~ "^" gene_field "=") {
                    split(info_fields[i], kv, "=")
                    gene = kv[2]
                    break
                }
            }
        }
        
        # Try SnpEff EFF field (older format)
        # EFF=missense_variant(MODERATE|MISSENSE|Gca/Aca|A142T|gene_name|...
        if (gene == "NA" && $8 ~ /EFF=/) {
            match($8, /EFF=[^;]*/)
            eff_str = substr($8, RSTART+4, RLENGTH-4)
            # Parse: effect(impact|type|aa_change|aa_pos|gene_name|...
            if (match(eff_str, /\([^)]+\)/)) {
                paren_content = substr(eff_str, RSTART+1, RLENGTH-2)
                split(paren_content, eff_parts, "|")
                if (length(eff_parts) >= 5) {
                    gene = eff_parts[5]
                }
            }
        }
        
        if (gene == "NA" || gene == "") next
        variants_with_gene++
        
        # Check polymorphism in population
        has_ref = 0
        has_alt = 0
        
        for (i = 1; i <= n_samples; i++) {
            col = sample_idx[i]
            gt = $col
            split(gt, gt_parts, ":")
            genotype = gt_parts[1]
            
            if (genotype ~ /0/) has_ref = 1
            if (genotype ~ /1/) has_alt = 1
        }
        
        is_polymorphic = (has_ref && has_alt) ? 1 : 0
        
        # Check divergence with outgroup
        # Divergent = variant present in population but absent in outgroup
        # (doesn not need to be fixed in population, just present)
        is_divergent = 0
        if (n_outgroup > 0 && has_alt) {
            out_has_alt = 0
            
            for (i = 1; i <= n_outgroup; i++) {
                col = outgroup_idx[i]
                gt = $col
                split(gt, gt_parts, ":")
                genotype = gt_parts[1]
                
                # Check if outgroup has the alternate allele
                if (genotype ~ /1/) out_has_alt = 1
            }
            
            # Divergent if alternate allele present in population but not in outgroup
            if (!out_has_alt) {
                is_divergent = 1
            }
        }
        
        # Output: population, gene, is_polymorphic, is_divergent
        print pop "\t" gene "\t" is_polymorphic "\t" is_divergent
    }
    END {
        print "DEBUG: Processed " variants_processed " variants, " variants_with_gene " had gene annotations" > "/dev/stderr"
    }' > "$output_file"
    
    # Debug: count output lines
    local output_lines=$(wc -l < "$output_file")
    echo "    Output lines: $output_lines" >&2
}

# Process each population
for pop in "${!POP_SAMPLES[@]}"; do
    echo "Processing population: $pop" >&2
    samples="${POP_SAMPLES[$pop]}"
    
    # Count nonsynonymous variants
    echo "  - Counting nonsynonymous variants..." >&2
    count_variants "$NONSYN_VCF" "$pop" "$samples" "$OUTGROUP" "$TMPDIR/${pop}_nonsyn.txt"
    
    # Count synonymous variants
    echo "  - Counting synonymous variants..." >&2
    count_variants "$SYN_VCF" "$pop" "$samples" "$OUTGROUP" "$TMPDIR/${pop}_syn.txt"
done

# Combine and calculate DoS
echo "" >&2
echo "Calculating DoS statistics..." >&2

{
    echo -e "Population\tGene\tDn\tDs\tPn\tPs\tInformative_Sites\tDoS"
    
    for pop in "${!POP_SAMPLES[@]}"; do
        # Aggregate counts by gene
        cat "$TMPDIR/${pop}_nonsyn.txt" "$TMPDIR/${pop}_syn.txt" | \
        awk -v pop="$pop" '
        {
            gene = $2
            is_poly = $3
            is_div = $4
            
            # Check if from nonsyn or syn based on file
            if (FILENAME ~ /nonsyn/) {
                if (is_poly) Pn[gene]++
                if (is_div) Dn[gene]++
            } else {
                if (is_poly) Ps[gene]++
                if (is_div) Ds[gene]++
            }
            genes[gene] = 1
        }
        END {
            for (gene in genes) {
                pn = (gene in Pn) ? Pn[gene] : 0
                dn = (gene in Dn) ? Dn[gene] : 0
                ps = (gene in Ps) ? Ps[gene] : 0
                ds = (gene in Ds) ? Ds[gene] : 0
                
                print pop "\t" gene "\t" dn "\t" ds "\t" pn "\t" ps
            }
        }' "$TMPDIR/${pop}_nonsyn.txt" "$TMPDIR/${pop}_syn.txt"
    done
} | awk -v min_sites="$MIN_SITES" '
NR == 1 {
    print
    next
}
{
    pop = $1
    gene = $2
    Dn = $3
    Ds = $4
    Pn = $5
    Ps = $6
    
    informative = Pn + Dn
    
    # Calculate DoS
    if ((Dn + Ds) > 0 && (Pn + Ps) > 0) {
        dos = (Dn / (Dn + Ds)) - (Pn / (Pn + Ps))
    } else {
        dos = "NA"
    }
    
    # Apply filter
    if (informative >= min_sites && dos != "NA") {
        print pop "\t" gene "\t" Dn "\t" Ds "\t" Pn "\t" Ps "\t" informative "\t" dos
    }
}' > "$OUTPUT"

echo "Results saved to: $OUTPUT" >&2

# Print summary statistics
echo "" >&2
echo "=== SUMMARY ===" >&2

for pop in "${!POP_SAMPLES[@]}"; do
    echo "" >&2
    echo "$pop:" >&2
    
    # Extract DoS values for this population and calculate stats
    grep "^${pop}" "$OUTPUT" | awk -F'\t' '
    {
        if ($8 != "NA") {
            print $8
        }
    }' | sort -n > "$TMPDIR/${pop}_dos_sorted.txt"
    
    if [[ -s "$TMPDIR/${pop}_dos_sorted.txt" ]]; then
        total=$(wc -l < "$TMPDIR/${pop}_dos_sorted.txt")
        mean=$(awk '{sum+=$1} END {print sum/NR}' "$TMPDIR/${pop}_dos_sorted.txt")
        median=$(awk '{a[NR]=$1} END {if(NR%2==1) print a[int(NR/2)+1]; else print (a[NR/2]+a[NR/2+1])/2}' "$TMPDIR/${pop}_dos_sorted.txt")
        pos_sel=$(awk '$1 > 0' "$TMPDIR/${pop}_dos_sorted.txt" | wc -l)
        neg_sel=$(awk '$1 < 0' "$TMPDIR/${pop}_dos_sorted.txt" | wc -l)
        upper_10=$(awk -v n=$total 'NR==int(n*0.9) {print}' "$TMPDIR/${pop}_dos_sorted.txt")
        lower_10=$(awk -v n=$total 'NR==int(n*0.1) {print}' "$TMPDIR/${pop}_dos_sorted.txt")
        
        echo "  Total genes analyzed: $total" >&2
        printf "  Mean DoS: %.4f\n" "$mean" >&2
        printf "  Median DoS: %.4f\n" "$median" >&2
        echo "  Genes under positive selection (DoS > 0): $pos_sel" >&2
        echo "  Genes under purifying selection (DoS < 0): $neg_sel" >&2
        printf "  Upper 10 percent DoS threshold: %.4f\n" "$upper_10" >&2
        printf "  Lower 10 percent DoS threshold: %.4f\n" "$lower_10" >&2
    else
        echo "  No data for this population" >&2
    fi
done

echo "" >&2
echo "Done!" >&2
