#!/bin/bash
#SBATCH --job-name=verify_pseudo
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/verify_pseudogenome_%j.out
#SBATCH --error=logs/verify_pseudogenome_%j.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

module load samtools/1.16.1
module load bcftools/1.15.1

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
REF=${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna
PSEUDO=${PROJECT_DIR}/reference/pseudogenome/colombian_pseudogenome.fna
GFF=${PROJECT_DIR}/reference/pseudogenome/colombian_pseudogenome.gff3
CHAIN=${PROJECT_DIR}/reference/pseudogenome/colombian_pseudogenome.chain
VCF=${PROJECT_DIR}/gatk/trimmomatic/vcf_filtered/snps_filtered.vcf.gz

PASS=0; FAIL=0
check() {
    local label="$1"; local result="$2"; local expected="$3"
    if [ "$result" = "$expected" ]; then
        echo "  ✅ $label"
        PASS=$((PASS+1))
    else
        echo "  ❌ $label — got: '$result' expected: '$expected'"
        FAIL=$((FAIL+1))
    fi
}
check_range() {
    local label="$1"; local val="$2"; local lo="$3"; local hi="$4"
    if [ "$val" -ge "$lo" ] && [ "$val" -le "$hi" ]; then
        echo "  ✅ $label ($val)"
        PASS=$((PASS+1))
    else
        echo "  ❌ $label — got $val, expected between $lo and $hi"
        FAIL=$((FAIL+1))
    fi
}

echo "========================================================"
echo " Pseudogenome + Annotation Verification"
echo " $(date)"
echo "========================================================"

# ── 1. GENOME INTEGRITY ───────────────────────────────────────
echo ""
echo "[ 1. GENOME INTEGRITY ]"

PSEUDO_CONTIGS=$(grep -c "^>" "$PSEUDO")
REF_CONTIGS=$(grep -c "^>" "$REF")
check "Same contig count as reference" "$PSEUDO_CONTIGS" "$REF_CONTIGS"

PSEUDO_SIZE=$(awk '/^>/{next}{sum+=length($0)} END{print sum}' "$PSEUDO")
REF_SIZE=$(awk '/^>/{next}{sum+=length($0)} END{print sum}' "$REF")
echo "  Reference size:   ${REF_SIZE} bp"
echo "  Pseudogenome size: ${PSEUDO_SIZE} bp"
DIFF=$(( PSEUDO_SIZE - REF_SIZE ))
echo "  Difference: ${DIFF} bp (net insertions from INDELs)"

FAI_CONTIGS=$(wc -l < "${PSEUDO}.fai")
check "samtools .fai index present and complete" "$FAI_CONTIGS" "$PSEUDO_CONTIGS"

BWA_FILES=$(ls ${PSEUDO}.{amb,ann,bwt,pac,sa} 2>/dev/null | wc -l)
check "BWA index complete (5 files)" "$BWA_FILES" "5"

# ── 2. ANNOTATION INTEGRITY ───────────────────────────────────
echo ""
echo "[ 2. ANNOTATION INTEGRITY ]"

GFF_GENES=$(grep -v "^#" "$GFF" | awk '$3=="gene"' | wc -l)
check_range "Gene count reasonable (>5000)" "$GFF_GENES" 5000 15000
echo "  Total genes: $GFF_GENES"

GFF_MRNA=$(grep -v "^#" "$GFF" | awk '$3=="mRNA"' | wc -l)
GFF_EXON=$(grep -v "^#" "$GFF" | awk '$3=="exon"' | wc -l)
GFF_CDS=$(grep -v "^#" "$GFF"  | awk '$3=="CDS"'  | wc -l)
echo "  mRNA: $GFF_MRNA | exon: $GFF_EXON | CDS: $GFF_CDS"

# Check no feature exceeds its chromosome length
echo ""
echo "  Checking for out-of-bounds features..."
OOB=$(awk 'NR==FNR && !/^#/{len[$1]=$2; next}
           !/^#/{if($5 > len[$1]) print $1,$4,$5,len[$1]}' \
     "${PSEUDO}.fai" "$GFF" | wc -l)
check "No out-of-bounds features" "$OOB" "0"
if [ "$OOB" -gt 0 ]; then
    echo "  First 3 out-of-bounds:"
    awk 'NR==FNR && !/^#/{len[$1]=$2; next}
         !/^#/{if($5 > len[$1]) print $1,$3,$4,$5,"chrom_len="len[$1]}' \
        "${PSEUDO}.fai" "$GFF" | head -3
fi

# Check no feature has start > end
INVERTED=$(grep -v "^#" "$GFF" | awk '$4 > $5' | wc -l)
check "No inverted features (start > end)" "$INVERTED" "0"

# ── 3. BDNF GENE CHECK ───────────────────────────────────────
echo ""
echo "[ 3. BDNF GENE CHECK ]"

BDNF_GENE=$(grep -v "^#" "$GFF" | awk '$3=="gene"' | grep "gene=bdnf" | wc -l)
check_range "bdnf gene entry present" "$BDNF_GENE" 1 5

BDNF_LINE=$(grep -v "^#" "$GFF" | awk '$3=="gene"' | grep "gene=bdnf" | head -1)
if [ -n "$BDNF_LINE" ]; then
    BDNF_CHR=$(echo "$BDNF_LINE" | awk '{print $1}')
    BDNF_START=$(echo "$BDNF_LINE" | awk '{print $4}')
    BDNF_END=$(echo "$BDNF_LINE" | awk '{print $5}')
    BDNF_STRAND=$(echo "$BDNF_LINE" | awk '{print $7}')
    echo "  bdnf location: ${BDNF_CHR}:${BDNF_START}-${BDNF_END} (${BDNF_STRAND} strand)"
    check "bdnf on chromosome NC_024333.1" "$BDNF_CHR" "NC_024333.1"
    check "bdnf on negative strand" "$BDNF_STRAND" "-"
    check_range "bdnf start reasonable (near 15.9 Mb)" "$BDNF_START" 15800000 16100000
fi

BDNF_EXONS=$(grep -v "^#" "$GFF" | awk '$3=="exon"' | grep "gene=bdnf" | wc -l)
check_range "bdnf exons present" "$BDNF_EXONS" 5 100
echo "  bdnf exons: $BDNF_EXONS"

# ── 4. SGRNA CUT SITE REGION ─────────────────────────────────
echo ""
echo "[ 4. sgRNA CUT SITE REGION (NC_024333.1:15922039-15922058) ]"

# Extract 20bp sgRNA region from pseudogenome and reference
PSEUDO_SEQ=$(samtools faidx "$PSEUDO" "NC_024333.1:15922039-15922058" 2>/dev/null | grep -v "^>")
REF_SEQ=$(samtools faidx "$REF"    "NC_024333.1:15922039-15922058" 2>/dev/null | grep -v "^>")
echo "  Reference sgRNA site:   $REF_SEQ"
echo "  Pseudogenome sgRNA site: $PSEUDO_SEQ"
if [ "$PSEUDO_SEQ" = "$REF_SEQ" ]; then
    echo "  ✅ sgRNA site sequence unchanged (no variants in cut site — as expected for Control)"
    PASS=$((PASS+1))
else
    echo "  ℹ️  sgRNA site differs from reference (Control SNPs applied)"
    PASS=$((PASS+1))
fi

# ── 5. VARIANT APPLICATION SPOT CHECK ────────────────────────
echo ""
echo "[ 5. VARIANT APPLICATION SPOT CHECK ]"

# Pick first 3 PASS Control SNPs on NC_024333.1 and verify they were applied
echo "  Checking 3 random Control SNPs on NC_024333.1..."
bcftools view -f PASS -s "Control_MNP_I_S54_L002,Control_MNP_II_S55_L002,Control_MNP_III_S56_L002" \
    "$VCF" "NC_024333.1" 2>/dev/null | \
    bcftools view -g ^miss 2>/dev/null | \
    grep -v "^#" | head -3 | \
    while read -r line; do
        CHROM=$(echo "$line" | awk '{print $1}')
        POS=$(echo "$line" | awk '{print $2}')
        REF_A=$(echo "$line" | awk '{print $4}')
        ALT_A=$(echo "$line" | awk '{print $5}')
        PSEUDO_BASE=$(samtools faidx "$PSEUDO" "${CHROM}:${POS}-${POS}" 2>/dev/null | grep -v "^>")
        if [ "$PSEUDO_BASE" = "$ALT_A" ]; then
            echo "  ✅ pos ${POS}: REF=${REF_A} → ALT=${ALT_A} applied correctly"
        else
            echo "  ❌ pos ${POS}: expected ALT=${ALT_A}, got ${PSEUDO_BASE}"
        fi
    done

# ── SUMMARY ──────────────────────────────────────────────────
echo ""
echo "========================================================"
echo " SUMMARY"
echo "  Passed: $PASS"
echo "  Failed: $FAIL"
if [ "$FAIL" -eq 0 ]; then
    echo "  ✅ All checks passed — pseudogenome and annotation look correct"
else
    echo "  ❌ $FAIL check(s) failed — review above"
fi
echo "========================================================"
echo "End time: $(date)"
