#!/bin/bash
#SBATCH --job-name=tgsgapcloser_genome
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --output=logs/tgsgapcloser_genome_%j.out
#SBATCH --error=logs/tgsgapcloser_genome_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Roadmap option 2: gap-filling with the existing Nanopore data (see
# "Genome Assembly — Further Improvement Options" in CLAUDE.md Pending
# Analyses). Unlike NextPolish short-read polishing (tried 2026-08,
# no improvement - see "NextPolish — No Improvement" Known Issue), this
# uses genuinely NEW information (long reads spanning repeats/gaps that
# short reads can't) rather than reprocessing the same Illumina data.
#
# Input: pooled telencephalon Nanopore data (10 fish, same Colombian
# population, NOT the same 3 WGS individuals) - ~2.4Gbp combined, N50
# ~4080bp, only ~3.2x depth against the ~700-750Mb genome. Too shallow for
# a hybrid reassembly (ruled out 2026-07-10) but a reasonable fit for
# targeted gap-filling, since it only needs reads that span a gap, not
# genome-wide coverage. See [[nanopore_epigenome_data]] memory for the full
# provenance/yield breakdown.
#
# TGS-GapCloser v1.2.1 (tgsgapcloser_env, see 00_setup_tgsgapcloser_env.sh).
# --racon (not --ne or the --ngs/--pilon path) for error correction of
# newly-filled gap sequence - Nanopore reads have a much higher raw error
# rate than Illumina, so correction matters here in a way it didn't for
# NextPolish (which was correcting already-consistent short-read-derived
# sequence with no new information - see that Known Issue for why it
# didn't help). racon is purpose-built for long-read consensus polishing.
# --tgstype ont since this is Nanopore, not PacBio.

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
GENOME=${PROJECT_DIR}/reference/colombian_scaffolded_genome/colombian_scaffolded.fna
OUT_DIR=${PROJECT_DIR}/assembly/tgsgapcloser_output
OUT_PREFIX=${OUT_DIR}/colombian_gapfilled
RACON=/hpcfs/apps/conda4.12.0/envs/racon-1.5.0/bin/racon

NANOPORE_RUN1=/hpcfs/home/ing_civil/da.martinez33/Guppy_epigenome/Control_Epigen_Guppy_Telencefalo_051224/no_sample_id/20241205_1452_MN47264_FBA87937_26c2ff3d/fastq_pass/barcode01
NANOPORE_RUN2=/hpcfs/home/ing_civil/da.martinez33/Guppy_epigenome/Control_Epigen_Guppy_Telencefalo_051224/no_sample_id/20241206_1024_MN47264_FBA87937_7302966d/fastq_pass/barcode01
COMBINED_READS=${OUT_DIR}/nanopore_combined.fastq.gz

mkdir -p "$OUT_DIR" logs/

if [ ! -f "$GENOME" ]; then
    echo "ERROR: genome not found at $GENOME"
    exit 1
fi

# TGS-GapCloser writes done_stepN_tag marker files into the CURRENT WORKING
# DIRECTORY (here, wherever this script is launched from), NOT under
# --output - completely undocumented and independent of the --output
# prefix. Job 707123 silently "resumed" from a previous failed run's stale
# tags (steps 1/2.1 skipped as "already done") and then fatally errored
# looking for an output file that didn't exist, because --output's actual
# directory HAD been cleaned but these tags, sitting elsewhere, hadn't.
# Remove them at the start of every run so this can never resume from
# stale state.
rm -f done_step*_tag

echo "Start time: $(date)"

echo "Combining Nanopore reads from both runs into ${COMBINED_READS}..."
# NOT a plain `cat *.gz > combined.gz` - that produces a valid but
# multi-member gzip stream (one member per source file, ~490 members
# here). Rebuilding as a single clean gzip stream via zcat|gzip ruled out
# multi-stream gzip as a cause of the crash below (it wasn't - see next
# comment) but is still correct practice, so kept.
zcat "${NANOPORE_RUN1}"/*.fastq.gz "${NANOPORE_RUN2}"/*.fastq.gz | gzip > "$COMBINED_READS"
echo "Combined reads: $(du -h "$COMBINED_READS" | cut -f1)"

LINE_COUNT=$(zcat "$COMBINED_READS" | wc -l)
READ_COUNT=$((LINE_COUNT / 4))
echo "Combined read count: ${READ_COUNT} (${LINE_COUNT} lines)"
if [ $((LINE_COUNT % 4)) -ne 0 ]; then
    echo "ERROR: combined fastq line count not divisible by 4 - malformed file, aborting"
    exit 1
fi
if [ "$READ_COUNT" -lt 1000000 ]; then
    echo "ERROR: expected ~1.1M reads combined, got ${READ_COUNT} - something went wrong"
    exit 1
fi

# Convert FASTQ -> FASTA. Root cause of jobs 705937/710309 both crashing
# identically (Assertion 'line.size() > 1' failed in
# BGIQD::FASTA::Id_Desc_Head::Init, LoadONTReads) even after the
# multi-stream gzip fix above: the tgsgapcloser WRAPPER SCRIPT
# unconditionally calls its internal tgsgapcandidate/gapcloser binaries
# with --ont_reads_a (grep confirms --ont_reads_q, the FASTQ-format flag,
# is never used anywhere in the wrapper) even though the top-level --reads
# flag's own --help just says "input TGS read file" with no format
# specified. The internal binary's own --help (run it directly, e.g.
# .../tgsgapcloserbin/tgsgapcandidate with no args) reveals the two
# separate flags: --ont_reads_q "in fastq format" vs --ont_reads_a "in
# fasta format" - the wrapper only ever uses the latter. Feeding it our
# raw FASTQ under the FASTA-only code path made its header parser choke on
# '@'/'+' lines it wasn't expecting. Fix: pre-convert to FASTA ourselves.
COMBINED_FASTA=${OUT_DIR}/nanopore_combined.fasta.gz
echo "Converting to FASTA (tgsgapcloser wrapper only ever uses --ont_reads_a, never --ont_reads_q - see comment above)..."
zcat "$COMBINED_READS" | awk 'NR%4==1 {print ">"substr($0,2)} NR%4==2 {print}' | gzip > "$COMBINED_FASTA"
echo "FASTA reads: $(du -h "$COMBINED_FASTA" | cut -f1)"
FASTA_SEQ_COUNT=$(zcat "$COMBINED_FASTA" | grep -c "^>")
echo "FASTA sequence count: ${FASTA_SEQ_COUNT}"
if [ "$FASTA_SEQ_COUNT" -ne "$READ_COUNT" ]; then
    echo "ERROR: FASTA seq count (${FASTA_SEQ_COUNT}) != FASTQ read count (${READ_COUNT}) - conversion failed"
    exit 1
fi

CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate tgsgapcloser_env

tgsgapcloser \
    --scaff "$GENOME" \
    --reads "$COMBINED_FASTA" \
    --output "$OUT_PREFIX" \
    --racon "$RACON" \
    --tgstype ont \
    --thread "${SLURM_CPUS_PER_TASK}" \
    >"${OUT_PREFIX}.pipe.log" 2>&1

echo "End time: $(date)"

echo "=== Output files ==="
ls -la "${OUT_PREFIX}"* 2>&1
