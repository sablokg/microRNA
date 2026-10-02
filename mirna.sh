#!/usr/bin/env bash
#
# ============================================================================
#  mirna_detection_pipeline.sh Gaurav Sablok gsablok@proton.me
#  Complete small RNA-seq -> miRNA detection pipeline
#
#  Stages:
#    1. QC                 -> FastQC / MultiQC
#    2. Adapter trimming    -> fastp (fast, modern replacement for cutadapt)
#    3. Length filtering     -> keep 18-25 nt reads (mature miRNA size range)
#    4. Known miRNA quant + isomiRs -> miRge3.0 (2021, actively maintained,
#         Python3, handles UMIs, much faster than miRDeep2, reports isomiRs
#         and tRFs in GFF3)
#    5. Novel miRNA discovery       -> miRDeep2 (still the gold-standard for
#         de novo/novel miRNA prediction via hairpin structure scoring;
#         miRge3.0 does not do full novel prediction)
#    6. Merge + summary report
#
#  Usage:
#    ./mirna_detection_pipeline.sh -i sample.fastq.gz -g hg38 -o results/
#
#  Requires: conda/mamba
# ============================================================================

set -euo pipefail

# ---------------------------- CONFIG ---------------------------------------
SAMPLE=""
GENOME_BUILD="hg38"          # hg38 | mm10 | etc. Adjust species DBs below.
OUTDIR="mirna_results"
THREADS=8
ADAPTER="AGATCGGAAGAGCACACGTCT"   # default Illumina small RNA adapter; edit if needed
MIN_LEN=18
MAX_LEN=25

usage() {
  echo "Usage: $0 -i <reads.fastq.gz> [-g hg38] [-o outdir] [-t threads] [-a adapter]"
  exit 1
}

while getopts "i:g:o:t:a:h" opt; do
  case $opt in
    i) SAMPLE="$OPTARG" ;;
    g) GENOME_BUILD="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    a) ADAPTER="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -z "$SAMPLE" ]] && usage
[[ ! -f "$SAMPLE" ]] && { echo "Input file not found: $SAMPLE"; exit 1; }

BASENAME=$(basename "$SAMPLE" | sed -E 's/\.(fastq|fq)(\.gz)?$//')
mkdir -p "$OUTDIR"/{qc_raw,trimmed,qc_trimmed,mirge3,mirdeep2,logs,final}

echo "=== miRNA detection pipeline: $BASENAME ==="

# ---------------------------- 0. ENV SETUP ----------------------------------
# One-time environment build. Comment out after first successful run.
setup_env() {
  echo "[setup] Creating conda environment 'mirna-env' ..."
  mamba create -y -n mirna-env -c bioconda -c conda-forge \
      fastqc multiqc fastp bowtie=1.3.1 samtools=1.19 \
      mirdeep2 mirge3 sra-tools
}
# Uncomment the line below on first run:
# setup_env

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mirna-env

# ---------------------------- 1. RAW QC --------------------------------------
echo "[1/6] Raw read QC (FastQC)"
fastqc -t "$THREADS" -o "$OUTDIR/qc_raw" "$SAMPLE" > "$OUTDIR/logs/fastqc_raw.log" 2>&1

# ---------------------------- 2. ADAPTER TRIMMING ----------------------------
echo "[2/6] Adapter trimming + quality filtering (fastp)"
TRIMMED="$OUTDIR/trimmed/${BASENAME}.trimmed.fastq.gz"
fastp \
  -i "$SAMPLE" \
  -o "$TRIMMED" \
  --adapter_sequence "$ADAPTER" \
  --length_required "$MIN_LEN" \
  --length_limit "$MAX_LEN" \
  --qualified_quality_phred 20 \
  --thread "$THREADS" \
  --json "$OUTDIR/logs/${BASENAME}.fastp.json" \
  --html "$OUTDIR/logs/${BASENAME}.fastp.html" \
  > "$OUTDIR/logs/fastp.log" 2>&1

echo "[2/6] Post-trim QC (FastQC)"
fastqc -t "$THREADS" -o "$OUTDIR/qc_trimmed" "$TRIMMED" >> "$OUTDIR/logs/fastqc_trimmed.log" 2>&1

# ---------------------------- 3. REFERENCE SETUP -----------------------------
# Adjust paths to your local reference / miRBase / miRge3.0 library installs.
GENOME_FA="refs/${GENOME_BUILD}/genome.fa"
BOWTIE_INDEX="refs/${GENOME_BUILD}/bowtie_index/${GENOME_BUILD}"
MIRBASE_HAIRPIN="refs/mirbase/hairpin.fa"
MIRBASE_MATURE="refs/mirbase/mature.fa"
MIRGE3_LIB="refs/mirge3_lib"          # downloaded via `miRge3.0 annotate --download`
SPECIES="human"                        # human | mouse | ... (miRge3.0 naming)

# ---------------------------- 4. KNOWN miRNA QUANT (miRge3.0) ---------------
echo "[3/6] Known miRNA quantification + isomiRs (miRge3.0)"
miRge3.0 \
  -s "$TRIMMED" \
  -lib "$MIRGE3_LIB" \
  -on "$SPECIES" \
  -db mirbase \
  -o "$OUTDIR/mirge3" \
  -gff \
  -a illumina \
  -cpu "$THREADS" \
  > "$OUTDIR/logs/mirge3.log" 2>&1

# ---------------------------- 5. NOVEL miRNA DISCOVERY (miRDeep2) -----------
echo "[4/6] Genome alignment for novel miRNA discovery (bowtie -> miRDeep2)"
COLLAPSED="$OUTDIR/mirdeep2/${BASENAME}_collapsed.fa"
mapper.pl "$TRIMMED" -e -h -i -j -l "$MIN_LEN" -m \
  -p "$BOWTIE_INDEX" \
  -s "$COLLAPSED" \
  -t "$OUTDIR/mirdeep2/${BASENAME}_reads_vs_genome.arf" \
  -v > "$OUTDIR/logs/mapper.log" 2>&1

echo "[5/6] Novel + known miRNA prediction (miRDeep2 core)"
cd "$OUTDIR/mirdeep2"
miRDeep2.pl \
  "../../${COLLAPSED}" \
  "../../${GENOME_FA}" \
  "../../${OUTDIR}/mirdeep2/${BASENAME}_reads_vs_genome.arf" \
  "../../${MIRBASE_MATURE}" \
  none \
  "../../${MIRBASE_HAIRPIN}" \
  -t human \
  > "../../${OUTDIR}/logs/mirdeep2.log" 2>&1
cd - > /dev/null

# ---------------------------- 6. SUMMARY -------------------------------------
echo "[6/6] Aggregating reports (MultiQC + summary)"
multiqc "$OUTDIR" -o "$OUTDIR/final" > "$OUTDIR/logs/multiqc.log" 2>&1

cat <<EOF > "$OUTDIR/final/README.txt"
miRNA detection pipeline complete for sample: $BASENAME

Outputs:
  - Known miRNA counts / isomiRs / GFF3 : $OUTDIR/mirge3/
  - Novel miRNA candidates (survey.csv, result*.html) : $OUTDIR/mirdeep2/
  - QC reports : $OUTDIR/final/multiqc_report.html
  - Raw logs : $OUTDIR/logs/

Next steps:
  - Filter miRDeep2 novel candidates by miRDeep score >= 4 and randfold p<0.05
  - Merge known (miRge3.0) + novel (miRDeep2) counts for a combined expression matrix
  - For differential expression, feed the count matrix to DESeq2 / edgeR in R
EOF

echo "=== Done. See $OUTDIR/final/README.txt ==="
