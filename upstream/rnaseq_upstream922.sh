#!/bin/bash
#$ -S /bin/bash
#$ -cwd
set -Eeuo pipefail
shopt -s nullglob

# Run from 01.RawData inside an allocated compute node.
# Confirm the library strandedness before final counting.
# Examples:
#   PREFLIGHT_ONLY=1 bash rnaseq_upstream.sh
#   THREADS=16 LIBRARY_TYPE=reverse MODE=both bash rnaseq_upstream.sh
# RAW/OUT/REF and other settings can be overridden at launch.

die() {
    echo "ERROR: $*" >&2
    exit 1
}
need_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "command not found: $1"
}
need_file() {
    [[ -f "$1" && -r "$1" && -s "$1" ]] || die "file missing, unreadable or empty: $1"
}
positive_int() {
    [[ "$2" =~ ^[1-9][0-9]*$ ]] || die "$1 must be a positive integer"
}
boolean_setting() {
    [[ "$2" == 0 || "$2" == 1 ]] || die "$1 must be 0 or 1"
}
trap 'rc=$?; printf "ERROR: line %s failed (exit %s)\n" "$LINENO" "$rc" >&2; exit "$rc"' ERR
(( BASH_VERSINFO[0] >= 4 )) || die "Bash 4 or newer is required; run with bash, not sh"

RAW="${RAW:-$PWD}"
[[ -d "$RAW" ]] || die "input directory not found: $RAW"
RAW="$(cd -- "$RAW" && pwd -P)"
OUT="${OUT:-$RAW/rnaseq_upstream}"
need_cmd realpath
OUT="$(realpath -m -- "$OUT")"
[[ "$OUT" != "$RAW" ]] || die "OUT must be a separate result directory, not RAW itself"

# Retained from the original script; these are paths on the Linux server.
REF="${REF:-/ref/hg38}"
FASTA="${FASTA:-${REF}/hg38.fa}"
GTF="${GTF:-${REF}/gencode.v48.annotation.gtf}"
HISAT_INDEX="${HISAT_INDEX:-${REF}/hg38.hisat_index}"
FASTA="$(realpath -m -- "$FASTA")"
GTF="$(realpath -m -- "$GTF")"
HISAT_INDEX="$(realpath -m -- "$HISAT_INDEX")"

THREADS="${THREADS:-16}"
positive_int THREADS "$THREADS"
(( THREADS >= 2 )) || die "THREADS must be >= 2 for the concurrent HISAT2/samtools pipe"
# Reserve one core for the samtools main thread, in addition to sort workers.
HISAT_THREADS="${HISAT_THREADS:-$(( THREADS / 2 ))}"
positive_int HISAT_THREADS "$HISAT_THREADS"
SORT_THREADS="${SORT_THREADS:-$(( THREADS - HISAT_THREADS - 1 ))}"
[[ "$SORT_THREADS" =~ ^(0|[1-9][0-9]*)$ ]] || die "SORT_THREADS must be a non-negative integer"
(( HISAT_THREADS + SORT_THREADS + 1 <= THREADS )) || die "HISAT_THREADS + SORT_THREADS + 1 exceeds THREADS"
FASTQC_THREADS="${FASTQC_THREADS:-$(( THREADS < 4 ? THREADS : 4 ))}"
positive_int FASTQC_THREADS "$FASTQC_THREADS"
(( FASTQC_THREADS <= THREADS )) || die "FASTQC_THREADS exceeds THREADS"
SAMTOOLS_THREADS=$(( THREADS - 1 ))
SORT_MEM="${SORT_MEM:-768M}"  # samtools sort memory limit per thread
LIBRARY_TYPE="${LIBRARY_TYPE:-unstranded}"  # unstranded | forward | reverse
MODE="${MODE:-gene}"                       # gene | circrna | both
SKIP_MD5="${SKIP_MD5:-0}"                  # 1 only if checksums were already checked
PREFLIGHT_ONLY="${PREFLIGHT_ONLY:-0}"      # inspect inputs/tools; no result files written
RESUME_FROM_ALIGNMENT="${RESUME_FROM_ALIGNMENT:-0}"
AUTO_ACTIVATE_CONDA="${AUTO_ACTIVATE_CONDA:-1}"
CONDA_SH="${CONDA_SH:-/anaconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-rnaseq_upstream}"
for setting in SKIP_MD5 PREFLIGHT_ONLY RESUME_FROM_ALIGNMENT AUTO_ACTIVATE_CONDA; do
    boolean_setting "$setting" "${!setting}"
done
case "$MODE" in
    gene|circrna|both) ;;
    *) die "MODE must be gene, circrna, or both" ;;
esac
case "$LIBRARY_TYPE" in
    unstranded) HISAT_STRAND=(); FC_STRAND=0; CIRI_STRAND=0 ;;
    forward) HISAT_STRAND=(--rna-strandness FR); FC_STRAND=1; CIRI_STRAND=1 ;;
    reverse) HISAT_STRAND=(--rna-strandness RF); FC_STRAND=2; CIRI_STRAND=2 ;;
    *) die "LIBRARY_TYPE must be unstranded, forward, or reverse" ;;
esac

if [[ "$AUTO_ACTIVATE_CONDA" == 1 ]]; then
    need_file "$CONDA_SH"
    # Some conda activation hooks reference unset variables.
    set +u
    source "$CONDA_SH"
    conda activate "$CONDA_ENV"
    set -u
fi
unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
hash -r
printf 'Conda environment: %s\nRAW: %s\nOUT: %s\nMODE: %s\nLIBRARY_TYPE: %s\n' \
    "${CONDA_PREFIX:-<externally managed>}" "$RAW" "$OUT" "$MODE" "$LIBRARY_TYPE"
echo "Library type is a user setting; the script does not infer it from FASTQ names."

need_file "$FASTA"
need_file "$GTF"
# HISAT2 accepts both small (.ht2) and large (.ht2l) index sets.
HISAT_INDEX_EXT=""
for ext in ht2 ht2l; do
    complete=1
    for i in {1..8}; do
        [[ -r "${HISAT_INDEX}.${i}.${ext}" && -s "${HISAT_INDEX}.${i}.${ext}" ]] || complete=0
    done
    if [[ "$complete" == 1 ]]; then
        HISAT_INDEX_EXT="$ext"
        break
    fi
done
[[ -n "$HISAT_INDEX_EXT" ]] || die "incomplete HISAT2 index: ${HISAT_INDEX}.[1-8].ht2 or .ht2l"

for cmd in hisat2 samtools multiqc sort awk; do
    need_cmd "$cmd"
done
if [[ "$RESUME_FROM_ALIGNMENT" == 0 ]]; then
    need_cmd fastqc
    need_cmd fastp
fi
if [[ "$SKIP_MD5" == 0 ]]; then
    need_cmd md5sum
fi
if [[ "$MODE" == gene || "$MODE" == both ]]; then
    need_cmd featureCounts
    # Help may return nonzero; avoid a pipefail/SIGPIPE false negative.
    FC_HELP="$(featureCounts -h 2>&1 || true)"
    [[ -n "$FC_HELP" ]] || die "featureCounts -h returned no help text"
    FC_PAIR_ARGS=(-p -B -C)
    if [[ "$FC_HELP" == *--countReadPairs* ]]; then
        FC_PAIR_ARGS+=(--countReadPairs)
    else
        # Only verified pre-2.0.2 versions may use -p alone for fragment counting.
        FC_VERSION="$(featureCounts -v 2>&1 || true)"
        if [[ "$FC_VERSION" =~ v([0-9]+)\.([0-9]+)\.([0-9]+) ]]; then
            major="${BASH_REMATCH[1]}"; minor="${BASH_REMATCH[2]}"; patch="${BASH_REMATCH[3]}"
            if ! (( major == 1 || (major == 2 && minor == 0 && patch < 2) )); then
                die "featureCounts version should support --countReadPairs but help did not advertise it"
            fi
        else
            die "cannot establish paired-fragment counting support in featureCounts"
        fi
    fi
    printf 'featureCounts pair options: %s\n' "${FC_PAIR_ARGS[*]}"
fi
if [[ "$MODE" == circrna || "$MODE" == both ]]; then
    for cmd in bwa stringtie CIRIquant java perl; do
        need_cmd "$cmd"
    done
    # Check that the installed CIRIquant entry point and its Python imports work.
    CIRIquant --version
fi

echo "[1/8] Discover sample directories and validate paired FASTQs"
# Confirmed layout: RAW/1N/1N_1.fq.gz, RAW/1N/1N_2.fq.gz, RAW/1N/MD5.txt.
# One biological sample per directory; no lane merge or raw-file symlinks needed.
SAMPLE_ORDER=()
for sample_dir in "$RAW"/*/; do
    sample="${sample_dir%/}"
    sample="${sample##*/}"
    [[ "$sample" =~ ^[0-9]+[NT]$ ]] || continue
    SAMPLE_ORDER+=("$sample")
done
(( ${#SAMPLE_ORDER[@]} > 0 )) || die "no sample directories such as 1N/ or 1T/ found in $RAW"
mapfile -t SAMPLE_ORDER < <(printf '%s\n' "${SAMPLE_ORDER[@]}" | sort -V)
R1_FILES=()
R2_FILES=()
declare -A SAMPLE_PRESENT
for sample in "${SAMPLE_ORDER[@]}"; do
    sample_dir="$RAW/$sample"
    sample_real="$(cd -- "$sample_dir" && pwd -P)"
    [[ "$OUT" != "$sample_real" && "$OUT" != "$sample_real/"* ]] || die "OUT must not be inside a sample directory: $sample"
    r1="$sample_dir/${sample}_1.fq.gz"
    r2="$sample_dir/${sample}_2.fq.gz"
    need_file "$r1"
    need_file "$r2"
    # Refuse extra FASTQs rather than silently omit lanes or another library.
    sample_fastqs=("$sample_dir"/*.fq.gz "$sample_dir"/*.fastq.gz "$sample_dir"/*.fq "$sample_dir"/*.fastq)
    (( ${#sample_fastqs[@]} == 2 )) || die "expected exactly ${sample}_1.fq.gz and ${sample}_2.fq.gz in $sample_dir; found ${#sample_fastqs[@]} FASTQ files"
    if [[ "$SKIP_MD5" == 0 ]]; then
        need_file "$sample_dir/MD5.txt"
    fi
    if [[ "$RESUME_FROM_ALIGNMENT" == 1 ]]; then
        need_file "$OUT/02_clean/${sample}_1.clean.fq.gz"
        need_file "$OUT/02_clean/${sample}_2.clean.fq.gz"
    fi
    SAMPLE_PRESENT[$sample]=1
    R1_FILES+=("$r1")
    R2_FILES+=("$r2")
done

emit_sample_sheet() {
    local sample patient condition pair_status
    printf 'sample\tpatient\tcondition\tpair_status\tfastq_1\tfastq_2\tmd5_file\n'
    for sample in "${SAMPLE_ORDER[@]}"; do
        patient="${sample%[NT]}"
        condition="${sample: -1}"
        pair_status=unpaired
        if [[ -n "${SAMPLE_PRESENT[${patient}N]:-}" && -n "${SAMPLE_PRESENT[${patient}T]:-}" ]]; then
            pair_status=paired
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$sample" "$patient" "$condition" "$pair_status" \
            "$RAW/$sample/${sample}_1.fq.gz" "$RAW/$sample/${sample}_2.fq.gz" "$RAW/$sample/MD5.txt"
    done
}
echo "Detected ${#SAMPLE_ORDER[@]} samples"
if [[ "$PREFLIGHT_ONLY" == 1 ]]; then
    emit_sample_sheet
    echo "Preflight passed: file presence and required commands checked; no analysis/output files created."
    echo "MD5 content, FASTQ records, reference consistency, output write access and alignment have NOT been tested."
    exit 0
fi

mkdir -p \
    "$OUT/00_metadata" \
    "$OUT/01_fastqc_raw" \
    "$OUT/02_clean" \
    "$OUT/03_bam" \
    "$OUT/04_counts" \
    "$OUT/05_multiqc" \
    "$OUT/06_circrna"
SAMPLE_SHEET="$OUT/00_metadata/samples.tsv"
emit_sample_sheet > "$SAMPLE_SHEET"
column -t -s $'\t' "$SAMPLE_SHEET" 2>/dev/null || cat "$SAMPLE_SHEET"
printf 'RAW=%s\nOUT=%s\nFASTA=%s\nGTF=%s\nHISAT_INDEX=%s\nMODE=%s\nLIBRARY_TYPE=%s\nTHREADS=%s\nRESUME_FROM_ALIGNMENT=%s\n' \
    "$RAW" "$OUT" "$FASTA" "$GTF" "$HISAT_INDEX" "$MODE" "$LIBRARY_TYPE" "$THREADS" "$RESUME_FROM_ALIGNMENT" \
    > "$OUT/00_metadata/run_settings.txt"

echo "[2/8] Check MD5 for each sample before QC/alignment"
if [[ "$SKIP_MD5" == 0 ]]; then
    : > "$OUT/00_metadata/md5check.log"
    for sample in "${SAMPLE_ORDER[@]}"; do
        (
            cd -- "$RAW/$sample"
            echo "Sample: $sample"
            # Standard md5sum format: HASH  filename (or HASH *filename).
            # Normalize CRLF in a pipe; never rewrite the vendor's MD5.txt.
            # Require both input files to be listed before claiming validation.
            for mate in 1 2; do
                expected="${sample}_${mate}.fq.gz"
                awk -v expected="$expected" '
                    { sub(/\r$/, ""); hash=$1; name=$0
                      sub(/^[^[:space:]]+[[:space:]]+\*?/, "", name)
                      sub(/^\.\//, "", name)
                      if (length(hash)==32 && hash !~ /[^[:xdigit:]]/ && name==expected) found=1 }
                    END { exit !found }
                ' MD5.txt || die "MD5.txt must contain a checksum for $expected (md5sum format, relative to this sample directory)"
            done
            sed 's/\r$//' MD5.txt | md5sum --strict -c -
        ) 2>&1 | tee -a "$OUT/00_metadata/md5check.log"
    done
else
    echo "MD5 verification explicitly skipped (SKIP_MD5=1)" | tee "$OUT/00_metadata/md5check.log"
fi

FASTP_THREADS="$THREADS"
((FASTP_THREADS > 16)) && FASTP_THREADS=16
CLEAN_R1_FILES=()
CLEAN_R2_FILES=()

if [[ "$RESUME_FROM_ALIGNMENT" == "1" ]]; then
    echo "[3-4/8] Resume mode: reuse existing clean FASTQs, skip FastQC/fastp"
else
    echo "[3/8] Raw FASTQ QC"
    fastqc \
        --threads "$FASTQC_THREADS" \
        --outdir "$OUT/01_fastqc_raw" \
        "${R1_FILES[@]}" "${R2_FILES[@]}"

    echo "[4/8] Adapter and quality trimming"
fi

for r1 in "${R1_FILES[@]}"; do
    sample="$(basename "$r1" _1.fq.gz)"
    r2="${r1%_1.fq.gz}_2.fq.gz"
    clean_r1="$OUT/02_clean/${sample}_1.clean.fq.gz"
    clean_r2="$OUT/02_clean/${sample}_2.clean.fq.gz"

    if [[ "$RESUME_FROM_ALIGNMENT" == "1" ]]; then
        need_file "$clean_r1"
        need_file "$clean_r2"
    else
        fastp \
            --in1 "$r1" \
            --in2 "$r2" \
            --out1 "$clean_r1" \
            --out2 "$clean_r2" \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --unqualified_percent_limit 40 \
            --length_required 35 \
            --thread "$FASTP_THREADS" \
            --html "$OUT/02_clean/${sample}.fastp.html" \
            --json "$OUT/02_clean/${sample}.fastp.json"
    fi

    CLEAN_R1_FILES+=("$clean_r1")
    CLEAN_R2_FILES+=("$clean_r2")
done

if [[ "$RESUME_FROM_ALIGNMENT" != "1" ]]; then
    fastqc \
        --threads "$FASTQC_THREADS" \
        --outdir "$OUT/02_clean" \
        "${CLEAN_R1_FILES[@]}" "${CLEAN_R2_FILES[@]}"
fi

if [[ "$MODE" == "gene" || "$MODE" == "both" ]]; then
    echo "[5/8] HISAT2 alignment and BAM QC"
    SPLICE_ARGS=()
    if command -v hisat2_extract_splice_sites.py >/dev/null 2>&1; then
        SPLICE_SITES="$OUT/00_metadata/annotation.splice_sites.txt"
        hisat2_extract_splice_sites.py "$GTF" > "$SPLICE_SITES"
        SPLICE_ARGS=(--known-splicesite-infile "$SPLICE_SITES")
    else
        echo "WARNING: hisat2_extract_splice_sites.py not found; continuing without an explicit splice-site file" >&2
    fi

    BAM_FILES=()
    for clean_r1 in "${CLEAN_R1_FILES[@]}"; do
        sample="$(basename "$clean_r1" _1.clean.fq.gz)"
        clean_r2="$OUT/02_clean/${sample}_2.clean.fq.gz"
        bam="$OUT/03_bam/${sample}.sorted.bam"

        hisat2 \
            --threads "$HISAT_THREADS" \
            -x "$HISAT_INDEX" \
            ${HISAT_STRAND[@]+"${HISAT_STRAND[@]}"} \
            ${SPLICE_ARGS[@]+"${SPLICE_ARGS[@]}"} \
            -1 "$clean_r1" \
            -2 "$clean_r2" \
            --summary-file "$OUT/03_bam/${sample}.hisat2.summary.txt" \
            2> "$OUT/03_bam/${sample}.hisat2.stderr.log" \
        | samtools sort \
            --threads "$SORT_THREADS" \
            -m "$SORT_MEM" \
            -o "$bam" -

        samtools quickcheck -v "$bam"
        samtools index -@ "$SAMTOOLS_THREADS" "$bam"
        samtools flagstat -@ "$SAMTOOLS_THREADS" "$bam" \
            > "$OUT/03_bam/${sample}.flagstat.txt"
        samtools stats -@ "$SAMTOOLS_THREADS" "$bam" \
            > "$OUT/03_bam/${sample}.samtools.stats.txt"
        BAM_FILES+=("$bam")
    done

    echo "[6/8] Gene-level raw counts with featureCounts"

    featureCounts \
        -T "$THREADS" \
        -a "$GTF" \
        -o "$OUT/04_counts/gene_counts.featureCounts.txt" \
        -t exon \
        -g gene_id \
        -s "$FC_STRAND" \
        "${FC_PAIR_ARGS[@]}" \
        "${BAM_FILES[@]}"

    awk '
        BEGIN { FS=OFS="\t" }
        $1 == "Geneid" {
            printf "gene_id"
            for (i=7; i<=NF; i++) {
                name=$i
                sub(/^.*\//, "", name)
                sub(/\.sorted\.bam$/, "", name)
                printf OFS name
            }
            print ""
            next
        }
        $1 !~ /^#/ {
            printf "%s", $1
            for (i=7; i<=NF; i++) printf OFS $i
            print ""
        }
    ' "$OUT/04_counts/gene_counts.featureCounts.txt" \
      > "$OUT/04_counts/gene_counts.matrix.tsv"
else
    echo "[5-6/8] Standard gene alignment/counting skipped because MODE=$MODE"
fi

if [[ "$MODE" == "circrna" || "$MODE" == "both" ]]; then
    echo "[7/8] circRNA back-splice junction detection and quantification"
    CIRI_REF_DIR="$OUT/06_circrna/reference"
    CIRI_RESULTS="$OUT/06_circrna/results"
    BWA_INDEX="${BWA_INDEX:-$CIRI_REF_DIR/hg38.bwa_index}"
    CIRI_CONFIG="$CIRI_REF_DIR/hg38.ciriquant.yml"
    mkdir -p "$CIRI_REF_DIR" "$CIRI_RESULTS"

    BWA_INDEX_COMPLETE=1
    for ext in amb ann bwt pac sa; do
        [[ -s "${BWA_INDEX}.${ext}" ]] || BWA_INDEX_COMPLETE=0
    done
    if [[ "$BWA_INDEX_COMPLETE" != "1" ]]; then
        bwa index -p "$BWA_INDEX" "$FASTA"
    fi

    if [[ ! -s "${FASTA}.fai" ]]; then
        [[ -w "$(dirname "$FASTA")" ]] || die "FASTA directory is not writable; prepare ${FASTA}.fai with samtools faidx first"
        samtools faidx "$FASTA"
    fi

    printf '%s\n' \
        'name: hg38' \
        'reference:' \
        "  fasta: $FASTA" \
        "  gtf: $GTF" \
        "  bwa_index: $BWA_INDEX" \
        "  hisat_index: $HISAT_INDEX" \
        > "$CIRI_CONFIG"

    for clean_r1 in "${CLEAN_R1_FILES[@]}"; do
        sample="$(basename "$clean_r1" _1.clean.fq.gz)"
        clean_r2="$OUT/02_clean/${sample}_2.clean.fq.gz"
        mkdir -p "$CIRI_RESULTS/$sample"

        CIRIquant \
            --threads "$THREADS" \
            --read1 "$clean_r1" \
            --read2 "$clean_r2" \
            --config "$CIRI_CONFIG" \
            --out "$CIRI_RESULTS/$sample" \
            --prefix "$sample" \
            --library-type "$CIRI_STRAND"
    done
else
    echo "[7/8] circRNA branch skipped; use MODE=both to enable it"
fi

echo "[8/8] Aggregate QC reports"
multiqc \
    --force \
    --outdir "$OUT/05_multiqc" \
    "$OUT/01_fastqc_raw" \
    "$OUT/02_clean" \
    "$OUT/03_bam" \
    "$OUT/04_counts" \
    "$OUT/06_circrna"

echo "Done"
echo "Sample metadata: $SAMPLE_SHEET"
if [[ "$MODE" == "gene" || "$MODE" == "both" ]]; then
    echo "Gene count matrix: $OUT/04_counts/gene_counts.matrix.tsv"
fi
if [[ "$MODE" == "circrna" || "$MODE" == "both" ]]; then
    echo "CIRIquant results: $OUT/06_circrna/results"
fi
echo "MultiQC report: $OUT/05_multiqc/multiqc_report.html"
