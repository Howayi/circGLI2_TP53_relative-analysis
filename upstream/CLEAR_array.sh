#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N CLEAR
#$ -pe smp 8
#$ -l h_vmem=10G
#$ -tc 6
# SGE does not expand variables in #$ lines, so these paths must be literal.
# The logs directory must exist before you run qsub.
#$ -o /X101SC26080554-Z01-J006/01.RawData/rnaseq_upstream/logs/
#$ -e/X101SC26080554-Z01-J006/01.RawData/rnaseq_upstream/logs/
#
# Submit as an array job, one task per line of samples.txt:
#   qsub -t 1-$(wc -l < samples.txt) CLEAR_array.sh

set -o pipefail

BASE_DIR="/X101SC26080554-Z01-J006/01.RawData/rnaseq_upstream"
INPUT_DIR="${BASE_DIR}/02_clean"
OUTPUT_DIR="${BASE_DIR}/CLEAR_output"
LOG_DIR="${BASE_DIR}/logs"
SAMPLE_LIST="${BASE_DIR}/samples.txt"
CONDA_SH="/anaconda3/etc/profile.d/conda.sh"
CONDA_ENV="clear_py27"
GENOME_FA="/ref/hg38/hg38.fa"
HISAT_INDEX="/ref/hg38/hg38.hisat_index"
BOWTIE_INDEX="/ref/hg38/hg38.bowtie_index"
GTF="/ref/hg38/gencode.v48.annotation.gtf"

# Match the tool's thread count to the slots SGE actually gave the job.
THREADS="${NSLOTS:-8}"

# SGE sets TMPDIR to a node-local per-task directory and removes it afterwards.
# Falls back to the shared directory if it is not set.
SCRATCH="${TMPDIR:-${BASE_DIR}/CLEAR_tmp}"

if [[ -z "${SGE_TASK_ID:-}" || "${SGE_TASK_ID}" == "undefined" ]]; then
    echo "ERROR: submit as an array job: qsub -t 1-\$(wc -l < samples.txt) CLEAR_array.sh" >&2
    exit 1
fi

cd "${BASE_DIR}" || exit 1
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${SCRATCH}"

sample="$(sed -n "${SGE_TASK_ID}p" "${SAMPLE_LIST}")"
if [[ -z "${sample}" ]]; then
    echo "ERROR: no sample on line ${SGE_TASK_ID} of ${SAMPLE_LIST}" >&2
    exit 1
fi

r1_gz="${INPUT_DIR}/${sample}_1.clean.fq.gz"
r2_gz="${INPUT_DIR}/${sample}_2.clean.fq.gz"
sample_out="${OUTPUT_DIR}/${sample}"
done_flag="${sample_out}/.clear_done"

# Lets you resubmit the whole array after failures; finished samples are skipped.
# quant/quant.txt also counts as finished, for samples done by the old script.
if [[ -f "${done_flag}" || -s "${sample_out}/quant/quant.txt" ]]; then
    echo "Already finished, skipping: ${sample}"
    exit 0
fi

# clear_quant calls os.mkdir() on its output directory and aborts if it exists,
# so move any leftover from an unfinished run aside (never deleted).
if [[ -e "${sample_out}" ]]; then
    stale="${sample_out}.incomplete_$(date +%Y%m%d_%H%M%S)"
    echo "Moving unfinished output aside: ${sample_out} -> ${stale}"
    mv "${sample_out}" "${stale}" || exit 1
fi

for f in "${r1_gz}" "${r2_gz}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: input file missing for ${sample}: ${f}" >&2
        exit 1
    fi
done

if [[ ! -f "${CONDA_SH}" ]]; then
    echo "ERROR: conda initialization script not found: ${CONDA_SH}" >&2
    exit 1
fi
source "${CONDA_SH}" || exit 1
conda activate "${CONDA_ENV}" || exit 1
unset PERL5LIB
unset PERLLIB
unset PERL5OPT

if ! command -v clear_quant >/dev/null 2>&1; then
    echo "ERROR: clear_quant was not found in conda environment ${CONDA_ENV}." >&2
    exit 1
fi

echo "======================================"
echo "Processing ${sample} (task ${SGE_TASK_ID}) on $(hostname), ${THREADS} threads"
echo "Scratch: ${SCRATCH}"
date
echo "======================================"

r1_fastq="${SCRATCH}/${sample}_1.clean.fastq"
r2_fastq="${SCRATCH}/${sample}_2.clean.fastq"
# Remove the uncompressed FASTQ however the task ends.
trap 'rm -f "${r1_fastq}" "${r2_fastq}"' EXIT

# Decompress both mates at the same time.
gzip -dc "${r1_gz}" > "${r1_fastq}" & pid1=$!
gzip -dc "${r2_gz}" > "${r2_fastq}" & pid2=$!
wait "${pid1}"; s1=$?
wait "${pid2}"; s2=$?
if (( s1 != 0 || s2 != 0 )); then
    echo "ERROR: failed to decompress input for ${sample} (check scratch space: df -h ${SCRATCH})" >&2
    exit 1
fi
ls -lh "${r1_fastq}" "${r2_fastq}"

clear_quant \
    -1 "${r1_fastq}" \
    -2 "${r2_fastq}" \
    -g "${GENOME_FA}" \
    -i "${HISAT_INDEX}" \
    -j "${BOWTIE_INDEX}" \
    -G "${GTF}" \
    -p "${THREADS}" \
    -o "${sample_out}" \
    > "${LOG_DIR}/${sample}_clear_quant.log" \
    2> "${LOG_DIR}/${sample}_clear_quant.err"
clear_status=$?

if (( clear_status == 0 )); then
    touch "${done_flag}"
    echo "Finished ${sample}"
else
    echo "ERROR: clear_quant failed for ${sample} (exit ${clear_status})." >&2
    echo "See ${LOG_DIR}/${sample}_clear_quant.err" >&2
fi
date
exit "${clear_status}"
