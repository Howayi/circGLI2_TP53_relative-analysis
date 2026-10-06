#!/bin/bash
#$ -S /bin/bash
#$ -cwd
set -Eeuo pipefail
shopt -s nullglob

# Exploratory screen for candidate SNVs/small indels in CDKN2A/CDKN2B (INK4)
# and CTNNB1 (beta-catenin) in the HISAT2 BAMs from rnaseq_upstream922.sh,
# followed by a paired tumor-normal comparison.
#
# Every output row is a CANDIDATE for manual review (IGV, strand/position
# checks, RNA-editing databases, DNA confirmation). Labels are heuristics;
# none of them establishes a somatic DNA mutation.
#
# Run from rnaseq_upstream/03_bam on a compute node:
#   bash ink4_ctnnb1_check.sh
#   GENES="CDKN2A CDKN2B CDKN2C CDKN2D CTNNB1" MIN_ALT=5 bash ink4_ctnnb1_check.sh
# The file must have Unix (LF) line endings: sed -i 's/\r$//' ink4_ctnnb1_check.sh

die() { echo "ERROR: $*" >&2; exit 1; }
trap 'rc=$?; printf "ERROR: line %s failed (exit %s)\n" "$LINENO" "$rc" >&2; exit "$rc"' ERR

BAM_DIR="${BAM_DIR:-$PWD}"
[[ -d "$BAM_DIR" ]] || die "BAM_DIR not found: $BAM_DIR"
BAM_DIR="$(cd -- "$BAM_DIR" && pwd -P)"
OUT="${OUT:-$BAM_DIR/../07_ink4_ctnnb1}"
REF="${REF:-/ref/hg38}"
FASTA="${FASTA:-$REF/hg38.fa}"
GTF="${GTF:-$REF/gencode.v48.annotation.gtf}"
GENES="${GENES:-CDKN2A CDKN2B CTNNB1}"
PAD="${PAD:-10}"                 # bp added around each exon (covers splice sites)
MIN_ALT="${MIN_ALT:-3}"          # alt reads required in at least one sample
MIN_MAPQ="${MIN_MAPQ:-20}"       # HISAT2: 60 = unique, 0/1 = multi-mapped
MIN_BASEQ="${MIN_BASEQ:-20}"
MAX_DEPTH="${MAX_DEPTH:-100000}" # used for both -d and -L (see step 3)
MAX_P0="${MAX_P0:-0.05}"         # see tumor_normal_pairs labels in step 5
AUTO_ACTIVATE_CONDA="${AUTO_ACTIVATE_CONDA:-1}"
CONDA_SH="${CONDA_SH:-/anaconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-rnaseq_upstream}"

if [[ "$AUTO_ACTIVATE_CONDA" == 1 ]]; then
    set +u; source "$CONDA_SH"; conda activate "$CONDA_ENV"; set -u
fi
for cmd in bcftools awk sort gzip; do
    command -v "$cmd" >/dev/null 2>&1 || die "command not found: $cmd"
done
[[ -s "$FASTA" && -s "$GTF" ]] || die "FASTA or GTF missing"
[[ -s "$FASTA.fai" ]] || die "missing $FASTA.fai; run: samtools faidx $FASTA"
mpileup_help="$(bcftools mpileup 2>&1 || true)"
[[ "$mpileup_help" == *--ignore-RG* ]] || die "this bcftools has no mpileup --ignore-RG; please update bcftools"
bcftools_version="$(bcftools --version)"
bcftools_version="${bcftools_version%%$'\n'*}"

mkdir -p "$OUT"
cat > "$OUT/run_settings.txt" <<EOF
$bcftools_version
BAM_DIR=$BAM_DIR
FASTA=$FASTA
GTF=$GTF
GENES=$GENES
PAD=$PAD MIN_ALT=$MIN_ALT MIN_MAPQ=$MIN_MAPQ MIN_BASEQ=$MIN_BASEQ MAX_DEPTH=$MAX_DEPTH MAX_P0=$MAX_P0
mpileup: --ignore-RG -d $MAX_DEPTH -L $MAX_DEPTH -p; default read filters (UNMAP,SECONDARY,QCFAIL,DUP; anomalous pairs dropped)
Indel candidates need default -m/-F thresholds (bcftools 1.19: >= 2 gapped reads and >= 5% of reads) within one sample.
Lower-fraction indels are not detected; their absence is not evidence of no indel.
EOF

echo "[1/5] Samples and tumor/normal pairing"
BAMS=("$BAM_DIR"/*.sorted.bam)
(( ${#BAMS[@]} > 0 )) || die "no *.sorted.bam in $BAM_DIR"
mapfile -t BAMS < <(printf '%s\n' "${BAMS[@]}" | sort -V)
: > "$OUT/sample_names.txt"
for bam in "${BAMS[@]}"; do
    [[ "$bam" != *[[:space:]]* ]] || die "BAM path contains whitespace: $bam"
    [[ -s "$bam.bai" || -s "${bam%.bam}.bai" ]] || die "missing index for $bam"
    basename "$bam" .sorted.bam >> "$OUT/sample_names.txt"
done
# Naming convention: <patient number><N|T>, e.g. 1N / 1T. Everything else is reported.
awk 'BEGIN { FS=OFS="\t" }
    { s[NR] = $1
      if ($1 ~ /^[0-9]+[NT]$/) has[substr($1, 1, length($1)-1), substr($1, length($1))] = 1 }
    END { print "sample", "patient", "condition", "pair_status"
          for (i=1; i<=NR; i++) {
              x = s[i]
              if (x !~ /^[0-9]+[NT]$/) { print x, "NA", "NA", "excluded_unrecognized_name"; continue }
              p = substr(x, 1, length(x)-1); c = substr(x, length(x)); o = (c == "N") ? "T" : "N"
              print x, p, c, (((p, o) in has) ? "paired" : "unpaired_no_" o)
          } }
' "$OUT/sample_names.txt" > "$OUT/sample_pairing.tsv"
awk -F'\t' 'NR > 1 && $4 != "paired" { print "WARNING: " $1 " not in tumor/normal comparison (" $4 ")" > "/dev/stderr" }' \
    "$OUT/sample_pairing.tsv"
awk -F'\t' 'NR > 1 { n[$4]++ } END { for (k in n) print "  " k ": " n[k] }' "$OUT/sample_pairing.tsv"

echo "[2/5] Target exons (+/- ${PAD} bp) from GTF: $GENES"
BED="$OUT/target_exons.bed"
awk -v genes="$GENES" -v pad="$PAD" '
    BEGIN { FS=OFS="\t"; n=split(genes, g, " "); for (i=1; i<=n; i++) want[g[i]]=1 }
    ($3 == "gene" || $3 == "exon") && match($9, /gene_name "[^"]+"/) {
        name = substr($9, RSTART+11, RLENGTH-12)
        if (!(name in want)) next
        if ($3 == "gene") { seen[name]++; next }
        s = $4 - 1 - pad; if (s < 0) s = 0
        print $1, s, $5 + pad, name
    }
    END { for (x in want) if (seen[x] != 1) { print "gene " x " found " seen[x]+0 " times in GTF" > "/dev/stderr"; bad=1 }
          exit bad }
' "$GTF" | sort -k1,1 -k2,2n | awk '
    BEGIN { FS=OFS="\t" }
    $1 == c && $2 <= e { if ($3 > e) e = $3; if (index("," g ",", "," $4 ",") == 0) g = g "," $4; next }
    { if (c != "") print c, s, e, g; c=$1; s=$2; e=$3; g=$4 }
    END { if (c != "") print c, s, e, g }
' > "$BED"
awk -F'\t' '{ b[$4] += $3 - $2 } END { for (k in b) print "  " k ": " b[k] " target bases" }' "$BED"

echo "[3/5] Joint pileup of all BAMs"
# --ignore-RG: one BAM = one sample, named by the path given below; checked afterwards.
# -L equals -d: at most MAX_DEPTH reads per file enter the pileup, so the average
# per-file depth can never exceed -L and indel detection is never skipped for depth.
# -p applies the indel thresholds (-m >= 2 gapped reads, -F >= 5% gapped reads) per
# sample; by default they are applied to reads pooled over all BAMs, which hides an
# indel present in only one of ~20 samples.
# ADT (total AD over all alleles) is computed before the multiallelic split, so
# VAF = alt / ADT stays correct for split records.
PILEUP="$OUT/pileup.split.bcf"
bcftools mpileup -Ou --ignore-RG -f "$FASTA" -R "$BED" \
        -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR \
        -q "$MIN_MAPQ" -Q "$MIN_BASEQ" -d "$MAX_DEPTH" -L "$MAX_DEPTH" -p \
        "${BAMS[@]}" \
    | bcftools +fill-tags -Ou -- -t 'FORMAT/ADT:1=int(smpl_sum(FORMAT/AD))' \
    | bcftools norm -Ou -f "$FASTA" -m -any \
    | bcftools view -Ob -o "$PILEUP.tmp"
vcf_samples="$(bcftools query -l "$PILEUP.tmp")"
[[ "$vcf_samples" == "$(printf '%s\n' "${BAMS[@]}")" ]] \
    || die "pileup sample columns do not match the BAM list one-to-one; see bcftools query -l $PILEUP.tmp"
bcftools reheader -s "$OUT/sample_names.txt" -o "$PILEUP" "$PILEUP.tmp"
rm -f "$PILEUP.tmp"
bcftools index -f "$PILEUP"

echo "[4/5] Per-site usable depth (same reads/filters as the variant counts)"
# depth = ADT at each SNV/reference site; target bases absent from the pileup have depth 0.
{
    printf 'chrom\tpos'; while read -r s; do printf '\t%s' "$s"; done < "$OUT/sample_names.txt"; echo
    bcftools query -i 'TYPE!="indel"' -f '%CHROM\t%POS[\t%ADT]\n' "$PILEUP" \
        | awk -F'\t' '!seen[$1 FS $2]++'
} > "$OUT/depth_per_site.tsv"
awk -v bed="$BED" -v cap="$MAX_DEPTH" '
    BEGIN { FS=OFS="\t"
            while ((getline line < bed) > 0) { split(line, b, "\t"); nb++; c[nb]=b[1]; s[nb]=b[2]; e[nb]=b[3]; g[nb]=b[4]
                                                size[b[4]] += b[3] - b[2]; if (!(b[4] in gi)) { gi[b[4]]=1; go[++ng]=b[4] } } }
    NR == 1 { for (j=3; j<=NF; j++) name[j] = $j; ns = NF; next }
    { gene = ""
      for (i=1; i<=nb; i++) if ($1 == c[i] && $2 > s[i] && $2 <= e[i]) { gene = g[i]; break }
      if (gene == "") next
      for (j=3; j<=ns; j++) { if ($j >= 10) ge10[gene, j]++; if ($j >= 30) ge30[gene, j]++; if ($j > mx[gene, j]) mx[gene, j] = $j } }
    END { print "target", "sample", "target_bases", "frac_depth_ge10", "frac_depth_ge30", "max_depth"
          for (k=1; k<=ng; k++) for (j=3; j<=ns; j++) {
              t = go[k]
              printf "%s\t%s\t%d\t%.3f\t%.3f\t%d\n", t, name[j], size[t], ge10[t, j]/size[t], ge30[t, j]/size[t], mx[t, j]
              if (mx[t, j] >= 0.9 * cap) printf "WARNING: %s %s depth %d is near MAX_DEPTH=%d; reads may be capped\n", t, name[j], mx[t, j], cap > "/dev/stderr"
          } }
' "$OUT/depth_per_site.tsv" > "$OUT/depth_summary.tsv"
gzip -f "$OUT/depth_per_site.tsv"

echo "[5/5] Candidate alleles and tumor/normal comparison"
bcftools view -Oz -o "$OUT/candidates.vcf.gz" \
    -i "ALT!=\"<*>\" && MAX(FMT/AD[*:1])>=$MIN_ALT" "$PILEUP"
bcftools index -f "$OUT/candidates.vcf.gz"
n_candidates="$(bcftools view -H "$OUT/candidates.vcf.gz" | wc -l)"
echo "Candidate alleles: $n_candidates"

# Long table: one row per candidate allele x sample. After the split, AD = REF,ALT
# (true REF count) and ADT = all counted alleles; VAF is NA without coverage.
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE\t%AD\t%ADT\t%ADF\t%ADR]\n' "$OUT/candidates.vcf.gz" \
| awk -v bed="$BED" '
    BEGIN { FS=OFS="\t"
            while ((getline line < bed) > 0) { split(line, b, "\t"); n++; c[n]=b[1]; s[n]=b[2]; e[n]=b[3]; g[n]=b[4] }
            print "chrom", "pos", "ref", "alt", "target", "sample", "depth", "ref_reads", "alt_reads", "vaf", "alt_fwd", "alt_rev" }
    { gene = "."
      for (i=1; i<=n; i++) if ($1 == c[i] && $2 > s[i] && $2 <= e[i]) { gene = g[i]; break }
      for (j=5; j<NF; j+=5) {
          split($(j+1), ad, ","); split($(j+3), f, ","); split($(j+4), r, ",")
          dp = $(j+2) + 0; vaf = (dp > 0) ? sprintf("%.3f", ad[2] / dp) : "NA"
          print $1, $2, $3, $4, gene, $j, dp, ad[1]+0, ad[2]+0, vaf, f[2]+0, r[2]+0
      } }
' > "$OUT/allele_counts.long.tsv"

# Pair table: patients marked "paired" in sample_pairing.tsv; rows where N or T has >= MIN_ALT alt reads.
# p_zero_normal_alt = (1 - T_vaf)^N_depth: chance of seeing 0 alt reads in the normal if the
# allele were present there at the tumor VAF (independent-sampling model, RNA VAF only).
# Labels:
#   candidate_tumor_only           T_alt >= MIN_ALT, N_alt = 0, p_zero_normal_alt <= MAX_P0
#   tumor_alt_normal_underpowered  T_alt >= MIN_ALT, N_alt = 0, but normal depth too low to say
#   tumor_alt_normal_no_coverage   T_alt >= MIN_ALT, normal depth 0
#   tumor_alt_normal_low_alt       T_alt >= MIN_ALT, 0 < N_alt < MIN_ALT
#   alt_in_both                    both >= MIN_ALT (germline, shared artifact, or editing)
#   normal_alt_only                only N >= MIN_ALT
awk -v min_alt="$MIN_ALT" -v max_p0="$MAX_P0" '
    BEGIN { FS=OFS="\t" }
    FNR == NR { if (FNR > 1 && $4 == "paired") { role[$1] = $3; pt_of[$1] = $2; if (!($2 in pi)) { pi[$2]=1; po[++np]=$2 } } next }
    FNR == 1 { next }
    ($6 in role) { site = $1 OFS $2 OFS $3 OFS $4 OFS $5
                   if (!(site in si)) { si[site]=1; so[++ns]=site }
                   k = site SUBSEP pt_of[$6] SUBSEP role[$6]
                   dp[k] = $7; alt[k] = $9; vaf[k] = $10; fw[k] = $11; rv[k] = $12 }
    END {
      print "chrom", "pos", "ref", "alt", "target", "patient", "N_depth", "N_alt", "N_vaf", "T_depth", "T_alt", "T_vaf", "T_alt_fwd", "T_alt_rev", "p_zero_normal_alt", "label"
      for (i=1; i<=ns; i++) for (q=1; q<=np; q++) {
          site = so[i]; p = po[q]; kn = site SUBSEP p SUBSEP "N"; kt = site SUBSEP p SUBSEP "T"
          na = alt[kn] + 0; ta = alt[kt] + 0; nd = dp[kn] + 0
          if (na < min_alt && ta < min_alt) continue
          p0 = (vaf[kt] != "NA" && vaf[kt] != "") ? sprintf("%.3g", (1 - vaf[kt]) ^ nd) : "NA"
          if (ta >= min_alt) {
              if (na == 0) label = (nd == 0) ? "tumor_alt_normal_no_coverage" : ((p0 + 0 <= max_p0) ? "candidate_tumor_only" : "tumor_alt_normal_underpowered")
              else label = (na < min_alt) ? "tumor_alt_normal_low_alt" : "alt_in_both"
          } else label = "normal_alt_only"
          print site, p, nd, na, (vaf[kn] == "" ? "NA" : vaf[kn]), dp[kt] + 0, ta, (vaf[kt] == "" ? "NA" : vaf[kt]), fw[kt] + 0, rv[kt] + 0, p0, label
      } }
' "$OUT/sample_pairing.tsv" "$OUT/allele_counts.long.tsv" > "$OUT/tumor_normal_pairs.tsv"

echo "Done. Outputs in $OUT:"
echo "  sample_pairing.tsv        which samples were paired, unpaired or excluded"
echo "  depth_summary.tsv         fraction of target bases with usable depth >= 10 / >= 30"
echo "  depth_per_site.tsv.gz     usable depth at every covered target base"
echo "  candidates.vcf.gz         candidate alleles (annotate with VEP)"
echo "  allele_counts.long.tsv    depth, ref/alt reads, VAF, alt strand counts per sample"
echo "  tumor_normal_pairs.tsv    N vs T per paired patient, with review labels"
echo "candidate_tumor_only rows (review each one; not confirmed somatic mutations):"
awk -F'\t' 'NR == 1 || $16 == "candidate_tumor_only"' "$OUT/tumor_normal_pairs.tsv"
