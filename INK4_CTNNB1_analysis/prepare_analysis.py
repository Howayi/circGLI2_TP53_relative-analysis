"""Inspect supplied RNA allele evidence and annotate against public gene references locally."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
import math
import struct
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
KEY = ["chrom", "pos", "ref", "alt"]
BASES = "TCAG"
CODONS = dict(zip(("".join(x) for x in itertools.product(BASES, repeat=3)),
                  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
AA3 = dict(zip("ACDEFGHIKLMNPQRSTVWY*", ["Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr", "Ter"]))


def complement(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def write_tsv(frame, filename):
    frame.to_csv(ROOT / "source_data" / filename, sep="\t", index=False, na_rep="NA")


def load_transcript(gene, transcript=None, sequence_name=None):
    obj = json.loads((ROOT / "references" / f"{gene}.lookup.json").read_text())
    tid = transcript or obj["canonical_transcript"].split(".")[0]
    tr = next(t for t in obj["Transcript"] if t["id"] == tid)
    assert obj["assembly_name"] == "GRCh38"
    strand = tr["strand"]
    exons = sorted(tr["Exon"], key=lambda e: e["start"], reverse=strand == -1)
    positions = []
    for e in exons:
        positions.extend(range(e["start"], e["end"] + 1) if strand == 1 else range(e["end"], e["start"] - 1, -1))
    translation = tr["Translation"]
    first = translation["start"] if strand == 1 else translation["end"]
    last = translation["end"] if strand == 1 else translation["start"]
    start, end = positions.index(first), positions.index(last)
    mapped_length = end - start + 1
    assert mapped_length in (translation["length"] * 3, translation["length"] * 3 + 3)
    seq_path = ROOT / "references" / f"{sequence_name or gene}.cds.json"
    seq = json.loads(seq_path.read_text())["seq"].upper() if seq_path.exists() else None
    if seq is not None:
        assert len(seq) in (mapped_length, mapped_length + 3), (gene, len(seq), mapped_length)
        assert seq[:3] == "ATG"
    stop_extension = 3 if mapped_length == translation["length"] * 3 and (seq is None or len(seq) == mapped_length + 3) else 0
    coding_positions = positions[start:end + 1 + stop_extension]
    return {"gene": gene, "object": obj, "tr": tr, "positions": positions,
            "cdna": {p: i + 1 for i, p in enumerate(positions)},
            "cds": {p: i + 1 for i, p in enumerate(coding_positions)},
            "cds_start_index": start, "cds_end_index": end,
            "seq": seq, "exons": exons}


def context(model, pos):
    if pos not in model["cdna"]:
        return "outside_canonical_exon"
    index = model["cdna"][pos] - 1
    if pos in model["cds"]:
        return "CDS"
    return "5_prime_UTR" if index < model["cds_start_index"] else "3_prime_UTR"


def annotate_snv(model, pos, ref, alt):
    out = {"canonical_context": context(model, pos), "canonical_consequence": "",
           "cdna_position": model["cdna"].get(pos), "cds_position": model["cds"].get(pos),
           "protein_position": None, "protein_change": "", "codon_change": "", "reference_check": "not_CDS_SNV"}
    if out["canonical_context"] != "CDS":
        out["canonical_consequence"] = out["canonical_context"]
        return out
    if len(ref) != 1 or len(alt) != 1:
        out["canonical_consequence"] = "frameshift_candidate" if (len(alt) - len(ref)) % 3 else "inframe_indel_candidate"
        return out
    cdspos = out["cds_position"]
    out["protein_position"] = (cdspos - 1) // 3 + 1
    if model["seq"] is None:
        out["canonical_consequence"] = "coding_SNV_unresolved"
        out["reference_check"] = "CDS_sequence_unavailable"
        return out
    ref_t, alt_t = (ref, alt) if model["tr"]["strand"] == 1 else (complement(ref), complement(alt))
    assert model["seq"][cdspos - 1] == ref_t, (model["gene"], pos, ref, model["seq"][cdspos - 1])
    codon_start = ((cdspos - 1) // 3) * 3
    codon = model["seq"][codon_start:codon_start + 3]
    changed = list(codon)
    changed[(cdspos - 1) % 3] = alt_t
    alt_codon = "".join(changed)
    aa, alt_aa = CODONS[codon], CODONS[alt_codon]
    consequence = "synonymous_variant" if aa == alt_aa else "stop_gained" if alt_aa == "*" else "stop_lost" if aa == "*" else "missense_variant"
    if cdspos <= 3 and aa != alt_aa:
        consequence = "start_lost_candidate"
    out.update(canonical_consequence=consequence,
               protein_change=f"p.{AA3[aa]}{out['protein_position']}{'=' if aa == alt_aa else AA3[alt_aa]}",
               codon_change=f"{codon}>{alt_codon}", reference_check="match")
    return out


def inspect_bcf(path):
    with gzip.open(path, "rb") as f:
        assert f.read(5) == b"BCF\x02\x02"
        hlen = struct.unpack("<I", f.read(4))[0]
        header = f.read(hlen).decode().rstrip("\x00")
        count = 0
        while True:
            lengths = f.read(8)
            if not lengths:
                break
            assert len(lengths) == 8
            shared, individual = struct.unpack("<II", lengths)
            assert len(f.read(shared + individual)) == shared + individual
            count += 1
    return count, header


def zero_probability(alt, depth, normal_depth):
    if normal_depth == 0 or depth == 0:
        return 1.0
    if alt == depth:
        return 0.0
    return math.exp(normal_depth * math.log1p(-alt / depth))


def main(input_dir):
    p = Path(input_dir)
    (ROOT / "source_data").mkdir(exist_ok=True)
    a = pd.read_csv(p / "allele_counts.long.tsv", sep="\t")
    pairs = pd.read_csv(p / "tumor_normal_pairs.tsv", sep="\t")
    pairing = pd.read_csv(p / "sample_pairing.tsv", sep="\t", dtype={"patient": str})
    samples = (p / "sample_names.txt").read_text().splitlines()
    bed = pd.read_csv(p / "target_exons.bed", sep="\t", header=None, names=["chrom", "start", "end", "target"])
    depth = pd.read_csv(p / "depth_per_site.tsv.gz", sep="\t")
    coverage = pd.read_csv(p / "depth_summary.tsv", sep="\t")
    assert pairing["sample"].tolist() == samples
    assert depth.columns.tolist() == ["chrom", "pos"] + samples
    assert not a.duplicated(KEY + ["sample"]).any()
    assert (a.alt_fwd + a.alt_rev == a.alt_reads).all()
    assert (a.ref_reads + a.alt_reads <= a.depth).all()
    expected_vaf = a.alt_reads / a.depth.replace(0, float("nan"))
    assert (expected_vaf[a.depth > 0] - a.vaf[a.depth > 0]).abs().max() <= 0.00050001
    assert a.loc[a.depth == 0, "vaf"].isna().all()
    variants = a[KEY + ["target"]].drop_duplicates().copy()
    a["vaf_exact"] = expected_vaf
    a["variant_id"] = a.chrom + ":" + a.pos.astype(str) + ":" + a.ref + ">" + a.alt
    variants["variant_id"] = variants.chrom + ":" + variants.pos.astype(str) + ":" + variants.ref + ">" + variants.alt

    vcf_rows, header = [], []
    with gzip.open(p / "candidates.vcf.gz", "rt") as f:
        for line in f:
            fields = line.rstrip().split("\t")
            if line.startswith("#"):
                header.append(line.rstrip())
                if line.startswith("#CHROM"):
                    assert fields[9:] == samples
                continue
            chrom, pos, _, ref, alt = fields[:5]
            tags = fields[8].split(":")
            info = dict(x.split("=", 1) for x in fields[7].split(";") if "=" in x)
            vcf_rows.append({"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt,
                             **{f"vcf_{tag}": info.get(tag) for tag in ["VDB", "RPBZ", "MQBZ", "BQBZ", "SCBZ"]}})
            observed = a[(a.chrom == chrom) & (a.pos == int(pos)) & (a.ref == ref) & (a.alt == alt)].set_index("sample")
            assert len(observed) == len(samples)
            for sample, encoded in zip(samples, fields[9:]):
                fmt = dict(zip(tags, encoded.split(":")))
                ad = [int(v) for v in fmt["AD"].split(",")]
                assert ad == [int(observed.loc[sample, "ref_reads"]), int(observed.loc[sample, "alt_reads"])]
                assert int(fmt["ADT"]) == int(observed.loc[sample, "depth"])
    assert len(vcf_rows) == len(variants)
    variants = variants.merge(pd.DataFrame(vcf_rows), on=KEY, validate="one_to_one")

    models = {g: load_transcript(g) for g in ["CTNNB1", "CDKN2A", "CDKN2B"]}
    arf = load_transcript("CDKN2A", "ENST00000579755", "CDKN2A_ARF")
    annotated = []
    blocks = []
    for gene, model in models.items():
        tr = model["tr"]
        for exon_number, exon in enumerate(model["exons"], 1):
            poslist = list(range(exon["start"], exon["end"] + 1))
            for ctx in ["5_prime_UTR", "CDS", "3_prime_UTR"]:
                selected = [x for x in poslist if context(model, x) == ctx]
                if selected:
                    blocks.append({"target": gene, "transcript": f"{tr['id']}.{tr['version']}", "exon_number": exon_number,
                                   "context": ctx, "genomic_start": min(selected), "genomic_end": max(selected),
                                   "cdna_start": min(model["cdna"][x] for x in selected),
                                   "cdna_end": max(model["cdna"][x] for x in selected)})
    for row in variants.to_dict("records"):
        model = models[row["target"]]
        ann = annotate_snv(model, row["pos"], row["ref"], row["alt"])
        coding_in = []
        for tr in model["object"]["Transcript"]:
            if "Translation" not in tr:
                continue
            bounds = tr["Translation"]
            if min(bounds["start"], bounds["end"]) <= row["pos"] <= max(bounds["start"], bounds["end"]) and any(e["start"] <= row["pos"] <= e["end"] for e in tr["Exon"]):
                coding_in.append(f"{tr['id']}.{tr['version']}")
        ann["coding_transcripts_at_anchor"] = ";".join(coding_in)
        ann["canonical_transcript"] = f"{model['tr']['id']}.{model['tr']['version']}"
        ann["ARF_context"] = ""
        ann["ARF_protein_change"] = ""
        if row["target"] == "CDKN2A":
            arf_ann = annotate_snv(arf, row["pos"], row["ref"], row["alt"])
            ann["ARF_context"] = arf_ann["canonical_context"]
            ann["ARF_protein_change"] = arf_ann["protein_change"]
        annotated.append({**row, **ann})
    annotations = pd.DataFrame(annotated)
    a = a.merge(annotations[KEY + ["canonical_context", "canonical_consequence", "protein_change"]], on=KEY, validate="many_to_one")
    pairs = pairs.merge(annotations, on=KEY + ["target"], validate="many_to_one")
    pairs["patient"] = pairs.patient.astype(str)
    for side in ["N", "T"]:
        pairs[f"{side}_vaf_exact"] = pairs[f"{side}_alt"] / pairs[f"{side}_depth"].replace(0, float("nan"))
        side_rows = a.copy()
        side_rows["patient"] = side_rows["sample"].str[:-1]
        side_rows = side_rows[side_rows["sample"].str.endswith(side)]
        merged = pairs.merge(side_rows[KEY + ["patient", "depth", "alt_reads"]], on=KEY + ["patient"], validate="many_to_one")
        assert len(merged) == len(pairs)
        assert (merged[f"{side}_depth"] == merged.depth).all()
        assert (merged[f"{side}_alt"] == merged.alt_reads).all()
    pairs["p_zero_exact"] = [zero_probability(r.T_alt, r.T_depth, r.N_depth) for r in pairs.itertuples()]
    pairs["p_zero_rounded_vaf"] = [zero_probability(float(r.T_vaf) * r.T_depth, r.T_depth, r.N_depth) if pd.notna(r.T_vaf) else 1.0 for r in pairs.itertuples()]
    pairs["label_original"] = pairs["label"]
    pairs["p_zero_normal_alt_original"] = pairs["p_zero_normal_alt"]
    eligible = (pairs.T_alt >= 3) & (pairs.N_alt == 0) & (pairs.N_depth > 0)
    assert pairs.loc[eligible, "label"].isin(["candidate_tumor_only", "tumor_alt_normal_underpowered"]).all()
    pairs.loc[eligible, "label"] = pairs.loc[eligible, "p_zero_exact"].map(
        lambda p0: "candidate_tumor_only" if p0 <= 0.05 else "tumor_alt_normal_underpowered")
    pairs["p_zero_normal_alt"] = pairs["p_zero_exact"]
    pairs["label_changed_by_exact_p0"] = pairs.label.ne(pairs.label_original)
    pairs["original_candidate_removed"] = pairs.label_original.eq("candidate_tumor_only") & ~pairs.label.eq("candidate_tumor_only")
    pairs["new_candidate_added"] = ~pairs.label_original.eq("candidate_tumor_only") & pairs.label.eq("candidate_tumor_only")
    pairs["one_direction_ALT"] = (pairs.T_alt > 0) & ((pairs.T_alt_fwd == 0) | (pairs.T_alt_rev == 0))
    pairs["normal_depth_lt20"] = pairs.N_depth < 20
    pairs["tumor_depth_lt20"] = pairs.T_depth < 20
    pairs["tumor_ALT_lt10"] = pairs.T_alt < 10
    pairs["review_high_support"] = ((pairs.label == "candidate_tumor_only") & (pairs.T_alt >= 10) & (pairs.N_depth >= 20) & (pairs.T_depth >= 20) & (pairs.T_alt_fwd >= 2) & (pairs.T_alt_rev >= 2))
    pairs["p0_exact_above_threshold"] = (pairs.label == "candidate_tumor_only") & (pairs.p_zero_exact > 0.05)
    candidates = pairs[pairs.label == "candidate_tumor_only"].sort_values(["patient", "T_alt"], ascending=[True, False])

    coverage_recomputed = []
    for gene in models:
        b = bed[bed.target == gene]
        selected = pd.Series(False, index=depth.index)
        for r in b.itertuples():
            selected |= (depth.chrom == r.chrom) & (depth.pos > r.start) & (depth.pos <= r.end)
        assert int(selected.sum()) == int((b.end - b.start).sum())
        for sample in samples:
            values = depth.loc[selected, sample]
            reference = coverage[(coverage.target == gene) & (coverage["sample"] == sample)].iloc[0]
            frac10, frac30 = (values >= 10).mean(), (values >= 30).mean()
            assert abs(frac10 - reference.frac_depth_ge10) <= 0.00050001
            assert abs(frac30 - reference.frac_depth_ge30) <= 0.00050001
            assert int(values.max()) == int(reference.max_depth)
            coverage_recomputed.append({"target": gene, "sample": sample, "target_bases": len(values),
                                        "frac_depth_ge10": frac10, "frac_depth_ge30": frac30, "max_depth": int(values.max()),
                                        "mean_depth": values.mean(), "median_depth": values.median()})
    coverage_clean = pd.DataFrame(coverage_recomputed).merge(pairing, on="sample", validate="many_to_one")
    bcf_count, bcf_header = inspect_bcf(p / "pileup.split.bcf")
    assert next(line for line in bcf_header.splitlines() if line.startswith("#CHROM")).split("\t")[9:] == samples

    write_tsv(a, "allele_evidence_annotated.tsv")
    write_tsv(annotations, "variant_annotations.tsv")
    write_tsv(pairs, "paired_evidence_annotated.tsv")
    write_tsv(candidates, "tumor_only_candidates_review.tsv")
    write_tsv(pairs[pairs.label_changed_by_exact_p0].sort_values(["patient", "pos"]), "P0_label_changes.tsv")
    write_tsv(coverage_clean, "coverage_recomputed.tsv")
    write_tsv(pairing, "sample_metadata.tsv")
    write_tsv(pd.DataFrame(blocks), "canonical_transcript_blocks.tsv")
    write_tsv(pairs.groupby(["target", "patient", "label"]).size().reset_index(name="records"), "status_counts.tsv")
    write_tsv(annotations.groupby(["target", "canonical_consequence"]).size().reset_index(name="variant_records"), "consequence_counts.tsv")
    write_tsv(pd.DataFrame([{"file": x.name, "bytes": x.stat().st_size, "sha256": hashlib.sha256(x.read_bytes()).hexdigest()} for x in sorted(p.iterdir()) if x.is_file()]), "input_manifest.tsv")
    summary = {
        "input_samples": len(samples), "paired_patients": len(pairing.loc[pairing.pair_status == "paired", "patient"].unique()),
        "unpaired_samples": pairing.loc[pairing.pair_status != "paired", "sample"].tolist(),
        "unique_variant_records": len(annotations), "unique_variants_by_gene": annotations.groupby("target").size().to_dict(),
        "paired_evidence_records": len(pairs), "tumor_only_pair_records": len(candidates),
        "tumor_only_unique_variants": len(candidates[KEY].drop_duplicates()),
        "original_tumor_only_pair_records": int(pairs.label_original.eq("candidate_tumor_only").sum()),
        "P0_labels_changed": int(pairs.label_changed_by_exact_p0.sum()),
        "P0_candidates_removed": int(pairs.original_candidate_removed.sum()),
        "P0_candidates_added": int(pairs.new_candidate_added.sum()),
        "tumor_only_missense_records": int(candidates.canonical_consequence.eq("missense_variant").sum()),
        "tumor_only_single_direction": int(candidates.one_direction_ALT.sum()),
        "tumor_only_ALT_lt10": int(candidates.tumor_ALT_lt10.sum()),
        "tumor_only_low_normal_depth": int(candidates.normal_depth_lt20.sum()),
        "tumor_only_low_tumor_depth": int(candidates.tumor_depth_lt20.sum()),
        "tumor_only_p0_exact_above_threshold": int(candidates.p0_exact_above_threshold.sum()),
        "high_support_pair_records": int(candidates.review_high_support.sum()),
        "allele_rows_with_other_ALT_counts": int((a.ref_reads + a.alt_reads < a.depth).sum()),
        "bcf_records": bcf_count, "candidate_vcf_records": len(vcf_rows), "depth_rows": len(depth),
        "CDS_reference_matches": int((annotations.reference_check == "match").sum()),
        "paired_label_counts": pairs.label.value_counts().to_dict(),
        "tumor_only_canonical_consequence_counts": candidates.canonical_consequence.value_counts().to_dict(),
        "reference_retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "annotation_scope": "Local canonical-transcript annotation against downloaded Ensembl GRCh38 references; CDKN2A ARF additionally inspected; not full VEP/clinical pathogenicity annotation.",
        "checks": "PASS: VCF sample identity, all allele counts/ADT, strand sums, paired counts, VAF rounding, all coverage-summary cells, BCF structure and samples; labels corrected in both directions using exact integer-count P0"
    }
    (ROOT / "validation.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print("\nTUMOR-ONLY ANNOTATIONS\n", candidates[["patient", "pos", "ref", "alt", "T_alt", "T_depth", "N_depth", "canonical_context", "protein_change", "p_zero_exact", "p0_exact_above_threshold"]].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Directory produced by ink4_ctnnb1_check.sh")
    main(parser.parse_args().input)
