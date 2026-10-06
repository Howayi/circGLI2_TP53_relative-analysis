"""Check evidence/export integrity and package the finished result; no graphics are drawn."""
import argparse
import hashlib
import json
import struct
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import pandas as pd
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parent


def main(input_dir):
    summary = json.loads((ROOT / "validation.json").read_text())
    ann = pd.read_csv(ROOT / "source_data/variant_annotations.tsv", sep="\t")
    candidates = pd.read_csv(ROOT / "source_data/tumor_only_candidates_review.tsv", sep="\t")
    assert len(ann) == summary["unique_variant_records"] == 628
    assert len(candidates) == summary["tumor_only_pair_records"] == 34
    assert candidates[["chrom", "pos", "ref", "alt"]].drop_duplicates().shape[0] == 31
    assert candidates.canonical_consequence.eq("missense_variant").sum() == 9
    assert candidates.p0_exact_above_threshold.sum() == 0
    assert candidates.review_high_support.sum() == 6
    assert ann[ann.canonical_context == "CDS"].query("ref.str.len() == 1 and alt.str.len() == 1", engine="python").reference_check.eq("match").all()
    original = Path(input_dir)
    manifest = pd.read_csv(ROOT / "source_data/input_manifest.tsv", sep="\t")
    for row in manifest.itertuples():
        assert hashlib.sha256((original / row.file).read_bytes()).hexdigest() == row.sha256

    paired = pd.read_csv(ROOT / "source_data/paired_evidence_annotated.tsv", sep="\t")
    raw_pair = pd.read_csv(original / "tumor_normal_pairs.tsv", sep="\t")
    key = ["chrom", "pos", "ref", "alt", "target", "patient"]
    joined = raw_pair.merge(paired, on=key, suffixes=("_raw", "_corrected"), validate="one_to_one")
    assert len(joined) == len(raw_pair) == len(paired) == 727
    for column in set(raw_pair.columns) - set(key):
        left = joined[column + "_raw"]
        right = joined[column + "_original"] if column in ["label", "p_zero_normal_alt"] else joined[column + "_corrected"]
        assert (left.eq(right) | (left.isna() & right.isna())).all(), column
    eligible = (paired.T_alt >= 3) & (paired.N_alt == 0) & (paired.N_depth > 0)
    independent_p0 = (1 - paired.T_alt / paired.T_depth.replace(0, float("nan"))) ** paired.N_depth
    assert (independent_p0[eligible] - paired.loc[eligible, "p_zero_exact"]).abs().max() < 1e-11
    expected_candidate = eligible & independent_p0.le(0.05)
    assert paired.label.eq("candidate_tumor_only").eq(expected_candidate).all()
    assert paired.label_changed_by_exact_p0.sum() == 7
    assert paired.original_candidate_removed.sum() == 6
    assert paired.new_candidate_added.sum() == 1
    added = paired[paired.new_candidate_added].iloc[0]
    assert added.patient == 1 and added.pos == 41224972 and added.protein_change == "p.Ala87Val"
    assert 0.049 < added.p_zero_exact < 0.05
    expected_per_patient = {1: 16, 3: 0, 6: 4, 28: 3, 33: 4, 34: 0, 42: 2, 57: 5}
    actual_per_patient = paired.assign(is_candidate=expected_candidate).groupby("patient").is_candidate.sum().to_dict()
    assert actual_per_patient == expected_per_patient
    plot_qa = json.loads((ROOT / "plot_data_QA.json").read_text())
    assert plot_qa == {"comparison_records": 727, "tumor_only_records": 34, "heatmap_cells": 496,
                       "heatmap_candidate_borders": 34, "missense_points": 9}

    exports = pd.read_csv(ROOT / "figure_exports.tsv", sep="\t")
    assert len(exports) == 6
    rows = []
    for row in exports.itertuples():
        prefix = ROOT / "figures" / row.figure
        png = prefix.with_suffix(".png").read_bytes()
        assert png[:8] == b"\x89PNG\r\n\x1a\n"
        width, height = struct.unpack(">II", png[16:24])
        assert abs(width - row.width_mm / 25.4 * 300) < 2
        assert abs(height - row.height_mm / 25.4 * 300) < 2
        svg = ET.parse(prefix.with_suffix(".svg")).getroot()
        svg_text = " ".join("".join(t.itertext()) for t in svg.iter() if t.tag.endswith("}text"))
        assert len(svg_text) > 100
        assert not any("\u4e00" <= c <= "\u9fff" for c in svg_text)
        pdf = PdfReader(prefix.with_suffix(".pdf"))
        assert len(pdf.pages) == 1
        pdf_text = pdf.pages[0].extract_text()
        assert len(pdf_text) > 100
        assert not any("\u4e00" <= c <= "\u9fff" for c in pdf_text)
        for extension in ["png", "pdf", "svg", "tiff"]:
            assert prefix.with_suffix("." + extension).stat().st_size > 1000
        rows.append({"figure": row.figure, "png_pixels": [width, height], "pdf_pages": 1,
                     "svg_editable_text": True, "pdf_extractable_English_text": True})
    qa = {"status": "PASS", "original_inputs_unchanged": True,
          "all_coding_SNV_references_match": True,
          "exact_P0_relabelling": {"changed": 7, "removed": 6, "added": 1,
                                   "candidate_records": 34, "candidate_sites": 31, "missense_records": 9},
          "plotted_data": plot_qa,
          "images": rows,
          "figure_language": "English",
          "visual_review": "Manual visual review required after regenerating figures.",
          "limitations": "Exploratory RNA allele evidence; no BAM-level visual validation or DNA somatic calling."}
    (ROOT / "export_QA.json").write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding="utf-8")
    files = [x for x in ROOT.iterdir() if x.is_file() and x.suffix in {".py", ".R", ".md", ".json", ".tsv", ".txt"}]
    for folder in ["figures", "source_data", "references"]:
        files.extend(x for x in (ROOT / folder).rglob("*") if x.is_file())
    destination = ROOT / "INK4_CTNNB1_corrected_EN_results.zip"
    with zipfile.ZipFile(destination, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for file in sorted(files):
            z.write(file, str(Path("INK4_CTNNB1_corrected_EN_results") / file.relative_to(ROOT)))
    with zipfile.ZipFile(destination) as z:
        assert z.testzip() is None
    print(json.dumps({"QA": "PASS", "figures": len(rows), "bundle_files": len(files),
                      "bundle_bytes": destination.stat().st_size}, ensure_ascii=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Original checkpoint result directory")
    main(parser.parse_args().input)
