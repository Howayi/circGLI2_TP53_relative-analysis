#!/usr/bin/env python3
"""Matched TCGA PanCancer Atlas TP53 mutation vs PROGENy-derived p53/DDR footprint.

Inputs are the patient-level TP53 mutation table from the companion mutation
analysis and the frozen top-100 p53 weights from Bioconductor progeny 1.34.0.
Expression is downloaded from the matching cBioPortal PanCancer Atlas studies.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
from http.client import IncompleteRead, RemoteDisconnected
import json
import math
import os
from pathlib import Path
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / "results" / ".mplconfig"))
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


API = "https://www.cbioportal.org/api"
EXPR_SUFFIX = "_rna_seq_v2_mrna_median_all_sample_Zscores"
RNA_LIST_SUFFIX = "_rna_seq_v2_mrna"
MODEL_URL = "https://bioconductor.org/packages/3.23/bioc/src/contrib/progeny_1.34.0.tar.gz"
MODEL_VERSION = "progeny 1.34.0 (Bioconductor 3.23)"
MODEL_TARBALL_SHA256 = "3e8fe926eef24e587e4a1089e4efe0b81a8aa8bd385c91018edeffa9757791d1"
MODEL_WEIGHTS_SHA256 = "c2b9aba4a7793d1ab7b5ebccda4001dd66a6e654b39a722f4bc8ec0b38c945d9"
GDC_EXPR_UUID = "3586c0da-64d0-4b74-a449-5ff4d9136611"
GDC_EXPR_MD5 = "02e72c33071307ff6570621480d3c90b"
GDC_EXPR_FILENAME = "EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"
USER_AGENT = "TP53-PROGENy-pan-cancer-analysis/1.0"
GENE_CHUNK_SIZE = 5
WT_COLOR = "#4C78A8"
MUT_COLOR = "#D7301F"
# Seven symbols in the 2018 PROGENy model were subsequently renamed. Entrez IDs
# preserve gene identity across the model and PanCanAtlas expression profile.
LEGACY_SYMBOL_TO_ENTREZ = {
    "WDR63": 126820,     # current symbol DNAI3
    "TMEM27": 57393,     # current symbol CLTRN
    "C17orf82": 388407,  # current symbol LINC02875
    "TEX37": 200523,     # current symbol SPMIP9
    "C2orf66": 401027,   # current cBioPortal capitalization C2ORF66
    "RGAG4": 340526,     # current symbol RTL5
    "LOC400710": 400710, # current cBioPortal symbol ZNF473CR
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=SCRIPT_DIR / "inputs")
    parser.add_argument("--out", type=Path, default=SCRIPT_DIR / "results")
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260921)
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Use existing compressed API caches only; do not access cBioPortal.",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def post_json(url: str, body: object, retries: int = 7) -> object:
    payload = json.dumps(body, separators=(",", ":")).encode("utf-8")
    request = Request(
        url,
        data=payload,
        headers={
            "Accept": "application/json",
            "Content-Type": "application/json",
            "User-Agent": USER_AGENT,
        },
        method="POST",
    )
    for attempt in range(retries):
        try:
            with urlopen(request, timeout=60) as response:
                return json.load(response)
        except (
            HTTPError,
            URLError,
            TimeoutError,
            ConnectionError,
            ConnectionResetError,
            IncompleteRead,
            RemoteDisconnected,
            json.JSONDecodeError,
        ) as exc:
            retryable = not isinstance(exc, HTTPError) or exc.code in {
                408,
                429,
                500,
                502,
                503,
                504,
            }
            if not retryable or attempt == retries - 1:
                raise
            wait_seconds = min(30, 2**attempt)
            print(f"  API retry {attempt + 1}/{retries - 1} in {wait_seconds}s: {exc}")
            time.sleep(wait_seconds)
    raise RuntimeError("unreachable")


def read_gzip_json(path: Path) -> object:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def write_gzip_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", compresslevel=6) as handle:
        json.dump(value, handle, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def chunks(values: list[int], size: int) -> list[list[int]]:
    return [values[start : start + size] for start in range(0, len(values), size)]


def expression_cache_request(
    study_id: str, sample_ids: list[str], entrez_ids: list[int]
) -> tuple[dict, str]:
    metadata = {
        "cache_schema": 1,
        "study_id": study_id,
        "expression_profile_id": study_id + EXPR_SUFFIX,
        "rna_sample_list_id": study_id + RNA_LIST_SUFFIX,
        "analysis_sample_ids": sorted(sample_ids),
        "entrez_gene_ids": sorted(int(value) for value in entrez_ids),
        "projection": "DETAILED",
    }
    serialized = json.dumps(metadata, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return metadata, hashlib.sha256(serialized).hexdigest()


def split_sample_ids(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return [item for item in str(value).split(";") if item]


def sample_type_code(sample_id: str) -> str:
    parts = sample_id.split("-")
    return parts[3][:2] if len(parts) >= 4 and len(parts[3]) >= 2 else "NA"


def sample_priority(code: str) -> int:
    order = {
        "01": 1,  # primary solid tumor
        "03": 2,  # primary blood-derived cancer, peripheral blood
        "09": 3,  # primary blood-derived cancer, bone marrow
        "02": 4,  # recurrent solid tumor
        "04": 5,  # recurrent blood-derived cancer
        "05": 6,  # additional new primary
        "06": 7,  # metastatic
        "07": 8,  # additional metastatic
        "08": 9,  # human tumor original cells
    }
    return order.get(code, 99)


def fetch_gene_mapping(
    genes: list[str], cache_path: Path, force: bool, offline: bool
) -> list[dict]:
    if cache_path.exists() and not force:
        return read_gzip_json(cache_path)
    if offline:
        raise FileNotFoundError(f"Offline cache missing: {cache_path}")
    result = post_json(
        f"{API}/genes/fetch?geneIdType=HUGO_GENE_SYMBOL&projection=SUMMARY",
        genes,
    )
    write_gzip_json(cache_path, result)
    return result


def fetch_study_expression(
    study_id: str,
    sample_ids: list[str],
    entrez_ids: list[int],
    cache_path: Path,
    force: bool,
    offline: bool,
) -> list[dict]:
    request_metadata, request_sha256 = expression_cache_request(
        study_id, sample_ids, entrez_ids
    )
    if cache_path.exists() and not force:
        cached = read_gzip_json(cache_path)
        if (
            isinstance(cached, dict)
            and cached.get("request_sha256") == request_sha256
            and isinstance(cached.get("records"), list)
        ):
            return cached["records"]
        if offline:
            raise ValueError(f"Offline cache fingerprint mismatch: {cache_path}")
    if offline:
        raise FileNotFoundError(f"Offline cache missing: {cache_path}")
    profile_id = study_id + EXPR_SUFFIX
    sample_list_id = study_id + RNA_LIST_SUFFIX
    endpoint = f"{API}/molecular-profiles/{profile_id}/molecular-data/fetch?projection=DETAILED"
    records: list[dict] = []
    gene_chunks = chunks(entrez_ids, GENE_CHUNK_SIZE)
    for chunk_index, gene_chunk in enumerate(gene_chunks, start=1):
        print(f"  expression chunk {chunk_index}/{len(gene_chunks)}", flush=True)
        response = post_json(
            endpoint,
            {"sampleListId": sample_list_id, "entrezGeneIds": gene_chunk},
        )
        records.extend(response)
        time.sleep(0.25)
    write_gzip_json(
        cache_path,
        {
            "request": request_metadata,
            "request_sha256": request_sha256,
            "records": records,
        },
    )
    return records


def choose_one_sample_per_patient(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["sample_type_code"] = frame["sample_id"].map(sample_type_code)
    frame["sample_priority"] = frame["sample_type_code"].map(sample_priority)
    frame = frame.sort_values(
        ["patient_id", "sample_priority", "n_signature_genes_observed", "sample_id"],
        ascending=[True, True, False, True],
        kind="mergesort",
    )
    return frame.drop_duplicates("patient_id", keep="first").reset_index(drop=True)


def fast_cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    y_sorted = np.sort(y)
    less = np.searchsorted(y_sorted, x, side="left")
    at_most = np.searchsorted(y_sorted, x, side="right")
    u = np.sum(less + 0.5 * (at_most - less))
    return float(2 * u / (len(x) * len(y)) - 1)


def mann_whitney_asymptotic(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    n_x, n_y = len(x), len(y)
    combined = np.concatenate([x, y])
    ranks = pd.Series(combined).rank(method="average").to_numpy()
    u_x = float(ranks[:n_x].sum() - n_x * (n_x + 1) / 2)
    delta = float(2 * u_x / (n_x * n_y) - 1)
    _, tie_counts = np.unique(combined, return_counts=True)
    n_total = n_x + n_y
    tie_term = float(np.sum(tie_counts**3 - tie_counts))
    variance = n_x * n_y / 12 * (
        (n_total + 1) - tie_term / (n_total * (n_total - 1))
    )
    if variance <= 0:
        return u_x, 1.0, delta
    numerator = max(0.0, abs(u_x - n_x * n_y / 2) - 0.5)
    z_value = numerator / math.sqrt(variance)
    p_value = math.erfc(z_value / math.sqrt(2))
    return u_x, p_value, delta


def bootstrap_delta_ci(
    x: np.ndarray, y: np.ndarray, n_boot: int, rng: np.random.Generator
) -> tuple[float, float]:
    estimates = np.empty(n_boot, dtype=float)
    for index in range(n_boot):
        x_boot = x[rng.integers(0, len(x), len(x))]
        y_boot = y[rng.integers(0, len(y), len(y))]
        estimates[index] = fast_cliffs_delta(x_boot, y_boot)
    low, high = np.percentile(estimates, [2.5, 97.5])
    return float(low), float(high)


def bh_adjust(p_values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid = p_values.dropna()
    if valid.empty:
        return result
    order = np.argsort(valid.to_numpy())
    ranked = valid.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    result.loc[valid.index[order]] = adjusted
    return result


def summarize_groups(
    scores: pd.DataFrame, n_boot: int, seed: int, min_group: int = 10
) -> pd.DataFrame:
    rows: list[dict] = []
    for cancer, group in scores.groupby("cancer_type", sort=True):
        mutant = group.loc[group["tp53_status"].eq("MUT"), "activity_z_within_cancer"].to_numpy()
        wildtype = group.loc[group["tp53_status"].eq("WT"), "activity_z_within_cancer"].to_numpy()
        row = {
            "cancer_type": cancer,
            "n_matched": len(group),
            "n_mut": len(mutant),
            "n_wt": len(wildtype),
            "median_mut": float(np.median(mutant)) if len(mutant) else np.nan,
            "q1_mut": float(np.quantile(mutant, 0.25)) if len(mutant) else np.nan,
            "q3_mut": float(np.quantile(mutant, 0.75)) if len(mutant) else np.nan,
            "median_wt": float(np.median(wildtype)) if len(wildtype) else np.nan,
            "q1_wt": float(np.quantile(wildtype, 0.25)) if len(wildtype) else np.nan,
            "q3_wt": float(np.quantile(wildtype, 0.75)) if len(wildtype) else np.nan,
            "median_difference_mut_minus_wt": (
                float(np.median(mutant) - np.median(wildtype))
                if len(mutant) and len(wildtype)
                else np.nan
            ),
            "mann_whitney_u": np.nan,
            "p_value": np.nan,
            "cliffs_delta": np.nan,
            "cliffs_delta_ci_low": np.nan,
            "cliffs_delta_ci_high": np.nan,
            "tested": len(mutant) >= min_group and len(wildtype) >= min_group,
        }
        if row["tested"]:
            u_value, p_value, delta = mann_whitney_asymptotic(mutant, wildtype)
            stable_cancer_seed = int.from_bytes(cancer.encode("ascii"), "little") % (2**32)
            rng = np.random.default_rng(seed + stable_cancer_seed)
            ci_low, ci_high = bootstrap_delta_ci(mutant, wildtype, n_boot, rng)
            row.update(
                {
                    "mann_whitney_u": u_value,
                    "p_value": p_value,
                    "cliffs_delta": delta,
                    "cliffs_delta_ci_low": ci_low,
                    "cliffs_delta_ci_high": ci_high,
                }
            )
        rows.append(row)
    result = pd.DataFrame(rows)
    result["q_value_bh"] = bh_adjust(result["p_value"])
    result["significant_fdr_0_05"] = result["q_value_bh"].lt(0.05).fillna(False)
    result["direction_mut_vs_wt"] = np.select(
        [result["cliffs_delta"].lt(0), result["cliffs_delta"].gt(0)],
        ["lower", "higher"],
        default="not_tested_or_equal",
    )
    return result


def plot_distributions(scores: pd.DataFrame, stats: pd.DataFrame, path: Path, seed: int) -> None:
    cancers = sorted(scores["cancer_type"].unique())
    stat_lookup = stats.set_index("cancer_type")
    fig, axes = plt.subplots(4, 8, figsize=(24, 13), sharey=True)
    rng = np.random.default_rng(seed)
    score_min = float(scores["activity_z_within_cancer"].min())
    score_max = float(scores["activity_z_within_cancer"].max())
    lower = min(-3.0, score_min - 0.15)
    upper = max(3.0, score_max + 0.15)
    for axis, cancer in zip(axes.flat, cancers):
        cohort = scores.loc[scores["cancer_type"].eq(cancer)]
        arrays = [
            cohort.loc[cohort["tp53_status"].eq(status), "activity_z_within_cancer"].to_numpy()
            for status in ["WT", "MUT"]
        ]
        for position, (values, color) in enumerate(zip(arrays, [WT_COLOR, MUT_COLOR])):
            if len(values) >= 2:
                violin = axis.violinplot(
                    values,
                    positions=[position],
                    widths=0.72,
                    showmeans=False,
                    showmedians=False,
                    showextrema=False,
                )
                violin["bodies"][0].set_facecolor(color)
                violin["bodies"][0].set_edgecolor("none")
                violin["bodies"][0].set_alpha(0.34)
            if len(values):
                box = axis.boxplot(
                    values,
                    positions=[position],
                    widths=0.28,
                    patch_artist=True,
                    showfliers=False,
                    medianprops={"color": "white", "linewidth": 1.5},
                    boxprops={"facecolor": color, "edgecolor": color, "linewidth": 1},
                    whiskerprops={"color": color, "linewidth": 1},
                    capprops={"color": color, "linewidth": 1},
                )
                _ = box
                plotted = values if len(values) <= 180 else rng.choice(values, 180, replace=False)
                jitter = rng.uniform(-0.13, 0.13, size=len(plotted))
                axis.scatter(
                    np.full(len(plotted), position) + jitter,
                    plotted,
                    s=5,
                    color=color,
                    alpha=0.24,
                    linewidths=0,
                    rasterized=True,
                )
        q_value = stat_lookup.loc[cancer, "q_value_bh"]
        title_color = "#8B0000" if pd.notna(q_value) and q_value < 0.05 else "#222222"
        q_text = f"q={q_value:.2g}" if pd.notna(q_value) else "not tested"
        axis.set_title(f"{cancer}  {q_text}", fontsize=9.5, color=title_color, fontweight="bold")
        axis.set_xticks([0, 1], [f"WT\nn={len(arrays[0])}", f"MUT\nn={len(arrays[1])}"], fontsize=8)
        axis.axhline(0, color="#999999", linewidth=0.6, linestyle="--", zorder=0)
        axis.set_xlim(-0.55, 1.55)
        axis.set_ylim(lower, upper)
        axis.grid(axis="y", color="#E7E7E7", linewidth=0.6)
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.tick_params(axis="y", labelsize=8, length=0)
    for axis in axes.flat[len(cancers) :]:
        axis.axis("off")
    fig.suptitle(
        "PROGENy-derived p53/DDR footprint by TP53 mutation status across TCGA PanCancer Atlas",
        x=0.06,
        y=0.995,
        ha="left",
        fontsize=19,
        fontweight="bold",
    )
    fig.text(
        0.06,
        0.968,
        "One matched RNA sample per sequenced patient; footprint standardized within cancer. Red titles: BH-FDR q < 0.05.",
        ha="left",
        fontsize=11,
        color="#444444",
    )
    fig.text(0.012, 0.5, "p53/DDR footprint score (within-cancer z)", rotation=90, va="center", fontsize=12)
    fig.tight_layout(rect=(0.035, 0.025, 0.995, 0.945), h_pad=1.1, w_pad=0.8)
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_forest(stats: pd.DataFrame, path: Path) -> None:
    plotted = stats.loc[stats["tested"]].sort_values("cliffs_delta").reset_index(drop=True)
    y_position = np.arange(len(plotted))
    colors = np.where(
        plotted["significant_fdr_0_05"] & plotted["cliffs_delta"].lt(0),
        "#2166AC",
        np.where(
            plotted["significant_fdr_0_05"] & plotted["cliffs_delta"].gt(0),
            "#B2182B",
            "#777777",
        ),
    )
    fig, axis = plt.subplots(figsize=(10.5, max(8.5, 0.34 * len(plotted) + 2.5)))
    left_error = plotted["cliffs_delta"] - plotted["cliffs_delta_ci_low"]
    right_error = plotted["cliffs_delta_ci_high"] - plotted["cliffs_delta"]
    for row_index, row in plotted.iterrows():
        axis.errorbar(
            row["cliffs_delta"],
            row_index,
            xerr=np.array([[left_error.iloc[row_index]], [right_error.iloc[row_index]]]),
            fmt="o",
            color=colors[row_index],
            ecolor=colors[row_index],
            markersize=5.5 + 0.012 * math.sqrt(row["n_matched"]),
            capsize=2,
            linewidth=1.35,
            zorder=3,
        )
        q_text = f"q={row['q_value_bh']:.2g}"
        axis.text(
            1.03,
            row_index,
            f"{int(row['n_mut'])}/{int(row['n_wt'])}   {q_text}",
            transform=axis.get_yaxis_transform(),
            va="center",
            fontsize=8.5,
            clip_on=False,
        )
    axis.axvline(0, color="#333333", linewidth=1, linestyle="--")
    axis.set_yticks(y_position, plotted["cancer_type"], fontsize=9)
    axis.set_xlim(-1.03, 1.03)
    axis.set_xlabel("Cliff's delta: P(MUT > WT) − P(MUT < WT)", fontsize=11)
    axis.set_title(
        "Effect of TP53 mutation status on the PROGENy-derived p53/DDR footprint",
        loc="left",
        fontsize=16,
        fontweight="bold",
        pad=28,
    )
    axis.text(
        0,
        1.015,
        "Points show Cliff's delta; bars show 95% stratified-bootstrap CI. Tested when both groups n ≥ 10.",
        transform=axis.transAxes,
        fontsize=10,
        color="#444444",
        va="bottom",
    )
    axis.text(1.03, 1.015, "n MUT/WT   BH q", transform=axis.transAxes, ha="left", va="bottom", fontsize=9)
    axis.text(0.01, -0.08, "← lower activity in TP53-mutant", transform=axis.transAxes, color="#2166AC", fontsize=9)
    axis.text(0.99, -0.08, "higher activity in TP53-mutant →", transform=axis.transAxes, ha="right", color="#B2182B", fontsize=9)
    axis.grid(axis="x", color="#E5E5E5", linewidth=0.7)
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.tick_params(axis="y", length=0)
    fig.subplots_adjust(left=0.12, right=0.78, top=0.9, bottom=0.1)
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(scores: pd.DataFrame, stats: pd.DataFrame, path: Path) -> None:
    order = (
        stats.sort_values("cliffs_delta", na_position="last")["cancer_type"].tolist()
    )
    medians = scores.pivot_table(
        index="cancer_type",
        columns="tp53_status",
        values="activity_z_within_cancer",
        aggfunc="median",
    ).reindex(index=order, columns=["WT", "MUT"])
    counts = scores.pivot_table(
        index="cancer_type",
        columns="tp53_status",
        values="patient_id",
        aggfunc="count",
    ).reindex(index=order, columns=["WT", "MUT"])
    values = medians.to_numpy(dtype=float)
    limit = max(0.5, float(np.nanmax(np.abs(values))))
    fig, axis = plt.subplots(figsize=(7.2, 12))
    image = axis.imshow(
        values,
        aspect="auto",
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
    )
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            value = values[row, column]
            count = counts.iloc[row, column]
            if np.isfinite(value):
                text_color = "white" if abs(value) > 0.58 * limit else "#222222"
                axis.text(
                    column,
                    row,
                    f"{value:+.2f}\n(n={int(count)})",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color=text_color,
                )
    axis.set_xticks([0, 1], ["TP53 WT", "TP53 MUT"], fontsize=11, fontweight="bold")
    axis.set_yticks(np.arange(len(order)), order, fontsize=9)
    axis.tick_params(length=0)
    axis.spines[:].set_visible(False)
    axis.set_title(
        "Median PROGENy-derived p53/DDR footprint by cancer and TP53 status",
        loc="left",
        fontsize=15,
        fontweight="bold",
        pad=25,
    )
    axis.text(
        0,
        1.01,
        "Rows ordered by Cliff's delta (MUT vs WT); values are within-cancer z-score medians.",
        transform=axis.transAxes,
        fontsize=9.5,
        color="#444444",
        va="bottom",
    )
    colorbar = fig.colorbar(image, ax=axis, fraction=0.045, pad=0.06)
    colorbar.set_label("Median footprint z-score", fontsize=10)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def build_manifest(
    studies: pd.DataFrame,
    inputs: Path,
    out: Path,
    retrieved_at: str,
) -> pd.DataFrame:
    rows = [
        {
            "role": "TP53 mutation status",
            "source": "cBioPortal TCGA PanCancer Atlas 2018",
            "dataset_or_profile": "Per-study mutation profile + *_sequenced sample list",
            "identifier": "See study_id, mutation_profile_id and sample_list_id in inputs/study_ids.tsv",
            "normalization_or_processing": "Patient roster and audit label; WT is restricted to sequenced patients",
            "retrieved_at_utc": retrieved_at,
            "url": "https://www.cbioportal.org/",
            "local_file": "inputs/tp53_patient_status.tsv",
            "sha256": sha256_file(inputs / "tp53_patient_status.tsv"),
        },
        {
            "role": "Selected-sample TP53 mutation status",
            "source": "cBioPortal TCGA PanCancer Atlas 2018",
            "dataset_or_profile": "Retained protein-altering TP53 mutation events",
            "identifier": "Exact sample_id from per-study mutation profile",
            "normalization_or_processing": "MUT when the selected RNA sample has >=1 retained event; otherwise WT within the sequenced roster",
            "retrieved_at_utc": retrieved_at,
            "url": "https://www.cbioportal.org/",
            "local_file": "inputs/tp53_mutations.tsv",
            "sha256": sha256_file(inputs / "tp53_mutations.tsv"),
        },
        {
            "role": "Underlying PanCanAtlas expression matrix",
            "source": "NCI GDC PanCanAtlas publication supplement",
            "dataset_or_profile": GDC_EXPR_FILENAME,
            "identifier": f"GDC UUID {GDC_EXPR_UUID}; MD5 {GDC_EXPR_MD5}",
            "normalization_or_processing": "Upper-quartile-normalized RSEM; EB++ batch correction",
            "retrieved_at_utc": retrieved_at,
            "url": f"https://api.gdc.cancer.gov/data/{GDC_EXPR_UUID}",
            "local_file": "Not downloaded in full; cBioPortal profile slices cached under results/raw/cbioportal_expression",
            "sha256": "NA (official MD5 shown in identifier)",
        },
        {
            "role": "PROGENy p53 model",
            "source": "Bioconductor",
            "dataset_or_profile": MODEL_VERSION,
            "identifier": "model_human_full; pathway=p53; top 100 by ascending p.value",
            "normalization_or_processing": "Frozen signed weights; score = sum(expression z-score × weight)",
            "retrieved_at_utc": retrieved_at,
            "url": MODEL_URL,
            "local_file": "inputs/progeny_p53_top100_weights.tsv",
            "sha256": sha256_file(inputs / "progeny_p53_top100_weights.tsv"),
        },
    ]
    for study in studies.itertuples(index=False):
        profile_id = study.study_id + EXPR_SUFFIX
        rows.append(
            {
                "role": "Matched expression",
                "source": "cBioPortal TCGA PanCancer Atlas 2018",
                "dataset_or_profile": profile_id,
                "identifier": f"study={study.study_id}; sample_list={study.study_id + RNA_LIST_SUFFIX}",
                "normalization_or_processing": "mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM)",
                "retrieved_at_utc": retrieved_at,
                "url": f"{API}/molecular-profiles/{profile_id}",
                "local_file": f"results/raw/cbioportal_expression/{study.cancer_type}_{profile_id}.json.gz",
                "sha256": sha256_file(
                    out / "raw" / "cbioportal_expression" / f"{study.cancer_type}_{profile_id}.json.gz"
                ),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    inputs = args.inputs.resolve()
    out = args.out.resolve()
    table_dir = out / "tables"
    figure_dir = out / "figures"
    raw_dir = out / "raw"
    cache_dir = raw_dir / "cbioportal_expression"
    for path in [table_dir, figure_dir, cache_dir, out / ".mplconfig"]:
        path.mkdir(parents=True, exist_ok=True)

    retrieved_at = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    status = pd.read_csv(inputs / "tp53_patient_status.tsv", sep="\t", dtype=str).fillna("")
    mutation_events = pd.read_csv(inputs / "tp53_mutations.tsv", sep="\t", dtype=str).fillna("")
    studies = pd.read_csv(inputs / "study_ids.tsv", sep="\t", dtype=str)
    weights = pd.read_csv(inputs / "progeny_p53_top100_weights.tsv", sep="\t")
    if sha256_file(inputs / "progeny_p53_top100_weights.tsv") != MODEL_WEIGHTS_SHA256:
        raise ValueError("Frozen PROGENy p53 weight-file SHA-256 does not match")
    expected_ranks = np.arange(1, 101)
    valid_weights = (
        len(weights) == 100
        and weights["pathway"].eq("p53").all()
        and weights["gene"].is_unique
        and np.array_equal(weights["rank"].to_numpy(), expected_ranks)
        and weights["p.value"].is_monotonic_increasing
        and np.isfinite(weights[["weight", "p.value"]].to_numpy(dtype=float)).all()
    )
    if not valid_weights:
        raise ValueError("Frozen PROGENy weight file is not the expected 100-row p53 model")

    mutation_events["dbd_core_bool"] = mutation_events["dbd_core_102_292"].str.lower().eq("true")
    mutation_events["sample_event_category"] = np.select(
        [
            mutation_events["dbd_core_bool"] & mutation_events["variant_group"].eq("truncating"),
            mutation_events["dbd_core_bool"] & mutation_events["variant_group"].eq("missense"),
            mutation_events["variant_group"].eq("truncating"),
            mutation_events["variant_group"].eq("missense"),
        ],
        ["DBD-truncating", "DBD-missense", "non-DBD-truncating", "non-DBD-missense"],
        default="other",
    )
    category_priority = {
        "DBD-truncating": 1,
        "DBD-missense": 2,
        "non-DBD-truncating": 3,
        "non-DBD-missense": 4,
        "other": 5,
    }
    mutation_events["sample_category_priority"] = mutation_events["sample_event_category"].map(
        category_priority
    )
    mutation_events = mutation_events.sort_values(
        ["cancer_type", "sample_id", "sample_category_priority", "protein_position"],
        kind="mergesort",
    )
    sample_mutation = (
        mutation_events.groupby(["cancer_type", "patient_id", "sample_id"], as_index=False)
        .agg(
            assigned_category=("sample_event_category", "first"),
            n_retained_events=("sample_id", "size"),
            protein_changes=("protein_change", lambda values: ";".join(sorted(set(values) - {""}))),
            mutation_types=("mutation_type", lambda values: ";".join(sorted(set(values) - {""}))),
        )
    )
    sample_mutation["tp53_status"] = "MUT"

    genes = weights["gene"].tolist()
    mapping_records = fetch_gene_mapping(
        genes,
        raw_dir / "cbioportal_gene_mapping.json.gz",
        args.force_download,
        args.offline,
    )
    direct_symbol_to_entrez = {
        record["hugoGeneSymbol"]: int(record["entrezGeneId"]) for record in mapping_records
    }
    symbol_to_entrez = dict(direct_symbol_to_entrez)
    symbol_to_entrez.update(
        {symbol: entrez for symbol, entrez in LEGACY_SYMBOL_TO_ENTREZ.items() if symbol in genes}
    )
    weights["entrez_gene_id"] = weights["gene"].map(symbol_to_entrez).astype("Int64")
    weights["mapping_basis"] = np.select(
        [
            weights["gene"].isin(direct_symbol_to_entrez),
            weights["gene"].isin(LEGACY_SYMBOL_TO_ENTREZ),
        ],
        ["direct current HUGO symbol", "legacy symbol mapped by stable Entrez ID"],
        default="unmapped/suppressed locus",
    )
    weights["mapped_in_cbioportal"] = weights["entrez_gene_id"].notna()
    weights.to_csv(table_dir / "progeny_p53_top100_weights_with_entrez.tsv", sep="\t", index=False)
    mapped = weights.loc[weights["mapped_in_cbioportal"]].copy()
    entrez_to_gene = dict(zip(mapped["entrez_gene_id"].astype(int), mapped["gene"]))
    model_genes = mapped["gene"].tolist()
    model_weights = mapped.set_index("gene").loc[model_genes, "weight"].to_numpy(dtype=float)
    minimum_coverage = math.ceil(0.80 * len(model_genes))

    selected_frames: list[pd.DataFrame] = []
    wide_frames: list[pd.DataFrame] = []
    coverage_rows: list[dict] = []
    for study in studies.itertuples(index=False):
        cancer = study.cancer_type
        study_id = study.study_id
        cohort_status = status.loc[status["cancer_type"].eq(cancer)].copy()
        requested_samples = sorted(
            {sample for value in cohort_status["sample_ids"] for sample in split_sample_ids(value)}
        )
        profile_id = study_id + EXPR_SUFFIX
        cache_path = cache_dir / f"{cancer}_{profile_id}.json.gz"
        print(f"[{cancer}] requesting/caching {len(requested_samples)} samples × {len(mapped)} genes")
        records = fetch_study_expression(
            study_id,
            requested_samples,
            mapped["entrez_gene_id"].astype(int).tolist(),
            cache_path,
            args.force_download,
            args.offline,
        )
        if not records:
            print(f"  WARNING: no expression returned for {cancer}")
            continue
        long = pd.DataFrame(
            {
                "sample_id": [record["sampleId"] for record in records],
                "patient_id": [record["patientId"] for record in records],
                "entrez_gene_id": [int(record["entrezGeneId"]) for record in records],
                "value": [float(record["value"]) for record in records],
            }
        )
        long["gene"] = long["entrez_gene_id"].map(entrez_to_gene)
        long = long.dropna(subset=["gene"])
        wide = long.pivot_table(
            index=["sample_id", "patient_id"], columns="gene", values="value", aggfunc="first"
        ).reset_index()
        # API retrieval uses the optimized RNA sample list; analysis remains a
        # strict intersection with the mutation-sequenced sample IDs.
        wide = wide.loc[wide["sample_id"].isin(requested_samples)].copy()
        for gene in model_genes:
            if gene not in wide.columns:
                wide[gene] = np.nan
        wide = wide[["sample_id", "patient_id", *model_genes]]
        observed = wide[model_genes].notna().sum(axis=1)
        values = wide[model_genes].fillna(0.0).to_numpy(dtype=float)
        wide["n_signature_genes_observed"] = observed
        wide["activity_raw_weighted_sum"] = values @ model_weights
        wide["eligible_gene_coverage"] = observed.ge(minimum_coverage)
        wide["cancer_type"] = cancer
        wide_frames.append(wide.copy())

        score_frame = wide.loc[wide["eligible_gene_coverage"]].copy()
        score_frame = choose_one_sample_per_patient(score_frame)
        if score_frame.empty:
            continue
        raw_mean = score_frame["activity_raw_weighted_sum"].mean()
        raw_sd = score_frame["activity_raw_weighted_sum"].std(ddof=1)
        if not np.isfinite(raw_sd) or raw_sd <= 0:
            raise ValueError(f"Non-positive or non-finite raw activity SD in {cancer}: {raw_sd}")
        score_frame["activity_z_within_cancer"] = (
            score_frame["activity_raw_weighted_sum"] - raw_mean
        ) / raw_sd
        score_frame["activity_percentile_within_cancer"] = (
            score_frame["activity_raw_weighted_sum"].rank(method="average") - 0.5
        ) / len(score_frame)
        cohort_sample_mutation = sample_mutation.loc[
            sample_mutation["cancer_type"].eq(cancer),
            [
                "patient_id",
                "sample_id",
                "tp53_status",
                "assigned_category",
                "n_retained_events",
                "protein_changes",
                "mutation_types",
            ],
        ]
        score_frame = score_frame.merge(
            cohort_sample_mutation,
            on=["patient_id", "sample_id"],
            how="left",
            validate="one_to_one",
        )
        score_frame["tp53_status"] = score_frame["tp53_status"].fillna("WT")
        score_frame["assigned_category"] = score_frame["assigned_category"].fillna("WT")
        score_frame["n_retained_events"] = score_frame["n_retained_events"].fillna(0).astype(int)
        score_frame[["protein_changes", "mutation_types"]] = score_frame[
            ["protein_changes", "mutation_types"]
        ].fillna("")
        patient_status_columns = [
            "patient_id",
            "tp53_status",
            "assigned_category",
        ]
        patient_status = cohort_status[patient_status_columns].rename(
            columns={
                "tp53_status": "tp53_patient_status",
                "assigned_category": "patient_assigned_category",
            }
        )
        score_frame = score_frame.merge(
            patient_status, on="patient_id", how="inner", validate="one_to_one"
        )
        score_frame["selected_sample_vs_patient_status_concordant"] = score_frame[
            "tp53_status"
        ].eq(score_frame["tp53_patient_status"])
        score_frame["expression_profile_id"] = profile_id
        selected_frames.append(score_frame)
        coverage_rows.append(
            {
                "cancer_type": cancer,
                "study_id": study_id,
                "expression_profile_id": profile_id,
                "n_sequenced_patients": len(cohort_status),
                "n_requested_sample_ids": len(requested_samples),
                "n_samples_with_any_expression": wide["sample_id"].nunique(),
                "n_patients_scored": score_frame["patient_id"].nunique(),
                "matched_patient_fraction": score_frame["patient_id"].nunique() / len(cohort_status),
                "min_genes_observed": int(wide["n_signature_genes_observed"].min()),
                "median_genes_observed": float(wide["n_signature_genes_observed"].median()),
                "max_genes_observed": int(wide["n_signature_genes_observed"].max()),
                "minimum_genes_required": minimum_coverage,
                "n_samples_below_coverage_threshold": int((~wide["eligible_gene_coverage"]).sum()),
            }
        )

    if not selected_frames:
        raise RuntimeError("No matched expression samples were scored")
    scores = pd.concat(selected_frames, ignore_index=True)
    scores = scores.sort_values(["cancer_type", "patient_id"], kind="mergesort").reset_index(drop=True)
    if scores.duplicated(["cancer_type", "patient_id"]).any():
        raise AssertionError("Duplicate patient after sample selection")
    coverage = pd.DataFrame(coverage_rows).sort_values("cancer_type")
    wide_all = pd.concat(wide_frames, ignore_index=True)
    wide_all = wide_all.sort_values(["cancer_type", "patient_id", "sample_id"], kind="mergesort")

    score_columns = [
        "cancer_type",
        "patient_id",
        "sample_id",
        "sample_type_code",
        "tp53_status",
        "tp53_patient_status",
        "selected_sample_vs_patient_status_concordant",
        "assigned_category",
        "patient_assigned_category",
        "n_retained_events",
        "protein_changes",
        "mutation_types",
        "n_signature_genes_observed",
        "activity_raw_weighted_sum",
        "activity_z_within_cancer",
        "activity_percentile_within_cancer",
        "expression_profile_id",
    ]
    scores[score_columns].to_csv(table_dir / "patient_activity_scores.tsv", sep="\t", index=False)
    coverage.to_csv(table_dir / "cohort_expression_coverage.tsv", sep="\t", index=False)
    wide_all.to_csv(
        raw_dir / "progeny_p53_expression_zscores_wide.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    gene_coverage = weights.copy()
    if model_genes:
        nonmissing = scores[model_genes].notna()
        gene_coverage["n_selected_samples_with_value"] = gene_coverage["gene"].map(
            nonmissing.sum(axis=0).to_dict()
        ).fillna(0).astype(int)
        gene_coverage["n_cancers_with_any_value"] = gene_coverage["gene"].map(
            nonmissing.groupby(scores["cancer_type"]).any().sum(axis=0).to_dict()
        ).fillna(0).astype(int)
    gene_coverage.to_csv(table_dir / "signature_gene_coverage.tsv", sep="\t", index=False)

    stats = summarize_groups(scores, args.bootstrap, args.seed)
    stats.to_csv(table_dir / "cancer_mut_vs_wt_statistics.tsv", sep="\t", index=False)
    plot_distributions(scores, stats, figure_dir / "Figure1_activity_distributions", args.seed)
    plot_forest(stats, figure_dir / "Figure2_effect_size_forest")
    plot_heatmap(scores, stats, figure_dir / "Figure3_status_median_heatmap")

    manifest = build_manifest(studies, inputs, out, retrieved_at)
    manifest.to_csv(table_dir / "data_manifest.tsv", sep="\t", index=False)
    significant = stats.loc[stats["significant_fdr_0_05"]].sort_values("cliffs_delta")
    metadata = {
        "analysis_name": "Selected-sample TP53 mutation status versus PROGENy-derived p53/DDR transcriptional footprint",
        "created_at_utc": retrieved_at,
        "expression_profile_suffix": EXPR_SUFFIX,
        "expression_normalization": (
            "cBioPortal all-sample z-score profile derived from log RNA Seq V2 RSEM; "
            "underlying PanCanAtlas matrix is upper-quartile normalized and EB++ batch corrected"
        ),
        "underlying_expression_gdc_uuid": GDC_EXPR_UUID,
        "underlying_expression_gdc_md5": GDC_EXPR_MD5,
        "model": MODEL_VERSION,
        "model_url": MODEL_URL,
        "model_tarball_sha256": MODEL_TARBALL_SHA256,
        "signature_selection": "p53 pathway, top 100 genes by ascending model p.value",
        "signature_gene_count": int(len(weights)),
        "signature_gene_count_mapped": int(len(mapped)),
        "minimum_per_sample_gene_coverage": minimum_coverage,
        "score_formula": "raw=sum(z_expression_gene * signed_PROGENy_weight); then z-score raw values within cancer",
        "progeny_preprocessing_note": (
            "This is a transparent gene-z-score adaptation of the PROGENy weighted model, not the original-paper "
            "TCGA DESeq2-VST preprocessing; numerical scores are not directly interchangeable with that implementation"
        ),
        "missing_expression_handling": (
            "A missing gene z-score contributes zero (the cohort reference mean); samples below 80% gene coverage are excluded"
        ),
        "patient_sample_rule": (
            "Exact mutation-profile sample IDs only; one RNA sample per patient, prioritized 01,03,09,02,04,05,06,07,08; "
            "MUT/WT is assigned from the selected sample's retained TP53 events"
        ),
        "statistical_test": (
            "Two-sided asymptotic Mann-Whitney U with tie correction; Cliff's delta; "
            f"{args.bootstrap} stratified bootstrap replicates; BH correction across cohorts with n_MUT,n_WT>=10"
        ),
        "random_seed": args.seed,
        "n_cancers": int(scores["cancer_type"].nunique()),
        "n_matched_patients": int(len(scores)),
        "n_mut": int(scores["tp53_status"].eq("MUT").sum()),
        "n_wt": int(scores["tp53_status"].eq("WT").sum()),
        "n_selected_sample_patient_status_discordant": int(
            (~scores["selected_sample_vs_patient_status_concordant"]).sum()
        ),
        "n_tested_cancers": int(stats["tested"].sum()),
        "n_fdr_significant_cancers": int(stats["significant_fdr_0_05"].sum()),
        "fdr_significant_cancers": significant[
            ["cancer_type", "cliffs_delta", "q_value_bh", "direction_mut_vs_wt"]
        ].to_dict("records"),
        "important_interpretation": (
            "This is an association between mutation status and a transcriptional footprint, not a causal or clinical score. "
            "GTEx is excluded because it has no comparable somatic TP53 mutation calls."
        ),
    }
    with (out / "analysis_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2)

    checks = [
        ("PASS" if len(weights) == 100 else "FAIL", "frozen PROGENy p53 model has 100 rows"),
        (
            "PASS" if sha256_file(inputs / "progeny_1.34.0.tar.gz") == MODEL_TARBALL_SHA256 else "FAIL",
            "bundled progeny source tarball SHA-256 matches recorded value",
        ),
        (
            "PASS"
            if sha256_file(inputs / "progeny_p53_top100_weights.tsv") == MODEL_WEIGHTS_SHA256
            else "FAIL",
            "frozen top-100 weight TSV SHA-256 matches recorded value",
        ),
        ("PASS" if len(mapped) >= 95 else "WARN", f"{len(mapped)}/100 signature genes mapped to Entrez"),
        (
            "PASS" if scores["cancer_type"].nunique() == 32 else "FAIL",
            "all 32 requested TCGA cancer cohorts are represented",
        ),
        (
            "PASS" if set(scores["tp53_status"]) <= {"MUT", "WT"} else "FAIL",
            "selected-sample TP53 status contains only MUT/WT",
        ),
        (
            "PASS"
            if scores.loc[scores["tp53_status"].eq("MUT"), "n_retained_events"].ge(1).all()
            and scores.loc[scores["tp53_status"].eq("WT"), "n_retained_events"].eq(0).all()
            else "FAIL",
            "selected-sample MUT/WT agrees with retained mutation-event counts",
        ),
        (
            "PASS" if not scores.duplicated(["cancer_type", "patient_id"]).any() else "FAIL",
            "one activity row per cancer/patient",
        ),
        (
            "PASS" if scores["n_signature_genes_observed"].ge(minimum_coverage).all() else "FAIL",
            f"all reported samples meet >= {minimum_coverage} observed signature genes",
        ),
        (
            "PASS"
            if scores.groupby("cancer_type")["activity_z_within_cancer"].mean().abs().max() < 1e-10
            else "FAIL",
            "within-cancer activity means are zero after standardization",
        ),
        (
            "PASS"
            if np.allclose(
                scores.groupby("cancer_type")["activity_z_within_cancer"].std(ddof=1).to_numpy(),
                1.0,
                atol=1e-10,
            )
            else "FAIL",
            "within-cancer activity sample SDs are one",
        ),
        (
            "PASS"
            if all((figure_dir / name).exists() for name in [
                "Figure1_activity_distributions.png",
                "Figure1_activity_distributions.pdf",
                "Figure2_effect_size_forest.png",
                "Figure2_effect_size_forest.pdf",
                "Figure3_status_median_heatmap.png",
                "Figure3_status_median_heatmap.pdf",
            ])
            else "FAIL",
            "all PNG and PDF figures exist",
        ),
    ]
    report_lines = [
        "TP53 / PROGENy-derived p53-DDR footprint analysis validation",
        f"Created UTC: {retrieved_at}",
        "",
        *[f"[{status_code}] {message}" for status_code, message in checks],
        "",
        f"Matched patients: {len(scores):,}",
        f"Cancer types: {scores['cancer_type'].nunique()}",
        f"MUT / WT: {scores['tp53_status'].eq('MUT').sum():,} / {scores['tp53_status'].eq('WT').sum():,}",
        f"Selected-sample vs patient-level status discordant: {(~scores['selected_sample_vs_patient_status_concordant']).sum():,}",
        f"Tested cancers: {stats['tested'].sum()}",
        f"BH-FDR significant cancers: {stats['significant_fdr_0_05'].sum()}",
    ]
    (out / "validation_report.txt").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    if any(code == "FAIL" for code, _ in checks):
        raise AssertionError("One or more validation checks failed; inspect validation_report.txt")
    print("Analysis complete")
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
