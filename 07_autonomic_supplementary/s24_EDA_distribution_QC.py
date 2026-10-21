"""Audit exported EDA distributions and metric-specific QC, without reprocessing signals.

Run using the physio Python environment. No automatic ECG flag exclusion and no
claim of manual visual review. Historical cleaning affects SCR only, not SCL.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import neurokit2 as nk

BASE = Path(__file__).resolve().parents[1] / "neurokit_out"
OUT = BASE / "eda_distribution_QC_audit"
OUT.mkdir(exist_ok=True)
raw = pd.read_csv(BASE / "physio_training_long_neurokit.csv")
clean = pd.read_csv(BASE / "physio_training_long_neurokit_scr_cleaned.csv")
qc = pd.read_csv(BASE / "physio_qc_flags.csv")
if qc.duplicated(["Subject", "SegmentType", "Segment"]).any():
    raise ValueError("Duplicate QC segment identifiers")
if raw.duplicated(["Subject", "Block"]).any() or clean.duplicated(["Subject", "Block"]).any():
    raise ValueError("Duplicate analysis observation identifiers")
raw = raw.set_index(["Subject", "Block"]).sort_index()
clean = clean.set_index(["Subject", "Block"]).sort_index()
if not raw.index.equals(clean.index):
    raise ValueError("Raw and cleaned Subject/Block keys differ")

rows = []
for outcome, baseline in [("SCL_mean", "SCL_baseline"), ("SCR_rate", "SCR_baseline")]:
    d = clean[[outcome, baseline]].replace([np.inf, -np.inf], np.nan).dropna()
    for scope, x in [("all_finite_training", clean[outcome]), ("model_complete_cases", d[outcome])]:
        x = pd.to_numeric(x, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna().to_numpy()
        for scale in ["raw", "ln"]:
            if scale == "ln" and not (x > 0).all():
                continue
            y = np.log(x) if scale == "ln" else x
            rows.append(dict(outcome=outcome, scope=scope, scale=scale, n=len(y), skew=stats.skew(y, bias=False), shapiro_p=stats.shapiro(y).pvalue))
            fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
            ax[0].hist(y, bins="auto", edgecolor="black")
            ax[0].set_title(f"{outcome}: {scope}, {scale}")
            stats.probplot(y, plot=ax[1])
            fig.tight_layout()
            fig.savefig(OUT / f"{outcome}_{scope}_{scale}.png", dpi=150)
            plt.close(fig)
pd.DataFrame(rows).to_csv(OUT / "distribution_summary.csv", index=False)

audit = []
for metric, training_col, baseline_col, segment_col in [
    ("HR", "HR", "HR_baseline", "ECG_HR"),
    ("RMSSD", "RMSSD", "RMSSD_baseline", "ECG_RMSSD"),
    ("SCL", "SCL_mean", "SCL_baseline", "EDA_SCL_mean"),
    ("SCR", "SCR_rate", "SCR_baseline", "EDA_SCR_rate")]:
    rest = qc[qc.SegmentType == "rest"].copy()
    invalid_rest = ~np.isfinite(pd.to_numeric(rest[segment_col], errors="coerce"))
    if metric == "SCR":
        invalid_rest |= pd.to_numeric(rest[segment_col], errors="coerce") > 30
    # Finite original values made missing by existing cleaning, not all QC flags.
    direct = np.isfinite(raw[training_col]) & ~np.isfinite(clean[training_col])
    missing_baseline = ~np.isfinite(clean[baseline_col])
    model_missing = ~np.isfinite(clean[training_col]) | missing_baseline
    affected = set(rest.loc[invalid_rest, "Subject"]) | set(clean.index.get_level_values("Subject")[direct])
    audit.append(dict(metric=metric, invalid_rest_segments=int(invalid_rest.sum()), direct_training_exclusions=int(direct.sum()), directly_affected_participants=len(affected), missing_baseline_participants=len(set(clean.index.get_level_values("Subject")[missing_baseline])), model_excluded_training_observations=int(model_missing.sum()), model_observations=int((~model_missing).sum()), model_participants=len(set(clean.index.get_level_values("Subject")[~model_missing]))))
    rest.loc[invalid_rest].to_csv(OUT / f"{metric}_invalid_rest_segments.csv", index=False)
pd.DataFrame(audit).to_csv(OUT / "QC_missingness_summary.csv", index=False)
# Export proposed waveform-review inventory; blank decisions do not imply review.
flags = [c for c in qc if c.startswith("Flag_") and not c.endswith("_bool")]
flagged = qc[qc[flags].apply(lambda c: c.astype(str).str.lower().isin(["true", "1"])).any(axis=1)].copy()
flagged["Visual_review_decision"] = ""
flagged["Visual_review_notes"] = ""
review_path = OUT / "visual_review_inventory.csv"
if not review_path.exists():
    flagged.to_csv(review_path, index=False)
metadata = dict(neurokit_version_current_environment=nk.__version__, original_extraction_version="not recorded in historical CSV", threshold="candidate phasic peak prominence / maximum candidate peak prominence in the same segment >= 0.1", visual_review="not documented in historical cleaning log", no_model_imputation=True)
(OUT / "implementation_and_provenance.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
print(pd.DataFrame(rows).to_string(index=False))
print(pd.DataFrame(audit).to_string(index=False))
