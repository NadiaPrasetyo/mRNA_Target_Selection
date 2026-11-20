#!/usr/bin/env python3
"""
Random Forest antigen-prediction pipeline - single leave-one-clade-out model.

DESIGN
  1. Load the labelled training bacteria (antigen = 1 / non-antigen = 0).
  2. Using the NCBI taxonomy below, pick ONE clade at random (default rank: phylum;
     seeded, so it is reproducible) and hold ALL species of that clade out of training.
  3. Train ONE Random Forest (default hyper-parameters, no tuning) on the remaining
     clades.
  4. Evaluate that model on the held-out clade (confusion matrix, sensitivity, specificity,
     precision, F1, MCC, ROC-AUC, PR-AUC with 95% CIs).
  5. Only then use the SAME trained model to predict the bacterium of interest
     (default S.aureus). That bacterium is never used for training, tuning, feature
     selection or evaluation; any labels it has are ignored.

  The feature space is defined by the TRAINING species only; the held-out clade and the
  bacterium of interest are aligned to it (extra features dropped, missing filled with 0).
  No homology filtering is performed (see notes in the report).

All details are printed to the console and written to <output_dir>/model_report.txt.
"""

import os
import sys
import glob
import json
import argparse
import logging
import platform
from datetime import datetime, timezone
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import requests
import scipy
import sklearn
from scipy.stats import norm
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.preprocessing import LabelEncoder

# ----------------------
# NCBI taxonomy (from the provided NCBI tree). Each entry: (rank, name, tax_id)
# ----------------------
RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

_BACILLATI = ("kingdom", "Bacillati", 1783272)
_PSEUDOMONADATI = ("kingdom", "Pseudomonadati", 3379134)
_BACILLOTA = ("phylum", "Bacillota", 1239)
_PSEUDOMONADOTA = ("phylum", "Pseudomonadota", 1224)
_GAMMA = ("class", "Gammaproteobacteria", 1236)

LINEAGE: Dict[str, List[Tuple[str, str, int]]] = {
    "S.pneumoniae": [_BACILLATI, _BACILLOTA, ("class", "Bacilli", 91061),
                     ("order", "Lactobacillales", 186826), ("family", "Streptococcaceae", 1300),
                     ("genus", "Streptococcus", 1301), ("species", "Streptococcus pneumoniae", 1313)],
    "S.pyogenes": [_BACILLATI, _BACILLOTA, ("class", "Bacilli", 91061),
                   ("order", "Lactobacillales", 186826), ("family", "Streptococcaceae", 1300),
                   ("genus", "Streptococcus", 1301), ("species", "Streptococcus pyogenes", 1314)],
    "H.pylori": [_PSEUDOMONADATI, ("phylum", "Campylobacterota", 29547),
                 ("class", "Epsilonproteobacteria", 3031852), ("order", "Campylobacterales", 213849),
                 ("family", "Helicobacteraceae", 72293), ("genus", "Helicobacter", 209),
                 ("species", "Helicobacter pylori", 210)],
    "B.melitensis": [_PSEUDOMONADATI, _PSEUDOMONADOTA, ("class", "Alphaproteobacteria", 28211),
                     ("order", "Hyphomicrobiales", 356), ("family", "Brucellaceae", 118882),
                     ("genus", "Brucella", 234), ("species", "Brucella melitensis", 29459)],
    "N.gonorrhoeae": [_PSEUDOMONADATI, _PSEUDOMONADOTA, ("class", "Betaproteobacteria", 28216),
                      ("order", "Neisseriales", 206351), ("family", "Neisseriaceae", 481),
                      ("genus", "Neisseria", 482), ("species", "Neisseria gonorrhoeae", 485)],
    "C.burnetii": [_PSEUDOMONADATI, _PSEUDOMONADOTA, _GAMMA, ("order", "Legionellales", 118969),
                   ("family", "Coxiellaceae", 118968), ("genus", "Coxiella", 776),
                   ("species", "Coxiella burnetii", 777)],
    "H.influenzae": [_PSEUDOMONADATI, _PSEUDOMONADOTA, _GAMMA, ("order", "Pasteurellales", 135625),
                     ("family", "Pasteurellaceae", 712), ("genus", "Haemophilus", 724),
                     ("species", "Haemophilus influenzae", 727)],
    "P.aeruginosa": [_PSEUDOMONADATI, _PSEUDOMONADOTA, _GAMMA, ("order", "Pseudomonadales", 72274),
                     ("family", "Pseudomonadaceae", 135621), ("genus", "Pseudomonas", 286),
                     ("species", "Pseudomonas aeruginosa", 287)],
    "V.cholerae": [_PSEUDOMONADATI, _PSEUDOMONADOTA, _GAMMA, ("order", "Vibrionales", 135623),
                   ("family", "Vibrionaceae", 641), ("genus", "Vibrio", 662),
                   ("species", "Vibrio cholerae", 666)],
    "T.pallidum": [_PSEUDOMONADATI, ("phylum", "Spirochaetota", 203691), ("class", "Spirochaetia", 203692),
                   ("order", "Treponematales", 3120456), ("family", "Treponemataceae", 2845253),
                   ("genus", "Treponema", 157), ("species", "Treponema pallidum", 160)],
    "C.trachomatis": [_PSEUDOMONADATI, ("phylum", "Chlamydiota", 204428), ("class", "Chlamydiia", 204429),
                      ("order", "Chlamydiales", 51291), ("family", "Chlamydiaceae", 809),
                      ("genus", "Chlamydia", 810), ("species", "Chlamydia trachomatis", 813)],
    # --- NOT in the supplied tree: added from NCBI Taxonomy knowledge - please verify IDs ---
    "S.enteritidis": [_PSEUDOMONADATI, _PSEUDOMONADOTA, _GAMMA, ("order", "Enterobacterales", 91347),
                      ("family", "Enterobacteriaceae", 543), ("genus", "Salmonella", 590),
                      ("species", "Salmonella enterica subsp. enterica serovar Enteritidis", 149539)],
    "S.aureus": [_BACILLATI, _BACILLOTA, ("class", "Bacilli", 91061), ("order", "Bacillales", 1385),
                 ("family", "Staphylococcaceae", 90964), ("genus", "Staphylococcus", 1279),
                 ("species", "Staphylococcus aureus", 1280)],
}

DEFAULT_TRAIN_BACTERIA = [
    "B.melitensis", "C.burnetii", "C.trachomatis", "H.influenzae", "H.pylori",
    "N.gonorrhoeae", "P.aeruginosa", "S.enteritidis", "S.pneumoniae", "S.pyogenes",
    "T.pallidum", "V.cholerae",
]

METRIC_ORDER = ["sensitivity_recall", "specificity", "precision", "f1", "mcc",
                "accuracy", "roc_auc", "pr_auc"]


def clade_at(species: str, rank: str) -> str:
    for r, name, _ in LINEAGE[species]:
        if r == rank:
            return name
    raise KeyError(f"{species} has no {rank} in LINEAGE")


def lineage_table(species_list: List[str]) -> pd.DataFrame:
    rows = []
    for sp in species_list:
        row = {"bacterium": sp}
        for r, name, tid in LINEAGE[sp]:
            if r != "species":
                row[r] = f"{name} ({tid})"
        rows.append(row)
    return pd.DataFrame(rows)


# ----------------------
# Output helper: print to console AND keep for model_report.txt
# ----------------------
REPORT: List[str] = []


def emit(text: str = "") -> None:
    print(text, flush=True)
    REPORT.append(text)


def section(title: str) -> None:
    emit("")
    emit("=" * 78)
    emit(title)
    emit("=" * 78)


# ----------------------
# Metrics with confidence intervals
# ----------------------
def wilson_ci(k: int, n: int, alpha: float = 0.05) -> Tuple[float, float]:
    if n == 0:
        return (np.nan, np.nan)
    z = norm.ppf(1 - alpha / 2)
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def _safe_div(a: float, b: float) -> float:
    return a / b if b > 0 else np.nan


def confusion_counts(y_true, y_pred) -> Tuple[int, int, int, int]:
    y_true = np.asarray(y_true).astype(int)
    y_pred = np.asarray(y_pred).astype(int)
    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())
    return tn, fp, fn, tp


def _point_metrics(y_true, probs, threshold) -> Dict[str, float]:
    tn, fp, fn, tp = confusion_counts(y_true, (probs >= threshold).astype(int))
    mcc_denom = np.sqrt(float(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if len(np.unique(y_true)) < 2:
        roc = pr = np.nan
    else:
        roc, pr = roc_auc_score(y_true, probs), average_precision_score(y_true, probs)
    return {"sensitivity_recall": _safe_div(tp, tp + fn), "specificity": _safe_div(tn, tn + fp),
            "precision": _safe_div(tp, tp + fp), "f1": _safe_div(2 * tp, 2 * tp + fp + fn),
            "mcc": (tp * tn - fp * fn) / mcc_denom if mcc_denom > 0 else 0.0,
            "accuracy": _safe_div(tp + tn, tp + tn + fp + fn), "roc_auc": roc, "pr_auc": pr}


def flat_metrics(y_true, probs, threshold, n_boot, seed, alpha=0.05) -> Dict[str, float]:
    """
    Counts, confusion matrix, metrics and 95% CIs in one flat dict.
    Wilson CIs: sensitivity, specificity, precision, accuracy.
    Stratified bootstrap percentile CIs: F1, MCC, ROC-AUC, PR-AUC.
    """
    y_true = np.asarray(y_true).astype(int)
    probs = np.asarray(probs, dtype=float)
    est = _point_metrics(y_true, probs, threshold)
    tn, fp, fn, tp = confusion_counts(y_true, (probs >= threshold).astype(int))

    ci = {"sensitivity_recall": wilson_ci(tp, tp + fn, alpha), "specificity": wilson_ci(tn, tn + fp, alpha),
          "precision": wilson_ci(tp, tp + fp, alpha), "accuracy": wilson_ci(tp + tn, len(y_true), alpha)}
    boot_keys = ["f1", "mcc", "roc_auc", "pr_auc"]
    pos, neg = np.where(y_true == 1)[0], np.where(y_true == 0)[0]
    if len(pos) and len(neg) and n_boot > 0:
        rng = np.random.default_rng(seed)
        boots = {k: [] for k in boot_keys}
        for _ in range(n_boot):
            idx = np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))])
            m = _point_metrics(y_true[idx], probs[idx], threshold)
            for k in boot_keys:
                boots[k].append(m[k])
        for k in boot_keys:
            a = np.asarray(boots[k], dtype=float)
            ci[k] = ((np.nanpercentile(a, 100 * alpha / 2), np.nanpercentile(a, 100 * (1 - alpha / 2)))
                     if not np.isnan(a).all() else (np.nan, np.nan))
    else:
        ci.update({k: (np.nan, np.nan) for k in boot_keys})

    n = len(y_true)
    out = {"n": n, "n_pos": int(y_true.sum()), "n_neg": n - int(y_true.sum()),
           "prevalence": y_true.mean() if n else np.nan, "threshold": threshold,
           "tn": tn, "fp": fp, "fn": fn, "tp": tp}
    for m in METRIC_ORDER:
        out[m], out[f"{m}_ci_lo"], out[f"{m}_ci_hi"] = est[m], ci[m][0], ci[m][1]
    return out


def print_metrics(title: str, r: Dict[str, float]) -> None:
    emit(f"{title}")
    emit(f"  n={r['n']} (pos={r['n_pos']}, neg={r['n_neg']}, prevalence={r['prevalence']:.3f}); "
         f"threshold={r['threshold']}")
    emit("  Confusion matrix (rows = true, cols = predicted):")
    emit(f"                    pred_0   pred_1")
    emit(f"    true_0 (non)    {r['tn']:>6}   {r['fp']:>6}")
    emit(f"    true_1 (antigen){r['fn']:>6}   {r['tp']:>6}")
    for m in METRIC_ORDER:
        emit(f"  {m:<20s} {r[m]:.3f}  (95% CI {r[m + '_ci_lo']:.3f}-{r[m + '_ci_hi']:.3f})")
    emit(f"  PR-AUC no-skill baseline (= prevalence): {r['prevalence']:.3f}")


# ----------------------
# Data loading
# ----------------------
def load_bacterium_data(base_dir: str, bacterium: str) -> pd.DataFrame:
    folder = os.path.join(base_dir, bacterium, "raw_data")
    files = glob.glob(os.path.join(folder, "*_raw_data.csv"))
    if not files:
        logging.warning(f"No files found for {bacterium} at {folder}")
        return pd.DataFrame()
    dfs = []
    for f in files:
        try:
            dfs.append(pd.read_csv(f))
        except Exception as e:
            logging.error(f"Failed to read {f}: {e}")
    if not dfs:
        return pd.DataFrame()
    out = pd.concat(dfs, ignore_index=True)
    out.columns = [c.strip() for c in out.columns]
    return out


def preprocess_and_pivot(df: pd.DataFrame):
    """accession x feature_subfeature matrix (mean over replicates, missing -> 0) + labels (or None)."""
    missing = {"accession", "feature", "subfeature", "value"} - set(df.columns)
    if missing:
        raise ValueError(f"Input dataframe missing required cols: {missing}")
    df = df.copy()
    df = df[df["accession"].astype(str).str.len() > 4]
    if df.empty:
        return pd.DataFrame(), None
    df["feature_subfeature"] = df["feature"].astype(str) + "_" + df["subfeature"].astype(str)
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df[df["value"].notna()]
    agg = df.groupby(["accession", "feature_subfeature"], observed=True)["value"].mean().reset_index()
    X = agg.pivot(index="accession", columns="feature_subfeature", values="value").fillna(0.0)
    labels = None
    if "label" in df.columns:
        labels = (df.groupby("accession", observed=True)["label"]
                    .agg(lambda x: x.mode().iat[0] if not x.mode().empty else x.iloc[0])).reindex(X.index)
    return X, labels


def make_binary_labels(labels: pd.Series) -> np.ndarray:
    """0/1 labels. Numeric {0,1}; or strings containing antigen/positive/... (not starting with
    non/not/neg/no/false) -> 1; or two classes (second sorted -> 1, with warning)."""
    labs = labels.astype(str)
    unique = sorted(set(labs))
    try:
        if set(pd.to_numeric(labels.dropna().unique()).tolist()).issubset({0, 1}):
            return labels.astype(int).to_numpy()
    except Exception:
        pass
    pos_tokens = ("antigen", "positive", "pos", "1", "true", "yes")
    neg_prefixes = ("non", "not", "neg", "no", "false")
    pos_vals = [u for u in unique if any(t in u.lower() for t in pos_tokens)
                and not u.lower().startswith(neg_prefixes)]
    if pos_vals and len(pos_vals) < len(unique):
        return labs.map({u: int(u in pos_vals) for u in unique}).astype(int).to_numpy()
    if len(unique) == 2:
        logging.warning(f"Ambiguous labels {unique}: mapping '{unique[1]}' -> 1.")
        return LabelEncoder().fit(unique).transform(labs)
    raise ValueError(f"Unable to convert labels to binary. Found classes: {unique}")


def load_labelled_species(base_dir: str, species: str):
    """Returns (X_df indexed by accession, y) for a labelled species, or None."""
    df = load_bacterium_data(base_dir, species)
    if df.empty:
        return None
    X, labels = preprocess_and_pivot(df)
    if X.empty or labels is None:
        logging.error(f"{species}: no usable features/labels; skipping.")
        return None
    keep = labels.notna()
    if (~keep).any():
        logging.warning(f"{species}: dropping {int((~keep).sum())} accessions with missing labels.")
    X, labels = X.loc[keep], labels.loc[keep]
    return X, make_binary_labels(labels)


# ----------------------
# UniProt names (optional)
# ----------------------
def fetch_names(accessions) -> pd.DataFrame:
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    rows = {}
    for acc in accessions:
        try:
            r = requests.get(base_url, headers={"accept": "application/json"}, timeout=10,
                             params={"query": f"accession:{acc}", "fields": "accession,protein_name,gene_names",
                                     "format": "json", "size": 50})
            r.raise_for_status()
            pn, gn = set(), set()
            for e in r.json().get("results", []):
                d = e.get("proteinDescription", {})
                if "recommendedName" in d:
                    v = d["recommendedName"].get("fullName", {}).get("value")
                    if v:
                        pn.add(v)
                for key in ("submissionNames", "alternativeNames"):
                    for s in d.get(key, []):
                        v = s.get("fullName", {}).get("value")
                        if v:
                            pn.add(v)
                for g in e.get("genes", []):
                    for v in [g.get("geneName", {}).get("value")] + \
                             [s.get("value") for s in g.get("synonyms", [])] + \
                             [o.get("value") for o in g.get("orfNames", [])]:
                        if v:
                            gn.add(v)
            rows[acc] = {"protein_names": ", ".join(sorted(pn)) or None, "gene_names": ", ".join(sorted(gn)) or None}
        except Exception as e:
            logging.warning(f"UniProt lookup failed for {acc}: {e}")
    return pd.DataFrame.from_dict(rows, orient="index").rename_axis("accession")


# ----------------------
# Main
# ----------------------
def main(args: argparse.Namespace) -> None:
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.WARNING,
                        format="%(asctime)s - %(levelname)s - %(message)s")
    os.makedirs(args.output_dir, exist_ok=True)
    target = args.target_bacterium
    if target in args.train_bacteria:
        raise ValueError("The bacterium of interest must not be among --train-bacteria.")
    for sp in args.train_bacteria:  # the bacterium of interest needs no lineage
        if sp not in LINEAGE:
            raise ValueError(f"{sp} is missing from LINEAGE; add it to the taxonomy table.")

    section("RANDOM FOREST ANTIGEN PREDICTION - SINGLE LEAVE-ONE-CLADE-OUT MODEL")
    emit(f"Run time (UTC): {datetime.now(timezone.utc).isoformat(timespec='seconds')}")
    emit(f"Command: {' '.join(sys.argv)}")
    emit(f"Versions: python {platform.python_version()}, numpy {np.__version__}, pandas {pd.__version__}, "
         f"scikit-learn {sklearn.__version__}, scipy {scipy.__version__}")

    section("1. NCBI TAXONOMY OF THE BACTERIA (name (tax id))")
    tbl = lineage_table(args.train_bacteria)
    emit(tbl.to_string(index=False))
    tbl.to_csv(os.path.join(args.output_dir, "taxonomy_lineages.csv"), index=False)

    # ---- load all labelled training species ----
    data = {}
    for sp in args.train_bacteria:
        res = load_labelled_species(args.base_dir, sp)
        if res is not None:
            data[sp] = res
    if len(data) < 2:
        raise RuntimeError("Need at least two labelled species.")

    # ---- randomly choose the clade to leave out ----
    rank = args.holdout_rank
    clades = sorted({clade_at(sp, rank) for sp in data})
    if args.holdout_clade:
        if args.holdout_clade not in clades:
            raise ValueError(f"--holdout-clade {args.holdout_clade} not among {clades}")
        held_clade, how = args.holdout_clade, "specified with --holdout-clade"
    else:
        held_clade = str(np.random.default_rng(args.seed).choice(clades))
        how = f"chosen at random (numpy default_rng, seed={args.seed})"
    test_species = [sp for sp in data if clade_at(sp, rank) == held_clade]
    train_species = [sp for sp in data if sp not in test_species]

    section(f"2. TRAIN / TEST SPLIT (rank = {rank})")
    emit(f"Clades available at rank '{rank}': {', '.join(clades)}")
    emit(f"Held-out clade: {held_clade} - {how}")
    emit(f"TRAIN species ({len(train_species)}): {', '.join(train_species)}")
    emit(f"TEST (held-out clade) species ({len(test_species)}): {', '.join(test_species)}")
    emit(f"Bacterium of interest: {target} - predict-only, NOT used for training/tuning/evaluation")
    if target in LINEAGE and clade_at(target, rank) == held_clade:  # optional note only
        emit(f"  NOTE: {target} belongs to the held-out clade, so the model has seen NO {held_clade} "
             "proteins when predicting it.")
    emit("")
    emit(f"{'species':<16}{'role':<8}{'clade':<20}{'n':>6}{'pos':>6}{'neg':>6}")
    for sp, (X_sp, y_sp) in data.items():
        emit(f"{sp:<16}{('TEST' if sp in test_species else 'TRAIN'):<8}{clade_at(sp, rank):<20}"
             f"{len(y_sp):>6}{int(y_sp.sum()):>6}{int((y_sp == 0).sum()):>6}")

    # ---- build matrices; feature space = TRAIN species only ----
    features = sorted(set().union(*(data[sp][0].columns for sp in train_species)))

    def stack(species_list):
        Xs, ys, sp_col, acc = [], [], [], []
        for sp in species_list:
            X_sp, y_sp = data[sp]
            Xs.append(X_sp.reindex(columns=features, fill_value=0.0))
            ys.append(y_sp)
            sp_col += [sp] * len(y_sp)
            acc += list(X_sp.index.astype(str))
        return pd.concat(Xs, axis=0), np.concatenate(ys), np.array(sp_col), np.array(acc)

    X_tr, y_tr, _, _ = stack(train_species)
    X_te, y_te, sp_te, acc_te = stack(test_species)
    n_dropped = len(set().union(*(data[sp][0].columns for sp in test_species)) - set(features))
    emit("")
    emit(f"TRAIN: {X_tr.shape[0]} proteins (pos={int(y_tr.sum())}, neg={int((y_tr == 0).sum())}), "
         f"{len(features)} features (union over training species).")
    emit(f"TEST : {X_te.shape[0]} proteins (pos={int(y_te.sum())}, neg={int((y_te == 0).sum())}); "
         f"{n_dropped} test-only features dropped.")
    if len(np.unique(y_tr)) < 2 or len(np.unique(y_te)) < 2:
        raise RuntimeError("Train or test set has a single class - check labels / chosen clade.")

    # ---- train ----
    section("3. RANDOM FOREST (default settings, no hyper-parameter tuning)")
    rf = RandomForestClassifier(n_estimators=args.n_estimators, class_weight="balanced",
                                random_state=args.seed, n_jobs=-1)
    rf.fit(X_tr.to_numpy(dtype=np.float32), y_tr)
    params = rf.get_params()
    emit("Hyper-parameters (everything not listed in the script is the scikit-learn default):")
    for k in sorted(params):
        emit(f"  {k:<26s} {params[k]}")
    emit(f"Classification threshold: {args.threshold} (fixed a priori, not tuned)")
    emit("Hyper-parameter tuning: none. Nothing was tuned on the held-out clade or the bacterium of interest.")

    # ---- evaluate on held-out clade ----
    section(f"4. EVALUATION ON HELD-OUT CLADE: {held_clade}")
    p_te = rf.predict_proba(X_te.to_numpy(dtype=np.float32))[:, 1]
    overall = flat_metrics(y_te, p_te, args.threshold, args.n_boot, args.seed)
    print_metrics(f"Held-out clade {held_clade} ({', '.join(test_species)}):", overall)
    rows = [{"group": held_clade, **overall}]
    if len(test_species) > 1:
        for sp in test_species:
            m = sp_te == sp
            r = flat_metrics(y_te[m], p_te[m], args.threshold, args.n_boot, args.seed)
            emit("")
            print_metrics(f"  Per species - {sp}:", r)
            rows.append({"group": sp, **r})
    pd.DataFrame(rows).to_csv(os.path.join(args.output_dir, "heldout_metrics.csv"), index=False)
    pd.DataFrame({"accession": acc_te, "species": sp_te, "true_label": y_te, "prob_antigen": p_te,
                  "pred_label": (p_te >= args.threshold).astype(int)}
                 ).to_csv(os.path.join(args.output_dir, "heldout_predictions.csv"), index=False)

    fi = pd.DataFrame({"feature": features, "importance": rf.feature_importances_})
    fi = fi.sort_values("importance", ascending=False).reset_index(drop=True)
    fi["rank"] = fi.index + 1
    fi.to_csv(os.path.join(args.output_dir, "feature_importances.csv"), index=False)
    emit("")
    emit("Top 10 features by impurity importance (full list in feature_importances.csv):")
    emit(fi.head(10).round(5).to_string(index=False))

    # ---- predict the bacterium of interest with the SAME model ----
    section(f"5. PREDICTION FOR BACTERIUM OF INTEREST: {target}")
    df_t = load_bacterium_data(args.base_dir, target)
    if df_t.empty:
        emit(f"No data found for {target}; no predictions made.")
    else:
        X_t, _ignored = preprocess_and_pivot(df_t)  # any label column is ignored
        only_t = len(set(X_t.columns) - set(features))
        absent = len(set(features) - set(X_t.columns))
        X_t = X_t.reindex(columns=features, fill_value=0.0)
        p_t = rf.predict_proba(X_t.to_numpy(dtype=np.float32))[:, 1]
        out = pd.DataFrame({"accession": X_t.index, "prob_antigen": p_t,
                            "pred_label": (p_t >= args.threshold).astype(int)}).set_index("accession")
        if not args.skip_annotation:
            out = out.merge(fetch_names(out.index), left_index=True, right_index=True, how="left")
        out = out.sort_values("prob_antigen", ascending=False)
        out.insert(0, "rank", np.arange(1, len(out) + 1))
        tag = target.replace(".", "_")
        out.to_csv(os.path.join(args.output_dir, f"{tag}_predictions.csv"))
        X_t.assign(prob_antigen=p_t).sort_values("prob_antigen", ascending=False).to_csv(
            os.path.join(args.output_dir, f"{tag}_features_with_probs.csv"))
        emit(f"{len(out)} proteins scored with the model trained on {len(train_species)} species "
             f"(held-out clade: {held_clade}).")
        emit(f"Feature alignment: {only_t} {target}-only features dropped; {absent} training features "
             "absent in the data (filled with 0).")
        emit(f"{int(out['pred_label'].sum())} proteins at/above threshold {args.threshold} "
             "(descriptive only - NOT a performance estimate; no labels were used).")
        emit(f"Probability quantiles (min/25/50/75/max): "
             f"{', '.join(f'{q:.3f}' for q in np.quantile(p_t, [0, .25, .5, .75, 1]))}")

    section("NOTES")
    emit("* No homology filtering was done: related proteins in the training clades may resemble "
         "held-out proteins, so the held-out estimate can still be optimistic.")
    emit("* A single random held-out clade gives ONE estimate; its value depends on which clade was drawn "
         "(change --seed or use --holdout-clade).")
    emit("* CIs: Wilson for sensitivity/specificity/precision/accuracy; stratified bootstrap "
         f"(n={args.n_boot}) for F1/MCC/ROC-AUC/PR-AUC. They resample proteins and ignore species/family structure.")

    with open(os.path.join(args.output_dir, "model_report.txt"), "w") as fh:
        fh.write("\n".join(REPORT) + "\n")
    with open(os.path.join(args.output_dir, "run_config.json"), "w") as fh:
        json.dump({"seed": args.seed, "holdout_rank": rank, "holdout_clade": held_clade,
                   "train_species": train_species, "test_species": test_species,
                   "bacterium_of_interest": target, "threshold": args.threshold,
                   "rf_params": {k: v for k, v in params.items()},
                   "n_features": len(features), "n_bootstrap": args.n_boot,
                   "command_line": " ".join(sys.argv)}, fh, indent=2, default=str)
    print(f"\nSaved outputs to {args.output_dir}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Single leave-one-clade-out Random Forest antigen pipeline")
    p.add_argument("--base-dir", default="./results", help="Folder containing <bacterium>/raw_data/*_raw_data.csv")
    p.add_argument("--output-dir", default="./results")
    p.add_argument("--train-bacteria", nargs="+", default=DEFAULT_TRAIN_BACTERIA)
    p.add_argument("--target-bacterium", default="S.aureus", help="Bacterium of interest (predict-only)")
    p.add_argument("--holdout-rank", default="phylum", choices=RANKS[:-1],
                   help="Taxonomic rank at which a clade is left out (default: phylum)")
    p.add_argument("--holdout-clade", default=None, help="Force a specific clade instead of a random one")
    p.add_argument("--seed", type=int, default=42, help="Seed for the random clade choice, RF and bootstrap")
    p.add_argument("--n-estimators", type=int, default=500)
    p.add_argument("--threshold", type=float, default=0.5)
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--skip-annotation", action="store_true", help="Skip UniProt name lookup")
    p.add_argument("--verbose", action="store_true")
    main(p.parse_args())