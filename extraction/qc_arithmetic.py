"""Arithmetic QC of the extraction table: flags only, nothing transcribed is changed.

Run from the project root:   python3 extraction/qc_arithmetic.py [--workbook path.xlsx]

Reads    IMD_genetics_analysis_bundle (1).zip    extraction_variants_round1.csv and extraction_studies_round1.csv (never edited)
         data/extraction_corrections.csv          rows already corrected by hand (marked, not re-applied)
         data/panel_candidates.csv                 the variants with three or more meta-analysable studies (the panel candidates)
         data/panel_allele_frequencies.csv         gnomAD v4 non-Finnish European allele frequencies for those variants
Writes   extraction/output/qc_rows.csv             one row per variant row: keys plus qc_ columns
         extraction/output/qc_studies.csv          one row per study: counts of flags and a priority (P1 to P3)
         extraction/output/qc_summary.csv          coverage: how many rows each check could see
         a reviewer workbook (xlsx, gitignored) when --workbook is given

What the checks are (see extraction/README.md):
  - orientation: which of six contrast forms (recessive, dominant, allelic, each for the effect and for the other allele) of a
    row's own six genotype cells reproduces its printed odds ratio within a factor of 1.25. "other_allele" is not an error.
  - interval: the estimate inside its own confidence interval, the interval not inverted, the log-scale arms roughly symmetric.
  - Hardy-Weinberg in controls (chi-square, 1 df; flagged at p < 0.001).
  - control allele frequency against gnomAD NFE, flagged only when exchanging the two homozygote cells would bring it markedly
    closer (a sign of swapped columns). Panel variants only, and only meaningful for European-ancestry cohorts.

A row that passes has passed a test against itself, not a comparison with its paper.
"""
import argparse, csv, io, json, math, os, sys, zipfile
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_kernel import (qc_num, qc_reconcile_row, qc_hwe_p, qc_ci_check, qc_control_eaf, qc_frequency_check,
                       qc_assert_transcription_intact, qc_write_workbook)

BUNDLE = "IMD_genetics_analysis_bundle (1).zip"
CELLS = ["n_cases_hom_effect", "n_cases_het", "n_cases_hom_other", "n_ctrl_hom_effect", "n_ctrl_het", "n_ctrl_hom_other"]
TRANSCRIBED = CELLS + ["study_id", "gene_symbol", "rsid", "variant_name_as_reported", "effect_allele", "other_allele", "genetic_model",
                       "comparison_as_reported", "subgroup", "effect_measure", "effect_estimate", "ci_lower", "ci_upper", "p_value",
                       "adjusted", "direction", "maf_cases", "maf_controls", "hwe_p_controls", "page_or_table", "notes"]
HWE_ALPHA = 1e-3


def read_bundle(name):
    with zipfile.ZipFile(BUNDLE) as z:
        with z.open("IMD_genetics_bundle/data/" + name) as f:
            return list(csv.DictReader(io.TextIOWrapper(f, encoding="utf-8")))


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def cells_of(row):
    v = [qc_num(row.get(c)) for c in CELLS]
    return None if any(x is None for x in v) else v


def panel_variants():
    return {r["rsid"]: r["gene"] for r in read_csv("data/panel_candidates.csv") if int(r["studies_inclusive"] or 0) >= 3}


def reference_freqs():
    out = {}
    for r in read_csv("data/panel_allele_frequencies.csv"):
        try:
            out[r["rsid"]] = json.loads(r["gnomad_nfe_freqs"])
        except (ValueError, TypeError):
            pass
    return out


def check_row(r, panel, ref, corrected):
    q = {"variant_row_id": r["variant_row_id"], "study_id": r["study_id"], "rsid": r["rsid"], "gene_symbol": r["gene_symbol"],
         "usable_tier": r.get("usable_tier", ""), "extraction_method": r.get("extraction_method", ""),
         "page_or_table": r.get("page_or_table", ""), "qc_panel_variant": "yes" if r["rsid"] in panel else "no",
         "qc_already_corrected": ";".join(sorted(corrected.get(r["variant_row_id"], [])))}
    cells = cells_of(r)
    est = qc_num(r.get("effect_estimate"))
    measure = r.get("effect_measure", "")
    adjusted = r.get("adjusted") == "yes"

    # Orientation: odds ratios only, against the row's own cells
    if measure == "OR" and cells and est:
        rec = qc_reconcile_row(cells, est, is_adjusted=adjusted)
        q["qc_orientation"] = rec["orientation"]
        q["qc_reproducing_contrast"] = rec["reproducing_contrast"]
        q["qc_all_matching_contrasts"] = ";".join(rec["matches"])
        q["qc_closest_recomputed_or"] = "" if not rec["ratios"] else "%.3g" % min(
            (v for v in rec["ratios"].values() if v), key=lambda v: abs(math.log(v / est)))
    else:
        q["qc_orientation"] = "not_checkable"
        q["qc_reproducing_contrast"] = q["qc_all_matching_contrasts"] = q["qc_closest_recomputed_or"] = ""
        q["qc_orientation_reason"] = ("no six-cell counts" if not cells else "no estimate" if not est else
                                      "measure is %s, not an odds ratio" % (measure or "missing"))

    # Interval
    ci = qc_ci_check(r.get("effect_estimate"), r.get("ci_lower"), r.get("ci_upper"))
    q["qc_interval_flag"] = ci["flag"] or ""
    q["qc_interval_detail"] = ci["detail"]

    # Hardy-Weinberg in controls
    hwe = qc_hwe_p(*(cells[3:] if cells else [None] * 3))
    q["qc_hwe_p_controls"] = "" if hwe is None else "%.3g" % hwe
    q["qc_hwe_flag"] = "hwe_violation_controls" if hwe is not None and hwe < HWE_ALPHA else ""

    # Control allele frequency against the reference (panel variants with named alleles)
    eaf = qc_control_eaf(cells) if cells else None
    q["qc_control_eaf"] = "" if eaf is None else "%.3f" % eaf
    freqs = ref.get(r["rsid"], {})
    ea = (r.get("effect_allele") or "").strip().upper()
    q["qc_reference_eaf"] = "%.3f" % freqs[ea] if ea in freqs else ""
    fc = qc_frequency_check(cells, freqs.get(ea)) if cells and ea in freqs else {"flag": None, "detail": "no reference frequency"}
    q["qc_frequency_flag"] = fc["flag"] or ""
    q["qc_frequency_detail"] = fc["detail"]

    q["qc_flags"] = ";".join(f for f in [
        "orientation_not_reproducible" if q["qc_orientation"] == "not_reproducible" else "",
        "orientation_other_allele" if q["qc_orientation"] == "other_allele" else "",
        q["qc_interval_flag"], q["qc_hwe_flag"], q["qc_frequency_flag"],
        "hand_transcribed" if r.get("extraction_method") == "agent_read_not_model_extracted" else ""] if f)
    return q


def priority(study, rows, no_pdf):
    if study in no_pdf:
        return "P0 no source"
    serious = sum(1 for q in rows if q["qc_orientation"] == "not_reproducible" or q["qc_interval_flag"] in
                  ("estimate_outside_own_interval", "interval_inverted") or q["qc_frequency_flag"])
    other = sum(1 for q in rows if q["qc_flags"]) - serious
    panel = any(q["qc_panel_variant"] == "yes" for q in rows)
    if serious:
        return "P1 panel" if panel else "P1"
    if other:
        return "P2 panel" if panel else "P2"
    return "P3 panel" if panel else "P3"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workbook", help="write a reviewer workbook (xlsx) to this path")
    args = ap.parse_args()

    rows = read_bundle("extraction_variants_round1.csv")
    original = [dict(r) for r in rows]
    studies = {s["study_id"]: s for s in read_bundle("extraction_studies_round1.csv")}
    corrected = defaultdict(set)
    for c in read_csv("data/extraction_corrections.csv"):
        corrected[c["variant_row_id"]].add(c["field"])
    panel, ref = panel_variants(), reference_freqs()
    no_pdf = {s for s, v in studies.items() if v.get("source_status", "").startswith(("VOID", "no_source"))}

    qc = [check_row(r, panel, ref, corrected) for r in rows]
    qc_assert_transcription_intact(rows, original, TRANSCRIBED)

    os.makedirs("extraction/output", exist_ok=True)
    cols = list(qc[0].keys()) + ["qc_orientation_reason"]
    cols = list(dict.fromkeys(cols))
    with open("extraction/output/qc_rows.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(qc)

    by_study = defaultdict(list)
    for q in qc:
        by_study[q["study_id"]].append(q)
    srows = []
    for sid in sorted(set(studies) | set(by_study)):
        rs = by_study.get(sid, [])
        s = studies.get(sid, {})
        srows.append({"study_id": sid, "first_author": s.get("first_author", ""), "year": s.get("year", ""),
                      "source_status": s.get("source_status", ""), "rows": len(rs),
                      "panel_rows": sum(1 for q in rs if q["qc_panel_variant"] == "yes"),
                      "rows_checkable_orientation": sum(1 for q in rs if q["qc_orientation"] != "not_checkable"),
                      "not_reproducible": sum(1 for q in rs if q["qc_orientation"] == "not_reproducible"),
                      "other_allele": sum(1 for q in rs if q["qc_orientation"] == "other_allele"),
                      "interval_flags": sum(1 for q in rs if q["qc_interval_flag"]),
                      "hwe_flags": sum(1 for q in rs if q["qc_hwe_flag"]),
                      "frequency_flags": sum(1 for q in rs if q["qc_frequency_flag"]),
                      "hand_transcribed": sum(1 for q in rs if "hand_transcribed" in q["qc_flags"]),
                      "already_corrected_rows": sum(1 for q in rs if q["qc_already_corrected"]),
                      "priority": priority(sid, rs, no_pdf) if rs or sid in no_pdf else "no rows"})
    order = {"P0": 0, "P1": 1, "P2": 2, "P3": 3, "no": 4}
    srows.sort(key=lambda s: (order[s["priority"][:2]], "panel" not in s["priority"], -s["not_reproducible"], s["study_id"]))
    with open("extraction/output/qc_studies.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(srows[0].keys()))
        w.writeheader()
        w.writerows(srows)

    n = len(qc)
    cnt = Counter(q["qc_orientation"] for q in qc)
    summary = [("variant rows", n),
               ("rows with six cells and a printed odds ratio (orientation checkable)", n - cnt["not_checkable"]),
               ("  reproduced for the effect allele", cnt["effect_allele"]),
               ("  reproduced only for the other allele (not an error)", cnt["other_allele"]),
               ("  not reproducible, crude estimate", cnt["not_reproducible"]),
               ("  adjusted estimate, not recomputable", cnt["adjusted_estimate_not_recomputable"]),
               ("rows with a complete interval", sum(1 for q in qc if q["qc_interval_detail"] != "interval incomplete")),
               ("  interval flags", sum(1 for q in qc if q["qc_interval_flag"])),
               ("    estimate outside its own interval", sum(1 for q in qc if q["qc_interval_flag"] == "estimate_outside_own_interval")),
               ("    interval inverted", sum(1 for q in qc if q["qc_interval_flag"] == "interval_inverted")),
               ("    log-scale arms asymmetric", sum(1 for q in qc if q["qc_interval_flag"] == "interval_asymmetric_on_log_scale")),
               ("rows with a Hardy-Weinberg test in controls", sum(1 for q in qc if q["qc_hwe_p_controls"])),
               ("  p < 0.001", sum(1 for q in qc if q["qc_hwe_flag"])),
               ("rows with a reference frequency (panel variants)", sum(1 for q in qc if q["qc_reference_eaf"] and q["qc_control_eaf"])),
               ("  closer to reference if homozygotes exchanged", sum(1 for q in qc if q["qc_frequency_flag"])),
               ("hand-transcribed rows", sum(1 for q in qc if "hand_transcribed" in q["qc_flags"])),
               ("rows already in extraction_corrections.csv", sum(1 for q in qc if q["qc_already_corrected"])),
               ("rows with any flag", sum(1 for q in qc if q["qc_flags"])),
               ("rows no check could see", sum(1 for q in qc if q["qc_orientation"] == "not_checkable" and
                                               q["qc_interval_detail"] == "interval incomplete" and not q["qc_hwe_p_controls"]))]
    with open("extraction/output/qc_summary.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["measure", "rows"])
        w.writerows(summary)
    for k, v in summary:
        print("%-75s %5d" % (k, v))
    print("studies by priority:", dict(Counter(s["priority"] for s in srows)))

    if args.workbook:
        flagged = [dict(q, **{"MATCHES PDF?": "", "CORRECTED VALUE": "", "READ FROM (page/table)": "", "CHECKER": "", "DATE": ""})
                   for q in qc if q["qc_flags"] and "orientation_other_allele" != q["qc_flags"]]
        prio = {s["study_id"]: s["priority"] for s in srows}
        flagged.sort(key=lambda q: (order[prio[q["study_id"]][:2]], q["qc_panel_variant"] != "yes", q["study_id"], q["variant_row_id"]))
        for q in flagged:
            q["priority"] = prio[q["study_id"]]
        row_cols = ["priority", "variant_row_id", "study_id", "rsid", "gene_symbol", "qc_panel_variant", "page_or_table", "qc_flags",
                    "qc_closest_recomputed_or", "qc_interval_detail", "qc_hwe_p_controls", "qc_frequency_detail",
                    "qc_already_corrected", "MATCHES PDF?", "CORRECTED VALUE", "READ FROM (page/table)", "CHECKER", "DATE"]
        study_cols = list(srows[0].keys()) + ["CHECKED BY", "DATE CHECKED", "OUTCOME", "NOTES"]
        readme = [("Arithmetic QC of the PROSE extraction table", 14, True),
                  ("Internal consistency is not verification. Agreement with a published estimate is not evidence a value is right: "
                   "two compensating errors reproduce it exactly as well as a correct reading does. A row that passes every check here has "
                   "passed a test against itself and has not been compared with its paper.", 10, False),
                  ("Blank means not checked, never confirmed.", 10, True),
                  ("What was checked: for odds ratios with all six genotype cells, which of six contrasts reproduces the printed estimate; "
                   "each estimate against its own interval; Hardy-Weinberg in controls; control allele frequency against gnomAD NFE for panel "
                   "variants. 'Reproduced only for the other allele' is not an error and is left out of the rows sheet.", 10, False),
                  ("Order of work: P0 studies need the paper first (procurement, not checking). Then P1 panel, P1, P2 panel, P2. Record "
                   "corrections in data/extraction_corrections.csv with the page or table read from; never edit the bundle.", 10, False)]
        qc_write_workbook(args.workbook, [
            {"name": "studies", "columns": study_cols, "rows": srows,
             "dropdowns": {"OUTCOME": ["no change", "corrected", "could not check"]}},
            {"name": "rows to check", "columns": row_cols, "rows": flagged,
             "dropdowns": {"MATCHES PDF?": ["yes", "no", "not found"]}}], readme=readme, priority_column="priority")
        print("workbook:", args.workbook)


if __name__ == "__main__":
    main()
