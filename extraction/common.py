"""Shared helpers for the extraction check scripts.

Everything here reads. Nothing writes to the bundle or to the corrections log.
"""
import csv, io, os, re, zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUNDLE = os.path.join(ROOT, "IMD_genetics_analysis_bundle (1).zip")
PREFIX = "IMD_genetics_bundle/data/"
CORRECTIONS = os.path.join(ROOT, "data", "extraction_corrections.csv")
HERE = os.path.dirname(os.path.abspath(__file__))

# Values copied from a paper. The agents re-read these; nothing else is compared.
COUNT_FIELDS = ["n_cases_hom_effect", "n_cases_het", "n_cases_hom_other",
                "n_ctrl_hom_effect", "n_ctrl_het", "n_ctrl_hom_other"]
NUMBER_FIELDS = ["effect_estimate", "ci_lower", "ci_upper"]
TEXT_FIELDS = ["effect_allele", "other_allele", "genetic_model", "comparison_as_reported",
               "subgroup", "effect_measure", "p_value", "adjusted", "page_or_table"]
COMPARED_FIELDS = (["effect_allele", "other_allele", "genetic_model", "effect_measure", "adjusted"]
                   + COUNT_FIELDS + NUMBER_FIELDS + ["p_value"])


def read_bundle(name):
    with zipfile.ZipFile(BUNDLE) as z:
        with z.open(PREFIX + name) as f:
            return list(csv.DictReader(io.TextIOWrapper(f, encoding="utf-8-sig")))


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, columns=None):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    columns = columns or (list(rows[0].keys()) if rows else [])
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return path


def corrected_variants():
    """Round-1 variant rows with the corrections log applied field by field.

    This is the table the analysis actually uses, so it is what the agents are
    checked against. Returns (rows, set of corrected row ids).
    """
    rows = read_bundle("extraction_variants_round1.csv")
    by_id = {r["variant_row_id"]: r for r in rows}
    fixed = set()
    if os.path.exists(CORRECTIONS):
        for c in read_csv(CORRECTIONS):
            r = by_id.get(c["variant_row_id"])
            if r is not None and c["field"] in r:
                r[c["field"]] = c["corrected_value"]
                fixed.add(c["variant_row_id"])
    return rows, fixed


def panel_rsids(min_studies=3):
    return {r["rsid"].strip().lower() for r in read_bundle("panel_candidates.csv")
            if int(r["studies_inclusive"] or 0) >= min_studies}


def num(x):
    try:
        s = str(x).strip().replace("−", "-")
        return float(s) if s not in ("", "NA", "None", "null") else None
    except ValueError:
        return None


def norm_rsid(x):
    m = re.search(r"rs\s*(\d+)", str(x or ""), re.I)
    return "rs" + m.group(1) if m else ""


def norm_text(x):
    return re.sub(r"\s+", " ", str(x or "").strip().lower())


def sources(pdf_dir=None):
    """Map each study to its source PDF(s).

    extraction/sources.csv wins when present (one row per file, so a study can
    list a supplement too). Otherwise the bundle's pdf_source_audit is used and
    VOID / no_source studies are left out. Each entry is (file, pages), where
    pages is blank for the whole file or a range such as "322-323" for an
    abstract inside a conference supplement.
    """
    path = os.path.join(HERE, "sources.csv")
    out = {}
    if os.path.exists(path):
        for r in read_csv(path):
            if r.get("file"):
                out.setdefault(r["study_id"], []).append((r["file"], r.get("pages", "")))
    else:
        for r in read_bundle("pdf_source_audit.csv"):
            if r["pdf_used"] and r["source_status"] in ("verified", "re-extracted_round1b"):
                out.setdefault(r["study_id"], []).append((r["pdf_used"], ""))
    if pdf_dir:
        out = {s: [(os.path.join(pdf_dir, f), p) for f, p in fs] for s, fs in out.items()}
    return out
