#!/usr/bin/env python3
"""Compare the blind re-extraction with the corrected round-1 table, then build the reviewer workbook.

  python3 extraction/compare.py --run pilot              # writes disagreements.csv for the judge
  python3 extraction/compare.py --run pilot --workbook   # after the judge: workbook + proposed corrections

Rows are paired within each study and variant on genetic model, comparison,
subgroup and values. Where the two extractions name the alleles the other way
round, the re-extraction's genotype cells are reversed and its estimate inverted
before comparing, because that is the same result written from the other allele.

Agreement between two extractions is not proof that either is right: both can
misread the same table. The judge reads the paper for every disagreement, and a
human confirms each change before it goes into data/extraction_corrections.csv.
"""
import argparse, datetime, difflib, glob, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (HERE, COUNT_FIELDS, COMPARED_FIELDS, corrected_variants, panel_rsids,
                    norm_rsid, norm_text, num, read_csv, write_csv)

LOOSE = {"adjusted": {"not_stated": "unclear"}}


def same_number(a, b, rel=0.01, abs_=0.011):
    x, y = num(a), num(b)
    if x is None or y is None:
        return x is None and y is None
    return abs(x - y) <= max(abs_, rel * max(abs(x), abs(y)))


def norm_p(p):
    s = norm_text(p).replace(" ", "").replace("p=", "").replace("p", "")
    s = s.lstrip("=")
    if s.startswith("0.") or s.startswith("<0.") or s.startswith(">0."):
        s = s.replace("0.", ".", 1)
    try:
        return "%.3g" % float(s)
    except ValueError:
        return s


def flipped(r):
    """The same result expressed for the other allele."""
    f = dict(r)
    f["effect_allele"], f["other_allele"] = r.get("other_allele"), r.get("effect_allele")
    for grp in ("cases", "ctrl"):
        f["n_%s_hom_effect" % grp], f["n_%s_hom_other" % grp] = r.get("n_%s_hom_other" % grp), r.get("n_%s_hom_effect" % grp)
    est, lo, hi = num(r.get("effect_estimate")), num(r.get("ci_lower")), num(r.get("ci_upper"))
    ratio = str(r.get("effect_measure", "")).upper() in ("OR", "HR", "RR", "SUBDISTRIBUTION_HR")
    if ratio and est:
        f["effect_estimate"] = 1 / est
        f["ci_lower"], f["ci_upper"] = (1 / hi if hi else None), (1 / lo if lo else None)
    return f


def field_equal(field, o, r):
    a, b = o.get(field), r.get(field)
    if field in COUNT_FIELDS:
        x, y = num(a), num(b)
        return x == y
    if field in ("effect_estimate", "ci_lower", "ci_upper"):
        return same_number(a, b, rel=0.02 if r.get("_flipped") else 0.01)
    if field == "p_value":
        return norm_p(a) == norm_p(b)
    x, y = norm_text(a), norm_text(b)
    m = LOOSE.get(field, {})
    return m.get(x, x) == m.get(y, y)


def diff_fields(o, r):
    return [f for f in COMPARED_FIELDS if not field_equal(f, o, r)]


def oriented(o, r):
    """Return the re-extracted row in the original's allele orientation."""
    oe, oo = norm_text(o.get("effect_allele")), norm_text(o.get("other_allele"))
    re_, ro = norm_text(r.get("effect_allele")), norm_text(r.get("other_allele"))
    if oe and re_ and oe == ro and oo == re_ and oe != oo:
        f = flipped(r)
        f["_flipped"] = True
        return f
    return dict(r)


def pair_score(o, r):
    rr = oriented(o, r)
    s = 3.0 * (norm_text(o["genetic_model"]) == norm_text(r["genetic_model"]))
    s += 3.0 * difflib.SequenceMatcher(None, norm_text(o["comparison_as_reported"]),
                                       norm_text(r["comparison_as_reported"])).ratio()
    s += 2.0 * difflib.SequenceMatcher(None, norm_text(o["subgroup"]), norm_text(r["subgroup"])).ratio()
    s += 2.0 * (num(o["effect_estimate"]) is not None and same_number(o["effect_estimate"], rr["effect_estimate"]))
    s += 2.0 * all(field_equal(f, o, rr) for f in COUNT_FIELDS)
    return s


def variant_key(gene, rsid, name):
    rs = norm_rsid(rsid)
    return rs if rs else "%s|%s" % (norm_text(gene), norm_text(name))


def load_extractions(d):
    out = {}
    for p in sorted(glob.glob(os.path.join(d, "extract_*.json"))):
        with open(p, encoding="utf-8") as f:
            j = json.load(f)
        out[j["meta"]["study_id"]] = j
    return out


def compare(run):
    d = os.path.join(HERE, "runs", run)
    ext = load_extractions(d)
    if not ext:
        sys.exit("No extract_*.json files in %s. Run run_agents.py extract first." % d)
    rows, fixed = corrected_variants()
    panel = panel_rsids()
    pairs, disagreements = [], []
    n_item = 0
    for study, j in ext.items():
        scope = j["meta"].get("scope", "panel")
        orig = [r for r in rows if r["study_id"] == study
                and (scope == "all" or norm_rsid(r["rsid"]) in panel)]
        new = [dict(r, _i=i) for i, r in enumerate(j["rows"])]
        keys = {variant_key(o["gene_symbol"], o["rsid"], o["variant_name_as_reported"]) for o in orig}
        keys |= {variant_key(r["gene_symbol"], r["rsid"], r["variant_name_as_reported"]) for r in new}
        for k in sorted(keys):
            O = [o for o in orig if variant_key(o["gene_symbol"], o["rsid"], o["variant_name_as_reported"]) == k]
            R = [r for r in new if variant_key(r["gene_symbol"], r["rsid"], r["variant_name_as_reported"]) == k]
            cand = sorted(((pair_score(o, r), oi, ri) for oi, o in enumerate(O) for ri, r in enumerate(R)),
                          reverse=True)
            used_o, used_r, matched = set(), set(), []
            for s, oi, ri in cand:
                if s >= 3.0 and oi not in used_o and ri not in used_r:
                    used_o.add(oi); used_r.add(ri); matched.append((O[oi], R[ri], s))
            for o, r, s in matched:
                rr = oriented(o, r)
                diffs = diff_fields(o, rr)
                n_item += 1
                item = "I%04d" % n_item
                pairs.append(pair_row(item, study, o, r, "agree" if not diffs else "disagree", s, diffs, rr))
                for f in diffs:
                    disagreements.append(dis_row(item, study, o, f, o.get(f), r.get(f), rr.get("_flipped")))
            for oi, o in enumerate(O):
                if oi not in used_o:
                    n_item += 1
                    item = "I%04d" % n_item
                    pairs.append(pair_row(item, study, o, None, "original_only", 0, [], None))
                    disagreements.append(dis_row(item, study, o, "row_present", "reported in paper", "", False))
            for ri, r in enumerate(R):
                if ri not in used_r:
                    n_item += 1
                    item = "I%04d" % n_item
                    pairs.append(pair_row(item, study, None, r, "reextraction_only", 0, [], None))
                    disagreements.append(dis_row(item, study, r, "row_present", "", "reported in paper", False))
    for p in pairs:
        p["already_corrected"] = "yes" if p["variant_row_id"] in fixed else ""
    write_csv(os.path.join(d, "pairs.csv"), pairs)
    write_csv(os.path.join(d, "disagreements.csv"), disagreements,
              ["item_id", "study_id", "variant_row_id", "gene_symbol", "rsid", "genetic_model", "comparison",
               "subgroup", "field", "original", "reextraction", "alleles_named_other_way"])
    summary = summarise(pairs, disagreements)
    write_csv(os.path.join(d, "agreement_summary.csv"), summary, ["measure", "value"])
    for s in summary:
        print("%-46s %s" % (s["measure"], s["value"]))
    return pairs, disagreements


def pair_row(item, study, o, r, status, score, diffs, rr):
    base = o or r
    return {"item_id": item, "study_id": study, "variant_row_id": (o or {}).get("variant_row_id", ""),
            "gene_symbol": base.get("gene_symbol", ""), "rsid": base.get("rsid", ""),
            "genetic_model": base.get("genetic_model", ""),
            "comparison": base.get("comparison_as_reported", ""), "subgroup": base.get("subgroup", ""),
            "status": status, "match_score": "%.1f" % score, "fields_differing": ";".join(diffs),
            "alleles_named_other_way": "yes" if rr and rr.get("_flipped") else "",
            "reextraction_confidence": (r or {}).get("confidence", ""),
            "reextraction_page": (r or {}).get("page_or_table", ""),
            "reextraction_source_text": (r or {}).get("verbatim_source", ""),
            "reextraction_note": (r or {}).get("note", "")}


def dis_row(item, study, base, field, original, reextraction, flip):
    return {"item_id": item, "study_id": study, "variant_row_id": base.get("variant_row_id", ""),
            "gene_symbol": base.get("gene_symbol", ""), "rsid": base.get("rsid", ""),
            "genetic_model": base.get("genetic_model", ""), "comparison": base.get("comparison_as_reported", ""),
            "subgroup": base.get("subgroup", ""), "field": field,
            "original": "" if original is None else str(original),
            "reextraction": "" if reextraction is None else str(reextraction),
            "alleles_named_other_way": "yes" if flip else ""}


def summarise(pairs, dis):
    from collections import Counter
    st = Counter(p["status"] for p in pairs)
    matched = st["agree"] + st["disagree"]
    out = [{"measure": "studies compared", "value": len({p["study_id"] for p in pairs})},
           {"measure": "original rows in scope", "value": matched + st["original_only"]},
           {"measure": "rows paired with a re-extracted row", "value": matched},
           {"measure": "paired rows agreeing on every field", "value": st["agree"]},
           {"measure": "paired rows with at least one difference", "value": st["disagree"]},
           {"measure": "original rows the re-extraction did not find", "value": st["original_only"]},
           {"measure": "re-extracted rows not in the original", "value": st["reextraction_only"]}]
    by_field = Counter(x["field"] for x in dis if x["field"] != "row_present")
    for f, n in by_field.most_common():
        out.append({"measure": "differences in " + f, "value": n})
    return out


# ---------------------------------------------------------------- after the judge

def unblind(d):
    verdicts = {}
    for p in sorted(glob.glob(os.path.join(d, "judge_*.json"))):
        with open(p, encoding="utf-8") as f:
            j = json.load(f)
        for it in j["items"]:
            a_orig = j["blinding"].get(it["item_id"]) == "A=original"
            for fv in it["fields"]:
                v = fv["verdict"]
                who = {"A": "original" if a_orig else "reextraction",
                       "B": "reextraction" if a_orig else "original"}.get(v, v)
                verdicts[(it["item_id"], fv["field"])] = {
                    "judge_supports": who, "value_in_source": fv.get("value_in_source") or "",
                    "judge_quote": fv.get("quote", ""), "judge_location": fv.get("page_or_table", ""),
                    "judge_note": it.get("item_note", ""), "judge_model": j["meta"].get("model_served", "")}
    return verdicts


def workbook(run):
    d = os.path.join(HERE, "runs", run)
    pairs, dis = compare(run)
    verdicts = unblind(d)
    if not verdicts:
        print("No judge_*.json yet; the workbook will show differences without verdicts.")
    for x in dis:
        x.update(verdicts.get((x["item_id"], x["field"]), {"judge_supports": "not judged"}))
    rank = {"reextraction": 0, "neither": 1, "not_found": 2, "not judged": 3, "original": 4, "both": 5}
    dis.sort(key=lambda x: (rank.get(x["judge_supports"], 3), x["study_id"], x["item_id"]))
    write_csv(os.path.join(d, "disagreements_judged.csv"), dis)

    today = datetime.date.today().isoformat()
    proposed = [{"variant_row_id": x["variant_row_id"], "study_id": x["study_id"], "field": x["field"],
                 "original_value": x["original"], "corrected_value": x["value_in_source"] or x["reextraction"],
                 "evidence": "Blind re-extraction and judge (run %s) both read: %s" % (run, x.get("judge_quote", "")),
                 "source_location": x.get("judge_location", ""), "checked_on": "", "checked_by": ""}
                for x in dis if x["judge_supports"] == "reextraction" and x["variant_row_id"]
                and x["field"] != "row_present"]
    write_csv(os.path.join(d, "proposed_corrections.csv"), proposed,
              ["variant_row_id", "study_id", "field", "original_value", "corrected_value", "evidence",
               "source_location", "checked_on", "checked_by"])

    from collections import Counter
    jc = Counter(x["judge_supports"] for x in dis)
    readme = [
        ("Blind re-extraction check, run %s (%s)" % (run, today), 14, True),
        ("Internal consistency is not verification. Two extractions that agree can both be wrong, "
         "because they can misread the same table the same way. Agreement here lowers the chance of an "
         "error; it does not prove a value.", 10, True),
        ("What was done: an agent read each source PDF with no access to the round-1 values and extracted "
         "the listed variants. Its rows were paired with the corrected round-1 table. Every difference went "
         "to a second agent, which read the PDF and said which of two unlabelled extractions matched the paper, "
         "quoting the text it read.", 10, False),
        ("What it did not do: nothing was changed. The bundle and data/extraction_corrections.csv are "
         "untouched. proposed_corrections.csv in the run folder lists the changes the judge supports, in the "
         "corrections-log format, with checked_on and checked_by left empty. A person opens the PDF, confirms "
         "each one, fills those two columns and copies the confirmed lines into the corrections log.", 10, False),
        ("Order of work: the Differences sheet is sorted with the rows where the judge sided with the "
         "re-extraction first, then 'neither', then 'not found'. Rows where the judge sided with the original "
         "need no change but are listed so the sample can be audited.", 10, False),
        ("Blank in a reviewer column means not checked, never confirmed.", 10, True),
        ("Judge verdicts: " + ", ".join("%s %d" % (k, v) for k, v in jc.most_common()), 10, False),
    ]
    summary_rows = [{"measure": s["measure"], "value": s["value"]} for s in summarise(pairs, dis)]
    review = ["MATCHES PDF?", "CORRECTED VALUE", "READ FROM (page/table)", "CHECKER", "DATE"]
    dcols = ["judge_supports", "study_id", "variant_row_id", "gene_symbol", "rsid", "genetic_model",
             "comparison", "subgroup", "field", "original", "reextraction", "value_in_source", "judge_quote",
             "judge_location", "judge_note", "alleles_named_other_way"] + review
    pcols = ["status", "study_id", "variant_row_id", "gene_symbol", "rsid", "genetic_model", "comparison",
             "subgroup", "fields_differing", "alleles_named_other_way", "already_corrected",
             "reextraction_confidence", "reextraction_page", "reextraction_source_text", "reextraction_note"]
    from qc_kernel import qc_write_workbook
    out = os.path.join(d, "reextraction_review_%s.xlsx" % run)
    qc_write_workbook(out, [
        {"name": "Summary", "columns": ["measure", "value"], "rows": summary_rows, "widths": [52, 12]},
        {"name": "Differences", "columns": dcols, "rows": dis, "highlight_column": "judge_supports",
         "dropdowns": {"MATCHES PDF?": ["original", "re-extraction", "neither", "cannot tell"]}},
        {"name": "All paired rows", "columns": pcols, "rows": pairs},
    ], readme=readme)
    print("Workbook: %s" % out)
    print("Proposed corrections (unconfirmed): %d" % len(proposed))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", default="pilot")
    p.add_argument("--workbook", action="store_true")
    a = p.parse_args()
    workbook(a.run) if a.workbook else compare(a.run)


if __name__ == "__main__":
    main()
