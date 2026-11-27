#!/usr/bin/env python3
"""Run the two checking agents on Russ's Mac, where the PDFs are.

  extract  Blind re-extraction. The agent gets the PDF and the list of variants
           to look for, never the round-1 values.
  judge    For each disagreement that compare.py found, a second agent reads the
           PDF and says which of two unlabelled extractions matches it.

Every request and response is saved under extraction/runs/<run>/ together with
the model that answered, the prompt hash and the token usage, so the run can be
described in the methods and audited later. Nothing here touches the bundle or
data/extraction_corrections.csv.

  python3 extraction/run_agents.py extract --pdf-dir ~/PROSE_pdfs --scope panel --run pilot
  python3 extraction/compare.py --run pilot
  python3 extraction/run_agents.py judge --pdf-dir ~/PROSE_pdfs --run pilot
  python3 extraction/compare.py --run pilot --workbook

Needs `pip install anthropic pypdf` and an API key (ANTHROPIC_API_KEY or `ant auth login`).
--dry-run builds every request and prints its size without sending anything.
"""
import argparse, base64, datetime, hashlib, json, os, random, sys
from collections import OrderedDict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import HERE, corrected_variants, panel_rsids, sources, norm_rsid, read_csv
from schema import extraction_schema, judge_schema

DEFAULT_MODEL = "claude-opus-5-5"


def prompt(name):
    with open(os.path.join(HERE, "prompts", name), encoding="utf-8") as f:
        return f.read()


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()[:12]


def run_dir(run):
    d = os.path.join(HERE, "runs", run)
    os.makedirs(d, exist_ok=True)
    return d


def excluded_studies():
    """Studies excluded from the evidence base when read for extraction (eligible_confirmed = no)."""
    path = os.path.join(os.path.dirname(HERE), "data", "extraction_studies_corrected.csv")
    if not os.path.exists(path):
        return set()
    return {r["study_id"] for r in read_csv(path) if r.get("eligible_confirmed", "").strip().lower() == "no"}


def targets(scope, studies):
    """Variants to look for in each study: names only, no values. Excluded studies are left out."""
    rows, _ = corrected_variants()
    panel = panel_rsids()
    skip = excluded_studies()
    out = OrderedDict()
    for r in sorted(rows, key=lambda r: (r["study_id"], r["variant_row_id"])):
        if studies and r["study_id"] not in studies:
            continue
        if r["study_id"] in skip:
            continue
        rs = norm_rsid(r["rsid"])
        if scope == "panel" and rs not in panel:
            continue
        key = (r["gene_symbol"].strip().upper(), rs, r["variant_name_as_reported"].strip())
        out.setdefault(r["study_id"], OrderedDict())[key] = True
    return {s: list(v) for s, v in out.items()}


def pdf_bytes(path, pages):
    """The whole file, or only the listed pages (for abstracts inside a conference supplement)."""
    if not pages:
        with open(path, "rb") as f:
            return f.read()
    import io
    from pypdf import PdfReader, PdfWriter
    lo, _, hi = pages.partition("-")
    reader, writer = PdfReader(path), PdfWriter()
    for i in range(int(lo) - 1, int(hi or lo)):
        writer.add_page(reader.pages[i])
    buf = io.BytesIO()
    writer.write(buf)
    return buf.getvalue()


def pdf_blocks(files):
    blocks = []
    for p, pages in files:
        data = base64.standard_b64encode(pdf_bytes(p, pages)).decode()
        title = os.path.basename(p) + (" (pages %s)" % pages if pages else "")
        blocks.append({"type": "document", "title": title,
                       "source": {"type": "base64", "media_type": "application/pdf", "data": data}})
    return blocks


def call(client, model, system, content, schema, effort):
    """One structured-output request. Returns (parsed JSON or None, record)."""
    with client.beta.messages.stream(
        model=model, max_tokens=64000, system=system,
        messages=[{"role": "user", "content": content}],
        output_config={"effort": effort, "format": {"type": "json_schema", "schema": schema}},
        betas=["server-side-fallback-2026-07-01"], fallbacks="default",
    ) as stream:
        msg = stream.get_final_message()
    served = [msg.model]
    for b in msg.content:
        if b.type == "fallback":
            served.append(b.to.model)
    rec = {"request_id": getattr(msg, "_request_id", None), "model_requested": model,
           "model_served": served[-1], "fallback_used": len(served) > 1,
           "stop_reason": msg.stop_reason,
           "usage": {"input": msg.usage.input_tokens, "output": msg.usage.output_tokens},
           "finished_at": datetime.datetime.now().isoformat(timespec="seconds")}
    if msg.stop_reason in ("refusal", "max_tokens"):
        return None, rec
    text = next((b.text for b in msg.content if b.type == "text"), "")
    try:
        return json.loads(text), rec
    except json.JSONDecodeError:
        rec["error"] = "response was not valid JSON"
        return None, rec


def client_or_none(dry):
    if dry:
        return None
    import anthropic
    return anthropic.Anthropic()


def cmd_extract(a):
    studies = set(a.studies.split(",")) if a.studies else None
    want = targets(a.scope, studies)
    src = sources(os.path.expanduser(a.pdf_dir))
    d = run_dir(a.run)
    system = prompt("extract_system.md")
    schema = extraction_schema()
    client = client_or_none(a.dry_run)
    log = []
    for study, variants in want.items():
        out_path = os.path.join(d, "extract_%s.json" % study)
        if os.path.exists(out_path) and not a.redo:
            print("%s already done, skipping (use --redo)" % study)
            continue
        files = src.get(study, [])
        missing = [f for f, _ in files if not os.path.exists(f)]
        if not files or missing:
            print("%s: no source PDF%s" % (study, (" (missing: %s)" % missing[0]) if missing else ""))
            log.append({"study_id": study, "status": "no_source"})
            continue
        listing = "\n".join("- %s %s (%s)" % (g, rs or "no rsID", name) for g, rs, name in variants)
        text = ("Study %s. Extract every result the paper reports for these variants:\n%s\n"
                "Return the rows, the variants not found, and a one-line paper_note." % (study, listing))
        content = pdf_blocks(files) + [{"type": "text", "text": text}]
        size = sum(len(pdf_bytes(f, p)) for f, p in files)
        if a.dry_run:
            print("%s: %d variants, %d PDF(s), %.1f MB" % (study, len(variants), len(files), size / 1e6))
            continue
        print("%s: extracting %d variants ..." % (study, len(variants)), flush=True)
        try:
            result, rec = call(client, a.model, system, content, schema, a.effort)
        except Exception as e:  # one bad request (e.g. a PDF over the page limit) should not stop the run
            if type(e).__name__ in ("AuthenticationError", "PermissionDeniedError", "TypeError"):
                raise
            print("   request failed, skipped: %s" % str(e)[:300])
            log.append({"study_id": study, "status": "request_failed", "error": str(e)[:500]})
            continue
        rec.update({"study_id": study, "files": [os.path.basename(f) + (":" + p if p else "") for f, p in files],
                    "prompt_sha": sha(system), "scope": a.scope, "n_variants": len(variants)})
        log.append(rec)
        if result is not None:
            with open(out_path, "w", encoding="utf-8") as f:
                json.dump({"meta": rec, "targets": variants, **result}, f, indent=1)
            print("   %d rows, %d not found" % (len(result["rows"]), len(result["not_found"])))
        else:
            print("   no result: %s" % (rec.get("error") or rec["stop_reason"]))
    append_log(d, "extract", log)


def cmd_judge(a):
    d = run_dir(a.run)
    path = os.path.join(d, "disagreements.csv")
    if not os.path.exists(path):
        sys.exit("Run compare.py --run %s first." % a.run)
    rows = read_csv(path)
    src = sources(os.path.expanduser(a.pdf_dir))
    system = prompt("judge_system.md")
    schema = judge_schema()
    client = client_or_none(a.dry_run)
    rng = random.Random(a.seed)
    by_study = OrderedDict()
    for r in rows:
        by_study.setdefault(r["study_id"], OrderedDict()).setdefault(r["item_id"], []).append(r)
    log = []
    for study, items in by_study.items():
        out_path = os.path.join(d, "judge_%s.json" % study)
        if os.path.exists(out_path) and not a.redo:
            print("%s already judged, skipping (use --redo)" % study)
            continue
        files = [(f, p) for f, p in src.get(study, []) if os.path.exists(f)]
        if not files:
            print("%s: no source PDF" % study)
            continue
        blinding, blocks = {}, []
        for item_id, fields in items.items():
            first = fields[0]
            a_is_original = rng.random() < 0.5
            blinding[item_id] = "A=original" if a_is_original else "A=reextraction"
            lines = []
            for f in fields:
                va, vb = (f["original"], f["reextraction"]) if a_is_original else (f["reextraction"], f["original"])
                lines.append("  %s: A=%s | B=%s" % (f["field"], va or "(blank)", vb or "(blank)"))
            blocks.append("Item %s. paper_result: %s %s, %s, %s, subgroup %s\n%s" % (
                item_id, first["gene_symbol"], first["rsid"], first["genetic_model"],
                first["comparison"], first["subgroup"], "\n".join(lines)))
        text = "Study %s. Check these items against the paper:\n\n%s" % (study, "\n\n".join(blocks))
        if a.dry_run:
            print("%s: %d items" % (study, len(items)))
            continue
        print("%s: judging %d items ..." % (study, len(items)), flush=True)
        try:
            result, rec = call(client, a.model, system, pdf_blocks(files) + [{"type": "text", "text": text}],
                               schema, a.effort)
        except Exception as e:
            if type(e).__name__ in ("AuthenticationError", "PermissionDeniedError", "TypeError"):
                raise
            print("   request failed, skipped: %s" % str(e)[:300])
            log.append({"study_id": study, "status": "request_failed", "error": str(e)[:500]})
            continue
        rec.update({"study_id": study, "prompt_sha": sha(system), "n_items": len(items)})
        log.append(rec)
        if result is not None:
            with open(out_path, "w", encoding="utf-8") as f:
                json.dump({"meta": rec, "blinding": blinding, **result}, f, indent=1)
    append_log(d, "judge", log)


def append_log(d, step, log):
    if not log:
        return
    with open(os.path.join(d, "run_log.jsonl"), "a", encoding="utf-8") as f:
        for r in log:
            f.write(json.dumps({"step": step, **r}) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    for name in ("extract", "judge"):
        s = sub.add_parser(name)
        s.add_argument("--pdf-dir", required=True, help="folder holding the source PDFs")
        s.add_argument("--run", default="pilot", help="name of the run folder under extraction/runs/")
        s.add_argument("--model", default=DEFAULT_MODEL)
        s.add_argument("--effort", default="high", choices=["low", "medium", "high", "xhigh", "max"])
        s.add_argument("--dry-run", action="store_true")
        s.add_argument("--redo", action="store_true", help="repeat studies already done in this run")
        if name == "extract":
            s.add_argument("--scope", default="panel", choices=["panel", "all"])
            s.add_argument("--studies", help="comma-separated study ids, e.g. S026,S071")
        else:
            s.add_argument("--seed", type=int, default=2026, help="seed for the A/B blinding")
    a = p.parse_args()
    {"extract": cmd_extract, "judge": cmd_judge}[a.cmd](a)


if __name__ == "__main__":
    main()
