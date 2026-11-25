"""Helpers for arithmetic QC of an extraction table. Flag-only: nothing here writes."""
import math

QC_TOLERANCE = 1.25
QC_CONTRAST_NAMES = ("recessive_for_effect", "dominant_for_effect", "allelic_for_effect",
                     "recessive_for_other", "dominant_for_other", "allelic_for_other")


def qc_num(x):
    """Coerce to float, or None if not numeric."""
    try:
        return float(str(x).strip())
    except Exception:
        return None


def qc_odds(x1, x0, y1, y0):
    """Odds ratio (x1/x0)/(y1/y0) with a 0.5 continuity correction when any cell is zero."""
    vals = [x1, x0, y1, y0]
    if any(v is None for v in vals):
        return None
    if min(vals) == 0:
        x1, x0, y1, y0 = x1 + 0.5, x0 + 0.5, y1 + 0.5, y0 + 0.5
    if 0 in (x0, y1, y0):
        return None
    return (x1 / x0) / (y1 / y0)


def qc_contrast_set(cells):
    """Six contrast forms from [case_hom_eff, case_het, case_hom_oth, ctrl_hom_eff, ctrl_het, ctrl_hom_oth].

    Each genetic model is computed both as the paper may express it for the effect
    allele and as it may express it for the other allele, because the transcribed
    model word is relative to whichever allele the paper grouped.
    """
    if cells is None or len(cells) != 6 or any(c is None for c in cells):
        return {}
    a, b, c, d, e, f = cells
    return {"recessive_for_effect": qc_odds(a, b + c, d, e + f),
            "dominant_for_effect": qc_odds(a + b, c, d + e, f),
            "allelic_for_effect": qc_odds(2 * a + b, 2 * c + b, 2 * d + e, 2 * f + e),
            "recessive_for_other": qc_odds(c, a + b, f, d + e),
            "dominant_for_other": qc_odds(b + c, a, e + f, d),
            "allelic_for_other": qc_odds(2 * c + b, 2 * a + b, 2 * f + e, 2 * d + e)}


def qc_reconcile_row(cells, reported, tolerance=None, is_adjusted=False):
    """Which contrast forms reproduce `reported` from `cells`.

    Returns {orientation, reproducing_contrast, matches, ratios}. Orientation is
    effect_allele / other_allele / not_reproducible / adjusted_estimate_not_recomputable
    / not_checkable. Never decides a row is wrong — a mismatch is a flag for a human.
    """
    if tolerance is None:
        tolerance = QC_TOLERANCE
    rep = qc_num(reported)
    ratios = qc_contrast_set(cells)
    if rep is None or rep <= 0 or not ratios:
        return {"orientation": "not_checkable", "reproducing_contrast": "", "matches": [], "ratios": ratios}
    lim = math.log(tolerance)
    hits = [k for k in QC_CONTRAST_NAMES
            if ratios.get(k) and ratios[k] > 0 and abs(math.log(ratios[k] / rep)) <= lim]
    eff = [k for k in hits if k.endswith("_for_effect")]
    if eff:
        orient, which = "effect_allele", eff[0]
    elif hits:
        orient, which = "other_allele", hits[0]
    elif is_adjusted:
        orient, which = "adjusted_estimate_not_recomputable", ""
    else:
        orient, which = "not_reproducible", ""
    return {"orientation": orient, "reproducing_contrast": which, "matches": hits, "ratios": ratios}


def qc_hwe_p(hom_a, het, hom_b, min_n=20):
    """Chi-square (1 df) Hardy-Weinberg p for one genotype triple, or None if underpowered."""
    vals = [qc_num(hom_a), qc_num(het), qc_num(hom_b)]
    if any(v is None for v in vals):
        return None
    obs_aa, obs_ab, obs_bb = vals
    n = obs_aa + obs_ab + obs_bb
    if n < min_n:
        return None
    p = (2 * obs_aa + obs_ab) / (2 * n)
    exp = [n * p * p, 2 * n * p * (1 - p), n * (1 - p) * (1 - p)]
    if min(exp) < 1:
        return None
    chi = sum((o - e) ** 2 / e for o, e in zip([obs_aa, obs_ab, obs_bb], exp))
    return math.erfc(math.sqrt(chi / 2.0))


def qc_ci_check(estimate, lower, upper, asymmetry=1.2):
    """Flag an estimate outside its own interval, an inverted interval, or lopsided log-scale arms."""
    est, lo, hi = qc_num(estimate), qc_num(lower), qc_num(upper)
    if None in (est, lo, hi):
        return {"flag": None, "detail": "interval incomplete"}
    if lo > hi:
        return {"flag": "interval_inverted", "detail": "lower limit exceeds upper limit"}
    if not (lo <= est <= hi):
        return {"flag": "estimate_outside_own_interval", "detail": "estimate %g lies outside %g to %g" % (est, lo, hi)}
    if min(est, lo, hi) > 0:
        left, right = math.log(est) - math.log(lo), math.log(hi) - math.log(est)
        if min(left, right) > 0 and max(left, right) / min(left, right) > asymmetry:
            return {"flag": "interval_asymmetric_on_log_scale",
                    "detail": "log-scale arms differ by %.2fx" % (max(left, right) / min(left, right))}
    return {"flag": None, "detail": "consistent"}


def qc_control_eaf(cells):
    """Effect-allele frequency among controls, from the last three cells."""
    if cells is None or len(cells) != 6 or any(c is None for c in cells[3:]):
        return None
    d, e, f = cells[3:]
    n = d + e + f
    return None if n == 0 else (2 * d + e) / (2 * n)


def qc_frequency_check(cells, reference_freq, shift=0.15):
    """Flag when exchanging the homozygote cells moves the control frequency markedly closer to reference."""
    ref = qc_num(reference_freq)
    if ref is None or cells is None or len(cells) != 6:
        return {"flag": None, "detail": "no reference frequency"}
    here = qc_control_eaf(cells)
    there = qc_control_eaf([cells[2], cells[1], cells[0], cells[5], cells[4], cells[3]])
    if here is None or there is None:
        return {"flag": None, "detail": "control counts incomplete"}
    if abs(there - ref) + shift < abs(here - ref):
        return {"flag": "frequency_closer_if_homozygotes_exchanged",
                "detail": "control frequency %.3f vs reference %.3f; %.3f if exchanged" % (here, ref, there)}
    return {"flag": None, "detail": "control frequency %.3f vs reference %.3f" % (here, ref)}


def qc_assert_transcription_intact(rows, original_rows, transcribed_columns, key="variant_row_id"):
    """Raise unless every transcribed column still matches the loaded table. Run before saving."""
    base = {r[key]: r for r in original_rows}
    drift = []
    for r in rows:
        o = base.get(r[key])
        if o is None:
            drift.append((r[key], "row absent from the original table"))
            continue
        for c in transcribed_columns:
            if str(r.get(c)) != str(o.get(c)):
                drift.append((r[key], "%s: %r -> %r" % (c, o.get(c), r.get(c))))
    if drift:
        raise AssertionError("transcribed values changed in %d place(s): %s" % (len(drift), drift[:10]))
    return True


def qc_write_workbook(path, sheets, readme=None, priority_column=None, highlight_column=None):
    """Write a reviewer workbook.

    sheets: list of {name, columns, rows, widths?, freeze?, dropdowns?}. dropdowns maps
    a column name to its allowed values. readme: list of (text, size, bold) for a leading
    Read me sheet. priority_column tints rows by a value starting P1/P2/P3; highlight_column
    tints non-empty cells. Returns path.
    """
    import openpyxl
    from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
    from openpyxl.utils import get_column_letter
    from openpyxl.worksheet.datavalidation import DataValidation
    head_fill, head_font = PatternFill("solid", fgColor="1F3864"), Font(bold=True, color="FFFFFF", size=10)
    tints = {"P1": PatternFill("solid", fgColor="F8CBAD"), "P2": PatternFill("solid", fgColor="FFE699"),
             "P3": PatternFill("solid", fgColor="E2EFDA")}
    mark = PatternFill("solid", fgColor="D9E1F2")
    edge = Border(*[Side("thin", color="BFBFBF")] * 4)
    wb = openpyxl.Workbook()
    del wb["Sheet"]
    if readme:
        ws = wb.create_sheet("Read me")
        ws.column_dimensions["A"].width = 118
        for text, size, bold in readme:
            ws.append([text])
            cell = ws.cell(ws.max_row, 1)
            cell.font = Font(size=size, bold=bold, color="1F3864" if (bold and size >= 12) else "000000")
            cell.alignment = Alignment(wrap_text=True, vertical="top")
    for spec in sheets:
        cols = spec["columns"]
        ws = wb.create_sheet(spec["name"][:31])
        ws.append(cols)
        widths = spec.get("widths") or [max(9, min(44, len(str(c)) + 6)) for c in cols]
        for i, _ in enumerate(cols, 1):
            cell = ws.cell(1, i)
            cell.fill, cell.font = head_fill, head_font
            cell.alignment = Alignment(wrap_text=True, vertical="center")
            ws.column_dimensions[get_column_letter(i)].width = widths[i - 1]
        ws.row_dimensions[1].height = 42
        for r in spec["rows"]:
            ws.append([r.get(c, "") for c in cols])
            n = ws.max_row
            for i in range(1, len(cols) + 1):
                ws.cell(n, i).border = edge
                ws.cell(n, i).alignment = Alignment(wrap_text=True, vertical="top")
            if priority_column and priority_column in cols:
                tint = tints.get(str(r.get(priority_column, ""))[:2])
                if tint:
                    for i in range(1, min(5, len(cols) + 1)):
                        ws.cell(n, i).fill = tint
            if highlight_column and highlight_column in cols and r.get(highlight_column):
                ws.cell(n, cols.index(highlight_column) + 1).fill = mark
        ws.freeze_panes = spec.get("freeze", "A2")
        ws.auto_filter.ref = "A1:%s%d" % (get_column_letter(len(cols)), ws.max_row)
        for col, opts in (spec.get("dropdowns") or {}).items():
            if col not in cols or ws.max_row < 2:
                continue
            letter = get_column_letter(cols.index(col) + 1)
            dv = DataValidation(type="list", formula1='"' + ",".join(opts) + '"', allow_blank=True)
            ws.add_data_validation(dv)
            dv.add("%s2:%s%d" % (letter, letter, ws.max_row))
    wb.save(path)
    return path
