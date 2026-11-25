"""Output schemas for the two agents, with allowed values taken from the bundle's codebook.csv."""
from common import read_bundle

NUM = {"anyOf": [{"type": "number"}, {"type": "null"}]}
INT = {"anyOf": [{"type": "integer"}, {"type": "null"}]}
STR = {"anyOf": [{"type": "string"}, {"type": "null"}]}


def codebook_enum(field, table="2_variants"):
    for r in read_bundle("codebook.csv"):
        if r["table"] == table and r["field"] == field:
            return [v.strip() for v in r["allowed values / type"].split("|")]
    raise KeyError(field)


def _obj(props, required=None):
    return {"type": "object", "properties": props,
            "required": list(required or props), "additionalProperties": False}


def extraction_schema():
    row = _obj({
        "listed_variant": {"type": "boolean"},
        "gene_symbol": {"type": "string"},
        "rsid": STR,
        "variant_name_as_reported": {"type": "string"},
        "effect_allele": STR,
        "other_allele": STR,
        "genetic_model": {"type": "string", "enum": codebook_enum("genetic_model")},
        "comparison_as_reported": {"type": "string"},
        "subgroup": {"type": "string"},
        "n_cases_hom_effect": INT, "n_cases_het": INT, "n_cases_hom_other": INT,
        "n_ctrl_hom_effect": INT, "n_ctrl_het": INT, "n_ctrl_hom_other": INT,
        "effect_measure": {"type": "string", "enum": codebook_enum("effect_measure")},
        "effect_estimate": NUM, "ci_lower": NUM, "ci_upper": NUM,
        "p_value": STR,
        "adjusted": {"type": "string", "enum": codebook_enum("adjusted")},
        "page_or_table": {"type": "string"},
        "verbatim_source": {"type": "string"},
        "confidence": {"type": "string", "enum": ["high", "medium", "low"]},
        "note": {"type": "string"},
    })
    missing = _obj({"variant": {"type": "string"}, "reason": {"type": "string"}})
    return _obj({"rows": {"type": "array", "items": row},
                 "not_found": {"type": "array", "items": missing},
                 "paper_note": {"type": "string"}})


def judge_schema():
    field = _obj({
        "field": {"type": "string"},
        "value_in_source": STR,
        "verdict": {"type": "string", "enum": ["A", "B", "both", "neither", "not_found"]},
        "quote": {"type": "string"},
        "page_or_table": {"type": "string"},
    })
    item = _obj({"item_id": {"type": "string"},
                 "fields": {"type": "array", "items": field},
                 "item_note": {"type": "string"}})
    return _obj({"items": {"type": "array", "items": item}})
