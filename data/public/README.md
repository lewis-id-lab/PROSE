# Public cohort table

`profile_admissions.csv` is the only information from the patient cohort in this repository. It holds, for each of the 128
combinations of the seven binary predictors of the Stanzani 2019 score, the number of the 4,078 base-cohort admissions
with that combination.

| column | meaning |
|---|---|
| `profile` | 1 to 128; profile *i* has predictor *j* equal to bit *j* of *i* − 1 (`all_profiles()`, `R/baseline.R`) |
| the seven predictors | 0 or 1 |
| `admissions` | admissions with the profile; empty when suppressed |
| `suppressed` | TRUE for a profile with 1 to 4 admissions |

`profile_admissions_meta.csv` gives the total (4,078), the minimum cell size (5) and the number of suppressed profiles.

**Disclosure rule (25 November 2026).** No published cell has fewer than five admissions: 39 profiles are published
exactly (3,996 admissions), 36 are suppressed (82 admissions in all) and 53 are empty. No outcome, patient, date, ward or
transplant information is included. Event counts per profile are not published, because 55 of the 133 events sit in
profiles with fewer than five events.

**Use.** When the patient-level files are absent, `load_base_cohort()` rebuilds the cohort from this table, sharing the
suppressed admissions equally over the suppressed profiles (2 or 3 each). That moves 13 admissions (0.3% of the cohort)
away from their true profile. `reproduce.R` reruns every step that does not need outcomes or patient identity and reports
how far each output is from the committed one.

**Regenerating.** `scripts/56_export_public_cohort.R`, on the machine that holds the patient data. It refuses to write a
table that breaks the rule (`check_public_table()`, `R/public_cohort.R`).
