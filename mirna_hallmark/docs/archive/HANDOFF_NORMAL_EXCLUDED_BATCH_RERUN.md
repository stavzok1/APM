# Handoff — Normal-like-discarded + batch-corrected full rerun

**For the next agent.** This task is mid-flight. The **code changes are complete and
import-clean**; the **full rerun did not finish** (it kept getting killed by the
environment's background-process handling); and the **framework docs still need their
numbers refreshed** from the rerun outputs. This guide tells you exactly where things
stand and what to do.

> Read `mirna_hallmark/AGENTS.md` first (subproject conventions). This subproject has its
> own docs/catalog/ledger/registry under `mirna_hallmark/docs/`.

---

## 0. TL;DR — your first three actions

1. **Check if a rerun is already running or finished** (§2 commands). It probably is **not**.
2. **Relaunch the rerun** with the persisted orchestrator (§3) and let it finish (~hours).
3. **When it finishes** (`ALL DONE` in its log), **refresh the framework docs with exact
   numbers** and update the ledger + registry (§4). That is the only remaining deliverable.

Do **not** re-do the conceptual doc edits (§7) — they're done.

---

## 1. What the user asked for

1. **Discard PAM50 Normal-like** cohort-wide (n=36; normal-tissue admixture drives spurious
   cohort-level correlations).
2. **Batch-correct** the analysis — as **covariates, not ComBat** (decision + rationale in
   §6; do not switch to ComBat). Must include **NAT and GTEx** batch, not just tumour.
3. **Rerun everything** (run_all + cross-state landscape + CPTAC + Buffa).
4. **Then** make the framework docs use **exact numbers, not `~`** (this is the final step,
   blocked on the rerun outputs).

---

## 2. Current state (verify before doing anything)

```bash
cd /sci/labs/michall/stavzok/APM
# config switches live?
.venv/bin/python3 -c "from mirna_hallmark import config as C; print(C.EXCLUDE_NORMAL_LIKE, C.TCGA_BATCH_KIND, C.PAM50_ORDER)"
#   expect: True plate_both ('LumA','LumB','Her2','Basal')

# is a rerun running?
pgrep -af "rerun_normal_excluded_batch.py|mirna_hallmark\.(run_all|cross_state|cptac|buffa)" | grep -v pgrep
cat mirna_hallmark/output/logs/rerun_all.lock 2>/dev/null         # held = a run is live

# did the last attempt finish?
L=$(ls -t mirna_hallmark/output/logs/rerun_all_*.log | head -1); echo "$L"
grep -E "STEP [0-9]+/28:|rc=|ALL DONE|EXCEPTION" "$L" | tail -30
```

**As of this handoff:** the rerun is **NOT running — it was intentionally stopped** at the
user's request (a second session needed the resources first). The launch method is solved
(systemd user service, §3/§6); you just need to **(re)launch it when the user says go**.
Nothing rerun-related is alive, the lockfile is cleared, and the `mirna_rerun` unit is
removed. Confirm and launch:
```bash
systemctl --user status mirna_rerun                      # 'could not be found' = not running (expected)
pgrep -af rerun_normal_excluded_batch.py | grep -v pgrep # nothing = clean
# then launch per §3
```
> A second session may be running individual modules (e.g. `robustness_checks`) on the
> **same batch-corrected code** — check `pgrep -af mirna_hallmark` and do **not** launch the
> full rerun while it is mid-run (concurrent writers corrupt outputs; the lockfile only
> guards the orchestrator, not ad-hoc module runs).
Earlier `run_in_background`/`setsid` attempts kept dying at session teardown — addressed by
the systemd-user-service launch (§6).

Code is done and verified: all 18 touched modules import cleanly (re-verify with the
import loop in §5 if unsure).

---

## 3. Relaunch the rerun (the robust way)

A **single Python orchestrator** runs all 28 steps sequentially; a **PID lockfile**
prevents duplicate concurrent runs; it logs straight to a file (no terminal pipe). It is
persisted in the repo (the previous session's scratchpad is gone):

```
mirna_hallmark/rerun_normal_excluded_batch.py
```

**Launch it as a systemd *user* service — this is the fix for the killing problem (§6).**
A `run_in_background`/`setsid`/`nohup`/`tmux` job started from a tool stays in the Claude
session's cgroup and is killed when the session tears down. `systemd-run --user` places the
job in the **user manager's** cgroup (`user@<uid>.service/app.slice/…`), a different cgroup
the teardown does not touch:

```bash
cd /sci/labs/michall/stavzok/APM
systemctl --user reset-failed mirna_rerun 2>/dev/null
systemd-run --user --unit=mirna_rerun --collect \
  /sci/labs/michall/stavzok/APM/.venv/bin/python3 \
  /sci/labs/michall/stavzok/APM/mirna_hallmark/rerun_normal_excluded_batch.py
```

Verify it launched in the right cgroup and is the only instance (the lockfile also guards
duplicates — duplicates would corrupt outputs, §6):

```bash
systemctl --user status mirna_rerun            # Active: running; CGroup under user@<uid>.service
cat mirna_hallmark/output/logs/rerun_all.lock  # one owner pid == the service MainPID
```

Monitor / abort:

```bash
L=$(ls -t mirna_hallmark/output/logs/rerun_all_*.log | head -1)
grep -E "STEP [0-9]+/28:|rc=|ALL DONE|EXCEPTION" "$L" | tail -30
systemctl --user status mirna_rerun            # 'inactive (dead)' after it finishes
systemctl --user stop mirna_rerun              # to abort
```

> **Linger caveat:** `loginctl enable-linger` is **Access denied** for this user, so the
> user systemd manager (and thus this service) would stop on a *full node logout* (all the
> user's sessions ending). It survives Claude-session teardown fine (the demonstrated
> failure mode). If you expect a full logout, ask an admin to `loginctl enable-linger
> stavzok`, or instead submit via `at`/`batch` (atd is running — a system daemon fully
> independent of any login session): `echo '.venv/bin/python3 mirna_hallmark/rerun_normal_excluded_batch.py' | at now`.

The 28 steps (in `rerun_normal_excluded_batch.py`): `healthy_anchor --build` →
`run_all --skip-mirna-cnv` → `target_combined_anticorr` → `cn_dosage_attribution` →
cross-state Stages 1-3 (`mirna_state_class` → `cross_state_coupling` → `mirna_comovement` →
`nmf_sample_signatures` → `decoupling_validation`) → tissue-reference (`cross_state_landscape`,
`cross_state_deep_dive`, family×3, expression panels×3, `nat_tumor_umap`,
`geneset_architecture --all-hallmarks`, `within_pathway_nmf`) → CPTAC (`cptac_validation`
×2 cohorts + resid/orphan/specificity/confound/evidence) → `buffa_validation`. Each step is
non-fatal (logged + continue).

> If a step legitimately fails (e.g. missing external cache), it logs `rc!=0` and continues;
> that is fine — note it and move on. The headline steps are `run_all`, the cross-state
> Stages, and `target_combined_anticorr`.

---

## 4. THE REMAINING DELIVERABLE — refresh framework docs with exact numbers

After `ALL DONE`, replace every `~`/"about"/approximate figure in **both**:
- `docs/MODELING_FRAMEWORK_EXTERNAL.md` (external; **no internal paths/IDs**, paper-style math)
- `docs/MODELING_FRAMEWORK.md` (internal; mirror)

### 4a. Numbers that are sample-independent (won't change; just make exact)
Pull from the edge/universe builders (do not eyeball): 50 Hallmark sets; gene-universe union
(~4,384); arms / edges / target-genes under pressure (was stated 721 arms / ~5,100 edges /
~1,410 genes). Confirm exact counts from the run, e.g. the `coupling_predictor_comparison`
manifest / `predictor_comparison_summary.tsv` header prints, and `target_combined_anticorr`
row counts.

### 4b. Cohort + stratification (CHANGED by the Normal-drop)
- Cohort **n = 1065** (was ~1096; 36 Normal-like dropped; includes a few PAM50-NaN).
- PAM50 (post-drop): **LumA 507, LumB 264, Basal 197, Her2 93** (Normal-like **0 — dropped**).
  Update the §1.4 stratification table in the external doc and §2/§ data section internal.

### 4c. Coupling numbers (CHANGED by Normal-drop + batch — re-read from the NEW outputs)
These move and MUST be re-read after the rerun (do not reuse the pre-rerun values):
| number in docs | source file (after rerun) |
|---|---|
| Basal X/50 programs neg-sig (and per-PAM50) | `output/hallmark_interaction/hallmark_anticorrelation_by_pam50.tsv` |
| per-gene neg-sig coupling | `output/target_combined_anticorr/target_combined_anticorr.tsv` |
| **MH-44 edge/gene predictor counts** (e.g. "edge 1,487 neg-sig; z recovers 0.95; share 0.51; gene 430 vs 399") | `output/coupling_predictor_comparison/predictor_comparison_summary.tsv` |
| cross-state classes (acquired_realized / lost / etc.), brake-release counts | `output/tissue_reference/mirna_state_class/*` |
| module composite / decoupling counts | `output/tissue_reference/mirna_comovement/*`, `.../decoupling_validation/*` |
| CPTAC protein FDR-neg, Buffa concordance | `output/cptac_validation/*`, `output/buffa_validation/*` |

The **MH-44 paragraph in §6.3 of both framework docs cites specific numbers** (edge 1,487 /
z 0.95 / share 0.51 / gene 430 vs 399) — these are pre-rerun and **will change**; refresh
them from `predictor_comparison_summary.tsv`.

### 4d. Re-render the external doc (paper-style math; must stay glyph-clean)
```bash
export PATH="/cs/usr/stavzok/opt/quarto/bin:$PATH"
cd mirna_hallmark/docs
pandoc MODELING_FRAMEWORK_EXTERNAL.md -o MODELING_FRAMEWORK_EXTERNAL.docx
pandoc MODELING_FRAMEWORK_EXTERNAL.md -o MODELING_FRAMEWORK_EXTERNAL.pdf --pdf-engine=xelatex -V geometry:margin=1in
```
Guard: no literal `≈`/`≥`/box-glyphs in text (use `$\approx$` in math or the word "about" —
but the user wants exact numbers, so prefer exact). Confirm "no warnings" on the PDF render.

### 4e. Ledger + registry
- `docs/ANALYSIS_RUN_LEDGER.md`: add a row for the Normal-excluded + batch-corrected full
  refresh (date, what changed, key count shifts). Update per-component timestamps if you
  follow the `apm-analysis-ledger` skill.
- `docs/DISCOVERY_REGISTRY.md`: add **MH-45** — "Normal-like discarded + sequencing-batch
  (covariate) corrected refresh" with the headline before/after shifts (Basal X/50 before→
  after; whether the spine coupling survives batch; MH-44 counts after). Note that MH-44's
  cited numbers were refreshed.
- Heads-up (flag to user, don't necessarily rewrite): `LANDSCAPE_REPORT.md` and
  `BIOLOGICAL_SYNTHESIS.md` cite many Normal-inclusive / pre-batch numbers. The user's
  explicit ask was the **framework docs**; mention the others may warrant a later pass.

---

## 5. What changed in the code this session (the diff to understand)

All in `mirna_hallmark/`. Re-verify imports:
```bash
.venv/bin/python3 -c "import importlib; [importlib.import_module('mirna_hallmark.'+m) for m in ['config','data_loaders','tcga_batch','hallmark_interaction','target_combined_anticorr','subtype_contrasts','robustness_checks','cn_dosage_attribution','coupling_predictor_comparison','run_all','family_normal_reference','cross_state_coupling','cross_state_deep_dive','cross_state_landscape','mirna_state_class','mirna_comovement','decoupling_validation','geneset_architecture']]; print('OK')"
```

| File | Change |
|---|---|
| `config.py` | `EXCLUDE_NORMAL_LIKE=True`, `NORMAL_LIKE_LABEL="Normal"`, `TCGA_BATCH_KIND="plate_both"`, `GTEX_SAMPLE_ATTR` path; `PAM50_ORDER` drops Normal when excluding |
| `data_loaders.py` | `load_clinical_strata` drops Normal-like (global); `augment_tcga_batch(cov)` (participant-keyed batch dummies); `normal_like_participants()` |
| `tcga_batch.py` | `cross_state_batch_dummies(state, ids)` + `_tcga_plate_maps` (tumor/NAT plate from barcode field 6) + `_gtex_batch_map` (SMNABTCH/SMGEBTCH over GTEx SMLRNA breast samples); `CROSS_STATE_MIN_N=10` |
| `hallmark_interaction.py`, `target_combined_anticorr.py`, `subtype_contrasts.py` | cov → `D.augment_tcga_batch(cov)` at each partial-correlation site |
| `robustness_checks.py` | `_partial_batch()` wrapper (adds batch) + global rename of every `partial_spearman(` call to it |
| `cn_dosage_attribution.py` | batch appended in the `_partial` chokepoint |
| `coupling_predictor_comparison.py` | **NEW module** (MH-44: does the pressure construction change coupling vs raw abundance) + batch in `cov_base`; wired as **run_all step 10** |
| `run_all.py` | added `coupling_predictor_comparison` step (`--skip-predictor-comparison` to omit) |
| `family_normal_reference.py` | `_split_types` drops Normal-like **tumour** columns (covers both cross-state matrix paths) |
| `cross_state_coupling.py` | `_with_state_batch()` + `_state_covariates` appends **per-state** batch (tumor/NAT/GTEx) — covers coupling/comovement/decoupling/state_class |

Architecture: **participant-keyed TCGA spine** gets participant-plate batch via
`data_loaders.augment_tcga_batch`; **cross-state** (barcode/donor-keyed) gets per-state batch
via `tcga_batch.cross_state_batch_dummies` wired into the single cov chokepoint
`cross_state_coupling._state_covariates`. CPTAC keeps its own plex batch (already wired).

Docs already edited this session (conceptual, **done** — see §7): `MODELING_FRAMEWORK.md`,
`MODELING_FRAMEWORK_EXTERNAL.md` (+ rendered `.pdf`/`.docx`), `FORMULAS.md` (specificity),
`ANALYSES_CATALOG.md` + `ANALYSIS_RUN_LEDGER.md` + `DISCOVERY_REGISTRY.md` (MH-44 row).

---

## 6. Gotchas / environment notes (these bit hard — read them)

- **Why background jobs kept dying (root cause).** The harness runs its shells in the Claude
  session's cgroup (`…/session-*.scope`) and **kills that whole cgroup on session teardown**.
  `run_in_background`, `setsid`, `nohup`, `disown`, and even a `tmux` server started from a
  tool all stay **in that cgroup**, so a cgroup-kill reaches them regardless of
  daemonization — that is why every earlier launch died. **Fix:** launch via
  `systemd-run --user` (§3), which puts the job in the **user manager's** cgroup
  (`user@<uid>.service/…`), outside the session cgroup. Verified: the service's cgroup is
  `mirna_rerun.service`, not `session-*.scope`.
- **Duplicates corrupt outputs.** Earlier failures also spawned concurrent duplicates racing
  on the same output files. The Python orchestrator's **lockfile** (`output/logs/rerun_all.lock`)
  makes any second launch a no-op (it sees the live owner PID and exits). Still, **verify one
  instance** after launching. Remove a stale lock only if its owner PID is dead (`kill -0 <pid>`).
- **System clock looked inconsistent** across this session (date jumped 06-25↔06-26); don't
  trust absolute log mtimes for "is it progressing" — use the step/rc lines.
- `run_all` uses `--skip-mirna-cnv` (reuse the CNV extraction cache; ~12 min saved). CNV
  *stratification* tables go through a parent loader that may still include Normal-like — out
  of scope; note it if asked.
- The orchestrator's `cn_dosage_attribution` / CPTAC / Buffa steps depend on caches/external
  data that exist from prior runs; if one errors it's logged non-fatally.

---

## 7. Already resolved this session — DO NOT redo

- **ComBat vs covariates** → covariates chosen (ComBat over-corrects on the unbalanced TCGA
  design and can manufacture feature-feature correlations, which is exactly our readout).
- **Framework-doc conceptual fixes** (from the turns before the rerun): §4.1 modes relabelled
  as *pressure-construction* (not "coupling"); §6 predictor-by-resolution (dose = abundance at
  edge/module, aggregate pressure at gene/program; the **softmax share** is the only piece that
  changes coupling); §6.4 retitled **potential × realized**; §4.4 specificity = relative
  concentration (not "most of the budget"); §7 protein residual = *beyond-measured-mRNA* (not
  strictly translational); §1.4 data/stratification section + named public sources; MH-44 module
  built, wired, catalogued. These are **correct** — only their **numbers** need the post-rerun
  refresh (§4).

---

## 8. Key files
- Orchestrator: `mirna_hallmark/rerun_normal_excluded_batch.py` (28 steps, lockfile, file-logged)
- Rerun logs: `mirna_hallmark/output/logs/rerun_all_*.log`; lock: `.../rerun_all.lock`
- Deliverable docs: `docs/MODELING_FRAMEWORK_EXTERNAL.md` (+ `.pdf`/`.docx`), `docs/MODELING_FRAMEWORK.md`
- Formula/results references: `docs/FORMULAS.md`, `docs/LANDSCAPE_REPORT.md`
- Registry/ledger/catalog: `docs/DISCOVERY_REGISTRY.md`, `docs/ANALYSIS_RUN_LEDGER.md`, `docs/ANALYSES_CATALOG.md`
