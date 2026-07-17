# Design — a new miRNA-pressure / decoupling prognostic (parsimonious + functional)

**Goal.** A *minimal* set of miRNAs or genes whose **pressure state** (potential, realized, or decoupled)
captures as much of the tumour's overall miRNA-repression biology as possible — and, ideally, prognoses.
Motivated by the MH-60→69 verdict: coarse pressure *magnitude* is not prognostic, but **realized
repression (anti-corr) > abundance**, **escape (decoupling) opens a TCGA crack**, and the signal is
**functional**, not level-based. So the construct must be built on the *functional/realized/acquired*
axes, with **target-aware weighting**, and selected for **parsimony + generalization**.

---

## 1. The three design axes (choose along each)

**(a) Unit** — what carries a feature:
- **miRNA-centric** — features = miRNAs; each summarises its action over *its target set*.
- **gene-centric** — features = target genes; each summarises the *pressure it receives* / its *escape*.
- Build **both** (they answer different questions: "which regulators" vs "which controlled nodes").

**(b) Readout axis** — what the feature measures (per patient):
| axis | miRNA-centric | gene-centric | rationale |
|---|---|---|---|
| **potential** | abundance × edges × AGO-gate (input) | Σ over miRNAs of input pressure | what *could* be repressed |
| **realized** | −Σ_g w(g)·z(expr_g) over targets (anti-corr) | −z(expr_g) given its pressure | what *is* repressed — **beats potential (MH-67)** |
| **decoupling / escape** | resid(realized ~ potential) | resid(expr ~ received-pressure) | dysfunction — **opened the only TCGA crack (MH-68)** |
| **acquired** | potential/realized(tumour) − (NAT or GTEx) | same | cancer-*specific* change; NAT-referenced = within-platform |

Lead with **realized + decoupling**; add **acquired** (now that matched NAT exists) as a cancer-specific
third axis. Potential/abundance kept only as the baseline to beat.

**(c) Target-aware aggregation** — the key lever the coarse analyses lacked. When a miRNA's readout is a
sum over its targets, weight each target g by `w(g) =` product of:
1. **edge evidence** (miRTarBase HE / TargetScan context score) — confidence the edge is real;
2. **target expressed-ness** — `z(expr_g)` floor: *you can only repress what is on*; un-expressed targets contribute noise;
3. **cancer-relevance** `r(g)` — DepMap `dep_role_weight` (dependency), driver status, malignancy_sign — repression of a *driver/dependency* matters more than a bystander;
4. **(optional) program membership** — restrict/weight to the program of interest.
- **Dropout** (sparsify the target set): keep the **top-K targets by `w(g)·|coupling|`**, drop the long tail of weak/noisy edges — a miRNA's prognostic signal may live in 2–5 targets (cf. miR-200→**LOX**, MH-69), so the *mean over all targets dilutes it*. Dropout is the fix for the MH-69 finding that role *subsetting* lost power (subsetting threw away genes; **weighting + top-K keeps all, emphasises the relevant**).

---

## 2. Minimal-set selection — "max information in fewest features"

Two objectives, two methods, then intersect:

**(i) Information-capture (unsupervised) — a faithful *summary* of the pressure landscape.**
- Build the per-patient feature matrix (chosen unit × axis, target-weighted).
- **mRMR** (max-relevance / min-redundancy) or **sparse-NMF / sparse-PCA** → the minimal panel that
  *reconstructs* the global pressure/decoupling structure (spans its independent axes). This is the
  "captures as much information about overall cancer miRNA pressure as possible" objective — valuable
  **independent of outcome** (a parsimonious read of the tumour's repression state).

**(ii) Prognostic (supervised) — generalising, not overfit.**
- **Sparse Cox** (LASSO / elastic-net) with **stability selection** (bootstrap; keep features chosen in
  ≥X % of resamples), **nested CV** for honest concordance, and **cross-cohort** (train TCGA → test
  Buffa-DRFS, and vice-versa). Only features that survive *across* cohorts count.

**(iii) Curated vs learned (offer both).**
- **Curated** (biology-first): the famous oncomiR/TSmiR set (MH-66) restricted to those whose targets are
  **drivers/dependencies**, on the realized/decoupling axis — a small interpretable panel.
- **Learned**: the stability-selected panel from (ii), intersected with the info-panel from (i).
- Compare: does learned beat curated? Does either beat clinical out-of-sample?

Target panel size: **~10–30 features** (events ≈100 → must stay parsimonious; rule-of-thumb ≤ events/10
free parameters, enforced by the L1 penalty + stability selection).

---

## 3. Concrete pipeline (buildable on existing machinery)

1. **Tensors** — per-(miRNA, patient) and per-(gene, patient) for {potential, realized, decoupling,
   acquired}, using `compute_gene_pressure_contributions` (per-edge), `load_gene_dependency` / roles for
   `r(g)`, RNA for expressed-ness, and matched-NAT (or GTEx/L3) for acquired. Target-weighted + top-K dropout.
2. **Feature pool** — miRNA-centric (~500 arms × 4 axes) + gene-centric (~2.5k genes × 3 axes), pre-filtered
   to expressed/edge-supported.
3. **Select** — (i) mRMR/sparse-NMF info-panel; (ii) stability-selected sparse Cox (nested CV, cross-cohort);
   intersect → the minimal panel(s).
4. **Score** — per-patient composite (panel's weighted realized/decoupling/acquired).
5. **Validate** — nested-CV concordance **vs clinical-only** (the honest gate); cross-cohort (TCGA↔Buffa);
   per-subtype where powered; **METABRIC-full** when the EGA miRNA lands (OS + n≈1,302 — the real test).
6. **Interpret** — which miRNAs/genes/axes/targets; tie to hypoxia/ECM (miR-210, miR-200→LOX) biology.

---

## 4. Honest expectations & guardrails

- **Prior is sober**: pressure magnitude was exhaustively non-prognostic; the realized/decoupling/acquired
  axes + target-relevance weighting + dropout are the levers most likely to *improve* on the coarse means,
  but signal is **Buffa/METABRIC-concentrated**, so a learned signature will likely be METABRIC-trained and
  must cross-validate to TCGA (thin). Overfitting risk is real (~100 events) → **stability selection +
  nested CV + cross-cohort are mandatory**, and "beats clinical out-of-sample" is the only acceptance bar.
- **The information-capture objective may be the more durable payoff** than prognosis: a parsimonious panel
  that *summarises the tumour's miRNA-repression state* (potential/realized/decoupled) is a useful descriptor
  of regulatory biology regardless of whether it predicts survival — and it is the framework-native deliverable.
- **Powered confirmation** of any prognostic claim = **METABRIC-full (EGA-pending)**; the coherence→OS lead
  (MH-65) is the existing thread this construct would formalise.

---

## Build-1 result (MH-73)
The realized/target-weighted+dropout construct **beats full clinical in METABRIC out-of-sample (+0.056 concordance, fair baseline)** — design principles validated — but **does not generalize to TCGA** (pressure < strong TNM+age clinical). Acquired (NAT/GTEx) no better than realized. **Build 1b (gene-centric mirror, MH-74) CONFIRMS unit-robustly:** realized beats full clinical in METABRIC (+0.060), decoupling +0.052, **magnitude does NOT (−0.016)** — the prognostic content is the *functional repression state*, not magnitude. TCGA is sparse for outcome (overfit) and cannot arbitrate. ⇒ **the construct is validated and ready; METABRIC-full (EGA) is the powered/cross-cohort confirmation.** **Builds 1c/1d (MH-75) COMPLETE the design:** curated 2-composite beats full clinical in METABRIC (+0.028); the repression landscape is **low-dimensional (20-miRNA panel = 90% variance)** and that unsupervised info-panel also beats clinical (+0.019). Design done; all functional variants beat clinical in METABRIC, magnitude never does. **BUT build 1e (MH-76): the gain does NOT generalize — a METABRIC-trained frozen panel adds ~nothing on TCGA (+0.002/+0.006); combine/RSF give no lift. The signal is real but METABRIC-SPECIFIC (likely platform: array vs seq).** A generalizable biomarker needs METABRIC-full + a same-platform external set — not achievable on current data.

## 5. Recommended first build (scoped)
Start **miRNA-centric, realized + decoupling axes, target-relevance-weighted + top-K dropout**, both
**curated** (driver-targeting TS/oncomiRs) and **learned** (stability-selected sparse Cox), TCGA+Buffa,
nested-CV vs clinical. If it clears the bar, add the **acquired (NAT-referenced)** axis and the
**gene-centric** mirror. One module (`pressure_prognostic_signature.py`) + a registry MH-id.
