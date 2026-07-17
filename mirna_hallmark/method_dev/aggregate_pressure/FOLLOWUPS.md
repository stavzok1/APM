# Method-dev follow-ups — aggregate pressure arc

Tracked open tasks from the force-vs-abundance / acquired-pressure arc.
Findings: `AGGREGATE_PRESSURE_FINDINGS.md`. Design log: `AGGREGATE_FORCE_VS_ABUNDANCE_DESIGN.md`.

## Priority — pursue

1. **Burden axis (THE key gap).** Coupling (correlation) is structurally **blind to constitutive
   repression** — a miRNA uniformly elevated in all tumors holds its target down everywhere but adds
   ~0 cross-sample variance → ρ≈0. Test a **level statistic**: per-gene *mean acquired pressure*
   `Σ_m mean_s[max(a_m−h_m,0)]·w_eff`, **no correlation**, vs the ± sets. Q: does burden separate
   positives (possibly better than coupling), and do burden+coupling together beat either alone?
2. **Constitutive vs modulated decomposition.** For the positives, split each gene's acquired pressure
   into constitutive (mean) vs modulated (cross-sample variance). If mostly constitutive ⇒ proves
   coupling is the wrong instrument and burden (1) is the right one.
3. **Within-sample compositional / CLR normalization** layered **before** the acquired transform.
   miRNA-seq is compositional; RPM doesn't remove hub-domination. Normalize each sample across the
   miRNA pool (CLR-style) to control per-sample global miRNA level (depth/DICER/proliferation), then
   apply `max(a−h,0)`. The only z-variant expected to do real work (across **all** miRNAs = a global
   normalization; across only g's regulators = the softmax/identity axis, do NOT use for detection).
4. **Protein (CPTAC) readout** for the positives. miRNA acts partly on translation; protein
   anti-correlation may exceed mRNA. CPTAC is available; use as a sensitivity on `acq`.
5. **Widen the positive set + TarBase v9.** n=18 limits power; TarBase v9 (6M validated, cell-type
   contexts) recovers the 4 dropped positives (HOXD10/DICER1/RECK/SMAD4) and expands coverage. See
   `[[mirna-edge-sources-landscape]]`. Ingest as its own evidence tier (HTP-CLIP vs LTP), test edge-universe sensitivity.
6. **CIBERSORTx composition sensitivity (Test 1, tumor-only).** Swap the marker metagenes for CIBERSORTx
   fractions (`annotations/cibersortx_brca_participant_features.tsv`; prior retest at
   `output/core_coupling_deconv_retest_cibersortx/`) — a more principled composition control. **Tumor-only**:
   for the symmetric Δρ it would require running CIBERSORTx on **GTEx** with the same signature (not done).

## Do NOT pursue (settled negatives / wrong axis)

- Edge weighting (validated evidence OR TargetScan context++) — **inert** (paired test).
- Promiscuity normalization (budget-split, any α>0) — **harmful**.
- Cross-sample z of acquired (centering inert; /σ marginal); softmax-for-detection (degenerate at aggregate,
  destroys magnitude — identity axis, keep for attribution only); within-regulator within-sample z (= softmax).
- Linear distance-from-healthy (un-rectified) — provably inert for a correlation (additive constant).
