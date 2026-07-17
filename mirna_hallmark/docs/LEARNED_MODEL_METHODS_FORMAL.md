---
title: "Learned model — METHODS"
subtitle: "Formula-level spec of every estimator"
---

The formal math for each estimator in `mirna_hallmark/learned/`. This document is the **spec** —
notation, objective, formula, and one line of "what it computes" — kept clean of prose. **Which
estimator for which job, and why over the alternatives → §18 (in this doc).** Companions:
`LEARNED_MODEL_RATIONALE.md` (the *why* behind every section here — section numbers mirror this
doc), `LEARNED_MODEL_DISCOVERY_SYNTHESIS.md` §6b (the converged model in one line), and
`LEARNED_MODEL_VALIDATION.md` (does it work). The precursor **pressure** model's math is separate,
in `METHODS.md`.

# §0 — Notation & identification

Per gene $g$, aligned on tumour participants (PAM50 Normal-like dropped):

- $Y$ — target mRNA, $\log_2(\text{TPM}+1)$, vector over participants.
- $X$ — regulator-arm abundance, $\log_2(\text{RPM}+1)$, participant $\times$ arm (undetected arm $\to 0$).
- $C$ — confounder block (participant $\times$ conf), see §1.
- $w$ — adaptive-penalty prior weight per arm, see §2.

**Candidate regulators** (inclusion, migrated 2026-07-06): $\text{POOLED-HE} = \text{miRTarBase-HE}\ \cup\ \text{TarBase-v9}$
low-throughput functional (reporter / western / proteomics, non-weak) — `ledger.pooled_he_edges()`,
$\approx 5{,}940$ edges over $1{,}571$ Hallmark genes; `assemble_gene(he_only=True)` returns this.
Replaced miRTarBase-HE-only; gate-stable on anchors ($\Delta\rho \le 0.003$).

**Residualisation** (matrix-view identification, Design §2/§B): for a vector $v$ and matrix $C$,

$$\operatorname{resid}(v \mid C) \;=\; v - C\,(C^{+} v)$$

(OLS fit removed; intercept included). Every coupling / weight below is estimated in the residual
space $(I - P_C)$ — i.e. $C$ is partialled out of **both** the predictor and the target.

# §1 — Confounder block $C$

*Modules: `data.assemble_gene`.*

Core (always): $C = [\,\text{CPE},\ \text{target\_cn},\ \text{mal\_prolif}\,]$.

- $\text{CPE}$ — consensus tumour purity (per participant).
- $\text{target\_cn}$ — gene $g$'s own ASCAT3 copy number (`_target_cn`), the gene's genomic dosage.
- $\text{mal\_prolif}$ (`_malignant_prolif`) — malignant-compartment proliferation:

$$\text{score} = \operatorname*{mean}_{g' \in (\text{E2F} \cup \text{G2M})} z_{\text{row}}(Y_{g'}),
\qquad
\text{mal\_prolif} = \operatorname{resid}\!\big(\text{score} \mid \text{CancerEpithelial\_fraction}\big)$$

(removes the composition-driven part, leaving per-malignant-cell proliferation).

Optional add-ons (joined as extra $C$ columns):

- `deconv` — the 8 non-malignant CIBERSORTx fractions $[\text{CAFs},\text{T},\text{Myeloid},\text{B},\text{Endothelial},\text{PVL},\text{NormalEpi},\text{Plasmablasts}]$ (Cancer-Epithelial excluded — conditioning on the target's own compartment over-controls).
- `mycaf` (`_mycaf`) $= \operatorname{resid}\big(\operatorname{mean} z(\text{ACTA2,TAGLN,POSTN,FAP,COL11A1,THBS2,MYH11,CNN1}) \mid \text{CAF\_fraction}\big)$ — myCAF-vs-iCAF axis.
- `mal_emt` (`_mal_emt`) $= \operatorname{resid}\big(\operatorname{mean} z(\text{VIM,FN1,CDH2,SNAI2,ZEB1,ZEB2,TWIST1,SPARC}) \mid \text{non-malignant fractions}\big)$ — within-cancer EMT state.
- HRD off; naive batch off (post-fit conditioning only, never in-fit).

**Defaults (which $C$ is actually used).** The core $[\text{CPE},\text{target\_cn},\text{mal\_prolif}]$
is **always on**. The non-malignant CIBERSORTx fractions (`deconv`) are **opt-in**: OFF for the
gate / coupling (`assemble_gene` `deconv=False`), but **forced ON for the canonical attribution**
(`canonical_M` / `budget_shift` / `decompose` call `deconv=True`) — those weights are *cell-intrinsic*.
The metagene-$z$ composition (`_state_metagene_cov`, a per-program $z$-metagene) is only the
cross-state fallback for GTEx (no CIBERSORTx there); NAT / tumour use the CIBERSORTx fractions.
So: tumour attribution = fractions; gate / coupling = core only unless `deconv` requested; healthy =
metagene-$z$ until the pending GTEx CIBERSORTx run lands.

**Post-fit conditioning.** `batch` (TCGA plate dummies) is **never in the lasso** — an in-fit
ComBat-style correction removes biological heterogeneity along with the batch — so batch is appended
to $C$ **only** for the partial-Spearman $\rho$ significance computation (`oof_gate`), conditioning
the *test*, not the *model*.

**Within-cancer vs non-malignant.** $\text{mal\_prolif}$ residualises the proliferation metagene on
the Cancer-Epithelial fraction $\to$ per-malignant-cell proliferation; the `deconv` block adds the
non-malignant fractions (CAFs / immune / …) to strip stromal / immune composition. Cancer-Epithelial
is deliberately excluded from `deconv` — conditioning on the compartment the target is expressed in
over-controls the signal.

# §2 — The prior $w$

*Modules: `evidence/ledger.py`, `evidence/methods.py`.*

Evidence enters as **magnitude**, non-circularly. Per (edge, assay class) count distinct PMIDs; the
fused edge weight is

$$w_m \;=\; \sum_{\text{class}} w_{\text{class}} \cdot \log\!\big(1 + \#\,\text{distinct\_PMID}(m,\text{class})\big).$$

- $w_{\text{class}}$ (assay directness) is either **hand-set** (`CLASS_WEIGHT`: reporter $3 >$ western/proteomics $2.5 >$ qpcr $1.5 >$ degradome/chimeric $1 >$ ago-clip $0.5 >$ other $0$) or **learned** (`CLASS_WEIGHT_MRNA`): non-negative least squares of measured transfection repression (McGeary GSE140217) on per-class $\log$-PMID, within-family-centred.
- $w = \operatorname{clip}(w_m,\ \ge 10^{-3})$. Alternatives: `scanmir` ($K_D$ affinity), `fused` (existence $\times$ magnitude).

**Default.** The `assemble_gene` / `oof_gate` signature default is `evidence_score` (hand-weighted),
but the canonical / reported fits force `ledger` (the weighted, PMID-deduped, method-centric fusion)
— so we use **weighted evidence, not binary inclusion**. Empirically the source is non-critical (hub
panel, family gate): mean $\rho_{\text{model}}$ is ledger $-0.345 \approx$ evidence $-0.340 \approx$
fused $-0.342 \approx$ scanmir $-0.331$; stability $0.86/0.86/0.85/0.84$; ledger passes $5/5$ gate,
scanmir/fused $4/5$ (scanmir breaks GATA3, fused breaks PTEN vs-curated). $\Rightarrow$ ledger is the
marginal best and is canonical; the prior sets *ordering*, the data sets *magnitude* (§3).

# §3 — The model: gene-focused non-negative adaptive lasso

*Modules: `regression.fit_gene`.*

Objective (one gene), on $C$-residualised $\tilde{y}=\operatorname{resid}(Y\mid C)$, $\tilde{X}=\operatorname{resid}(X\mid C)$:

$$\hat{M} \;=\; \operatorname*{arg\,min}_{M \ge 0}\ \big\lVert\, \tilde{y} - \tilde{X} M \,\big\rVert^2 \;+\; \alpha \sum_m \frac{\lvert M_m \rvert}{w_m}.$$

- Non-negativity encodes "miRNAs repress" (Design §Decision C).
- **Adaptive** penalty $\lambda_m = \alpha / w_m$: strong-evidence edges resist shrinkage, weak edges shrink out. The prior is an **inverse penalty** (resistance to being zeroed), **not** a location — it sets ordering / selection, not magnitude.
- Implementation (Zou reparametrisation): scale column $m$ by $w_m \to$ plain non-negative lasso (uniform $\alpha$) $\to$ unscale $M_m = b_m w_m$. Default $\alpha = 0.005$ (fixed; non-critical — global grid $0.001$–$0.02$ gives flat mean $\rho_{\text{model}} \approx -0.340$, gate $8/8$; $\alpha$ is a **sparsity** knob [nonzero $9 \to 7$], not a coupling knob; only $0.05$ over-shrinks). $\alpha$ touches **only** the lasso — the canonical attribution (§5) is unpenalized bagged NNLS.
- Aggregate repression pressure: $\rho_{\text{agg}} = X M$ (higher $\Rightarrow$ more repression $\Rightarrow$ lower $Y$).
- This is the model used for coupling / selection / validation (§6) and discovery. Contribution $\approx a \cdot M$ (unsaturated linear occupancy; the saturating $a/(a+K)$ link was tried and failed the gate).

# §3b — Observation likelihood: Gaussian vs Student-t &nbsp;(`spike_slab._gibbs_posterior(nu=)`, MH-92, 2026-07-10)

The dense §6b posterior fits $r = X\beta + \varepsilon$ with default **Gaussian** $\varepsilon_i \sim \mathcal{N}(0,\sigma^2)$ (voom-style, on log-CPM residuals). A residual diagnostic (85 genes $\times$ 1041 patients) found heavy tails — pooled excess-kurtosis $3.5$ (median $1.2$ across well-expressed genes), Student-t MLE $\nu \approx 7$, $P(|z|>4)=33\times$ Gaussian — from amplified-subset tumours ($+$skew) and a few near-floor genes ($-$skew). **NB is rejected** ($\rho_{\text{Spearman}}(\text{expr},\text{kurtosis})=+0.19$, wrong sign for a low-count structure). The upgrade is **Student-t** via a scale-mixture of normals (Gibbs-conjugate, no HMC):

$$\varepsilon_i \sim \mathcal{N}\!\big(0,\ \sigma^2/\lambda_i\big),\qquad \lambda_i \sim \text{Gamma}(\nu/2,\ \nu/2)\ \Longrightarrow\ \varepsilon_i \sim t_\nu(0,\sigma^2)$$
$$\lambda_i \mid \varepsilon_i,\sigma^2 \sim \text{Gamma}\!\Big(\tfrac{\nu+1}{2},\ \tfrac{\nu + \varepsilon_i^2/\sigma^2}{2}\Big)$$

The per-coordinate conditional gains weights $A_m = \sum_i \lambda_i x_{im}^2/\sigma^2$, $B_m = \sum_i \lambda_i x_{im} r_{m,i}/\sigma^2$, and $\sigma^2 \sim \text{InvGamma}(a_0+n/2,\ b_0+\tfrac12\sum_i \lambda_i \varepsilon_i^2)$; a large-residual patient gets small $\lambda_i$ (down-weighted). **$\nu=\text{None} \Rightarrow$ Gaussian, bit-identical to §6b.** Verdict: correctness/robustness upgrade, not a performance lever at $n\approx1000$ — nests Gaussian ($\Delta\beta\approx0.001$), robust $1.76\times$ under $5\%$ contamination, OOF coupling $\Delta\rho=+0.001$ (n.s.), attribution unchanged (corr $0.999$). The same mixture applies to the CN channel term ($t$ on the $\hat\pi_s$ pseudo-obs, `t['nu']`; $s^2_\pi$ estimated $\Rightarrow$ $t$ guards collider-outlier segments).

# §4 — Seed-family collapse

*Modules: `families.collapse_by_family`, Design §F.*

Same-seed arms are near-collinear $\Rightarrow$ the identified estimand is family $\to$ gene (arm =
nomination).

$$X_{\text{fam}} = \log_2\!\Big(1 + \sum_{m \in \text{family}} (2^{x_m} - 1)\Big) \quad\text{(true RPM pool, re-logged)},
\qquad
w_{\text{fam}} = \max_{m \in \text{family}} w_m.$$

# §5 — Canonical attribution: bagged family NNLS

> ✅ **SETTLED 2026-07-17 (MH-138) — no conflict; one word, two jobs.**
> **Magnitude** $=$ the dense Gibbs posterior mean $\beta$ ⇒ **§6b**. **Identity** $=$
> `attribution.shapley_identity` $=$ LMG on $R^2$, which needs **fixed weights** $m$, supplied by
> $\texttt{canonical\_M}$ = **the bagged NNLS below** ⇒ **§5, for that job**. Code already coherent
> (`identity_vs_magnitude:131`, `bayes_shapley_identity:91` both call $\texttt{ST.canonical\_M}$) and matches
> §2e's deployable form: **NNLS fixes support + signs; Gibbs draws give the width**.
> **⭐ Why identity cannot take the Gibbs mean (§2e's mechanism, MEASURED, $n=5{,}117$):** the half-normal slab
> $\mathcal{N}^{+}(0,\nu^{2})$ has mean $\nu\sqrt{2/\pi} > 0$, so an un-informed family **cannot be zeroed** —
> it relaxes to the prior. $\beta = 0$ in $\mathbf{0/5{,}117}$; all $\beta>0$, where NNLS returns exact zeros.
> **Pooled $\mathbf{73.0\%}$ of $\beta$-mass sits on unidentified edges** ($\lvert z\rvert \le 2$); per-gene
> median $55.5\%$; $42.7\%$ of genes at $100\%$ (452 single-family, where $\beta \equiv$ uniform); $n_{\text{fam}}\ge3$
> median $32.0\%$.
> Also settled: **bagged NNLS is not retired** — correct for the **cross-cohort gauge** ($a=4.1$ vs split-half
> truth $1.0$ under Gibbs's heavy-tailed SDs) ⇒ ***NNLS for the GAUGE, Gibbs for the MODEL***; on the model
> Gibbs reproduces better ($\rho$ $0.822$ vs $0.729$). ⚠ The single-*lasso* premise below is pre-convergence
> framing; the collinearity instability it describes is real regardless. See §18 and MH-138.

*Modules: `states._bagged_nnls`, `canonical_M`.*

The stable readout of per-edge weight (single-lasso coefficients are unstable under collinearity,
corr $0.03$). On $C$-residualised, $z$-scored family predictors:

$$
\begin{aligned}
y_r &= -\operatorname{resid}(Y \mid C) &&\text{(sign: repression} \to +)\\
X_z &= \operatorname{zscore\_cols}(X_{\text{fam}}) &&= (X_{\text{fam}} - \operatorname{mean})/\operatorname{sd}\\
X_z[:,\ \operatorname{sd} < 0.1] &= 0 &&\text{(variance floor)}\\
\text{for } b = 1..B:\ \ \text{idx}_b &= \operatorname{bootstrap}(n),\ \ c_b = \operatorname{NNLS}(X_z[\text{idx}_b],\, y_r[\text{idx}_b])\\
M_{\text{fam}} &= \tfrac{1}{B}\textstyle\sum_b c_b &&\text{(bagged, } B=40)
\end{aligned}
$$

with the sign chosen so repression maps to a positive weight, a **variance floor** that zeroes any
arm undetectable in the state ($\operatorname{sd}<0.1$), and $B=40$ bootstrap replicates whose
bagging drives coefficient variance down by $\sim 1/B$.

- Reproducible (corr $\sim 0.99$ across seeds) and recovers curated drivers (PTEN miR-21/103a/181/182/96). $z$-scored because level-scale NNLS is ill-conditioned (recovers noise arms).
- **arm_level broadcast:** $M_{\text{arm}} = M_{\text{fam}}(\text{family}(\text{arm}))$ — family weight shared to member arms (per-arm split = nomination). Used everywhere coefficients are read: budget, $\Delta M$/wiring, decompose, card.
- **$\Delta M$ / cross-state wiring:** fit $M_{\text{fam}}$ separately per state (tumour "01" / NAT "11" / GTEx) on the **same fixed** family support $\to \Delta M = M_{\text{state1}} - M_{\text{state2}}$ (fixed support $\Rightarrow$ comparable estimand).

## §5a — Budget contribution (family-size-correct) &nbsp;(`states.budget_shift`)

Realized within-gene budget, apportioned so a multi-arm family isn't over-counted:

$$\operatorname{contribution}(\text{arm},\text{state}) = M_{\text{fam}} \cdot X_{\text{fam}}(\text{state}) \cdot \frac{2^{a_{\text{arm}}(\text{state})}}{\sum_{m \in \text{fam}} 2^{a_m(\text{state})}},$$

$$\operatorname{share}(\text{arm}) = \frac{\operatorname{contribution}}{\sum_{\text{arms}} \operatorname{contribution}}, \qquad \operatorname{rank} = \text{descending rank of contribution},$$

so $\sum_{\text{arm}\in\text{fam}} \operatorname{contribution} = M_{\text{fam}} \cdot X_{\text{fam}}$ (no family-size inflation). A singleton's $X_{\text{fam}}$ = its own abundance, so contribution $= M a$ (unchanged).

# §6 — Validation gate: out-of-fold aggregate coupling

*Modules: `mvp.oof_gate`, `mvp.gate_fdr`.*

5-fold over participants; fit $M$ on train, score held-out:

$$
\begin{aligned}
\text{oof\_model}[te] &= X[te]\, M_{\text{train}}, \qquad \text{oof\_abund}[te] = \operatorname*{mean}_m X[te] && \text{(per participant)}\\
\rho_{\text{model}} &= \operatorname{spearman}\!\big(\operatorname{resid}(\text{oof\_model}\mid C),\, \operatorname{resid}(Y \mid C)\big) && \text{(expect} < 0)\\
\rho_{\text{abund}} &= \operatorname{spearman}\!\big(\operatorname{resid}(\text{oof\_abund}\mid C),\, \operatorname{resid}(Y \mid C)\big) && \text{Bar-1: equal-weight inclusion}\\
\rho_{\text{curated}} &= \operatorname{spearman}\!\big(\operatorname{resid}(X w_{\text{prior}}\mid C),\, \operatorname{resid}(Y \mid C)\big) && \text{Bar-2: frozen evidence weights (no fit)}
\end{aligned}
$$

Gate passes if $\rho_{\text{model}} < \rho_{\text{abund}}$ (learning beats abundance) **and**
$\rho_{\text{model}} \le \rho_{\text{curated}}$ ($\ge$ frozen prior). Genome-wide FDR (`gate_fdr`): per
gene one-sided partial-$t$ $p$ for $\rho_{\text{model}}<0$ on $\mathrm{df}\approx n-8$, then BH and
Benjamini–Yekutieli across genes (BY for shared-regulator dependence). Only meaningful for
$\ge 2$-regulator genes (singletons: $X M \equiv$ abundance up to scaling).

# §7 — Per-edge coupling: family-aware permutation partial-Spearman

*Modules: `coupling.edge_coupling` → `coupling_inference.permutation_pvalues`.*

Per gene, predictor rows = arms, target = gene mRNA (repeated), covariates = $C$, families = seed
family:

$$\operatorname{obs}_m = \sum_s \operatorname{unit}\!\big(\operatorname{resid}(\operatorname{rank}(x_m)\mid \operatorname{rank}(C))\big)_s \cdot \operatorname{unit}\!\big(\operatorname{resid}(\operatorname{rank}(Y)\mid \operatorname{rank}(C))\big)_s \quad \text{(partial Spearman)}$$

- **null:** the same permutation $\pi$ of the predictor-residual sample order applied to **every** row (preserves family dependence).

$$p_{\text{perm}} = \frac{1 + \#\{\text{null} \le \operatorname{obs}\}}{n_{\text{perm}}+1}, \qquad p_z = \Phi\!\Big(\frac{\operatorname{obs} - \mu_{\text{null}}}{\sigma_{\text{null}}}\Big) \quad \text{(smooth tail from null moments)}$$

$$q_{\text{bh}} = \operatorname{BH}(p_z), \quad q_{\text{by}} = \operatorname{BY}(p_z), \quad q_{\text{family}} = \text{Simes-per-family } \operatorname{BH}(p_z) \quad \text{(§14)}$$

Called **per gene** ($C$ carries per-gene $\text{target\_cn}$). Genome-wide: second-level BY across
all edges.

# §8 — Causality: CN-locus instrument

> ⛔ **THERE IS NO LIVE "CAUSALITY" JOB — BOTH CN INSTRUMENTS ARE RETRACTED (banner added 2026-07-17).**
> The machinery below still runs; **its causal interpretation is dead**, and the retraction is in the production
> code (`learned/instrument.py:1128–1148`).
> **(1) $\pi_{\text{causal}}$ is NOT an IV** (MH-124r/126): it is $\gamma_s\cdot b_{\text{fam}}$, a
> **product-of-coefficients mediation** estimator whose $b$ is the **observational OLS slope** (the
> instrument-orthogonal, endogenous variation) — so it re-derives the anti-correlation it was meant to validate.
> A 2SLS ratio would be $\pi/\gamma$, not $\gamma\cdot b$. $\gamma_s$ is **site-blind** (HE $+0.19931$ vs decoy
> $+0.19844$, $p=0.20$); the clean reduced form **clustered at the arm** gives $p=0.115$ (n.s.).
> **(2) The within-arm two-way-FE replacement is also retracted** (MH-133): control class **71% real binders**;
> against a genuinely site-free control $\tau=-0.0007$, $p=0.84$; fails the site-type efficacy ladder
> ($p=0.26$). **Refutations at $n=216$k–$235$k pairs, not power failures.**
> ⇒ **Edge existence rests on ONE observational line.** The $F>10 \wedge \text{T1-clean}$ rule remains a
> **validity** rule (exclusion fails for ~73% of well-instrumented edges), never a reality test.
> Current state: [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) Axis 2.

*Modules: `instrument.edge_instrument`, `run_clean`.*

miRNA-locus copy number $\text{CN}(m)$ instruments arm dose:

$$
\begin{aligned}
\text{first stage:}\quad & x_m = \gamma\,\text{CN}(m) + C\delta + u, \qquad F = \frac{(\text{RSS}_{\text{reduced}} - \text{RSS}_{\text{full}})/1}{\text{RSS}_{\text{full}}/\mathrm{df}} \quad \text{usable if } F > 10\\
\text{reduced form:}\quad & \rho_{\text{reduced}} = \operatorname{spearman}\!\big(\operatorname{resid}(\text{CN}\mid C),\, \operatorname{resid}(Y\mid C)\big) \quad \text{expect} < 0 \text{ (genetic dose} \downarrow \text{target)}\\
\text{observational:}\quad & \rho_{\text{obs}} = \operatorname{spearman}\!\big(\operatorname{resid}(x_m\mid C),\, \operatorname{resid}(Y\mid C)\big)\\
& \text{causal\_concordant} = (\rho_{\text{reduced}} < 0) \wedge (\rho_{\text{obs}} < 0) \wedge (F > 10)\\
& p_{\text{reduced}} = \text{one-sided Spearman-}p(\rho_{\text{reduced}}<0), \quad q_{\text{reduced,by}} = \operatorname{BY}(p_{\text{reduced}})
\end{aligned}
$$

**Cluster-exclusion** (`cluster_attribution`): CN moves the whole genomic cluster. For co-located
miRNAs that also target $g$ (set $S = \{m\} \cup \text{co-targeters}$), condition each on the others $+ C$:

$$\operatorname{part}(k) = \operatorname{spearman}\!\big(\operatorname{resid}(x_k \mid \{x_{j\neq k}\} \cup C),\ \operatorname{resid}(Y \mid \{x_{j\neq k}\} \cup C)\big).$$

$\text{arm\_unique} = \big(\operatorname{part}(m) < -0.1\big) \wedge (\text{no co-targeter survives})$;
$\ \text{strong\_causal} = \text{concordant} \wedge (\text{cluster\_clean} \vee \text{arm\_unique})$;
$\ \text{strong\_causal\_fdr}$ adds $q_{\text{reduced,by}} < 0.05$.

# §9 — Driver resolution inside a family

*Modules: `families.resolution_report`.*

For a model-selected multi-member family, each member's anti-corr is measured *net of its
family-mates* $(+C)$:

$$\operatorname{part}(\text{arm}) = \operatorname{spearman}\!\big(\operatorname{resid}(x_{\text{arm}} \mid \{\text{mates}\}\cup C),\ \operatorname{resid}(Y \mid \{\text{mates}\}\cup C)\big).$$

A member **survives** if $\operatorname{part} < -0.1$ (real conditioned repression). The nomination
reflects how the surviving anti-corr distributes (not a forced single winner):

- **0 survivors $\to$ "family-only"** — the split is unidentified (collinear members lose their residual variance $\Rightarrow \operatorname{part}\to 0$, auto-fail; e.g. ESR1 miR-221/222).
- **1 survivor $\to$ "single-driver"** (e.g. PTEN miR-181b, miR-429).
- **$\ge 2$ survivors $\to$ "co-drivers"** — the repression is carried by several members, reported as a set (e.g. PTEN miR-17$\sim$92 $=$ miR-106b $-0.22$ $+$ miR-20a $-0.13$; ZEB1 $=$ miR-429 $+$ miR-200c).

Columns: `resolution` (the class), `n_drivers`, `drivers` (the surviving set, strongest first),
`driver` (strongest survivor, compat), `member_partials` (every member's part). Empirically
$\approx 13\!:\!7\!:\!11$ single : co : family-only.

# §10 — Uncertainty: hierarchical Bayes

*Modules: `hierarchical._gibbs`.*

Per program (genes sharing regulators), on $C$-residualised $z$-scored $X$:

$$
\begin{aligned}
y_g &= X_g \beta_g + \varepsilon_g, & \varepsilon_g &\sim \mathcal{N}(0, \sigma_g^2)\\
\beta(m,g) &\sim \mathcal{N}(\mu_m, \tau^2) & &\text{(miRNA-level pooling across the genes } m \text{ targets)}\\
\mu_m &\sim \mathcal{N}(0, \omega^2), & \tau^2, \sigma_g^2 &\sim \text{InvGamma}
\end{aligned}
$$

Gibbs (conjugate): $\beta_g \mid \cdot \sim \mathcal{N}\!\big(V(X^\top X/\sigma^2 + \mu/\tau^2),\, V\big)$,
$V = (X^\top X/\sigma^2 + I/\tau^2)^{-1}$; $\mu_m \mid \cdot$ pooled across its genes;
$\tau^2, \sigma_g^2 \mid \cdot$ inverse-gamma. Output: posterior mean $\pm$ sd of $\beta$; **identified
edge** $= \lvert \beta/\text{sd}\rvert > 2$. (Cross-gene magnitude pooling does not improve OOF
coupling here — used as the uncertainty object; coupling stays per-gene, §6.)

# §11 — Subtype tests

*Modules: `subtype.py`.*

- **Level** (Kruskal): $\operatorname{kruskal\_across\_strata}(X M,\ \text{PAM50})$ — does pressure magnitude differ across subtypes.
- **Slope** (interaction): nested rank-OLS $F$-test $\operatorname{rank}(Y) \sim \operatorname{rank}(\rho_{\text{agg}})\cdot\text{subtype} + C$ vs common-slope — does the coupling differ; BH across genes.
- **Wiring, independent** (`subtype_wiring`): fit canonical $M_{\text{fam}}$ (§5) separately per subtype; mean pairwise cross-subtype $\operatorname{spearman}(M)$ vs an $n$-matched noise floor (random cohort subsamples). High variance (differences of two small-$n$ fits) — *over-calls* remodeling.
- **Wiring, pooled** (Decision H, `subtype_wiring_pooled`): whole-cohort $M_{\text{all}}$ (§5) is the **prior mean**; per subtype $M_s = M_{\text{all}} + \delta_s$ via ridge-toward-$M_{\text{all}}$,

$$\delta_s = \operatorname*{arg\,min}_{\delta}\ \big\lVert (y_r - X_z M_{\text{all}})[s] - X_z[s]\,\delta \big\rVert^2 + \lambda \lVert \delta \rVert^2 \quad \text{(closed form; } \lambda = 10).$$

$\Delta$ is regularised $\to$ small-$n$ noise shrinks to $0$; remodeling $= \lVert \delta_s \rVert_1$ vs
an $n$-matched null. **Corrects the independent version:** PTEN borderline ($\sim 1.35\times$ across
$\lambda = 1$–$30$, not robust), only RB1-Her2 modest ($\sim 1.5\times$) $\Rightarrow$ subtype wiring
mostly conserved. This is the design's $M_T = M_H + \Delta$ nesting, subtype flavour (state flavour:
same with $M_H$ from GTEx as the prior mean).

# §12 — Discovery

*Modules: `discovery.py`.*

Orphan edge $(m,g)$ beyond the curated POOLED-HE regulators:

$$\operatorname{spearman}\!\big(\operatorname{resid}(x_m \mid C \cup \text{pooledHE-aggregate}),\ \operatorname{resid}(Y \mid \cdot)\big) < \min\_\text{partial},$$

gated by scanMiR biochemical support and (`deconv`) composition retention. Edge-level FDR from a
permutation null-scan (real vs permuted hits). Family-level: collapse to seed-family $\to$ gene
hypotheses (§14).

# §13 — Enrichment

*Modules: `enrichment.py`, `stats.hypergeom_enrichment`.*

For a target set (discoveries / a hub miRNA's targets) vs Hallmark program $P$ over the HE universe $U$:

$$\text{fold} = \frac{k/\lvert P\rvert}{K/\lvert U\rvert}, \qquad p = P(X \ge k)\ \text{hypergeometric with } (\lvert U\rvert, K, \lvert P\rvert); \quad \text{BH across programs}.$$

# §14 — Multiple-testing machinery

*Modules: `stats.bh_fdr`, `coupling_inference.{benjamini_yekutieli, simes_p, family_simes_fdr}`.*

- **BH:** $q_{(i)} = \min_{j \ge i}\, p_{(j)}\, n/j$ (step-up).
- **Benjamini–Yekutieli** (arbitrary dependence): BH with $\alpha/H_n$, $H_n = \sum_{i=1}^{n} 1/i$.
- **Simes per family:** family $p = \min_i\, p_{(i)}\, k/i$ over the family's $k$ tests; BH across families $\to$ per-test $q$.
- **Family-wise min-statistic null:** per draw take the family's most-extreme null statistic; $p_{\text{family}} = \big(1 + \#\{\text{null}_{\min} \le \text{obs}_{\min}\}\big)/(n_{\text{perm}}+1)$.

Hypothesis unit = seed family (Design §F): family-mate arms are **one** test, not $N$.

# §15 — Identity vs magnitude attribution

*Modules: `attribution.identity_vs_magnitude`, Design §Decision I.*

Over the gene's nonzero seed-families, two objects — never conflated:

- **Magnitude** (budget) $= M_{\text{fam}} X_{\text{fam}}$ normalised (§5a) — realized force (abundance-included).
- **Identity** $=$ **Shapley** credit for the fixed-weight aggregate's explained variance:

$$\varphi_f = \operatorname*{mean}_{\text{orderings } \pi}\big[\, R^2(S \cup \{f\}) - R^2(S) \,\big], \quad S = f\text{'s predecessors in } \pi,$$

$$\text{pred}_S = \sum_{f' \in S} M_{f'}\, \tilde{X}_{f'}, \qquad R^2(S) = \operatorname{corr}(\text{pred}_S, \tilde{y})^2,$$

with $\tilde{X}/\tilde{y}$ $C$-residualised ($\tilde{X}$ $z$-scored). Sampled orderings ($n_{\text{perm}}=400$;
$=$ **LMG** for the linear aggregate — exact $= 2^p$). $\text{identity\_share} = \varphi_f / \sum \varphi$.

- $\text{gap} = \text{identity} - \text{magnitude}$: **$+$** = explains more of $Y$ than its budget (quiet on-target owner, e.g. ESR1 miR-18: identity $0.54$ / budget $0.20$); **$-$** = high budget but doesn't explain $Y$ (loud passenger, e.g. ESR1 miR-22: budget $0.26$ / identity $0.01$). Budget rewards abundance ($M\bar{X}$); identity rewards on-target coupling — the gap is *primarily* abundance-vs-coupling, and *secondarily* collinearity (co-varying families split the shared $R^2$).
- **Discipline:** identity $\perp$ coupling (§6/§7) — the share distributes the coupling that *exists*; it is never evidence a driver exists. This retires the precursor softmax-share identity.

# §16 — Structural identity: the loss lens

*Modules: `structural_identity.py`, Design §Decision I loss variant.*

Expression-free "who is **designed** to specifically repress $g$", so a *silenced* specialist
($M$ / budget / Shapley all $0$ — no variance, no abundance) still surfaces. Candidate families =
pooled-HE (§0) $\cap$ scanMiR-modelled arms, with a TargetScan-context++ fallback for arms scanMiR
lacks (recovers e.g. miR-137). Per seed-family $f$ of $g$:

- $\operatorname{potential}(f,g) = \text{biochem strength}$ (scanMiR $\lvert\text{repression}\rvert$ / TargetScan $\text{ts\_mag}$) $\times\ \big(1 + \log(1 + \text{evidence ledger})\big)$, direction-filtered (drop TarBase-v9 validated-activating). Family $=$ its strongest member.
- $\displaystyle \operatorname{specificity}(f,g) = \frac{\sum_{\text{members}} \operatorname{ev}(\cdot, g)}{\sum_{\text{members}} \sum_{g'} \operatorname{ev}(\cdot, g')}$ — the family's pooled evidence-mass fraction aimed at $g$ (pool the whole seed family, not one representative arm — else a thinly-studied member fakes concentration). There is no biochemical specialist (every arm hits $\sim 700$ targets $\sim$uniformly), so specificity is purely an evidence construct — no biochemical floor.
- $\operatorname{confidence}(f) = \min\!\big(1,\ \operatorname{pooled\_ev\_depth}(f)/10\big)$, $\operatorname{ev\_depth} = \#$ genes with any family evidence.
- $\operatorname{structural}(f,g) = \operatorname{potential} \cdot \operatorname{specificity} \cdot \operatorname{confidence}$ (confidence a **multiplier**: no evidence $\Rightarrow$ no ownership). $\text{structural\_share} =$ gene-local normalisation.

$\operatorname{lost\_specialist}(\text{gene})$ overlays baseline-active (GTEx/NAT $> \log_2(11)$)
$\wedge$ tumour-silenced ($<$ floor) $+$ §17.

## §16b — Sequence-designed ownership: affinity-percentile attribution-identity &nbsp;(`attribution_identity.py`, 2026-07-09)

A THIRD identity object, distinct from §15 (Shapley **realized** coupling ownership, from TCGA) and §16
(`structural_identity` = potential·specificity·confidence): the **abundance-removed, expression-free
SEQUENCE-designed owner** — *who g is designed to be owned by*. Non-circular (sequence + AGO chimeric, never TCGA
$X/Y$).

**L1 — family ownership** (seed-family unit, §4; `seq_specificity`):
$$\text{aff}(F,g)=\operatorname*{mean}_{m\in F}\ \textstyle\sum_{\text{sites}}\lvert\text{weighted context++}\rvert(m\!\to\! g),\qquad
\text{pct}(F,g)=\operatorname{rank}_{\%}\big[\text{aff}(F,\cdot)\big]_g,\qquad
\text{own}(F\mid g)=\operatorname{softmax}_F\!\big(z(\text{pct}(F,g))\big).$$
Co-seed matures are **pooled** by the canonical seed family (`seed_family_map`, from `miR_Family_Info.txt`);
$\text{pct}$ = "is $g$ among $F$'s strongest genome targets" — the **percentile**, NOT concentration-share (share is
confounded by targetome size; it crowned miR-96 over canonical miR-200 on ZEB1). context++ beats scanMiR $K_D$
(checked; $K_D$ scores weak/6mer sites). $\;$⚠ **SUPERSEDED (MH-87):** HE-restricted only — GENOME-WIDE all-site
$K_D$ beats context++ at per-arm specialist recovery (0.89 vs 0.79); at family-L1 the swap is mixed. L1 recovers canonical specialists (miR-19/26→PTEN, miR-15/16→CCND1,
miR-200/141→ZEB1), abundance-removed.

**L2 — arm nomination within $F$** (`chimeric_evidence`): a distribution over member arms from **pooled chimeric
provenance** (Manakov + TarBase qCLASH/CLASH, overlap-resolved `mean_norm`), else expression$\times$loading.
Combined arm ranking $=\text{own}(F\mid g)\cdot\text{share}(\text{arm}\mid F)$.

**Verdict:** a **descriptive** model (breast-binding-enriched $p=7.5\mathrm{e}{-7}$, seed-caveated; arm-level 91%
cross-source) that earns **NO regression value** — inert on §6 coupling (p0$\times$spec\_gain sweep), §15 magnitude,
§12 discovery, and identifiability (`identity_payoff_genome`: 0% top-driver flip genome-wide; evidence-π saturates).
Identity $\perp$ what the model estimates: a hypothesis-generator, not a prior. See VALIDATION §1, SYNTHESIS §6c.

# §17 — Loss/gain mechanism: promoter methylation gate

*Modules: `methylation.py`.*

Tumour $-$ normal promoter $\Delta\beta$ over CpG probes directly overlapping the miRNA promoter (full
sesame HM450 manifest, $485{,}577$ probes; window $[\text{TSS}-2000, \text{TSS}+500]$ around the
pri-miRNA hairpin; **no cCRE hop** — the gene-centric probe annotation misses miRNA loci, which is why
the pipeline fell back to cCRE). Group means over $802$ tumour vs $99$ normal raw betas (no barcode map
needed). **Bidirectional:**

- $\text{hyper} = (\Delta\beta \ge +0.15) \wedge (\beta_{\text{tumour}} \ge 0.20) \to$ arm **silenced** $\to$ repression **lost** (target $\uparrow$).
- $\text{hypo} = (\Delta\beta \le -0.15) \wedge (\beta_{\text{normal}} \ge 0.20) \to$ arm **de-repressed** $\to$ repression **gained** (target $\downarrow$).

Validated: TSG-miRNAs mean $\Delta\beta +0.061$ vs oncomiRs $-0.070$ (miR-124 $+0.193$, 48 CGI probes;
miR-21 $-0.176$). Promoter-probe-union betas cached at `output/learned/mirna_promoter_betas.parquet`.

# §18 — Estimator selection: what estimator for what job

*Modules: `spike_slab._gibbs_posterior`, `attribution_eb`, `attribution.shapley_identity`, `discovery.py`, `gauge.py`, `mvp.oof_gate`.*

There is **one model** — a **dense learned-$\tau^2$ non-negative Bayesian posterior per gene** — read **two**
ways:

$$
\begin{aligned}
r &= X\beta + \varepsilon, & \varepsilon &\sim \mathcal{N}(0,\sigma^2) && r = -\operatorname{resid}(Y\mid C)\\
\beta_m &= z_m\,\theta_m, & \theta_m &\sim \mathcal{N}^{+}(0,\nu^2) && \text{half-normal slab: } \beta \ge 0\ \text{("miRNAs repress")}\\
& & z_m &\sim \text{Bernoulli}(\pi_m), & \nu^2,\sigma^2 &\sim \text{InvGamma} \quad (\nu^2 \text{ LEARNED, not guessed})
\end{aligned}
$$

- $\pi \equiv 1$ (**dense**) $\to$ learned-$\tau^2$ ridge: **coupling + attribution + identifiability**.
- **evidence-$\pi$** (**inclusion**) $\to$ spike-and-slab with a genuine PIP: **discovery**.

The learned $\nu^2$ — not $L_2$ — is the edge over the lasso ($\texttt{nnridge\_cv} \approx \text{lasso}$). **The
adaptive lasso (§3) is a BASELINE, not the model.** Canonical statement: `LEARNED_MODEL_DISCOVERY_SYNTHESIS.md`
§6b; the *why* per estimator: `RATIONALE.md` §2a–§2g.

**Inclusion policy (section-wide).** Every estimator below draws a gene's candidate regulators from
$\text{POOLED-HE}$ (§0, `ledger.pooled_he_edges()`). Scope is the **Hallmark universe by design**;
`he_only=False` $\to$ the full hallmark set.

| # | Job / estimand | Estimator | Chosen over |
|---|----------------|-----------|-------------|
| 1 | **Coupling** | aggregate $X\,\mathbb{E}[\beta]$ of the **dense** ($\pi\equiv1$) posterior, on $C$-residualised $z$-scored family predictors (`_gibbs_posterior`) | **adaptive lasso** (§3): OOF mean $\rho$ $\mathbf{-0.168}$ vs $-0.152$, Wilcoxon $\mathbf{p=9\mathrm{e}{-16}}$ — at $n \gg p$ sparse selection over-shrinks the aggregate, dense learned-$\tau^2$ keeps the diffuse many-arm signal; **CV-tuned NN-ridge** ($\texttt{nnridge\_cv}\approx$ lasso $\Rightarrow$ the win is the *learned* $\nu^2$, not $L_2$); **inclusion mode** (selection discards the diffuse signal) |
| 2 | **Attribution — magnitude** | posterior **mean $\pm$ sd** of the *same* dense fit (`attribution_eb`) | **bagged NNLS** (§5): equal *and* unified — agree $0.80$, reproducibility $0.97 > 0.98\cdot\text{NNLS}$, split-half $\rho$ $0.822$ vs $0.729$; the calibrated sd comes free (bagging only approximates it) |
| 3 | **Attribution — identity / ownership** | **Shapley/LMG** decomposition of the fixed-weight aggregate's $R^2$ across families (§15; `attribution.shapley_identity`) | **budget share** (§5a — that is *force*, abundance-included: a different question); raw per-family $\rho$ (double-counts shared variance under collinearity). **Discipline:** identity $\perp$ coupling — it says *who*, not *whether* |
| 4 | **Identifiability** | $\lvert z\rvert = \lvert \text{mean}/\text{sd}\rvert > 2$ off the same dense posterior | the **§9 conditioned-partial**: equal (precision $0.86$, recall $0.89$); $\lvert z\rvert>2$ **is** the full-conditional partial-Spearman, shrinkage-stabilised (Spearman $0.88$) $\Rightarrow$ no separate estimator needed. Cross-gene hierarchical pooling (§10) does not improve coupling |
| 5 | **Discovery** | **PIP** from the **inclusion** mode (evidence-$\pi$ / learned base-rate) $+$ permutation-FDR $+$ scanMiR duplex gate $+$ deconv composition-robustness gate (`_gibbs_ss`, `priors.inclusion_prior`, §12) | the dense $\lvert z\rvert$ (over-conditions: discovery must nominate a *representative* of a collinear orphan group — model-averaging — which the full-conditional $\lvert z\rvert$ drops; the two are **complementary**, SYNTHESIS §4); NNLS / fixed support (structurally *cannot* select); the lasso (baseline) |
| 6 | **Cross-cohort GAUGE** $\beta_{\text{src}} = a\,\beta_{\text{tgt}} + c + \delta$ | **bagged NNLS** $\beta/\text{sd}$ — **not** the Gibbs posterior (`gauge.py`, `_bagged_nnls_meansd`) | the dense Gibbs posterior — see the **carve-out** below |
| 7 | **OOF prediction / the gate** | out-of-fold **partial-Spearman** of the aggregate vs held-out mRNA $\mid C$ (§6), 5-fold over participants (`mvp.oof_gate`) | any in-sample fit statistic (the estimand is held-out prediction). *Which comparator the gate is scored against is a validation question, not an estimator one* $\to$ `STATE_OF_PLAY.md` Axis 4 |
| 8 | **Per-edge coupling readout** | **partial-Spearman** of arm vs gene $\mid C$ (§7) — estimator-free | any fitted $M$ (a single edge needs no model; direct partial-$\rho$ is assumption-light and robust) |

**Engine: Gibbs, not HMC — measured, not aesthetic.** Every live channel is **Gaussian-conjugate**; the
Student-t likelihood is a **scale mixture of normals** and stays conjugate (§3b); the protein discordance link
is **additive** $\Rightarrow$ conjugate. HMC was warranted only for the A3 program-wise field, which the A4
probe **gated out** (pooling is a wash, $\Delta\rho = -0.0006$). Only a future **non-conjugate** channel
reopens this — binding $=$ Poisson/NB, methylation $=$ Beta.

**⚠ The gauge carve-out (job 6) — fit it with bagged NNLS, keep Gibbs for the model.** Fed the dense Gibbs
posterior's $\beta/\text{sd}$, the gauge returns $a = 4.1$ where the truth is $1.0$. **Gibbs is not the worse
estimator** (split-half reproducibility $\rho = 0.822$ vs bagged NNLS's $0.729$) — the **errors-in-variables
correction** is what breaks: it divides by $\operatorname{Var}(\hat\beta) - \operatorname{mean}(se^2)$, and for
Gibbs $\operatorname{mean}(se^2)$ is dominated by a heavy tail of a few enormous posterior SDs
($\sqrt{\operatorname{mean}(se^2)} = 0.055$ against a *typical* $se$ of $0.015$), collapsing the denominator
(reliability $0.17$ vs NNLS's $0.72$). `gauge._MIN_RELIABILITY` **refuses** that gauge rather than silently
returning $4.1$.

**Where the prior acts.** The prior enters as **inclusion** ($\pi$, evidence-graded, §2) and as slab-**scale**
channels ($\mu =$ scanMiR biochemical magnitude, $\tau =$ evidence depth); the slab **location stays at $0$**,
so the prior sets *ordering*, the data sets *magnitude*. The **dense ($\pi\equiv1$) fit is prior-inert by
design** — the prior's work is in inclusion/discovery, not in coupling.

---

*Every formula above has a runnable implementation at the cited module.function; see the reproduce
index in `LEARNED_MODEL_VALIDATION.md`.*
