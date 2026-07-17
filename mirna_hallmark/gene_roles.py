"""Per-gene cancer-role (oncogene / tumor-suppressor) annotation.

A curated `malignancy_sign` per gene so that *repressing* a node can be scored for its
malignancy direction independently of the gene-gene OmniPath topology (which is blind to
non-transcriptional effectors such as PTEN's lipid-phosphatase activity). Convention:

    malignancy_sign =  +1   canonical oncogene  (repressing it is ANTI-tumor)
                       -1   tumor suppressor    (repressing it is PRO-tumor)
                        0   dual-role / context-dependent / unknown

The repression-direction credit used downstream is therefore `-malignancy_sign`:
silencing a TSG (sign -1) contributes +1 toward malignancy; silencing an oncogene (+1)
contributes -1.

The built-in dictionary covers high-confidence drivers (COSMIC Cancer Gene Census /
OncoKB consensus, with breast-cancer relevance prioritised). It is intentionally
conservative: genes not listed get sign 0 (unknown) rather than a guessed polarity.
Drop a two-column TSV (`gene`, `role` in {oncogene,tsg,oncogene/tsg}) at
`config.GENE_ROLES_OVERRIDE` to extend/override without editing code.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

# --- canonical oncogenes (repression = anti-tumor) ------------------------- #
_ONCOGENES = {
    "ABL1", "AKT1", "AKT2", "AKT3", "ALK", "AR", "AURKA", "AURKB", "BCL2", "BCL2L1",
    "BCL6", "BRAF", "CCND1", "CCND2", "CCND3", "CCNE1", "CDK4", "CDK6", "CTNNB1",
    "EGFR", "ERBB2", "ERBB3", "ERBB4", "ETV1", "ETV4", "ETV5", "EZH2", "FGFR1",
    "FGFR2", "FGFR3", "FOSL1", "FOXM1", "GAB2", "HRAS", "IGF1R", "JAK2", "JUN",
    "KIT", "KRAS", "MAPK1", "MAP3K8", "MDM2", "MDM4", "MET", "MITF", "MYB", "MYC",
    "MYCN", "NRAS", "NTRK1", "PDGFRA", "PDGFRB", "PIK3CA", "PIK3CB", "PIK3CD",
    "PIM1", "PLK1", "PPM1D", "RAF1", "REL", "RELA", "RET", "ROS1", "SKP2", "SRC",
    "STAT3", "TERT", "WNT1", "YAP1", "YES1", "MTOR", "RHEB", "RAC1", "CDC25A",
    "E2F1", "E2F3", "FOXA1", "GATA3", "SOX2", "RICTOR", "HMGA1", "HMGA2", "MELK",
    "BIRC5", "TWIST1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "PTK2", "PTK2B", "PRKCA",
    "S1PR1", "SPHK1", "TGFB1", "VEGFA", "HIF1A", "EPAS1", "SLC2A1", "LDHA", "FASN",
    "ACLY", "GLS", "MTDH", "MUC1", "CKS1B", "CKS2", "MELK", "TOP2A", "DLGAP5",
    "KIF11", "MYBL2", "CDC20", "BUB1", "TTK", "CENPF", "RRM2", "TYMS",
}

# --- canonical tumor suppressors (repression = pro-tumor) ------------------ #
_TSGS = {
    "APC", "ARID1A", "ARID1B", "ARID2", "ATM", "ATR", "ATRX", "BAP1", "BRCA1",
    "BRCA2", "CDH1", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B", "CHEK1", "CHEK2",
    "CREBBP", "CYLD", "DICER1", "EP300", "FAT1", "FBXW7", "FOXO1", "FOXO3",
    "GATA4", "KEAP1", "KMT2C", "KMT2D", "LATS1", "LATS2", "LKB1", "STK11", "MAP2K4",
    "MEN1", "MLH1", "MSH2", "MSH6", "NF1", "NF2", "NOTCH2", "PALB2", "PBRM1",
    "PDCD4", "PHLPP1", "PHLPP2", "PIK3R1", "PTEN", "PTPN14", "RASA1", "RB1",
    "RBL2", "SETD2", "SFRP1", "SMAD2", "SMAD3", "SMAD4", "SMARCA4", "SMARCB1",
    "SOCS1", "SOCS3", "SOCS5", "SPRED1", "SPRY2", "STK11", "TGFBR1", "TGFBR2",
    "TP53", "TP53BP1", "TP63", "TSC1", "TSC2", "VHL", "WT1", "WWOX", "DLC1",
    "DAB2", "DAB2IP", "MEOX2", "RUNX3", "TIMP3", "MAP3K7", "INPP4B", "DEPTOR",
    "TRIB3", "AKT1S1", "THEM4", "DUSP1", "DUSP4", "DUSP6", "GADD45A", "GADD45B",
    "BTG2", "KLF4", "KLF6", "NDRG1", "NDRG2", "IRF1", "STAT1", "BNIP3", "BNIP3L",
    "NPRL2", "TXNIP", "CEBPA", "MXD1", "MNT", "TCF21", "ID4",
}

# --- well-known dual / context-dependent (assign 0 explicitly) ------------- #
_DUAL = {"NOTCH1", "TGFB1", "MYB", "WT1", "RUNX1", "FOXP1", "BCL6", "ETS1",
         "KLF5", "ATF3", "NFKB1", "CEBPB"}


def _builtin_table() -> pd.DataFrame:
    rows = []
    for g in sorted(_ONCOGENES - _DUAL):
        rows.append((g, "oncogene", 1))
    for g in sorted(_TSGS - _DUAL - _ONCOGENES):
        rows.append((g, "tsg", -1))
    for g in sorted(_DUAL):
        rows.append((g, "oncogene/tsg", 0))
    return pd.DataFrame(rows, columns=["gene", "role", "malignancy_sign"])


def load_gene_roles(genes: Optional[Iterable[str]] = None,
                    override_path: Optional[Path] = None) -> pd.DataFrame:
    """Return gene -> {role, malignancy_sign}. Optional external TSV (`gene`,`role`)
    is merged on top of the built-in dictionary (override wins). Genes not annotated
    are returned with sign 0 only when an explicit gene list is supplied."""
    tab = _builtin_table()
    path = override_path
    if path is None:
        try:
            from mirna_hallmark import config as C
            path = getattr(C, "GENE_ROLES_OVERRIDE", None)
        except Exception:
            path = None
    if path is not None and Path(path).exists():
        ext = pd.read_csv(path, sep="\t")
        ext["gene"] = ext["gene"].astype(str)
        role_to_sign = {"oncogene": 1, "tsg": -1, "tumor_suppressor": -1,
                        "oncogene/tsg": 0, "dual": 0}
        ext["malignancy_sign"] = ext["role"].str.lower().map(role_to_sign).fillna(0).astype(int)
        tab = pd.concat([tab[~tab.gene.isin(ext.gene)], ext[["gene", "role", "malignancy_sign"]]],
                        ignore_index=True)
    tab = tab.drop_duplicates("gene")
    if genes is not None:
        idx = pd.Index(sorted(set(genes)), name="gene")
        tab = tab.set_index("gene").reindex(idx).reset_index()
        tab["role"] = tab["role"].fillna("unknown")
        tab["malignancy_sign"] = tab["malignancy_sign"].fillna(0).astype(int)
    return tab


def load_gene_dependency(genes: Optional[Iterable[str]] = None, *,
                         panel: str = "breast") -> pd.DataFrame:
    """Continuous, ALL-GENE role weight from DepMap CRISPR (Chronos) gene-effect — the graded
    counterpart to the binary `malignancy_sign`. Returns gene -> `dep_role_weight` oriented so it reads
    like `-malignancy_sign` (the repression-direction credit):

        dep_role_weight(g) = +gene_effect(g)   (breast-panel mean by default)

    gene-effect << 0 ⇒ essential / tumour *dependency* (oncogene-like) ⇒ weight < 0 (repressing it is
    ANTI-tumor); gene-effect ≳ 0 ⇒ TSG-like ⇒ weight > 0 (repressing it is pro-tumor). Defined for ~18.5k
    genes, not just the 232 curated drivers. **One-sided caveat**: CRISPR essentiality resolves the
    dependency (negative) tail well but the TSG (positive) tail weakly, so this score is the all-gene
    graded version of the *onco-collateral / anti-tumor* axis, complementing (not replacing) the curated
    TSG-credit. Built by `_build_depmap_dependency.py` → `data/CCLE/depmap_gene_effect_summary.tsv`."""
    from mirna_hallmark import config as C
    src = C.REPO_ROOT / "data/CCLE/depmap_gene_effect_summary.tsv"
    if not src.exists():
        raise FileNotFoundError(f"{src} missing — run `python -m mirna_hallmark.analyses.builders._build_depmap_dependency`")
    dep = pd.read_csv(src, sep="\t")
    col = "dep_effect_breast" if panel == "breast" else "dep_effect_pan"
    dep["dep_role_weight"] = pd.to_numeric(dep[col], errors="coerce")
    out = dep[["gene", "dep_role_weight", "dep_effect_breast", "dep_effect_pan"]].copy()
    if genes is not None:
        idx = pd.Index(sorted(set(genes)), name="gene")
        out = out.set_index("gene").reindex(idx).reset_index()
        # genes with no CRISPR readout stay neutral (0) rather than NaN-propagating
        out["dep_role_weight"] = out["dep_role_weight"].fillna(0.0)
    return out


def load_tf_census(genes: Optional[Iterable[str]] = None) -> pd.DataFrame:
    """Curated human transcription-factor census (Lambert et al. 2018, 1,639 TFs) — an identity-based
    hierarchy annotation that complements the OmniPath-induced `w_arch` (which only weakly tracks TF
    identity: ~48% of curated TFs fall below median w_arch because transcriptional regulation is
    under-covered in the signed-PPI network). Returns gene -> `is_tf` (bool) + `dbd` (DNA-binding domain
    family). Built by `_build_tf_census.py` → `annotations/humantfs_lambert2018_tf.tsv`."""
    from mirna_hallmark import config as C
    src = C.REPO_ROOT / "annotations/humantfs_lambert2018_tf.tsv"
    if not src.exists():
        raise FileNotFoundError(f"{src} missing — run `python -m mirna_hallmark.analyses.builders._build_tf_census`")
    tf = pd.read_csv(src, sep="\t")[["gene", "is_tf", "dbd"]]
    if genes is not None:
        idx = pd.Index(sorted(set(genes)), name="gene")
        tf = tf.set_index("gene").reindex(idx).reset_index()
        tf["is_tf"] = tf["is_tf"].fillna(False).astype(bool)
        tf["dbd"] = tf["dbd"].fillna("")
    return tf


if __name__ == "__main__":
    t = load_gene_roles()
    print(f"built-in roles: {len(t)} genes "
          f"({(t.malignancy_sign == 1).sum()} oncogene, "
          f"{(t.malignancy_sign == -1).sum()} TSG, "
          f"{(t.malignancy_sign == 0).sum()} dual)")
