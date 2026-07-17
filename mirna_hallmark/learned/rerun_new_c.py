"""Detached re-run of the key learned-model outputs under the CORRECTED confounder block
(design-core C = CPE + target-CN + malignant-proliferation, HRD off; CIBERSORTx tumour+NAT / metagene GTEx
cross-state). Everything picks up the new default C automatically via data.assemble_gene. 2026-07-05."""
import time, traceback
def step(name, fn):
    t = time.time(); print(f"\n{'='*70}\n[rerun] START {name}  {time.strftime('%H:%M:%S')}\n{'='*70}", flush=True)
    try:
        fn(); print(f"[rerun] DONE {name}  ({(time.time()-t)/60:.1f} min)", flush=True)
    except Exception as e:
        print(f"[rerun] FAILED {name}: {e}\n{traceback.format_exc()[:800]}", flush=True)
if __name__ == "__main__":
    from mirna_hallmark.learned import programs, discovery, states
    step("PROGRAMS cell-intrinsic core C (5 Hallmarks)", lambda: programs.run())
    step("PROGRAMS +deconv (composition pass)", lambda: programs.run(deconv=True))
    step("DISCOVERY genome-wide (FDR + deconv top-150)", lambda: discovery.run_all(top=150))
    step("BUDGET_SHIFT hub (GTEx→NAT→tumour)",
         lambda: [states.budget_shift(g) for g in ["PTEN","ESR1","ZEB1","BCL2","GATA3","CDKN1A"]])
    step("REALIZATION hub (paired within-patient)",
         lambda: states.realization(["PTEN","ZEB1","ESR1","GATA3","BCL2","CDKN1A","PDCD4","RB1"]))
    print("\n[rerun] ===== ALL DONE =====", flush=True)
