"""Within-wave mediation X_t -> M_t -> AUDIT_t, estimated separately at T5, T6, T7 (PROCESS Model 4,
HC3, bootstrap) and pooled (GEE a- and b-paths, Monte Carlo CI).  Logged to output/audit/."""
import os, sys, itertools, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from stats_engine import *
warnings.filterwarnings("ignore")
DATA = sys.argv[1]
OUT = os.path.join(os.path.dirname(__file__), "..", "output", "audit")
long = pd.read_pickle(os.path.join(DATA, "long.pkl"))
PRED = ["DEP", "ANX", "PSS", "SOCSUP", "EMOSUP", "TANSUP", "AFFSUP", "POSINT", "MEIM", "MARB", "MBSS",
        "FAMSUP", "FAMOBL", "RELIG", "HISP", "AMER", "NSAFE"]
SUP = {"SOCSUP", "EMOSUP", "TANSUP", "AFFSUP", "POSINT"}
L = prep_long(long, PRED + ["AUDIT"])
COV = covariate_terms()
waves = {}
for wv in [0, 1, 2]:
    d = long[long.WAVE == wv].copy()
    for v in PRED + ["AUDIT"] + COVS_CONT:
        d["z" + v] = zscore(d[v])
    waves[wv] = d.reset_index(drop=True)
rng = np.random.default_rng(1)
rows = []
for x, m in itertools.product(PRED, PRED):
    if x == m or {x, m} <= SUP or {x, m} <= {"MARB", "MBSS"}:
        continue
    for adj in (False, True):
        cov = COV if adj else []
        rec = {"X": x, "M": m, "adjusted": adj}
        allsig = True
        for wv in [0, 1, 2]:
            res = process4(waves[wv], "zAUDIT", "z" + x, "z" + m, cov, nboot=500, seed=wv + 7)
            rec[f"ind_T{wv+5}"] = res["indirect"]; rec[f"LL_T{wv+5}"] = res["boot_LL"]; rec[f"UL_T{wv+5}"] = res["boot_UL"]
            rec[f"sig_T{wv+5}"] = res["sig"]
        cv = "" if not cov else " + " + " + ".join(cov)
        ra = gee_fit(L, f"z{m} ~ z{x} + W6 + W7{cv}")
        rb = gee_fit(L, f"zAUDIT ~ z{x} + z{m} + W6 + W7{cv}")
        a, sa, b, sb = ra.params["z" + x], ra.bse["z" + x], rb.params["z" + m], rb.bse["z" + m]
        mc = rng.normal(a, sa, 20000) * rng.normal(b, sb, 20000)
        rec.update({"a_pooled": a, "a_p": ra.pvalues["z" + x], "b_pooled": b, "b_p": rb.pvalues["z" + m],
                    "ab_pooled": a * b, "MC_LL": np.percentile(mc, 2.5), "MC_UL": np.percentile(mc, 97.5)})
        rec["pooled_sig"] = rec["MC_LL"] > 0 or rec["MC_UL"] < 0
        rec["n_waves_sig"] = sum(rec[f"sig_T{k}"] for k in (5, 6, 7))
        rows.append(rec)
R = pd.DataFrame(rows)
R.to_csv(os.path.join(OUT, "audit_screen_E_withinwave_mediation.csv"), index=False)
print("E", len(R), "pooled sig adj", R[R.adjusted].pooled_sig.sum(), "3-wave sig adj", (R[R.adjusted].n_waves_sig == 3).sum())
