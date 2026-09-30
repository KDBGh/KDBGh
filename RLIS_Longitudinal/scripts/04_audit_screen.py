"""
AUDIT-as-outcome screen (all designs).  Every test is logged to output/audit_screen_*.csv.

A  three-wave GEE moderation      zAUDIT_t ~ X_t * W_t + wave (+ covs)
B  wave interaction               zAUDIT_t ~ X_t * wave (+ covs)
C  three-wave mediation           X5 -> M6 -> AUDIT7 | M5, AUDIT5 (+ covs), HC3, bootstrap
D  prospective lagged moderation  AUDIT7 ~ X5 * W5 + AUDIT5 (+ covs), HC3
E  three-way                      zAUDIT_t ~ X_t * W_t * wave (+ covs)  [for A hits]
"""
import os, sys, itertools, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from stats_engine import *
warnings.filterwarnings("ignore")
DATA = sys.argv[1]
OUT = os.path.join(os.path.dirname(__file__), "..", "output", "audit")
os.makedirs(OUT, exist_ok=True)
long = pd.read_pickle(os.path.join(DATA, "long.pkl"))
wide = pd.read_pickle(os.path.join(DATA, "wide.pkl"))
CONT = ["AUDIT", "DEP", "ANX", "PSS", "SOCSUP", "EMOSUP", "TANSUP", "AFFSUP", "POSINT", "MEIM", "MARB",
        "MBSS", "FAMSUP", "FAMOBL", "RELIG", "HISP", "AMER", "NSAFE"]
L = prep_long(long, CONT)
COV = covariate_terms()
PRED = [c for c in CONT if c != "AUDIT"]
MODS = PRED + ["gender", "EDUHI", "EMPLOY", "maristat", "ORIG_SA", "ORIG_CA"]
SUP = {"SOCSUP", "EMOSUP", "TANSUP", "AFFSUP", "POSINT"}


def skip(x, w):
    return x == w or ({x, w} <= SUP) or ({x, w} <= {"MARB", "MBSS"})


def wald(r, names):
    R = np.zeros((len(names), len(r.params)))
    for i, n in enumerate(names):
        R[i, list(r.params.index).index(n)] = 1
    t = r.wald_test(R, scalar=True)
    return float(t.statistic), float(t.pvalue)


# A ---------------------------------------------------------------------------
rows = []
for x in PRED:
    for w in MODS:
        if skip(x, w):
            continue
        zx = "z" + x
        zw = w if w in COVS_BIN else "z" + w
        if (w, x) in [(r["X"], r["W"]) for r in rows[-400:]]:
            pass
        for adj in (False, True):
            cov = [c for c in (COV if adj else []) if c != zw]
            f = f"zAUDIT ~ {zx} * {zw} + W6 + W7" + ("" if not cov else " + " + " + ".join(cov))
            r = gee_fit(L, f)
            t = f"{zx}:{zw}"
            rows.append({"X": x, "W": w, "adjusted": adj, "B_int": r.params[t], "SE": r.bse[t], "p_int": r.pvalues[t],
                         "B_X": r.params[zx], "p_X": r.pvalues[zx], "B_W": r.params[zw], "p_W": r.pvalues[zw],
                         "N_obs": int(r.nobs)})
A = pd.DataFrame(rows)
A.to_csv(os.path.join(OUT, "audit_screen_A_gee_moderation.csv"), index=False)
print("A", len(A), "sig adj", (A[A.adjusted].p_int < .05).sum(), "/", A.adjusted.sum(), flush=True)

# B ---------------------------------------------------------------------------
rows = []
for x in PRED + ["gender", "EDUHI", "EMPLOY", "maristat", "ORIG_SA", "ORIG_CA"]:
    zx = x if x in COVS_BIN else "z" + x
    for adj in (False, True):
        cov = [c for c in (COV if adj else []) if c != zx]
        f = f"zAUDIT ~ {zx} * (W6 + W7)" + ("" if not cov else " + " + " + ".join(cov))
        r = gee_fit(L, f)
        c, p = wald(r, [f"{zx}:W6", f"{zx}:W7"])
        rows.append({"X": x, "adjusted": adj, "B_X_T5": r.params[zx], "p_X_T5": r.pvalues[zx],
                     "B_XxT6": r.params[f"{zx}:W6"], "p_XxT6": r.pvalues[f"{zx}:W6"],
                     "B_XxT7": r.params[f"{zx}:W7"], "p_XxT7": r.pvalues[f"{zx}:W7"], "Wald_chi2": c, "Wald_p": p})
B = pd.DataFrame(rows)
B.to_csv(os.path.join(OUT, "audit_screen_B_wave_interaction.csv"), index=False)
print("B", len(B), flush=True)

# C ---------------------------------------------------------------------------
W = wide.copy()
for v in CONT:
    for k in "567":
        if f"{v}{k}" in W:
            W[f"z{v}{k}"] = zscore(W[f"{v}{k}"])
for v in COVS_CONT:
    W["z" + v] = zscore(W[v])
rows = []
for x, m in itertools.product(PRED, PRED):
    if skip(x, m):
        continue
    for adj in (False, True):
        covs = [f"z{m}5", "zAUDIT5"] + (COV if adj else [])
        res = process4(W, "zAUDIT7", f"z{x}5", f"z{m}6", covs, nboot=1000)
        rows.append({"X": x, "M": m, "adjusted": adj, **{k: v for k, v in res.items() if k not in ("ra", "rb", "rc")}})
C = pd.DataFrame(rows)
C.to_csv(os.path.join(OUT, "audit_screen_C_mediation.csv"), index=False)
print("C", len(C), "sig adj", C[C.adjusted].sig.sum(), flush=True)

# D ---------------------------------------------------------------------------
rows = []
for x in PRED:
    for w in MODS:
        if skip(x, w):
            continue
        zw = w if w in COVS_BIN else f"z{w}5"
        for adj in (False, True):
            cov = ["zAUDIT5"] + [c for c in (COV if adj else []) if c != zw]
            r, sl, jn, n = process1(W, "zAUDIT7", f"z{x}5", zw, cov, w_is_binary=w in COVS_BIN, jn=False)
            t = f"z{x}5:{zw}"
            rows.append({"X": x, "W": w, "adjusted": adj, "B_int": r.params[t], "SE_HC3": r.bse[t],
                         "p_int": r.pvalues[t], "B_X": r.params[f"z{x}5"], "p_X": r.pvalues[f"z{x}5"], "n": n})
    for adj in (False, True):  # lagged main effect
        cov = ["zAUDIT5"] + (COV if adj else [])
        r = ols_hc3(W, f"zAUDIT7 ~ z{x}5 + " + " + ".join(cov))
        rows.append({"X": x, "W": "(main effect)", "adjusted": adj, "B_X": r.params[f"z{x}5"], "p_X": r.pvalues[f"z{x}5"],
                     "n": int(r.nobs)})
D = pd.DataFrame(rows)
D.to_csv(os.path.join(OUT, "audit_screen_D_lagged.csv"), index=False)
print("D", len(D), flush=True)
