"""
Systematic, theory-constrained screen of candidate longitudinal models.
Every test is logged (output/screen_*.csv) so the full search is transparent.
"""
import os, sys, itertools, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from stats_engine import *
warnings.filterwarnings("ignore")
DATA = sys.argv[1]
OUT = os.path.join(os.path.dirname(__file__), "..", "output")
long = pd.read_pickle(os.path.join(DATA, "long.pkl"))
wide = pd.read_pickle(os.path.join(DATA, "wide.pkl"))
CONT = ["DEP", "ANX", "PSS", "SOCSUP", "EMOSUP", "MEIM", "MARB", "MBSS", "FAMSUP", "FAMOBL",
        "RELIG", "HISP", "AMER", "NCOH", "NSAFE", "POSAOE", "AUDIT"]
L = prep_long(long, CONT)
COV = covariate_terms()
COV_NOG = covariate_terms(exclude=("gender",))

# ---------------- A. Pooled 3-wave GEE moderation  Y_t ~ X_t * W_t + wave (+ covs)
rows = []
outcomes = ["DEP", "ANX", "AUDIT", "PSS"]
preds = ["PSS", "SOCSUP", "EMOSUP", "MEIM", "FAMSUP", "FAMOBL", "RELIG", "HISP", "AMER",
         "NCOH", "NSAFE", "MARB", "MBSS", "DEP", "ANX"]
mods = ["SOCSUP", "EMOSUP", "MEIM", "FAMSUP", "FAMOBL", "RELIG", "HISP", "AMER", "NCOH",
        "NSAFE", "MARB", "MBSS", "gender", "DOCSTAT", "EDUHI", "ORIG_SA", "ORIG_CA"]
for y in outcomes:
    for x in preds:
        for w in mods:
            if len({x, y, w}) < 3: continue
            if {x, w} <= {"SOCSUP", "EMOSUP"} or {x, w} <= {"MARB", "MBSS"}: continue
            zx, zy = "z" + x, "z" + y
            zw = w if w in COVS_BIN else "z" + w
            for adj in (False, True):
                cov = [c for c in (COV if adj else []) if c != zw]
                f = f"{zy} ~ {zx} * {zw} + W6 + W7" + ("" if not cov else " + " + " + ".join(cov))
                try:
                    r = gee_fit(L, f)
                except Exception as e:
                    continue
                t = f"{zx}:{zw}"
                rows.append({"design": "GEE_3wave", "Y": y, "X": x, "W": w, "adjusted": adj,
                             "B_int": r.params[t], "SE": r.bse[t], "p_int": r.pvalues[t],
                             "B_X": r.params[zx], "p_X": r.pvalues[zx], "N_obs": int(r.nobs)})
modA = pd.DataFrame(rows)
modA.to_csv(os.path.join(OUT, "screen_A_gee_moderation.csv"), index=False)
print("A done", len(modA))

# ---------------- B. Wave interaction  Y_t ~ X_t * WAVE (+ covs)
rows = []
for y in outcomes:
    for x in preds:
        if x == y: continue
        zx, zy = "z" + x, "z" + y
        for adj in (False, True):
            cov = COV if adj else []
            f = f"{zy} ~ {zx} * (W6 + W7)" + ("" if not cov else " + " + " + ".join(cov))
            r = gee_fit(L, f)
            # joint Wald test of the two interaction terms
            names = [f"{zx}:W6", f"{zx}:W7"]
            R = np.zeros((2, len(r.params)))
            for i, nme in enumerate(names): R[i, list(r.params.index).index(nme)] = 1
            wt = r.wald_test(R, scalar=True)
            rows.append({"design": "GEE_wave", "Y": y, "X": x, "adjusted": adj,
                         "B_X_T5": r.params[zx], "p_X": r.pvalues[zx],
                         "B_XxT6": r.params[names[0]], "p_XxT6": r.pvalues[names[0]],
                         "B_XxT7": r.params[names[1]], "p_XxT7": r.pvalues[names[1]],
                         "Wald_chi2": float(wt.statistic), "Wald_p": float(wt.pvalue)})
modB = pd.DataFrame(rows)
modB.to_csv(os.path.join(OUT, "screen_B_wave_interaction.csv"), index=False)
print("B done", len(modB))

# ---------------- C. Three-wave half-longitudinal mediation X5 -> M6 -> Y7 | M5, Y5
W = wide.copy()
for v in CONT:
    for w_ in "567":
        if f"{v}{w_}" in W: W[f"z{v}{w_}"] = zscore(W[f"{v}{w_}"])
for v in ["AGE", "LNHINC"]:
    W["z" + v] = zscore(W[v])
rows = []
Xs = ["MEIM", "FAMSUP", "FAMOBL", "RELIG", "HISP", "AMER", "NCOH", "NSAFE", "MARB", "MBSS",
      "SOCSUP", "EMOSUP", "PSS"]
Ms = ["PSS", "SOCSUP", "EMOSUP", "DEP", "ANX", "MEIM", "NCOH", "FAMSUP"]
Ys = ["DEP", "ANX", "AUDIT", "PSS"]
for x, m, y in itertools.product(Xs, Ms, Ys):
    if len({x, m, y}) < 3: continue
    if {x, m} <= {"SOCSUP", "EMOSUP"}: continue
    base = [f"z{m}5", f"z{y}5"]
    for adj in (False, True):
        covs = base + (COV if adj else [])
        res = process4(W, f"z{y}7", f"z{x}5", f"z{m}6", covs, nboot=1000)
        rows.append({"design": "MED_T5T6T7", "X": x, "M": m, "Y": y, "adjusted": adj,
                     **{k: v for k, v in res.items() if k not in ("ra", "rb", "rc")}})
modC = pd.DataFrame(rows)
modC.to_csv(os.path.join(OUT, "screen_C_mediation.csv"), index=False)
print("C done", len(modC))

# ---------------- D. Prospective lagged moderation Y7 ~ X5 * W5 + Y5 (+covs), HC3
rows = []
for y in outcomes:
    for x in preds:
        for w in mods:
            if len({x, y, w}) < 3: continue
            if {x, w} <= {"SOCSUP", "EMOSUP"} or {x, w} <= {"MARB", "MBSS"}: continue
            zw = w if w in COVS_BIN else f"z{w}5"
            for adj in (False, True):
                cov = [f"z{y}5"] + [c for c in (COV if adj else []) if c != zw]
                try:
                    r, sl, jn, n = process1(W, f"z{y}7", f"z{x}5", zw, cov, w_is_binary=w in COVS_BIN, jn=False)
                except Exception:
                    continue
                t = f"z{x}5:{zw}"
                rows.append({"design": "LAG_T5toT7", "Y": y, "X": x, "W": w, "adjusted": adj,
                             "B_int": r.params[t], "SE_HC3": r.bse[t], "p_int": r.pvalues[t], "n": n})
modD = pd.DataFrame(rows)
modD.to_csv(os.path.join(OUT, "screen_D_lagged_moderation.csv"), index=False)
print("D done", len(modD))
