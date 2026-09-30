"""
Final five longitudinal models (RLIS T5-T6-T7), unadjusted and adjusted.

All continuous variables are z-standardized.  Pooled three-wave models use GEE
(SPSS GENLIN: normal/identity, exchangeable, robust sandwich SEs).  Every model is
also re-estimated separately at T5, T6 and T7 with OLS + HC3 SEs (SPSS PROCESS hc=3),
and the mediation model uses PROCESS Model 4 logic with HC3 SEs and 5,000 percentile
bootstrap resamples.
"""
import os, sys, warnings
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from stats_engine import *  # noqa
warnings.filterwarnings("ignore")

DATA = sys.argv[1]
OUT = os.path.join(os.path.dirname(__file__), "..", "output", "final")
os.makedirs(OUT, exist_ok=True)
long = pd.read_pickle(os.path.join(DATA, "long.pkl"))
wide = pd.read_pickle(os.path.join(DATA, "wide.pkl"))
CONT = ["DEP", "ANX", "PSS", "SOCSUP", "EMOSUP", "RELIG", "HISP", "AMER", "NSAFE"]
L = prep_long(long, CONT)
COV = covariate_terms()
LABEL = {"zDEP": "Depressive symptoms", "zANX": "Anxiety symptoms", "zPSS": "Perceived stress",
         "zRELIG": "Religiosity", "zAMER": "Americanism", "zHISP": "Hispanicism",
         "zEMOSUP": "Emotional support", "zSOCSUP": "Social support", "zNSAFE": "Neighborhood safety",
         "W6": "Wave T6 (vs T5)", "W7": "Wave T7 (vs T5)", "gender": "Male (vs female)",
         "maristat": "Married/partnered", "DOCSTAT": "Documented", "EDUHI": "Education > high school",
         "EMPLOY": "Employed", "ORIG_SA": "South American (vs Cuban)", "ORIG_CA": "Central American (vs Cuban)",
         "zAGE": "Age", "zLNHINC": "Household income (log)"}
report = []


def lincom(r, weights):
    """Linear combination of GEE coefficients: estimate, SE, z, p, 95% CI (SPSS EMMEANS/TEST equivalent)."""
    names = list(r.params.index)
    w = np.zeros(len(names))
    for k, v in weights.items():
        w[names.index(k)] = v
    est = float(w @ r.params.values)
    se = float(np.sqrt(w @ r.cov_params().values @ w))
    z = est / se
    p = 2 * stats.norm.sf(abs(z))
    return {"estimate": est, "SE": se, "z": z, "p": p, "LL": est - 1.96 * se, "UL": est + 1.96 * se}


def jn_gee(r, x, xw, wname, wvals):
    b1, b3 = r.params[x], r.params[xw]
    V = r.cov_params()
    a = b3 ** 2 - 1.96 ** 2 * V.loc[xw, xw]
    b = 2 * (b1 * b3 - 1.96 ** 2 * V.loc[x, xw])
    c = b1 ** 2 - 1.96 ** 2 * V.loc[x, x]
    disc = b ** 2 - 4 * a * c
    if disc < 0 or a == 0:
        return []
    roots = sorted([(-b - np.sqrt(disc)) / (2 * a), (-b + np.sqrt(disc)) / (2 * a)])
    return [(round(rt, 3), round(float((wvals > rt).mean() * 100), 1)) for rt in roots if wvals.min() <= rt <= wvals.max()]


def gee_table(r, model, adj):
    t = tidy(r)
    t.insert(0, "term", [LABEL.get(i, i) for i in t.index])
    t["Wald_chi2"] = (t["B"] / t["SE"]) ** 2
    t.insert(0, "adjusted", adj)
    t.insert(0, "model", model)
    return t


def wave_ols(y, x, w, adj, w_binary=False):
    rows = []
    for wv in [0, 1, 2]:
        d = long[long.WAVE == wv].copy()
        for v in {y, x, w} - set(COVS_BIN) | set(COVS_CONT):
            d["z" + v] = zscore(d[v])
        zw = w if w in COVS_BIN else "z" + w
        cov = [c for c in (COV if adj else []) if c != zw]
        r, sl, jn, n = process1(d, "z" + y, "z" + x, zw, cov, w_is_binary=w_binary)
        t = f"z{x}:{zw}"
        rows.append({"wave": f"T{wv + 5}", "adjusted": adj, "n": n, "B_X": r.params["z" + x], "p_X": r.pvalues["z" + x],
                     "B_W": r.params[zw], "p_W": r.pvalues[zw], "B_int": r.params[t], "SE_HC3": r.bse[t],
                     "t": r.tvalues[t], "p_int": r.pvalues[t], "R2": r.rsquared,
                     "slope_lowW": sl.effect[0], "p_lowW": sl.p[0], "slope_meanW": sl.effect[1], "p_meanW": sl.p[1],
                     "slope_highW": sl.effect[2], "p_highW": sl.p[2], "JN": jn})
    return pd.DataFrame(rows)


def moderation_model(tag, y, x, w):
    zy, zx, zw = "z" + y, "z" + x, "z" + w
    tabs, slopes = [], []
    for adj in (False, True):
        f = f"{zy} ~ {zx} * {zw} + W6 + W7" + (" + " + " + ".join(COV) if adj else "")
        r = gee_fit(L, f)
        tabs.append(gee_table(r, tag, adj))
        for lvl, lab in [(-1, "-1 SD"), (0, "Mean"), (1, "+1 SD")]:
            s = lincom(r, {zx: 1, f"{zx}:{zw}": lvl})
            slopes.append({"model": tag, "adjusted": adj, "moderator_level": f"{w} {lab}", **s})
        jn = jn_gee(r, zx, f"{zx}:{zw}", zw, L[zw].dropna())
        slopes.append({"model": tag, "adjusted": adj, "moderator_level": f"J-N boundary (z, % above): {jn}"})
        report.append(f"{tag} {'ADJ' if adj else 'UNADJ'}: N_obs={int(r.nobs)} N_persons={len(np.unique(r.model.groups))} "
                      f"{x}x{w} B={r.params[f'{zx}:{zw}']:.3f} SE={r.bse[f'{zx}:{zw}']:.3f} p={r.pvalues[f'{zx}:{zw}']:.2e}")
    pd.concat(tabs).round(4).to_csv(os.path.join(OUT, f"{tag}_GEE_coefficients.csv"))
    pd.DataFrame(slopes).round(4).to_csv(os.path.join(OUT, f"{tag}_simple_slopes.csv"), index=False)
    pw = pd.concat([wave_ols(y, x, w, False), wave_ols(y, x, w, True)])
    pw.round(4).to_csv(os.path.join(OUT, f"{tag}_per_wave_HC3.csv"), index=False)


def wave_model(tag, y, x):
    zy, zx = "z" + y, "z" + x
    tabs, slopes = [], []
    for adj in (False, True):
        f = f"{zy} ~ {zx} * (W6 + W7)" + (" + " + " + ".join(COV) if adj else "")
        r = gee_fit(L, f)
        tabs.append(gee_table(r, tag, adj))
        R = np.zeros((2, len(r.params)))
        for i, nme in enumerate([f"{zx}:W6", f"{zx}:W7"]):
            R[i, list(r.params.index).index(nme)] = 1
        wt = r.wald_test(R, scalar=True)
        for wv, wts in [("T5", {zx: 1}), ("T6", {zx: 1, f"{zx}:W6": 1}), ("T7", {zx: 1, f"{zx}:W7": 1})]:
            slopes.append({"model": tag, "adjusted": adj, "wave": wv, **lincom(r, wts)})
        slopes.append({"model": tag, "adjusted": adj, "wave": f"Omnibus {x} x Wave Wald chi2(2)={float(wt.statistic):.2f}, p={float(wt.pvalue):.2e}"})
        report.append(f"{tag} {'ADJ' if adj else 'UNADJ'}: {x} x Wave Wald chi2(2)={float(wt.statistic):.2f} p={float(wt.pvalue):.2e}")
    pd.concat(tabs).round(4).to_csv(os.path.join(OUT, f"{tag}_GEE_coefficients.csv"))
    pd.DataFrame(slopes).round(4).to_csv(os.path.join(OUT, f"{tag}_wave_specific_effects.csv"), index=False)


def three_way_model(tag, y, x, w):
    zy, zx, zw = "z" + y, "z" + x, "z" + w
    tabs, eff = [], []
    for adj in (False, True):
        f = f"{zy} ~ {zx} * {zw} * (W6 + W7)" + (" + " + " + ".join(COV) if adj else "")
        r = gee_fit(L, f)
        tabs.append(gee_table(r, tag, adj))
        R = np.zeros((2, len(r.params)))
        for i, nme in enumerate([f"{zx}:{zw}:W6", f"{zx}:{zw}:W7"]):
            R[i, list(r.params.index).index(nme)] = 1
        wt = r.wald_test(R, scalar=True)
        for wv, extra in [("T5", ""), ("T6", "W6"), ("T7", "W7")]:
            # X x W interaction at each wave
            wts = {f"{zx}:{zw}": 1}
            if extra:
                wts[f"{zx}:{zw}:{extra}"] = 1
            eff.append({"model": tag, "adjusted": adj, "wave": wv, "quantity": f"{x} x {w} interaction", **lincom(r, wts)})
            for lvl, lab in [(-1, "-1 SD"), (1, "+1 SD")]:
                ws = {zx: 1, f"{zx}:{zw}": lvl}
                if extra:
                    ws[f"{zx}:{extra}"] = 1
                    ws[f"{zx}:{zw}:{extra}"] = lvl
                eff.append({"model": tag, "adjusted": adj, "wave": wv, "quantity": f"slope of {x} at {w} {lab}", **lincom(r, ws)})
        eff.append({"model": tag, "adjusted": adj, "wave": f"Omnibus {x} x {w} x Wave Wald chi2(2)={float(wt.statistic):.2f}, p={float(wt.pvalue):.2e}"})
        report.append(f"{tag} {'ADJ' if adj else 'UNADJ'}: {x} x {w} x Wave Wald chi2(2)={float(wt.statistic):.2f} p={float(wt.pvalue):.2e}")
    pd.concat(tabs).round(4).to_csv(os.path.join(OUT, f"{tag}_GEE_coefficients.csv"))
    pd.DataFrame(eff).round(4).to_csv(os.path.join(OUT, f"{tag}_wave_specific_effects.csv"), index=False)


def mediation_model(tag, x, m, ys):
    W = wide.copy()
    for v in CONT:
        for w_ in "567":
            if f"{v}{w_}" in W:
                W[f"z{v}{w_}"] = zscore(W[f"{v}{w_}"])
    for v in COVS_CONT:
        W["z" + v] = zscore(W[v])
    rows = []
    for y in ys:
        for adj in (False, True):
            covs = [f"z{m}5", f"z{y}5"] + (COV if adj else [])
            res = process4(W, f"z{y}7", f"z{x}5", f"z{m}6", covs, nboot=5000)
            for path, rr in [("a: X->M", res["ra"]), ("b & c': X,M->Y", res["rb"]), ("c: total", res["rc"])]:
                t = tidy(rr)
                t.insert(0, "term", [LABEL.get(i.rstrip("567"), i) + (f" (T{i[-1]})" if i[-1] in "567" else "") for i in t.index])
                t.insert(0, "equation", path)
                t.insert(0, "Y", y)
                t.insert(0, "adjusted", adj)
                t.to_csv(os.path.join(OUT, f"{tag}_{y}_{'adj' if adj else 'unadj'}_{path[0]}_HC3.csv"))
            rows.append({"model": tag, "Y": y, "adjusted": adj, **{k: v for k, v in res.items() if k not in ("ra", "rb", "rc")}})
            report.append(f"{tag} Y={y} {'ADJ' if adj else 'UNADJ'}: n={res['n']} a={res['a']:.3f}(p={res['a_p']:.3f}) "
                          f"b={res['b']:.3f}(p={res['b_p']:.3f}) ab={res['indirect']:.4f} 95%CI[{res['boot_LL']:.4f},{res['boot_UL']:.4f}]")
    pd.DataFrame(rows).round(4).to_csv(os.path.join(OUT, f"{tag}_indirect_effects.csv"), index=False)


def descriptives():
    rows = []
    for v in CONT + ["AUDIT"]:
        for wv in [0, 1, 2]:
            s = long.loc[long.WAVE == wv, v]
            rows.append({"variable": v, "wave": f"T{wv + 5}", "n": int(s.notna().sum()), "M": s.mean(), "SD": s.std(),
                         "min": s.min(), "max": s.max(), "skew": s.skew(), "kurtosis": s.kurt()})
    pd.DataFrame(rows).round(3).to_csv(os.path.join(OUT, "Table1b_descriptives_by_wave.csv"), index=False)
    base = wide.copy()
    t1 = []
    for c in COVS_BIN:
        t1.append({"variable": LABEL.get(c, c), "n": int(base[c].notna().sum()), "value": f"{base[c].mean() * 100:.1f}%"})
    for c in COVS_CONT:
        t1.append({"variable": c, "n": int(base[c].notna().sum()), "value": f"M={base[c].mean():.2f}, SD={base[c].std():.2f}"})
    pd.DataFrame(t1).to_csv(os.path.join(OUT, "Table1a_baseline_sample.csv"), index=False)
    # attrition: retained at T7 vs lost, on T5 study variables
    base["retainedT7"] = base["DEP7"].notna()
    att = []
    for v in [f"{c}5" for c in CONT] + COVS_CONT + COVS_BIN:
        a, b = base.loc[base.retainedT7, v].dropna(), base.loc[~base.retainedT7, v].dropna()
        tt = stats.ttest_ind(a, b, equal_var=False)
        att.append({"variable": v, "M_retained": a.mean(), "M_lost": b.mean(), "n_retained": len(a), "n_lost": len(b),
                    "Welch_t": tt.statistic, "p": tt.pvalue})
    pd.DataFrame(att).round(3).to_csv(os.path.join(OUT, "Table_attrition.csv"), index=False)
    for wv in [0, 1, 2]:
        long.loc[long.WAVE == wv, CONT].corr().round(2).to_csv(os.path.join(OUT, f"Table2_correlations_T{wv + 5}.csv"))


if __name__ == "__main__":
    descriptives()
    moderation_model("M1a_Religiosity_x_Stress_on_Depression", "DEP", "PSS", "RELIG")
    moderation_model("M1b_Religiosity_x_Stress_on_Anxiety", "ANX", "PSS", "RELIG")
    moderation_model("M2_EmotionalSupport_x_Americanism_on_Anxiety", "ANX", "AMER", "EMOSUP")
    wave_model("M3a_Hispanicism_x_Wave_on_Depression", "DEP", "HISP")
    wave_model("M3b_Hispanicism_x_Wave_on_Anxiety", "ANX", "HISP")
    wave_model("M3c_Hispanicism_x_Wave_on_PerceivedStress", "PSS", "HISP")
    three_way_model("M3d_Stress_x_Hispanicism_x_Wave_on_Depression", "DEP", "PSS", "HISP")
    three_way_model("M4_Stress_x_EmotionalSupport_x_Wave_on_Anxiety", "ANX", "PSS", "EMOSUP")
    moderation_model("M4s_Stress_x_EmotionalSupport_on_Anxiety_pooled", "ANX", "PSS", "EMOSUP")
    mediation_model("M5_NeighborhoodSafety_SocialSupport_Stress", "NSAFE", "SOCSUP", ["PSS", "DEP", "ANX"])
    open(os.path.join(OUT, "RESULTS_SUMMARY.txt"), "w").write("\n".join(report) + "\n")
    print("\n".join(report))
