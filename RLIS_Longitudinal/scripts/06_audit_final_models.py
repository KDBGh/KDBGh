"""
Final five AUDIT-outcome longitudinal models (RLIS T5-T6-T7), unadjusted and adjusted.

Primary estimator: GEE (SPSS GENLIN) normal/identity on z-AUDIT, exchangeable, robust SEs.
Per-wave replication: OLS + HC3 (SPSS PROCESS hc=3).
Sensitivity: (i) negative-binomial GEE on raw AUDIT (log link, ancillary = 1, the SPSS GENLIN
NEGBIN(1) default), (ii) linear GEE on z[ln(AUDIT + 1)].
"""
import os, sys, warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.genmod.cov_struct import Exchangeable
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from stats_engine import *  # noqa
warnings.filterwarnings("ignore")

DATA = sys.argv[1]
OUT = os.path.join(os.path.dirname(__file__), "..", "output", "audit", "final")
os.makedirs(OUT, exist_ok=True)
long = pd.read_pickle(os.path.join(DATA, "long.pkl"))
CONT = ["AUDIT", "ANX", "DEP", "PSS", "FAMSUP", "EMOSUP", "MEIM", "HISP"]
L = prep_long(long, CONT)
L["zlnAUDIT"] = zscore(np.log1p(L["AUDIT"]))
COV = covariate_terms()
LABEL = {"zAUDIT": "Hazardous drinking (AUDIT)", "zANX": "Anxiety symptoms", "zDEP": "Depressive symptoms",
         "zPSS": "Perceived stress", "zFAMSUP": "Familism support values", "zEMOSUP": "Emotional support",
         "zMEIM": "Ethnic identity (MEIM)", "zHISP": "Hispanicism", "maristat": "Married/partnered",
         "W6": "Wave T6 (vs T5)", "W7": "Wave T7 (vs T5)", "gender": "Male (vs female)", "DOCSTAT": "Documented",
         "EDUHI": "Education > high school", "EMPLOY": "Employed", "ORIG_SA": "South American (vs Cuban)",
         "ORIG_CA": "Central American (vs Cuban)", "zAGE": "Age", "zLNHINC": "Household income (log)"}
report = []


def cov_for(adj, excl=()):
    return [c for c in COV if c not in excl] if adj else []


def rhs(terms, adj, excl=()):
    c = cov_for(adj, excl)
    return terms + ("" if not c else " + " + " + ".join(c))


def nb_gee(formula):
    return smf.gee(formula, groups="caseid", data=L, cov_struct=Exchangeable(),
                   family=sm.families.NegativeBinomial(alpha=1.0), missing="drop").fit(cov_type="robust")


def lincom(r, weights):
    names = list(r.params.index)
    w = np.zeros(len(names))
    for k, v in weights.items():
        w[names.index(k)] = v
    est = float(w @ r.params.values)
    se = float(np.sqrt(w @ r.cov_params().values @ w))
    z = est / se
    return {"estimate": est, "SE": se, "z": z, "p": 2 * stats.norm.sf(abs(z)), "LL": est - 1.96 * se, "UL": est + 1.96 * se}


def wald(r, names):
    R = np.zeros((len(names), len(r.params)))
    for i, n in enumerate(names):
        R[i, list(r.params.index).index(n)] = 1
    t = r.wald_test(R, scalar=True)
    return float(t.statistic), float(t.pvalue)


def jn(r, x, xw, wv):
    b1, b3, V = r.params[x], r.params[xw], r.cov_params()
    a = b3 ** 2 - 1.96 ** 2 * V.loc[xw, xw]
    b = 2 * (b1 * b3 - 1.96 ** 2 * V.loc[x, xw])
    c = b1 ** 2 - 1.96 ** 2 * V.loc[x, x]
    disc = b ** 2 - 4 * a * c
    if disc < 0 or a == 0:
        return []
    roots = sorted([(-b - np.sqrt(disc)) / (2 * a), (-b + np.sqrt(disc)) / (2 * a)])
    return [(round(rt, 3), round(float((wv > rt).mean() * 100), 1)) for rt in roots if wv.min() <= rt <= wv.max()]


def table(r, tag, adj, est="GEE-normal (z-AUDIT)"):
    t = tidy(r)
    t.insert(0, "term", [LABEL.get(i, i) for i in t.index])
    if "NB" in est:
        t["IRR"], t["IRR_LL"], t["IRR_UL"] = np.exp(t.B), np.exp(t.LL), np.exp(t.UL)
    t.insert(0, "estimator", est)
    t.insert(0, "adjusted", adj)
    t.insert(0, "model", tag)
    return t


def per_wave(x, w, adj):
    rows = []
    wbin = w in COVS_BIN
    for k in [0, 1, 2]:
        d = long[long.WAVE == k].copy()
        for v in {"AUDIT", x, w} - set(COVS_BIN) | set(COVS_CONT):
            d["z" + v] = zscore(d[v])
        zw = w if wbin else "z" + w
        r, sl, jnp, n = process1(d, "zAUDIT", "z" + x, zw, cov_for(adj, (zw,)), w_is_binary=wbin)
        t = f"z{x}:{zw}"
        rows.append({"wave": f"T{k + 5}", "adjusted": adj, "n": n, "B_int": r.params[t], "SE_HC3": r.bse[t],
                     "t": r.tvalues[t], "p_int": r.pvalues[t], "R2": r.rsquared,
                     "slope_W_low": sl.effect[0], "p_low": sl.p[0], "slope_W_high": sl.effect[len(sl) - 1],
                     "p_high": sl.p[len(sl) - 1]})
    return rows


def moderation(tag, x, w):
    wbin = w in COVS_BIN
    zx, zw = "z" + x, (w if wbin else "z" + w)
    xw = f"{zx}:{zw}"
    tabs, slopes, sens, pw = [], [], [], []
    for adj in (False, True):
        r = gee_fit(L, "zAUDIT ~ " + rhs(f"{zx} * {zw} + W6 + W7", adj, (zw,)))
        tabs.append(table(r, tag, adj))
        levels = [(0, f"{w} = 0"), (1, f"{w} = 1")] if wbin else [(-1, "-1 SD"), (0, "Mean"), (1, "+1 SD")]
        for lv, lab in levels:
            slopes.append({"model": tag, "adjusted": adj, "level": lab, **lincom(r, {zx: 1, xw: lv})})
        if not wbin:
            slopes.append({"model": tag, "adjusted": adj, "level": f"J-N (z, % above) {jn(r, zx, xw, L[zw].dropna())}"})
        rn = nb_gee("AUDIT ~ " + rhs(f"{zx} * {zw} + W6 + W7", adj, (zw,)))
        tabs.append(table(rn, tag, adj, "GEE-NB(1) raw AUDIT"))
        rl = gee_fit(L, "zlnAUDIT ~ " + rhs(f"{zx} * {zw} + W6 + W7", adj, (zw,)))
        sens.append({"model": tag, "adjusted": adj, "linear_B": r.params[xw], "linear_p": r.pvalues[xw],
                     "NB_IRR": np.exp(rn.params[xw]), "NB_p": rn.pvalues[xw], "lnAUDIT_B": rl.params[xw],
                     "lnAUDIT_p": rl.pvalues[xw]})
        pw += per_wave(x, w, adj)
        report.append(f"{tag} {'ADJ' if adj else 'UNADJ'}: N_obs={int(r.nobs)} {x}x{w} B={r.params[xw]:.3f} "
                      f"SE={r.bse[xw]:.3f} p={r.pvalues[xw]:.2e} | NB IRR={np.exp(rn.params[xw]):.3f} p={rn.pvalues[xw]:.2e}"
                      f" | ln p={rl.pvalues[xw]:.2e}")
    pd.concat(tabs).round(4).to_csv(os.path.join(OUT, f"{tag}_coefficients.csv"))
    pd.DataFrame(slopes).round(4).to_csv(os.path.join(OUT, f"{tag}_simple_slopes.csv"), index=False)
    pd.DataFrame(sens).round(4).to_csv(os.path.join(OUT, f"{tag}_sensitivity.csv"), index=False)
    pd.DataFrame(pw).round(4).to_csv(os.path.join(OUT, f"{tag}_per_wave_HC3.csv"), index=False)


def wave_model(tag, x):
    zx = "z" + x
    names = [f"{zx}:W6", f"{zx}:W7"]
    tabs, eff, sens = [], [], []
    for adj in (False, True):
        r = gee_fit(L, "zAUDIT ~ " + rhs(f"{zx} * (W6 + W7)", adj))
        tabs.append(table(r, tag, adj))
        c, p = wald(r, names)
        for wv, wts in [("T5", {zx: 1}), ("T6", {zx: 1, names[0]: 1}), ("T7", {zx: 1, names[1]: 1})]:
            eff.append({"model": tag, "adjusted": adj, "wave": wv, **lincom(r, wts)})
        eff.append({"model": tag, "adjusted": adj, "wave": f"Omnibus {x} x Wave Wald chi2(2)={c:.2f}, p={p:.2e}"})
        rn = nb_gee("AUDIT ~ " + rhs(f"{zx} * (W6 + W7)", adj))
        tabs.append(table(rn, tag, adj, "GEE-NB(1) raw AUDIT"))
        cn, pn = wald(rn, names)
        rl = gee_fit(L, "zlnAUDIT ~ " + rhs(f"{zx} * (W6 + W7)", adj))
        cl, pl = wald(rl, names)
        for wv, wts in [("T5", {zx: 1}), ("T6", {zx: 1, names[0]: 1}), ("T7", {zx: 1, names[1]: 1})]:
            e = lincom(rn, wts)
            eff.append({"model": tag, "adjusted": adj, "wave": f"{wv} (NB, IRR={np.exp(e['estimate']):.3f})", **e})
        sens.append({"model": tag, "adjusted": adj, "linear_chi2": c, "linear_p": p, "NB_chi2": cn, "NB_p": pn,
                     "ln_chi2": cl, "ln_p": pl})
        report.append(f"{tag} {'ADJ' if adj else 'UNADJ'}: {x} x Wave Wald chi2(2)={c:.2f} p={p:.2e} | NB p={pn:.2e} | ln p={pl:.2e}")
    pd.concat(tabs).round(4).to_csv(os.path.join(OUT, f"{tag}_coefficients.csv"))
    pd.DataFrame(eff).round(4).to_csv(os.path.join(OUT, f"{tag}_wave_specific_effects.csv"), index=False)
    pd.DataFrame(sens).round(4).to_csv(os.path.join(OUT, f"{tag}_sensitivity.csv"), index=False)
    # per-wave HC3 simple regressions of AUDIT on X at each wave
    rows = []
    for adj in (False, True):
        for k in [0, 1, 2]:
            d = long[long.WAVE == k].copy()
            for v in ["AUDIT", x] + COVS_CONT:
                d["z" + v] = zscore(d[v])
            rr = ols_hc3(d, "zAUDIT ~ " + rhs(zx, adj))
            rows.append({"wave": f"T{k + 5}", "adjusted": adj, "n": int(rr.nobs), "B": rr.params[zx], "SE_HC3": rr.bse[zx],
                         "p": rr.pvalues[zx]})
    pd.DataFrame(rows).round(4).to_csv(os.path.join(OUT, f"{tag}_per_wave_HC3.csv"), index=False)


def descriptives():
    rows = []
    for v in CONT:
        for k in [0, 1, 2]:
            s = long.loc[long.WAVE == k, v]
            rows.append({"variable": v, "wave": f"T{k + 5}", "n": int(s.notna().sum()), "M": s.mean(), "SD": s.std(),
                         "min": s.min(), "max": s.max(), "skew": s.skew()})
    for k in [0, 1, 2]:
        s = long.loc[long.WAVE == k, "AUDIT"]
        rows.append({"variable": "AUDIT >= 8 (%)", "wave": f"T{k + 5}", "n": int(s.notna().sum()), "M": (s >= 8).mean() * 100})
        rows.append({"variable": "AUDIT = 0 (%)", "wave": f"T{k + 5}", "n": int(s.notna().sum()), "M": (s == 0).mean() * 100})
    pd.DataFrame(rows).round(3).to_csv(os.path.join(OUT, "Table1_AUDIT_study_variables_by_wave.csv"), index=False)
    for k in [0, 1, 2]:
        long.loc[long.WAVE == k, CONT].corr().round(2).to_csv(os.path.join(OUT, f"Table2_correlations_T{k + 5}.csv"))


if __name__ == "__main__":
    descriptives()
    moderation("A1_Anxiety_x_FamilismSupport_on_AUDIT", "ANX", "FAMSUP")
    moderation("A1s_Depression_x_FamilismSupport_on_AUDIT", "DEP", "FAMSUP")
    moderation("A2_Stress_x_Partnered_on_AUDIT", "PSS", "maristat")
    wave_model("A3_Anxiety_x_Wave_on_AUDIT", "ANX")
    wave_model("A3s_Stress_x_Wave_on_AUDIT", "PSS")
    wave_model("A4_EmotionalSupport_x_Wave_on_AUDIT", "EMOSUP")
    moderation("A5_EthnicIdentity_x_Hispanicism_on_AUDIT", "MEIM", "HISP")
    open(os.path.join(OUT, "RESULTS_SUMMARY_AUDIT.txt"), "w").write("\n".join(report) + "\n")
    print("\n".join(report))
