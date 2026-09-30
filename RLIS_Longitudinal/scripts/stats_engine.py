"""
Estimation engine that reproduces the SPSS procedures used in this project.

  gee_fit()      -> SPSS GENLIN (GEE): identity link, normal distribution,
                    subject = caseid, within-subject = WAVE, exchangeable working
                    correlation, ROBUST (sandwich) covariance matrix.
  ols_hc3()      -> SPSS REGRESSION / UNIANOVA with /ROBUST=HC3, or PROCESS v4 with HC=3.
  process4()     -> PROCESS v4 Model 4 (simple mediation), HC3 SEs for a and b paths,
                    percentile bootstrap CI for the indirect effect (5,000 resamples).
  process1()     -> PROCESS v4 Model 1 (moderation) with HC3, simple slopes at -1/0/+1 SD
                    and Johnson-Neyman region.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.genmod.cov_struct import Exchangeable

COVS_BIN = ["gender", "maristat", "DOCSTAT", "EDUHI", "EMPLOY", "ORIG_SA", "ORIG_CA"]
COVS_CONT = ["AGE", "LNHINC"]


def zscore(s: pd.Series) -> pd.Series:
    """SPSS DESCRIPTIVES /SAVE: (x - mean) / SD (n-1)."""
    return (s - s.mean()) / s.std(ddof=1)


def covariate_terms(exclude=()):
    terms = [c for c in COVS_BIN if c not in exclude] + ["z" + c for c in COVS_CONT if c not in exclude]
    return terms


def prep_long(long: pd.DataFrame, cont_vars):
    d = long.copy()
    for v in cont_vars + COVS_CONT:
        d["z" + v] = zscore(d[v])  # grand-mean standardization in the person-period file
    d["W6"] = (d["WAVE"] == 1).astype(float)
    d["W7"] = (d["WAVE"] == 2).astype(float)
    return d


def prep_wide(wide: pd.DataFrame, cont_vars):
    d = wide.copy()
    for v in cont_vars + COVS_CONT:
        if v in d:
            d["z" + v] = zscore(d[v])
    return d


def gee_fit(d: pd.DataFrame, formula: str):
    m = smf.gee(formula, groups="caseid", data=d, cov_struct=Exchangeable(),
                family=sm.families.Gaussian(), missing="drop")
    r = m.fit(cov_type="robust")
    return r


def tidy(r, terms=None):
    ci = r.conf_int()
    out = pd.DataFrame({"B": r.params, "SE": r.bse, "z_or_t": r.tvalues if hasattr(r, "tvalues") else r.params / r.bse,
                        "p": r.pvalues, "LL": ci[0], "UL": ci[1]})
    if terms is not None:
        out = out.loc[[t for t in terms if t in out.index]]
    return out


def ols_hc3(d: pd.DataFrame, formula: str):
    return smf.ols(formula, data=d, missing="drop").fit(cov_type="HC3")


def process1(d, y, x, w, covs, w_is_binary=False, jn=True):
    """PROCESS Model 1 (HC3). Returns fitted model, interaction row, simple slopes, J-N bounds."""
    f = f"{y} ~ {x} * {w}" + ("" if not covs else " + " + " + ".join(covs))
    cols = [y, x, w] + covs
    dd = d[cols].dropna()
    r = smf.ols(f, data=dd).fit(cov_type="HC3")
    V = r.cov_params()
    b1, b3 = r.params[x], r.params[f"{x}:{w}"]
    levels = [0, 1] if w_is_binary else [-1.0 * dd[w].std(ddof=1) + dd[w].mean(), dd[w].mean(),
                                         dd[w].std(ddof=1) + dd[w].mean()]
    slopes = []
    for lv in levels:
        est = b1 + b3 * lv
        se = np.sqrt(V.loc[x, x] + lv ** 2 * V.loc[f"{x}:{w}", f"{x}:{w}"] + 2 * lv * V.loc[x, f"{x}:{w}"])
        tval = est / se
        from scipy import stats
        p = 2 * stats.t.sf(abs(tval), r.df_resid)
        tc = stats.t.ppf(0.975, r.df_resid)
        slopes.append({"W": lv, "effect": est, "SE": se, "t": tval, "p": p, "LL": est - tc * se, "UL": est + tc * se})
    jn_pts = None
    if jn and not w_is_binary:
        from scipy import stats
        tc = stats.t.ppf(0.975, r.df_resid)
        a = b3 ** 2 - tc ** 2 * V.loc[f"{x}:{w}", f"{x}:{w}"]
        bq = 2 * (b1 * b3 - tc ** 2 * V.loc[x, f"{x}:{w}"])
        c = b1 ** 2 - tc ** 2 * V.loc[x, x]
        disc = bq ** 2 - 4 * a * c
        if disc >= 0 and a != 0:
            roots = sorted([(-bq - np.sqrt(disc)) / (2 * a), (-bq + np.sqrt(disc)) / (2 * a)])
            wmin, wmax = dd[w].min(), dd[w].max()
            jn_pts = [(rt, float((dd[w] > rt).mean() * 100)) for rt in roots if wmin <= rt <= wmax]
    return r, pd.DataFrame(slopes), jn_pts, len(dd)


def process4(d, y, x, m, covs, nboot=5000, seed=20260924):
    """PROCESS Model 4 with HC3 SEs; covariates enter both equations (PROCESS default)."""
    cols = list(dict.fromkeys([y, x, m] + covs))
    dd = d[cols].dropna().reset_index(drop=True)
    cv = "" if not covs else " + " + " + ".join(covs)
    ra = smf.ols(f"{m} ~ {x}{cv}", dd).fit(cov_type="HC3")
    rb = smf.ols(f"{y} ~ {x} + {m}{cv}", dd).fit(cov_type="HC3")
    rc = smf.ols(f"{y} ~ {x}{cv}", dd).fit(cov_type="HC3")
    a, b = ra.params[x], rb.params[m]
    ind = a * b
    rng = np.random.default_rng(seed)
    Xa = sm.add_constant(dd[[x] + covs].values, has_constant="add")
    Xb = sm.add_constant(dd[[x, m] + covs].values, has_constant="add")
    ym, yy = dd[m].values, dd[y].values
    n = len(dd)
    boots = np.empty(nboot)
    for i in range(nboot):
        idx = rng.integers(0, n, n)
        ba = np.linalg.lstsq(Xa[idx], ym[idx], rcond=None)[0][1]
        bb = np.linalg.lstsq(Xb[idx], yy[idx], rcond=None)[0][2]
        boots[i] = ba * bb
    ll, ul = np.percentile(boots, [2.5, 97.5])
    return {"n": n, "a": a, "a_se": ra.bse[x], "a_p": ra.pvalues[x],
            "b": b, "b_se": rb.bse[m], "b_p": rb.pvalues[m],
            "c_prime": rb.params[x], "c_prime_se": rb.bse[x], "c_prime_p": rb.pvalues[x],
            "c_total": rc.params[x], "c_total_se": rc.bse[x], "c_total_p": rc.pvalues[x],
            "indirect": ind, "boot_se": boots.std(ddof=1), "boot_LL": ll, "boot_UL": ul,
            "sig": (ll > 0) or (ul < 0), "ra": ra, "rb": rb, "rc": rc}
