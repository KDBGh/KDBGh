"""
RLIS T5-T6-T7 longitudinal dataset builder.

Builds every psychometric scale from its ITEMS, using the identical item set,
identical reverse-scoring and identical scoring rule (item mean) at T5, T6 and T7,
so that scores are directly comparable across waves.  Computes Cronbach's alpha
and McDonald's omega (one-factor) per wave, merges baseline (T5) covariates from
the master file, and writes wide and long (person-period) files.

Mirrors SPSS:  RECODE (reverse items) -> COMPUTE scale = MEAN.k(items) ->
RELIABILITY /MODEL=ALPHA -> MATCH FILES /BY caseid -> VARSTOCASES.
"""
import os
import sys
import numpy as np
import pandas as pd
import pyreadstat

UP = sys.argv[1] if len(sys.argv) > 1 else "/root/.claude/uploads/640cb044-2bb0-5f80-946e-f1e1f6e27ea1"
OUT = os.path.join(os.path.dirname(__file__), "..", "output")
DATA_OUT = sys.argv[2] if len(sys.argv) > 2 else OUT  # person-level data stays out of git

FILES = {
    "T5": "e0652465-RLIS_T5_Longitudinal_Model_WORKING.sav",
    "T6": "30270963-RLIS_T6_Longitudinal_Model_WORKING.sav",
    "T7": "9fbd6525-RLIS_T7_Longitudinal_Model_WORKING1.sav",
}
MASTER = "37be42f2-RLIS_T5T6T7_Master_Longitudinal_ANALYSIS.sav"

# ------------------------------------------------------------------ scale map
# (items, reverse-scored items, (min,max) of the response scale, minimum items answered)
SCALES = {
    # ---- retained: alpha >= .70 at T5, T6 AND T7 (see output/reliability_by_wave.csv) ----
    # CES-D negative-affect items (8 of CES-D-10; the two positively worded items form a
    # method factor in Latino samples and were unreliable here: CITC .00/.28/-.36).  1-4
    "DEP": (["ds_bothered", "ds_focus", "ds_depressed", "ds_effort", "ds_fearful",
             "ds_restless", "ds_lonely", "ds_ngoing"], [], (1, 4), 6),
    # GAD-7 (Spitzer et al., 2006); 0-3
    "ANX": (["gad_nervous", "gad_nstopworry", "gad_worry", "gad_nrelax", "gad_restless",
             "gad_irritable", "gad_afraid"], [], (0, 3), 6),
    # Perceived Stress Scale - Perceived Helplessness factor (6 negatively worded items;
    # Taylor, 2015; Perera et al., 2017); 0-4
    "PSS": (["pss_unexpected", "pss_ncontrol", "pss_nervous", "pss_ncope", "pss_angered",
             "pss_difficulties"], [], (0, 4), 5),
    # MOS Social Support Survey (Sherbourne & Stewart, 1991) 1-5
    "SOCSUP": (["ssi_listen", "ssi_info", "ssi_advice", "ssi_confide", "ssi_want", "ssi_share",
                "ssi_suggestions", "ssi_understands", "ssi_bed", "ssi_doctor", "ssi_meals",
                "ssi_chores", "ssi_affection", "ssi_wanted", "ssi_hugs", "ssi_goodtime",
                "ssi_relax", "ssi_enjoyable", "ssi_mindoff"], [], (1, 5), 15),
    "EMOSUP": (["ssi_listen", "ssi_info", "ssi_advice", "ssi_confide", "ssi_want", "ssi_share",
                "ssi_suggestions", "ssi_understands"], [], (1, 5), 6),
    "TANSUP": (["ssi_bed", "ssi_doctor", "ssi_meals", "ssi_chores"], [], (1, 5), 3),
    "AFFSUP": (["ssi_affection", "ssi_wanted", "ssi_hugs"], [], (1, 5), 2),
    "POSINT": (["ssi_goodtime", "ssi_relax", "ssi_enjoyable", "ssi_mindoff"], [], (1, 5), 3),
    # Multigroup Ethnic Identity Measure (Phinney, 1992) 0-4
    "MEIM": (["mgei_learnhistory", "mgei_activeethnicgrp", "mgei_ethnicbckgrnd", "mgei_lifeaffected",
              "mgei_happy", "mgei_belonging", "mgei_understand", "mgei_talk2others", "mgei_pride",
              "mgei_cultprac", "mgei_feelgood", "mgei_attachment"], [], (0, 4), 10),
    # Marianismo Beliefs Scale (Castillo et al., 2010) 1-4, total and Silencing-Self subscale
    "MARB": (["marbs_pure", "marbs_strength", "marbs_mainstrength", "marbs_mthrunify",
              "marbs_teachloyal", "marbs_famhappy", "marbs_virgin", "marbs_wait4child",
              "marbs_religionval", "marbs_faithful2prtnr", "marbs_satisfyprtnrneeds",
              "marbs_nspkagainstmen", "marbs_respectmenop", "marbs_avoidno", "marbs_obeymale",
              "marbs_nbirthcntrl", "marbs_nneeds", "marbs_guilt", "marbs_ntalksex", "marbs_forgive",
              "marbs_agreeable", "marbs_spiritlead", "marbs_relservices", "marbs_spiritgrow"], [], (1, 4), 20),
    "MBSS": (["marbs_nbirthcntrl", "marbs_nneeds", "marbs_guilt", "marbs_ntalksex", "marbs_forgive",
              "marbs_agreeable"], [], (1, 4), 5),
    # Mexican American Cultural Values Scale (Knight et al., 2010) 1-5
    "FAMSUP": (["macv_first", "macv_security", "macv_united", "macv_close", "macv_holiday",
                "macv_affect"], [], (1, 5), 5),
    "FAMOBL": (["macv_financial", "macv_share", "macv_rmodel", "macv_sacrifice", "macv_duty"], [], (1, 5), 4),
    "RELIG": (["macv_strength", "macv_god1", "macv_faith", "macv_word", "macv_religion1",
               "macv_thank", "macv_teachpray"], [], (1, 5), 6),
    # Bicultural Involvement Questionnaire (Szapocznik et al., 1980): Hispanicism / Americanism
    "HISP": (["biq_spanhome", "biq_spanschool", "biq_spanwork", "biq_spanfriend", "biq_spangen",
              "biq_hispmusic", "biq_hispdance", "biq_hispplaces", "biq_hisprec", "biq_hisptv",
              "biq_hipradio", "biq_hispbookmag"], [], (0, 4), 10),
    "AMER": (["biq_enghome", "biq_engschool", "biq_engwork", "biq_engfriend", "biq_enggen",
              "biq_amermusic", "biq_amerdance", "biq_amerplaces", "biq_amerrec", "biq_amertv",
              "biq_amerradio", "biq_amerbookmag"], [], (0, 4), 10),
    # Neighborhood safety 1-5; raw coding 1 = "totalmente de acuerdo" (per T6 value labels), so all
    # three items are reversed so that HIGHER = SAFER neighborhood
    "NSAFE": (["ns_safewalk", "ns_nviolence", "ns_safe"], ["ns_safewalk", "ns_nviolence", "ns_safe"], (1, 5), 2),
    # Alcohol outcome expectancies 1-6
    "POSAOE": (["aoes_accepted", "aoes_outgoing", "aoes_easy2sclze", "aoes_freelytalk",
                "aoes_friendlier", "aoes_morescl", "aoes_buzz", "aoes_happy", "aoes_goodtime",
                "aoes_fun", "aoes_pleasant", "aoes_good", "aoes_mdesire", "aoes_mactive",
                "aoes_mresponsive", "aoes_massertive", "aoes_lnegative", "aoes_lstress",
                "aoes_nprob"], [], (1, 6), 15),
}

# ---- screened out: alpha < .70 at >= 1 wave (reported for transparency only) ----
EXCLUDED = {
    "CESD10": (["ds_bothered", "ds_focus", "ds_depressed", "ds_effort", "ds_hopeful", "ds_fearful",
                "ds_restless", "ds_happy", "ds_lonely", "ds_ngoing"], ["ds_hopeful", "ds_happy"], (1, 4), 8),
    "PSS10": (["pss_unexpected", "pss_ncontrol", "pss_nervous", "pss_confident", "pss_goyourway",
               "pss_ncope", "pss_irritations", "pss_ontop", "pss_angered", "pss_difficulties"],
              ["pss_confident", "pss_goyourway", "pss_irritations", "pss_ontop"], (0, 4), 8),
    "FCOH": (["ff_help", "ff_together", "ff_ntogether", "ff_getalong", "ff_avoid"],
             ["ff_ntogether", "ff_avoid"], (1, 4), 4),
    "FCONF": (["ff_fight", "ff_throw", "ff_ntemper", "ff_hit", "ff_ncriticize"],
              ["ff_ntemper", "ff_ncriticize"], (1, 4), 4),
    "MACH": (["macs_mensuperior", "macs_fatherlaw", "macs_malebirth", "macs_nweakest",
              "macs_nguarddown", "macs_cryshame", "macs_controlwife", "macs_fight",
              "macs_beautifulwomen", "macs_mansbills"], [], (1, 7), 8),
    "CABA": (["macs_manners", "macs_affectionchild", "macs_respectelder", "macs_womanloyal",
              "macs_fairness", "macs_defendfam", "macs_famimportant", "macs_mothers",
              "macs_nbragsex", "macs_childbetterlife"], [], (1, 7), 8),
    "INDEP": (["macv_selfcare", "macv_independent", "macv_decisions", "macv_onlyself",
               "macv_solve"], [], (1, 5), 4),
    # NCES cohesion: alpha adequate but validity unstable (r with CES-D = .09 / .47 / -.12 across waves)
    "NCOH": (["nces_help", "nces_closeknit", "nces_trust"], [], (1, 5), 2),
    # SRPB: alpha high but T7 distribution collapses (59% = 1) and r with religiosity ~ 0
    "SRPB": (["srpb_hardtimes", "srpb_stress", "srpb_understandothers", "srpb_comfort"], [], (1, 5), 3),
}

AUDIT_T5 = ["frqdrink", "numdrink", "frq5more", "nostopdk", "failnorm", "drinkam", "guiltdrk",
            "unabremb", "injurdrk", "concerdk"]
AUDIT_T67 = ["audit_often", "audit_howmany", "audit_5ormore", "audit_cantstop", "audit_didntdo",
             "audit_morningdrink", "audit_guilt", "audit_cantremember", "audit_injured",
             "audit_dctrconcern"]


def cronbach_alpha(items: pd.DataFrame) -> float:
    x = items.dropna()
    k = x.shape[1]
    if k < 2 or len(x) < 10:
        return np.nan
    return k / (k - 1) * (1 - x.var(ddof=1).sum() / x.sum(axis=1).var(ddof=1))


def omega_total(items: pd.DataFrame) -> float:
    """McDonald's omega from a one-factor principal-axis solution (iterated)."""
    x = items.dropna()
    x = x.loc[:, x.std() > 0]
    if x.shape[1] < 3 or len(x) < 10:
        return np.nan
    r = np.corrcoef(x.values, rowvar=False)
    h2 = 1 - 1 / np.diag(np.linalg.pinv(r))  # SMC start
    for _ in range(200):
        rr = r.copy()
        np.fill_diagonal(rr, h2)
        w, v = np.linalg.eigh(rr)
        lam = v[:, -1] * np.sqrt(max(w[-1], 1e-9))
        new = np.clip(lam ** 2, 1e-4, 0.999)
        if np.max(np.abs(new - h2)) < 1e-6:
            break
        h2 = new
    lam = np.abs(lam)
    return lam.sum() ** 2 / (lam.sum() ** 2 + (1 - lam ** 2).sum())


def load(t):
    df, _ = pyreadstat.read_sav(os.path.join(UP, FILES[t]), encoding="latin1")
    return df


def clean_items(df, items, lo, hi):
    x = df[items].apply(pd.to_numeric, errors="coerce")
    return x.where((x >= lo) & (x <= hi))  # -88/-99/out-of-range (e.g., 93) -> system missing


def main():
    rel_rows = []
    wide = None
    for t in ["T5", "T6", "T7"]:
        df = load(t)
        w = t[-1]
        out = pd.DataFrame({"caseid": df["caseid"]})
        for name, (items, rev, (lo, hi), kmin) in {**SCALES, **EXCLUDED}.items():
            x = clean_items(df, items, lo, hi)
            for r in rev:
                x[r] = lo + hi - x[r]
            score = x.mean(axis=1).where(x.notna().sum(axis=1) >= kmin)
            if name in SCALES:
                out[f"{name}{w}"] = score
            rel_rows.append({"scale": name, "retained": name in SCALES or name == "AUDIT", "wave": t, "k": len(items), "n": int(score.notna().sum()),
                             "alpha": cronbach_alpha(x), "omega": omega_total(x),
                             "mean": score.mean(), "sd": score.std()})
        # AUDIT-10 total (0-40): items 1-8 scored 0-4, items 9-10 scored 0/2/4
        a_items = AUDIT_T5 if t == "T5" else AUDIT_T67
        a = df[a_items].apply(pd.to_numeric, errors="coerce")
        a = a.where((a >= 0) & (a <= 4))
        # non-drinkers (frequency = never) skip items 2-10 -> scored 0
        never = a[a_items[0]] == 0
        a.loc[never, a_items[1:]] = a.loc[never, a_items[1:]].fillna(0)
        out[f"AUDIT{w}"] = a.sum(axis=1).where(a.notna().sum(axis=1) >= 8)
        rel_rows.append({"scale": "AUDIT", "retained": True, "wave": t, "k": 10, "n": int(out[f"AUDIT{w}"].notna().sum()),
                         "alpha": cronbach_alpha(a), "omega": omega_total(a),
                         "mean": out[f"AUDIT{w}"].mean(), "sd": out[f"AUDIT{w}"].std()})
        wide = out if wide is None else wide.merge(out, on="caseid", how="outer")

    master, _ = pyreadstat.read_sav(os.path.join(UP, MASTER))
    t5 = load("T5")
    cov = master[["caseid", "gender", "maristat", "DOCSTAT", "EDUHI", "EMPLOY", "ORIG_SA", "ORIG_CA"]].copy()
    cov = cov.merge(t5[["caseid", "demo_age", "demo_hhinc"]], on="caseid", how="left")
    cov["AGE"] = cov["demo_age"].where(cov["demo_age"] > 0)
    inc = cov["demo_hhinc"].where(cov["demo_hhinc"] >= 0)
    cov["LNHINC"] = np.log(inc + 1)  # identical to master LNHINC5 (log household income)
    cov = cov.drop(columns=["demo_age", "demo_hhinc"])
    wide = cov.merge(wide, on="caseid", how="right")
    wide = wide[wide["caseid"].notna()]

    rel = pd.DataFrame(rel_rows)
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(DATA_OUT, exist_ok=True)
    rel.round(3).to_csv(os.path.join(OUT, "reliability_by_wave.csv"), index=False)
    wide.to_pickle(os.path.join(DATA_OUT, "wide.pkl"))
    wide.to_csv(os.path.join(DATA_OUT, "RLIS_T5T6T7_wide_scales.csv"), index=False)

    # person-period (long) file: SPSS VARSTOCASES equivalent
    stems = list(SCALES) + ["AUDIT"]
    rows = []
    for w, wave in zip("567", [0, 1, 2]):
        part = wide[["caseid", "gender", "maristat", "DOCSTAT", "EDUHI", "EMPLOY", "ORIG_SA", "ORIG_CA",
                     "AGE", "LNHINC"] + [f"{s}{w}" for s in stems]].copy()
        part.columns = [c[:-1] if c[:-1] in stems and c[-1] == w else c for c in part.columns]
        part["WAVE"] = wave
        rows.append(part)
    long = pd.concat(rows, ignore_index=True).sort_values(["caseid", "WAVE"])
    long.to_pickle(os.path.join(DATA_OUT, "long.pkl"))
    print(rel.pivot(index="scale", columns="wave", values="alpha").round(3))
    print("wide", wide.shape, "long", long.shape)


if __name__ == "__main__":
    main()
