"""Figure 1 for AUDIT Model 1: (A) predicted AUDIT across anxiety at low/high familism support,
(B) Johnson-Neyman plot of the conditional anxiety effect across familism support. Adjusted GEE."""
import os, sys, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from stats_engine import *
warnings.filterwarnings("ignore")
DATA = sys.argv[1]
OUT = os.path.join(os.path.dirname(__file__), "..", "output", "audit", "final", "figures")
L = prep_long(pd.read_pickle(os.path.join(DATA, "long.pkl")), ["AUDIT", "ANX", "FAMSUP"])
COV = covariate_terms()
r = gee_fit(L, "zAUDIT ~ zANX*zFAMSUP + W6 + W7 + " + " + ".join(COV))
b, V = r.params, r.cov_params()
dd = L.dropna(subset=["zAUDIT", "zANX", "zFAMSUP"] + COV)
m, s = L.AUDIT.mean(), L.AUDIT.std()
anx_sd, anx_m = L.ANX.std(), L.ANX.mean()

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
C = {-1: "#2a78d6", 1: "#eb6834"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": INK2,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.2), facecolor=SURF)
for ax in (a1, a2):
    ax.set_facecolor(SURF); ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)

# Panel A: predictions (covariates and wave at sample means)
zgrid = np.linspace(dd.zANX.min(), 2.5, 60)
base = np.array([1.0 if k == "Intercept" else dd[k].mean() if k in dd else 0.0 for k in b.index])
for f, lab in [(-1, "Low familism support (−1 SD)"), (1, "High familism support (+1 SD)")]:
    pred, lo, hi = [], [], []
    for z in zgrid:
        x = base.copy()
        x[list(b.index).index("zANX")] = z
        x[list(b.index).index("zFAMSUP")] = f
        x[list(b.index).index("zANX:zFAMSUP")] = z * f
        e = x @ b.values; se = np.sqrt(x @ V.values @ x)
        pred.append(m + e * s); lo.append(m + (e - 1.96 * se) * s); hi.append(m + (e + 1.96 * se) * s)
    raw = anx_m + zgrid * anx_sd
    a1.fill_between(raw, lo, hi, color=C[f], alpha=0.12, lw=0)
    a1.plot(raw, pred, color=C[f], lw=2)
    a1.annotate(lab.split(" (")[0], (raw[-1], pred[-1]), xytext=(4, 0), textcoords="offset points",
                color=INK, fontsize=9, va="center")
a1.set_xlabel("Anxiety symptoms (GAD-7 item mean, 0–3)")
a1.set_ylabel("Predicted AUDIT score")
a1.set_title("A. Predicted hazardous drinking by anxiety", loc="left", color=INK, fontsize=11)
a1.set_xlim(0, (anx_m + 2.5 * anx_sd) * 1.45)

# Panel B: Johnson-Neyman conditional effect
fg = np.linspace(dd.zFAMSUP.min(), dd.zFAMSUP.max(), 200)
eff = b["zANX"] + b["zANX:zFAMSUP"] * fg
se = np.sqrt(V.loc["zANX", "zANX"] + fg ** 2 * V.loc["zANX:zFAMSUP", "zANX:zFAMSUP"] + 2 * fg * V.loc["zANX", "zANX:zFAMSUP"])
lo, hi = eff - 1.96 * se, eff + 1.96 * se
sig = (lo > 0) | (hi < 0)
a2.fill_between(fg, lo, hi, color="#2a78d6", alpha=0.12, lw=0)
a2.plot(fg, eff, color="#2a78d6", lw=2)
a2.axhline(0, color=INK2, lw=1)
roots = fg[np.where(np.diff(sig.astype(int)) != 0)[0]]
for rt in roots:
    a2.axvline(rt, color=INK2, lw=1, ls=(0, (3, 3)))
    a2.annotate(f"z = {rt:.2f}", (rt, hi.max()), xytext=(3, -2), textcoords="offset points", fontsize=8.5, color=INK2, va="top")
a2.set_xlabel("Familism support values (z)")
a2.set_ylabel("Effect of anxiety on AUDIT (β, 95% CI)")
a2.set_title("B. Johnson–Neyman regions of significance", loc="left", color=INK, fontsize=11)
fig.text(0.01, -0.03, "Adjusted three-wave GEE (robust SEs); covariates held at sample means. "
         "Dashed lines mark where the anxiety effect becomes significant (p < .05).", fontsize=8.5, color=INK2)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"Figure1_Model1_Anxiety_x_Familism_AUDIT.{ext}"), dpi=300, bbox_inches="tight", facecolor=SURF)
print("roots", roots)
