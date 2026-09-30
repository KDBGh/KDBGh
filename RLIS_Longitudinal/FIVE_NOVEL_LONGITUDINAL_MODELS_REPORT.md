# RLIS T5–T6–T7: Five Novel Longitudinal Models

**Data.** Recent Latino Immigrant Study, three annual waves (T5 n = 518, T6 n = 505, T7 n = 493; 1,554 person-waves; 518 persons). Baseline covariates come from the master file.

**Covariates (adjusted models).** Gender (male = 1), age (z), household income (ln, z), education (> high school = 1), marital status (married/partnered = 1), employment (employed = 1), documentation status (documented = 1), South American origin, Central American origin (Cuban = reference).

**Standardization.** All continuous variables are z-scored. Pooled models use grand-mean z-scores in the long file. Per-wave models use within-wave z-scores. Binary covariates are left 0/1.

**Standard errors.**
- Pooled three-wave models use GEE (SPSS GENLIN: normal/identity link, exchangeable working correlation, **robust sandwich** covariance). This is the cluster-robust extension; plain HC3 assumes independent rows and is invalid for repeated measures.
- Every model is re-estimated at T5, T6 and T7 separately with **OLS + HC3** (SPSS PROCESS `/hc=3`).
- The mediation model uses **PROCESS Model 4 with HC3** and 5,000 percentile bootstrap resamples.

**Excluded by request.** Fatalism (your Paper 3).

---

## 1. Measures retained (α ≥ .70 at *every* wave)

| Construct | Instrument (items) | α T5 | α T6 | α T7 |
|---|---|---|---|---|
| Depressive symptoms | CES-D negative-affect items (8) | .80 | .86 | .89 |
| Anxiety symptoms | GAD-7 (7) | .89 | .88 | .93 |
| Perceived stress | PSS perceived-helplessness factor (6) | .79 | .81 | .76 |
| Social support (total) | MOS-SSS (19) | .98 | .96 | .97 |
| Emotional/informational support | MOS-SSS subscale (8) | .96 | .97 | .96 |
| Religiosity | MACVS Religion (7) | .94 | .91 | .94 |
| Hispanicism | BIQ (12) | .87 | .78 | .85 |
| Americanism | BIQ (12) | .95 | .90 | .91 |
| Neighborhood safety (higher = safer) | 3 items | .77 | .90 | .74 |

**Screened out: α < .70 at ≥ 1 wave or invalid scoring.**
- Full CES-D-10 (.68 / .86 / .70). The two positively worded items had corrected item–total r of .00 / .28 / −.36.
- Full PSS-10 (.70 / .79 / .60).
- Family cohesion (.59 / .81 / .62).
- Family conflict (.27 / .34 / .46).
- Machismo (.83 / .68 / .69).
- Caballerismo (.63 / .60 / .62).
- MACVS Independence (.76 / .52 / .86).
- SRPB spirituality: at T7, 59% of scores = 1, and r with religiosity ≈ 0.
- NCES cohesion: correlation with depression flips across waves (.09 / .47 / −.12).
- The T7 immigration-stress instrument differs from T5/T6, so it can't be compared across waves.

Full table: `output/reliability_by_wave.csv`.

---

## 2. How the five models were chosen (transparent search)

The search had four blocks. Every test is logged in `output/screen_*.csv`.

- **A:** 1,810 pooled three-wave GEE moderation tests.
- **B:** 112 wave-interaction tests.
- **C:** 678 three-wave (T5→T6→T7) mediation tests.
- **D:** 1,810 prospective T5→T7 lagged moderation tests.

A model was retained only if it met all six criteria:
1. It is anchored in an established theory.
2. It is not in the RLIS portfolio you sent: published papers, Cabrera Tineo projects, medical-student projects, or Castillo abstracts.
3. It is significant both **unadjusted and adjusted**.
4. It survives **Benjamini–Hochberg FDR** correction within its search block.
5. It **replicates within each wave separately** (moderation models) or across the three-wave temporal sequence (mediation).
6. It is robust to log-transforming the skewed outcomes.

---

## MODEL 1 — Religiosity buffers the stress → depressive-symptom (and anxiety) link

**Type:** moderation (three-wave GEE, with replication at each wave).

**Theory.** Transactional stress-and-coping theory (Lazarus & Folkman, 1984) and Pargament's (1997) religious coping theory: religious meaning systems change how stressors are appraised, which weakens their emotional impact. The stress-buffering hypothesis (Cohen & Wills, 1985) predicts a negative stress × resource interaction.

**Gap / novelty.**
- In RLIS, religion has only been studied as religious coping → acculturative stress (Sanchez et al., 2012) and religious coping → acculturative stress and alcohol (Sanchez et al., 2015), plus religious social capital (Sanchez et al., 2019).
- None of these tested religiosity as a moderator of *general* perceived stress on depressive/anxiety symptoms, or did so across three waves.
- Outside RLIS, the evidence is mediation-only and in Spain ([Revollo et al., 2010](https://pubmed.ncbi.nlm.nih.gov/20619614/)), or cross-sectional ([HCHS/SOL, PLOS One 2017](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0185661)).
- The only RLIS-cohort stress → depression moderation paper tests emotion regulation and mindfulness, not religiosity ([PubMed 38015128](https://pubmed.ncbi.nlm.nih.gov/38015128/)).

**H1.** Religiosity will weaken the positive association between perceived stress and (a) depressive symptoms and (b) anxiety.

**1a. Outcome: depressive symptoms** (pooled GEE, 1,508 / 1,415 observations; 518 / 485 persons)

| Term | Unadjusted B (SE) | p | Adjusted B (SE) | p | Adjusted 95% CI |
|---|---|---|---|---|---|
| Perceived stress | .538 (.027) | < .001 | .521 (.028) | < .001 | .467, .575 |
| Religiosity | −.038 (.020) | .056 | −.047 (.022) | .034 | −.090, −.004 |
| **Stress × Religiosity** | **−.164 (.018)** | **< .001** | **−.177 (.018)** | **< .001** | −.212, −.142 |

**Simple slopes of stress, adjusted:**
- Low religiosity (−1 SD): β = .70 [.64, .76]
- Mean religiosity: β = .52
- High religiosity (+1 SD): β = .34 [.28, .41]

The stress → depression slope is **51% weaker** at high religiosity.

**1b. Outcome: anxiety.**
- Stress × Religiosity: unadjusted B = −.282 (SE .023), p < .001; adjusted B = −.289 (SE .024), p < .001.
- Simple slopes of stress: β = .74 at low religiosity vs .16 at high religiosity.

**Per-wave HC3 replication** (interaction B, adjusted):

| Outcome | T5 | T6 | T7 |
|---|---|---|---|
| Depression | −.157 (p = .006) | −.093 (p < .001) | −.179 (p < .001) |
| Anxiety | −.301 (p < .001) | −.234 (p < .001) | −.273 (p < .001) |

**Robustness.**
- FDR q < 10⁻¹⁹.
- Log-transformed outcomes: p < 10⁻¹⁵.

**Suggested title.** *"Religiosity as a Longitudinal Stress Buffer: Perceived Stress, Depressive and Anxiety Symptoms Among Latino/a Immigrants Across Three Years."*

---

## MODEL 2 — Emotional support moderates the Americanism → anxiety link

**Type:** moderation (three-wave GEE, with replication at each wave).

**Theory.**
- Berry's (1997) acculturation framework: adopting the receiving culture (Americanism) brings acculturative demands, including new norms, language pressure and identity negotiation. Those demands become psychologically costly when coping resources are thin.
- The stress-buffering hypothesis and the social convoy model: emotional/informational support is the resource that turns acculturative demands into manageable challenges.

**Gap / novelty.**
- RLIS papers have used Americanism only as a predictor of resilience (Cabrera Tineo; women only), of physical activity (ecological paper), or of alcohol (Westrick et al., 2021).
- Wuyke et al. (2023) examined *Hispanicism* as a moderator of context of reception → depression.
- No study has tested *emotional support × Americanism* on *anxiety*. The literature search found only main-effect or acculturative-stress–buffering studies (e.g., [Acculturative stress & resilience, PMC10299814](https://pmc.ncbi.nlm.nih.gov/articles/PMC10299814/)).

**H2.** Americanism will be positively associated with anxiety only when emotional support is low.

| Term | Unadjusted B (SE) | p | Adjusted B (SE) | p | Adjusted 95% CI |
|---|---|---|---|---|---|
| Americanism | .091 (.024) | < .001 | .084 (.025) | .001 | .035, .134 |
| Emotional support | −.286 (.041) | < .001 | −.307 (.044) | < .001 | −.393, −.220 |
| **Americanism × Emotional support** | **−.150 (.034)** | **< .001** | **−.165 (.035)** | **< .001** | −.233, −.097 |

**Simple slopes of Americanism, adjusted:**
- Low emotional support: β = .25 [.16, .34], p < .001
- Mean support: β = .08, p = .001
- High support: β = −.08 [−.16, −.01], p = .035

**Johnson–Neyman:** Americanism → anxiety is significantly positive when emotional support is below z = 0.22 (63.5% of observations). It is significantly negative above z = 0.95 (26.3%).

**Per-wave HC3 replication** (interaction B, adjusted): T5 −.170 (p = .044); T6 −.159 (p = .048); T7 −.164 (p < .001).

**Robustness.** FDR q < 10⁻⁴; log-transformed anxiety p < 10⁻⁵.

**Suggested title.** *"When Americanization Becomes Anxiety-Provoking: Emotional Support as a Moderator Among Latino/a Immigrants Over Three Years."*

---

## MODEL 3 — Hispanicism × Wave: from protective to risk factor (wave interaction)

**Type:** wave interaction (three-wave GEE).

**Theory.**
- Acculturation as a *time-dependent* process (Schwartz, Unger, Zamboanga & Szapocznik, 2010).
- Person–environment (cultural) fit: heritage-culture retention protects at arrival because it supplies co-ethnic networks and familiar meaning. As years in the U.S. accumulate, exclusive heritage orientation becomes a mismatch with host-society demands, restricting opportunity and increasing strain.
- This is the mechanism proposed for erosion of the immigrant health paradox.

**Gap / novelty.**
- RLIS treats Hispanicism as a static predictor or moderator (Wuyke et al., 2023; resilience and physical-activity projects).
- No RLIS or external study has tested whether the *effect* of Hispanicism on mental health *changes direction across waves*. The search found cross-sectional and adolescent studies only (e.g., [PMC11864884](https://pmc.ncbi.nlm.nih.gov/articles/PMC11864884/)).

**H3.** The association of Hispanicism with distress will change across waves: protective at T5, attenuating or reversing by T7.

**3a. Depressive symptoms** (pooled GEE)

| Term | Unadjusted B (SE) | p | Adjusted B (SE) | p |
|---|---|---|---|---|
| Hispanicism (effect at T5) | −.304 (.055) | < .001 | −.270 (.058) | < .001 |
| Hispanicism × T6 | .293 (.074) | < .001 | .227 (.077) | .003 |
| Hispanicism × T7 | .639 (.082) | < .001 | .609 (.085) | < .001 |
| **Omnibus Hispanicism × Wave, Wald χ²(2)** | **60.54** | **< .001** | **52.10** | **< .001** |

**Wave-specific Hispanicism effects, adjusted:**
- T5: β = −.27 [−.38, −.16], p < .001
- T6: β = −.04, p = .38
- T7: β = **+.34** [.21, .47], p < .001

**Replicated on other outcomes (adjusted):**
- **3b. Anxiety:** Wald χ²(2) = 31.86, p < .001. Hispanicism effect: T5 −.26 → T6 +.07 → T7 +.23.
- **3c. Perceived stress:** Wald χ²(2) = 27.40, p < .001. Hispanicism effect: T5 −.33 → T6 −.05 → T7 +.05.

**3d. Mechanism, stress × Hispanicism × Wave → depression (adjusted):** Wald χ²(2) = 48.62, p < .001.
- T5: Hispanicism *buffers* stress (interaction −.167, p < .001).
- T6: Hispanicism *amplifies* stress (+.116, p < .001).
- T7: amplification grows (+.308, p < .001). The stress slope is 1.03 at high Hispanicism vs .41 at low.

**Robustness.**
- FDR q < 10⁻¹⁰ (depression) and < 10⁻⁵ (anxiety, stress).
- Log-transformed outcomes: p < 10⁻¹⁰.

**Suggested title.** *"The Shifting Mental-Health Role of Heritage Culture Retention: A Three-Wave Study of Hispanicism Among Latino/a Immigrants."*

---

## MODEL 4 — Emotional-support buffering of stress → anxiety *emerges* over time (stress × support × wave)

**Type:** moderation with a wave interaction (three-way GEE).

**Theory.**
- The stress-buffering hypothesis (Cohen & Wills, 1985) combined with the social-convoy model (Kahn & Antonucci, 1980) and the social-support deterioration/deterrence model (Kaniasty & Norris, 1993).
- Migration disrupts support networks. Buffering capacity should grow as convoys are rebuilt, predicting that the stress × support interaction becomes stronger at later waves.

**Gap / novelty.**
- In the RLIS portfolio, Castillo Abstract 5 models *main effects* of social support and stress on anxiety at the latest wave, and Cano et al. (2016) used social support as a moderator of immigration stress → *alcohol*.
- No study has tested the stress × support **buffering interaction on anxiety** or **whether buffering changes across waves**.

**H4.** Emotional support will increasingly buffer the stress → anxiety association across waves.

| Wave | Stress × Emotional support B (95% CI), adjusted | p | Stress slope at low support | Stress slope at high support |
|---|---|---|---|---|
| T5 | −.075 (−.191, .042) | .21 | .41 | .26 |
| T6 | −.166 (−.274, −.058) | .003 | .62 | .29 |
| T7 | **−.446 (−.577, −.316)** | **< .001** | **1.16** | **.27** |

- **Omnibus three-way Wald χ²(2):** unadjusted 22.55, p < .001; adjusted 21.11, p < .001.
- **Pooled two-way (averaged over waves):** B = −.186 (unadjusted, p < .001); −.181 (adjusted, p < .001).
- **Per-wave HC3, adjusted:** T5 −.091 (p = .36); T6 −.144 (p = .005); T7 −.295 (p < .001).
- **Robustness:** log-transformed anxiety three-way p < 10⁻⁵.

> **Partial-overlap flag.** Castillo Abstract 5 uses the same constructs (social support, stress, anxiety) as *main effects*. Frame this paper around the **interaction and its emergence over time**, and cite Abstract 5. Clear it with Dr. Castillo before submission.

**Suggested title.** *"Rebuilding the Convoy: Emotional Support Increasingly Buffers Stress-Related Anxiety Among Latino/a Immigrants Across Three Years."*

---

## MODEL 5 — Neighborhood safety (T5) → social support (T6) → perceived stress (T7)

**Type:** three-wave (half-longitudinal) mediation with temporal precedence. X, M and Y are measured at successive waves, with prior M and Y controlled (Cole & Maxwell, 2003).

**Theory.**
- Social disorganization theory and collective efficacy (Sampson, Raudenbush & Earls, 1997).
- Neighborhood stress process (Aneshensel, 1992).
- Safe neighborhoods allow outdoor contact, trust and network formation, which build perceived support. Support in turn lowers stress appraisals.

**Gap / novelty.**
- RLIS papers have not examined neighborhood safety as a longitudinal antecedent of social support. Castillo Abstract 4 has T7 group differences; Castillo Abstract 5 has collective efficacy as a main effect.
- Externally, similar paths exist cross-sectionally or in rural mothers ([PubMed 41396502](https://pubmed.ncbi.nlm.nih.gov/41396502/); [Kwon, 2020](https://onlinelibrary.wiley.com/doi/10.1002/jcop.22288)), but not in a three-wave adult immigrant cohort.

**H5.** Greater T5 neighborhood safety will predict lower T7 perceived stress indirectly, through higher T6 social support.

**Primary outcome: T7 perceived stress** (PROCESS Model 4, HC3, 5,000 bootstraps; controls: T5 social support and T5 stress)

| Path | Unadjusted (n = 491) B (SE) | p | Adjusted (n = 462) B (SE) | p |
|---|---|---|---|---|
| a: Safety T5 → Support T6 | .116 (.052) | .026 | .166 (.057) | .004 |
| b: Support T6 → Stress T7 | −.140 (.051) | .006 | −.173 (.053) | .001 |
| c′: direct | .079 (.052) | .131 | .051 (.057) | .372 |
| c: total | .062 (.051) | .223 | .022 (.057) | .694 |
| **Indirect a × b [95% boot CI]** | **−.016 [−.037, −.001]** | sig | **−.029 [−.057, −.008]** | sig |

**Secondary outcomes (adjusted only):**
- T7 depression: indirect −.019 [−.045, −.001]
- T7 anxiety: indirect −.020 [−.046, −.002]

Both are non-significant unadjusted, so report them as secondary.

The total effect is not significant. Report this as an *indirect-only* (full) mediation (Zhao, Lynch & Chen, 2010), which is accepted practice.

**Suggested title.** *"Safe Streets, Stronger Ties: A Three-Wave Mediation of Neighborhood Safety, Social Support, and Perceived Stress Among Latino/a Immigrants."*

---

## 3. Novelty check against the RLIS portfolio you provided

| Model | Closest existing RLIS work | What is new |
|---|---|---|
| 1 Stress × Religiosity → Dep/Anx | Sanchez 2012, 2015 (religious coping → acculturative stress/alcohol) | Religiosity as moderator of general stress; mental-health outcomes; three waves |
| 2 Americanism × Emotional support → Anx | Westrick 2021 (alcohol); resilience project | Acculturation × support interaction; anxiety outcome |
| 3 Hispanicism × Wave | Wuyke 2023 (Hispanicism moderator, men, cross-sectional) | Time-varying *reversal* of the Hispanicism effect across 3 waves |
| 4 Stress × Emotional support × Wave → Anx | Castillo Abstract 5 (main effects); Cano 2016 (alcohol) | Buffering interaction and its temporal emergence (**flag**) |
| 5 Safety → Support → Stress | Castillo Abstracts 4 and 5 (cross-sectional, main effects) | Temporal three-wave mediation |

None of the five uses fatalism, alcohol outcomes (avoiding the Cabrera Tineo sex × depression → alcohol paper and the AOE papers), marianismo or machismo (avoiding the resilience and physical-activity projects), physical activity or diet.

---

## 4. Limitations to state up front (reviewers will check these)

1. **Low rank-order stability.** One-year test–retest correlations of the psychosocial scales are low (r ≈ .05–.28; AUDIT .56–.72). Record linkage was verified: age r = .96–.98 across waves, and gender is consistent. The pooled moderation models therefore capture **replicated concurrent (within-wave) associations**, not prospective prediction. Word them as "across three waves, the association between X and Y was consistently weaker when…". Avoid "X predicts later Y".
2. **Prospective lagged moderation (T5 → T7) was null for these interactions** (24 / 905 significant, which is chance level). Say so honestly. Model 5 is the model that uses temporal ordering.
3. **Skewed outcomes.** GAD-7 skew is ≈ 3; all results replicate with ln-transformed outcomes. Report this as a sensitivity analysis.
4. **Documentation status is highly unbalanced** (28 undocumented). It is kept as a covariate but not tested as a moderator.
5. **GEE vs HC3.** SPSS GENLIN robust (sandwich) SEs are the correct HC-type estimator for clustered data. HC3 is applied wherever observations are independent (per-wave models and mediation).
6. **Attrition was minimal.** 25 of 518 participants (4.8%) were lost by T7. Those lost did not differ on T5 depression, anxiety, stress, support, religiosity or acculturation (all p > .05). They were more often South American (50% vs 27%) and less often Central American, and they reported slightly lower T5 neighborhood safety (p = .04). GEE uses all available person-waves. See `output/final/Table_attrition.csv`.
7. **Exploratory search.** Disclose that candidate models were screened and that retained effects survived FDR correction and per-wave replication.

---

## 5. Files

| File | Contents |
|---|---|
| `spss_syntax/RLIS_Five_Longitudinal_Models.sps` | Complete SPSS syntax: scale building, reliability, merge, restructure, z-scores, GENLIN, PROCESS |
| `scripts/01_build_dataset.py` | Scale construction and reliability (Python mirror of SPSS STEP 1–4) |
| `scripts/stats_engine.py` | GEE, OLS-HC3, PROCESS Models 1 and 4 re-implementations |
| `scripts/02_screen.py` | Full transparent model screen |
| `scripts/03_final_models.py` | Final five models (unadjusted and adjusted) |
| `output/final/*` | Coefficient tables, simple slopes, J-N regions, per-wave HC3, indirect effects, descriptives, attrition, correlations |

Expected SPSS-vs-Python agreement: coefficients and SEs should match to about the 2nd–3rd decimal. Small differences come from GENLIN's scale-parameter and working-correlation estimation, and from the bootstrap seed, which moves CI limits by about ±.002.
