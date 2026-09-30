# RLIS T5–T6–T7: Five Longitudinal Models with AUDIT (Hazardous Drinking) as the Outcome

**Outcome.** AUDIT total score (0–40), which matches the stored `AUDIT_SS` / `audit_ss` exactly (r = 1.00).

| | T5 | T6 | T7 |
|---|---|---|---|
| α | .78 | .71 | .83 |
| M (SD) | 3.99 (3.76) | 3.54 (2.95) | 4.20 (3.97) |
| AUDIT ≥ 8 | 16.0% | 8.9% | 15.8% |

**Sample.** 518 persons, 1,554 person-waves.

**Covariates (adjusted models).** Gender, age, household income (ln), education, marital status, employment, documentation status, South American origin, Central American origin (Cuban = reference). All continuous variables are z-standardized.

**Estimation.**
- **Primary:** three-wave GEE (SPSS GENLIN: normal/identity on z-AUDIT, exchangeable working correlation, robust sandwich SEs). This is the cluster-robust analogue of HC3, required because each person contributes three rows.
- **Per-wave replication:** OLS with **HC3** (PROCESS `/hc=3`; UNIANOVA `/ROBUST=HC3`).
- **Sensitivity 1:** negative-binomial GEE on raw AUDIT (log link, SPSS `NEGBIN(1)`). Exp(B) is an incidence-rate ratio (IRR).
- **Sensitivity 2:** linear GEE on ln(AUDIT + 1). AUDIT is right-skewed (skew ≈ 1.3–1.7) with many zeros (13–15%).

**Fatalism was not used.**

---

## Predictor measures (α ≥ .70 at every wave)

| Construct | Instrument (items) | α T5 / T6 / T7 |
|---|---|---|
| Anxiety symptoms | GAD-7 (7) | .89 / .88 / .93 |
| Perceived stress | PSS perceived-helplessness factor (6) | .79 / .81 / .76 |
| Emotional/informational support | MOS-SSS (8) | .96 / .97 / .96 |
| Familism support values | MACVS Familism-Support (6) | .88 / .89 / .81 |
| Depressive symptoms (secondary only) | CES-D negative-affect items (8) | .80 / .86 / .89 |
| Married/partnered | Baseline binary | — |

---

## How the models were chosen

Every test is logged in `output/audit/`.

| Block | Tests | Result |
|---|---|---|
| A. Three-wave GEE moderation | 704 | 58 / 352 adjusted significant (16%, ≈ 3× chance) |
| B. Predictor × wave | 46 | Strongest block |
| C. Three-wave mediation, X T5 → M T6 → AUDIT T7 (controlling M T5 and AUDIT T5) | 500 | **0 significant adjusted** |
| D. Prospective lagged moderation, T5 → T7 | 738 | Chance level |
| E. Within-wave mediation replicated at T5/T6/T7, plus pooled | 500 | None significant both unadjusted and adjusted |

**Why mediation failed.** AUDIT is highly stable across waves (r = .56–.72), while the psychosocial predictors are not (r ≈ .05–.28). After controlling prior AUDIT there is little change left to explain.

I therefore selected **moderation and wave-interaction models**, as your instructions allowed. A model was retained only if it:
1. rests on an established theory;
2. does not duplicate your RLIS portfolio;
3. is significant **unadjusted and adjusted** in the primary GEE;
4. holds in the **negative-binomial and log-AUDIT** sensitivity models;
5. shows the same direction when estimated at each wave separately.

---

## MODEL 1 — Familism support values buffer the anxiety → hazardous-drinking link

**Theory.**
- Tension-reduction theory and the self-medication hypothesis: anxious people drink to relieve negative affect.
- Cultural-values-as-protective-resource perspective (familismo). Endorsing family closeness, loyalty and mutual support supplies alternative coping resources and normative constraints on drinking-to-cope.

**Gap.**
- RLIS papers studied *pre-immigration* familismo → alcohol misuse (Dillon et al., 2013) and family *cohesion* (a family-functioning behavior) → AUDIT (Cano et al., 2017; Castillo abstracts).
- No RLIS paper tests familism *values* as a **moderator of anxiety** → AUDIT. Cabrera Tineo's paper uses *depression* × *sex*.
- The only external test, in college students, found no moderation ([PubMed 40138516](https://pubmed.ncbi.nlm.nih.gov/40138516/)). This would be the first adult-immigrant, three-wave test.

**H1.** Familism support values weaken the positive association between anxiety and hazardous drinking.

| Term | Unadjusted B (SE) | p | Adjusted B (SE) | p | Adjusted 95% CI |
|---|---|---|---|---|---|
| Anxiety | .043 (.032) | .181 | .050 (.036) | .163 | −.020, .120 |
| Familism support | −.032 (.021) | .127 | −.043 (.022) | .048 | −.086, −.000 |
| **Anxiety × Familism support** | **−.115 (.039)** | **.003** | **−.127 (.044)** | **.004** | −.214, −.040 |

**Simple slopes of anxiety, adjusted:**
- Low familism (−1 SD): β = **+.18** [.03, .32], p = .018
- Mean familism: β = .05, ns
- High familism (+1 SD): β = **−.08** [−.14, −.02], p = .010

**Johnson–Neyman:** anxiety → AUDIT is significantly positive when familism support is below z = −0.37 (33.6% of observations). It is significantly negative above z = 0.78 (42.0%).

**Sensitivity:**
- Negative binomial: IRR = 0.925, p = .0006 (adjusted); IRR = 0.938, p = .004 (unadjusted).
- ln(AUDIT + 1): p = .0003 (adjusted); p = .0004 (unadjusted).

**Per-wave HC3 interaction B, adjusted:** T5 −.178 (p = .004); T6 −.094 (p = .070); T7 −.230 (p = .13). The sign is the same at every wave.

**FDR:** q = .039 (unadjusted) and .077 (adjusted) within the distress × buffer family.

**Secondary analysis:** depression × familism support → AUDIT. Unadjusted B = −.082 (p = .030); adjusted B = −.090 (p = .030); NB p = .013; log p = .038. Present it as a supplementary analysis, citing Cabrera Tineo Paper 1.

*Suggested title:* "Familismo as a Cultural Buffer: Anxiety and Hazardous Drinking Among Latino/a Immigrants Across Three Years."

---

## MODEL 2 — Being partnered buffers the stress → hazardous-drinking link

**Theory.**
- Stress-and-coping and stress-buffering (Cohen & Wills, 1985).
- Social control theory of marriage (Umberson, 1987): a partner regulates health behavior and supplies daily support, weakening drinking-to-cope under stress.

**Gap.**
- Marital status appears in RLIS only as a *covariate* or *main effect*. Cabrera Tineo found partnership protective for binge drinking.
- Moderators tested in Cano et al. (2016) were gender, immigration status and social support, with *immigration* stress as the predictor.
- The one Hispanic-immigrant AUDIT–stress study found no stress effect and did not test partnership ([PMC9738040](https://pmc.ncbi.nlm.nih.gov/articles/PMC9738040/)).

**H2.** Perceived stress is associated with hazardous drinking among unpartnered, but not partnered, immigrants.

| Term | Unadjusted B (SE) | p | Adjusted B (SE) | p | Adjusted 95% CI |
|---|---|---|---|---|---|
| Perceived stress (among unpartnered) | .196 (.060) | .001 | .182 (.062) | .003 | .061, .304 |
| Married/partnered | −.228 (.087) | .009 | −.156 (.085) | .068 | −.322, .011 |
| **Stress × Partnered** | **−.218 (.062)** | **< .001** | **−.209 (.064)** | **.001** | −.334, −.085 |

**Simple slopes of stress, adjusted:**
- Unpartnered: β = **.18** [.06, .30], p = .003
- Partnered: β = −.03 [−.07, .01], p = .18

**Sensitivity:**
- Negative binomial: IRR = 0.850, p = .0004 (adjusted); IRR = 0.853, p < .0001 (unadjusted).
- ln(AUDIT + 1): p = .007 (adjusted); p = .003 (unadjusted).

**Per-wave HC3 interaction B:**
- Unadjusted: T5 −.343 (p = .008); T6 −.269 (p = .020); T7 −.073 (ns).
- Adjusted: T5 −.375 (p = .006); T6 −.194 (p = .075); T7 −.041 (ns).

The buffering is strongest in the earlier waves.

**FDR:** q = .015 (unadjusted) and .035 (adjusted). Significant.

*Suggested title:* "Partnership as a Stress Buffer: Perceived Stress and Hazardous Drinking Among Latino/a Immigrants."

---

## MODEL 3 — Anxiety × Wave: the anxiety → hazardous-drinking link fades over time

**Theory.**
- Tension-reduction / self-medication.
- Acculturation-as-process: as immigrants settle, coping repertoires widen (services, networks, routines), and drinking becomes less of a response to anxiety.

**Gap.**
- The Cabrera Tineo paper models *depression* × *sex*, with wave only as time.
- No RLIS or external study tests whether the **anxiety–drinking association itself changes across waves** in adult immigrants.

**H3.** The anxiety–AUDIT association will weaken across waves.

| Term | Unadjusted B (SE) | p | Adjusted B (SE) | p |
|---|---|---|---|---|
| Anxiety (effect at T5) | .223 (.070) | .002 | .236 (.077) | .002 |
| Anxiety × T6 | −.198 (.078) | .011 | −.216 (.082) | .008 |
| Anxiety × T7 | −.294 (.071) | < .001 | −.318 (.081) | < .001 |
| **Omnibus Anxiety × Wave, Wald χ²(2)** | **17.19** | **< .001** | **15.42** | **< .001** |

**Wave-specific anxiety effect, adjusted:**
- T5: β = **.24** [.09, .39], p = .002
- T6: β = .02, ns
- T7: β = −.08, p = .08

**Sensitivity:**
- Negative binomial: Wald p = .039 (adjusted), .019 (unadjusted). IRR = 1.14 at T5 → 0.92 at T7.
- ln(AUDIT + 1): p < .0001.

**Per-wave HC3 simple regression, adjusted:** T5 B = .299 (p = .003); T6 .104 (p = .07); T7 −.063 (ns).

**FDR:** q = .002 (unadjusted) and .003 (adjusted).

*Suggested title:* "Does Anxiety-Related Drinking Fade With Settlement? A Three-Wave Study of Latino/a Immigrants."

---

## MODEL 4 — Emotional support × Wave: protective at T5, not by T6–T7

**Theory.**
- Social-network / social-influence theory of drinking (Rosenquist et al., 2010).
- Social-convoy model. Early on, emotional support is embedded in family/co-ethnic networks that constrain drinking. As immigrants integrate, support networks increasingly include peers and coworkers whose drinking norms are more permissive ([assimilation and alcohol evidence](https://pubmed.ncbi.nlm.nih.gov/40237107)). The protective association therefore disappears.

**Gap.**
- RLIS papers estimate social support → alcohol as a single *time-invariant* effect: Cano et al. (2017); Castillo Abstracts 3 and 5, where support predicts fewer binge days.
- None tests whether the support–drinking association **changes across waves**.

> **Overlap flag:** same constructs as Cano 2017 and Castillo Abstract 3. The contribution is the *time-varying* effect. Cite both and confirm with Dr. Castillo.

**H4.** Emotional support will be protective at T5 and weaken over time.

| Term | Unadjusted B (SE) | p | Adjusted B (SE) | p |
|---|---|---|---|---|
| Emotional support (effect at T5) | −.123 (.041) | .002 | −.121 (.042) | .004 |
| Support × T6 | .145 (.048) | .003 | .175 (.049) | < .001 |
| Support × T7 | .195 (.050) | < .001 | .198 (.053) | < .001 |
| **Omnibus Support × Wave, Wald χ²(2)** | **16.85** | **< .001** | **18.75** | **< .001** |

**Wave-specific support effect, adjusted:**
- T5: β = **−.12** [−.20, −.04], p = .004
- T6: β = +.05, p = .050
- T7: β = +.08, p = .082

**Sensitivity:** negative binomial p = .006 (adjusted), .003 (unadjusted); ln(AUDIT + 1) p < .001.

**Per-wave HC3:**
- Unadjusted: T5 −.270 (p < .001); T6 −.115 (p = .005); T7 +.039 (ns).
- Adjusted: T5 −.166 (p = .016); T6 .000; T7 +.058.

**FDR:** q = .002.

*Suggested title:* "When Support Stops Protecting: The Changing Link Between Emotional Support and Hazardous Drinking Among Latino/a Immigrants."

---

## MODEL 5 — Perceived stress × Wave: the stress → hazardous-drinking link reverses by T7

**Theory.**
- Stress-response-dampening and drinking-to-cope models, which are distinct from anxiety-based tension reduction.
- Life-course and role-accumulation perspectives: later in settlement, stress increasingly comes from work and family *role demands*, which constrain drinking time and opportunity rather than encouraging it.

**Gap.** Perceived stress → AUDIT has not been tested in RLIS as a time-varying association. Existing RLIS work uses *immigration* or *acculturative* stress, and treats the effect as constant.

**H5.** The stress–AUDIT association will change across waves.

| Term | Unadjusted B (SE) | p | Adjusted B (SE) | p |
|---|---|---|---|---|
| Stress (effect at T5) | .094 (.040) | .021 | .084 (.042) | .044 |
| Stress × T6 | −.018 (.046) | .689 | −.007 (.048) | .881 |
| Stress × T7 | −.178 (.059) | .003 | −.184 (.062) | .003 |
| **Omnibus Stress × Wave, Wald χ²(2)** | **11.14** | **.004** | **11.73** | **.003** |

**Wave-specific stress effect, adjusted:**
- T5: β = .08, p = .044
- T6: β = .08, p = .004
- T7: β = **−.10**, p = .054

With negative binomial, the T7 effect is IRR = 0.887, p = .020.

**Sensitivity:** negative binomial p = .001 (adjusted), .002 (unadjusted); ln(AUDIT + 1) p = .0003.

**FDR:** q = .016–.022.

> Models 3 and 5 both concern internalizing distress and drinking over time. They can be written as two papers (anxiety/tension-reduction vs. stress/role-demand frameworks) or merged into one stronger paper. **Alternate Model 5** if you prefer a different construct: ethnic identity × Hispanicism → AUDIT (adjusted B = .047, p = .003; NB p = .0004). This is weaker: unadjusted log-model p = .07, and it replicates only at T7.

---

## Novelty check against your RLIS portfolio

| Model | Closest existing work | Distinct contribution |
|---|---|---|
| 1 Anxiety × Familism support → AUDIT | Dillon 2013 (pre-immigration familismo); Cano 2017 (family cohesion); Cabrera Tineo (depression × sex) | Familism *values* as moderator of *anxiety*-related drinking; three waves |
| 2 Stress × Partnered → AUDIT | Cabrera Tineo (married → binge, main effect); Cano 2016 (immigration stress; gender/status/support moderators) | Partnership as stress buffer; perceived (general) stress |
| 3 Anxiety × Wave → AUDIT | Cabrera Tineo (depression → AUDIT) | Time-varying anxiety–drinking coupling |
| 4 Emotional support × Wave → AUDIT | Cano 2017; Castillo Abstracts 3 and 5 (**flag**) | Time-varying direction of the support effect |
| 5 Stress × Wave → AUDIT | Cano 2016 (immigration stress) | Time-varying general-stress–drinking association |

None of the five uses fatalism, alcohol expectancies (AOE papers; Cabrera Tineo 2026), machismo/caballerismo, marianismo, neighborhood → support → drinking (Kanamori et al., 2016, [PMC5129265](https://pmc.ncbi.nlm.nih.gov/articles/PMC5129265/)), or immigration stress.

## Limitations to disclose

1. **Moderation effects are concurrent (within-wave) associations** that repeat across three waves, not prospective prediction. The prospective T5 → T7 tests were at chance level. Avoid causal wording ("predicts later").
2. **The Model 1 and 2 interactions are strongest at T5–T6** and weaker at T7 (same direction). Report the per-wave table.
3. **Effect sizes are small** (interaction β ≈ .1–.2 SD). This is typical for AUDIT moderation studies; report IRRs from the NB models for clinical interpretation.
4. **No mediation model met your criteria** (0 of 500 adjusted three-wave mediations). The reason is AUDIT's high stability combined with low stability of the predictors.
5. **Exploratory search with full logs.** Disclose the screen and the FDR results:
   - Models 2–5 survive FDR within their families.
   - Model 1 is borderline (adjusted q = .077) but holds in every sensitivity analysis.
6. **Attrition was 4.8%** by T7 (see `output/final/Table_attrition.csv`).

## Files

| File | Contents |
|---|---|
| `spss_syntax/RLIS_AUDIT_Five_Longitudinal_Models.sps` | Complete SPSS syntax (GENLIN primary + NEGBIN sensitivity, PROCESS/UNIANOVA HC3 per wave) |
| `scripts/04_audit_screen.py`, `scripts/05_audit_wave_mediation_screen.py` | Full transparent AUDIT screens |
| `scripts/06_audit_final_models.py` | Final five models, unadjusted and adjusted, with sensitivity analyses |
| `output/audit/final/*` | Coefficients (normal and NB), simple slopes, J-N, wave-specific effects, per-wave HC3, sensitivity tables, descriptives, correlations |
