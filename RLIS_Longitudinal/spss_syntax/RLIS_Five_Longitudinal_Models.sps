* Encoding: UTF-8.
* =============================================================================.
* RLIS T5-T6-T7: FIVE NOVEL LONGITUDINAL MODELS  (SPSS 27+; PROCESS v4.x macro).
* Reproduces scripts/01_build_dataset.py + scripts/03_final_models.py step by step.
* Run blocks IN ORDER.  Edit the FILE HANDLE paths in STEP 0 first.
* Requirements: PROCESS v4 macro (Hayes) defined via  Extensions > Utilities > Custom Dialog
*   or by running process.sps once per session (INSERT FILE='...\process.sps').
* =============================================================================.

* -----------------------------------------------------------------------------.
* STEP 0. File locations (EDIT THESE).
* -----------------------------------------------------------------------------.
FILE HANDLE T5  /NAME='C:\RLIS\RLIS_T5_Longitudinal_Model_WORKING.sav'.
FILE HANDLE T6  /NAME='C:\RLIS\RLIS_T6_Longitudinal_Model_WORKING.sav'.
FILE HANDLE T7  /NAME='C:\RLIS\RLIS_T7_Longitudinal_Model_WORKING1.sav'.
FILE HANDLE MST /NAME='C:\RLIS\RLIS_T5T6T7_Master_Longitudinal_ANALYSIS.sav'.
FILE HANDLE OUT /NAME='C:\RLIS\analysis'.

* -----------------------------------------------------------------------------.
* STEP 1. Build scales identically in each wave.
*   Macro !BUILD w=5|6|7 : cleans out-of-range codes (-88,-99, 93 ...), reverse-scores,
*   computes item-mean scores, and runs RELIABILITY (Cronbach alpha) for each scale.
* -----------------------------------------------------------------------------.
DEFINE !BUILD (w = !TOKENS(1))

* --- set impossible / missing codes to system-missing (item response ranges) ---.
RECODE ds_bothered ds_focus ds_depressed ds_effort ds_fearful ds_restless ds_lonely ds_ngoing
   (LOWEST THRU 0.5=SYSMIS) (4.5 THRU HIGHEST=SYSMIS).
RECODE gad_nervous gad_nstopworry gad_worry gad_nrelax gad_restless gad_irritable gad_afraid
   (LOWEST THRU -0.5=SYSMIS) (3.5 THRU HIGHEST=SYSMIS).
RECODE pss_unexpected pss_ncontrol pss_nervous pss_ncope pss_angered pss_difficulties
   (LOWEST THRU -0.5=SYSMIS) (4.5 THRU HIGHEST=SYSMIS).
RECODE ssi_listen ssi_info ssi_advice ssi_confide ssi_want ssi_share ssi_suggestions ssi_understands
   ssi_bed ssi_doctor ssi_meals ssi_chores ssi_affection ssi_wanted ssi_hugs ssi_goodtime ssi_relax ssi_enjoyable ssi_mindoff
   (LOWEST THRU 0.5=SYSMIS) (5.5 THRU HIGHEST=SYSMIS).
RECODE macv_strength macv_god1 macv_faith macv_word macv_religion1 macv_thank macv_teachpray
   (LOWEST THRU 0.5=SYSMIS) (5.5 THRU HIGHEST=SYSMIS).
RECODE biq_spanhome biq_spanschool biq_spanwork biq_spanfriend biq_spangen biq_hispmusic
   biq_hispdance biq_hispplaces biq_hisprec biq_hisptv biq_hipradio biq_hispbookmag
   biq_enghome biq_engschool biq_engwork biq_engfriend biq_enggen biq_amermusic biq_amerdance
   biq_amerplaces biq_amerrec biq_amertv biq_amerradio biq_amerbookmag
   (LOWEST THRU -0.5=SYSMIS) (4.5 THRU HIGHEST=SYSMIS).
RECODE ns_safewalk ns_nviolence ns_safe (LOWEST THRU 0.5=SYSMIS) (5.5 THRU HIGHEST=SYSMIS).
EXECUTE.

* --- Neighborhood safety: raw 1 = strongly agree, so reverse all 3 items (higher = safer) ---.
COMPUTE ns_safewalk_rv  = 6 - ns_safewalk.
COMPUTE ns_nviolence_rv = 6 - ns_nviolence.
COMPUTE ns_safe_rv      = 6 - ns_safe.

* --- Scale scores (item means; minimum items required as in the Python build) ---.
* Depressive symptoms: 8 negative-affect CES-D items (1-4).
COMPUTE DEP!w    = MEAN.6(ds_bothered, ds_focus, ds_depressed, ds_effort, ds_fearful, ds_restless, ds_lonely, ds_ngoing).
* Anxiety: GAD-7 (0-3).
COMPUTE ANX!w    = MEAN.6(gad_nervous, gad_nstopworry, gad_worry, gad_nrelax, gad_restless, gad_irritable, gad_afraid).
* Perceived stress: PSS perceived-helplessness factor, 6 items (0-4).
COMPUTE PSS!w    = MEAN.5(pss_unexpected, pss_ncontrol, pss_nervous, pss_ncope, pss_angered, pss_difficulties).
* MOS social support: total (19 items) and emotional/informational (8 items) (1-5).
COMPUTE SOCSUP!w = MEAN.15(ssi_listen, ssi_info, ssi_advice, ssi_confide, ssi_want, ssi_share, ssi_suggestions,
                   ssi_understands, ssi_bed, ssi_doctor, ssi_meals, ssi_chores, ssi_affection, ssi_wanted, ssi_hugs,
                   ssi_goodtime, ssi_relax, ssi_enjoyable, ssi_mindoff).
COMPUTE EMOSUP!w = MEAN.6(ssi_listen, ssi_info, ssi_advice, ssi_confide, ssi_want, ssi_share, ssi_suggestions, ssi_understands).
* Religiosity: MACVS religion subscale, 7 items (1-5).
COMPUTE RELIG!w  = MEAN.6(macv_strength, macv_god1, macv_faith, macv_word, macv_religion1, macv_thank, macv_teachpray).
* Hispanicism / Americanism: BIQ, 12 items each (0-4).
COMPUTE HISP!w   = MEAN.10(biq_spanhome, biq_spanschool, biq_spanwork, biq_spanfriend, biq_spangen, biq_hispmusic,
                   biq_hispdance, biq_hispplaces, biq_hisprec, biq_hisptv, biq_hipradio, biq_hispbookmag).
COMPUTE AMER!w   = MEAN.10(biq_enghome, biq_engschool, biq_engwork, biq_engfriend, biq_enggen, biq_amermusic,
                   biq_amerdance, biq_amerplaces, biq_amerrec, biq_amertv, biq_amerradio, biq_amerbookmag).
* Neighborhood safety (1-5, higher = safer).
COMPUTE NSAFE!w  = MEAN.2(ns_safewalk_rv, ns_nviolence_rv, ns_safe_rv).
EXECUTE.

* --- Internal consistency (report alpha per wave in Table 1) ---.
RELIABILITY /VARIABLES=ds_bothered ds_focus ds_depressed ds_effort ds_fearful ds_restless ds_lonely ds_ngoing
  /SCALE('DEP') ALL /MODEL=ALPHA /SUMMARY=TOTAL.
RELIABILITY /VARIABLES=gad_nervous gad_nstopworry gad_worry gad_nrelax gad_restless gad_irritable gad_afraid
  /SCALE('ANX') ALL /MODEL=ALPHA /SUMMARY=TOTAL.
RELIABILITY /VARIABLES=pss_unexpected pss_ncontrol pss_nervous pss_ncope pss_angered pss_difficulties
  /SCALE('PSS') ALL /MODEL=ALPHA /SUMMARY=TOTAL.
RELIABILITY /VARIABLES=ssi_listen ssi_info ssi_advice ssi_confide ssi_want ssi_share ssi_suggestions ssi_understands
  ssi_bed ssi_doctor ssi_meals ssi_chores ssi_affection ssi_wanted ssi_hugs ssi_goodtime ssi_relax ssi_enjoyable ssi_mindoff
  /SCALE('SOCSUP') ALL /MODEL=ALPHA.
RELIABILITY /VARIABLES=ssi_listen ssi_info ssi_advice ssi_confide ssi_want ssi_share ssi_suggestions ssi_understands
  /SCALE('EMOSUP') ALL /MODEL=ALPHA.
RELIABILITY /VARIABLES=macv_strength macv_god1 macv_faith macv_word macv_religion1 macv_thank macv_teachpray
  /SCALE('RELIG') ALL /MODEL=ALPHA.
RELIABILITY /VARIABLES=biq_spanhome biq_spanschool biq_spanwork biq_spanfriend biq_spangen biq_hispmusic
  biq_hispdance biq_hispplaces biq_hisprec biq_hisptv biq_hipradio biq_hispbookmag
  /SCALE('HISP') ALL /MODEL=ALPHA.
RELIABILITY /VARIABLES=biq_enghome biq_engschool biq_engwork biq_engfriend biq_enggen biq_amermusic
  biq_amerdance biq_amerplaces biq_amerrec biq_amertv biq_amerradio biq_amerbookmag
  /SCALE('AMER') ALL /MODEL=ALPHA.
RELIABILITY /VARIABLES=ns_safewalk_rv ns_nviolence_rv ns_safe_rv /SCALE('NSAFE') ALL /MODEL=ALPHA.

SORT CASES BY caseid.
SAVE OUTFILE=!QUOTE(!CONCAT(!UNQUOTE('OUT/scales_T'), !w, !UNQUOTE('.sav')))
  /KEEP=caseid DEP!w ANX!w PSS!w SOCSUP!w EMOSUP!w RELIG!w HISP!w AMER!w NSAFE!w.
!ENDDEFINE.

GET FILE=T5.
!BUILD w=5.
GET FILE=T6.
!BUILD w=6.
GET FILE=T7.
!BUILD w=7.

* -----------------------------------------------------------------------------.
* STEP 2. Baseline covariates (from master + T5) and WIDE file.
*   gender 1=male; maristat 1=married/partnered; DOCSTAT 1=documented;
*   EDUHI 1=more than high school; EMPLOY 1=employed; ORIG_SA/ORIG_CA dummies (Cuban = ref).
*   AGE = T5 age; LNHINC5 = ln(household income + 1)  (already in the T5 file).
* -----------------------------------------------------------------------------.
GET FILE=T5 /KEEP=caseid demo_age LNHINC5.
RECODE demo_age (LOWEST THRU 0=SYSMIS).
RENAME VARIABLES (demo_age=AGE) (LNHINC5=LNHINC).
SORT CASES BY caseid.
SAVE OUTFILE='OUT/t5cov.sav'.

GET FILE=MST /KEEP=caseid gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA.
SORT CASES BY caseid.
MATCH FILES /TABLE=* /FILE='OUT/t5cov.sav' /BY caseid.
SAVE OUTFILE='OUT/cov.sav'.

MATCH FILES /FILE='OUT/scales_T5.sav' /FILE='OUT/scales_T6.sav' /FILE='OUT/scales_T7.sav'
  /TABLE='OUT/cov.sav' /BY caseid.
SELECT IF NOT MISSING(DEP5) OR NOT MISSING(DEP6) OR NOT MISSING(DEP7).
EXECUTE.
SAVE OUTFILE='OUT/RLIS_wide.sav'.

* -----------------------------------------------------------------------------.
* STEP 3. Descriptives, attrition and correlations (Table 1, Table 2).
* -----------------------------------------------------------------------------.
FREQUENCIES gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA.
DESCRIPTIVES AGE LNHINC DEP5 DEP6 DEP7 ANX5 ANX6 ANX7 PSS5 PSS6 PSS7 SOCSUP5 SOCSUP6 SOCSUP7
  EMOSUP5 EMOSUP6 EMOSUP7 RELIG5 RELIG6 RELIG7 HISP5 HISP6 HISP7 AMER5 AMER6 AMER7 NSAFE5 NSAFE6 NSAFE7
  /STATISTICS=MEAN STDDEV MIN MAX SKEWNESS KURTOSIS.
* Attrition: retained at T7 vs lost.
COMPUTE retainT7 = NOT MISSING(DEP7).
T-TEST GROUPS=retainT7(0 1) /VARIABLES=DEP5 ANX5 PSS5 SOCSUP5 EMOSUP5 RELIG5 HISP5 AMER5 NSAFE5 AGE LNHINC gender.

* -----------------------------------------------------------------------------.
* STEP 4. Restructure to LONG (person-period) file and standardize.
*   WAVE: 0 = T5, 1 = T6, 2 = T7.  Grand-mean z-scores (DESCRIPTIVES /SAVE) in the long file.
* -----------------------------------------------------------------------------.
GET FILE='OUT/RLIS_wide.sav'.
VARSTOCASES
  /MAKE DEP FROM DEP5 DEP6 DEP7
  /MAKE ANX FROM ANX5 ANX6 ANX7
  /MAKE PSS FROM PSS5 PSS6 PSS7
  /MAKE SOCSUP FROM SOCSUP5 SOCSUP6 SOCSUP7
  /MAKE EMOSUP FROM EMOSUP5 EMOSUP6 EMOSUP7
  /MAKE RELIG FROM RELIG5 RELIG6 RELIG7
  /MAKE HISP FROM HISP5 HISP6 HISP7
  /MAKE AMER FROM AMER5 AMER6 AMER7
  /MAKE NSAFE FROM NSAFE5 NSAFE6 NSAFE7
  /INDEX=WAVE(3)
  /KEEP=caseid gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA AGE LNHINC
  /NULL=KEEP.
COMPUTE WAVE = WAVE - 1.
VALUE LABELS WAVE 0 'T5' 1 'T6' 2 'T7'.
COMPUTE W6 = (WAVE = 1).
COMPUTE W7 = (WAVE = 2).
COMPUTE W5 = (WAVE = 0).
EXECUTE.
DESCRIPTIVES DEP ANX PSS SOCSUP EMOSUP RELIG HISP AMER NSAFE AGE LNHINC /SAVE.
* Creates zDEP zANX zPSS zSOCSUP zEMOSUP zRELIG zHISP zAMER zNSAFE zAGE zLNHINC.

* Product terms (computed so that GENLIN, PROCESS and simple-slope re-centering all use the same terms).
COMPUTE PSSxREL   = zPSS * zRELIG.
COMPUTE AMERxEMO  = zAMER * zEMOSUP.
COMPUTE PSSxEMO   = zPSS * zEMOSUP.
COMPUTE HISPxW6   = zHISP * W6.
COMPUTE HISPxW7   = zHISP * W7.
COMPUTE PSSxHISP  = zPSS * zHISP.
COMPUTE PSSxW6    = zPSS * W6.
COMPUTE PSSxW7    = zPSS * W7.
COMPUTE HISPW6    = zHISP * W6.
COMPUTE HISPW7    = zHISP * W7.
COMPUTE PHxW6     = zPSS * zHISP * W6.
COMPUTE PHxW7     = zPSS * zHISP * W7.
COMPUTE EMOxW6    = zEMOSUP * W6.
COMPUTE EMOxW7    = zEMOSUP * W7.
COMPUTE PExW6     = zPSS * zEMOSUP * W6.
COMPUTE PExW7     = zPSS * zEMOSUP * W7.
EXECUTE.
SORT CASES BY caseid WAVE.
SAVE OUTFILE='OUT/RLIS_long.sav'.

* -----------------------------------------------------------------------------.
* GEE MACRO  (SPSS GENLIN = Analyze > Generalized Linear Models > Generalized Estimating Equations)
*   Normal / identity, subject = caseid, within-subject = WAVE, exchangeable working correlation,
*   ROBUST (sandwich) covariance.  !cov is empty for UNADJUSTED models.
* -----------------------------------------------------------------------------.
DEFINE !GEE (dv = !TOKENS(1) / terms = !ENCLOSE('[',']') / cov = !DEFAULT() !ENCLOSE('[',']'))
GENLIN !dv WITH !terms !cov
  /MODEL !terms !cov INTERCEPT=YES DISTRIBUTION=NORMAL LINK=IDENTITY
  /CRITERIA SCALE=MLE PCONVERGE=1E-006(ABSOLUTE) SINGULAR=1E-012 ANALYSISTYPE=3(WALD) CILEVEL=95 LIKELIHOOD=FULL
  /REPEATED SUBJECT=caseid WITHINSUBJECT=WAVE SORT=YES CORRTYPE=EXCHANGEABLE ADJUSTCORR=NO
     COVB=ROBUST MAXITERATIONS=100 PCONVERGE=1E-006(ABSOLUTE) UPDATECORR=1
  /MISSING CLASSMISSING=EXCLUDE
  /PRINT CPS DESCRIPTIVES MODELINFO FIT SUMMARY SOLUTION WORKINGCORR.
!ENDDEFINE.

* Adjusted models add: gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC.

* =============================================================================.
* MODEL 1. Religiosity x Perceived stress -> Depressive symptoms (1a) and Anxiety (1b).
* =============================================================================.
GET FILE='OUT/RLIS_long.sav'.
* 1a unadjusted / adjusted.
!GEE dv=zDEP terms=[zPSS zRELIG PSSxREL W6 W7].
!GEE dv=zDEP terms=[zPSS zRELIG PSSxREL W6 W7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
* 1b unadjusted / adjusted.
!GEE dv=zANX terms=[zPSS zRELIG PSSxREL W6 W7].
!GEE dv=zANX terms=[zPSS zRELIG PSSxREL W6 W7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].

* Simple slopes of stress at LOW (-1 SD) and HIGH (+1 SD) religiosity (re-centering method):
*   the coefficient of zPSS in each run below IS the simple slope at that moderator level.
COMPUTE RELlo = zRELIG + 1.
COMPUTE RELhi = zRELIG - 1.
COMPUTE PSSxRELlo = zPSS * RELlo.
COMPUTE PSSxRELhi = zPSS * RELhi.
EXECUTE.
!GEE dv=zDEP terms=[zPSS RELlo PSSxRELlo W6 W7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
!GEE dv=zDEP terms=[zPSS RELhi PSSxRELhi W6 W7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
!GEE dv=zANX terms=[zPSS RELlo PSSxRELlo W6 W7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
!GEE dv=zANX terms=[zPSS RELhi PSSxRELhi W6 W7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].

* Wave-by-wave replication with HC3 (PROCESS Model 1). Within-wave z-scores via SPLIT FILE.
SORT CASES BY WAVE caseid.
SPLIT FILE LAYERED BY WAVE.
DESCRIPTIVES VARIABLES=DEP (wzDEP) ANX (wzANX) PSS (wzPSS) RELIG (wzRELIG) AMER (wzAMER)
  EMOSUP (wzEMOSUP) HISP (wzHISP) AGE (wzAGE) LNHINC (wzLNHINC) /SAVE.
* (wz* = z-scores computed WITHIN each wave, for the per-wave HC3 replications).
SPLIT FILE OFF.
TEMPORARY.
SELECT IF WAVE = 0.
PROCESS y=wzDEP /x=wzPSS /w=wzRELIG /cov=gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /model=1 /hc=3 /jn=1 /moments=1 /plot=1.
TEMPORARY.
SELECT IF WAVE = 1.
PROCESS y=wzDEP /x=wzPSS /w=wzRELIG /cov=gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /model=1 /hc=3 /jn=1 /moments=1 /plot=1.
TEMPORARY.
SELECT IF WAVE = 2.
PROCESS y=wzDEP /x=wzPSS /w=wzRELIG /cov=gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /model=1 /hc=3 /jn=1 /moments=1 /plot=1.
* Unadjusted per-wave: delete the /cov= line.  Repeat with y=wzANX for Model 1b.
SAVE OUTFILE='OUT/RLIS_long.sav'.

* =============================================================================.
* MODEL 2. Emotional support x Americanism -> Anxiety.
* =============================================================================.
!GEE dv=zANX terms=[zAMER zEMOSUP AMERxEMO W6 W7].
!GEE dv=zANX terms=[zAMER zEMOSUP AMERxEMO W6 W7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
COMPUTE EMOlo = zEMOSUP + 1.
COMPUTE EMOhi = zEMOSUP - 1.
COMPUTE AMERxEMOlo = zAMER * EMOlo.
COMPUTE AMERxEMOhi = zAMER * EMOhi.
EXECUTE.
!GEE dv=zANX terms=[zAMER EMOlo AMERxEMOlo W6 W7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
!GEE dv=zANX terms=[zAMER EMOhi AMERxEMOhi W6 W7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
* Per-wave HC3 replication.
TEMPORARY.
SELECT IF WAVE = 0.
PROCESS y=wzANX /x=wzAMER /w=wzEMOSUP /cov=gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /model=1 /hc=3 /jn=1 /moments=1 /plot=1.
TEMPORARY.
SELECT IF WAVE = 1.
PROCESS y=wzANX /x=wzAMER /w=wzEMOSUP /cov=gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /model=1 /hc=3 /jn=1 /moments=1 /plot=1.
TEMPORARY.
SELECT IF WAVE = 2.
PROCESS y=wzANX /x=wzAMER /w=wzEMOSUP /cov=gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /model=1 /hc=3 /jn=1 /moments=1 /plot=1.

* =============================================================================.
* MODEL 3. Hispanicism x Wave -> Depressive symptoms (3a), Anxiety (3b), Perceived stress (3c).
*   Omnibus Hispanicism x Wave test = joint Wald test of HISPxW6 and HISPxW7:
*   request it with the /TEST subcommand is not available in GENLIN, so use the
*   Type III (ANALYSISTYPE=3) table after re-running with WAVE as a FACTOR (second call).
* =============================================================================.
!GEE dv=zDEP terms=[zHISP W6 W7 HISPxW6 HISPxW7].
!GEE dv=zDEP terms=[zHISP W6 W7 HISPxW6 HISPxW7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
* Omnibus Wald chi-square (2 df) for HISP x WAVE (see "Tests of Model Effects" row WAVE*zHISP).
GENLIN zDEP BY WAVE (ORDER=ASCENDING) WITH zHISP gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
  /MODEL WAVE zHISP WAVE*zHISP gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
     INTERCEPT=YES DISTRIBUTION=NORMAL LINK=IDENTITY
  /CRITERIA SCALE=MLE ANALYSISTYPE=3(WALD) CILEVEL=95
  /REPEATED SUBJECT=caseid WITHINSUBJECT=WAVE SORT=YES CORRTYPE=EXCHANGEABLE ADJUSTCORR=NO COVB=ROBUST
  /PRINT MODELINFO FIT SUMMARY SOLUTION.
* Wave-specific Hispanicism effects: the zHISP coefficient is the T5 effect in the first call.
*   For T6 / T7 effects, change the reference wave (coefficient of zHISP = effect at the reference wave):.
COMPUTE HISPxW5 = zHISP * W5.
EXECUTE.
!GEE dv=zDEP terms=[zHISP W5 W7 HISPxW5 HISPxW7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
!GEE dv=zDEP terms=[zHISP W5 W6 HISPxW5 HISPxW6] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
* Repeat the four !GEE calls above with dv=zANX (3b) and dv=zPSS (3c).

* 3d. Three-way: Perceived stress x Hispanicism x Wave -> Depressive symptoms.
!GEE dv=zDEP terms=[zPSS zHISP PSSxHISP W6 W7 PSSxW6 PSSxW7 HISPW6 HISPW7 PHxW6 PHxW7].
!GEE dv=zDEP terms=[zPSS zHISP PSSxHISP W6 W7 PSSxW6 PSSxW7 HISPW6 HISPW7 PHxW6 PHxW7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
* PSSxHISP = stress x Hispanicism interaction at T5; PSSxHISP + PHxW6 = at T6; PSSxHISP + PHxW7 = at T7
*   (obtain T6/T7 directly by re-running with W5 and W7 [or W5 and W6] as the wave dummies, as above).

* =============================================================================.
* MODEL 4. Perceived stress x Emotional support x Wave -> Anxiety (emergent stress-buffering).
* =============================================================================.
!GEE dv=zANX terms=[zPSS zEMOSUP PSSxEMO W6 W7 PSSxW6 PSSxW7 EMOxW6 EMOxW7 PExW6 PExW7].
!GEE dv=zANX terms=[zPSS zEMOSUP PSSxEMO W6 W7 PSSxW6 PSSxW7 EMOxW6 EMOxW7 PExW6 PExW7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
* Pooled two-way (averaged over waves):.
!GEE dv=zANX terms=[zPSS zEMOSUP PSSxEMO W6 W7].
!GEE dv=zANX terms=[zPSS zEMOSUP PSSxEMO W6 W7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
* Per-wave HC3 replication (buffering absent at T5, present at T6, strongest at T7).
TEMPORARY.
SELECT IF WAVE = 0.
PROCESS y=wzANX /x=wzPSS /w=wzEMOSUP /cov=gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /model=1 /hc=3 /jn=1 /moments=1 /plot=1.
TEMPORARY.
SELECT IF WAVE = 1.
PROCESS y=wzANX /x=wzPSS /w=wzEMOSUP /cov=gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /model=1 /hc=3 /jn=1 /moments=1 /plot=1.
TEMPORARY.
SELECT IF WAVE = 2.
PROCESS y=wzANX /x=wzPSS /w=wzEMOSUP /cov=gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /model=1 /hc=3 /jn=1 /moments=1 /plot=1.

* =============================================================================.
* MODEL 5. Three-wave mediation (WIDE file):
*   Neighborhood safety (T5) -> Social support (T6) -> Perceived stress (T7),
*   controlling prior levels of the mediator (SOCSUP5) and outcome (PSS5).
*   PROCESS Model 4, HC3 SEs, 5,000 percentile bootstrap resamples.
* =============================================================================.
GET FILE='OUT/RLIS_wide.sav'.
DESCRIPTIVES NSAFE5 SOCSUP5 SOCSUP6 PSS5 PSS7 DEP5 DEP7 ANX5 ANX7 AGE LNHINC /SAVE.
* Unadjusted (autoregressive controls only).
PROCESS y=zPSS7 /x=zNSAFE5 /m=zSOCSUP6 /cov=zSOCSUP5 zPSS5
  /model=4 /hc=3 /boot=5000 /seed=20260924 /total=1 /effsize=1.
* Adjusted.
PROCESS y=zPSS7 /x=zNSAFE5 /m=zSOCSUP6 /cov=zSOCSUP5 zPSS5 gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
  /model=4 /hc=3 /boot=5000 /seed=20260924 /total=1 /effsize=1.
* Secondary outcomes (report as secondary; significant when adjusted only).
PROCESS y=zDEP7 /x=zNSAFE5 /m=zSOCSUP6 /cov=zSOCSUP5 zDEP5 gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
  /model=4 /hc=3 /boot=5000 /seed=20260924 /total=1.
PROCESS y=zANX7 /x=zNSAFE5 /m=zSOCSUP6 /cov=zSOCSUP5 zANX5 gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
  /model=4 /hc=3 /boot=5000 /seed=20260924 /total=1.
* END OF SYNTAX.
