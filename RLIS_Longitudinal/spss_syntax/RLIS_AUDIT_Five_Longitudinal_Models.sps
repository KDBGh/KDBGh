* Encoding: UTF-8.
* =============================================================================.
* RLIS T5-T6-T7: FIVE LONGITUDINAL MODELS WITH AUDIT (HAZARDOUS DRINKING) AS OUTCOME.
* Reproduces scripts/06_audit_final_models.py.  SPSS 27+, PROCESS v4.x macro loaded.
* Run blocks IN ORDER.  Edit the FILE HANDLE paths in STEP 0 first.
* =============================================================================.

* STEP 0. File locations (EDIT).
FILE HANDLE T5  /NAME='C:\RLIS\RLIS_T5_Longitudinal_Model_WORKING.sav'.
FILE HANDLE T6  /NAME='C:\RLIS\RLIS_T6_Longitudinal_Model_WORKING.sav'.
FILE HANDLE T7  /NAME='C:\RLIS\RLIS_T7_Longitudinal_Model_WORKING1.sav'.
FILE HANDLE MST /NAME='C:\RLIS\RLIS_T5T6T7_Master_Longitudinal_ANALYSIS.sav'.
FILE HANDLE OUT /NAME='C:\RLIS\analysis_audit'.

* -----------------------------------------------------------------------------.
* STEP 1. Build scales identically at each wave (item means) + Cronbach alpha.
*   AUDIT: the stored total score is used (AUDIT_SS at T5/T6, audit_ss at T7). It was verified to be
*   identical (r = 1.00) to a fresh 10-item AUDIT score (items 1-8 = 0-4; items 9-10 = 0/2/4;
*   non-drinkers scored 0 on items 2-10).
* -----------------------------------------------------------------------------.
DEFINE !BUILD (w = !TOKENS(1) / aud = !TOKENS(1))
RECODE ds_bothered ds_focus ds_depressed ds_effort ds_fearful ds_restless ds_lonely ds_ngoing
   (LOWEST THRU 0.5=SYSMIS) (4.5 THRU HIGHEST=SYSMIS).
RECODE gad_nervous gad_nstopworry gad_worry gad_nrelax gad_restless gad_irritable gad_afraid
   (LOWEST THRU -0.5=SYSMIS) (3.5 THRU HIGHEST=SYSMIS).
RECODE pss_unexpected pss_ncontrol pss_nervous pss_ncope pss_angered pss_difficulties
   (LOWEST THRU -0.5=SYSMIS) (4.5 THRU HIGHEST=SYSMIS).
RECODE ssi_listen ssi_info ssi_advice ssi_confide ssi_want ssi_share ssi_suggestions ssi_understands
   (LOWEST THRU 0.5=SYSMIS) (5.5 THRU HIGHEST=SYSMIS).
RECODE macv_first macv_security macv_united macv_close macv_holiday macv_affect
   (LOWEST THRU 0.5=SYSMIS) (5.5 THRU HIGHEST=SYSMIS).
EXECUTE.
* Anxiety: GAD-7 (0-3).
COMPUTE ANX!w    = MEAN.6(gad_nervous, gad_nstopworry, gad_worry, gad_nrelax, gad_restless, gad_irritable, gad_afraid).
* Depressive symptoms: CES-D negative-affect items (1-4) [secondary analysis only].
COMPUTE DEP!w    = MEAN.6(ds_bothered, ds_focus, ds_depressed, ds_effort, ds_fearful, ds_restless, ds_lonely, ds_ngoing).
* Perceived stress: PSS perceived-helplessness factor (0-4).
COMPUTE PSS!w    = MEAN.5(pss_unexpected, pss_ncontrol, pss_nervous, pss_ncope, pss_angered, pss_difficulties).
* Emotional/informational support: MOS-SSS (1-5).
COMPUTE EMOSUP!w = MEAN.6(ssi_listen, ssi_info, ssi_advice, ssi_confide, ssi_want, ssi_share, ssi_suggestions, ssi_understands).
* Familism support values: MACVS (1-5).
COMPUTE FAMSUP!w = MEAN.5(macv_first, macv_security, macv_united, macv_close, macv_holiday, macv_affect).
COMPUTE AUDIT!w  = !aud.
EXECUTE.
RELIABILITY /VARIABLES=gad_nervous gad_nstopworry gad_worry gad_nrelax gad_restless gad_irritable gad_afraid
  /SCALE('ANX') ALL /MODEL=ALPHA /SUMMARY=TOTAL.
RELIABILITY /VARIABLES=pss_unexpected pss_ncontrol pss_nervous pss_ncope pss_angered pss_difficulties
  /SCALE('PSS') ALL /MODEL=ALPHA /SUMMARY=TOTAL.
RELIABILITY /VARIABLES=ssi_listen ssi_info ssi_advice ssi_confide ssi_want ssi_share ssi_suggestions ssi_understands
  /SCALE('EMOSUP') ALL /MODEL=ALPHA.
RELIABILITY /VARIABLES=macv_first macv_security macv_united macv_close macv_holiday macv_affect
  /SCALE('FAMSUP') ALL /MODEL=ALPHA.
RELIABILITY /VARIABLES=ds_bothered ds_focus ds_depressed ds_effort ds_fearful ds_restless ds_lonely ds_ngoing
  /SCALE('DEP') ALL /MODEL=ALPHA.
SORT CASES BY caseid.
SAVE OUTFILE=!QUOTE(!CONCAT(!UNQUOTE('OUT/audit_scales_T'), !w, !UNQUOTE('.sav')))
  /KEEP=caseid ANX!w DEP!w PSS!w EMOSUP!w FAMSUP!w AUDIT!w.
!ENDDEFINE.

GET FILE=T5.
!BUILD w=5 aud=AUDIT_SS.
* AUDIT item alpha T5.
RELIABILITY /VARIABLES=frqdrink numdrink frq5more nostopdk failnorm drinkam guiltdrk unabremb injurdrk concerdk
  /SCALE('AUDIT T5') ALL /MODEL=ALPHA.
GET FILE=T6.
!BUILD w=6 aud=AUDIT_SS.
RELIABILITY /VARIABLES=audit_often audit_howmany audit_5ormore audit_cantstop audit_didntdo audit_morningdrink
  audit_guilt audit_cantremember audit_injured audit_dctrconcern /SCALE('AUDIT T6') ALL /MODEL=ALPHA.
GET FILE=T7.
!BUILD w=7 aud=audit_ss.
RELIABILITY /VARIABLES=audit_often audit_howmany audit_5ormore audit_cantstop audit_didntdo audit_morningdrink
  audit_guilt audit_cantremember audit_injured audit_dctrconcern /SCALE('AUDIT T7') ALL /MODEL=ALPHA.
* NOTE: for AUDIT alpha, non-drinkers have system-missing items 2-10; to match the reported alpha,
*   first RECODE those items (SYSMIS=0) for cases with frequency item = 0.

* -----------------------------------------------------------------------------.
* STEP 2. Baseline covariates + WIDE file.
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
MATCH FILES /FILE='OUT/audit_scales_T5.sav' /FILE='OUT/audit_scales_T6.sav' /FILE='OUT/audit_scales_T7.sav'
  /TABLE='OUT/cov.sav' /BY caseid.
EXECUTE.
SAVE OUTFILE='OUT/RLIS_AUDIT_wide.sav'.
DESCRIPTIVES AUDIT5 AUDIT6 AUDIT7 ANX5 ANX6 ANX7 PSS5 PSS6 PSS7 EMOSUP5 EMOSUP6 EMOSUP7 FAMSUP5 FAMSUP6 FAMSUP7
  /STATISTICS=MEAN STDDEV MIN MAX SKEWNESS KURTOSIS.
CORRELATIONS /VARIABLES=AUDIT5 ANX5 PSS5 EMOSUP5 FAMSUP5.
CORRELATIONS /VARIABLES=AUDIT6 ANX6 PSS6 EMOSUP6 FAMSUP6.
CORRELATIONS /VARIABLES=AUDIT7 ANX7 PSS7 EMOSUP7 FAMSUP7.

* -----------------------------------------------------------------------------.
* STEP 3. LONG file, dummies, grand-mean z-scores, product terms.
* -----------------------------------------------------------------------------.
VARSTOCASES
  /MAKE AUDIT FROM AUDIT5 AUDIT6 AUDIT7
  /MAKE ANX FROM ANX5 ANX6 ANX7
  /MAKE DEP FROM DEP5 DEP6 DEP7
  /MAKE PSS FROM PSS5 PSS6 PSS7
  /MAKE EMOSUP FROM EMOSUP5 EMOSUP6 EMOSUP7
  /MAKE FAMSUP FROM FAMSUP5 FAMSUP6 FAMSUP7
  /INDEX=WAVE(3)
  /KEEP=caseid gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA AGE LNHINC
  /NULL=KEEP.
COMPUTE WAVE = WAVE - 1.
VALUE LABELS WAVE 0 'T5' 1 'T6' 2 'T7'.
COMPUTE W5 = (WAVE = 0).
COMPUTE W6 = (WAVE = 1).
COMPUTE W7 = (WAVE = 2).
COMPUTE lnAUDIT = LN(AUDIT + 1).
EXECUTE.
DESCRIPTIVES AUDIT lnAUDIT ANX DEP PSS EMOSUP FAMSUP AGE LNHINC /SAVE.
COMPUTE ANXxFAM  = zANX * zFAMSUP.
COMPUTE DEPxFAM  = zDEP * zFAMSUP.
COMPUTE PSSxMAR  = zPSS * maristat.
COMPUTE ANXxW6   = zANX * W6.
COMPUTE ANXxW7   = zANX * W7.
COMPUTE ANXxW5   = zANX * W5.
COMPUTE PSSxW6   = zPSS * W6.
COMPUTE PSSxW7   = zPSS * W7.
COMPUTE PSSxW5   = zPSS * W5.
COMPUTE EMOxW6   = zEMOSUP * W6.
COMPUTE EMOxW7   = zEMOSUP * W7.
COMPUTE EMOxW5   = zEMOSUP * W5.
EXECUTE.
SORT CASES BY caseid WAVE.
SAVE OUTFILE='OUT/RLIS_AUDIT_long.sav'.

* -----------------------------------------------------------------------------.
* GEE macros.  !GEE  = normal/identity on z-scored DV (PRIMARY).
*              !GEENB = negative binomial (ancillary 1) / log link on raw AUDIT (SENSITIVITY; Exp(B) = IRR).
*   Adjusted models: add  gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC  to cov=[ ].
*   (For Model 2 maristat is the moderator, so it is not repeated as a covariate.)
* -----------------------------------------------------------------------------.
DEFINE !GEE (dv = !TOKENS(1) / terms = !ENCLOSE('[',']') / cov = !DEFAULT() !ENCLOSE('[',']'))
GENLIN !dv WITH !terms !cov
  /MODEL !terms !cov INTERCEPT=YES DISTRIBUTION=NORMAL LINK=IDENTITY
  /CRITERIA SCALE=MLE ANALYSISTYPE=3(WALD) CILEVEL=95 LIKELIHOOD=FULL
  /REPEATED SUBJECT=caseid WITHINSUBJECT=WAVE SORT=YES CORRTYPE=EXCHANGEABLE ADJUSTCORR=NO COVB=ROBUST
     MAXITERATIONS=100 PCONVERGE=1E-006(ABSOLUTE) UPDATECORR=1
  /MISSING CLASSMISSING=EXCLUDE
  /PRINT CPS DESCRIPTIVES MODELINFO FIT SUMMARY SOLUTION WORKINGCORR.
!ENDDEFINE.
DEFINE !GEENB (terms = !ENCLOSE('[',']') / cov = !DEFAULT() !ENCLOSE('[',']'))
GENLIN AUDIT WITH !terms !cov
  /MODEL !terms !cov INTERCEPT=YES DISTRIBUTION=NEGBIN(1) LINK=LOG
  /CRITERIA ANALYSISTYPE=3(WALD) CILEVEL=95 LIKELIHOOD=FULL
  /REPEATED SUBJECT=caseid WITHINSUBJECT=WAVE SORT=YES CORRTYPE=EXCHANGEABLE ADJUSTCORR=NO COVB=ROBUST
     MAXITERATIONS=100 PCONVERGE=1E-006(ABSOLUTE) UPDATECORR=1
  /PRINT MODELINFO FIT SUMMARY SOLUTION (EXPONENTIATED).
!ENDDEFINE.

* =============================================================================.
* MODEL 1. Anxiety x Familism support values -> AUDIT.
* =============================================================================.
GET FILE='OUT/RLIS_AUDIT_long.sav'.
!GEE dv=zAUDIT terms=[zANX zFAMSUP ANXxFAM W6 W7].
!GEE dv=zAUDIT terms=[zANX zFAMSUP ANXxFAM W6 W7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
* Simple slopes of anxiety at low (-1 SD) / high (+1 SD) familism: coefficient of zANX in each run.
COMPUTE FAMlo = zFAMSUP + 1.
COMPUTE FAMhi = zFAMSUP - 1.
COMPUTE ANXxFAMlo = zANX * FAMlo.
COMPUTE ANXxFAMhi = zANX * FAMhi.
EXECUTE.
!GEE dv=zAUDIT terms=[zANX FAMlo ANXxFAMlo W6 W7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
!GEE dv=zAUDIT terms=[zANX FAMhi ANXxFAMhi W6 W7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
* Sensitivity: negative binomial and ln(AUDIT+1).
!GEENB terms=[zANX zFAMSUP ANXxFAM W6 W7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
!GEE dv=zlnAUDIT terms=[zANX zFAMSUP ANXxFAM W6 W7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
* Secondary: depressive symptoms x familism support.
!GEE dv=zAUDIT terms=[zDEP zFAMSUP DEPxFAM W6 W7].
!GEE dv=zAUDIT terms=[zDEP zFAMSUP DEPxFAM W6 W7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].

* Per-wave HC3 replication (PROCESS Model 1; within-wave z-scores).
SORT CASES BY WAVE caseid.
SPLIT FILE LAYERED BY WAVE.
DESCRIPTIVES VARIABLES=AUDIT (wzAUDIT) ANX (wzANX) DEP (wzDEP) PSS (wzPSS) EMOSUP (wzEMOSUP) FAMSUP (wzFAMSUP)
  AGE (wzAGE) LNHINC (wzLNHINC) /SAVE.
SPLIT FILE OFF.
SAVE OUTFILE='OUT/RLIS_AUDIT_long.sav'.
TEMPORARY.
SELECT IF WAVE = 0.
PROCESS y=wzAUDIT /x=wzANX /w=wzFAMSUP /cov=gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /model=1 /hc=3 /jn=1 /moments=1 /plot=1.
TEMPORARY.
SELECT IF WAVE = 1.
PROCESS y=wzAUDIT /x=wzANX /w=wzFAMSUP /cov=gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /model=1 /hc=3 /jn=1 /moments=1 /plot=1.
TEMPORARY.
SELECT IF WAVE = 2.
PROCESS y=wzAUDIT /x=wzANX /w=wzFAMSUP /cov=gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /model=1 /hc=3 /jn=1 /moments=1 /plot=1.
* Unadjusted per-wave: delete the /cov= list.

* =============================================================================.
* MODEL 2. Perceived stress x Married/partnered -> AUDIT.
* =============================================================================.
!GEE dv=zAUDIT terms=[zPSS maristat PSSxMAR W6 W7].
!GEE dv=zAUDIT terms=[zPSS maristat PSSxMAR W6 W7] cov=[gender DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
* Simple slope for PARTNERED participants: re-code the moderator so partnered = 0.
COMPUTE unpart = 1 - maristat.
COMPUTE PSSxUNP = zPSS * unpart.
EXECUTE.
!GEE dv=zAUDIT terms=[zPSS unpart PSSxUNP W6 W7] cov=[gender DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
* (zPSS coefficient: first adjusted run = slope among unpartnered; this run = slope among partnered.)
!GEENB terms=[zPSS maristat PSSxMAR W6 W7] cov=[gender DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
!GEE dv=zlnAUDIT terms=[zPSS maristat PSSxMAR W6 W7] cov=[gender DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
TEMPORARY.
SELECT IF WAVE = 0.
PROCESS y=wzAUDIT /x=wzPSS /w=maristat /cov=gender DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /model=1 /hc=3 /plot=1.
TEMPORARY.
SELECT IF WAVE = 1.
PROCESS y=wzAUDIT /x=wzPSS /w=maristat /cov=gender DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /model=1 /hc=3 /plot=1.
TEMPORARY.
SELECT IF WAVE = 2.
PROCESS y=wzAUDIT /x=wzPSS /w=maristat /cov=gender DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /model=1 /hc=3 /plot=1.

* =============================================================================.
* MODEL 3. Anxiety x Wave -> AUDIT (tension-reduction coupling over time).
* =============================================================================.
!GEE dv=zAUDIT terms=[zANX W6 W7 ANXxW6 ANXxW7].
!GEE dv=zAUDIT terms=[zANX W6 W7 ANXxW6 ANXxW7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
* Omnibus 2-df Wald test: Type III row "WAVE * zANX".
GENLIN zAUDIT BY WAVE (ORDER=ASCENDING) WITH zANX gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
  /MODEL WAVE zANX WAVE*zANX gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
     INTERCEPT=YES DISTRIBUTION=NORMAL LINK=IDENTITY
  /CRITERIA SCALE=MLE ANALYSISTYPE=3(WALD) CILEVEL=95
  /REPEATED SUBJECT=caseid WITHINSUBJECT=WAVE SORT=YES CORRTYPE=EXCHANGEABLE ADJUSTCORR=NO COVB=ROBUST
  /PRINT MODELINFO FIT SUMMARY SOLUTION.
* Wave-specific anxiety effects (zANX coefficient = effect at the omitted/reference wave):.
!GEE dv=zAUDIT terms=[zANX W5 W7 ANXxW5 ANXxW7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
!GEE dv=zAUDIT terms=[zANX W5 W6 ANXxW5 ANXxW6] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
!GEENB terms=[zANX W6 W7 ANXxW6 ANXxW7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
* Per-wave HC3: AUDIT regressed on anxiety at each wave.
* SPSS REGRESSION has no HC3 option, so UNIANOVA /ROBUST=HC3 is used (Univariate GLM > Options > robust SEs).
TEMPORARY.
SELECT IF WAVE = 0.
UNIANOVA wzAUDIT WITH wzANX gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /ROBUST=HC3 /PRINT=PARAMETER /DESIGN=wzANX gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC.
TEMPORARY.
SELECT IF WAVE = 1.
UNIANOVA wzAUDIT WITH wzANX gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /ROBUST=HC3 /PRINT=PARAMETER /DESIGN=wzANX gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC.
TEMPORARY.
SELECT IF WAVE = 2.
UNIANOVA wzAUDIT WITH wzANX gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC
  /ROBUST=HC3 /PRINT=PARAMETER /DESIGN=wzANX gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA wzAGE wzLNHINC.

* =============================================================================.
* MODEL 4. Emotional support x Wave -> AUDIT.
* =============================================================================.
!GEE dv=zAUDIT terms=[zEMOSUP W6 W7 EMOxW6 EMOxW7].
!GEE dv=zAUDIT terms=[zEMOSUP W6 W7 EMOxW6 EMOxW7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
!GEE dv=zAUDIT terms=[zEMOSUP W5 W7 EMOxW5 EMOxW7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
!GEE dv=zAUDIT terms=[zEMOSUP W5 W6 EMOxW5 EMOxW6] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
!GEENB terms=[zEMOSUP W6 W7 EMOxW6 EMOxW7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
* Omnibus: rerun the "GENLIN zAUDIT BY WAVE ..." block above replacing zANX by zEMOSUP.

* =============================================================================.
* MODEL 5. Perceived stress x Wave -> AUDIT.
* =============================================================================.
!GEE dv=zAUDIT terms=[zPSS W6 W7 PSSxW6 PSSxW7].
!GEE dv=zAUDIT terms=[zPSS W6 W7 PSSxW6 PSSxW7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
!GEE dv=zAUDIT terms=[zPSS W5 W7 PSSxW5 PSSxW7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
!GEE dv=zAUDIT terms=[zPSS W5 W6 PSSxW5 PSSxW6] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
!GEENB terms=[zPSS W6 W7 PSSxW6 PSSxW7] cov=[gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC].
* Omnibus: rerun the "GENLIN zAUDIT BY WAVE ..." block replacing zANX by zPSS.
* END.
