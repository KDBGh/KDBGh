* Encoding: UTF-8.
* =============================================================================.
* MODEL 1: Familism support values moderate the Anxiety -> Hazardous drinking (AUDIT) link.
* RLIS T5-T6-T7, three-wave GEE (robust SEs) + per-wave PROCESS Model 1 (HC3).
* Run PART A -> PART G in order.  Compare each CHECKPOINT with your output before moving on.
* =============================================================================.

* ---- EDIT THESE 5 PATHS (folder for OUT must already exist) ----.
FILE HANDLE T5  /NAME='C:\RLIS\RLIS_T5_Longitudinal_Model_WORKING.sav'.
FILE HANDLE T6  /NAME='C:\RLIS\RLIS_T6_Longitudinal_Model_WORKING.sav'.
FILE HANDLE T7  /NAME='C:\RLIS\RLIS_T7_Longitudinal_Model_WORKING1.sav'.
FILE HANDLE MST /NAME='C:\RLIS\RLIS_T5T6T7_Master_Longitudinal_ANALYSIS.sav'.
FILE HANDLE OUT /NAME='C:\RLIS\model1'.

* =============================================================================.
* PART A. Build Anxiety (GAD-7), Familism support (MACVS) and AUDIT at each wave.
* CHECKPOINT A (Cronbach alpha):  ANX  .892 / .879 / .933
*                                 FAMSUP .881 / .893 / .814
* =============================================================================.
DEFINE !M1BUILD (w = !TOKENS(1) / aud = !TOKENS(1))
RECODE gad_nervous gad_nstopworry gad_worry gad_nrelax gad_restless gad_irritable gad_afraid
   (LOWEST THRU -0.5=SYSMIS) (3.5 THRU HIGHEST=SYSMIS).
RECODE macv_first macv_security macv_united macv_close macv_holiday macv_affect
   (LOWEST THRU 0.5=SYSMIS) (5.5 THRU HIGHEST=SYSMIS).
EXECUTE.
COMPUTE ANX!w    = MEAN.6(gad_nervous, gad_nstopworry, gad_worry, gad_nrelax, gad_restless, gad_irritable, gad_afraid).
COMPUTE FAMSUP!w = MEAN.5(macv_first, macv_security, macv_united, macv_close, macv_holiday, macv_affect).
COMPUTE AUDIT!w  = !aud.
EXECUTE.
RELIABILITY /VARIABLES=gad_nervous gad_nstopworry gad_worry gad_nrelax gad_restless gad_irritable gad_afraid
  /SCALE('ANX') ALL /MODEL=ALPHA /SUMMARY=TOTAL.
RELIABILITY /VARIABLES=macv_first macv_security macv_united macv_close macv_holiday macv_affect
  /SCALE('FAMSUP') ALL /MODEL=ALPHA /SUMMARY=TOTAL.
SORT CASES BY caseid.
SAVE OUTFILE=!QUOTE(!CONCAT(!UNQUOTE('OUT/m1_T'), !w, !UNQUOTE('.sav'))) /KEEP=caseid ANX!w FAMSUP!w AUDIT!w.
!ENDDEFINE.

GET FILE=T5.
!M1BUILD w=5 aud=AUDIT_SS.
GET FILE=T6.
!M1BUILD w=6 aud=AUDIT_SS.
GET FILE=T7.
!M1BUILD w=7 aud=audit_ss.

* =============================================================================.
* PART B. Baseline covariates + wide file.
* CHECKPOINT B:  518 cases.
*   ANX5 M=.178 SD=.378 | ANX6 .161/.389 | ANX7 .206/.493
*   FAMSUP5 4.589/.451  | FAMSUP6 4.577/.481 | FAMSUP7 4.504/.544
*   AUDIT5 3.988/3.760  | AUDIT6 3.540/2.951 | AUDIT7 4.197/3.965
* =============================================================================.
GET FILE=T5 /KEEP=caseid demo_age LNHINC5.
RECODE demo_age (LOWEST THRU 0=SYSMIS).
RENAME VARIABLES (demo_age=AGE) (LNHINC5=LNHINC).
SORT CASES BY caseid.
SAVE OUTFILE='OUT/m1_t5cov.sav'.
GET FILE=MST /KEEP=caseid gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA.
SORT CASES BY caseid.
MATCH FILES /TABLE=* /FILE='OUT/m1_t5cov.sav' /BY caseid.
SAVE OUTFILE='OUT/m1_cov.sav'.
MATCH FILES /FILE='OUT/m1_T5.sav' /FILE='OUT/m1_T6.sav' /FILE='OUT/m1_T7.sav'
  /TABLE='OUT/m1_cov.sav' /BY caseid.
EXECUTE.
SAVE OUTFILE='OUT/m1_wide.sav'.
DESCRIPTIVES ANX5 ANX6 ANX7 FAMSUP5 FAMSUP6 FAMSUP7 AUDIT5 AUDIT6 AUDIT7
  /STATISTICS=MEAN STDDEV MIN MAX SKEWNESS.
FREQUENCIES gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA.

* =============================================================================.
* PART C. Long (person-period) file, wave dummies, z-scores, product term.
* CHECKPOINT C:  1,554 rows.  AUDIT n=1514 M=3.907 SD=3.592
*                ANX n=1513 M=.181 SD=.422 | FAMSUP n=1514 M=4.557 SD=.494
*                AGE M=35.232 SD=5.081 | LNHINC M=9.633 SD=2.892
* =============================================================================.
VARSTOCASES
  /MAKE AUDIT FROM AUDIT5 AUDIT6 AUDIT7
  /MAKE ANX FROM ANX5 ANX6 ANX7
  /MAKE FAMSUP FROM FAMSUP5 FAMSUP6 FAMSUP7
  /INDEX=WAVE(3)
  /KEEP=caseid gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA AGE LNHINC
  /NULL=KEEP.
COMPUTE WAVE = WAVE - 1.
VALUE LABELS WAVE 0 'T5' 1 'T6' 2 'T7'.
COMPUTE W6 = (WAVE = 1).
COMPUTE W7 = (WAVE = 2).
COMPUTE lnAUDIT = LN(AUDIT + 1).
EXECUTE.
DESCRIPTIVES AUDIT lnAUDIT ANX FAMSUP AGE LNHINC /SAVE.
COMPUTE ANXxFAM = zANX * zFAMSUP.
EXECUTE.
SORT CASES BY caseid WAVE.
SAVE OUTFILE='OUT/m1_long.sav'.

* =============================================================================.
* PART D. Three-wave GEE (PRIMARY).  Analyze > Generalized Linear Models > GEE.
* CHECKPOINT D1 (unadjusted): N = 1,511 obs / 518 subjects.
*    zANX*zFAMSUP (ANXxFAM)  B = -.115  SE = .039  Wald chi2 = 8.68  p = .003
* CHECKPOINT D2 (adjusted):   N = 1,418 obs / 485 subjects.
*    ANXxFAM  B = -.127  SE = .044  Wald chi2 = 8.16  p = .004
*    zANX .050 (p=.163) | zFAMSUP -.043 (p=.048) | male .756 | income .100 | educ -.153
* =============================================================================.
* D1 UNADJUSTED.
GENLIN zAUDIT WITH zANX zFAMSUP ANXxFAM W6 W7
  /MODEL zANX zFAMSUP ANXxFAM W6 W7 INTERCEPT=YES DISTRIBUTION=NORMAL LINK=IDENTITY
  /CRITERIA SCALE=MLE ANALYSISTYPE=3(WALD) CILEVEL=95 LIKELIHOOD=FULL
  /REPEATED SUBJECT=caseid WITHINSUBJECT=WAVE SORT=YES CORRTYPE=EXCHANGEABLE ADJUSTCORR=NO COVB=ROBUST
     MAXITERATIONS=100 PCONVERGE=1E-006(ABSOLUTE) UPDATECORR=1
  /MISSING CLASSMISSING=EXCLUDE
  /PRINT CPS DESCRIPTIVES MODELINFO FIT SUMMARY SOLUTION WORKINGCORR.
* D2 ADJUSTED.
GENLIN zAUDIT WITH zANX zFAMSUP ANXxFAM W6 W7 gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
  /MODEL zANX zFAMSUP ANXxFAM W6 W7 gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
     INTERCEPT=YES DISTRIBUTION=NORMAL LINK=IDENTITY
  /CRITERIA SCALE=MLE ANALYSISTYPE=3(WALD) CILEVEL=95 LIKELIHOOD=FULL
  /REPEATED SUBJECT=caseid WITHINSUBJECT=WAVE SORT=YES CORRTYPE=EXCHANGEABLE ADJUSTCORR=NO COVB=ROBUST
     MAXITERATIONS=100 PCONVERGE=1E-006(ABSOLUTE) UPDATECORR=1
  /MISSING CLASSMISSING=EXCLUDE
  /PRINT CPS DESCRIPTIVES MODELINFO FIT SUMMARY SOLUTION WORKINGCORR.

* =============================================================================.
* PART E. Simple slopes of anxiety at LOW (-1 SD) and HIGH (+1 SD) familism (re-centering).
*   In each run, the zANX coefficient IS the simple slope.
* CHECKPOINT E (adjusted): low familism  B = .177  SE = .075  p = .018
*                          high familism B = -.077 SE = .030  p = .010
*   (unadjusted: low .158 p=.015 ; high -.071 p=.018)
* =============================================================================.
COMPUTE FAMlo = zFAMSUP + 1.
COMPUTE FAMhi = zFAMSUP - 1.
COMPUTE ANXxFAMlo = zANX * FAMlo.
COMPUTE ANXxFAMhi = zANX * FAMhi.
EXECUTE.
GENLIN zAUDIT WITH zANX FAMlo ANXxFAMlo W6 W7 gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
  /MODEL zANX FAMlo ANXxFAMlo W6 W7 gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
     INTERCEPT=YES DISTRIBUTION=NORMAL LINK=IDENTITY
  /CRITERIA SCALE=MLE ANALYSISTYPE=3(WALD) CILEVEL=95
  /REPEATED SUBJECT=caseid WITHINSUBJECT=WAVE SORT=YES CORRTYPE=EXCHANGEABLE ADJUSTCORR=NO COVB=ROBUST
  /PRINT MODELINFO SOLUTION.
GENLIN zAUDIT WITH zANX FAMhi ANXxFAMhi W6 W7 gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
  /MODEL zANX FAMhi ANXxFAMhi W6 W7 gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
     INTERCEPT=YES DISTRIBUTION=NORMAL LINK=IDENTITY
  /CRITERIA SCALE=MLE ANALYSISTYPE=3(WALD) CILEVEL=95
  /REPEATED SUBJECT=caseid WITHINSUBJECT=WAVE SORT=YES CORRTYPE=EXCHANGEABLE ADJUSTCORR=NO COVB=ROBUST
  /PRINT MODELINFO SOLUTION.
* For UNADJUSTED simple slopes, delete the 9 covariates from WITH and /MODEL in both runs.

* =============================================================================.
* PART F. Per-wave replication with HC3 (PROCESS v4 Model 1).  Requires PROCESS macro loaded.
* CHECKPOINT F (adjusted interaction, "Int_1" row):
*    T5  B = -.178  SE(HC3) = .063  p = .004   (n = 484)
*    T6  B = -.094  SE(HC3) = .052  p = .070   (n = 472)
*    T7  B = -.230  SE(HC3) = .153  p = .132   (n = 462)
*  (unadjusted: T5 -.142 p=.048 | T6 -.166 p=.006 | T7 -.232 p=.113)
* =============================================================================.
SORT CASES BY WAVE caseid.
SPLIT FILE LAYERED BY WAVE.
DESCRIPTIVES VARIABLES=AUDIT (wzAUDIT) ANX (wzANX) FAMSUP (wzFAMSUP) AGE (wzAGE) LNHINC (wzLNHINC) /SAVE.
SPLIT FILE OFF.
SAVE OUTFILE='OUT/m1_long.sav'.
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
* Unadjusted per wave: remove the /cov= list.

* =============================================================================.
* PART G. Sensitivity analyses (skewed AUDIT).
* CHECKPOINT G1 negative binomial (adjusted): ANXxFAM  B = -.078  SE = .023  p = .001  Exp(B) = IRR = .925
*   (unadjusted: B = -.064, p = .004, IRR = .938)
* CHECKPOINT G2 ln(AUDIT+1) (adjusted): ANXxFAM p = .0003  (unadjusted p = .0004)
* =============================================================================.
SORT CASES BY caseid WAVE.
GENLIN AUDIT WITH zANX zFAMSUP ANXxFAM W6 W7 gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
  /MODEL zANX zFAMSUP ANXxFAM W6 W7 gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
     INTERCEPT=YES DISTRIBUTION=NEGBIN(1) LINK=LOG
  /CRITERIA ANALYSISTYPE=3(WALD) CILEVEL=95
  /REPEATED SUBJECT=caseid WITHINSUBJECT=WAVE SORT=YES CORRTYPE=EXCHANGEABLE ADJUSTCORR=NO COVB=ROBUST
  /PRINT MODELINFO FIT SUMMARY SOLUTION (EXPONENTIATED).
GENLIN zlnAUDIT WITH zANX zFAMSUP ANXxFAM W6 W7 gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
  /MODEL zANX zFAMSUP ANXxFAM W6 W7 gender maristat DOCSTAT EDUHI EMPLOY ORIG_SA ORIG_CA zAGE zLNHINC
     INTERCEPT=YES DISTRIBUTION=NORMAL LINK=IDENTITY
  /CRITERIA SCALE=MLE ANALYSISTYPE=3(WALD) CILEVEL=95
  /REPEATED SUBJECT=caseid WITHINSUBJECT=WAVE SORT=YES CORRTYPE=EXCHANGEABLE ADJUSTCORR=NO COVB=ROBUST
  /PRINT MODELINFO SOLUTION.
* For unadjusted sensitivity runs, delete the 9 covariates.
* END MODEL 1.
