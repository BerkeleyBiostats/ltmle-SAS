********************************************************************************;
* %TMLE_RCS                                                                    *;
*                                                                              *;
* SAS MACRO FOR FOR TARGETED MARXIMUM LIKELIHOOD ESTIMATION OF THE             *;
* FIXED INTERVENTION-SPECIFIC MARGINAL CUMULATIVE EVENT PROBABILITY IN         *;
* SURVIVAL DATA WITH TIME-DEPENDENT COVARIATES AND WITH INFORMATIVE RIGHT-     *;
* CENSORING.                                                                   *;
*                                                                              *;
* REQUIRES %SUPERLEARNER SAS MACRO FOR THE INITIAL Q-FACTOR ESTIMATORS         *;
*                                                                              *;
* JORDAN BROOKS                                                                *;
*                                                                              *;
********************************************************************************;
* TESTED ON SAS 9.2 WINDOWS 64-BIT                                             *;
********************************************************************************;
* INPUTS TO THE %TMLE_RCS MACRO:                                               *;
*                                                                              *;
* DSN                = Input SAS dataset containing observations in long       *;
*                      format, i.e., one row per subject per time-point.       *;
*                                                                              *;
* ID                 = Name of the variable that uniquely identifies           *;
*                      independent observations/subjects.                      *;
*                                                                              *;
* t                  = Name of the time-stamp variable. The first row for each *;
*                      subject has t=0.                                        *;
*                                                                              *;
* L1                 = The event counting process. This is really L1(t+1), a   *;
*                      random variable that takes value 0 until the time-point *;
*                      just before the event occurs. The first row for each    *;
*                      subject is then L1(0+1), the indicator that the event   *;
*                      process jumps at time t=1. Once the event process jumps *;
*                      It stays at 1 for all remaining time-points.            *;
*                      Must be NUMERIC coded 0/1.                              *;
*                      Missing values are not allowed.                         *;
*                                                                              *;
* L2                 = The time-dependent covariate process history. This is   *;
*                      is really L2(t), a vector of random variables that      *;
*                      represent the covariate history up to and including     *;
*                      time-point t. All variables here must be NUMERIC.       *;
*                      Missing values are not allowed.                         *;
*                                                                              *;
* A                  = The fixed treatment. Must be a NUMERIC coded 0/1.       *;
*                      Missing values are not allowed.                         *;
*                                                                              *;
* C                  = The censoring process. Really this is C(t), a random    *;
*                      variable that takes value 0 at every time-point before  *;
*                      censoring occurs. When censoring occurs, the process    *;
*                      jumps to C(t)=1, and stays at 1 for all reamining time- *;
*                      points. Must be NUMERIC coded 0/1.                      *;
*                      Missing values are not allowed.                         *;
*                                                                              *;
* K                  = The pen-ultimate time-point of interest. That is, the   *;
*                      parameter of interest is the marginal intervention-     *;
*                      specific cumulative event probability by time t=K+1.    *;
*                                                                              *;
* setA               = The fixed treatment intervention of interest.           *;
*                                                                              *;
* gn_A               = A .sas program with DATA STEP code to compute the fixed *;
*                      treatment probabilities. This is a function of L2(0).   *;
*                                                                              *;
* gn_dC              = A .sas program file with DATA STEP code to compute the  *;
*                      censoring intensity esimtates.                          *;
*                      This is a function of L2, A, and t.                     *;
*                                                                              *;
* g_LOWER_BOUND      = Lower bound for g estimates. The default is 0.001       *;
*                                                                              *;
* g_LOWER_BOUND      = Lower bound for g estimates. The default is 0.999       *;
*                                                                              *;
* Qn_dL1             = A .sas program file containing the estimator mapping    *;
*                      for the event intensity/hazard. This is a function of   *;
*                      as a function L2, A, and t. This is used as the initial *;
*                      estimator of the ultimate Q-factor,                     *;
*                      Qbar_(K) = E[L1(K+1)|setA, L2(t=K)], in the TMLE.       *;
*                                                                              *;
* Q_LOWER_BOUND      = Lower bound for Qbar_(k) = E[L1(K+1)|setA, L2(t)]       *;
*                      The default is 0.00000001                               *;
*                                                                              *;
* Q_UPPER_BOUND      = Upper bound for Qbar_(k) = E[L1(K+1)|setA, L2(t)]       *;
*                      The default is 0.99999999                               *;
*                                                                              *;
* SL_LIBRARY         = Names of SAS macros for candidate estimators in the     *;
*                      Super Learner library. The macros must be saved under   *;
*                      the exact name of candidate estimator as a .sas program *;
*                      in the filepath given in LIBRARY_DIR.                   *;
*                                                                              *;
* LIBRARY_DIR        = Filepath for directory containing the candidate         *;
*                      estimator macros, which are saved as .sas programs.     *;
*                                                                              *;
* EnterpriseMiner    = Enter T/F to indicate use of SAS/EnterpriseMiner. If    *;
*                      set to T, a SAS datamining database catalog is          *;
*                      constructed.                                            *;
*                                                                              *;
* V                  = The number of folds for V-fold crossvalidation.         *;
*                                                                              *;
* SEED               = The SEED for the random number generator.               *;
*                                                                              *;
* FOLD (optional)    = Variable containing the fold number for cross           *;
*                      validation in the Super Learner.                        *;
*                                                                              *;
* VERBOSE (optional) = A T/F that indicates whether all output normally        *;
*                      produced by SAS data and procedure steps should be      *;
*                      sent to the listing device. The default is F.           *;
*                                                                              *;
* WD                 = The working directory. Initial estimator and TMLE       *;
*                      updated fits will be saved to this directory as .sas    *;
*                      programs. This directory will be assigned a SAS libname *;
*                      "SL" during each fitting process.                       *;
********************************************************************************;

%MACRO TMLE_RCS(DSN, ID, t, L1, L2, A, C, K, setA,
                gn_A, gn_dC, g_LOWER_BOUND, g_UPPER_BOUND,
                Qn_dL1, Q_LOWER_BOUND, Q_UPPER_BOUND,
                SL_LIBRARY, LIBRARY_DIR, EnterpriseMiner, V, SEED, FOLD, VERBOSE,
                WD);

ods listing close;

***************************************************************************;
* AUGMENT DATASET TO INCLUDE ONE ROW PER SUBJECT PER TIME-POINT UNTIL t=K *;
***************************************************************************;
data tmle;
 set &DSN (where=(&t <= &K) keep = &ID &t &L1 &L2 &A &C);
 by &ID &t;
 dead=0;
 output;
 if last.&ID and &t lt &K then do until (&t=&K);
  &t = &t+1;
  dead = &L1;
  output;
 end;
 drop _:;
run;

******************************************************************************;
* COMPUTE g (treatment and censoring) AND Qbar_K (event intensity) ESTIMATES *;
******************************************************************************;
data tmle;
 set tmle;
 by &ID &t;
 retain gn_A;
 if first.&ID then do;
  %include "&gn_A";
 end;
 keep &ID &t &L1 &L2 &A &C gn_A dead;
run;
data tmle;
 set tmle;
 by &ID &t;
 %include "&gn_dC";
 keep &ID &t &L1 &L2 &A &C gn_A gn_dC dead;
run;
data tmle;
 set tmle;
 by &ID &t;
 %include "&Qn_dL1";
 retain g;
 if first.&ID then do;
  if &A=1 then g=gn_A;
  if &A=0 then g=(1-gn_A);
 end;
 g=g*(1-gn_dC);
 g = max(g, &g_LOWER_BOUND);
 g = min(g, &g_UPPER_BOUND);
 H = (&A=&setA and &C=0)/g;
 Qn_dL1 = Qn_dL1;
 keep &ID &t &L1 &L2 &A &C g H Qn_dL1 dead;
run;

*******************************************;
* SORT DESCENDING BY t FOR TMLE ALGORITHM *;
*******************************************;
proc sort data=tmle;
 by &ID descending &t;
run;

*************************************************************************;
* DEFINE ACTIVE SET OF UNITS THAT HAD OBSERVED INTERVENTION OF INTEREST *;
*************************************************************************;
%LET step = &K;
%LET active = ((&t=&K) and (&A=&setA) and (&C=0)) ;
title2 "TMLE ALGORITHM STEP &K";

************************************************************;
* COMPUTE ESTIMATE OF Qbar_(K) = E[Y=L(K+1)|setA, L2(t=K)] *;
* IF EVENT HAS OCCURRED PRIOR CURRENT DAY THEN Qbar_(K)=1  *;
************************************************************;
data tmle;
 set tmle;
 by &ID descending &t;
 if ((&t=&K) and (&A=&setA)) then do;
  Qbar=(Qn_dL1);
  if (dead) then Qbar=1;
 end;
run;

*******************;
* FIT epsilon_(K) *;
*******************;
proc nlmixed data=tmle(where=((&active) and (dead=0)));
 parms eps=0;
 Qbar_star = 1/(1+exp(-1*(log(Qbar/(1-Qbar)) + eps)));
 logL = 2*H*(&L1*log(Qbar_star)+(1-&L1)*log(1-Qbar_star));
 model &L1 ~ general(logL);
 ods output ParameterEstimates=_eps;
run;
data _null_;
 set _eps;
 call symputx("eps_&K", Estimate);
run;
proc datasets lib=work; delete _: ; quit;
data eps;
 t = %EVAL(&K+1);
 &A = &setA;
 eps = &&eps_&K;
 output;
run;

*******************;
* UPDATE Qbar_(K) *;
*******************;
data tmle;
 set tmle;
 by &ID descending &t;
 if ((&t=&K) and (&A=&setA)) then do;
  if (dead=1) then Qbar_Star=1;
  else do;
   linpart = log(Qbar/(1-Qbar)) + &&eps_&K;
   if linpart < -700 then Qbar_star = 0;
   else if linpart > 700 then Qbar_star = 1;
   else Qbar_star = 1/(1+exp(-1*(linpart)));
  end;
 end;
run;

*******************************************;
* COMPUTE INFLUENCE CURVE COMPONENT D_(K) *;
*******************************************;
data tmle;
 set tmle;
 if &t=&K then do;
  D=0;
  if (H ne 0) then do;
   D = H*(&L1-Qbar_star);
  end;
 end;
run;

***************;
* NOW ITERATE *;
***************;
%DO i = 1 %to %sysevalf(&K);

*************************************************************************;
* DEFINE ACTIVE SET OF UNITS THAT HAD OBSERVED INTERVENTION OF INTEREST *;
*************************************************************************;
%LET step = %EVAL(&K-&i);
%LET active = ((&t=&K-&i) and (&A=&setA) and (&C=0)) ;
title2 "TMLE ALGORITHM STEP &step" ;

******************************************************;
* DEFINE PREVIOUS CONDITIONAL EXPECTATION AS OUTCOME *;
******************************************************;
data tmle;
 set tmle;
 by &ID descending &t;
 lagQbar_star=lag(Qbar_star);
run;

*************************************************************************************;
* SUPER LEARNER FOR INITIAL ESTIMATOR Qbar_(K-&i) = E[Qbar_(K-&i+1)|setA, L2(K-&i)] *;
*************************************************************************************;
data temp; set tmle(where=((&active) and (dead=0))); run;
%SUPERLEARNER(TRAIN=temp, TEST=,
              Y=lagQbar_star, Y_TYPE=CTS, X=&L2, ID=&ID, T=, WEIGHTS=,
              SL_LIBRARY=&SL_LIBRARY,
              LIBRARY_DIR=&LIBRARY_DIR,
              EnterpriseMiner=&EnterpriseMiner,
              LOWER_BOUND=&Q_LOWER_BOUND, UPPER_BOUND=&Q_UPPER_BOUND,
              LOSS=LOG, V=&V, SEED=&SEED, FOLD=&FOLD,
              WD=&WD, VERBOSE=F);
proc datasets kill lib=sl; run; quit;
data _null_;
  fname="fname";
  rc=filename(fname,"&WD\F_SuperLearner_&setA._&step..sas");
  if fexist(fname) then rc=fdelete(fname);
  rc=rename("&WD\F_SuperLearner.sas","&WD\F_SuperLearner_&setA._&step..sas","file");
run;

**************************************************************************;
* COMPUTE INITIAL ESTIMATE Qbar_(K-&i) = E[Qbar_(K-&i+1)|setA, L2(K-&i)] *;
**************************************************************************;
data tmle;
 set tmle;
 by &ID descending &t;
 %include "&WD\F_SuperLearner_&setA._&step..sas";
 if ((&t=&K-&i) and (&A=&setA)) then do;
  if (dead) then Qbar=1;
  else Qbar = p_SL;
 end;
run;

**********************;
* FIT epsilon_(K-&i) *;
**********************;
proc nlmixed data=tmle(where=((&active) and (dead=0)));
 parms eps=0;
 Qbar_star = 1/(1+exp(-1*(log(Qbar/(1-Qbar)) + eps)));
 logL = 2*H*(lagQbar_star*log(Qbar_star)+(1-lagQbar_star)*log(1-Qbar_star));
 model lagQbar_star ~ general(logL);
 ods output ParameterEstimates=_eps;
run;
data _null_;
 set _eps;
 call symputx("eps_&step", Estimate);
run;
data _eps;
 t = &step;
 &A = &setA;
 eps = &&eps_&step;
 output;
run;
proc append base=eps data=_eps; run;
proc datasets lib=work; delete _:; quit;

**********************;
* UPDATE Qbar_(K-&i) *;
**********************;
data tmle;
 set tmle;
 if ((&t=&K-&i) and (&A=&setA)) then do;
  if (dead) then Qbar_Star=1;
  else do;
   linpart = log(Qbar/(1-Qbar)) + &&eps_&step;
   if linpart < -700 then Qbar_star = 0;
   else if linpart > 700 then Qbar_star = 1;
   else Qbar_star = 1/(1+exp(-1*(linpart)));
  end;
 end;
run;

**********************************************;
* COMPUTE INFLUENCE CURVE COMPONENT D_(K-&i) *;
**********************************************;
data tmle;
 set tmle;
 lagD = lag(D);
 if &t=&K-&i then do;
  D=lagD;
  if (H ne 0) and (lagQbar_star ne .) then do;
   D = D + H*(lagQbar_star-Qbar_star);
  end;
 end;
run;

%END;

****************************************************************************;
* SET A and C FOR ALL SUBJECTS AT t=0 AND COMPUTE COUNTERFACTUAL ESTIMATES *;
****************************************************************************;
data cf;
 set tmle (keep = &ID &t &L1 &L2 D where=(&t=0));
 &A=&setA;
 &C=0;
 %include "&gn_A";
 keep &ID &t &L1 &L2 &A &C gn_A D;
run;
data cf;
 set cf;
 %include "&gn_dC";
 keep &ID &t &L1 &L2 &A &C gn_A gn_dC D;
run;
data cf;
 set cf;
 %include "&WD\f_SuperLearner_&setA._0.sas";
 %IF &setA = 1 %THEN %DO; H = 1/ ((gn_A)*(1-gn_dC)); %END;
 %IF &setA = 0 %THEN %DO; H = 1/ ((1-gn_A)*(1-gn_dC)); %END;
 Qbar = p_SL;
 linpart = log(Qbar/(1-Qbar)) + &eps_0;
 if linpart < -700 then Qbar_star = 0;
 else if linpart > 700 then Qbar_star = 1;
 else Qbar_star = 1/(1+exp(-1*(linpart)));
 keep &ID &t &L1 &L2 &A &C gn_A gn_dC Qbar H Qbar_star D; 
run;

*********************************************************************;
* COMPUTE TMLE AS THE EMPIRCAL MEAN OF THE COUNTERFACTUAL ESTIMATES *;
*********************************************************************;
proc means data=cf noprint;
 var Qbar_star;
 output out=cf_stats;
run;
data _null_;
 set cf_stats;
 if _STAT_="MEAN" then do;
  call symputx("n", _FREQ_);
  call symputx("tmle", Qbar_star);
 end;
run;

*******************************************************;
* COMPUTE STANDARD ERROR BASED ON THE INFLUENCE CURVE *;
*******************************************************;
data cf;
 set cf (keep = &ID &t &L2 &A D Qbar H Qbar_star);
 by &ID;
 Dtmle = D + Qbar_star - &tmle;
run;
proc means data=cf noprint;
 var Dtmle;
 output out=Dtmle_stats;
run;
data _null_;
 set Dtmle_stats;
 if _STAT_="STD" then do;
  call symputx("se_tmle", Dtmle/sqrt(&n));
 end;
run;

*******************************************;
* PRINT TMLE AND STANDARD ERROR ESTIMATES *;
*******************************************;
data _null_;
 file print;
 put "Marginal cumulative probability of &L1=1 by time &K";
 put " under fixed treatment: &A = &setA";
 put "Sample size = &n";
 put "g-factors bounded = [&g_LOWER_BOUND, &g_UPPER_BOUND]";
 put "TMLE (SE)  = &tmle (&se_tmle)";
run;

%MEND TMLE_RCS;