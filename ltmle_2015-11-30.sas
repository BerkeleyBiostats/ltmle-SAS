*****************************************************************************************;
* CAUSAL EFFECT CUMULATIVE EVENT PROBABILITY AT TIME K+1 								*;
*																						*;
* JORDAN BROOKS 																		*;
*																						*;
* General TMLE for longidutinal data structures from Laan & Gruber 2011					*;
* Adapted to reflect model constraints in right-censored survival data					*;
*																						*;
* TMLE with SL for INITIAL ESTIMATOR, \bar{Q}_{n,Y|\bar{L(k)}^{0}						*;
*****************************************************************************************;
options nofmtErr threads cpucount=6;
title 'TMLE FOR CUMULATIVE DEATH PROBABILITY AT TIME K+1: SIMULATION';

%LET home =C:\Users\Jordan\Desktop\CodeMark;
%include "&home\LTMLE\TMLE_RCS_DISSERTATION_042113.sas";
%include "&home\SL\SUPERLEARNER_JSS_040312.sas";
%include "&home\SL\library\LOGIT.sas";
%include "&home\SL\library\LOGIT_CTS01.sas";



****************************;
* SETUP SAS DATA LIBRARIES *;
****************************;
libname f0 "&home\data"; 

proc datasets lib=work kill; run; quit;
proc datasets lib=f0 kill; run; quit;

%LET ncounter = 1000;

*********************;
* SIMULATE DATA 	*;
*********************;
%LET myseed=715;
%LET n=1000;
%LET id = ID;
%LET t = k;
%LET A = A;
%LET K=10;
%LET eta = 2;

%LET g_lower_bound=0.01;
%LET g_upper_bound=0.99;

%LET X = W1 W2 W3 W4;
%LET nX = %sysfunc(COUNTW(&X, ' '));


*****************************************************;
* NOTE: PARAMETER IS E[Y^{a}=I(L^{a}(K+1))]			*;
* 	where a= \bar{a} = (A(0)=a, C(0)=...=C(K)=0)	*;
*****************************************************;

**************************************;
* BASELINE FIXED TREATMENT MECHANISM *;
**************************************;
%LET g0_A = 1/(1+exp(-1*&eta*(0.001*W1 + 0.01*W2 - 0.5*W3 + 0.5*W4 - 0.2*W3*W1))) ;
%LET g_Afitform = W1 W4; * <- FUNCTIONAL FORM FITTED BY GLM *;

**************************************;
* TIME-DEPENDENT COVARIATE MECHANISM *;
**************************************;
%LET Q0_dL = 1/(1+exp(-1*(-3 + 0.01*W1 -0.5*A + 0.01*&t))) ;

*********************************;
* HAZARD OF CENSORING MECHANISM *;
*********************************;
%LET g0_dC = 1/(1+exp(-1*(-7 + 0.001*W1 -0.002*W2 + 0.5*W3 - 3*W4 + 0.5*A + 4*L - 0.01*&t))) ;
%LET g_dCfitform = W1 W2 W3 W4 A L &t; * <- FUNCTIONAL FORM FITTED BY GLM *;

*********************************;
* HAZARD OF EVENT MECHANISM 	*;
*********************************;
%LET Q0_dN = 1/(1+exp(-1*(-5 + 0.01*W1 - 0.002*W2 + 2*W3 - 3*W4 - A + A*W3*W4 + 2*L + 0.01*&t))) ;



*********************************************************;
* COMPUTE TRUE PARAMETERS (INTERVENTION-SPECIFIC MEANS) *;
*********************************************************;
%MACRO TRUTH(setA);
data f0.truth;
 call streaminit(715); 
 length W1 W2 W3 W4 A L C N 3;
 do _i = 1 to 1000000;
  W1=.; W2=.; W3=.; W4=.; A=.; N=.; C=.; 
  ID = _i;
  k = 0;
  W1 = rand("Normal", 0, 1);
  W2 = rand("Normal", 0, 1);
  W3 = rand("Bernoulli", 0.5);
  W4 = rand("Bernoulli", 0.5);
  g0_A = &g0_A;
  A = &setA;
  L = rand("Bernoulli", &Q0_dL); 
  g0_dC = &g0_dC;
  C = 0;
  Q0_dN = &Q0_dN;	
  N = rand("Bernoulli", Q0_dN);
  if C=1 then N=0;
  output;
  if (N=0 and C=0) then do until (N=1 or C=1 or k=&K+1);
   k = k+1;
   Q0_dL = &Q0_dL;
   if (L=0) then L = rand("Bernoulli", Q0_dL); 
   g0_dC = &g0_dC;
   C = 0;
   Q0_dN = &Q0_dN;	
   N = rand("Bernoulli", Q0_dN);
   if C=1 then N=0;
   output;
  end;
 end; 
 drop _: ;
run;
data f0.psi0;
 set f0.truth;
 by ID k;
 if last.ID;
 Y_&setA = (k le &K);
run;
proc means noprint;
 var Y_&setA;
 output out=f0.truth_Y&setA;
run;
%MEND TRUTH;


title2 "COMPUTE TRUE INTERVENTION SPECIFIC MEANS: psi_0";
%TRUTH(setA=1);
data _null_;
 set f0.truth_Y1;
 if _STAT_ = "MEAN" then do; call symputx("Y_1", Y_1); end; 
run;

%TRUTH(setA=0);
data _null_;
 set f0.truth_Y0;
 if _STAT_ = "MEAN" then do; call symputx("Y_0", Y_0); end; 
run;

%LET diff = %sysevalf(&Y_1 - &Y_0);



%MACRO G;
%DO counter=1 %TO &ncounter;
dm output 'clear';
dm log 'clear';

*****************************;
* SIMULATE EMPIRICAL DATA 	*;
*****************************;
data f0.long;
 call streaminit(&counter); 
 length W1 W2 W3 W4 A C N 3;
 length g0_A g0_dC Q0_dL Q0_dN 8;
 format g0_A g0_dC Q0_dL Q0_dN 12.11;
 do _i = 1 to &n;
  W1=.; W2=.; W3=.; W4=.; A=.; N=.; C=.; 
  ID = _i;
  k = 0;
  W1 = rand("Normal", 0, 1);
  W2 = rand("Normal", 0, 1);
  W3 = rand("Bernoulli", 0.5);
  W4 = rand("Bernoulli", 0.5);
  g0_A = &g0_A;
  A = rand("Bernoulli", g0_A);
  Q0_dL = &Q0_dL;
  L = rand("Bernoulli", Q0_dL); 
  g0_dC = &g0_dC;
  C = rand("Bernoulli", g0_dC);
  Q0_dN = &Q0_dN;
  N = rand("Bernoulli", Q0_dN);
  if C=1 then N=0;
  output;
  if (N=0 and C=0) then do until (N=1 or C=1 or k=&K+1);
   k = k+1; 
   if (L=0) then Q0_dL = &Q0_dL; else Q0_dL = 1;
   L = rand("Bernoulli", Q0_dL); 
   g0_dC = &g0_dC; 
   C = rand("Bernoulli", g0_dC);
   Q0_dN = &Q0_dN;
   N = rand("Bernoulli", Q0_dN);
   if C=1 then N=0;
   output;
  end;
 end; 
 drop _: ;
run;

*********************************************;
* ESTIMATE TREATMENT ASSIGNMENT PROBABILITY *;
*********************************************;
data temp; set f0.long(where=(&t=0)); run;
%LOGIT(TRAIN=temp, 
	   Y=A, Y_TYPE=BIN, 
       X=&g_Afitform, 
	   ID=&ID, T=, WEIGHTS=, 
	   SEED=, 
	   WD=&home\fits);
data _null_;
  fname="fname";
  rc=filename(fname,"&home\fits\F_gn_A.sas");
  if fexist(fname) then rc=fdelete(fname);
  rc=rename("&home\fits\F_LOGIT.sas","&home\fits\F_gn_A.sas","file");
run;
data _null_;
 file "&home\fits\F_gn_A.sas" MOD;
 put "gn_A = p_LOGIT;";
 put "drop p_LOGIT ;";
run;

************************************;
* ESTIMATE THE HAZARD OF CENSORING *;
************************************;
%LOGIT(TRAIN=f0.long, 
	   Y=C, Y_TYPE=BIN, 
       X=&g_dCfitform, 
	   ID=&ID, T=&t, WEIGHTS=, 
	   SEED=, 
	   WD=&home\fits);
data _null_;
  fname="fname";
  rc=filename(fname,"&home\fits\F_gn_dC.sas");
  if fexist(fname) then rc=fdelete(fname);
  rc=rename("&home\fits\F_LOGIT.sas","&home\fits\F_gn_dC.sas","file");
run;
data _null_;
 file "&home\fits\F_gn_dC.sas" MOD;
 put "gn_dC = p_LOGIT;";
 put "drop p_LOGIT ;";
run;

********************************;
* ESTIMATE THE HAZARD OF EVENT *;
********************************;
%LOGIT(TRAIN=f0.long, 
	   Y=N, Y_TYPE=BIN, 
       X=W1 W2 W3 W4 A A*W3*W4 L &t, 
	   ID=&ID, T=&t, WEIGHTS=, 
	   SEED=, 
	   WD=&home\fits);
data _null_;
  fname="fname";
  rc=filename(fname,"&home\fits\F_Qn_dL1.sas");
  if fexist(fname) then rc=fdelete(fname);
  rc=rename("&home\fits\F_LOGIT.sas","&home\fits\F_Qn_dL1.sas","file");
run;
data _null_;
 file "&home\fits\F_Qn_dL1.sas" MOD;
 put "Qn_dL1 = p_LOGIT;";
 put "drop p_LOGIT ;";
run;

data f0.long;
 set f0.long;
 format gn_A gn_dC Qn_dL1 12.11;
 %include "&home\fits\F_gn_A.sas";
 %include "&home\fits\F_gn_dC.sas";
 %include "&home\fits\F_Qn_dL1.sas";
run;

*********;
* LTMLE *;
*********;
title2 'Targeted Maximum Likelihood Substitution Estimator';
%TMLE_RCS(DSN=f0.long, ID=ID, t=k, 
		  L1=N, L2=W1 W2 W3 W4 L, C=C, A=A, 
		  setA=0, K=10, 
		  gn_A=&home\fits\f_gn_A.sas, 
		  gn_dC=&home\fits\f_gn_dC.sas, 
		  g_LOWER_BOUND=0.01, g_UPPER_BOUND=0.99,
		  Qn_dL1=&home\fits\f_Qn_dL1.sas, 
		  Q_LOWER_BOUND=0.00000001, Q_UPPER_BOUND=0.99999999, 
		  SL_LIBRARY=MEAN OLS LOGIT_CTS01 NN2 TREE,
		  LIBRARY_DIR=&home\SL\library, 
		  EnterpriseMiner=T, V=5, SEED=715, FOLD=, VERBOSE=F,
	      WD=&home\fits);
data tmle_0;
 set cf (keep = &ID Qbar_star Dtmle
		 rename=(Qbar_star=Qbar_star_0
				 Dtmle=Dtmle_0));
run;

%TMLE_RCS(DSN=f0.long, ID=ID, t=k, 
		  L1=N, L2=W1 W2 W3 W4 L, C=C, A=A, 
		  setA=1, K=10, 
		  gn_A=&home\fits\f_gn_A.sas, 
		  gn_dC=&home\fits\f_gn_dC.sas, 
		  g_LOWER_BOUND=0.01, g_UPPER_BOUND=0.99,
		  Qn_dL1=&home\fits\f_Qn_dL1.sas, 
		  Q_LOWER_BOUND=0.00000001, Q_UPPER_BOUND=0.99999999, 
		  SL_LIBRARY=MEAN OLS LOGIT_CTS01 NN2 TREE,
		  LIBRARY_DIR=&home\SL\library, 
		  EnterpriseMiner=T, V=5, SEED=715, FOLD=, VERBOSE=F,
	      WD=&home\fits);
data tmle_1;
 set cf (keep = &ID Qbar_star Dtmle
		 rename=(Qbar_star=Qbar_star_1
				 Dtmle=Dtmle_1));
run;

data f0.tmle;
 merge tmle_1 (keep = &id Qbar_star_1 Dtmle_1) tmle_0 (keep = &id Qbar_star_0 Dtmle_0);
 by &id;
 diff = Qbar_star_1 - Qbar_star_0;
 Dtmle_diff = Dtmle_1 - Dtmle_0;
run;
proc means data=f0.tmle noprint;
 var Qbar_star_1 Qbar_star_0 diff;
 output out=diff;
run;
data _null_;
 set diff;
 if _STAT_="MEAN" then do; 
  call symputx("n", _FREQ_); 
  call symputx("tmle_1", Qbar_star_1); 
  call symputx("tmle_0", Qbar_star_0); 
  call symputx("tmle_diff", diff); 
 end;
run;
proc means data=f0.tmle noprint;
 var Dtmle_1 Dtmle_0 Dtmle_diff;
 output out=Dtmle_stats; 
run;
data _null_;
 set Dtmle_stats;
 if _STAT_="STD" then do; 
  call symputx("se_tmle_1", Dtmle_1/sqrt(&n)); 
  call symputx("se_tmle_0", Dtmle_0/sqrt(&n)); 
  call symputx("se_tmle_diff", Dtmle_diff/sqrt(&n)); 
 end;
 if _STAT_="MIN" then do; 
  call symputx("min_Dtmle_0", Dtmle_0); 
  call symputx("min_Dtmle_1", Dtmle_1); 
  call symputx("min_Dtmle_diff", Dtmle_diff); 
 end;
 if _STAT_="MAX" then do; 
  call symputx("max_Dtmle_0", Dtmle_0); 
  call symputx("max_Dtmle_1", Dtmle_1); 
  call symputx("max_Dtmle_diff", Dtmle_diff); 
 end;
run;
data tmle_est;
 tmle_1    = &tmle_1;
 tmle_0    = &tmle_0;
 tmle_diff = &tmle_diff;
 se_tmle_1 	  = &se_tmle_1;
 se_tmle_0 	  = &se_tmle_0;
 se_tmle_diff = &se_tmle_diff;
 min_Dtmle_1    = &min_Dtmle_1;
 min_Dtmle_0    = &min_Dtmle_0;
 min_Dtmle_diff = &min_Dtmle_diff;
 max_Dtmle_1    = &max_Dtmle_1;
 max_Dtmle_0    = &max_Dtmle_0;
 max_Dtmle_diff = &max_Dtmle_diff;
 tmle_1_lo = tmle_1 - 1.96*se_tmle_1;
 tmle_1_hi = tmle_1 + 1.96*se_tmle_1;
 tmle_0_lo = tmle_0 - 1.96*se_tmle_0;
 tmle_0_hi = tmle_0 + 1.96*se_tmle_0;
 tmle_diff_lo = tmle_diff - 1.96*se_tmle_diff;
 tmle_diff_hi = tmle_diff + 1.96*se_tmle_diff;
 tmle_cover_1 = (tmle_1_lo <= &Y_1 <= tmle_1_hi);
 tmle_cover_0 = (tmle_0_lo <= &Y_0 <= tmle_0_hi);
 tmle_cover_diff = (tmle_diff_lo <= &diff <= tmle_diff_hi);
run;
proc append base=f0.tmle_est data=tmle_est; run;

%END;
%MEND;
%G;


ods listing;
title2 'SIMULATION SETTINGS'; title3;
data _null_;
 file print;
 put "Sample size = &n";
 put "Treatment probability = &g0_A";
 put "Time-dependent covariate hazard = &Q0_dL";
 put "Censoring hazard = &g0_dC";
 put "Event hazard = &Q0_dN";
 put "K+1 = &K";
 put "Y_1  = &Y_1";
 put "Y_0  = &Y_0";
 put "diff = &diff";
 put "g lower_bound = &g_lower_bound";
 put "g upper_bound = &g_upper_bound";
run;


title2 'ESTIMATES'; title3;
proc means data=f0.tmle_est; var tmle_1 tmle_0 tmle_diff; run;

title2 '95% CI COVERAGE'; title3;
proc means data=f0.tmle_est mean; var tmle_cover_1 tmle_cover_0 tmle_cover_diff; run;

title2 'MEAN SQUARED ERROR'; title3;
data tmle_SQERR;
 set f0.tmle_est;
 tmle_1_SQERR = (tmle_1 - &Y_1)**2;
 tmle_0_SQERR = (tmle_0 - &Y_0)**2;
 tmle_diff_SQERR = (tmle_diff - &diff)**2; 
run;
proc means data=tmle_SQERR; var tmle_1_SQERR tmle_0_SQERR tmle_diff_SQERR; run;

title2 'INFLUENCE CURVES'; title3;
proc means data=f0.tmle_est; var min_Dtmle_0 max_Dtmle_0; run;
proc means data=f0.tmle_est; var min_Dtmle_1 max_Dtmle_1; run;
proc means data=f0.tmle_est; var min_Dtmle_diff max_Dtmle_diff; run;

