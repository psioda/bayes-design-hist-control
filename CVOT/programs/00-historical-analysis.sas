%macro setup;
proc datasets library=work noprint kill; run; quit;
%global root;

ods _all_ close;
option nomprint nosymbolgen nomlogic;

%* for Windows;
%if &SYSSCP = WIN %then %do;
  %let root = C:\Users\psioda\Documents\Research Papers Temp\bayesDesignCVOT\sas_program_GitHub\CVOT;

  %* libraries;
  libname raw "&root.\data\raw_data";
  libname sp  "&root.\data\sampling_priors";
  libname sc  "&root.\data\simulation_controls";
  libname res "&root.\data\simulation_results";

  %* macro variables;
  %let sysparm = 1;

  %* macros;
  %include "&root.\..\macros\collapse.sas";
  %include "&root.\..\macros\fit_approx.sas";
  %include "&root.\..\macros\fit_mcmc.sas";
  %include "&root.\..\macros\simulateStudy.sas";
  %include "&root.\..\macros\subsetInterim.sas";

%end;
%else %do; %* for Linux;
 %let root = /nas/longleaf/home/psioda/stat-projects/bayesDesignCVOT/CVOT;

  %* libraries;
  libname raw "&root./data/raw_data";
  libname sp  "&root./data/sampling_priors";
  libname sc  "&root./data/simulation_controls";
  libname res "&root./data/simulation_results";

  %* macros;
  %include "&root./../macros/collapse.sas";
  %include "&root./../macros/fit_approx.sas";
  %include "&root./../macros/fit_mcmc.sas";
  %include "&root./../macros/simulateStudy.sas";
  %include "&root./../macros/subsetInterim.sas";

%end;

%mend setup;
%setup;

data hist_condensed; set raw.hist_condensed;        run;
data empty;          set raw.hist_condensed(obs=0); run;


%let dsNew        = work.hist_condensed;
%let dsHist       = work.empty;
%let a0           = 0;
%let dsOutParmEst = sp.histEstGenMod;
%let dsOutModFit  = sp.histModGenMod;
%let shareParms   = 1;
%fit_approx;

%let dsNew        = work.hist_condensed;
%let dsHist       = work.empty;
%let a0           = 0;
%let dsOutParmEst = sp.histEstMCMC;
%let noTreatment  = 1;

%let lambdaPrior  = gamma(0.0001,iscale=0.0001);
%let regPrior     = normal(0,sd=1000);
%let sampMethod   = normal;
%let seed         = 4823892;
%let mcmcOptions  = nmc=1000000 nbi=500 nthin=5 seed=&seed. postout=sp.histPostSamp;
%fit_mcmc;


/*ods html newfile=proc;*/
/**/
/*data x;*/
/* merge raw.hist_condensed(where=(riskstratum=1) rename=(numEvents = nu1 riskTime=r1))*/
/*       raw.hist_condensed(where=(riskstratum=2) rename=(numEvents = nu2 riskTime=r2));*/
/* by x1 x2 x3;*/
/* drop simStudy riskStratum interval x1;*/
/* format x2 5.1 x3 3. r1 r2 6.2 nu1 nu2 3.;*/
/**/
/**/
/*run;*/
/*quit; proc sort; by x3 x2; run;*/
/**/
/*data _null_;*/
/* set x;*/
/*  put x2 " & " x3 " & " nu1 " & " r1 " & " nu2 " & " r2 " //";*/
/*run;*/
/**/
/*proc print data = raw.hist_condensed; run; quit;*/
