%macro setup;
proc datasets library=work noprint kill; run; quit;
%global root;

ods _all_ close;
option nomprint nosymbolgen nomlogic;

%* for Windows;
%if &SYSSCP = WIN %then %do;
  %let root = C:\Users\psioda\Documents\Research Papers Temp\bayesDesignCVOT\sas_program_GitHub\SIM;

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
 %let root = /nas/longleaf/home/psioda/stat-projects/bayesDesignCVOT/SIM;

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

** macro variables for historical study;
%let N0         = 6000;
%let sProb      = 0.5; ** two strata;
%let Haz0       = 0.020 0.035; ** hazard in each stratum;
%let beta0      = 0.4 0.4;
%let seed0      = 75959;
%let nuTarg     = 605;
%let censorDist = rand('exponential')/0.05;

data historical;
 call streaminit(&seed0.);
 array h[2] (&Haz0.);
 array b[2] (&beta0.);
 array x[3] x1 x2 x3;

 z        = 0;
 simStudy = 1;

 do subject = 1 to &N0.;
   riskStratum = rand('bernoulli',&sProb.)+1;
   x[1]        = 0;
   x[2]        = round(rand('uniform')*1.2-0.6,0.2);
   x[3]        = rand('bernoulli',0.5);
   y           = round(rand('exponential')/(h[riskStratum]*exp(x[2]*b[1] + x[3]*b[2])),1e-4);
   c           = round(&censorDist.,1e-4);

   e = 1;
   if y > c then do;
    y = c;
	e = 0;
   end;

   output;
 end;

run; proc sort; by descending e y; run;

%let targetTime = 9999999;
data _null_;
 set historical;
 if e = 1 and _n_ = &nuTarg. then do;
  call symput('targetTime',strip(put(y,best.)));
 end;
run; %put &=targetTime;

data historical;
 set historical;
  if y > &targetTime then do;
    y =  &targetTime.;
	e = 0;
  end;

  censor = 1-e;
  event  = e;
  obsTime = y;
  
  keep SimStudy subject obsTime censor event riskStratum z x:;
run;



data raw.hist;
 retain SimStudy subject obsTime censor event riskStratum z;
 set historical;
run;

%let fitHazComp   = 1 | 1;
%let fitHazBreaks = . | .;
%let stratlevels  = 2; 
%let dsFull       = raw.hist;
%let dsCollapsed  = raw.hist_condensed;
%let writeBreaks  = 1;
%let dsBreaks     = raw.haz_breaks;
%collapse;


data hist_condensed; set raw.hist_condensed; run;
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
%let sampMethod   = slice;
%let seed         = 4823892;
%let mcmcOptions  = nmc=1000000 nbi=500 nthin=5 seed=&seed. postout=sp.histPostSamp;
%fit_mcmc;
