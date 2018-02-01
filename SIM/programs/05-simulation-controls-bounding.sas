%macro setup;
proc datasets library=work noprint kill; run; quit;
%global root;

ods _all_ close;
option nomprint nosymbolgen nomlogic;

%if &SYSSCP = WIN %then %do;
  %let root = C:\Users\psioda\Documents\Research Papers Temp\bayesDesignCVOT\sas_program_GitHub\SIM;

  %* libraries;
  libname raw "&root.\data\raw_data";
  libname sp  "&root.\data\sampling_priors";
  libname sc  "&root.\data\simulation_controls";
  libname res "&root.\data\simulation_results";

  %* macro variables;
  %let sysparm = 2;

  %* macros;
  %include "&root.\..\macros\collapse.sas";
  %include "&root.\..\macros\fit_approx.sas";
  %include "&root.\..\macros\fit_mcmc.sas";
  %include "&root.\..\macros\simulateStudy.sas";
  %include "&root.\..\macros\subsetInterim.sas";

%end;
%else %do;
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


%let number_of_nodes = 800;
%let nPerLoop        = 250;
%let nLoops          = 600;

%let numSubjects     = 5000;
%let borrowAmt       = 0.50;
%let nuTarget        = 612;
%let histNumEvents   = 605;

%let a0        = %sysevalf(&nuTarget.*&borrowAmt./&histNumEvents.    ); %put &=a0;
%let nuInterim = %sysfunc(ceil(&nuTarget.-&nuTarget.*&borrowAmt.  )); %put &=nuInterim;
%let w0        = 0.45;

** compute posterior means;
proc means data = sp.histPostSamp noprint nway;
 var lambda1 lambda2 beta1 beta2;
 output out = parameter_perturbations(drop=_:) mean=haz11 haz21 b1 b2;
run;

data parameter_perturbations;
 set parameter_perturbations;
  gamma = 0;
run;

data parameter_perturbations_bounding(keep=perturbation h11 h21 beta1 beta2 gamma );
 length perturbation $50 ;
 set parameter_perturbations;

  array h(*) haz11 haz21;
  array b(*) b1 b2;


  do h11_shift = 0.65 to 1.35 by 0.050;
  do h21_shift = 0.65 to 1.35 by 0.050;
  do b1_shift  = 0.65 to 1.35 by 0.050;
  do b2_shift  = 0.65 to 1.35 by 0.050;

  perturbation = catx('-',put(h11_shift,5.2),put(h21_shift,5.2),put(b1_shift,5.2),put(b2_shift,5.2));

  h11 = haz11*h11_shift;
  h21 = haz21*h21_shift;
  beta1 = b1 + log(b1_shift);
  beta2 = b2 + log(b2_shift);

   output;
  end;
  end;
  end;
  end;


  drop h11_: h21_: b1_: b2_: ;
run;

proc IML;

 use raw.hist_condensed;             
   read all var {x1 x2 x3} into x; 
   read all var {riskStratum} into strata;
   read all var {numEvents} into event;
   read all var {riskTime} into risk;
 close raw.hist_condensed;

 use parameter_perturbations_bounding; 
     read all var {h11 h21} into haz;
     read all var {gamma beta1 beta2} into beta;
	 read all var {perturbation} into perturbation;
 close parameter_perturbations_bounding;
 

 use parameter_perturbations;
     read all var {haz11 haz21} into haz_opt;
     read all var {gamma b1 b2} into beta_opt;
 close parameter_perturbations;


 rParms = nrow(beta);
 cParms = ncol(beta);

  a0         = &a0.;
  b          = t(beta_opt);
  xb         = x*b;
  h          = haz_opt[strata];
 logLike_opt = (a0 * (  event # log(h) + event # xb + event # log(risk) - risk#h#exp(xb)   ))[+];

 w = J(rParms,2,0);
 do j = 1 to rParms;
  b  = t(beta[j,]);
  xb = x*b;
  h  = t(haz[j,]);
  h  = h[strata];

  logLike   = (a0 * (  event # log(h) + event # xb + event # log(risk) - risk#h#exp(xb)   ))[+];
  w[j,1]    = logLike_opt-logLike;
  w[j,2]    = logLike_opt;

 end;

 create w from w[c={"W" "logLike_opt"}];
  append from w;
 close w;
quit;

data sp.parameter_perturbations_bound;
 merge parameter_perturbations_bounding w;

 if w <= 8*&w0.;
 rename w = w_hist;
run;


data parameter_perturbations_bounding;
 set sp.parameter_perturbations_bound;
 hypothesis = 0; output;
 hypothesis = 1; output;
run; proc sort; by hypothesis; run;

data parameter_perturbations_bounding;
 length  node_idx 4.;
 set parameter_perturbations_bounding;

   node_idx+1;
   if node_idx > &number_of_nodes. then node_idx=1;
run;

data sc.simulation_controls_bound;
  length  node_idx inner_idx nPerLoop nLoops nSims numSubjects nuTotal 5.;
 set parameter_perturbations_bounding(keep=node_idx hypothesis perturbation);

 nPerLoop    = &nPerLoop.;
 nLoops      = &nLoops.;
 nSims       = nPerLoop*nLoops;
 numSubjects = &numSubjects.;
 nuTotal     = &nuTarget.;

 extra_events = 0;

 do inner_idx = 1 to nLoops;
    seed           = round(1 + rand('uniform')*2**30);
    covSampleSeed  = round(1 + rand('uniform')*2**30);
    parmSampleSeed = round(1 + rand('uniform')*2**30);
    output;
 end;
 drop nLoops;
run;

data parameter_perturbations_bounding(keep=perturbation h11 h21 beta1 beta2 gamma );
 length perturbation $50 ;
 set parameter_perturbations;

  array h(*) haz11 haz21;
  array b(*) b1 b2;


  do h11_shift = 0.65 to 1.35 by 0.01;
  h21_shift = h11_shift;
  do b1_shift  = 1.00 to 1.00 by 0.01;
  do b2_shift  = 1.00 to 1.00 by 0.01;

  perturbation = catx('-',put(h11_shift,5.2),put(h21_shift,5.2),put(b1_shift,5.2),put(b2_shift,5.2));

  h11 = haz11*h11_shift;
  h21 = haz21*h21_shift;
  beta1 = b1 + log(b1_shift);
  beta2 = b2 + log(b2_shift);

   output;
  end;
  end;
  end;


  drop h11_: h21_: b1_: b2_: ;
run;



proc IML;

 use raw.hist_condensed;             
   read all var {x1 x2 x3} into x; 
   read all var {riskStratum} into strata;
   read all var {numEvents} into event;
   read all var {riskTime} into risk;
 close raw.hist_condensed;

 use parameter_perturbations_bounding; 
     read all var {h11 h21} into haz;
     read all var {gamma beta1 beta2} into beta;
	 read all var {perturbation} into perturbation;
 close parameter_perturbations_bounding;
 

 use parameter_perturbations;
     read all var {haz11 haz21} into haz_opt;
     read all var {gamma b1 b2} into beta_opt;
 close parameter_perturbations;


 rParms = nrow(beta);
 cParms = ncol(beta);

  a0         = &a0.;
  b          = t(beta_opt);
  xb         = x*b;
  h          = haz_opt[strata];
 logLike_opt = (a0 * (  event # log(h) + event # xb + event # log(risk) - risk#h#exp(xb)   ))[+];

 w = J(rParms,2,0);
 do j = 1 to rParms;
  b  = t(beta[j,]);
  xb = x*b;
  h  = t(haz[j,]);
  h  = h[strata];

  logLike   = (a0 * (  event # log(h) + event # xb + event # log(risk) - risk#h#exp(xb)   ))[+];
  w[j,1]    = logLike_opt-logLike;
  w[j,2]    = logLike_opt;

 end;

 create w from w[c={"W" "logLike_opt"}];
  append from w;
 close w;
quit;

data sp.parameter_perturbations_bound2;
 merge parameter_perturbations_bounding w;
 if w <= 4*&w0.;
 rename w = w_hist;
run;


data parameter_perturbations_bound2;
 set sp.parameter_perturbations_bound2;
 hypothesis = 0; output;
 hypothesis = 1; output;
run; proc sort; by hypothesis; run;

data parameter_perturbations_bound2;
 length  node_idx 4.;
 set parameter_perturbations_bound2;

   node_idx+1;
   if node_idx > &number_of_nodes. then node_idx=1;
run;

data sc.simulation_controls_bound2;
  length  node_idx inner_idx nPerLoop nLoops nSims numSubjects nuTotal 5.;
 set parameter_perturbations_bound2(keep=node_idx hypothesis perturbation);

 nPerLoop    = &nPerLoop.;
 nLoops      = &nLoops.;
 nSims       = nPerLoop*nLoops;
 numSubjects = &numSubjects.;
 nuTotal     = &nuTarget.;

 extra_events = 0;

 do inner_idx = 1 to nLoops;
    seed           = round(1 + rand('uniform')*2**30);
    covSampleSeed  = round(1 + rand('uniform')*2**30);
    parmSampleSeed = round(1 + rand('uniform')*2**30);
    output;
 end;
 drop nLoops;
run;

