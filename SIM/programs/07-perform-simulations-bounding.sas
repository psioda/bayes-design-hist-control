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


%let node_idx    = &sysparm;
%let histDS      = raw.hist_condensed;
%let breakDS     = raw.haz_breaks;

%let stratLevels  = 2;
%let fitHazComp   = 1 | 1;
%let fitHazBreaks = . | .;
%let rProb        = 0.5;
%let enrollDist   = rand('beta',1,2)*3;
%let censorDist   = 1000;

** possible control event fractions for borrowing (must always include 0.00 last);
%let borrowAmt = 0.00 | 0.50 ; 
%let w0        = 0.45;

%macro loop;

 %** create suffix for name of the results dataset;
 data _null_; call symput('dsLabel',put(%scan(&node_idx.,1),z5.)); run;

 %** extract relevant simulation settings;
 %** determine how many simulation loops are needed;
 data work.__controls__;
  set sc.simulation_controls_bound2;
  where node_idx in ( &node_idx.);

  row_idx = _n_;
  call symput('nSimSettings',strip(put(_n_,best.)));
 run;

 %** determine the number of events from the historical trial;
 proc sql noprint;
  select sum(numEvents) into :histNumEvents from &histDS.;
 quit;
 %let histNumEvents = %left(&histNumEvents);

 %** extract number of covariates;
 proc contents data = raw.hist_condensed out = nCov(where=(   substr(upcase(name),1,1)='X'  )) noprint; run; quit;
 data _null_;
  set nCov;
  call symput('nCov',strip(put(_n_,best.)));
 run;


%do loop     = 1 %to &nSimSettings.;
%put iteration &loop of &nSimSettings.;


 data start; x = time(); run;

  %** create macro variables for current iteration;
  data _null_;
   set work.__controls__;
   where row_idx = &loop.;

   
   call symput('numSimulations',strip(put(nPerLoop,best.)));
   call symput('numSubjects',   strip(put(numSubjects,best.)));
   call symput('nuTarget',      strip(put(nuTotal,best.)));
   call symput('extra_events',  strip(put(extra_events,best.)));

   call symput('seed',          strip(put(seed,30.)));
   call symput('parmSampleSeed',strip(put(parmSampleSeed,30.)));
   call symput('covSampleSeed', strip(put(covSampleSeed,30.)));

   call symput('hypothesis',    strip(put(hypothesis,best.)));
   call symput('perturbation',  strip(perturbation));
  run;

  data tempsp;
   set sp.parameter_perturbations_bound2;
   where perturbation="&perturbation";
  run;

  %** simulate final data for new CVOT;
  %let dsOut       = work.final_data;
  %let trueHazComp = &fitHazComp.;
  %let dsCovDist   = &histDS.;
  %let sampPrior   = tempsp;
  %let dsBreaks    = &breakDS.;
  %simulateStudy;




  %** loop over borrowing amounts and back calculate interim data for analysis;
  %let borrowAmtNum = %eval(%sysfunc(count(&borrowAmt,|))+1); 
  %do z = 1 %to &borrowAmtNum.;

   %let nuInterim = %sysfunc(ceil(&nuTarget.*(1-%scan(&borrowAmt,&z,|))));
   %let nuFill    = %sysevalf(&nuTarget.-&nuInterim.);
   %let nuBorrow  = %sysfunc(min(%sysevalf(&extra_events. + &nuFill.),&histNumEvents));
   %let a0        = %sysfunc(min(%sysevalf(&nuBorrow./&histNumEvents.),1));
   %if &nuFill.   = 0 %then %do; %let a0 = 0; %let nuBorrow = 0; %end;
   %let totalInfo = %sysevalf(&nuBorrow. + &nuInterim.);

   %*%put &=nuInterim &=a0 &=nuFill &=nuBorrow &=totalInfo;

    %** back calculate interim data;
    %let dsFinal   = work.final_data;
	%let dsInterim = work.int_data;
    %subsetInterim;

    %** collapse interim data;
	%let dsFull      = work.int_data;
	%let dsCollapsed = work.int_condensed;
	%let writeBreaks = 0;
	%let dsBreaks    = ;
    %collapse;

    %** analyze data with asymptotic approximation / shared parameters;
	%let dsNew        = work.int_condensed;
	%let dsHist       = &histDS.;
	%let shareParms   = 1;
	%let dsOutParmEst = est_shared;
    %let dsOutModFit  = modFit_shared;
    %fit_Approx;

    data est_shared;
     set est_shared;
     where Parameter = 'POSTPROB';

       length sig&z. 3.;
       sig&z. = (Estimate >=0.975);

     keep simStudy Estimate sig&z.;
     rename Estimate = pprob&z.;
    run;

    data modFit_Shared;
     set modFit_shared;
     where Criterion = 'Log Likelihood';

     keep simStudy value;
     rename value=value1_a0;
    run;

    %** analyze data with asymptotic approximation / unshared parameters;
	%let dsNew        = work.int_condensed;
	%let dsHist       = &histDS.;
	%let shareParms   = 0;
	%let dsOutParmEst = est_unshared;
    %let dsOutModFit  = modFit_unshared;
    %fit_Approx;


    data modFit_unshared;
     set modFit_unshared;
     where Criterion = 'Log Likelihood';
     keep simStudy value;
      rename value=value2_a0;
    run;

    data work.__Results&z.;
     length simStudy row_idx 8.;
     merge est_shared interimTime(rename=(interimTime=interimTime&z.)) 
           modFit_Shared modFit_unshared; 
     by simStudy;
       row_idx = &loop.;
       LRStat_a0&z. = value2_a0-value1_a0;
       drop value:;
    run;

  proc datasets library=work noprint;
   delete ModFit_Unshared: ModFit_Shared: est_shared est_unshared
          int_condensed int_data interimTime;
  quit;

  %end;

  data work.__Results;
   merge %do z=1 %to &borrowAmtNum.; work.__Results&z. %end;;
   by simStudy;

	 sig_ovr = sig1*(LRStat_a02 > &w0.) + sig2*(LRStat_a02 <= &w0.);
	 keep simStudy sig_ovr sig2 row_idx LRStat_a02;
	 rename LRStat_a02 = w;
  run;

  proc means data = work.__Results noprint nway;
   class row_idx;
   var sig_ovr sig2 w;
   output out = work.__Results(keep=row_idx weight rejRate rejRateInt w) n(sig_ovr)=weight mean=rejRate rejRateInt w;
  run;

  proc append data = work.__Results base=work.Results force; run; quit;

  proc datasets library=work noprint;
   delete __Results: final_data: tempsp mem;
  quit;

  ods output Members=work.mem;
  proc contents data=work._all_ details; run; quit;

  proc means data=work.mem noprint;
   var filesize;
   output out = work.mem sum=TFS;
  run;

  data _null_;
   set work.mem;
   call symput('TFS',strip(put(TFS,best.)));
  run;

data stop;
 set start;
  y = time();
  h = (y-x)/60/60;
  elapsed_min = round(h*60,0.001);

  node_idx = &node_idx.; 
  loop=&loop;
 
  TFS = put(&TFS.,SIZEKMG12.);
  put node_idx= loop= elapsed_min= TFS=;
run;


%end;
proc sort data = work.Results; by row_idx; run;

data Results2;
 merge work.__controls__(in=a drop=nSims seed covSampleSeed parmSampleSeed: inner_idx numSubjects nuTotal) 
       work.Results(in=b) ;
 by row_idx;
  if a and b;
  drop row_idx ;
run;

proc means data = Results2 noprint nway;
 class hypothesis perturbation;
 freq weight;
 var rejRate rejRateInt w;
 output out = Results2(keep=hypothesis perturbation nSims rejRate rejRateInt w)
    n(w)=nSims mean=rejRate rejRateInt w;
run;

data res.sim_bounding2_&dsLabel.;
 set Results2;
run;


%mend;
option nonotes;
**option nomprint nomlogic nosymbolgen;
%let debug=0;
%loop;
option notes;

