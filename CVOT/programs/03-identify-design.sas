%macro setup;
proc datasets library=work noprint kill; run; quit;
%global root;

ods _all_ close;
option nomprint nosymbolgen nomlogic;

%if &SYSSCP = WIN %then %do;
  %let root = C:\Users\psioda\Documents\Research Papers Temp\bayesDesignCVOT\sas_program_GitHub\CVOT;

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

libname readonly (res) access=read;

proc contents data = readonly._all_ out = contents noprint; run; quit;
proc sort data = contents(keep=memname) nodupkey; by memname; run;

 ** create suffix for name of the dataset;
 data _null_;
  set contents;
  where upcase(memname) =: 'SIM_RESULTS';
  call symput('ds'||strip(_n_),strip(memname));
  call symput('numDS',strip(put(_n_,best.)));
 run; %put &=numDS.;

/*%macro build_list;*/
/*  libname resread (res) access=read;*/
/**/
/*  proc datasets library=work noprint;*/
/*   delete all_settings;*/
/*  run;*/
/*  quit;*/
/**/
/*  %do i = 1 %to &numDS.;*/
/*   proc sort data =  res.&&ds&i..(keep = perturbation) out = temp nodupkey;*/
/*    by perturbation;*/
/*   run;*/
/**/
/*   data temp;*/
/*    length ds $100.;*/
/*    set temp;*/
/*	 ds = "res.&&ds&i..";*/
/*   run;*/
/**/
/*   proc append data = temp base = all_settings force; run; quit;*/
/*  %end;*/
/**/
/*%mend;*/
/*%build_list;*/


%macro stack(dsout=results,timing=2,extra_events=306);
 data results;
  set %do i = 1 %to &numDS.;
            readonly.&&ds&i..(keep=hypothesis perturbation extra_events sig1 interimTime1
			              sig&timing. interimTime&timing. LRStat_a0&timing.)
      %end;;
  where extra_events=&extra_events;

  rename sig1                = sig_second
         interimTime1        = iTime_second
         sig&timing.         = sig_first
		 interimTime&timing. = iTime_first
		 LRStat_a0&timing.   = w;
 run;
%mend;

data design_parameters;
alpha = 0.025;
power = 0.900;

array t1e_bound[3] _temporary_ (0.025 0.050 0.075);
array pwr_bound[3] _temporary_ (0.050 0.100 0.150);

 do extra_events = 0,153,306;
 do timing  = 2 to 5;
 do bounds  = 1 to 3;

 node_idx + 1;
 delta_e = t1e_bound[bounds];
 delta_p = pwr_bound[bounds];
  output;

 end;
 end;
 end;
run; 

/*data design_parameters;*/
/* set design_parameters;*/
/* where node_idx = 2;*/
/* call symput('timing', strip(put(timing, best.)));*/
/* call symput('delta_e',strip(put(delta_e,best.)));*/
/* call symput('delta_p',strip(put(delta_p,best.)));*/
/* call symput('alpha',  strip(put(alpha,best.)));*/
/* call symput('power',  strip(put(power,best.)));*/
/**/
/*run;*/
/**/
/*%stack(dsout=results,timing=&timing.);*/


%macro identify_design(dsOut=Design,w0=0.0,delta_w=0.1);

data design_parameters;
 set design_parameters;
 where node_idx = &sysparm.;

 call symput('extra_events',strip(put(extra_events, best.)));
 call symput('timing', strip(put(timing, best.)));
 call symput('delta_e',strip(put(delta_e,best.)));
 call symput('delta_p',strip(put(delta_p,best.)));
 call symput('alpha',  strip(put(alpha,best.)));
 call symput('power',  strip(put(power,best.)));

run;

%stack(dsout=results,timing=&timing.,extra_events=&extra_events.);

%let var =
 prob_stoppage w marg_dur cond_dur marg_dur_perc cond_dur_perc rejRate;

 %let terminate = 0;
 %let last      = 0;

 %do %until(%eval(&terminate=1 & &last=1) );

  data temp;
   set results;

   prob_stoppage = (w <= &w0.);
   rejRate       = prob_stoppage*sig_first   + (1-prob_stoppage)*sig_second;
   marg_dur      = prob_stoppage*iTime_first + (1-prob_stoppage)*iTime_second;
   marg_dur_perc = (iTime_second - marg_dur) / iTime_second * 100;

   cond_dur      = iTime_first;
   cond_dur_perc = (iTime_second - cond_dur) / iTime_second * 100;

   keep extra_events hypothesis perturbation prob_stoppage w 
        marg_dur cond_dur marg_dur_perc cond_dur_perc rejRate;
  run;

  proc means data = temp noprint nway;
   class extra_events hypothesis perturbation;
   var &var.;
   output out = temp n(w)=n mean=&var.;
  run;

  %let rejRate_null = 0.00;
  %let rejRate_alt  = 1.00;
  data null alt;
   set temp;
    if hypothesis = 0 then output null;
    if hypothesis = 1 then output alt;
  run;
  proc sql noprint;
   select max(max(rejRate),0.00) into :rejRate_null from null;
   select min(min(rejRate),1.00) into :rejRate_alt  from alt;
  quit; 
  %let rejRate_null = %left(&rejRate_null.);
  %let rejRate_alt  = %left(&rejRate_alt.);




  %if &terminate = 0 %then %do;
   data _null_;
    terminate = &terminate.;
     if &rejRate_null. > &alpha. + &delta_e. then terminate = 1;
     if &rejRate_alt.  < &power. - &delta_p. then terminate = 1;
     call symput('terminate',strip(put(terminate,best.)));

   run;
   %put &=rejRate_null  - &=rejRate_alt - &=terminate - &=last - &w0.;

   data _null_;
     if &terminate = 0 then call symput('w0',strip(put(&w0.+&delta_w.,best.)));
	 else call symput('w0',strip(put(&w0.-&delta_w.,best.)));
   run;

  %end;
  %else %do;
   %let last = 1;
  %end;

  

  %if &terminate = 1 & &last = 1 %then %do;
  %put &=rejRate_null  - &=rejRate_alt - &=terminate - &=last - &w0.;
   data res.&dsout._%sysfunc(putn(&sysparm.,z5.));
    set temp;
	 w0     = &w0.;
	 exp_w0 = round(exp(w0),1e-4);

	 target_alpha   = &alpha.;
	 target_power   = &power.;
	 delta_e        = &delta_e.;
	 delta_p        = &delta_p.;
	 timing         = &timing.;
	 infoReduction  = (timing-1)/8*100;

	 drop _type_ _freq_;
   run;
  %end;

 
  %end;

  proc datasets library = work noprint; 
   delete alt null temp;
  quit;

%mend;

option nonotes;
%identify_design(dsOut=Design,w0=0.0,delta_w=0.05);
option notes;


