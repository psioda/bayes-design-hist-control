
%macro setup;
proc datasets library=work noprint kill; run; quit;
%global root outpath;

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
  %let sysparm = 3;
  %let outPath = &root.\data\simulation_results\;

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

  %* macro variables;
  %let outPath = &root./data/simulation_results/;

  %* macros;
  %include "&root./../macros/collapse.sas";
  %include "&root./../macros/fit_approx.sas";
  %include "&root./../macros/fit_mcmc.sas";
  %include "&root./../macros/simulateStudy.sas";
  %include "&root./../macros/subsetInterim.sas";
%end;

%mend setup;
%setup;

ods noresults;

data res;
 set res.sim_bound:;
 by hypothesis perturbation;
run;

proc means data = res noprint nway;
 class hypothesis perturbation;
 weight nSims;
 var rejRate rejRateInt w;
 output out = res n(w) = n mean=rejRate rejRateInt w;
run;

data res2;
 set res;

 pHaz1  = scan(perturbation,1,'-');
 pHaz2  = scan(perturbation,2,'-');
 pBeta1 = scan(perturbation,3,'-');
 pBeta2 = scan(perturbation,4,'-');


 if pHaz1 = pHaz2 and pHaz1 ^= '1.00' and pBeta1 = pBeta2 and pBeta2 = '1.00' then perturbCat = 1;
 else perturbCat = 2;

  wStat = log(w);

  if perturbCat = 1 then wStat1 = wStat;
  if perturbCat = 2 then wStat2 = wStat;


       if perturbCat = 1 and pHaz1 > '1.00' then biasCat = 1;
  else if perturbCat = 1 and pHaz1 = '1.00' then biasCat = 2;
  else if perturbCat = 1 and pHaz1 < '1.00' then biasCat = 3;


  scenerio = 'Always Stop Early';
  r = rejRateInt; r1=r;
   output;

  scenerio = 'Adaptive';
  r = rejRate; r1=.; r2=r;
   output;
run;

ods results;
ods html close;
options Papersize=("6in","8.3in") nodate nonumber;
ODS PDF file = "&outPath.SIM_bounding_t1e.pdf" dpi=500 startpage=no;

ods graphics / height=4in width=6in noborder;
proc sgplot data = res2;
where hypothesis = 0 and biasCat in (3,.);
 scatter y=r1 x=wStat2 / markerattrs=(color=Gray  size=6 symbol=circlefilled) legendlabel='Other Perturbation';
 loess y=r1 x=wStat1  / CLM="Confidence Limits" clmtransparency=0.8 datalabel=pHaz1 lineattrs=(color=veryDarkGray) markerattrs=(color=veryDarkGray size=8 symbol=diamondFilled) legendlabel='Uniform Perturbation';
  yaxis label='Type I Error Rate (Always Stop at Interim)';
  xaxis label='E[W|D]' values=(0 to 1 by 0.1);
run;

ods graphics / height=4in width=6in noborder;
proc sgplot data = res2;
where hypothesis = 0 and biasCat in (3,.);
 scatter y=r2 x=wStat2 / markerattrs=(color=Gray  size=6 symbol=circlefilled) legendlabel='Other Perturbation';
 loess y=r2 x=wStat1  / CLM="Confidence Limits" clmtransparency=0.8 datalabel=pHaz1 lineattrs=(color=veryDarkGray) markerattrs=(color=veryDarkGray size=8 symbol=diamondFilled) legendlabel='Uniform Perturbation';
  yaxis label='Type I Error Rate (Adaptive Design)';
  xaxis label='E[W|D]' values=(0 to 1 by 0.1);
run;

ods pdf close;


ods html close;
options Papersize=("6in","8.3in") nodate nonumber;
ODS PDF file = "&outPath.SIM_bounding_pwr.pdf" dpi=500 startpage=no;

ods graphics / height=4in width=6in noborder;
proc sgplot data = res2;
where hypothesis = 1 and biasCat in (1,.);
 scatter y=r1 x = wStat2 / markerattrs=(color=Gray  size=6 symbol=circlefilled) legendlabel='Other Perturbation';
 loess y=r1 x = wStat1  / CLM="Confidence Limits" clmtransparency=0.8 datalabel=pHaz1 lineattrs=(color=veryDarkGray) markerattrs=(color=veryDarkGray size=8 symbol=diamondFilled) legendlabel='Uniform Perturbation';
  yaxis label='Power (Always Stop at Interim)';
  xaxis label='E[W|D]' values=(0 to 1 by 0.1);
run;

ods graphics / height=4in width=6in noborder;
proc sgplot data = res2;
where hypothesis = 1 and biasCat in (1,.);
 scatter y=r2 x = wStat2 / markerattrs=(color=Gray  size=6 symbol=circlefilled) legendlabel='Other Perturbation';
 loess y=r2 x = wStat1  / CLM="Confidence Limits" clmtransparency=0.8 datalabel=pHaz1 lineattrs=(color=veryDarkGray) markerattrs=(color=veryDarkGray size=8 symbol=diamondFilled) legendlabel='Uniform Perturbation';
  yaxis label='Power (Adaptive Design)';
  xaxis label='E[W|D]' values=(0 to 1 by 0.1);
run;

ods pdf close;


