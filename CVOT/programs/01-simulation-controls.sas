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


%let number_of_nodes = 546;
%let nPerLoop        = 125;
%let nLoops          = 800;
%let numSubjects     = 5000;
%let nuTarget        = 612;
%let modStart        = 0.55;
%let modStop         = 1.45;
%let modBy           = 0.01;

** compute posterior means;
proc means data = sp.histPostSamp noprint nway;
 var lambda1-lambda4 beta1-beta6;
 output out = parameter_perturbations mean=h11 h21 h31 h32 beta1-beta6;
run;


data sp.parameter_perturbations;
 length perturbation $20. modFact 8.;
 set parameter_perturbations;

  gamma = 0;

  array h(*) h:;
  array b(*) beta:;

  do modFact = &modStart. to &modStop. by &modBy.;
  idx+1;
    perturbation = 'haz ('||strip(put(modFact,5.3))||')';
    do j = 1 to dim(h);
      h(j) = h(j)*modFact;
    end;
    output;
    do j = 1 to dim(h);
      h(j) = h(j)/modFact;
    end;
  end;
  keep perturbation modFact h: beta: gamma;
run;



data controls2;
 length hypothesis nPerLoop numSubjects 4.;
 set sp.parameter_perturbations(keep=perturbation);
 
 nPerLoop    = &nPerLoop.;
 nLoops      = &nLoops.;
 nSims       = nPerLoop*nLoops;
 numSubjects = &numSubjects.;

 do nuTotal          = &nuTarget.;
  do extra_events    = 0,153,306;
   do hypothesis     = 0,1;
     do parmSet      = 1 to 1;
      output; 
     end; 
   end;
  end;
 end;

run;
proc sort data = controls2 sortseq=linguistic(numeric_collation=on); by hypothesis perturbation parmSet; run;

data controls3;
 set controls2;
 node_idx+1;
 if node_idx > &number_of_nodes. then node_idx=1;
run;

data sc.simulation_controls;
 retain node_idx inner_idx;
 set controls3;

 do inner_idx = 1 to nLoops;
    seed           = round(1 + rand('uniform')*2**30);
    covSampleSeed  = round(1 + rand('uniform')*2**30);
    parmSampleSeed = round(1 + rand('uniform')*2**30);
    output;
 end;
 drop nLoops;
run;
