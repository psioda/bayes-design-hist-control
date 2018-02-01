
%macro setup;
proc datasets library=work noprint kill; run; quit;
%global root outpath;

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
  %let root = /nas/longleaf/home/psioda/stat-projects/bayesDesignCVOT/CVOT;

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

data designs;
 set res.design_:;

  AmountBorrowed = strip(put(infoReduction,5.2))||'%';
  modFactor      = input(scan(perturbation,2,'()'),best.);
  bounds         = put(target_alpha+delta_e,5.3)||'/'||put(target_power-delta_p,5.3);

run;



data myattrmap;
length id value $20 MARKERCOLOR linecolor $ 20 linepattern $ 9;
infile datalines dlm=' ';
input id $ value $ linecolor $ linepattern $;

MARKERCOLOR = linecolor;

datalines;
id  12.50%    black     1
id  25.00%    darkgray  2
id  37.50%    gray      1
id  50.00%    lightgray 2
;
run;

data myattrmap2;
length id value $20 MARKERCOLOR linecolor $ 20 linepattern $ 9;
infile datalines dlm=' ';
input id $ value $ linecolor $ linepattern $;

MARKERCOLOR = linecolor;

datalines;
id  0      black     1
id  153    veryDarkGray      2
id  306    gray 1
;
run;


** Type 1 Error Figure;
options Papersize=("8in","8in") nodate nonumber;
ods pdf file = "&outPath.t1e_p2A.pdf" dpi=500;

ods html close;
ods graphics / height=8in width=8in noborder;

proc sgpanel data = designs dattrmap=myattrmap;
 panelby bounds  extra_events / onepanel layout=lattice;
 where hypothesis = 0 and 0.6 <= modFactor <= 1.4;
 refline 0.00 0.025 0.05 0.075 0.10 / axis=y lineattrs=(pattern=3 color=verylightgray);
 loess x=modFactor y=rejRate / group = AmountBorrowed lineattrs=(thickness=2) markerattrs=(symbol=circleFilled size=6) attrid=id;
  label extra_events = 'Additional Events';
  label AmountBorrowed = 'Info. Reduction  ' 
        modFactor = 'Baseline Hazard Perturbation' 
        rejRate = 'Estimated Type I Error Rate' bounds = 'T1E/Power Bounds';
 rowaxis values=(0.000 to 0.10 by 0.02) ;
run;
quit;
run;

ods pdf close;


options Papersize=("8in","8in") nodate nonumber;
ods pdf file = "&outPath.t1e_p2B.pdf" dpi=500;

ods html close;
ods graphics / height=8in width=8in noborder;
proc sgpanel data = designs dattrmap=myattrmap2;
 panelby bounds  AmountBorrowed / onepanel layout=lattice;
 where hypothesis = 0 and 0.6 <= modFactor <= 1.4;
 refline 0.00 0.025 0.05 0.075 0.10 / axis=y lineattrs=(pattern=3 color=verylightgray);
 refline 1.00 / axis=x  lineattrs=(pattern=3 color=verylightgray);
 loess x=modFactor y=rejRate / nomarkers group = extra_events 
                                lineattrs=(thickness=2 ) /*markerattrs=(symbol=circleFilled size=6 )*/ attrid=id;
  label extra_events = 'Additional Events';
  label AmountBorrowed = 'Info. Reduction  ' 
        modFactor = 'Baseline Hazard Perturbation' 
        rejRate = 'Estimated Type I Error Rate' bounds = 'T1E/Power Bounds';
 rowaxis values=(0.000 to 0.10 by 0.02) ;
  keylegend / title='Additional Events ' ;
run;
quit;
run;
ods pdf close;




** Power Figure;
options Papersize=("8in","8in") nodate nonumber;
ods pdf file = "&outPath.power_p2A.pdf" dpi=500;

ods html close;
ods graphics / height=8in width=8in noborder;

proc sgpanel data = designs dattrmap=myattrmap;
 panelby bounds extra_events / onepanel layout=lattice;
 where hypothesis = 1 and 0.6 <= modFactor <= 1.4;
 refline 0.75 0.80 0.85 0.90 0.95 / axis=y lineattrs=(pattern=3 color=verylightgray);
 refline 1.00 / axis=x  lineattrs=(pattern=3 color=verylightgray);
 series x=modFactor y=rejRate / markers group = AmountBorrowed lineattrs=(thickness=2) markerattrs=(symbol=circleFilled size=6) attrid=id;
  label extra_events = 'Additional Events';
  label AmountBorrowed = 'Info. Reduction  ' 
        modFactor = 'Baseline Hazard Perturbation' 
        rejRate = 'Estimated Power' bounds = 'T1E/Power Bounds';
 rowaxis values=(0.75 to 0.95 by 0.05) ;
run;
quit;
run;
ods pdf close;

** Power Figure;
options Papersize=("8in","8in") nodate nonumber;
ods pdf file = "&outPath.power_p2B.pdf" dpi=500;

ods html close;
ods graphics / height=8in width=8in noborder;
proc sgpanel data = designs dattrmap=myattrmap2;
 panelby bounds  AmountBorrowed / onepanel layout=lattice;
 where hypothesis = 1 and 0.6 <= modFactor <= 1.4;
 refline 0.75 0.80 0.85 0.90 0.95 / axis=y lineattrs=(pattern=3 color=verylightgray);
 refline 1.00 / axis=x  lineattrs=(pattern=3 color=verylightgray);
 loess x=modFactor y=rejRate / nomarkers group = extra_events 
                                lineattrs=(thickness=2 ) /*markerattrs=(symbol=circleFilled size=6 )*/ attrid=id;
  label extra_events = 'Additional Events';
  label AmountBorrowed = 'Info. Reduction  ' 
        modFactor = 'Baseline Hazard Perturbation' 
        rejRate = 'Estimated Power' bounds = 'T1E/Power Bounds';
 rowaxis values=(0.75 to 0.95 by 0.05) ;
  keylegend / title='Additional Events ' ;

run;
quit;
run;
ods pdf close;

** Probability of Early Stoppage Figure;
options Papersize=("8in","8in") nodate nonumber;
ods pdf file = "&outPath.early_stoppage_p2A.pdf" dpi=500;

ods html close;
ods graphics / height=8in width=8in noborder;


proc sgpanel data = designs dattrmap=myattrmap;
 panelby bounds extra_events / onepanel layout=lattice;
 where hypothesis = 1 and 0.6 <= modFactor <= 1.5;
 refline 1.00 / axis=x  lineattrs=(pattern=3 color=verylightgray);
 series x=modFactor y=prob_stoppage / markers group = AmountBorrowed lineattrs=(thickness=2) markerattrs=(symbol=circleFilled size=6) attrid=id;
  label extra_events = 'Additional Events';
  label AmountBorrowed = 'Info. Reduction  ' 
        modFactor = 'Baseline Hazard Perturbation' 
        prob_stoppage = 'Probabiliity of Early Stoppage' bounds = 'T1E/Power Bounds';
 rowaxis values=(0.000 to 1.00 by 0.1) grid;
run;
quit;
run;
ods pdf close;

data temp1;
 set designs;
run; 
proc sort data = temp1; by bounds AmountBorrowed extra_events; run;
proc means data = temp1 noprint nway;
 class bounds AmountBorrowed extra_events;
 var prob_stoppage;
 output out = temp2 max=pMax;
run;

data temp2;
 set temp2;
 length pMaxc  $15;
 pMaxc = put(extra_events,3.)||'-'||put(pMax,4.2)||'   ';
 if pMaxc =: '  0' then pMaxc = '  '||pMaxc;
run;


proc transpose data = temp2 out = temp2B prefix=B;
 by bounds AmountBorrowed;
 id extra_events;
 var pMaxc;
run;

data temp3;
 merge temp1
             temp2B(keep=bounds AmountBorrowed B0 B153 B306);
  by bounds AmountBorrowed; 
run;


** Probability of Early Stoppage Figure;
options Papersize=("8in","8in") nodate nonumber;
ods pdf file = "&outPath.early_stoppage_p2B.pdf" dpi=500;

ods html close;
ods graphics / height=8in width=8in noborder;
proc sgpanel data = temp3 dattrmap=myattrmap2;
 panelby bounds AmountBorrowed  / onepanel layout=lattice;
 where hypothesis = 1 and 0.5 <= modFactor <= 1.5;
 refline 1.00 / axis=x  lineattrs=(pattern=3 color=verylightgray);
 loess x=modFactor y=prob_stoppage / nomarkers group = extra_events lineattrs=(thickness=2) /*markerattrs=(symbol=circleFilled size=6)*/ attrid=id;
  label extra_events = 'Additional Events';
  label AmountBorrowed = 'Info. Reduction  ' 
        modFactor = 'Baseline Hazard Perturbation' 
        prob_stoppage = 'Probabiliity of Early Stoppage' bounds = 'T1E/Power Bounds';
 rowaxis values=(0.000 to 1.00 by 0.2) grid;
  keylegend / title='Additional Events  ' ;
    inset  B0 B153 B306 / nolabel title="Maximum" position=topleft
                          titleattrs=(size=7);
run;
quit;
run;

ods pdf close;


** Marginal Mean Percent Reduction in Trial Length;
options Papersize=("8in","8in") nodate nonumber;
ods pdf file = "&outPath.trial_length_p2A.pdf" dpi=500;

ods html close;
ods graphics / height=8in width=8in noborder;


proc sgpanel data = designs dattrmap=myattrmap;
 panelby bounds  extra_events / onepanel layout=lattice;
 where hypothesis = 1 and 0.6 <= modFactor <= 1.4;
 refline 1.00 / axis=x  lineattrs=(pattern=3 color=verylightgray);
 series x=modFactor y=marg_dur_perc / markers group = AmountBorrowed lineattrs=(thickness=2) markerattrs=(symbol=circleFilled size=6) attrid=id;
  label extra_events = 'Additional Events';
  label AmountBorrowed = 'Info. Reduction  ' 
        modFactor = 'Baseline Hazard Perturbation' 
        marg_dur_perc = 'Average % Reduction in Trial Length' bounds = 'T1E/Power Bounds';
 rowaxis values=(0.000 to 25.00 by 5) grid;
run;
quit;
run;

ods pdf close;


data temp1;
 set designs;
run; 
proc sort data = temp1; by bounds AmountBorrowed extra_events; run;
proc means data = temp1 noprint nway;
 class bounds AmountBorrowed extra_events;
 var marg_dur_perc;
 output out = temp2 max=pMax;
run;

data temp2;
 set temp2;
 length pMaxc  $15;
 pMaxc = put(extra_events,3.)||'-'||strip(put(pMax,5.1))||'   ';
 if pMaxc =: '  0' then pMaxc = '  '||pMaxc;
run;


proc transpose data = temp2 out = temp2B prefix=B;
 by bounds AmountBorrowed;
 id extra_events;
 var pMaxc;
run;

data temp3;
 merge temp1
             temp2B(keep=bounds AmountBorrowed B0 B153 B306);
  by bounds AmountBorrowed; 
run;


** Marginal Mean Percent Reduction in Trial Length;
options Papersize=("8in","8in") nodate nonumber;
ods pdf file = "&outPath.trial_length_p2B.pdf" dpi=500;

ods html close;
ods graphics / height=8in width=8in noborder;


proc sgpanel data = temp3 dattrmap=myattrmap2;
 panelby bounds AmountBorrowed / onepanel layout=lattice;
 where hypothesis = 1 and 0.6 <= modFactor <= 1.4;
 refline 1.00 / axis=x  lineattrs=(pattern=3 color=verylightgray);
 loess x=modFactor y=marg_dur_perc / nomarkers group = extra_events  lineattrs=(thickness=2) /*markerattrs=(symbol=circleFilled size=6)*/ attrid=id;
  label extra_events = 'Additional Events';
  label AmountBorrowed = 'Info. Reduction  ' 
        modFactor = 'Baseline Hazard Perturbation' 
        marg_dur_perc = 'Average % Reduction in Trial Length' bounds = 'T1E/Power Bounds';
 rowaxis values=(0.000 to 25.00 by 5) grid;
 keylegend / title='Additional Events ' ;
    inset  B0 B153 B306 / nolabel title="Maximum" position=topleft
                          titleattrs=(size=7);
run;
quit;
run;

ods pdf close;
