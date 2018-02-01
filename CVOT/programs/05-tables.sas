
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
ods csv file="&outPath.Table1.csv";

proc print data = sp.histEstMCMC noobs;
format Mean StdDev HPDLower HPDUpper 6.4;
var Parameter Mean StdDev HPDLower HPDUpper;
run;
ods csv close;


/*

Parameter Mean StdDev HPDLower HPDUpper 
Parameter Mean StdDev HPDLower HPDUpper 
BETA1      0.3860 0.0935  0.2008  0.5664 
BETA2      0.6465 0.1040  0.4465  0.8530 
BETA3      0.4798 0.0839  0.3167  0.6453 
BETA4      0.0823 0.0425  0.0007  0.1667 
BETA5      1.4804 0.2452  1.0082  1.9651 
BETA6     -0.0141 0.0020 -0.0181 -0.0101 
LAMBDA_1_1 0.0164 0.0017  0.0131  0.0198 
LAMBDA_2_1 0.0218 0.0023  0.0174  0.0262 
LAMBDA_3_1 0.0223 0.0038  0.0151  0.0300 
LAMBDA_3_2 0.0361 0.0054  0.0259  0.0468 



\begin{table}
\centering
\caption{Posterior summaries for SAVOR trial} \label{tab:PostSumSAVOR}
\vspace{0.1cm}
\begin{tabular}{llccc}
    \hline
     Parameter & Characteristic                        & Mean   & SD     & HPD                  \\
    \hline
    $\beta_1$  & Male                                  & ~0.3860 & 0.0935 & (~0.2008,~0.5664)   \\ 
    $\beta_2$  & History of Stroke                     & ~0.6465 & 0.1040 & (~0.4465,~0.8530)   \\ 
    $\beta_3$  & History of MI                         & ~0.4798 & 0.0839 & (~0.3167,~0.6453)   \\ 
	$\beta_4$  & log(Duration of Diabetes (yrs))       & ~0.0823 & 0.0425 & (~0.0007,~0.1667)   \\ 				
	$\beta_5$  & log(HBA1c \%)                         & ~1.4804 & 0.2452 & (~1.0082,~1.9651)   \\ 
    $\beta_6$  & Estimated GFR (mL/min/1.73m$^2$)      & -0.0141 & 0.0020 & (-0.0181,-0.0101)   \\ 	\hline		                       
    $\lambda_{1,1}$ & Age $\le$ 65                     & ~0.0164 & 0.0017 & (~0.0131,~0.0198)  \\ 
    $\lambda_{2,1}$ & 65 $<$ Age $\le$ 75              & ~0.0218 & 0.0023 & (~0.0174,~0.0262)  \\
    $\lambda_{3,1}$ & Age $>$ 75                       & ~0.0223 & 0.0038 & (~0.0151,~0.0300)   \\
    $\lambda_{3,2}$ &                                  & ~0.0361 & 0.0054 & (~0.0259,~0.0468)   \\ \hline
\end{tabular}
\end{table}


\begin{table}[h!]
\centering
\caption{Posterior summaries for the SAVOR trial} \label{tab:PostSumSAVOR}
\vspace{0.1cm}
\begin{tabular}{llccc}
    \hline
     Parameter & Characteristic                               & Mean    & SD     & HPD                  \\
    \hline
    $\beta_1$  & Male                                         & ~0.3860 & 0.0935 & (~0.2008,~0.5664)   \\ 
    $\beta_2$  & History of Stroke                            & ~0.6465 & 0.1040 & (~0.4465,~0.8530)   \\ 
    $\beta_3$  & History of MI                                & ~0.4798 & 0.0839 & (~0.3167,~0.6453)   \\ 
    $\beta_4$  & log[Duration of Diabetes (yrs)]              & ~0.0823 & 0.0425 & (~0.0007,~0.1667)   \\ 				
	$\beta_5$  & log[HBA1c (\%)]                              & ~1.4804 & 0.2452 & (~1.0082,~1.9651)   \\ 
    $\beta_6$  & eGFR (mL/min/1.73m$^2$)                      & -0.0141 & 0.0020 & (-0.0181,-0.0101)   \\ 	\hline		                       
    $\lambda_{1,1}$ - $[0,\infty)$     & Age $\le$ 65         & ~0.0164 & 0.0017 & (~0.0131,~0.0198)  \\ 
    $\lambda_{2,1}$ - $[0,\infty)$     & 65 $<$ Age $\le$ 75  & ~0.0218 & 0.0023 & (~0.0174,~0.0262)  \\
    $\lambda_{3,1}$ - $[0,1.04)$      & Age $>$ 75            & ~0.0223 & 0.0038 & (~0.0151,~0.0300)   \\
    $\lambda_{3,2}$ - $[1.04,\infty)$ &                       & ~0.0361 & 0.0054 & (~0.0259,~0.0468)   \\ \hline
\end{tabular}
\end{table}

*/

data temp;
 set res.design_:;
  AmountBorrowed = strip(put(infoReduction,5.2))||'%';
  modFactor      = input(scan(perturbation,2,'()'),best.);
  bounds         = put(target_alpha+delta_e,5.3)||'/'||put(target_power-delta_p,5.3);
run; proc sort nodupkey; by Extra_Events bounds AmountBorrowed; run;

proc transpose data = temp out = temp2 prefix=p;
 by Extra_Events bounds;
 id AmountBorrowed;
 var exp_w0;
run;

data _null_;
 set temp2;

 ec = put(Extra_Events,3.);
 ec = tranwrd(ec," ","~");

 p12_50 = put(p12_50_,4.1);
 p12_50 = tranwrd(p12_50," ","~");

 p25_00 = put(p25_00_,4.1);
 p25_00 = tranwrd(p25_00," ","~");

 p37_50 = put(p37_50_,4.1);
 p37_50 = tranwrd(p37_50," ","~");

 p50_00 = put(p50_00_,4.1);
 p50_00 = tranwrd(p50_00," ","~");

  put "$" ec  "$ & $" bounds "$ & $" p12_50 "$ & $"  
                                              p25_00 "$ & $" 
                                              p37_50 "$ & $"
											  p50_00 "$ \\";
run;


