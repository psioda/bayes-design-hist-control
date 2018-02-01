/*
 This SAS macro is used simulate data for the new trial. The macro simulates data from a proportional hazards model with piecewise constant
 baseline hazard.
 The macro requires the following parameters:
  
  MACRO VARIABLES REQUIRED:
  
     [1] stratLevels     = number of strata;
   [2] trueHazComp     = list of the number of baseline hazard components in each stratum - each stratum separated by a pipe;
   [3] covSampleSeed   = seed used by PROC SURVEYSELECT when sampling covaraite vectors using unrestricted sampling;
   [4] dsCovDist       = dataset of covariate vectors to be sampled from;
   [5] numSimulations  = number of datasets to be simulated;
   [6] numSubjects     = number of subjects to enroll;
   [7] parmSampleSeed  = seed used by PROC SURVEYSELECT when sampling parameter values from sampling prior dataset using 
                         unrestricted random sampling;
   [8] sampPrior       = dataset that contains the parameter values that define the generative model; This code assumes that the
                         dataset has one observation as it is taken from a more general macro;
   [9] seed            = seed value used in random number generator when simulating new trial data;
  [10] dsBreaks        = dataset that stores the baseline hazard change points - required as input;
  [11] hypothesis      = an indicator for whether the null hypothesis (HR = 1.3) or alternative (HR=1.0) is true; (i.e., this code is not general)
  [12] rProb           = probability of randomization to treatment;
  [13] enrollDist     = distribution for enrollment - should be valid call to RAND function or a constant;
  [14] censorDist     = distribution for censorship - should be valid call to RAND function or a constant;
  [15] nuTarget       = targeted number of events for new trial;
  [16] dsOut          = dataset to store simulated new trial data;
  [17] DEBUG
*/

%macro simulateStudy;
%let maxComp = 0;
%do s = 1 %to &stratLevels.;
  %if %eval(&maxComp < %scan(&trueHazComp.,&s,|)) %then %let maxcomp = %scan(&trueHazComp.,&s,|);
%end; 


proc surveyselect data = &dsCovDist.(keep= riskStratum x:) out = __covariates(drop=NumberHits) seed=&covSampleSeed. rep=&numSimulations.
          outhits method=urs  n=&numSubjects. noprint; 	   
run; quit;

proc surveyselect data = &sampPrior.(keep= h: beta: gamma) out = __parameters(drop=NumberHits) seed=&parmSampleSeed. rep=&numSimulations.
          outhits method=urs n=1 noprint; 		   
run; quit;

data __simulated_subjects;
 merge __covariates __parameters;
 by replicate;
run;

data __simulated_subjects;
 set __simulated_subjects;
 by replicate;
  if first.replicate then subID = 0;
  subID+1;

 rename replicate = simStudy;
run;

data __temp__;
 call streaminit(&seed.);
 set __simulated_subjects;
 
 
   ** read in hazard break points;
   if _n_ = 1 then set &dsBreaks.;
  
   ** put parameters and covariates in array for processing;
   array h(&stratLevels.,&maxComp.) %do s = 1 %to &stratLevels.; %do c = 1 %to &maxComp.; h&s.&c. %end; %end; ;
   array t(&stratLevels.,&maxComp.) %do s = 1 %to &stratLevels.; %do c = 1 %to &maxComp.; t&s.&c. %end; %end;;  
  
   array x(*) x:;
   array beta(*) gamma beta: ;

   
   ** null or alternative value of treatment effect;
   if &hypothesis = 0 then beta(1) = log(1.3);
   else beta(1) = log(1.0);   
   

   ** simulate treatment allocation;
   x1 = rand('Bernoulli',&rProb.);
   
   ** calculate subject hazard ratio regression function;
   logHR = 0;
   do j = 1 to dim(beta);
      logHR = logHR + x(j)*beta(j);
   end;
   HR = exp(logHR);

   ** simulate enrollment time;
   enrollTime = round(&enrollDist.,0.0001);

   ** simulate censorship time;
   censTime = round(&censorDist.,0.0001);

   ** Simulate data from Piecewise-exponential model;
   k       = 1;
   obsTime = 0;
   stop    = 0;
   if t(riskStratum,k) > . then do while(stop=0);
     haz = h(riskStratum,k)*HR;
     obsTime   = obsTime + rand('exponential')/haz;

     if obsTime > t(riskStratum,k) and t(riskStratum,k) > . then do; obsTime = t(riskStratum,k); k+1; end;
     else stop = 1;
   end;
   else do;
     haz = h(riskStratum,k)*HR;
     obsTime   = obsTime + rand('exponential')/haz;
   end;
   
    obsTime = round(obsTime,0.0001);
  
    if obsTime > censTime then do; obsTime = censTime; censor = 1; end;
    else censor = 0;

    elapsedTime = enrollTime + obsTime;

   
 keep obsTime enrollTime elapsedTime x: censor riskStratum subID simstudy;
run;

proc sort data = __temp__ out = __events__;
 by simstudy elapsedTime;
 where censor = 0;
run;

data __events1__;
 set __events__ ;
 by simstudy;
 
  if first.simstudy then eventCounter = 0; eventCounter +1;
  
  if eventCounter = &nuTarget. or (last.simstudy and eventCounter < &nuTarget.) then output;
run;

data __events2__;
 set __events1__;
  targetTime = elapsedTime + 0.00005;
  keep simstudy targetTime;
run;

data &dsOut.;
 merge __temp__ __events2__;
 by simstudy;
 
 if enrollTime >= targetTime then delete;
 else if elapsedTime > targetTime then do;
   obsTime = targetTime - enrollTime;
   censor = 1;
   elapsedTime = enrollTime + obsTime;
 end;

run;

 %if &debug ^= 1 %then %do;
 proc datasets library=work noprint;
  delete __temp__ __events__ __events1__ __events2__ __covariates __parameters __simulated_subjects;
 run;
 quit; 
 %end;
%mend simulateStudy;
