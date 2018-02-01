/*
 This SAS macro is used back-calculate the interim analysis data using final data that was simulated with the %simulateStudy SAS macro.
 The macro requires the following parameters:
  
  MACRO VARIABLES REQUIRED:
  
   [1] dsFinal     = input subject-level dataset for final analysis opportunity;
   [2] nuInterim   = number of events at which to perform the interim analysis opportunity;
   [3] dsInterim   = output subject-level dataset for interim analysis opportunity;

*/

%macro subsetInterim;
proc sort data = &dsFinal.(drop=targetTime) out = __events__;
 by simStudy elapsedTime;
 where censor = 0;
run;

data __events1__;
 set __events__ end=last;
 by simStudy;
 if first.simStudy then nObs = 0;
 nObs+1;
 
  if nObs = &nuInterim. or (last.simStudy and nObs < &nuInterim.) then output;
   keep simStudy elapsedTime nObs;
    rename elapsedTime = targetTime;
run;


data __events2__;
 merge __events__ __events1__;
 by simStudy;
   keep simStudy elapsedTime targetTime;
    rename elapsedTime = nextTime;
run;

data __events3__;
 set __events2__;
 by simStudy;
 where nextTime > targetTime;
  if first.simStudy;
  keep simStudy nextTime;
 run;


data &dsInterim. interimTime(keep=simStudy modTargetTime rename=(modTargetTime=interimTime));
 merge &dsFinal.(drop=targetTime) __events1__ __events3__;
 by simStudy;

 if nextTime > . then modTargetTime = (targetTime + nextTime)/2;
 else modTargetTime = targetTime;
 
 if first.simStudy then output interimTime;
 
 if enrollTime >= modTargetTime  then delete;
 else if elapsedTime > modTargetTime then do;
   obsTime = modTargetTime - enrollTime;
   censor = 1;
   elapsedTime = enrollTime + obsTime;
 end;
 output &dsInterim.;

run;

 
 proc datasets library=work noprint;
  delete __event:;
 run;
 quit; 
 


%mend;