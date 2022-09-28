
/*************************************************************************************************************/
/*                                       MACRO PSPMCM  version 1.1 (feb. 2006)                               */ 
/*************************************************************************************************************/
/*  This program is freely distributed.                                                                      */
/*  It is proposed under the only responsability of its user.                                                */
/*  In particular, the authors cannot be considered as responsible for any dysfunction due to the program use*/
/*                                                                                                           */
/* Authors : Fabien CORBIERE, Pierre JOLY (p.joly@isped.u-bordeaux2.fr)                                      */                                                                  */                                                  */ 
/*                                                                                                           */
/*************************************************************************************************************/
/* Cure models have been developped to analyse failure time data with a cured fraction. For such data,       */
/* standard survival models are usually not appropriate because they do not account for the possibility of   */
/* cure. Mixture cure models assume that the studied population is a mixture of susceptible individuals, who */
/* may undergo the event of interest, and non susceptible individuals that will never undergo it.            */
/* The aim of the PSPMCM SAS macro is to fit parametric and semiparametric mixture cure models for indididual*/
/* data with covariates. The cure fraction can be modelled by various binary regression models. Parametric   */
/* and semiparametric Cox proportional hazard models can be used to model the survival of uncured individuals*/
/* The maximization of the likelihood function is performed using SAS PROC NLMIXED for parametric models and */
/* through an EM algorithm for the Cox PH mixture cure model. The variance of parameters estimates is        */
/* obtained by inverting the Hessian matrix or by non parametric bootstrap methods.                          */
/*                                                                                                           */
/*************************************************************************************************************/
/*                                      REQUIRED ARGUMENTS                                                   */
/*************************************************************************************************************/
/* %PSPMCM                                                                                                   */
/*				 (DATA=													                                     */
/*                ID=,													                                     */
/*                CENSCOD=,											                                         */		
/*                TIME=,												                                     */
/*				  VAR=,												                                         */	 
/*                INCPART=,	                                                                                 */
/*                SURVPART=,                                                                                 */
/*                ALPHA=,                                                                                    */
/*                BASELINE=,                                                                                 */ 
/*                TAIL=,                                                                                     */
/*                SU0MET=                                                                                    */
/*				  MAXITER=,									                                      			 */
/*                CONVCRIT=,																	             */ 
/*                FAST=,           								                                             */
/*				  BOOTSTRAP=,																     			 */
/*                STRATA=,                                  								                 */
/*                NSAMPLE=,																	    			 */ 
/*                BOOTMET=,														        					 */
/*                JACKDATA=,                   								                                 */
/*                GESTIMATE=,													        					 */
/*                SPLOT= )       								                                             */
/*                                                                                                           */
/* DATA           specifies the data set you are analysing.	    						                         */
/* ID             is a data set variable which identify individuals.                                      */
/* CENSCOD        is a data set variable which identify censored records (=0) and observed failure time (=1) */
/* TIME           is a data set variable, in numeric format, which identify failure or censoring time        */  
/* VAR            Contains the name of covariates, separated by blancks. Each variable is followed, into     */
/*                brackets, by the option I, S or IS, which indicates whether it has to be included in the   */ 
/*                covariates vector in the incidence part (I), survival part (S), or both of them (IS)       */
/*                In addition, when plots of survival functions are requested  the value of covariates at    */
/*                which the survival estimates are plotted must be specified after a comma.                  */                                                                            */
/*                example : GENDER(I S,1) STUDY(I,2) AGE(S, .)                                               */
/*                The effect of GENDER will be estimated in the incidence part (effect on the cure fraction) */
/*                and on the latency part (effect on survival of uncured individuals).                       */
/*                The effect of STUDY will be estimated only in the incidence part.                          */
/*                Conversely, the effect of AGE will be estimated only in the latency part.                  */ 
/*                The marginal and conditional survivals functions estimates will be plotted for individuals */
/*                with GENDER=1 and STUDY=2. To specify that an explanatory variable in not taken into account*/
/*                for these plots, just insert a dot (.) (as in the example for AGE)                          */
/* INCPART        Models the incidence part. the logit(LOGIT), probit(PROBIT) or complementary loglog(CLOGLOG)*/
/*                links are available. By default the logistic regression model (LOGIT) is assumed.           */
/* SURVPART       Specifies the form of the conditional baseline survival function. Parametric models, i.e.   */
/*                exponential(EXP), Weibull(WEIB), loglogistic(LLOG) and lognormal(LOGN, and semiparametric   */
/*                Cox's proportional hazards model (COX) are available.                                       */
/* ALPHA          Sets the significance level used for the confidence limits for the hazards ratios and odd's */
/*                ratio. The default value is 0.05.                                                           */
/* BASELINE       Y/N. If set Y indicates that the conditional baseline survival function estimate and        */
/*                parameters estimates are outputted in the BASELINE dataset. If bootstrap resampling is      */
/*                requested for the Cox model, the BASELINE\_T dataset will moreover contain the estimates    */
/*                for all replicates. By default it is assumed to be N (no ouput).                            */
/* SPLOT          if set to Y, the marginal and conditional survival functions estimates are plotted for      */
/*                values given in the VAR parameter. The Kaplan Meier estimates based on the observed data are*/
/*                also displayed for categorical variables.                                                   */
/* PLOTFIT        Y/N. When set to Y, the macro computes for each stratum defined by the covariate vectors Xi */
/*                and Zi the observed (empirical) marginal survival curve which is the kaplan-Meier estimate  */
/*                for the whole stratum (including right censored subjects). The plots of observed and fitted */
/*                values versus is a visual tool to examine the goodness of the model prediction.             */
/*                The correlation coefficient between observed and fitted values is also computed for each    */
/*                stratum and P-P plots are also plotted.                                                     */
/*                                                                                                            */    
/*                                                                                                            */
/* If the Cox PH mixture cure model is requested, additional options are available    :                       */
/*                                                                                                            */
/* SU0MET        indicates whether the breslow-type method (CH) or the product limit estimator (PL) is used to*/
/*               estimate the conditional baseline survival function.                                         */
/* TAIL          indicates if a constraint or a tail completion method is used for the estimate of the        */
/*               baseline function. The option TAIL=ZERO specifies that the zero tail constraint is used      */
/*               ETAIL or WTAIL specify that the exponential and Weibull tail completion methods are used     */
/*               NONE indicates that no tail constraint is used (may however yield some convergence problems. */
/*               The default value is ZERO.                                                                   */
/* MAXITER       is the maximum number of iterations to perform. By default, MAXITER=200. If convergence is   */
/*               not attained the displayed output and all output data sets created by the procedure contain  */
/*               results that are based on the last maximum likelihood iteration.                             */
/* CONVCRIT      sets the convergence criterion. The default value is 10E-5. The iterations are considered to */
/*               have converged when the maximum relative change in the parameters and likelihood estimates   */
/*               between iteration steps is less than the value specified.                                    */
/* FAST          if set to Y, outputs the results when convergence is attained.Parameters estimates and their */
/*               standard errors (computed by inverting the matrix of second derivates) are displayed.        */
/* BOOTSTRAP     if set to Y, indicates that non-parametric resampling with replacement from the original data*/
/*               set is performed. The default value is N.                                                    */
/* NSAMPLE       is the number of replicates produced for bootstrap confidence intervals computation.         */
/* STRATA        specifies the variables used to performe stratified non parametric resampling.               */
/* BOOTMET       specifies which type of bootstrap confidence intervals are computed : percentile (PTCL)      */
/*               hybrib method (HYB), normalized bias corrected (BOOTN), bias corrected (BC) and accelerated  */
/*               bias corrected (BCA) confidence intervals. Jackknife after bootstrap (JACK) is also available*/
/*               The option BOOTMET=ALL specifies that all confidence intervals will be computed              */
/* JACKDATA      specifies the name of the SAS data set containing the previously computed jacknife estimates */
/*               may be useful to fasten the computation of BCA confidence interval.                          */   
/* GESTIMATE     if set to Y, the macro displays Q-Q plots and bar charts of the distribution of parameters   */
/*               estimates for the bootstrap replicates.                                                      */
/*                                                                                                            */
/* OUTPUT DATASETS                                                                                            */
/* The following dataset are created by the macro:                                                            */
/*                                                                                                            */
/* For parametric mixture cure models :                                   								      */
/*	ESTIMATES      contains the estimates and standard errors for parameters. The prefix 'L' indicates that   */
/*                 the estimate is for the incidence part, while the prefix "S" indicates that the estimate is*/ 
/*                 for the survival part.                                                                     */ 
/*                                                                                                             */
/* For semiparametric mixture cure models :                                                                    */
/*	LIKELIHOOD     stores information about convergence and likelihood.                                        */
/*  FAST_INCI      and FAST_INCI_OD are created if the parameter FAST=Y is specified. They store the estimates */
/*                 and standard errors and the corresponding odds ratios estimates convergence was attained.   */
/*	FAST_SURV      does the same for the latency part.                                                         */                    
/*  PROBCURE       stores the estimated probability of been uncured.                                           */
/*  BASELINE       is created if the parameter BASELINE=Y is specified. Stores the estimate for the conditional*/
/*                 baseline survival function, as well as parameters estimates.                                */
/*	BASELINE_T     is created if the parameters BASELINE=Y and BOOTSTRAP=Y are specified. Stores the estimates */
/*                 for parameters and the conditional baseline survival function from the bootstrap replicates.*/
/*  BOOTCI         is created if the parameter BOOTSTRAP=Y is specified. Stores the bootstrap confidence       */
/*                 intervals for parameters estimates for the incidence. The BOOT_OR and BOOT_RR dataset store */
/*                 the corresponding bootstrap confidence intervals for the Odds Ratios and the Hazards Ratios.*/
/*  BOOTDIST_T     is created if the parameter BOOTSTRAP=Y is specified. Stores the parameters estimates and   */
/*                 some information about convergence for each bootstrap replicate.                            */
/*  JACKDIST_T     is created if the parameter BOOTMET=JACK or BOOTMET=ALL is specified. Stores the parameters */
/*                 estimates and some information about convergence for each jackknife replicate.              */
/**************************************************************************************************************/;

/*************************************************************************/
/*                   MODEL AND COMPUTATION SPECIFICATIONS                */ 
/*************************************************************************/

%macro prepare;
%global nobs nfail list_err err nlvar nsvar nvar s totvar exiterr SURVPART2 lvars svars KM_var phpar Kmvar;
%local i i1 i2 i3 i4 dsid var_dat u exiterr_0 exiterr_1 exiterr_2 exiterr_3 exiterr_4;

%let list_err=;
%let exiterr_0=0;
%let exiterr_1=0;
%let exiterr_2=0;
%let exiterr_3=0;
%let exiterr_4=0;
%let exit_bas=0;
%let exit_first=0;
%let exit_spec=0;
%let SURVPART2=;
%let SURVPART3=;
%let _first=;
%let _bas=;
%let _spec=;
%let phparl=;
%let phpars=;
%let kmpar=;
%let kmvar=;
%let lvars=;
%let svars=;
%let nlvar=0;
%let nsvar=0;
%let nvar=0;
%let totvar=;
%let boot_err= ;

/* Check for missing information */

%if &data= %then %do;
	%let exiterr_0=1;
	%put ** No data set specified : the macro will stop.; 
%end;
%if %index(&SURVPART,COX) %then %do; 
 %let SURVPART=COX;
 %let SURVPART3= PH COX;
 %let SURVPART2= ;
%end;

%if %index(&SURVPART,COX) ne 1 %then %do;
	%let SURVPART2=PARA;
	%if %index(&SURVPART,EX) %then %let SURVPART=EXP;
	%if %index(&SURVPART,W) %then %let SURVPART=WEIBULL;
	%if %index(&SURVPART,LLOG) %then %let SURVPART=LLOGISTIC;
	%if %index(&SURVPART,LOGN) %then %let SURVPART=LOGNORMAL;
	%if &AFT=Y %then %let SURVPART3=AFT &SURVPART;
	%if &AFT ne Y %then %let SURVPART3=PH &SURVPART;
%end;

%if &SURVPART3= %then %do;
	%let exiterr_0=1;
	%put ** The model is misspecified : The macro will stop;
	%put ** Available models are : WEIBULL, EXP,LLOG,LOGN,COX.;
%end;
%if &CENSCOD= %then %do;
	%let exiterr_0=1;
	%put ** No censure variable specified : the macro will stop.; 
%end;
%if &TIME= %then %do;
	%let exiterr_0=1;
	%put ** No time variable specified : the macro will stop.; 
%end;
%if ((&LINK ne LOGIT) and (&LINK ne CLOGLOG) and (&LINK ne PROBIT)) %then %do;
	%let exiterr_0=1;
	%put ** The link function &LINK for the incidence part is misspecified: the macro will stop.;
	%put NOTE : link function can be : LOGIT, CLOGLOG, PROBIT.;
%end; 

%if (&SURVPART=COX and (&TAIL NE NONE) and (&TAIL ne WTAIL) and (&TAIL ne ETAIL) and (&TAIL NE ZERO)) %then %do;
	%let exiterr_0=1;
	%put ** Tail completion asked = &TAIL (allowed :  NONE, ZERO,ETAIL or WTAIL);
%end;

%if &exiterr_0=1 %then %goto final;

/*Check for the presence of data sets */

 %if %sysfunc(exist(&data))=0 %then %do ;
	%let exiterr_1=1; 
	%goto final;
%end;

/*  Number and list of specified covariates, specification in incidence and/or latency */
%let i=1;
%do %while (%length(%scan(&VAR,&i,(,) )));
	%let _kmd=0;
	%let _bas2=;
	%let _first=%scan(&VAR,&i, (,) ); 
	%if %datatyp(&_first)=NUMERIC %then %do ;
		%let exit_first=1; 
		%put *** Error in covariates specification ****;
		%put *** No baseline value should be specified if the SPLOT or PLOTFIT macro options are not set to Y ****;
		%goto final;
	%end;
	%let _spec=%qcmpres(%scan(&VAR,%eval(&i+1),(,) ));
	%if %datatyp(&_spec)=NUMERIC %then %do; 
		%let exit_spec=1;	
		%put *** Error in variables specifications ;
		%put *** Please check if information for incidence (I) or survival (S) are well specified into brackets*****;
		%put *** Please check that value for survival plot are well specified after a comma ***;
		%goto final;
	%end;
	%if &splot=Y %then %do;
		%let _bas=%qcmpres(%scan(&VAR,%eval(&i+2), (,) )); 
		%let i2=%eval(&i+3);
		%let i=%eval(&i+3);
	%end;
	%else %do;
		%let i2=%eval(&i+2);
		%let i=%eval(&i+2);
	%end;
	%do %while(%length(%scan(&VAR,&i2,(,) )));
		%let _next= %scan(&VAR,&i2,(,)); 
		%if %quote(&_first) = %quote(&_next) %then 	%let i2=10000; 
        %else %do;
         %if &splot=Y %then %let i2=%eval(&i2+3); 
 		 %else %let i2=%eval(&i2+2);
		%end;              
	%end;
	%if (&i2<10000) %then %do;
		%let exiterr_4=1;
		 %if %index(&_spec,I)  %then %do;
	   	 	%let nlvar=%eval(&nlvar+1); 
		 	%let lvars= &lvars &_first;
			%if &splot=Y %then %do;
				%if %datatyp(&_bas)=NUMERIC %then %do;
					%let phparl=&phparl &_bas;
					%let kmpar=&kmpar &_bas;
					%let kmvar=&kmvar &_first;
					%let _kmd=1;
				%end;
				%else %do;
					%let _bas2=0;
					%let phparl=&phparl &_bas2;
				%end;
			%end;
			%let exiterr_4=0;
		%end;
		%if %index(&_spec,S) %then %do;
			%let nsvar=%eval(&nsvar+1); 
			%let svars= &svars &_first;
			%if &splot=Y %then %do ;
				%if %datatyp(&_bas)=NUMERIC %then %do;
					%let phpars= &phpars &_bas;
					%if &_kmd=0 %then %do;
						%let kmpar=&kmpar &_bas;
						%let kmvar=&kmvar &_first;
					%end;
				%end;
				%else %do;
					%let _bas2=0;
					%let phpars= &phpars &_bas2;
				%end;
			%end;
			%let exiterr_4=0;
		%end;
		%let totvar = &totvar &_first;
	%end;
%end;

%if &exiterr_4=1 %then %goto final;

%if &splot=Y %then %do;
	%let phpar= &phparl &phpars;
	%let i=1;
	%let KM_var=;
	%let _v=;
	%do %while (%length(%scan(&kmvar,&i,%str( ))));
		%let _v= %scan(&kmvar,&i,%str( ))=%scan(&kmpar, &i, %str( ));
		%if &i=1 %then %let KM_var=  &_v;
		%if &i>1 %then %let KM_var= &KM_var and &_v;
 		%let i=%eval(&i+1);
	%end;
%end;

%let VAR=&lvars &svars;
%let nvar=%eval(&nlvar + &nsvar);
%let s=%eval(&nlvar+1);

%put KMvar = &Kmvar ;

/* Check for the existence and type of specified covariates in the date set */
%let totvar2= &totvar &TIME &CENSCOD &ID;
%if &BOOTSTRAP=Y %then %let totvar2=&totvar2 &strata; ;
%let dsid=%sysfunc(open(&data,i));
%let i=1;
%do %while (%length(%scan(&totvar2,&i,%str( ))));
  %let var_list=%scan(&totvar2,&i,%str( )) ;
  %let i2=0;
  %do %while(&i2=0);
    %do u=1 %to %sysfunc(attrn(&dsid,nvars));
  		%let var_dat= %upcase(%sysfunc(varname(&dsid,&u)));
 		%if %quote(&var_list) = %quote(&var_dat) %then %do ;
				%if ((%sysfunc(vartype(&dsid,&u)))=N or (%quote(&var_list)= &ID)
                 or (%quote(&var_list)= &strata)) %then %let i2=1;
				%else %let i2=2;
		%end;
 	%end;
	%if &i2=0 %then %let i2=2;
  %end;
  %if &i2=2 %then %do ; 
	%let exiterr_2=1;
	%let list_err= &list_err &var_list;
  %end;
  %let i=%eval(&i+1);
%end;

%let dsid=%sysfunc(close(&dsid));

%if &exiterr_2=1 %then %goto final;

%if &SURVPART=COX %then %do;
	/* Check for bootstap method */
	 %if &BOOTMET=ALL %then %let BOOTMET=PCTL BOOTN HYB JACK BC BCA;
	 %if &BOOTMET NE N %then %do;
		%let i=1;
		%do %while(%length(%scan(&BOOTMET,&i,%str( ) )));
 			%if %scan(&BOOTMET,&i,%str( ))=HYB %then %let HYB=Y;
 			%else %if %scan(&BOOTMET,&i,%str( ))=PCTL %then %let PCTL=Y;
  			%else %if %scan(&BOOTMET,&i,%str( ))=BOOTN %then %let BOOTN=Y;
 			%else %if %scan(&BOOTMET,&i,%str( ))=JACK %then %let JACK=Y;
 			%else %if %scan(&BOOTMET,&i,%str( ))=BC %then %let BC=Y;
 			%else %if %scan(&BOOTMET,&i,%str( ))=BCA %then %let BCA=Y;
   			%else %do ; 
				%let boot_err=&boot_err %scan(&BOOTMET,&i,%str( ));
				%let exiterr_3=1;
			%end;
			%let i=%eval(&i+1);
		%end;
	%end;
	%if &exiterr_3=1 %then %goto final;

	%if ((&BCA=Y or &JACK=Y) and &JACKDATA ne ) %then %do;
 		%if %sysfunc(exist(&JACKDATA))=0 %then %do;
			%let exiterr_1=2;
			%goto final;
		%end;
	%end;

	/** Covariates values set to O for the BASELINE statement (PROC PHREG)*/
	
	data _cov; 
	lpi=0;
	%if (&svars ne) %then %do;
 		%do i=1 %to &nsvar;
 			%scan(&svars,&i,%str(  ))=0;
 		%end;
	%end;
	run;

%end;

/* Missing covariates, number of observations and number of failed individuals */

%let _mis=;
proc sort data= &data; by descending &CENSCOD;
data sample0 (drop=_mis _mfail); set &data nobs=nobs end=last; by descending &CENSCOD;
%if &SURVPART=COX %then %do;
	pin=1; 
	pi=&CENSCOD;
	_obs_=_N_;
	_sample_=0;
%end;
if _n_=1 then do;_mis=0; _mfail=0; end;
retain _mis _mfail;
if nmiss(of &VAR &CENSCOD &TIME ) >0 then do;
	_mis=_mis+1; _mfail=(&CENSCOD=1)*_mfail+1;
	delete ;
end;
if last then do;
	if _mis>0 then do;
	call symput('_mis','('||trim(left(put(_mis,8.)))||' '||'observation(s) deleted due to missing information)');
	end;
end;
if &CENSCOD=1 and last.&CENSCOD then do;
call symput('nfail',trim(left(put(_N_-_mfail,8.))));
end;
call symput('nobs',trim(left(put(nobs-_mis,8.))));
run;


%final:;
	%let exiterr=%eval(&exiterr_0+&exiterr_1+&exiterr_2+&exiterr_3+&exiterr_4
                       +&exit_first+&exit_spec+&exit_bas);

	%if &exiterr=0 %then %do;
		%put ************************************************************************************;
		%put ;
		%put ***   Data set : &data;
        %put ***   Id idenfifier : &ID ** Censoring Status : &CENSCOD ** Failure/censoring time : &TIME;    
		%put ***   Number of observations : &NOBS &_mis;
        %put ***   Number of event times : &NFAIL ;    
		%put ***   Number of variables : &NVAR;                                                   
		%put ***   Incidence variables are : &nlvar : &LVARS  ;                               
		%put ***   Latency   variables are : &nsvar : &SVARS ;                                
        %put ;
		%put *** MODEL SPECIFICATIONS; 
		%put ***   Incidence part : link = &LINK;
		%put ***   Latency  part : &survpart2 &SURVPART3; 
	 %if &SURVPART=COX %then %do;
		%put ***   Baseline survival function estimation : &SU0MET; 
		%put ***   Tail completion : &TAIL; 
		%put ***   Number of iterations allowed : &MAXITER ;
		%put ***   Convergence criteria : &CONVCRIT ;
		%put;
		%put *** CONFIDENCE  INTERVAL COMPUTATION  ;
	    %if &jackdata ne %then %do;
		 %put ***   Jackknife estimates data set : &JACKDATA;
        %end; 
 		%put ***   Number of bootstrap samples : &NSAMPLE   ;
		%put ***   Jackknife method (JACK): &JACK;
		%put ***   Normal Theory Bias corrected CI (BOOTN) : &BOOTN    ; 
        %put ***   Hybrid  method (HYB): &HYB; 
		%put ***   Percentile Method (PCTL): &PCTL      ; 
        %put ***   Bias Corrected Bootstrap Method (BC) = &BC;  
		%put ***   Bias Corrected and Accelerated Method (BCA) : &BCA   ;                                   
     %end;
		%put;
        %put *** OTHER OUTPUTS;
	 %if &SURVPART=COX %then %do;
		%put ***   Graph of estimates : &GESTIMATE        ;  
	 %end; 
		%put ***   Baseline survival function : &BASELINE        ; 
		%put ***   Baseline Survival functions plots : &SPLOT ; 
		%if &splot=Y %then %do;
			%put ***   Value(s) for SPLOT : &KM_var;
		%end;
		%put ***   Plot fit and correlation statistics : &PLOTFIT; 
		%put; 
 		%put ************************************************************************************;
    %end;
	%else %do;
		%if &exiterr_1=1 %then %put ** ERROR ** Data set &data does not exist; 
		%if &exiterr_1=2 %then %put ** ERROR ** Date set &JACKDATA does not exist;
		%if &exiterr_2=1 %then %put ** ERROR ** Variable(s) &list_err not found in the data set &data or in character format;
		%if &exiterr_3=1 %then %put ** ERROR ** BOOTSTRAP method &boot_err is not available;
        %if &exiterr_4=1 %then %put ** ERROR ** Misspecification for Incidence (I) or Latency (S) ;

		%put ** the macro will stop. **; 
	%end;

%mend prepare;

/*************************************************************************/
/*                    PARAMETRIC  MIXTURE CURE MODELS                    */ 
/*************************************************************************/

%macro parametric;
/*initilisation of parameters estimates using classic logistic regression 
		and parametric survival analysis */

proc logistic data=&data noprint outest=lp (keep= &lvars intercept) ;
model &CENSCOD(event='1')= &lvars/link=&LINK;
run;

proc lifereg data=&data noprint outest=sp (keep=intercept &svars _SCALE_ );
model &TIME*&CENSCOD(0)= &svars/dist=&SURVPART;
where &CENSCOD=1;
run;

data lp; set lp;
%do i=1 %to &nlvar;
     rename %scan(&lvars,&i,%str( )) = L_%scan(&lvars,&i,%str( ));
%end;
rename intercept=L_int;
run;

data sp; set sp ;
%do i=1 %to &nsvar;
		%scan(&svars,&i,%str( ))=-(%scan(&svars,&i,%str( ))/_scale_);
         rename %scan(&svars,&i,%str( )) = S_%scan(&svars,&i,%str( ));
%end;
rename intercept=_scale;
%if &SURVPART=EXP %then %do;
	drop _scale_;
%end;
%else %do;
	rename _SCALE_=_shape;
%end;
run;

data param; merge lp sp;
run;

/* PROC NLMIXED */

ods output ParameterEstimates=estimates;
proc nlmixed data=&data ;
title "results for &data";
title2 "distribution &SURVPART, link=&LINK";
parms /data=param;

/* Incidence part of the model */
Leta=L_int;
%local i;
%do i=1 %to &nlvar;
	Leta&i= L_%scan(&lvars,&i,%str( ))* %scan(&lvars,&i,%str( ));
	Leta=Leta+Leta&i;
%end;
%if &LINK=LOGIT %then %do;
	expeta = exp(Leta);
	p_i= expeta/(1+expeta);
%end;
%if &LINK=CLOGLOG %then %do;
    expeta=exp(leta);
	p_i=1-exp(-expeta);
%end;
%if &LINK=PROBIT %then %do;
   p_i= cdf('Normal',leta);
%end;

/* Survival part of the model */

linpsurv=0;
%do n=1 %to &nsvar;
	Seta&n = S_%scan(&svars,&n,%str(  ))* %scan(&svars,&n,%str(  ));
	linpsurv = linpsurv + Seta&n;
%end;

%if &SURVPART=EXP %then %do;
	Su_t=exp(-exp(log(&TIME)-_scale+linpsurv));
	hu_t=exp(linpsurv-_scale);
	fu_t = exp(linpsurv-_scale)*Su_t;
%end;

/*%if &SURVPART=WEIBULL %then %do;
alpha= exp(linpsurv);
Su_t=exp(-alpha*(&TIME**_shape));
hu_t=_shape*&TIME**(_shape-1)*(alpha);
fu_t= hu_t*Su_t;
%end;*/

%if &SURVPART=WEIBULL %then %do;
	Su_t=exp(-exp((log(&TIME)-_scale+linpsurv*_shape)/_shape));
	hu_t=1/_shape*&TIME**(1/_shape-1)*exp((_shape*linpsurv-_scale)/_shape);
	fu_t= hu_t*Su_t;
%end;

%if  &SURVPART=LOGNORMAL %then %do;
	Su_t=1-probnorm((log(&TIME)-_scale+linpsurv*_shape)/_shape);
	fu_t=(1/(SQRT(2*constant('PI'))*_shape*&TIME))*exp(-0.5*((log(&TIME)-_scale+linpsurv*_shape)/_shape)**2);
	hu_t=fu_t/Su_t;
%end;

/*%if  &SURVPART=LOGNORMAL %then %do;
Su_t=1-probnorm((log(&TIME)-linpsurv)/_shape);
fu_t=(1/(SQRT(2*constant('PI'))*_shape*&TIME))*exp(-0.5*((log(&TIME)-linpsurv)/_shape)**2);
hu_t=fu_t/Su_t;
%end;*/

%if &SURVPART=LLOGISTIC %then %do;
	Su_t=(1+exp((log(&TIME)-_scale+linpsurv*_shape)/_shape))**(-1);
	fu_t= exp((log(&TIME)-_scale+linpsurv*_shape)/_shape)/(_shape*(1+exp((log(&TIME)-_scale+linpsurv*_shape)/_shape))**2);
	hu_t=fu_t/Su_t;
%end;

llkl = (&CENSCOD=1)*(log(p_i)+ log(fu_t)) + (&CENSCOD=0)*(log(1-p_i+p_i*Su_t));

model &TIME ~ general(llkl);

run;

%if &LINK=LOGIT %then %do;

	data OR (keep= parameter OR upper lower probt) 
  	     RR (keep= parameter RR upper lower probt); set estimates;
		if substr(parameter,1,2)='L_' and parameter ne 'L_int' then do;
			OR = exp(estimate);
			Upper=exp(upper);
			lower=exp(lower);
			output OR;
		end;
		if substr(parameter,1,2)='S_' and parameter ne 'S_int' then do;
			RR = exp(estimate);
			Upper=exp(upper);
			lower=exp(lower);
			output RR;
		end;
	run;

	proc print data=OR noobs l ;
	var parameter OR lower upper probt;
	title "Odd's ratio estimates for data=&data";
	title2 "distribution &SURVPART, alpha=&ALPHA";
	proc print data=RR noobs l ;
	var parameter RR lower upper probt;
	title "Relative risk estimates for data=&data";
	title2 "distribution &SURVPART,alpha=&ALPHA ";
	run;

%end;


%mend parametric;

/*************************************************************************/
/*                           BOOTSTRAP                                   */ 
/*************************************************************************/

%macro bootstrap ;

%put ** SAS is resampling the original sample : &run over &nr time(s), please wait...;

proc multtest data=sample0 outsamp=sample (keep=_sample_ _obs_) 
noprint bootstrap  nsample=&NSAMPLE_b seed=&seed; 
class _obs_;
test mean(&CENSCOD);
%if (&strata ne) %then %do;
	strata &strata;
%end;
run;

proc sort data=sample; by _obs_;
proc sort data=sample0; by _obs_;

data sample ; merge sample (in=k) sample0 (in=m drop=_sample_); by _obs_; if k;
_sample_=&nb*(&run-1)+_sample_;
run;


%if &run=1 %then %do;
	proc append base=sample data=sample0;
	run;
%end;

%mend bootstrap;


/*************************************************************************/
/*               GET THE MAX FAILURE TIME FOR EACH REPLICATE            */ 
/*************************************************************************/

%macro max;

proc sort data= sample; by _sample_ &time descending &CENSCOD;
run;

data min_max (keep=_sample_ min max); set sample (keep=_sample_ &CENSCOD &TIME) ; by _sample_ descending &CENSCOD;
where &CENSCOD=1;
retain min max;
if first.&CENSCOD then min=&TIME;
if last.&CENSCOD then max=&TIME;
if last._sample_ then output;
run;

data sample; merge sample min_max; by _sample_;
run;


%mend max;

/*************************************************************************/
/*                   %COMPARE                                            */
/*    Tests for convergence in parameters estimates and likelihood       */
/*************************************************************************/

%macro compare;

%let iv=%eval(&it-1);
data compare  ; merge comp&iv (in=i rename=(est=old_est)) comp&it (in=j) end=last;
by _sample_; if j;
	retain crit_u crit_g conv_1;
	if first._sample_ then do; crit_u=0; end;
	diff=abs(old_est)-abs(est);
	denom=(abs(old_est)+abs(est))/2;
	if (denom>&CONVCRIT) then do;
		reldiff=abs(old_est-est)/denom;
		crit_u=max(crit_u,reldiff);
		crit_g=max(crit_g,reldiff);
	end;
	if last._sample_ then do;
		if crit_u<&CONVCRIT then do; 
			conv_u=1; conv_1=1; 
		end; 
		else do; 
			conv_u=0; 
		end; 
		output;
    end;
    if last then do;
		call symput('crit_g',trim(left(crit_g)));
	end;
	run;

	/*** Convergence indicators ***/

	data _null_; set compare end=last;
	if last then do;
		call symput('nb_nc',trim(left(_n_)));
		call symput('conv_1',trim(left(conv_1)));
	end;
	run;

	data _null_; 
		crit_g=&crit_g;
		if (crit_g<&CONVCRIT) then conv_g=1;
		else conv_g=0;
		call symput('conv_g', left(conv_g));
	run;

	/*** Convergence indicator for initial sample and selection of non converged samples ***/
     
        proc sort data=compare ; by descending conv_u  _sample_;
    	data
    	%if (&conv_g ne 1) %then n_conv (keep=_sample_) ; 
        %if (&conv_1=1) %then  conv (keep=_sample_ crit_u it);;
        set compare(keep=_sample_ conv_u crit_u) ; by descending conv_u;
		%if (&conv_0 ne 1) %then %do;
			if _sample_=0 then do;
				if (conv_u=1 or &it>=&MAXITER) then do; 
                	conv_0=1;
					call symput('crit_0',trim(left(crit_u)));
				end; 
				else conv_0=0;
				call symput('conv_0',trim(left(conv_0)));
			end;
		%end;
		%if (&conv_1=1) %then %do;
			if conv_u=1 then do; 
				it=&it;
				if last.conv_u then do;
					nb_c=_n_;
					call symput('nb_c',trim(left(nb_c))); 
				end;
				output conv;
			end;
		%end;
		%if (&conv_g ne 1) %then %do;
		if conv_u=0 then do;
			output n_conv;
		end;
		%end;
   		run;

	proc datasets lib=work nolist;
	delete comp&iv;
	run;


%mend compare;


/*************************************************************************/
/*                          WTAIL tail completion method                 */ 
/*************************************************************************/

%macro wtail;


	proc nlmixed data=sortie1;
	by _sample_;
	%if &it>1 %then %do;
		parms /data=_weib;
	%end;

	lambda=((-log(Su0_der))**(1/rho))/max;
	hu0=lambda**rho * rho *&TIME**(rho-1);
	Su0=exp(-(lambda*&TIME)**rho);

	llkl = (&CENSCOD=1)*(log(hu0)+est)+(pi*exp(est))*log(Su0);

	model &TIME ~ general(llkl);
	ods output ParameterEstimates=weib(keep= _sample_ parameter estimate);
	ods listing close;
	run;

	data _weib (drop=_sample_) weib1 (drop=parameter rename=(estimate=rho)) ; set weib;
	if _sample_=0 then output _weib weib1;
	else output weib1;
	run;

	data sortie1; merge sortie1 weib1; by _sample_; 
	lambda=((-log(Su0_der))**(1/rho))/max;
	if &TIME>=max then do;
	Su0=exp(-(lambda*&TIME)**(rho));
	pi=&CENSCOD+(1-&CENSCOD)*(IP_1*(su0**exp(est)))/((1-IP_1)+(IP_1*(su0**exp(est))));
	end;
	run;


%mend wtail;

/*************************************************************************/
/*                      EM ALGO : MAXIMISATION PART                      */ 
/*************************************************************************/

%macro maximisation;

	/** Incidence part **/

	proc logistic data=sample  noprint outest=lp(drop= _TYPE_ _LINK_ _STATUS_ _NAME_ )
    %if (&it>1) %then inest=lp;;
	model pi/pin=&lvars/ link=&LINK;
	output out=l prob=ip_1;
	by _sample_;
	run;

	** Rename incidence estimates and detection of sample(s) with quasi complete separation ;

    data lp2 (keep=var0-var&nlvar _sample_ _lnlike_ rename=(_lnlike_=lkl)) 
	%if &it=1 %then lp(keep=_sample_ intercept &lvars) sep(keep=_sample_);;
    set lp ;
    by _sample_;
		array old_var(&s) intercept &lvars;
		array new_var(&s) var0-var&nlvar;
		retain var0-var&nlvar ;
		if first._sample_ then do; i=1;end;
		do i=1 to &s;
			new_var(i)=old_var(i);
			if i>1 then do;
				if abs(new_var(i))>6 then sep=1;
			end;
		end;
		if last._sample_ then do; 
			if sep ne 1 then output lp2;	
			%if &it=1 %then %do;
				if sep ne 1 then output lp;
				if sep=1 then  output sep;
			%end;
		end;
	run;


%if &it=1 %then %do;
	data _null_; set sep nobs=nobs;
	call symput('nb_sep', trim(left(put(nobs,best8.))));
	run;
%put Number of sample(s) with quasi complete separation : &nb_sep ;
%end;

/** Survival part : cox PH ***/

** Sélection of uncured indivuduals i.e pi ne 0;

%if &it=1 %then %do;
	%if &nb_sep>0 %then %do;
		data surv l; merge l(in=i) lp(in=j keep=_sample_); by _sample_; if j; 
	%end;
	%else %do;
		data surv ; set l;
	%end;
 	%if ((&TAIL=ETAIL) or (&TAIL=WTAIL)) %then %do;
		if pi=0 then pi=1e-10;
	 %end;
	if pi ne 0 then do; lpi=log(pi) ; output surv; end;
	%if &nb_sep>0 %then %do;
	output l;
	%end;
	run;
%end;

%else %do;
	data surv; set l; 
	if pi ne 0;
	lpi=log(pi);
%end;


** Cox PH;

	proc phreg data=surv noprint outest=sp (drop= _TYPE_ _TIES_ _STATUS_ _NAME_  lpi) ;
	model &TIME*&CENSCOD(0)=&svars/offset=lpi;
	baseline out=bas (keep=_sample_ &TIME Su0 ) covariates=_cov  survival=su0/nomean method=&SU0MET;
	by _sample_;
	run;

** Rename survival estimates;

	%if (&svars ne )%then data sp2 (keep=var&s-var&nvar _sample_ _lnlike_ rename=(_lnlike_=lks));
    %if (&svars =) %then  data sp2 (keep= _sample_ _lnlike_ rename=(_lnlike_=lks));;
    set sp ; by _sample_;
	%if (&svars ne ) %then %do;
		array old_var(&nsvar) &svars ;
		array new_var(&nsvar) var&s-var&nvar;
		if first._sample_ then do; i=1;end;
			do i=1 to &nsvar;
			new_var(i)=old_var(i);
		end;
	%end;
	run;

/*** Transpose and rename Parapeters estimates  ***/

	data comp&it ;  merge lp2 sp2 ; by _sample_; 
	run;
	proc transpose data=comp&it out=comp&it(drop=_label_ rename=(col1=Est)) name=Variable;
	by _sample_;
	run;

%mend maximisation;

/*************************************************************************/
/*                           EM ALGO : EXPECTATION  PART                 */ 
/*************************************************************************/
%macro expectation;

/** Constraint : if time<min then  Su0=1
                 if time>maxi then Su0=0 (ZERO) SU0~weibull  (Wtail) Su0~exponentiel (ETAIL)(Su0(tk)->0)**/

data sample; merge l bas ; by _sample_ &TIME;
if &CENSCOD ne .;
run;

%if (&svars ne ) %then %do;
	data sample; merge sample sp2 (drop=lks); by _sample_;
	run;
%end;

%if ((&TAIL=ETAIL) or (&TAIL=WTAIL) or (&TAIL=NONE)) %then %do; 
	data Su0_der (keep= _sample_ Su0_der) ; set sample (keep= _sample_ &TIME max Su0); by _sample_;
	where Su0 ne . and 0<(1-&TIME/max)<0.5;
	if last._sample_ then do ; Su0_der=su0; output; end;
	run;
	data sample; merge sample Su0_der; by _sample_;
	run;
%end;

data sortie1 ; set sample end=last; by  _sample_ &TIME;
if &TIME<min then su0=1;
retain survi;
if su0 ne . then do; survi=su0; end;
su0=survi;
if &TIME>=max then do;
    %if (&TAIL=NONE) %then %do;
		Su0=Su0_der;
	%end;
	%if (&TAIL=ZERO) %then %do;
		Su0=0;
	%end;
	%if (&TAIL=ETAIL) %then %do;
   	  lambda=-log(Su0_der)/max;
	  Su0=exp(-lambda*&TIME);
	%end;
	%if ((&TAIL=WTAIL) and (&it=1)) %then %do;
      	 Su0=Su0_der;
	%end;
end;
est=0;
%if (&svars ne ) %then %do;
	array surv_var(&nsvar) &svars;
	array surv_par(&nsvar) var&s-var&nvar;
	do i=1 to &nsvar;
		est= est+ surv_var(i)*surv_par(i);
	end;
	drop i;
%end;
pi=&CENSCOD+(1-&CENSCOD)*(IP_1*(su0**exp(est)))/((1-IP_1)+(IP_1*(su0**exp(est))));
run;


%if (&TAIL=WTAIL) %then %do;
%wtail;
%end;

%mend expectation;


/*************************************************************************/
/*     END OF ITERATION SELECTION OF NON CONVERGED SAMPLES               */ 
/*************************************************************************/

%macro end_it;
%if (&it=1 or &conv_g=1) %then %do;
	data sample; set sortie1 (keep= _sample_ &ID &VAR &CENSCOD &TIME pi pin Su0 min max);
%end;
%else %if (&conv_g ne 1) %then %do;
	data sample; merge sortie1 (in=n keep= _sample_ &ID &VAR &CENSCOD &TIME pi pin Su0 min max) 
                             n_conv (in=m); by _sample_; if m;
%end;
	%if ((&TAIL=ZERO) or (&TAIL=NONE)) %then %do;
	Su0=.;
	%end;
	%if ((&TAIL=ETAIL or &TAIL=WTAIL)) %then %do;
	if &TIME<max then Su0=.;
	%end;
	run;

%if &it>1 %then %do;
	%if &conv_g ne 1 %then %do;
		%if &BOOTSTRAP=Y %then %do;
	       %put  * &nb_c replicate(s)(over &nb_nc left) converged in &it iterations (max. rel.diff. =&crit_g) ;
		   %if &conv_0=1 and &done_0 ne 1 %then %put * NOTE: Convergence attained for initial sample ;
		%end;
		%else %do;
		   %put * &msg1 : iteration &it : max rel. diff. in estimates =&crit_g **;
		%end;
	%end;
%end;

%if &conv_g=1 %then %do;
	%put  ;
	%put    *** Run &run : &msg1 have converged after &it iterations, M.R.D=&crit_g **; 
	%put  ;
%end;

%if &it>=&MAXITER %then %do;
	%if &BOOTSTRAP=Y %then %do;
	   %put;	
	   %put ********************************************************************************;
	   %put *** WARNING : &nb_nc replicate(s) did not converge after &MAXITER itérations ***;
	   %put ********************************************************************************;
	%end;
	%else %do;
	   %put ;
	   %put ****************************************************************;
	   %put *** WARNING : No convergence attained after &MAXITER itérations *;
	   %put *****************************************************************;
	%end;
%end;


%mend end_it;

/*************************************************************************/
/*                        RESULTS FOR INITIAL SAMPLE                     */ 
/*************************************************************************/
%macro result_0;

data estimates (rename=(est=estimate)); set comp&it;
if _sample_=0;
if variable='lkl' then do; 
	call symput('lkl',trim(left(est)));
	delete;
end;
if variable='lks' then do; 
	call symput('lks',trim(left(est)));
	delete;
end;
run;

data probcure ; 
attrib _sample_ &ID  &CENSCOD &TIME &VAR label='' pi label="Estimated prob. of being uncured" 
cured label="Estimated prob. of being cured" Su0 label="Conditionnal Baseline survival function" var0-var&nlvar label='' %if &svars ne %then var&s-var&nvar label='';;
merge lp2 (keep=_sample_ var0-var&nlvar)
 %if (&svars ne ) %then sortie1 (keep= _sample_ &ID &TIME &CENSCOD &VAR pi Su0 var&s-var&nvar);
 %if (&svars=)    %then sortie1 (keep= _sample_ &ID &TIME &CENSCOD &VAR pi Su0 );;
by _sample_;
where _sample_=0;
rename pi=uncured;
cured=1-pi;
keep &ID &VAR &TIME &CENSCOD var0-var&nvar pi cured Su0;
run;



data likelihood;
	N_Obs = &nobs;
	N_failed = &nfail;
	N_iter = &it;
	format Converg e12.;
	Converg=&crit_0;
	LogL_&LINK = &Lkl;
	LogL_coxPH = &lKs;
	LogTotal=(&lkl+&lks);
	label n_obs='Number of individuals' 
	      n_failed ='Number of individuals that failed'
		  n_iter='Number of iteration until convergence'
		  converg='Convergence indicator'
		  LogL_&LINK="loglikelihood for the &LINK part"	 
          logL_coxPH='loglikelihood for the survival (cox PH) part'
		  logTotal='log likelihood for the complete data';
	%if &BOOTSTRAP=Y %then %do;
		Bootstrap= &NSAMPLE;
		label boostrap='Number of boostrap re_samples';
	%end;
	run;

	Proc print data=likelihood noobs;
	Title1 "RESULTS FOR DATA &DATA";
	Title2 "Convergence and Log_likelihood";
    run;

%if &fast=Y %then %do;
   
	ods listing close;
	ods output ParameterEstimates=Fast_Inci ;
	ods output oddsratios=Fast_OR  ;

	proc logistic data=l ;
    model pi/pin=&lvars/ link=&LINK alpha=&ALPHA;
	where _sample_=0;
	/*title  "FAST RESULTS FOR THE &LINK PART (&data)";
	title2 "tail completion method : &TAIL";
	title3 "Su0 estimation : &SU0MET";*/
	run;

	data b; set l; if pi ne 0 and _sample_= 0; lpi=log(pi);

	proc phreg data=b  ; 
	%if (&svars ne ) %then ods output parameterestimates=Fast_Surv;;
	model &TIME*&CENSCOD(0)=&svars/offset=lpi rl alpha=&ALPHA;
	/*title  "FAST RESULTS FOR THE SURVIVAL PART (&data)";
	title2 "tail completion method : &TAIL";
	title3 "Su0 estimation : &SU0MET";*/
	run;  
 
	ods listing;

	proc print data=Fast_inci l noobs;
	title  "FAST RESULTS FOR THE &LINK PART (&data)";
	title2 'Analysis of Maximum Likelihood Estimates';
	run;

    proc print data=Fast_inci l noobs;
	title 'Odds Ratio Estimates';
	run;

	%if (&svars ne ) %then %do;
		proc print data=Fast_surv l noobs ;
		title  "FAST RESULTS FOR THE SURVIVAL PART (&data)";
		title2 'Analysis of Maximum Likelihood Estimates';
		run;
	%end;

%end;

%mend result_0;


/*************************************************************************/
/*                        GLOBAL RESULTS FOR BOOSTRAP                    */ 
/*************************************************************************/

%macro result;

data _boot_; merge lp2 (in=i drop=lkl)  conv (in=j) %if (&svars ne ) %then sp2 (in=i drop=lks);; 
by _sample_; if j;
run;

proc append base=&boot_res data=_boot_;
run;

%mend result;


/*************************************************************************/
/*              BASELINE SURVIVAL FONCTION                               */ 
/*************************************************************************/

%macro baseline;

%if &SURVPART=COX %then %do;
	%if ((&conv_0=1 or &it>&MAXITER) and &done_0 ne 1) %then %do;
		data _baseline ;
        attrib _sample_ &lvars %if &svars ne %then &svars; var0-var&nlvar %if &svars ne %then var&s-var&nvar;
        &TIME label='' pi label="prob_&CENSCOD=1" Su0 label="Conditionnal Baseline survival function" ;;

        merge lp2 (keep=_sample_ var0-var&nlvar)
	    %if (&svars ne ) %then sortie1 (keep= _sample_ &TIME &lvars  pi Su0 &svars var&s-var&nvar);
		%if (&svars=) %then sortie1 (keep= _sample_ &TIME &lvars pi Su0 );;
        by _sample_;
		where _sample_=0;
		it=&it; 
		crit_u=&crit_0;
	    run;

		%if &BOOTSTRAP=Y %then %do;
			proc append base=_baseline_t data=_baseline;
		%end;
	%end;

	%if (&BOOTSTRAP=Y and (&conv_1=1 or &it>=MAXITER)) %then %do;
		data bas_boot ; 
        attrib _sample_ &lvars %if &svars ne %then &svars; var0-var&nlvar %if &svars ne %then var&s-var&nvar;
        &TIME label='' pi label="prob_&CENSCOD=1" Su0 label="Conditionnal Baseline survival function" ;;
		merge conv(in=j) lp2 (in=i keep=_sample_ var0-var&nlvar)
	    %if (&svars ne ) %then sortie1 (in=i keep= _sample_ &TIME &lvars  pi Su0 &svars var&s-var&nvar);
		%if (svars =) %then sortie1 (in=i keep= _sample_ &TIME &lvars  pi Su0 );;
        by _sample_; if j;
		*label pi="prob_&CENSCOD=1";
	    run;

	    proc  append base=_baseline_t data=bas_boot ;
	    run;
	%end;
%end;


%if &SURVPART2=PARA %then %do;
	%let par=%eval(&s+1);
	%let par2=%eval(&nvar+3);
	%let par3=%eval(&s+2);
	%let par4=%eval(&nvar+2);
	
	proc sort data=&data; by &TIME;
	data _null_; set &data end=last;
	if _n_=1 then do;
		call symput('min',trim(left(&TIME)));
	end;
	if last then do;
		call symput('max',trim(left(&TIME)));
	end;
	run;

	proc transpose data=estimates  out=graph name=parameter;
	var estimate ;
	run;

	data _baseline; set graph ;
	rename col1=var0 col&par=_scale %if &SURVPART=EXP %then; 
           %if &SURVPART ne EXP %then col&par2=_shape;;
	label col1=L_int;
	%do i=2 %to &s;
		 %let i2=%eval(&i-1);
		 rename col&i=var&i2;*L_%scan(&lvars,&i2,%str( ));
		 label col&i="L_%scan(&lvars,&i2,%str( ))";
	%end;
	%do i=&par3 %to &par4;
		%let i2=%eval(&i-&par);
		%let i3=%eval(&i-2);
		 rename col&i=var&i3;* S_%scan(&svars,&i2,%str( ));
		 label col&i="S_%scan(&svars,&i2,%str( ))";
	%end;
	label parameter='';
	unit=(&max-&min)/50;
	min=&min;
	max=&max;
	do t=0.01 to max by unit;
	%if &SURVPART=WEIBULL %then %do;
		*Su0=exp(-exp(col&par)*(t**col&par2)); 
		Su0=exp(-exp((log(t)-col&par)/col&par2)); 
	%end;
	%if &SURVPART=LLOGISTIC %then %do;
		Su0=(1+exp((log(t)-col&par)/col&par2))**(-1); 
	%end;
	%if &SURVPART=LOGNORMAL %then %do;
		Su0=1-probnorm((log(t)-col&par)/col&par2);
	%end;
	%if &SURVPART=EXP %then %do;
		Su0 = exp(-exp(log(t)-col&par));
	%end;
	%if &LINK=LOGIT %then %do;
		p=exp(col1)/(1+exp(col1));
	%end;
	%if &LINK=PROBIT %then %do;
		p=cdf('Normal',col1);
	%end;
	%if &LINK=CLOGLOG %then %do;
		p=1-exp(-exp(col1));
	%end;
		S0=(1-p)+p*Su0;
		output;
	end;
	drop unit min max ;*col1 col&par col&par2 ;	
	run;
%end;

%mend baseline;


/*************************************************************************/
/*              SURVIVAL FUCNTIONS  PLOTS                                */ 
/*************************************************************************/

%macro survplot;

proc lifetest data=&data outsurv=_graphKMu (rename=(survival=KMu)keep= &TIME survival) noprint ;
time &TIME*&CENSCOD(0);
where &CENSCOD=1 and &KM_var; 
run;

proc lifetest data=&data outsurv=_graphKM (rename=(survival=KMp) keep=&TIME survival) noprint;
time &TIME*&CENSCOD(0);
where &KM_var;
run;

%if &syserr>4 %then %let nokm=1 ;

%if &nokm=1 %then %do;
%put;
%put *** SURVPLOT : ;
%put *** WARNING : The Kaplan-Meier estimates for &KM_VAR will not be displayed ****;
%put *** One of the covariates may be a continous variable or its value is not observed ****;
%put ;
%end;


%if &SURVPART=COX %then %do;
   	proc sort data=_baseline out=graph (keep=var0-var&nvar &TIME Su0); by &TIME;
        		
	data graph; 
	%if &nokm=0 %then  %do;
		merge graph _graphKMu _graphKM; by &TIME;
	%end;
	%if &nokm=1 %then %do;
   		set graph;
	%end;
		%let i=1;
		%do %while (%length(%scan(&VAR, &i, %str( ))));
			basva&i=%scan(&phpar,&i,%str( ));
			label basva&i="%scan(&VAR,&i,%str( ))";
 			%let i=%eval(&i+1);
		%end;
		%if &nokm=0 %then %do;
			retain sp su;
			if KMp ne . then do; sp=Kmp; end;
			Kmp=sp;
			if KMu ne . then do; su=Kmu; end;
			kmu=su;
			label Kmu="Product Limit estimate :uncensored"
	      	  	Kmp="Product limit estimate : population";
			drop su sp;
		%end;
		est_l=var0;
		est_s=0;
		%if &lvars ne %then %do;
			array basval(&nlvar) basva1-basva&nlvar;
			array inc_par(&nlvar) var1-var&nlvar;
			do n=1 to &nlvar;
				est_l=est_l+basval(n)*inc_par(n);
			end;
			drop n;
		%end;
		%if &svars ne %then %do; 
	    	array basvas(&nsvar) basva&s-basva&nvar;
	   		array surv_par(&nsvar) var&s-var&nvar;
       		 do m=1 to &nsvar;
  				est_s=est_s+basvas(m)*surv_par(m);
			end;
			drop m;
		%end;
		p=exp(est_l)/(1+exp(est_l));
		S_u=Su0**exp(est_s);
    	S_p=(1-p)+p*S_u;
		if _n_=1 then do;S_u=1;S_p=1; end;
		label S_u="PH cure Model : cured fraction"
			  
		  	  S_p="PH cure model : population";
		drop  var0-var&nvar  ;
		run;
	%end;

	
	%if &SURVPART2=PARA %then %do;
	     
		data graph; 
		%if &nokm=0 %then %do; 
			merge _baseline (rename=(t=&TIME)) _graphKMu _graphkm; by &TIME;
		%end;
		%if &nokm=1 %then %do;
			set _baseline (rename=(t=&TIME)); 
		%end;
		%let i=1;
		%do %while (%length(%scan(&VAR, &i, %str( ))));
			basva&i=%scan(&phpar,&i,%str( ));
			label basva&i="%scan(&VAR,&i,%str( ))";
 			%let i=%eval(&i+1);
		%end;
		%if &nokm=0 %then %do;
			retain sp su;
			if KMp ne . then do; sp=Kmp; end;
			Kmp=sp;
			if KMu ne . then do; su=Kmu; end;
			kmu=su;
			Label Kmu="Product Limit estimate :uncensored"
	     	 	  Kmp="Product limit estimate : population";
			drop su sp;
		%end;
		est_l=var0;
		est_s=0;
		%if &lvars ne %then %do;
			array basval(&nlvar) basva1-basva&nlvar;
			array inc_par(&nlvar) var1-var&nlvar;
			do n=1 to &nlvar;
				est_l=est_l+basval(n)*inc_par(n);
			end;
			drop n;
		%end;
		%if &svars ne %then %do; 
	    	array basvas(&nsvar) basva&s-basva&nvar;
	   		array surv_par(&nsvar) var&s-var&nvar;
       		 do m=1 to &nsvar;
  				est_s=est_s+basvas(m)*surv_par(m);
			end;
			drop m;
		%end;
		p=exp(est_l)/(1+exp(est_l));
		S_u=Su0**exp(est_s);
    	S_p=(1-p)+p*S_u;
		if _n_=1 then do;S_u=1;S_p=1; end;
		label S_u="&SURVPART : cured fraction"		  
		  	  S_p="&SURVPART : population";
        drop var0-var&nvar;
	    run;

	%end;


	goptions reset=global gunit=pct border cback=white
        	 ftitle=swissb ftext=swiss htitle=3 htext=2;

	symbol1 c=red w=1 i=join;
	symbol2 c=black  i=steplj;
	axis1 label=(h=2 angle=90 rotate=0 "Conditional Survival Function Estimate" )
	  	  order=(0 to 1 by 0.1)
   	      major=(height=1 width=1)
          minor=(number=4 color=black height=0.5 width=1);

    axis3 label=(h=2 angle=90 rotate=0 "Marginal Survival Function Estimate" )
     	  order=(0 to 1 by 0.1)
	 	  major=(height=1 width=1)
       	  minor=(number=4 color=black height=0.5 width=1);

    legend1 position=(bottom center outside) 
            value=(tick=1 "&LINK &SURVPART cure Model" tick=2 "Product Limit estimate ");
    legend2 position=(bottom center outside) 
            value=(tick=1 "&LINK &SURVPART cure Model" tick=2 "Product Limit estimate ");
       
	proc gplot data=graph; 
	plot S_u*&TIME=1 
	     %if &noKM=0 %then 
	 	 KMu*&TIME=2; /overlay legend=legend1 vaxis=axis1 cframe=ligr;;
		 title " Conditionnal Survival Function estimate for &data";
		 title2 "(uncured fraction : &KM_var)";
		 %if &SURVPART=COX %then %do;
		 title3 "Tail completion method : &TAIL Su0 estimation : &SU0MET";
		 %end;
		 %if &SURVPART2=PARA %then %do;
		 title3 "Survival Distribution function=&SURVPART";
		 %end;;
	run;
	quit;
	proc gplot data=graph; 
	plot S_p*&TIME=1
		%if &nokm=0 %then
	 	 KMp*&TIME=2; /overlay legend=legend2 vaxis=axis3 cframe=ligr;;
	     title " Marginal Survival Function estimate for &data";
	     title2 "(overall population : &KM_var)";
		 %if &SURVPART=COX %then %do;
		 title3 "Tail completion method : &TAIL Su0 estimation : &SU0MET";
		 %end;
		 %if &SURVPART2=PARA %then %do;
		 title3 "Survival Distribution function=&SURVPART";
		 %end;;
	run;
	quit;

%mend survplot;


/*************************************************************************/
/*                      PLOT FIT                                         */
/*************************************************************************/

%macro plotfit;

%if &nokm=1 %then %do ;
	%goto endp ;
%end;
 
/*  KME estimation for each group */

proc lifetest data=&data noprint outsurv=plotfit; 
time &TIME*&CENSCOD(0);
strata &KMvar ;
run;
data _leg; set plotfit; by stratum;
if first.stratum;
rename stratum=_stratum_;
keep stratum &KMvar;
run;

%if &syserr>4 %then %do ;
	%let &nokm =1;
    %goto endp;
%end;


data plotfit (drop=sp sdf_ucl sdf_lcl _censor_ rename=(stratum=_stratum_)) ;
set plotfit  end=last; by stratum ;
    retain sp ;
	if survival ne . then sp=survival;
    survival=sp;
	if last then do ;
		call symput('nstratum',trim(left(put(stratum,8.))));
	end;
run;

/* 2 So(t) Estimation for cure model */

%if &SURVPART=COX %then %do;
	proc sort data=_baseline out=cm (keep= var0-var&nvar su0 &TIME); by descending &TIME;
	run;
	proc sort data=plotfit; by descending &TIME;
	run;
%end;

%if &SURVPART2=PARA %then %do;
	data cm; set _baseline;
	keep var0-var&nvar _scale %if &SURVPART ne EXP %then _shape ;;
%end;

	data plotfit; merge plotfit cm;
	%if &SURVPART=COX %then by descending &TIME;;
	retain _v0-_v&nvar ;
    %if &SURVPART=COX %then %do;
    	retain _su0;
	%end;
	%if &SURVPART2=PARA %then %do; 
		retain _scal %if &SURVPART ne EXP %then _shap ;;
		if _scale ne . then _scal=_scale; _scale=_scal;
		if _shape ne . then _shap=_shape; _shape=_shap;
	%end;
	%if &SURVPART=COX %then %do;
		if su0 ne . then _su0=su0; su0=_su0;
	%end;
	%do nv=0 %to &nvar;
		if var&nv ne . then _v&nv=var&nv;
		var&nv=_v&nv;
	%end;
	est_l=var0;
	est_s=0;

	%if &lvars ne %then %do;
		array basval(&nlvar) &lvars ;
		array inc_par(&nlvar) var1-var&nlvar;
		do n=1 to &nlvar;
		    xb = basval(n)*inc_par(n);
			if xb ne . then do;
		 		est_l=est_l+basval(n)*inc_par(n);
			end;
		end;
		drop n;
	%end;
	%if &svars ne %then %do; 
	   	array basvas(&nsvar) &svars ;
	 	array surv_par(&nsvar) var&s-var&nvar;
     	do m=1 to &nsvar;
			Zg = basvas(m)*surv_par(m);
			if zg ne . then do;
  				est_s=est_s+basvas(m)*surv_par(m);
			end;
		end;
		drop m;
	%end;
	%if &SURVPART2=PARA %then %do;
		if &TIME ne 0 then do;
			%if &SURVPART=WEIBULL %then %do;
				*Su0=exp(-exp(_scale)*(&TIME**_shape)); 
				Su0=exp(-exp((log(&TIME)-_scale)/_shape)); 
			%end;
			%if &SURVPART=LLOGISTIC %then %do;
				Su0=(1+exp((log(&TIME)-_scale)/_shape))**(-1); 
			%end;
			%if &SURVPART=LOGNORMAL %then %do;
				Su0=1-probnorm((log(&TIME)-_scale)/_shape);
			%end;
			%if &SURVPART=EXP %then %do;
				Su0 = exp(-exp(log(&TIME)-_scale));
			%end;
	   end;
	%end;
	%if &LINK=LOGIT %then %do;
			p=exp(est_l)/(1+exp(est_l));
	%end;
	%if &LINK=PROBIT %then %do;
			p=cdf('Normal',est_l);
	%end;
	%if &LINK=CLOGLOG %then %do;
			p=1-exp(-exp(est_l));
	%end;	
	S_u=Su0**exp(est_s);
   	S_p=(1-p)+p*S_u;
	if &TIME=0 then do;
	S_u=1; S_p=1; Su0=1;
	end;
	drop _v0-_v&nvar ;
	%if &SURVPART=COX %then %do;
    drop _su0;
	%end;
    %if &SURVPART=PARA %then %do;
	drop _scal %if &SURVPART2 ne EXP %then _shap;;
	%end;
	label p  ="Proportion of uncured individuals"
		  s_p="Marginal Survival Distribution function"
	      s_u="Conditionnal Survival Distribution function";
	drop est_l est_s %if &lvars ne %then xb; %if &svars ne %then zg  ;;
	run;

	proc sort data=plotfit ;  
	by _stratum_ &TIME;
	run;

/* 3 plot KME(Xi,Yi)_pop and estimated S(t,Xi,Yi)_pop */

	goptions reset=all;
	%if &nstratum<5 %then %do;
		%let colour= red black blue green ;
		%do nst2=1 %to &nstratum;
			%let cs=%scan(&colour,&nst2, %str( ));
			symbol&nst2 i=join c=&cs ;
		%end;
		%let nst4=%eval(&nstratum+&nstratum);
		%do nst3=%eval(&nstratum+1) %to &nst4;
			%let cs2=%scan(&colour,%eval(&nst3-&nstratum), %str( ));
			symbol&nst3 i=steplj c=&cs2 l=4 w=1.8 ;
		%end;
	%end;
	%if &nstratum>4 %then %do;
		symbol1 c=red i=join;
		symbol2 c=black i=steplj l=4 w=1.8;
	%end;

	legend1 label=(justify=center "   Cure Model estimate ")
    shape=symbol(4,2)
	origin=(,0.2 cm);

	legend2 label=(justify=center " Kaplan-Meier estimate " )
        shape=symbol(4,2)
		origin=(,0.5 cm);

	axis1 label=(h=1 angle=90 rotate=0 " Survival Function Estimate" )
	  	  order=(0 to 1 by 0.1)
   	      major=(height=1 width=1)
          minor=(number=2 color=black height=0.5 width=1);
	axis2 label=none
	      order=(0 to 1 by 0.1)
          major=(height=1 width=1)
          minor=(number=2 color=black height=0.5 width=1);


	proc gplot data=plotfit;
	plot  s_p*&TIME=_stratum_ / vaxis=axis1 legend=legend1 cframe=ligr;
	plot2 survival*&TIME=_stratum_/vaxis=axis2 legend=legend2 cframe=ligr;
	%if &nstratum>4 %then %do;
		by _stratum_;
	%end;
	run;
	quit;

/* 4 Correlation coefficient between KME(Xi,Yi)_pop and estimated S(t,Xi,Yi)_pop */

	proc reg data=plotfit noprint outest=R_corr (keep=_stratum_ _RSQ_) rsquare;
	model s_p=survival;
	by _stratum_;
	run;
	quit;

	data R_corr; 
	attrib _stratum_ &KMvar label='' ;
	merge _leg R_corr;
	_Rcorr_=sqrt(_rsq_);
	label _stratum_='Stratum Number'
	      _rsq_='R-squared'
	      _Rcorr_="Pearson's correlation statistic";
	run;

	proc print data=R_corr noobs l;
	title " Correlation statistic between estimated and observed marginal survival functions";
	run;
	
/* 5  Plot estimated S(t,Xi,Yi)_pop against KME(Xi,Yi)_pop */

	goptions reset=all;
	axis1 label=(h=1 angle=90 rotate=0 " Cure model estimate" )
	      order=1 to 0 by -0.1
	  	  major=(height=1 width=1)
          minor=(number=2 color=black height=0.5 width=1);
	axis2 label=(h=1 rotate=0 " KME estimate" )
		  order=1 to 0 by -0.1
	  	  major=(height=1 width=1)
          minor=(number=2 color=black height=0.5 width=1);
	axis3 label=none
		  order=1 to 0 by -0.1
          major=(height=1 width=1)
          minor=(number=2 color=black height=0.5 width=1);
	%if &nstratum<5 %then %do;
	%do nst2=1 %to &nstratum;
		%let cs=%scan(&colour,&nst2, %str( ));
		symbol&nst2 i=join c=&cs;
	%end;
	%end;
	%if &nstratum > 4 %then %do;
		symbol1 i=join c=black;
	%end;
	%let nst2=%eval(&nstratum+1);
         symbol&nst2 i=join w=2 c=yellow;

	proc gplot data=plotfit;
	plot s_p*survival=_stratum_  / vaxis=axis1 haxis=axis2 cframe=ligr;
	plot2 survival*survival /  vaxis=axis3 haxis=axis2 cframe=ligr;
	%if &nstratum>4 %then %do;
		by _stratum_;
	%end;
	run;
	quit;

	%endp :;
	%if &nokm=1 %then %do;
		%put *** PLOTFIT : ;
		%put ***  Warning :  One of the covariates may be a continous variable or its value is not observed ****; 
		%put ***   No plot or statistic will be displayed ****;
	%end;


%mend plotfit;

/*************************************************************************/
/*                            JACKNIFE                                   */ 
/*************************************************************************/

%macro jacknife;

data JACKDATA/view=JACKDATA;
	do _sample_=1 to &nobs;
		do _i=1 to &nobs;
			if _i ne _sample_ then do;
			 _obs_=_i;
			 set sample_j point=_i;
			 output;
			end;
		end;
	end;
	stop;
run;

%mend jacknife;

/*************************************************************************/
/*                         SELECTION OF CONVERGED REPLICATES             */ 
/*************************************************************************/
%macro convrep; 

%global der;

	data BOOTDIST(drop=it crit_u); set bootdist_T  end=last;
	where _sample_ ne 0 and crit_u<&CONVCRIT;
	if last then do;
		call symput('der', trim(left(put(_n_,best8.))));
	end;
	run;

	%put   ;
	%put ** Number of sample with crit < &CONVCRIT = &der;
	%put   ;

%if (&JACK=Y or &BCA=Y) %then %do;
	%if &jackdata= %then %let jackdata=jackdist_t;

	data jackdist(drop=it crit_u); set &JACKDATA  end=last;
	where  crit_u<&CONVCRIT;
	if last then do;
		call symput('der_j', trim(left(put(_n_,best8.))));
	end;
	run;

	%put ** Number of Jacknife samples with crit < &CONVCRIT : &der_j **;
%end;

%mend convrep;

/*************************************************************************/
/*                         BOOTSTRAP CONFIDENCE INTERVALS                */ 
/*************************************************************************/

%macro bootci;
   data _ACTTR_ ; set estimates (drop=_sample_ rename=(estimate=value1));
   proc sort data=_ACTTR_; by variable;
   run;  

   *** transpose resampling distributions;

   proc sort data=BOOTDIST; by _sample_;
   run;
   proc transpose data=BOOTDIST prefix=col out=BOOTTRAN(rename=(col1=value)) name=variable;
   var var0-var&nvar;
   by _sample_;
   run;

   proc sort data=BOOTTRAN; by variable %if (&BC=Y or &BCA=Y) %then value;;
   run;
 %if (&BC=Y or &BCA=Y) %then %do;

      %if &BCA=Y %then %do;
       
          *** estimate acceleration for BCa;

         proc means data=JACKDIST noprint vardef=df;
               var var0-var&nvar; 
			   output out=JACKSKEW (drop= _type_ _freq_ ) skewness=;
         run;

         *** transpose skewness;
         proc transpose data=JACKSKEW prefix=col out=JACKSKEW (rename=(col1=skewness))name=variable;
               var var0-var&nvar;
         run;

         proc sort data=JACKSKEW; by variable ; 
         run;

      %end;

      *** estimate median bias for BC and BCA;

      data _BC_;
         retain _alpha _conf;
         drop value value1 _sample_ ;
         if _n_=1 then do;
            _alpha=&ALPHA;
            _conf=100*(1-_alpha);
            call symput('conf',trim(left(put(_conf,best8.))));
         end;
         merge _ACTTR_ BOOTTRAN;
         by variable;
         if first.variable then do; n=0; _z0=0; end;
         n+1;
         _z0+(value<value1)+.5*(value=value1);
         if last.variable then do;
            _z0=probit(_z0/n);
            output;
         end;
      run;

      *** compute percentiles;
      data BOOTPCTL_BC  ;
         retain _i 
			   %if &BC=Y %then  _loBC _upBC _nploBC _jloBC _gloBC _npupBC _jupBC _gupBC alcl_BC aucl_BC ; 
               %if &BCA=Y %then  _loBCA _upBCA _nploBCA _jloBCA _gloBCA _npupBCA _jupBCA _gupBCA alcl_BCA aucl_BCA ;;
         drop _alpha _sample_ _conf _i value
              %if &BC=Y %then  _nploBC _jloBC _gloBC _npupBC _jupBC _gupBC;
              %if &BCA=Y %then  _nploBCA _jloBCA _gloBCA _npupBCA _jupBCA _gupBCA ;;
         merge BOOTTRAN _BC_ %if &BCA=Y %then JACKSKEW;;
         by variable;
         label _z0='Bias Correction (Z0)'
		 		%if &BC=Y %then
			    _loBC='Lower BC Percentile Point'
                _upBC ='Upper BC Percentile Point'
			    alcl_BC="Method=BC Lower Confidence Limit"
			    aucl_BC="Method=BC Upper Confidence Limit";
				%if &BCA=Y %then
			    _loBCA='Lower BCA Percentile Point'
                _upBCA ='Upper BCA Percentile Point'
			    alcl_BCA='Method=BCA Lower Confidence Limit'
			    aucl_BCA='Method=BCA Upper Confidence Limit' ;;
              
         if first.variable then do;
            %if &BC=Y %then %do;
               _loBC=probnorm(_z0+(_z0+probit(_alpha/2)));
               _upBC=probnorm(_z0+(_z0+probit(1-_alpha/2)));

      	       _nploBC=min(n-.5,max(.5,fuzz(n*_loBC)));
               _jloBC=floor(_nploBC); _gloBC=_nploBC-_jloBC;
               _npupBC=min(n-.5,max(.5,fuzz(n*_upBC)));
               _jupBC=floor(_npupBC); _gupBC=_npupBC-_jupBC;
               _i=0;
            %end;

            %if &BCA=Y %then %do;
               drop skewness;
               retain _accel;
               label _accel='Acceleration';
               _accel=skewness/(-6*sqrt(&nobs))*
                      (&nobs-2)/&nobs/sqrt((&nobs-1)/&nobs);
               _i=_z0+probit(_alpha/2);
               _loBCA=probnorm(_z0+_i/(1-_i*_accel));
               _i=_z0+probit(1-_alpha/2);
               _upBCA=probnorm(_z0+_i/(1-_i*_accel));

			   _nploBCA=min(n-.5,max(.5,fuzz(n*_loBCA)));
               _jloBCA=floor(_nploBCA); _gloBCA=_nploBCA-_jloBCA;
               _npupBCA=min(n-.5,max(.5,fuzz(n*_upBCA)));
               _jupBCA=floor(_npupBCA); _gupBCA=_npupBCA-_jupBCA;
               _i=0;
           %end;  
         end;

         _i+1;

         %if &BCA=Y %then %do;
         	if _gloBCA then do;
           		if _i=_jloBCA+1 then alcl_BCA=value;
        	end;
        	else do;
           		if _i=_jloBCA then alcl_BCA=value;
            	else if _i=_jloBCA+1 then alcl_BCA=(alcl_BCA+value)/2;
         	end;
         	if _gupBCA then do;
            	if _i=_jupBCA+1 then aucl_BCA=value;
        	 end;
        	 else do;
            	if _i=_jupBCA then aucl_BCA=value;
                else if _i=_jupBCA+1 then aucl_BCA=(aucl_BCA+value)/2;
         	end;
		 %end;

        %if &BC=Y %then %do;
			if _gloBC then do;
            	if _i=_jloBC+1 then alcl_BC=value;
         	end;
         	else do;
            	if _i=_jloBC then alcl_BC=value;
            	else if _i=_jloBC+1 then alcl_BC=(alcl_BC+value)/2;
         	end;
         	if _gupBC then do;
            	if _i=_jupBC+1 then aucl_BC=value;
         	end;
         	else do;
            	if _i=_jupBC then aucl_BC=value;
            	else if _i=_jupBC+1 then aucl_BC=(aucl_BC+value)/2;
         	end;
		%end;

         if last.variable then do;
            output;
         end;
      run;

%end;

%if (&PCTL=Y or &HYB=Y) %then %do;
      %let pctlpre=a;
      %if &PCTl=Y %then %do;
      	%let pctlname_PCTL=lcl_PCTL ucl_PCTL;
	  %end;
	  %if &HYB=Y %then %do;
	  	%let pctlname_HYB=lcl_HYB ucl_HYB;
	  %end;

      data _null_;
         _alpha=&ALPHA;
         _conf=100*(1-_alpha);
         call symput('conf',trim(left(put(_conf,best8.))));
         %if &PCTL=Y %then %do;
            _loPCTL=100*_alpha/2;
            _upPCTL=(100-_loPCTL);
         call symput('pctlpts_PCTL',trim(left(put(_loPCTL,best8.)))||' '||
                               trim(left(put(_upPCTL,best8.))));

         %end;
         %if &HYB=Y %then %do;
            _upHYB=100*_alpha/2;
            _loHYB=(100-_upHYB);
         call symput('pctlpts_HYB',trim(left(put(_loHYB,best8.)))||' '||
                               trim(left(put(_upHYB,best8.))));

         %end;
         
      run;

  %if &PCTL=Y %then %do;
      proc univariate data=BOOTTRAN noprint pctldef=5;
         var value;
         output out=BOOTPCTL_PCTL n=n
            pctlpts=&pctlpts_PCTL pctlpre=&pctlpre pctlname=&pctlname_PCTL;
         by variable;
      run;
  %end;
  %if &HYB=Y %then %do;
      proc univariate data=BOOTTRAN noprint pctldef=5;
         var value;
         output out=BOOTPCTL_HYB n=n
            pctlpts=&pctlpts_HYB pctlpre=&pctlpre pctlname=&pctlname_HYB;
         by variable;
      run;
  %end;

%end;

%if (&JACK=Y or &BOOTN=Y) %then %do;
   %if &JACK=Y %then %let dist=JACKDIST;
   %if &BOOTN=Y %then %let dist=BOOTDIST;


  %compute:;

   *** compute mean, std, min, max of resampling distribution;
   proc means data=&dist (drop=_sample_) noprint %if &dist=JACKDIST %then vardef=n; %else vardef=df;;
    var var0-var&nvar;
      output out=_TMP2_(drop=_type_ _freq_);
   run;

   *** transpose statistics for resampling distribution;
   proc transpose data=_TMP2_  out=_TMP2_ name=variable;
   var var0-var&nvar;
   id _stat_;
   run;

   proc sort data=_TMP2_; by variable;
   run;

   data CI_&dist (rename=(%if &dist=BOOTDIST %then mean=bootmean aucl=aucl_BOOTN alcl=alcl_BOOTN));
                            %if &dist=JACKDIST %then mean=Jackmean aucl=aucl_JACK  alcl=alcl_JACK));;
      retain variable value mean stderr alcl biasco aucl confid method min max n;
      merge _ACTTR_(rename=(value1=value)) _TMP2_(rename=(std=stderr));	
      by variable;
         length method $20;
         retain z; drop z;
         if _n_=1 then do;
            z=probit(1-&ALPHA/2); 
            confid=&conf;
            %if &dist=BOOTDIST %then method='Bootstrap Normal'; 
			%if &dist=JACKDIST %then method='Jackknife';;
         end;
	  %if &dist=JACKDIST %then %do;
      	stderr=stderr*sqrt(&nobs-1);
     	bias=(mean-value)*(&nobs-1);
	  %end;
	  %if &dist=BOOTDIST %then %do;
       	bias=mean-value;
	  %end;
       biasco=value-bias;
          alcl=biasco-z*stderr;
          aucl=biasco+z*stderr;
      label variable  ='variable'
            value ='Observed Statistic'
            biasco='Bias-Corrected Statistic'
			%if &dist=BOOTDIST %then %do;
              alcl  ='Method=BOOTN  Lower Confidence Limit'
              aucl  ='Method=BOOTN  Upper Confidence Limit'
              mean='Bootstrap Mean'
              bias  ='Approximate Bias'
              stderr='Approximate Standard Error'
			%end;
			%if &dist=JACKDIST %then %do;
			  alcl  ='Method=JACK Lower Confidence Limit'
              aucl  ='Method=JACK  Upper Confidence Limit'
              mean='Jackknife Mean'
              bias  ='Estimated Bias'
              stderr='Estimated Standard Error'
			%end;
            confid='Confidence Level (%)'
            method='Method for Confidence Interval'
            min   ='Minimum Resampled Estimate'
            max   ='Maximum Resampled Estimate'
            n     ='Number of Resamples'
            ;
   run;

%final:;
   %if &dist=BOOTDIST %then %do;
   	%if (&JACK=Y) %then %do;
   	  %let dist=JACKDIST;
	  %goto compute;
	%end;
   %end;

%end;


%mend bootci;


/*************************************************************************/
/*             OUTPUT DATA SET                                          */ 
/*************************************************************************/

%macro output;

%let drop= _z0 n ;
%if &BC=Y %then %let drop=&drop _loBC _upBC ;
%if &BCA=Y %then %let drop=&drop _loBCA _upBCA _accel; 

   data BOOTCI;
      retain variable value confid n  
             %if &BC=Y %then alcl_BC aucl_BC ;
			 %if &BCA=Y %then alcl_BCA aucl_BCA;
			 %if &PCTL=Y %then alcl_PCTL aucl_PCTL;
			 %if &HYB=Y %then alcl_HYB aucl_HYB;
             %if &BOOTN=Y %then alcl_BOOTN aucl_BOOTN;
             %if &JACK=Y %then alcl_JACK alcl_JACK;;
      merge _ACTTR_ (rename=(value1=value))
	       %if (&BOOTN=Y) %then CI_bootdist (drop= bootmean method bias stderr biasco min max n confid );
		   %if (&JACK=Y) %then CI_jackdist  (drop= jackmean method bias stderr biasco min max n confid );
		   %if (&PCTL=Y) %then BOOTPCTL_PCTL ;
		   %if (&HYB=Y) %then BOOTPCTL_HYB ;
           %if (&BC=Y or &BCA=Y) %then BOOTPCTL_BC(drop= &drop);;
      by variable;

      %if &HYB=Y %then %do;
         aucl_HYB=2*value-aucl_HYB;
         alcl_HYB=2*value-alcl_HYB;
      %end;
      confid=&conf;
      label variable  ='Variable'
            value ='Observed Statistic'
            confid='Confidence Level (%)'
            n     ='Number of Resamples'
	        %if &PCTL=Y %then
	        alcl_PCTL="Method=PCTL Lower percentile"
			aucl_PCTL="Method=PCTL Upper percentile";
			%if &HYB=Y %then;
	        alcl_HYB="Method=HYB Lower percentile"
			aucl_HYB="Method=HYB Upper percentile";;
   run;


  data BOOTCI; set BOOTCI;
  format Variable $10.;
  if variable='var0' then variable='L_Int';
  %do i=1 %to &nlvar;
	if variable="var&i" then variable="L_%scan(&var,&i,%str( ))";
  %end;
  %do i=&s %to &nvar;
	if variable="var&i" then variable="S_%scan(&VAR,&i,%str( ))";
  %end;
  run;


%IF &LINK=LOGIT %then %do;
%local meth;
%let i=1;
data boot_OR (drop= value); set BOOTCI;
if substr(variable, 1,1)='L' and substr(variable,1,5) ne 'L_Int';
ODD=exp(value);
label ODD='Odds Ratio';
	%do %while(%length(%scan(&BOOTMET,&i,%str( ) )));
	    %let meth=%scan(&BOOTMET,&i,%str( ) );
    	&meth._LowerOR=exp(alcl_&meth);
		&meth._UpperOR=exp(aucl_&meth);
		label &meth._LowerOR="&meth Lower OR"
	    	  &meth._UpperOR="&meth Upper OR";
		drop aucl_&meth alcl_&meth;
		%let i=%eval(&i+1);
	%end;
run;
%let i=1;
data boot_RR (drop=value); set BOOTCI;
if substr(variable, 1,1)='S' ;
RelRisk=exp(value);
label RelRisk='Hazard Ratio';
	%do %while(%length(%scan(&BOOTMET,&i,%str( ) )));
	    %let meth=%scan(&BOOTMET,&i,%str( ) );
    	&meth._LowerRR=exp(alcl_&meth);
		&meth._UpperRR=exp(aucl_&meth);
		label &meth._LowerRR="&meth Lower RR"
	    	  &meth._UpperRR="&meth Upper RR";
		drop aucl_&meth alcl_&meth;
		%let i=%eval(&i+1);
	%end;
run;
%end;

ods listing ;

	proc print data=BOOTCI (drop= confid n) noobs l;
	id variable value;
	title "BOOTSTRAP CONFIDENCE INTERVAL FOR PARAMETERS ESTIMATES";
	title2 "Data set= &data";
	title3 "(confidence level=&conf %, &der bootstrap resamples)";
	run;


	%IF &LINK=LOGIT %then %do;
		proc print data=boot_OR (drop= confid n) noobs l;
	    id variable ODD;
		title "ODDS RATIO FOR THE LOGISTIC PART";
		title2 "Data set= &data";
	    title3 "(confidence level=&conf %, &der bootstrap resamples)";

		proc print data=boot_RR (drop= confid n) noobs l;
		id variable relrisk;
		title "HAZARD RATIO FOR THE SURVIVAL PART";
		title2 "Data set= &data";
	    title3 "(confidence level=&conf %, &der bootstrap resamples)";
		run;
	%end;


%mend output;


/*************************************************************************/
/*     QQ PLOTS AND DISTRIBUTION OF PARAMETETS ESTIMATES                 */ 
/*************************************************************************/
%macro graph;
	goptions reset=global Transparency NoBorder NoPrompt
             gunit=pct 
             ftitle=swissb ftext=swiss htitle=3 htext=2;

%let i=1;
%let _vec = L_int;
%do %while(%length(%scan(&var, &i, %str( ))));
	%let _first=%scan(&var, &i, %str( ));
	%if &i<&s %then %do;
	%let _vec = &_vec L_&_first;
	%end;
	%if &i>&nlvar %then %do;
	%let _vec = &_vec S_&_first;
	%end;
	%let i=%eval(&i+1);
%end;


Proc Univariate Data=bootdist noprint;
 Var &_vec;
 Histogram / Normal(mu=est sigma=est color=blue noprint) cframe=ligr ;
 title "distribution of estimates";
 title2 "over &der bootstap replicates";
 QQPlot / Normal(mu=est sigma=est color=yellow l=1 w=2 )cframe=ligr ;
Run;
Quit;



%mend graph;

/*************************************************************************/
/*     RENAME VARIABLE NAMES IN OUTPUT DATASETS                          */ 
/*************************************************************************/

%macro outset;

data estimates; set estimates;
	format variable $15.; label variable='variable';
    if variable='var0' then variable='Inc.Intercept'; 
	%if &lvars ne %then %do;
    	%do i=1 %to &nlvar;
			if variable="var&i" then variable="L_%scan(&VAR, &i, %str( ) )";
		%end;
	%end;
	%if &svars ne %then %do;
		%do i=&s %to &nvar;
			if variable="var&i" then variable="S_%scan(&VAR, &i, %str( ) )";
		%end;
	%end;
run;

%let datname=probcure; 
%if &BASELINE =Y %then %let datname= &datname _baseline;
%if &PLOTFIT=Y %then %let datname =&datname plotfit;
%if &BOOTSTRAP=Y %then %do;
	%let datname= &datname bootdist bootdist_t; 
	%if &BASELINE=Y %then %let datname=&datname _baseline_t;
	%if &JACK=Y and &JACKDATA= %then %let datname=&datnane Jackdist_t;
%end;

%let u=1;
%do %while(%length(%scan(&datname, &u, %str( ) ))); 
	%let dname= %scan(&datname, &u, %str ( ));
	data &dname; set &dname;
   	 rename var0=L_int;
	 %if &lvars ne %then %do;
   	 	%do i=1 %to &nlvar;
    	    rename var&i=L_%scan(&VAR,&i,%str( ));
			label var&i="L_%scan(&VAR,&i,%str( ))";
		%end;
	 %end;
	 %if (&svars ne ) %then %do;
   		%do i=&s %to &nvar;
	  		rename var&i=S_%scan(&VAR,&i,%str( ));
			label var&i="S_%scan(&VAR,&i,%str( ))";
		%end;
	 %end;
    run;
	%let u=%eval(&u+1);
%end;

%if &BASELINE=Y %then %do;
	data baseline; set _baseline;
	run;
%end;
%if &BOOTSTRAP=Y %then %do;
	%if &BASELINE=Y %then %do;
		data baseline_t; set _baseline_t;
		run;
	%end;
%end;



%mend outset;

/*************************************************************************/
/*                   %PSPMCM                                             */ 
/*                   put it all together                                 */
/*************************************************************************/


%macro pspmcm(DATA=,SURVPART=, AFT= ,ID=,CENSCOD=,TIME=,
			  VAR=, 
              INCPART=, 
			  TAIL=, SU0MET=,
			  FAST= ,BOOTSTRAP=,
			  NSAMPLE=, STRATA=,
			  MAXITER=,CONVCRIT=,ALPHA= , 
			  BOOTMET=, JACKDATA=,
              GESTIMATE=,
			  BASELINE=,
			  SPLOT= ,
              PLOTFIT= );

option  nonotes nomlogic nomprint nosymbolgen nosource;
*options notes  mprint source symbolgen mlogic  ;
options formdlim=' ' nodate nonumber  ls=100;

ods listing;

/* Change to uppercase */
%let SURVPART       = %qupcase(&SURVPART);
%let AFT=			= %qupcase(&AFT);
%let ID             = %qupcase(&ID);
%let CENSCOD        = %qupcase(&CENSCOD);
%let TIME           = %qupcase(&TIME);  
%let VAR            = %qupcase(&VAR);
%let LINK           = %upcase(&INCPART);
%let BOOTSTRAP      = %qupcase(&BOOTSTRAP);
%let strata         = %qupcase(&strata);
%let BASELINE       = %qupcase(&BASELINE);
%let GESTIMATE      = %qupcase(&GESTIMATE);
%let TAIL           = %qupcase(&TAIL);
%let FAST           = %qupcase(&FAST);
%let SU0MET         = %qupcase(&SU0MET);
%let BOOTMET        = %qupcase(&BOOTMET);
%let SPLOT          = %qupcase(&SPLOT);
%let PLOTFIT        = %qupcase(&PLOTFIT);


/* Default values  */
%let HYB=N ; %let PCTL=N; %let BOOTN=N; %let JACK=N ; %let BCA=N ; %let BC=N;
%let jack_n=0;
%let conv_0=0;
%let crit_0=0;
%let done_0=0;
%let nokm=0;


%if &INCPART=  %then  %let LINK=LOGIT;
%if &MAXITER= %then %let maxiter=200;
%if &TAIL= %then %let tail=ZERO;
%if &SU0MET= %then %let SU0MET=PL;
%if &CONVCRIT= %then %let CONVCRIT=1e-5;
%if &ALPHA= %then %let ALPHA=0.05;
%if &GESTIMATE ne Y %then %let GESTIMATE = N;
%if &BASELINE ne Y %then %let BASELINE= N;
%if &BOOTMET= %then %let BOOTMET=N;
%if (&SURVPART ne COX) %then %do;
	%let BOOTSTRAP=N; 
	%let tail=;
%end;
%if &BOOTSTRAP ne Y %then %do;
	%let nsample=0;
	%let BOOTMET=N;
	%let GESTIMATE= N; 
%end;
%if ((&BCA=Y or &JACK=Y) AND &JACKDATA= ) %then %let JACKDATA=Jackdist_t;
%if &SPLOT = %then %let SPLOT=N;
%if &SPLOT=Y %then %let BASELINE=Y;
%if &PLOTFIT=Y %then %do;
	%let BASELINE=Y;
	%let SPLOT=Y;
%end;
%if &PLOTFIT= %then %let PLOTFIT=N;

/*Time at the begenning */
data _null_;
time1=time();
call symput('time1',trim(left(time1)));
run;

*** compute confidence level;
data _null_;
conf=100*(1-&ALPHA);
call symput('conf',trim(left(put(conf,best8.))));
run;

/* Ckeck information*/
%prepare
%if (&exiterr ne 0) %then %goto exit;

/* Parametric mixture cure models */
%if (&SURVPART2=PARA) %then %do;
	 %parametric
	 %if (&BASELINE=Y) %then %do;
	   %baseline
	 %end;
     %if (&SPLOT=Y) %then %do;
       %survplot 
	%end;
	%if (&PLOTFIT=Y) %then %do;
       %plotfit 
	%end;
	%goto exit;
%end;

/* Error messsage if the specified number of replicates is not high enough to compute bootstrap CI */
%if (&BOOTSTRAP=Y and &BOOTMET ne N) %then %do;
	data _null_;
	n=1/(&ALPHA/2);
	call symput('nes',left(n));
	run;

	%if (&nes>&NSAMPLE) AND (&NSAMPLE>=0) %then %do;
		%put ** The number of bootstrap sreplicates specified is too low to allow computation of bootstrap Confidence Intervals **;
		%let exiterr=1;
		%goto exit;
	%end;
%end;

/* If Bootstrap=N then number of sample=0 ans sample0= intitial sample */
%if &BOOTSTRAP NE Y %then %do;
	%let nr=1;
	%let nb=0;
	%let msg1=original sample;
	data sample; set sample0; 
	run;
%end;


/*If bootstrap=Y */
%if &BOOTSTRAP=Y %then %do;

 /* Number of loop needed for the specified number of replicates (max size of data set set to 5.10e5 records) */
	data _null_;
	n=&NSAMPLE*&nobs;
	nr=ceil(n/5e5);
	nb=ceil(&NSAMPLE/nr);
	call symput('nr', trim(left(nr)));
	call symput('nb', trim(left(nb)));
	run;

	%if &nr>1 %then %do;
		%put **The dataset containing &NSAMPLE bootstrap replicates would by very large   **;
		%put **SAS will resampling in &nr times of &nb replicates each **;
	%end;
%end;

/* EM algorithm */
%let run=1;

%do %until(&run>&nr);

	%if &BOOTSTRAP=Y %then %do;
		%let nsample_b=&nb;
		%let seed=%eval(1234475+&run);
		%let msg1= &nb bootstrap replicates;
		%let boot_res=bootdist_t;
		%bootstrap
	%end;

	%put ;
	%put ** SAS will now compute estimates for &msg1;
	%max

	%let it=0;
	%let crit_g=0;
	%let nb_sep=0;
	%let conv_g=0;

	%do %until (&conv_g=1 or &it>=&MAXITER);
		%let it=%eval(&it+1);
		%let crit_u=0;
		%let conv_u=0;
		%let conv_1=0;
		%let nb_nc=0;
		%let nb_c=0;
		%maximisation
		%expectation
		%if &it>1 %then %do;
			%compare
		%end;
		%end_it
		%if ((&conv_0=1  or &it>=&MAXITER) and &done_0 ne 1)  %then %do;
			%result_0
			%if (&BASELINE=Y) %then %do;
				%baseline
			%end;
			%if (&SPLOT=Y) %then %do;
				%survplot
			%end;
			%if (&PLOTFIT=Y) %then %do;
				%plotfit
			%end;
			%let done_0=1;
		%end;
		%if ((&conv_1=1 or &it>=MAXITER) and &BOOTSTRAP=Y ) %then %do;
			%result
			%if (&BASELINE=Y) %then %do;
				%baseline
			%end;
		%end;
	%end;
    proc datasets lib=work nolist; delete comp&it;
    run;
	
	%let run=%eval(&run+1);

%end;


/* Itération for BCA ou JACKKNIFE*/
%if (&JACKDATA= and (&BCA=Y or &JACK=Y)) %then %do; 
		%put ** Creating the jackknifed replicates.....;
		%put ;

		proc sort data=sample0 out=sample_j(drop=_obs_ _sample_); by &TIME descending &CENSCOD;
		run;

	%jacknife;
	
	/* Calculating the number of run(s) needed to compute jackknife estimates */
	
	data _null_;
 	nb_j=ceil(1e6/(&nobs-1));
 	nr_j=ceil(&nobs/nb_j);
 	call symput('nb_j',trim(left(nb_j)));
 	call symput('nr_j',trim(left(nr_j)));
	run;

	%if &nr_j>1 %then %do;
		%put ** The number of jackknife replicates is very large **;
		%put ** Then macro will perform &nr_j loops of &nb_j replicates each**;
	%end;
	%else %do ; %let nb_j=&nobs; %let nr_j=1; %end;

	%let msg1=Jackknife re_samples;
	%let boot_res=jackdist_t;
    %let run_j=1;

	%do %until (&run_j>&nr_j);
	%put ** Calculating Jacknife estimates : loop &run_j/&nr_j (&nb_j/&nobs) replicates) **;
    
 		data sample ; set jackdata;
		%if &nr_j>1 %then %do;
		   %if (&run_j<&nr_j) %then %do;
 			  if &nb_j*(&run_j-1)+1=<_sample_=<(&nb_j*&run_j);
 		   %end;
 		   %if (&run_j=&nr_j) %then %do;
			  if _sample_>=&nb_j*(&run_j-1)+1;
 		   %end;
		%end;
 		run;
		
		%max

		%let it=0;
		%let crit_g=0;
		%let nb_sep=0;
		%let conv_g=0;

		%do %until (&conv_g=1 or &it>=&MAXITER);
			%let it=%eval(&it+1);
			%let conv_1=0;
			%let nb_nc=0;
			%let nb_c=0;
			%maximisation
			%expectation
			%if &it>1 %then %do;
				%compare
			%end;
			%end_it
			%if (&conv_1=1) %then %do;
				%result
			%end;
		%end;
	
		%let run_j=%eval(&run_j+1);
	%end; 
%end;


/* Selection of converged replicates for CI intervals compuation and graph */

%if &BOOTSTRAP=Y %then %do;
	%convrep
%end;

/* Bootstrap confidence interval computation*/
%if &BOOTMET ne N %then %do;
 %bootci
 %output
%end;

/*remane variable in output datasets */
%if (&SURVPART=COX) %then %do;
 %outset
%end;

/* QQ plot and distribution of parameters estimates*/
%if &GESTIMATE=Y %then %do;
 %graph
%end;



data _null_ ;
format time_tot time.;
time2=time();
time_tot= time2-&time1;
put 'Total time:' time_tot;
run;

%exit :;
%if (&exiterr ne 0) %then %do;
	%put ** The macro exited due to errors.**;
%end;

%else %do;
proc datasets lib=work nolist;
delete sample0 
%if &SURVPART=COX %then bas compare conv b lp l lp2 min_max n_conv sample sep sortie1 sp sp2   _cov 
%if (&BASELINE=Y and &BOOTSTRAP=Y) %then bas_boot _baseline_t; 
%if &BOOTSTRAP=Y %then fvar _acttr_ _boot_  ci_bootdist bootdist_t boot boottran;
%if &TAIL=WTAIL %then weib weib1;
%if &BOOTN=Y or &jack=Y %then _tmp2_ ;
%if &BCA=Y %then jackskew _jack_ sample_j bootpctl_bc;
%if &BC=Y %then boottran _bc_ bootpctl_bc;
%if &HYB=Y %then bootpctl_hyb;
%if &PCTL=Y %then bootpctl_pctl;
%if &JACK=Y %then jackdist_t ;
%if (&SPLOT=Y) or (&BASELINE=Y) %then  _baseline ;
%if (&SPLOT=Y) %then _graphkmu _graphkm graph;
%if &SURVPART2=PARA %then param ;
%if &PLOTFIT=Y and &nokm=0 %then cm _leg;;
run;
quit;
%end;
options notes source;

%mend pspmcm;

/*******************
   Macro 
********************

%pspmcm(DATA=simultex2,ID=,CENSCOD=status,TIME=survtime,
				  VAR= varX(IS, 0) varY1(IS, 0) varY2(IS,0) ,
				  INCPART=logit,
				  SURVPART=cox, 
				  TAIL=zero , SU0MET=pl,
				  FAST=Y,BOOTSTRAP=Y,
				  NSAMPLE=2000, STRATA=,
				  MAXITER=200,CONVCRIT=1e-5, ALPHA=0.05, 
				  BASELINE=Y, 
				  BOOTMET=ALL,
				  JACKDATA=,
				  GESTIMATE=Y,
                  SPLOT=Y,
                  PLOTFIT=Y);				  
run;



