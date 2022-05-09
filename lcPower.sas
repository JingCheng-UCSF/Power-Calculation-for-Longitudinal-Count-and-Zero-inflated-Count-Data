*option symbolgen;
**********************************************************************************
Macro lcPower to calculate power for longitudinal correlated count
and zero-inflated count outcome data.
Parameters:
ds= specifies the input simulation data set generated from macro simulateLC
ydist= Specifies outcome distribution. Choose from poisson, negbin, zip and zinb
mzi= specifies GLMM or Marginal ZI model for zip or zinb distributions. mzi=1 for
Marginalized ZI and mzi=0 for GLMM. The default is 0 (GLMM). Not applicable 
if ydist = poisson or negbin
gconv= specifies the relative gradient convergence criterion, 
required if ydist=zip or zinb. The default value is 1E–8. 
outds= specifies the output data set that saves the power(s) for interaction of 
       treatment by time
Written by Nancy Cheng
Reviewed by Jing Cheng
Date: 01/28/2021 
**********************************************************************************;
%macro lcPower(ds=, ydist=, mzi=0, gconv=1E-8, outds=);
proc datasets;
	delete _pe_ _pv_;
quit;

%let ydist=%sysfunc(upcase(&ydist));

%if &ydist=POISSON or &ydist= NEGBIN %then %do;
	proc glimmix method=quad data=&ds;
		title "Fit GLMM for &ydist outcome";
		by iSim;
		class iSite iSubject;
		model y = group time group*time x1 logx2/ dist=&ydist link=log cl;
		random int / subject=iSubject(iSite) type=vc ;
		ods output ParameterEstimates=_pe_;
	run;
	data _pv_;
		set _pe_;
		if Probt=. then rej=.;
		else rej=(Probt<=0.05);
		if effect='group*time' then output;
		keep isim Probt rej;
	run;
	*get mean rejection number;
	title '';
	proc means data=_pv_ n mean maxdec=5;
		var rej;
		ods output summary=&outds;
	run;

	data &outds;
		set &outds;
		label rej_mean='Power';
		rename rej_mean=power;
		keep rej_mean
	run;
	title 'Power Estimate';
	proc print data=&outds label noobs;run;
	title ' ';
%end;
%else %if &ydist=ZIP or &ydist=ZINB %then %do;
	%if &ydist=ZIP %then %do;
		%if &mzi=0 %then %do;
		title "Fit GLMM ZIP model";
		proc nlmixed data=&ds gconv=&gconv nthreads=-1;
			by iSim;
			parms gamma0=0 gamma1=0 gamma2=0 gamma3=0 gamma4=0 gamma5=0
			   	  alpha0=0 alpha1=0 alpha2=0 alpha3=0 alpha4=0 alpha5=0
				  sigma1=0.3 sigma12=-0.1 sigma2=0.3;
			bounds sigma1 sigma2 >= 0;
			/* linear predictor for the zero-inflation probability */
			logit_psi = gamma0 + gamma1*time + gamma2*group*time + gamma3*group + gamma4*x1 
					    + gamma5*logx2 + c1;

			psi = exp(logit_psi)/(1+exp(logit_psi));
			/*  mean \nu for Poisson Dist*/
			log_nu = alpha0 + alpha1*time + alpha2*group*time + alpha3*group + alpha4*x1 
					+ alpha5*logx2 + d1;
			
			lambda = exp(log_nu);

			if y=0 then loglike = log(psi + (1-psi)*exp(-lambda));
			else        loglike = log(1-psi) + y*log(lambda) - lambda - lgamma(y+1);

			model y ~ general(loglike);
			random c1 d1 ~ normal([0,0], [sigma1,sigma12,sigma2]) subject=iSubject;
			ods output ParameterEstimates=_pe_; 
		run;
		%end;
		%else %do; /*fit Marginalized ZIP*/
		title "Fit Marginalized ZIP";
		proc nlmixed data=&ds gconv=&gconv nthreads=-1;
			by iSim;
			parms gamma0= 0 gamma1= 0  gamma2=0 gamma3=0 gamma4=0 gamma5=0
			   	  alpha0= 0 alpha1= 0 alpha2= 0 alpha3= 0 alpha4= 0 alpha5= 0
				  sigma1=0.3 sigma12=-0.1 sigma2=0.3;
			bounds sigma1 sigma2 >= 0;
			/* linear predictor for the zero-inflation probability */
			logit_psi = gamma0 + gamma1*time + gamma2*group*time + gamma3*group + gamma4*x1 
					+ gamma5*logx2+c1;
			psi = exp(logit_psi)/(1+exp(logit_psi));
			/* Overall mean \nu */
			log_nu = alpha0 + alpha1*time + alpha2*group*time + alpha3*group + alpha4*x1 
					+ alpha5*logx2+d1;
			delta = log((1-psi)**(-1)) + log_nu;
		    /* Build the mZIP + RE log likelihood */
		    if y=0 then
		        loglike = log(psi + (1-psi)*(exp(-exp(delta))));
		    else loglike = log(1-psi) - exp(delta) + y*(delta) - lgamma(y + 1);

			model y ~ general(loglike);

			random c1 d1 ~ normal([0,0], [sigma1,sigma12,sigma2]) subject=iSubject;
			ods output ParameterEstimates=_pe_; 
		run;
		%end;
    %end;    
	%else %do; /*ZINB outcome*/
		%if &mzi=0 %then %do;
		title "Fit GLMM ZINB model";
		proc nlmixed data=&ds gconv=&gconv nthreads=-1;
			by iSim;
			parms gamma0=0 gamma1=0 gamma2=0 gamma3=0 gamma4=0 gamma5=0
			   	  alpha0=0 alpha1=0 alpha2=0 alpha3=0 alpha4=0 alpha5=0
				  sigma1=0.3 sigma12=-0.1 sigma2=0.3 r=1;
			bounds sigma1 sigma2 >= 0;
			/* linear predictor for the zero-inflation probability */
			logit_psi = gamma0 + gamma1*time + gamma2*group*time + gamma3*group + gamma4*x1 
					    + gamma5*logx2 + c1;
			psi = exp(logit_psi)/(1+exp(logit_psi));
			/*  mean \nu for Poisson Dist*/
			log_nu = alpha0 + alpha1*time + alpha2*group*time + alpha3*group + alpha4*x1 
					+ alpha5*logx2 + d1;	

			/* Model for the expectation under the negative binomial dist;*/
			lambda = exp(log_nu);
			/* for NB; */ 
			p = r/(r+lambda);
			
			/* Build the ZINB log likelihood */
			/* log-likelihood is log of appropriate density function */
		    if y=0 then loglike = log(psi + (1-psi)*(p**r));
		    else loglike = log(1-psi) + lgamma( y + r) - lgamma(y + 1)
		   			      - lgamma(r) + r*log(p) + y*log(1-p);

			model y ~ general(loglike);
			random c1 d1 ~ normal([0,0], [sigma1,sigma12,sigma2]) subject=iSubject;
			ods output ParameterEstimates=_pe_; 
		run;
		%end;
		%else %do; /*fit Marginalized ZINB*/
		proc nlmixed data=&ds  gconv=&gconv nthreads=-1; 
			title "Fit marginalized ZINB model";

			by iSim;
			parms gamma0= 0 gamma1= 0  gamma2=0 gamma3=0 gamma4=0 gamma5=0
			   	  alpha0= 0 alpha1= 0 alpha2= 0 alpha3= 0 alpha4= 0 alpha5= 0
				  phi=1   sigma1=0.3 sigma12=-0.1 sigma2=0.3;
			bounds sigma1 sigma2 >= 0;
			/* linear predictor for the zero-inflation probability */
			linpinfl = gamma0 + gamma1*time + gamma2*group*time + gamma3*group + gamma4*x1 
					+ gamma5*logx2+c1;
			psi = 1 / (1 + exp(-linpinfl));
			/* Overall mean \nu */
			nu = exp(alpha0 + alpha1*time + alpha2*group*time + alpha3*group + alpha4*x1 
					 + alpha5*logx2+d1);

			mu = nu/(1-psi);
			alpha = 1/phi;
		    theta = 1/(1+(mu/alpha));

			/* Build the mZINB + RE log likelihood */
			if y=0 then
			     loglike = log(psi + (1-psi)*(theta**alpha));
			else loglike = log(1-psi) + lgamma(y+alpha) - lgamma(alpha)
		                   + y*log(1-theta)+alpha*log(theta) - lgamma(y+1);

			model y ~ general(loglike);
			random c1 d1 ~ normal([0,0], [sigma1,sigma12,sigma2]) subject= iSubject;
			ods output ParameterEstimates=_pe_; 
		run;
		%end;

	%end;
	*compute rejection numbers;
	data _pe_;
		set _pe_;
		if Probt=. then rej=.;
		else rej=(Probt<=0.05);
	run;
	proc sql;
	create table _pv_ as
	select a.iSim, a.Probt as p_zi label='Pr > |t| from the logit part of the model', a.rej as rej_zi, 
	b.Probt as p_cnt label='Pr > |t| from the count part of the model', 
		b.rej as rej_cnt, 
		case when (a.rej + b.rej>=1) then 1
			when (a.rej + b.rej=0) then 0
		    else . end as rej_zi_or_cnt,
		case when (a.rej + b.rej=2) then 1
			when (a.rej + b.rej=0) or (a.rej + b.rej=1)  then 0
		    else . end as rej_zi_and_cnt
	from _pe_(where=(parameter='gamma2')) as a
	left join _pe_(where=(parameter='alpha2')) as b
	on a.iSim=b.iSim;
	quit;
	title ' ';

	*get mean rejection number;
	proc means data= _pv_ n mean maxdec=5;
		var rej_zi rej_cnt rej_zi_or_cnt rej_zi_and_cnt;
		ods output summary=&outds;
	run;
	*compute powers;
	data &outds;
		set &outds;
		*ask Jing if these labels need change!;
		label rej_zi_mean='Power to reject from the logit part of the model'
			rej_cnt_mean='Power to reject from the count part of the model'
			rej_zi_or_cnt_mean='Power to reject from the logit part or from the count part of the model '
			rej_zi_and_cnt_mean='Power to reject from the logit part and from the count part of the model ';

		rename rej_zi_mean=power_zi
			   rej_cnt_mean=power_cnt
			   rej_zi_or_cnt_mean=power_zi_or_cnt
			   rej_zi_and_cnt_mean=power_zi_and_cnt;
		keep rej_zi_mean rej_cnt_mean rej_zi_or_cnt_mean rej_zi_and_cnt_mean;
	run;
	title 'Power Estimate';
	proc print data=&outds label noobs;run;
	title ' ';
%end;
%else %do;
 put 'Warning: Outcome distribution &ydist is not supported.';
%end;
%mend;


%lcPower(ds=poi_data(where=(isim<=20)), ydist=Poisson, outds=Poi_power);
%lcPower(ds=nb_data, ydist=Negbin, outds=nb_power);
%lcPower(ds=zip_data(where=(isim<=3)), ydist=zip, gconv=0.0000001, mzi=0, outds=glmm_zip_power);
%lcPower(ds=zip_data(where=(isim<=3)), ydist=zip, gconv=0.0000001, mzi=1, outds=mzip_power);
%lcPower(ds=zinb_data(where=(isim<=3)), ydist=zinb, gconv=0.0000001, mzi=0, outds=glmm_zinb_power);
%lcPower(ds=zinb_data(where=(isim<=3)), ydist=zinb, gconv=0.0000001, mzi=1, outds=mzinb_power);



