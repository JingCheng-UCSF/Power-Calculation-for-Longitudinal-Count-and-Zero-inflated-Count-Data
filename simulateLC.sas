**************************************************************************************
macro simulateLC generates a time variable, binary treatment variable (group), a 
normal covariate x1, and log of Poisson covariate (logx2), and then generates 
longitudinal count and zero-inflated count. 
Parameters:
seed= specifies a random seed.
nSim= specifies the number of simulations.
nSite= specifies the number of study site. The default is 2.
nSubj_perSite= specifies the number of subjects within a study site.
ydist= specifies the outcome distribution, which can be poisson, negbin, zip and zinb
size= specifies the dispersion parameter for negative binomial outcome distribution. 
The default is 1
ymax= specifies the maximum value that the outcome is allowed.
x2m= specifies the Poisson mean to generate x2. The default is 87.
time_s= specifies the first study timepoint. The default is 0. 
time_e= specifies the last study timepoint. The default is 3.
sitecov= specifies a string with 6 numeric values separated by comma, eg, 
%str(1.5,0,0,0,0,0) as the upper triangular covariance matrix for site random effects
in the outcome model or in the count part of ZIP/ZINB outcome.
betas= specifies a string of 6 numeric values for the fixed effect coefficients in the
outcome model or in the count part of ZIP/ZINB outcome model. eg, 
%str(-12, 0.2, -0.1, -0.01, 0.1, 2.95), with the numbers match the order of 
intercept, time, time*group, group, x1, and logx2
gcov= specifies a string with 6 numeric values separated by comma, eg, %str(1.8,0.2,0,0,0,0)
as the upper triangular g-side covariance matrix in the outcome model or in the count part 
of ZIP/ZINB outcome model.
sitecov_zi= specifies a string with 6 numeric values separated by comma, eg, %str(1.5,0,0,0,0,0)
as the upper triangular covariance matrix for site random effects in the logit part of 
ZIP/ZINB outcome model. Required if ydist=zip/zinb
betas_zi= specifies a string of 6 numeric values for the fixed effect coefficients in the 
logit part of ZIP/ZINB outcome model. eg, %str(20, -2.04, -0.1, -0.01, -0.64, -5.73), 
with the numbers match the order of intercept, time, time*group, group, x1, and logx2
 Required if ydist=zip/zinb.
gcov_zi= specifies a string with 6 numeric values separated by comma, eg, %str(1.5,0.2,0,0,0,0)
as the upper triangular g-side covariance matrix in the logit part of ZIP/ZINB outcome model. 
Required if ydist=zip/zinb.
gzig_cov= specifies a string with 6 numeric values separated by comma, eg, 
%str(-4.17, 0.0,0.0, 0.0, 0.0, 0.0), as the upper triangular covariance matrix between 
random effects from the logit part and random effects from the count part of ZIP/ZINB outcome. 
Required if ydist=zip/zinb
outds= dataset to save the simulated data.
Written by Nancy Cheng
Date: 05/10/2019
Mod Date: 1/28/2021
**************************************************************************************;
%macro simulateLC(seed=, nSim=, nSite=2, nSubj_perSite=100, ydist=, size=1,
				  ymax=88, x2m=87, time_s=0, time_e=3, sitecov=, betas=, gcov=,
				  sitecov_zi=%str(), betas_zi=%str(), gcov_zi=%str(), gzig_cov=%str(), outds=);
	%let ydist=%sysfunc(upcase(&ydist));
    proc iml;
        call randseed(&seed);
        iSim = 1:&nSim;
        nSim = &nSim; 
        iSite = 1:&nSite;
        site = expandgrid(iSim, iSite);
        mattrib site colname = {"iSim", "iSite"};
        beta = {&betas};

        /* get site/site random effects */
        acov = sqrvech({&sitecov});
        a = randnormal(nrow(site), repeat(0, 3), acov);

        /* generate subject level random effects.  Use different subject ids for different sites in same simulation */
        iSubj = 1:(&nSite*&nSubj_perSite);

        nSubj = nSim * ncol(iSubj);
     
    	*generate binary treatment variable group;
		group = J(nSubj, 1, 0);
		do i=1 to &nSim;
			groupj_site = repeat({0 1}, 1, &nSubj_perSite/2);
			do j=1 to &nSite;
				
				groupj_ran=sample(groupj_site, &nSubj_perSite, 'WOR');
				group_site=t(groupj_ran);
				do k=1 to &nSubj_perSite;
					xpos=(i-1)*(&nSite*&nSubj_perSite) + (j-1)*&nSubj_perSite +k;
					group[xpos]=group_site[k]; *group - order by iSim, site;
				end;
			end;
		end;

		*generate normal covariate x1;
        x1 = randfun(nSubj, "Normal");
     	*generate Poisson distributed covariate x2- tooth surface count (max is 88 for a small child);
		x2 = J(nSubj, 1, 0);
		do i=1 to nSubj;
			x2[i]= randfun(1, 'Poisson', &x2m);
			if x2[i]>&ymax then i=i-1; 
		end;
	
		logx2=log(x2);
		*output matrix;
		subj = J(nSubj, ncol(site)+ncol(a)+15, 0);
		bG = randnormal(nSubj, repeat(0, 3), sqrvech({&gcov}));
	
		%if &ydist=ZIP or &ydist=ZINB %then %do;
			subj = J(nSubj, ncol(site)+ncol(a)+29, 0);
			*for zero inflation part;
		    beta_zi = {&betas_zi};

	  		acov_zi = sqrvech({&sitecov_zi});
	        a_zi = randnormal(nrow(site), repeat(0, 3), acov_zi);
			*gzig_cov - covariance between random effects (RE) from ZI part and REs from count part;
			zi_cnt_cov=sqrvech({&gcov_zi})||sqrvech({&gzig_cov})
					// (sqrvech({&gzig_cov})||sqrvech({&gcov}));
			*bG_zi= randnormal(nSubj, repeat(0, &g_ncov_zi), sqrvech({&gcov_zi}));
			bG_zi_bG= randnormal(nSubj, repeat(0, 6), zi_cnt_cov);
			bG_zi=bG_zi_bG[, 1:3];

			bG=bG_zi_bG[, 4:6];

		%end;
		
	    iOut = 1;
        iIn = 1;
        iIn2 = 1; *within simulation subject number, unique across sites;
    
		do i1 = 1 to nrow(site);
            do i2 = 1 to &nSubj_perSite;
                subj[iOut, 1:2] = site[i1,];
                subj[iOut, 3:5] = a[i1,];
                subj[iOut, 6] = iIn2;  *iSubject;
                subj[iOut, 7:9] = group[iIn] || x1[iIn] || logx2[iIn];
				subj[iOut, 10:12] = bG[iIn,];
				subj[iOut, 13:18] = t(beta);
				%if &ydist=ZIP or &ydist=ZINB %then %do;
					subj[iOut, 21:23] = a_zi[i1,];
					subj[iOut, 24:26] = bG_zi[iIn,];
					subj[iOut, 27:32] = t(beta_zi);
				%end;

                iOut = iOut+1;
                iIn = iIn+1;
                if iIn2 = ncol(iSubj) then  iIn2 =1;
                else iIn2 = iIn2+1;
             end;
        end;

        subj[,19] = beta[1] + subj[,3] + subj[,10] + subj[,7:9]*beta[4:6]; /*eta0*/
		do s=1 to nSubj;
			subj[s,20] = beta[2] + subj[s,4] + subj[s,11]+ (beta[3] + subj[s,5] + subj[s,12])*subj[s,7];
		end;
        subjColNames = {"iSim" "iSite" "a0" "a1" "a2" "iSubject" "group" "x1" "logx2" "b0" "b1" "b2" 
						'beta0' 'beta1' 'beta2' 'beta3' 'beta4' 'beta5' "eta0" "ab1"};
		%if %sysfunc(upcase(&ydist)) =ZIP or %sysfunc(upcase(&ydist))= ZINB %then %do;
			subj[,33] = beta_zi[1] + subj[,21] + subj[,24] + subj[,7:9]*beta_zi[4:6];
			do s=1 to nSubj;
				subj[s,34] = beta_zi[2] + subj[s,22] + subj[s,25]+ (beta_zi[3] + subj[s,23] + subj[s,26])*subj[s,7];
			end;
	        subjColNames = {"iSim" "iSite" "a0" "a1" "a2" "iSubject" "group" "x1" "logx2" "b0" "b1" "b2"
							'beta0' 'beta1' 'beta2' 'beta3' 'beta4' 'beta5' "eta0" "ab1"
							"a0_zi" "a1_zi" "a2_zi" "b0_zi" "b1_zi" "b2_zi" 
							'beta0_zi' 'beta1_zi' 'beta2_zi' 'beta3_zi' 'beta4_zi' 'beta5_zi'
							"eta0_zi" "ab1_zi"};
		%end;
		mattrib subj colname=subjColNames;        
		*output matrix to a sas dataset;
        create subjects from subj [colname=subjColNames];
        append from subj;
        close subjects;
    quit;

	*generate outcome data;
    data &outds;
        set subjects;
        call streaminit(&seed);
        do time = &time_s to &time_e;
            eta = eta0 + ab1*time;
            meanY = exp(eta);
			%if %sysfunc(upcase(&ydist))=POISSON %then %do;
				y = rand('POISSON', meanY);
            %end;
            %else %if %sysfunc(upcase(&ydist))=NEGBIN %then %do;
                p = (&size)/(meanY + &size);
                y=&ymax+1;
				do while (y>&ymax);
					y = rand('NEGBINOMIAL', p, &size);	
				end;	
             %end;
			%else %if %sysfunc(upcase(&ydist)) = ZIP %then %do;
				*- probability of zero inflation;
				eta_zi=eta0_zi + ab1_zi*time;;
			    p0_zi = 1 / (1 + exp(-eta_zi));

				*call ranbin(seed, 1, p0_zi, zero_inflate);
				zero_inflate=rand('BINOMIAL', p0_zi,1 );
				*-generate response;
				if zero_inflate then y=0;
				else do;	
					y = rand('POISSON', meanY);
				end;
                
             %end;
			 %else %if %sysfunc(upcase(&ydist)) = ZINB %then %do;
				*- probability of zero inflation;
				eta_zi=eta0_zi + ab1_zi*time;;
			    p0_zi = 1 / (1 + exp(-eta_zi));

				*call ranbin(seed, 1, p0_zi, zero_inflate);
				zero_inflate=rand('BINOMIAL', p0_zi,1 );

				*-generate response;
				if zero_inflate then y=0;
				else do;	
					size=&size;
					p = (size)/(meanY + size);
					y=&ymax+1;
					do while (y>&ymax);
						y = rand('NEGBINOMIAL', p, size);	
					end;	
				end;
            %end;
            
            %else %abort abend;
			y_zero=(y=0);
			if y<=&ymax then
            output;
        end;
	run;
	proc sort data=&outds; by iSim;run;

	title 'Descriptive statistics of observed zeros in outcome by group at each assessment time and treatment';
	proc means data=&outds n mean;
	class time group;
	var y_zero;
	run;
	title 'Descriptive statistics of outcome by group at each assessment time';
	proc means data=&outds n mean std min q1 median q3 max;
	class time group;
	var y;
	run;

	title 'Distribution of outcome by group at each assess sment time';
	proc freq data=&outds;
	table time*group*y/nocum nopercent norow;
	run;
	title ' ';
%mend;
***************************************************************************
*simulate POISSON data;
*use the values of parameters so that the simulated data are close to real caries data;
				*order: int, time, time*group, group, x1(age_stdz), and logx2(logSurfaceCount);
**************************************************************************;
%simulateLC(seed=5719, nSim=1000, nSubj_perSite=250, ydist=Poisson, 
			sitecov=%str(0.007, 0.0, 0, 0.0, 0, 0), 
			betas=%str(-12, 0.2, -0.1, -0.01, 0.1, 2.95), 
			gcov=%str(1.5,  0.0,0.0, 0.0, 0.0, 0.0), outds=Poi_data);

***********************************************************************
*simulate NB data;
*use the values of parameters so that the simulated data are close to real caries data;
**************************************************************************;
%simulateLC(seed=4191912, nSim=1000, nSubj_perSite=250, ydist=NegBin, size=0.8,
			sitecov=%str(0.007, 0.0, 0, 0.0, 0, 0), 
			betas=%str(-9, 0.2, -0.1, -0.01, 0.2, 2.2), 
			gcov=%str(1.5, 0.0,0.0, 0.0, 0.0, 0.0), 
			outds=NB_data);

************************************************************************************
*simulate ZIP data;
*use the values of parameters so that the simulated data are close to real caries data;
*order: int, time, time*group, group (group), x1(age_stdz), and logx2(logSurfaceCount);
**************************************************************************************;
%simulateLC(seed=20192403, nSim=1000, nSubj_perSite=250, ydist=ZIP,
			sitecov=%str(0.007, 0.0, 0, 0.0, 0, 0), 
			betas=%str(1.8, 0.4, -0.1, -0.01, 0.087, -0.605), 
			gcov=%str(1.6, 0.0,0.0, 0.0, 0.0, 0.0), 
			sitecov_zi=%str(0.007, 0.0, 0, 0.0, 0, 0), 
			betas_zi=%str(20, -2.04, -0.1, -0.01, -0.64, -5.73), 
			gcov_zi=%str(18.73, 0.0,0.0, 0.0, 0.0, 0.0), 
			gzig_cov=%str(-4.17, 0.0,0.0, 0.0, 0.0, 0.0),
			outds=ZIP_data);

*************************************************************************************
*simulate ZINB data;
*use the values of parameters so that the simulated data are close to real caries data;
***************************************************************************************;
%simulateLC(seed=20192503, nSim=1000, nSubj_perSite=250, ydist=ZINB, size=0.5,
			sitecov=%str(0.007, 0.0, 0, 0.0, 0, 0), 
			betas=%str(1.8, 0.4, -0.1, -0.01, 0.087, -0.605), 
			gcov=%str(1.6, 0.0,0.0, 0.0, 0.0, 0.0), 
			sitecov_zi=%str(0.007, 0.0, 0, 0.0, 0, 0), 
			betas_zi=%str(20, -2.04, -0.1, -0.01, -0.64, -5.73), 
			gcov_zi=%str(18.73, 0.0,0.0, 0.0, 0.0, 0.0), 
			gzig_cov=%str(-4.17,0.0,0.0, 0.0, 0.0, 0.0),
			outds=ZINB_data);



