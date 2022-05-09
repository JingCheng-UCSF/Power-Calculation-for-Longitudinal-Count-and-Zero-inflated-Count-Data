# Power-Calculation-for-Longitudinal-Count-and-Zero-inflated-Count-Data

A SAS macro lcPower has been developed to calculate power for longitudinal correlated count and zero-inflated count outcomes, which commonly occur in multi-center health science studies. The macro models longitudinal outcomes which follow Poisson, Negative Binomial, Zero-inflated Poisson (ZIP) or Zero-inflated Negative Binomial (ZINB) distributions. Such longitudinal outcomes have brought challenges to power analysis. A separate SAS macro simulateLC we developed allows users to generate longitudinal count and zero-inflated count data expected to be observed in their multi-center longitudinal studies, and then users can use the macro lcPower to choose generalized linear mixed models or marginal ZIP/ZINB models ( Long et al., 2015, Preisser 2016) to fit the simulated data for the power calculation. For the longitudinal Poisson and Negative Binomial outcomes, this macro lcPower fits generalized linear mixed effect models. For the longitudinal ZIP and ZINB outcomes, the macro lcPower allows users to choose the generalized linear mixed models or marginal ZIP/ZINB models (Long et al., 2015, Preisser 2016) to fit the data. Then the rejection rate under the alternative hypothesis of true group difference will be used to estimate the power of testing the parameter of treatment effect. Corresponding sample size n will be the estimated sample size to reach the power in testing the parameter of treatment effect. The lcPower macro provides a helpful tool for researchers to estimate the power and sample size for longitudinal count and zero-inflated count data in their studies. 

*************************************************************************************************************************************************************
%simulateLC SAS macro to simulate longitudinal correlated count and zero-inflated count outcome data

PURPOSE:

The %simulateLC macro creates a numeric time variable (time), a binary treatment variable (group), a normal covariate (x1) and log of Poisson covariate (logx2), and then generates longitudinal correlated count and zero-inflated count outcome data.

REQUIREMENTS:

Base SAS and SAS/IML software are required.

USAGE:

Save the %simulateLC macro definition to your computer system. In your SAS program, add a %inc statement to specify the physical name of an external file where the %simulateLC macro is stored and make it available for use: 
	%inc "<location of simulateLC macro>";
	
Example:  %inc 'c:\mysas\simulateLC.sas';
Following this statement, invoke the %simulateLC macro.
	
The %simulateLC macro has below parameters: 

seed=number
specifies a random seed. 

nSim=number >0
specifies the number of simulations.

nSite=number ≥2
specifies the number of study site. The default is 2.

nSubj_perSite=number >0
specifies the number of subjects within a study site. The default is 100.

ydist=’text-string’
specifies the outcome distribution, which can be “Poisson”, “Negbin”, “ZIP”, or “ZINB”.

size=number
          specifies the dispersion parameter for negative binomial outcome distribution. 
The default is 1.

ymax=number >0
specifies the maximum value that the outcome is allowed. The default is 88.

x2m=number >0
specifies the Poisson mean to generate the second covariate x2. The default is 87.

time_s=number ≥0
specifies the first study timepoint to be included in the outcome evaluation. The default is 0.

time_e=number >0
specifies the last study timepoint to be included in the outcome evaluation. The default is 3.

sitecov=%str(number1,number2, number3, number4, number5, number6)
specifies a string with 6 numeric values separated by commas, eg, %str(1.5,0,0,0,0,0) as the upper triangular covariance matrix for site random effects in the Poisson or NB outcome model or in the count part of ZIP/ZINB outcome model, in the order of variance of random intercept for sites, covariance between random intercept and random slope over time for sites, covariance between random intercept and random time*group interaction for sites, variance of random slope over time for sites, covariance between random slope over time and random time*group interaction for sites, and variance of random time*group interaction for sites.

betas=%str(number1,number2, number3, number4, number5, number6) 
  specifies a string of 6 numeric values, separated by commas for the fixed effect coefficients in the Poisson or NB outcome model or in the count part of ZIP/ZINB outcome model. eg, %str(-12, 0.2, -0.1, -0.01, 0.1, 2.95), with the fixed effects in the order of intercept, time, time*group, group, x1, and logx2.

gov=%str(number1,number2, number3, number4, number5, number6) 
  specifies a string of 6 numeric values separated by commas, eg, %str(1.8,0.2,0,0,0,0)
as the upper triangular g-side covariance matrix in the Poisson or NB outcome model or in the count part of ZIP/ZINB outcome model, in the order of variance of random intercept for individuals, covariance between random intercept and random slope over time for individuals, covariance between random intercept and random time*group interaction for individuals, variance of random slope over time for individuals, covariance between random slope over time and random time*group interaction for individuals, and variance of random time*group interaction for individuals.

sitecov_zi=%str(number1,number2, number3, number4, number5, number6)
specifies a string with 6 numeric values separated by commas, eg, %str(1.5,0,0,0,0,0) as the upper triangular covariance matrix for site random effects in the logit part of ZIP/ZINB outcome model, in the order of variance of random intercept for sites, covariance between random intercept and random slope over time for sites, covariance between random intercept and random time*group interaction for sites, variance of random slope over time for sites, covariance between random slope over time and random time*group interaction for sites, and variance of random time*group interaction for sites. Required if ydist=ZIP/ZINB.

betas_zi=%str(number1,number2, number3, number4, number5, number6) 
  specifies a string of 6 numeric values, separated by commas for the fixed effect coefficients in the logit part of ZIP/ZINB outcome model. eg, %str(-12, 0.2, -0.1, -0.01, 0.1, 2.95), with the fixed effects in the order of intercept, time, time*group, group, x1, and logx2. Required if ydist=ZIP/ZINB.

gov_zi=%str(number1,number2, number3, number4, number5, number6) 
  specifies a string of 6 numeric values separated by commas, eg, %str(1.8,0.2,0,0,0,0)
as the upper triangular g-side covariance matrix in the logit part of ZIP/ZINB outcome model, in the order of variance of random intercept for individuals, covariance between random intercept and random slope over time for individuals, covariance between random intercept and random time*group interaction for individuals, variance of random slope over time for individuals, covariance between random slope over time and random time*group interaction for individuals, and variance of random time*group interaction for individuals.
Required if ydist=ZIP/ZINB.

gzig_cov=%str(number1,number2, number3, number4, number5, number6) 
specifies a string of 6 numeric values separated by commas, eg, %str(-4.17, 0, 0, 0, 0, 0), as the upper triangular covariance matrix between random effects in the logit part and random effects in the count part of ZIP/ZINB outcome model, in the order of covariance between random intercepts in the logit and count parts for individuals, covariance between random intercept and random slope over time in the logit and count parts for individuals, covariance between random intercept and random time*group interaction in the logit and count parts for individuals, covariance of random slope over time in the logit and count parts for individuals, covariance between random slope over time and random time*group interaction in the logit and count parts for individuals, and covariance between random time*group interaction in the logit and count parts for individuals.
 Required if ydist=ZIP/ZINB

outds=output-SAS-data-set 
specifies a SAS dataset to save the simulated data.


EXAMPLES:
	
1.	 Generate 1000 simulations of Poisson data with 2 study sites over 3 timepoints (0, 1, and 2 time points), and 500 subjects for each study site.
 %simulateLC(seed=70439, nSim=1000, nSite=2, nSubj_perSite=500, ydist=poisson,
                      ymax=88, time_s=0, time_e=2, sitecov=%str(0.007, 0.0, 0, 0.0, 0, 0), 
                      betas=%str(-12, 0.2, -0.1, -0.01, 0.1, 2.95), gcov=%str(1.5,  0.0, 0.0, 0.0, 0.0, 0.0),    
                      outds=Poi_data)
2.	Generate 1000 simulations of ZINB data with 2 study sites over 4 timepoints (0, 1, 2, and 3 time points), and 250 subjects for each study site 
%simulateLC(seed=26417, nSim=1000, nSite=2, nSubj_perSite=250, ydist=zinb, size=0.5,
 	         ymax=88, time_s=0, time_e=3, sitecov=%str(0.007, 0.0, 0, 0.0, 0, 0), 
                      betas=%str(1.8, 0.4, 0, -0.01, 0.087, -0.605), gcov=%str(1.5,  0.0, 0.0, 0.0, 0.0, 0.0),
                     sitecov_zi=%str(0.007, 0.0, 0, 0.0, 0, 0), betas_zi=%str(20, -2.04, 0, -0.01, -0.64, -5.73),      
                     gcov_zi=%str(18.73,  0.0, 0.0, 0.0, 0.0, 0.0), 
                     gzig_cov=%str(-4.17, 0.0,0.0, 0.0, 0.0, 0.0), outds=zinb_data)

CITATION:
	
Cheng NF, Cheng J: A SAS macro to simulate longitudinal count and zero-inflated count data

CONTACTS:
	
Nancy.Cheng@ucsf.edu
	
Jing.Cheng@ucsf.edu



*************************************************************************************************************************************************************
%lcPower SAS macro to calculate power for longitudinal count and zero-inflated count outcome data

PURPOSE:
The %lcPower macro calculates power for longitudinal correlated count and zero-inflated count outcomes in multi-center studies with randomization stratified by study site. For the longitudinal Poisson and Negative Binomial outcome data the macro lcPower fits generalized linear mixed effect models (GLMM). For the longitudinal ZIP and ZINB outcomes the macro lcPower allows users to choose the generalized linear mixed effect models (GLMM) or marginal ZIP/ZINB models (Long et al., 2015, Preisser 2016) to fit the data. Then the rejection rate under the alternative hypothesis of true group difference over time will be used to estimate the power of testing the parameter of treatment effect over time (group*time interaction). Corresponding sample size n will be the estimated sample size to reach the power in detecting the treatment effect over time.

REQUIREMENTS:
	
Base SAS and SAS/STAT software are required.

USAGE:
	
Save the %lcPower macro definition to your computer system. In your SAS program, add a %inc statement to specify the physical name of an external file where the %lcPower macro is stored and make it available for use: 
	%inc "<location of lcPower macro>";
	
Example:  %inc 'c:\mysas\lcPower.sas';
Following this statement, invoke the %lcPower macro.
	
The %lcPower macro has parameters below: 

ds=input-SAS-data-set
specifies the input data set generated using macro %simulateLc. 

ydist=’text-string’
specifies the outcome distribution, which can be “Poisson”, “Negbin”, “ZIP”, or “ZINB”.

mzi=number 0 or 1
specifies the GLMM or Marginal ZI model for zip or zinb distributions, mzi=1 for Marginal ZI model and mzi=0 for GLMM. mzi=0 if ydist = Poisson or Negbin. The default is 0.

gconv=number
specifies the relative gradient convergence criterion. Required if ydist=ZIP/ZINB. The default is 1E-8.

outds=output-SAS-data-set 
specifies the name of output SAS dataset to save the power(s).

EXAMPLES:
	
1.	Calculate power for longitudinal Poisson data, which is generated by using macro %simulateLC (see example 1 in the documentation of %simulateLC)
%lcPower(ds=poi_data, ydist=Poisson, outds=Poi_power);

2.	Calculate power for longitudinal ZINB data, which is generated by using macro %simulateLC (see example 2 in the documentation of %simlateLC) by fitting marginal ZINB model
%lcPower(ds=zinb_data, ydist=zinb, gconv=0.0000001, mzi=1, outds=mzinb_power);



REFERENCES:
	
Long DL, Preisser JS, Herring AH, Golin CE. A Marginalized Zero-Inflated Poisson Regression Model with Random Effects. Journal of the Royal Statistical Society, Series C (Applied Statistics). 2015 November. DOI: 10.1111/rssc.12104.

Preisser JS, Das K, Long DL, Divaris K. Marginalized zero-inflated negative binomial regression with application to dental caries. Stat Med. 2016 May; 35(10): 1722–1735. doi:10.1002/sim.6804.


CITATION:
	
Cheng NF, Cheng J: A SAS macro to calculate power for longitudinal count and zero-inflated count outcome data.

CONTACTS:
	
Nancy.Cheng@ucsf.edu 
	
Jing.Cheng@ucsf.edu
