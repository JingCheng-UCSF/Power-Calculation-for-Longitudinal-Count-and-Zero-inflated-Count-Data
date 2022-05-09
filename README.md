# Power-Calculation-for-Longitudinal-Count-and-Zero-inflated-Count-Data

%lcPower SAS macro to calculate power for longitudinal count and zero-inflated count outcome data
PURPOSE
The %lcPower macro calculates power for longitudinal correlated count and zero-inflated count outcomes in multi-center studies with randomization stratified by study site. For the longitudinal Poisson and Negative Binomial outcome data the macro lcPower fits generalized linear mixed effect models (GLMM). For the longitudinal ZIP and ZINB outcomes the macro lcPower allows users to choose the generalized linear mixed effect models (GLMM) or marginal ZIP/ZINB models (Long et al., 2015, Preisser 2016) to fit the data. Then the rejection rate under the alternative hypothesis of true group difference over time will be used to estimate the power of testing the parameter of treatment effect over time (group*time interaction). Corresponding sample size n will be the estimated sample size to reach the power in detecting the treatment effect over time.

REQUIREMENTS
Base SAS and SAS/STAT software are required.

USAGE
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

EXAMPLES
1.	Calculate power for longitudinal Poisson data, which is generated by using macro %simulateLC (see example 1 in the documentation of %simulateLC)
%lcPower(ds=poi_data, ydist=Poisson, outds=Poi_power);

2.	Calculate power for longitudinal ZINB data, which is generated by using macro %simulateLC (see example 2 in the documentation of %simlateLC) by fitting marginal ZINB model
%lcPower(ds=zinb_data, ydist=zinb, gconv=0.0000001, mzi=1, outds=mzinb_power);



REFERENCES
Long DL, Preisser JS, Herring AH, Golin CE. A Marginalized Zero-Inflated Poisson Regression Model with Random Effects. Journal of the Royal Statistical Society, Series C (Applied Statistics). 2015 November. DOI: 10.1111/rssc.12104.

Preisser JS, Das K, Long DL, Divaris K. Marginalized zero-inflated negative binomial regression with application to dental caries. Stat Med. 2016 May; 35(10): 1722–1735. doi:10.1002/sim.6804.


CITATION
Cheng NF, Cheng J: A SAS macro to calculate power for longitudinal count and zero-inflated count outcome data.

CONTACTS
Nancy.Cheng@ucsf.edu
Jing.Cheng@ucsf.edu
