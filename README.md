# Power-Calculation-for-Longitudinal-Count-and-Zero-inflated-Count-Data

A SAS macro lcPower has been developed to calculate power for longitudinal correlated count and zero-inflated count outcomes, which commonly occur in multi-center health science studies. The macro models longitudinal outcomes which follow Poisson, Negative Binomial, Zero-inflated Poisson (ZIP) or Zero-inflated Negative Binomial (ZINB) distributions. Such longitudinal outcomes have brought challenges to power analysis. A separate SAS macro simulateLC we developed allows users to generate longitudinal count and zero-inflated count data expected to be observed in their multi-center longitudinal studies, and then users can use the macro lcPower to choose generalized linear mixed models or marginal ZIP/ZINB models ( Long et al., 2015, Preisser 2016) to fit the simulated data for the power calculation. For the longitudinal Poisson and Negative Binomial outcomes, this macro lcPower fits generalized linear mixed effect models. For the longitudinal ZIP and ZINB outcomes, the macro lcPower allows users to choose the generalized linear mixed models or marginal ZIP/ZINB models (Long et al., 2015, Preisser 2016) to fit the data. Then the rejection rate under the alternative hypothesis of true group difference will be used to estimate the power of testing the parameter of treatment effect. Corresponding sample size n will be the estimated sample size to reach the power in testing the parameter of treatment effect. The lcPower macro provides a helpful tool for researchers to estimate the power and sample size for longitudinal count and zero-inflated count data in their studies. 

%simulateLC SAS macro to simulate longitudinal correlated count and zero-inflated count outcome data

PURPOSE:
The %simulateLC macro creates a numeric time variable (time), a binary treatment variable (group), a normal covariate (x1) and log of Poisson covariate (logx2), and then generates longitudinal correlated count and zero-inflated count outcome data.


%lcPower SAS macro to calculate power for longitudinal count and zero-inflated count outcome data

PURPOSE:
The %lcPower macro calculates power for longitudinal correlated count and zero-inflated count outcomes in multi-center studies with randomization stratified by study site. For the longitudinal Poisson and Negative Binomial outcome data the macro lcPower fits generalized linear mixed effect models (GLMM). For the longitudinal ZIP and ZINB outcomes the macro lcPower allows users to choose the generalized linear mixed effect models (GLMM) or marginal ZIP/ZINB models (Long et al., 2015, Preisser 2016) to fit the data. Then the rejection rate under the alternative hypothesis of true group difference over time will be used to estimate the power of testing the parameter of treatment effect over time (group*time interaction). Corresponding sample size n will be the estimated sample size to reach the power in detecting the treatment effect over time.

REQUIREMENTS:

Base SAS and SAS/STAT software are required.


REFERENCES:
	
Long DL, Preisser JS, Herring AH, Golin CE. A Marginalized Zero-Inflated Poisson Regression Model with Random Effects. Journal of the Royal Statistical Society, Series C (Applied Statistics). 2015 November. DOI: 10.1111/rssc.12104.

Preisser JS, Das K, Long DL, Divaris K. Marginalized zero-inflated negative binomial regression with application to dental caries. Stat Med. 2016 May; 35(10): 1722–1735. doi:10.1002/sim.6804.


CITATION:
	
Cheng NF, Cheng J: A SAS macro to calculate power for longitudinal count and zero-inflated count outcome data.

CONTACTS:
	
Nancy.Cheng@ucsf.edu and Jing.Cheng@ucsf.edu
