SMNlmec R package.

The use of mixed-effect models to understand the evolution of the human immunodeficiency virus (HIV) and the progression of acquired immune deficiency syndrome (AIDS) has been the cornerstone of longitudinal data analysis in recent years. However, data from HIV/AIDS clinical trials have several complexities. Some of the most common recurrences are related to the situation where the HIV viral load can be undetectable, and the measures of the patient can be registered irregularly due to some problems in the data collection. Although censored mixed-effects models assuming conditionally independent normal random errors are commonly used to analyze this data type, this model may need to be more appropriate for accommodating outlying observations and responses recorded at irregular intervals.
Consequently, in this paper, we propose a Bayesian analysis of censored linear mixed-effects models that replace Gaussian assumptions with a flexible class of distributions, such as the scale mixture of normal family distributions, considering a damped exponential correlation structure which was employed to account for within-subject autocorrelation among irregularly observed measures. For this complex structure, Stan's default No-U-Turn sampler is utilized to obtain posterior simulations. The feasibility of the proposed methods was demonstrated through several simulation studies and their application to two AIDS case studies.

To install this package from this repo, please input the following command in R.

```r
library(devtools)

devtools::install_github("KelinZhong/SMNlmec")
```
