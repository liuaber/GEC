### Explaining why not using Mathematica, instead use R in this part.
require("expint")
library(expint)
f_gamma = function(alpha){
  gammainc(alpha+1,alpha)/gamma(alpha)+0.5*alpha^2-alpha }
uniroot(f_gamma,interval = c(0.001,10))

f_k = function(k){
  alpha = 0.308289
  k^alpha*exp(-k)/gamma(alpha)-alpha+k 
}
uniroot(f_k,interval = c(0.001,10))
library(invgamma)
1/(gammainc(0.3082908,0.144351)/gamma(alpha)+1)

   