### Explaining why not using Mathematica, instead use R in this part.
require("expint")
library(expint)
f_gamma = function(alpha){
  gammainc(alpha+1,alpha)/gamma(alpha)+0.5*alpha^2-alpha }
uniroot(f_gamma,interval = c(0.001,10))

