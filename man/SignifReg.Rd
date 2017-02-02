\name{SignifReg}
\alias{SignifReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Significant Variable Selection in Linear Regression
}
\description{
  Provide a significant variable selection procedure with different directions (forward, backward, stepwise) based on diverse criteria (Mallows' Cp, AIC, BIC, adjusted r-square, p-value). 
  The algorithm selects a final model with only significant variables based on a correction choice of False Discovery Rate, Bonferroni, or no correction.
}
\usage{
SignifReg(scope, data, alpha = 0.05, direction = "forward", 
criterion = "p-value", correction = "FDR")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scope}{
  The range of models examined in regression. It should be a formula containing a response variable and predictors.
  For the full model, \code{<response var>~.} can be used. 
  }
  \item{data}{
  Data frame, examined in regression.
}
  \item{alpha}{
  Significance level. Default value is 0.05.
}
  \item{direction}{
  Direction in variable selection: \code{direction = "step_null"}, 
 
  \code{direction = "step_full"}, \code{direction = "forward"},
   and 
   
   \code{direction = "backward"} are available. 
  \code{direction = "step_null"} is a stepwise selection starting from a null model. 
  \code{direction = "step_full"} is a stepwise selection starting from a full model. Default is \code{direction = "forward"}. 
  When criterion is p-value, \code{direction = "step_null"} and 
  
  \code{direction = "step_full"} are unavailable.
}
  \item{criterion}{
  Criterion to select predictor variables. \code{criterion = "AIC"}, 
  
  \code{criterion = "BIC"}, \code{criterion = "Cp"}(Mallows'Cp), \code{criterion = "r-adj"} (adjusted r-square), and \code{criterion = "p-value"} are available. Default is p-value.
}
  \item{correction}{
  Correction criterion to reduce multiple testing error. \code{correction = "FDR"} (False Discovery Rate), \code{correction = "Bonferroni"}, 
  and 
  
  \code{correction = "None"} (no correction) are available. Default is 
  
  \code{correction = "FDR"} . For Bonferroni correction, 
  
  either \code{correction = "Bonferroni"}   or \code{correction = "Bonf"} can be used.
}
}
\details{
SignifReg selects only significant predictors according to a designated correction. Although a model has the best feature for a chosen criterion, for example, the smallest AIC, this model will be excluded from the selection if it includes insignificant predictors based on the correction criterium. When the criterion is "p-value" and direction is "forward", at each step the variable to be added is the one that generates a model having the smallest maximum p-value of all variables included. This step is repeated as long as every predictor is significant according to the correction criterium. When the criterion is "p-value" and the direction is "backward", the process excludes a predictor having the biggest p-value one by one until the model includes only significant predictors according to the correction criterium. In the case that the criterion is "AIC", "BIC", or "Cp"(Mallows'Cp), SignifReg selects, at each step, the model having the smallest value of the criterion among models having only significant predictors according to the correction criterium. For "r-adj(adjusted r-square), SignifReg selects, at each step, the model having the largest value of adjusted r-square among models having only significant predictors according to the correction criterium."
}
\value{
SifnifReg returns an object of the class \code{lm} for a regression model.
}

\author{
Jongwook Kim <jongwook226@gmail.com> 

and Adriano Zanin Zambom <adriano.zambom@gmail.com>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
##mtcars data is used as an example.

data(mtcars)
scope1 <- mpg~.
model1 <- SignifReg(scope=scope1,data=mtcars,alpha=0.05,
direction="forward", criterion="p-value",
correction="FDR")

model2 <-SignifReg(scope=scope1,data=mtcars,alpha=0.05,
direction="step_full", criterion="Cp",
correction="Bonferroni")

scope2 <- mpg~cyl+disp+hp+drat+wt+qsec
model3 <-SignifReg(scope=scope2,data=mtcars,alpha=0.05,
direction="step_null", criterion="AIC",
correction="None")

summary(model1)
summary(model2)
summary(model3)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
