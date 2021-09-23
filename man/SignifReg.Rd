\name{SignifReg}
\alias{SignifReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Significance Controlled Variable Selection in (Generalized) Linear Regression
}
\description{
  Significance controlled variable selection selects variables in a generalized linear regression model with different directions of the algorithm (forward, backward, stepwise) based on a chosen criterion (AIC, BIC, adjusted r-square, PRESS or p-value). 
  The algorithm selects a final model with only significant variables based on a correction choice of False Discovery Rate, Bonferroni, etc from the p.adjust().
}
\usage{
SignifReg(fit, scope, alpha = 0.05, direction = "forward",
  criterion = "p-value", adjust.method = "fdr", trace=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
	an lm or glm object representing a model. It is an initial model for the variable selection.
}
  \item{scope}{
  defines the range of models examined in the stepwise search. This should be either a single formula, or a list containing components upper and lower, both formulae. See the details for how to specify the formulae and how they are used.
  }
  \item{alpha}{
  Significance level. Default value is 0.05.
}
  \item{direction}{
  Direction in variable selection: \code{direction = "both"}, 
 
  \code{direction = "forward"},
   and 
   
   \code{direction = "backward"} are available. 
  \code{direction = "both"} is a stepwise selection. Default is \code{direction = "forward"}. 
}
  \item{criterion}{
  Criterion to select predictor variables. \code{criterion = "AIC"}, 
  \code{criterion = "BIC"}, \code{criterion = "r-adj"} (adjusted r-square), \code{criterion = "PRESS"}, and \code{criterion = "p-value"} are available. Default is p-value.
}
  \item{adjust.method}{
  Correction for multiple testing accumulation of error. See \code{\link{p.adjust}}.
}
  \item{trace}{
  	If true, information is printed for each step of variable selection.
  	Default is \code{FALSE}.
  	Offers summaries of prospective models as each predictor in the scope is added to or removed from the model. \code{max_pvalue} indicates the maximum p-value from the multiple t-tests for each predictor in the model. 
  }
}
\details{
SignifReg selects only significant predictors according to a designated criterion. A model with the best criterion, for example, the smallest AIC, will not be considered if it includes insignificant predictors based on the chosen correction. When the criterion is "p-value", a predictor can be droped only if the current model has an insignificant pedictor, and a predictor can be added as long as the prospective model has all predictors significant (including the one to be added). The predictor to be added or removed is the one that generates a model having the smallest maximum p-value of the t-tests in the prospective models. This step is repeated as long as every predictor is significant according to the correction criterion. In the case that the criterion is "AIC", and "BIC", SignifReg selects, at each step, the model having the smallest value of the criterion among models having only significant predictors according to the chosen correction. 
}
\value{
SifnifReg returns an object of the class \code{lm} or \code{glm} for a generalized regression model with the additional component \code{steps.info}, which shows the steps taken during the variable selection and model metrics: Deviance, Resid.Df, Resid.Dev, AIC, BIC, adj.rsq, PRESS, max_pvalue,  max.VIF, and whether it passed the chosen p-value correction.
}
\references{
Zambom A Z, Kim J. Consistent significance controlled variable selection in high-dimensional regression. Stat.2018;7:e210. https://doi.org/10.1002/sta4.210
}
\author{
Jongwook Kim <jongwook226@gmail.com> 

Adriano Zanin Zambom <adriano.zambom@gmail.com>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{add1SignifReg}}, \code{\link{drop1SignifReg}}, \code{\link{add1summary}}, \code{\link{drop1summary}}
}
\examples{
##mtcars data is used as an example.

data(mtcars)
nullmodel = lm(mpg~1, mtcars)
fullmodel = lm(mpg~., mtcars)
scope = list(lower=formula(nullmodel),upper=formula(fullmodel))


fit1 <- lm(mpg~1, mtcars)
select.fit = SignifReg(fit1, scope = scope, direction = "forward", trace = TRUE)
select.fit$steps.info

fit = lm(mpg ~cyl + hp + am + gear, data = mtcars)
select.fit = SignifReg(fit,scope=scope, alpha = 0.05,direction = "backward",
  criterion = "p-value",adjust.method = "fdr",trace=TRUE)
select.fit$steps.info



fit = lm(mpg ~ cyl + hp + am + gear + disp, data = mtcars)
select.fit = SignifReg(fit,scope=scope, alpha = 0.5,direction = "both",
  criterion = "AIC",adjust.method = "fdr",trace=TRUE)
select.fit$steps.info


}

