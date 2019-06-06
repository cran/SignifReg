\name{SignifReg}
\alias{SignifReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Consistent Significance Controlled Variable Selection in Linear Regression
}
\description{
  Significance controlled variable selection selects variables in a linear regression model with different directions of th algorithm (forward, backward, stepwise) based on a chosen criterion (AIC, BIC, adjusted r-square, and p-value). 
  The algorithm selects a final model with only significant variables based on a correction choice of False Discovery Rate, Bonferroni, or fixed level.
}
\usage{
SignifReg(fit, scope, alpha = 0.05, direction = "forward", criterion = "p-value",
 correction = "FDR", trace=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
	an lm object representing a model. It is an initial model for the variable selection 
}
  \item{scope}{
  The range of models examined in regression. It should be either a data.frame or formula containing predictors. When scope is data.frame, all variables except the response variable in the data.frame are considered for the variable selection.
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
  \code{criterion = "BIC"}, \code{criterion = "r-adj"} (adjusted r-square), and \code{criterion = "p-value"} are available. Default is p-value.
}
  \item{correction}{
  Correction for multiple testing accumulation of error. \code{correction = "FDR"} (False Discovery Rate), \code{correction = "Bonferroni"}, 
  and 
  
  \code{correction = "None"} (no correction) are available. Default is 
  
  \code{correction = "FDR"} . For Bonferroni correction, 
  
  either \code{correction = "Bonferroni"}   or \code{correction = "Bonf"} can be used.
}
  \item{trace}{
  	If true, information is printed for each step of variable selection.
  	Default is \code{FALSE}.
  	Offers summaries of information as every single predictor in the scope is added to or removed from the model. \code{max_pvalue} indicates the maximum p-value from the multiple t-tests for each predictor in the model. This value can be used as a criterion in the case \code{criterion = "p-value"}. \code{alpha_cut_off}, \code{Bonferroni}, and \code{FDR} represent whether the model satisfies the significance by correction. \code{alpha_cut_off} means no correction and all p-values will be compared to it.
  }
}
\details{
SignifReg selects only significant predictors according to a designated creterion. Although a model has the best feature for a chosen criterion, for example, the smallest AIC, this model will be excluded from the selection if it includes insignificant predictors based on the chosen correction. When the criterion is "p-value", at each step, the variable to be added or substracted is the one that generates a model having the smallest maximum p-value of the t-tests. This step is repeated as long as every predictor is significant according to the correction criterion. In the case that the criterion is "AIC", and "BIC", SignifReg selects, at each step, the model having the smallest value of the criterion among models having only significant predictors according to the chosen correction. For "r-adj(adjusted r-square), SignifReg selects, at each step, the model having the largest value of adjusted r-square among models having only significant predictors according to the chosen correction.
}
\value{
SifnifReg returns an object of the class \code{lm} for a regression model.
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##mtcars data is used as an example.

data(mtcars)

fit1 <- lm(mpg~1, mtcars)
summary(SignifReg(fit1))
summary(SignifReg(fit1, scope=mtcars ,criterion="r-adj", direction="forward" 
,correction="None", trace=TRUE))
summary(SignifReg(fit1, scope=~.+disp+cyl+wt+qsec, criterion="p-value", direction="both",
correction="FDR", trace=TRUE))

fit2 <- lm(mpg~., mtcars)
summary(SignifReg(fit2, scope=~.-disp-cyl-wt-hp-drat-qset-carb , direction="backward", 
correction="Bonf", trace=TRUE))
summary(SignifReg(fit2, scope=mtcars ,criterion="AIC", direction="both",
correction="None", trace=TRUE))
}

