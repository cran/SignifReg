\name{drop1SignifReg}
\alias{drop1SignifReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Drop a predictor to a (generalized) linear regression model using the backward step in the Significance Controlled Variable Selection method
}
\description{
	drop1SignifReg removes from the model the predictor, out of the current predictors, which minimizes the criterion (AIC, BIC, r-ajd, PRESS, max p-value) when a) the p-values of the predictors in the current model do not pass the multiple testing correction (Bonferroni, FDR, None, etc) or b) when the p-values of both current and prospective models pass the correction but the criterion of the prospective model is smaller. 
	
	\code{max_pvalue} indicates the maximum p-value from the multiple t-tests for each predictor. More specifically, the algorithm computes the prospective models with each predictor included, and all p-values of this prospective model. Then, the predictor selected to be added to the model is the one whose generating model has the smallest p-values, in fact, the minimum of the maximum p-values in each prospective model.
}
\usage{
drop1SignifReg(fit, scope, alpha = 0.05, criterion = "p-value", 
   adjust.method = "fdr", override = FALSE, print.step = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
	an lm or glm object representing a model.
}
  \item{scope}{
  defines the range of models examined in the stepwise search. This should be either a single formula, or a list containing components upper and lower, both formulae. See the details for how to specify the formulae and how they are used.
}
  \item{alpha}{
  Significance level. Default value is 0.05.
}
  \item{criterion}{
  Criterion to select predictor variables. \code{criterion = "AIC"}, 
  \code{criterion = "BIC"}, \code{criterion = "r-adj"} (adjusted r-square), \code{criterion = "PRESS"}, and \code{criterion = "p-value"} are available. Default is p-value.
}
  \item{adjust.method}{
  Correction for multiple testing accumulation of error. See \code{\link{p.adjust}}.
}
  \item{override}{If \code{override = TRUE}, it returns a new \code{lm} or \code{glm} object that adds a new variable according to criterion even if the new model does not pass the multiple testing p-value correction.  }
  \item{print.step}{If true, information is printed for each step of variable selection.
  	Default is \code{FALSE}. }
  
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\value{
drop1SifnifReg returns an object of the class \code{lm} or \code{glm} for a generalized regression model with the additional component \code{steps.info}, which shows the steps taken during the variable selection and model metrics: Deviance, Resid.Df, Resid.Dev, AIC, BIC, adj.rsq, PRESS, max_pvalue,  max.VIF, and whether it passed the chosen p-value correction.
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
 \code{\link{SignifReg}},  \code{\link{add1summary}},  \code{\link{add1SignifReg}}, \code{\link{drop1summary}},  
}
\examples{
##mtcars data is used as an example.

data(mtcars)

fit <- lm(mpg~., mtcars)
drop1SignifReg(fit, print.step = TRUE)


}

