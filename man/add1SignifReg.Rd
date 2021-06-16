\name{add1SignifReg}
\alias{add1SignifReg}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add a predictor to a (generalized) linear regression model using the forward step in the Significance Controlled Variable Selection method
}
\description{
	add1SignifReg adds to the model the predictor, out of the available predictors, which minimizes the criterion (AIC, BIC, r-ajd, PRESS, max p-value) as long as all the p-values of the predictors in the prospective model (including the prospective predictor) are below the chosen correction method (Bonferroni, FDR, None, etc). The function returns the fitted model with the additional predictor if any. A summary table of the prospective models can be printed with \code{print.step = TRUE}. 
	
	\code{max_pvalue} indicates the maximum p-value from the multiple t-tests for each predictor. More specifically, the algorithm computes the prospective models with each predictor included, and all p-values of this prospective model. Then, the predictor selected to be added to the model is the one whose generating model has the smallest p-values, in fact, the minimum of the maximum p-values in each prospective model.
}
\usage{
add1SignifReg(fit, scope, alpha = 0.05, criterion = "p-value", 
   adjust.method = "fdr", override = FALSE, print.step = FALSE)
}
\arguments{
  \item{fit}{
	an lm or glm object representing a linear regression model. 
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
  \item{override}{If \code{override = TRUE}, it returns a new \code{lm} or \code{glm} object that adds a new variable according to criterion even if the new model does not pass the multiple testing p-value correction. 
  	}
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
add1SifnifReg returns an object of the class \code{lm} or \code{glm} for a generalized regression model with the additional component \code{steps.info}, which shows the steps taken during the variable selection and model metrics: Deviance, Resid.Df, Resid.Dev, AIC, BIC, adj.rsq, PRESS, max_pvalue,  max.VIF, and whether it passed the chosen p-value correction.
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
 \code{\link{SignifReg}},  \code{\link{add1summary}},  \code{\link{drop1summary}},  \code{\link{drop1SignifReg}}
}
\examples{
##mtcars data is used as an example.

data(mtcars)

nullmodel = lm(mpg~1, mtcars)
fullmodel = lm(mpg~., mtcars)
scope = list(lower=formula(nullmodel),upper=formula(fullmodel))
fit1 <- lm(mpg~1, data = mtcars)
add1SignifReg(fit1, scope = scope, print.step = TRUE)

fit2 <- lm(mpg~disp+cyl+wt+qsec, mtcars)
add1SignifReg(fit2, scope = scope, criterion="AIC", override="TRUE")
}

