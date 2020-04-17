\name{drop1SignifReg}
\alias{drop1SignifReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Drop a predictor to a linear regression model using the backward step in the Significance Controlled Variable Selection method
}
\description{
	drop1SignifReg removes from the model the predictor, out of the current predictors, which minimizes the criterion (AIC, BIC, r-ajd, PRESS, max p-value) when a) the p-values of the predictors in the current model do not pass the multiple testing correction (Bonferroni, FDR, None) or b) when the p-values of both current and prospective models pass the correction but the criterion of the prospective model is smaller. 
	
	\code{max_pvalue} indicates the maximum p-value from the multiple t-tests for each predictor. More specifically, the algorithm computes the prospective models with each predictor removed, and all p-values of this prospective model. Then, the predictor to be removed is the one whose generating model has the smallest p-values, in fact, the minimum of the maximum p-values in each prospective model.
}
\usage{
drop1SignifReg(fit, scope, alpha = 0.05, criterion = "p-value",
 correction = "FDR", override = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
	an lm object representing a model.
}
  \item{scope}{
  The range of models examined in regression. It should be either a data.frame or formula containing predictors. When scope is data.frame, all variables except the response variable in the data.frame are considered for the variable selection. 
}
  \item{alpha}{
  Significance level. Default value is 0.05.
}
  \item{criterion}{
  Criterion to select predictor variables. \code{criterion = "AIC"}, 
  \code{criterion = "BIC"}, \code{criterion = "r-adj"} (adjusted r-square), and \code{criterion = "p-value"} are available. Default is p-value.
}
  \item{correction}{
  Correction criterion to reduce multiple testing error. \code{correction = "FDR"} (False Discovery Rate), \code{correction = "Bonferroni"}, 
  and 
  
  \code{correction = "None"} (no correction) are available. Default is 
  
  \code{correction = "FDR"} . For Bonferroni correction, 
  
  either \code{correction = "Bonferroni"}   or \code{correction = "Bonf"} can be used.
}
  \item{override}{If \code{override = TRUE}, it returns a new lm object that excludes a predictor according to criterion even if the prospective model passes the correction. }
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

fit3 <- lm(mpg~., mtcars)
drop1SignifReg(fit3)
drop1SignifReg(fit3, scope=~.-disp-cyl-wt-hp-drat-qset-carb ,correction="Bonf")
}

