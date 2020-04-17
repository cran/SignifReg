\name{add1SignifReg}
\alias{add1SignifReg}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add a predictor to a linear regression model using the forward step in the Significance Controlled Variable Selection method
}
\description{
	add1SignifReg includes in the model the predictor, out of the available predictors, which minimizes the criterion (AIC, BIC, r-ajd, PRESS, max p-value) as long as all the p-values of the predictors in the prospective model (including the prospective predictor) are below the chosen correction method (Bonferroni, FDR, None). The function provides a summary table of the prospective models. 
	
	\code{max_pvalue} indicates the maximum p-value from the multiple t-tests for each predictor. More specifically, the algorithm computes the prospective models with each predictor included, and all p-values of this prospective model. Then, the predictor selected is the one whose generating model has the smallest p-values, in fact, the minimum of the maximum p-values in each prospective model.
}
\usage{
add1SignifReg(fit, scope, alpha = 0.05, criterion = "p-value",
 correction = "FDR", override = FALSE)
}
\arguments{
  \item{fit}{
	an lm object representing a linear regression model. 
}
  \item{scope}{
  The range of models examined in regression. It should be either a data.frame or formula containing predictors. When scope is data.frame, all variables except the response variable in the data.frame are considered.
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
  \item{override}{If \code{override = TRUE}, it returns a new lm object that adds a new variable according to criterion even if the new model does not pass the multiple testing p-value correction. }
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
 \code{\link{SignifReg}},  \code{\link{add1summary}},  \code{\link{drop1summary}},  \code{\link{drop1SignifReg}}
}
\examples{
##mtcars data is used as an example.

data(mtcars)

fit1 <- lm(mpg~1, mtcars)
add1SignifReg(fit1)

fit2 <- lm(mpg~disp+cyl+wt+qsec, mtcars)
add1SignifReg(fit2, criterion="AIC", override="TRUE")
}

