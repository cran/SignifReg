\name{drop1SignifReg}
\alias{drop1SignifReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Drop a possible single predictor to a regression model based on Consistent Significance Controlled Variable Selection
}
\description{
Drop a single term to the model using the method Consistent Significance Controlled Variable Selection.
The function also provides a summary table of the selection. \code{max_pvalue} indicates the maximum p-value from the multiple t-tests for each predictor. This value can be used as a criterion in the case \code{criterion = "p-value"}. \code{alpha_cut_off}, \code{Bonferroni}, and \code{FDR} represent whether the entire model satisfies the significance correction. \code{alpha_cut_off} means no correction.
 
}
\usage{
drop1SignifReg(fit, scope, alpha = 0.05, criterion = "p-value",
 correction = "FDR", override = FALSE)
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
  \item{override}{If \code{override = TRUE}, it returns a new lm object that excludes a new variable according to criterion even if the current model passes the correction. }
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##mtcars data is used as an example.

data(mtcars)

fit3 <- lm(mpg~., mtcars)
drop1SignifReg(fit3)
drop1SignifReg(fit3, scope=~.-disp-cyl-wt-hp-drat-qset-carb ,correction="Bonf")
}

