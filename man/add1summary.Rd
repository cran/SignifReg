\name{add1summary}
\alias{add1summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summaries of models when adding a predictor in (generalized) linear models
}
\description{
Offers summaries of prospective models as every available predictor in the scope is added to the model. 
}
\usage{
add1summary(fit, scope, alpha = 0.05, adjust.method = "fdr", sort.by = "p-value")
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
  \item{adjust.method}{
  Correction for multiple testing accumulation of error. See \code{\link{p.adjust}}.
}
  \item{sort.by}{ The criterion to use to sort the table of prospective models. Must be one of \code{criterion = "AIC"}, 
  \code{criterion = "BIC"}, \code{criterion = "r-adj"} (adjusted r-square), \code{criterion = "PRESS"}, and \code{criterion = "p-value"} are available. Default is p-value. }
}
\details{
\code{max_pvalue} indicates the maximum p-value from the multiple t-tests for each predictor.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\value{
a table with the possible inclusions and the metrics of the prospective models: AIC, BIC, adj.rsq, PRESS, max_pvalue,  max.VIF, and whether it passed the chosen p-value correction.
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
 \code{\link{SignifReg}},  \code{\link{add1SignifReg}},  \code{\link{drop1summary}},  \code{\link{drop1SignifReg}}
}
\examples{
##mtcars data is used as an example.
	
data(mtcars)

nullmodel = lm(mpg~1, mtcars)
fullmodel = lm(mpg~., mtcars)
scope = list(lower=formula(nullmodel),upper=formula(fullmodel))
fit1 <- lm(mpg~1, mtcars)
add1summary(fit1, scope = scope)

fit2 <- lm(mpg~disp+cyl+wt+qsec+cyl, data = mtcars)
add1summary(fit2, scope = scope)
}

