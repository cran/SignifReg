\name{drop1summary}
\alias{drop1summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summaries of add1 or drop1
}
\description{
Offers summaries of information as every single predictor in the scope is added to or removed from the model. \code{max_pvalue} indicates the maximum p-value from the multiple t-tests for each predictor. \code{alpha_cut_off}, \code{Bonferroni}, and \code{FDR} represent whether the model satisfies the significance by correction. \code{alpha_cut_off} means no correction.
}
\usage{
drop1summary(fit, scope, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
	an lm object representing a model. It is an initial model for the variable selection 

}
  \item{scope}{
  The range of models examined in regression. It should be either a data.frame of formula containing predictors. When scope is data.frame, all variables except the response variable in the data.frame are considered for the variable selection. See the examples how they can be used.
}
  \item{alpha}{
Significance level. Default value is 0.05.
}
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
drop1summary(fit3)
drop1summary(fit3, scope=~.-disp-cyl-wt-qsec-cyl, alpha=0.025)
}

