\name{How to use}
\alias{modeldiagnostic}
\title{GLM Diagnostic}

\description{This package generates the following plots for Generalized Linear Model diagnostics: Q-Q plot, Residuals vs Fitted, Histogram of Residuals and Cook's distance plot.}

\usage{
modeldiagnostic1(x)
modeldiagnostic2(x)
}

\arguments{
\item{x}{GLM model.}
}

\examples{

varY<-c(0,10,20,30,40,50,20,30,43,28)

varX1<-c(0,11,22,31,44,50,24,34,12,45)

varX2<-c("A","B","A","B","C","C","A","B","C","A")

data<-data.frame(varY,varX1,varX2)

fit <- glm(
varY ~ varX1 * varX2,
data = data
)

modeldiagnostic1(fit)
modeldiagnostic2(fit)

}
