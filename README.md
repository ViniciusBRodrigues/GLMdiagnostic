# GLMdiagnostic (biodata package)
This package generates the following plots for Generalized Linear Model diagnostics: Q-Q plot, Residuals vs Fitted, Histogram of Residuals and Cook's distance plot.

`Based on functions from RT4BIO package (Original authors: Ronaldo Reis Junior, Gladson Ramon Alves Borges and Magnel Lima de Oliveira. Link: https://www.researchgate.net/publication/282808626_RT4Bio_-_R_Tools_for_Biologists) and sjPlot package (Original authors: Daniel Lüdecke and Carsten Schwemmer).`


## Usage

`modeldiagnostic1(x)`

`modeldiagnostic2(x)`

## Examples

`varY<-c(0,10,20,30,40,50,20,30,43,28)`

`varX1<-c(0,11,22,31,44,50,24,34,12,45)`

`varX2<-c("A","B","A","B","C","C","A","B","C","A")`

`data<-data.frame(varY,varX1,varX2)`

`fit <- glm(
varY ~ varX1 * varX2,
data = data
)`

`modeldiagnostic1(fit)`

`modeldiagnostic2(fit)`

## Plots

![modeldiagnostic1](https://s13.postimg.org/854j9fshz/Rplot.png)

![modeldiagnostic1](https://s13.postimg.org/tg6oqx4zr/Rplot01.png)
