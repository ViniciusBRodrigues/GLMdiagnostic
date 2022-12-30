# GLMdiagnostic (biodata package)
This package generates plots for Generalized Linear Model diagnostics, such as Q-Q plot, Residuals vs Fitted, Histogram of Residuals and Cook's distance plot.

## Usage

The two main functions are `modeldiagnostic1(fit)` (Q-Q plot, Residuals vs Fitted) and `modeldiagnostic2(fit)` (Histogram of Residuals and Cook's distance plot).

## Examples

`varY<-c(0,10,20,30,40,50,20,30,43,28)
varX1<-c(0,11,22,31,44,50,24,34,12,45)
varX2<-c("A","B","A","B","C","C","A","B","C","A")

data<-data.frame(varY,varX1,varX2)

fit <- glm(
varY ~ varX1 * varX2,
data = data
)

modeldiagnostic1(fit)

modeldiagnostic2(fit)`

---

This package is a fork from [RT4BIO](https://www.researchgate.net/publication/282808626_RT4Bio_-_R_Tools_for_Biologists) package (Ronaldo Reis Junior, Gladson Ramon Alves Borges and Magnel Lima de Oliveira and [sjPlot](https://www.rdocumentation.org/packages/sjPlot/versions/2.8.12) package (Daniel Lüdecke and Carsten Schwemmer).
