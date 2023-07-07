# Compositional Regression Analysis 
R Script to estimate and impute cell type proportions from DNA methylation, and run compositional regression on example data set from planet R package.
In this example, we estimate the expected change in placenta cell type composition associated with 1 week increase in gestational age, along with the corresponding bootstrapped confidence intervals and corrected p-values.

# NOTE
One of the limitations of the 'compositions' package is that the p-value calculation is only correct for the last covariate in the model formula.  
In this example, the p-value for "gestation_wk" is correct, but not for infant sex. In order to calculate the correct p-value for "sex", 
you need to refit the same model with "sex" at the end of the formula: 
`ilr(outcome) ~ gestation_wk + sex`

# Reference
For further details, see Section 5.3.5 of "Analyzing compositional data with R" by K. Gerald van den Boogaart and Raimon Tolosana-Delgado (2013).
