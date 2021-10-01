# multiMaxEnt

An experimental R package to fit a multispecies extension of maxent (more specific maxnet implementation) using its equivalence with Poisson and Gibbs point process.

## The idea behind
[Fithian and Hastie 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4258396/) revealed the equivalence between [Phillips' MaxEnt](https://onlinelibrary.wiley.com/doi/full/10.1111/j.0906-7590.2008.5203.x) model, a Logistic regression with infinitely many "background points" a Poisson point process whose natural multi-class extension is Gibbs point process. In the same year, [Baddeley et al.](https://academic.oup.com/biomet/article/101/2/377/194771?login=true) proposed to use a similar equivalence between Logistic regression and a Gibbs point process which was implemented in the R package `spatstat.core`. 

This experimental package aim to combine the two paper and test the idea of using Gibbs point process to fit a multispecies version of MaxEnt.
