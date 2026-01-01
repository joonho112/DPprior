# DPprior: Principled Prior Elicitation for Dirichlet Process Mixture Models

Tools for eliciting Gamma hyperpriors on the concentration parameter
alpha in Dirichlet Process (DP) mixture models. Rather than requiring
specification of the abstract parameter alpha directly, 'DPprior' allows
researchers to express prior beliefs through intuitive quantities:
expected cluster counts and weight concentration. The package implements
the DORO 2.0 methodology with closed-form approximations (A1) and exact
Newton-based moment matching (A2), enabling principled prior
specification in low-information settings such as multisite trials and
meta-analyses. Features include dual-anchor elicitation for joint
control of cluster counts and weight behavior, comprehensive
diagnostics, and visualization tools. Based on methods described in Lee
et al. (2025)
[doi:10.3102/10769986241254286](https://doi.org/10.3102/10769986241254286)
.

## See also

Useful links:

- <https://joonho112.github.io/DPprior/>

- <https://github.com/joonho112/DPprior>

- Report bugs at <https://github.com/joonho112/DPprior/issues>

## Author

**Maintainer**: JoonHo Lee <jlee296@ua.edu>
([ORCID](https://orcid.org/0009-0006-4019-8703))
