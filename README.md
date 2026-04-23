
# sparseDGLM

<!-- badges: start -->
<!-- badges: end -->

The `sparseDGLM` package is designed to help users fit spatiotemporal dynamic generalized linear models to bicycle count data. It also supports high-dimensional spatial modeling using the SPDE approach introduced by Lindgren et al. (2011).

With this package, we may perform simulation-based MCMC inference for spatiotemporal prediction, including tasks like spatial interpolation and future forecasting. Additionally, packages also return an estimate of annual average bicycle counts, which is a key metric for understanding and improving bicycle infrastructure.

The package provides summaries of model parameter estimates for both sparse SPDE-based models and their dense counterparts, which use full Matérn covariance matrices. It also includes MCMC sampler for fitting Poisson generalized linear models (GLMs), offering flexibility depending on our needs.


## Installation

You can install the development version of sparseDGLM from [GitHub](https://github.com/) with:

``` r
remotes::install_github("yadavrishikesh/sparseDGLM") 
```

## Vignette File

For a detailed description of the package, including model specifications, implementation steps, and examples demonstrating how to generate results and reproduce outputs, please see the vignette available at the link below:

👉 [Vignette](https://yadavrishikesh.github.io/sparseDGLM/)
