
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AMEND

<!-- badges: start -->
<!-- badges: end -->

AMEND (Active Module identification with Experimental data and Network
Diffusion) is an algorithm designed to find a subset of connected nodes
in a molecular interaction network that have large experimental values.
It makes use of random walk with restart (RWR) to create node weights,
and a heuristic approach for solving the Maximum-weight Connected
Subgraph problem using these weights. This is performed iteratively
until an optimal subnetwork (i.e., module) is found. AMEND can now
accommodate multiplex and/or heterogeneous networks, making it a widely
applicable tool. These complex networks can include several node types,
edge types, and/or data types, which increases the types of biological
questions one can address with this method.

For a simple example, a user could have log fold changes from a
differential expression analysis of RNA-seq data comparing a treatment
group to a control group. Genes included in the experiment can be mapped
to a protein-protein interaction network, along with their log fold
changes. This network is the main input into AMEND, from which a single,
connected module is returned.

## Installation

You can install AMEND from [GitHub](https://github.com/samboyd0/AMEND)
with:

``` r
devtools::install_github("samboyd0/AMEND", build_vignettes = TRUE)
```

## Vignette

A vignette is available that illustrates how to use AMEND. It can be
accessed with the following code.

``` r
vignette("amend_vignette", package = "AMEND")
```

A short tutorial is available in the *vignettes/Tutorials* folder on the
AMEND GitHub page.

<!-- ## Example -->
<!-- AMEND contains three objects -->
<!-- ```{r example} -->
<!-- library(AMEND) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>. -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
