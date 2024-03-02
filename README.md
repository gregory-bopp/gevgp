# Metropolis-Hastings Samplers for Gaussian Processes and GEV-Gaussian Processes

This repository contains two R packages that contain Metropolis Hastings samplers for two different spatial processes. 

- **gpbayes**: is for a stationary [Gaussian processes](https://en.wikipedia.org/wiki/Gaussian_process) with a Matern covariance.

- **gevgpbayes**: an MCMC sampler for a hierarchical model, first proposed by [[1]](https://link.springer.com/article/10.1007/s13253-009-0010-1), where the marginal distributions at each point in space are Generalized Extreme Value (GEV), with spatially varying parameters, and the dependence copula is taken to be a Gaussian Process, in this case, again with a stationary Matern covariance.

References:
[Sang, Huiyan, and Alan E. Gelfand. "Continuous spatial process models for spatial extreme values." Journal of agricultural, biological, and environmental statistics 15.1 (2010): 49-65](https://link.springer.com/article/10.1007/s13253-009-0010-1)