# GenInvGaussian

[![Build Status](https://travis-ci.org/LMescheder/GenInvGaussian.jl.svg?branch=master)](https://travis-ci.org/LMescheder/GenInvGaussian.jl)

This package contains an implementation of the [generalized inverse gaussian distribution](https://en.wikipedia.org/wiki/Generalized_inverse_Gaussian_distribution). Sampling is done by a simple rejection sampling procedure from a gamma distributions. This works well, but naturally this is not the fastest method.

## Usage
```julia
using Distributions
using GenInvGaussian

d = GeneralizedInverseGaussian(1., 10., 10.)

println("Mean = ", mean(d))
println("Variance = ", var(d))
println("Mode = ", mode(d))

using PyPlot
s = rand(d, 1000000)
plt[:hist](s, bins=100)
println("Empirical mean = ", mean(s))
println("Empiricial variance = ", var(s))
println("Empiricial mode = ", mode(s))
```
Output is
```
Mean = 1.153417250747945
Variance = 0.13099554597623
Mode = 1.0
Empirical mean = 1.1538116272315135
Empiricial variance = 0.1308216411777261
Empiricial mode = 0.9096948416390549
```
