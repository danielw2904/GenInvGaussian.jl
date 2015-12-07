module GenInvGaussian

export GeneralizedInverseGaussian

import Base.mean
import Base.rand
using Distributions
import Distributions: mean, var, mode, mgf, cf, rand, params

immutable GeneralizedInverseGaussian <: ContinuousUnivariateDistribution
    p::Float64
    a::Float64
    b::Float64

    function GeneralizedInverseGaussian(p::Real, a::Real, b::Real)
        #@checkargs(GeneralizedInverseGaussian, a > zero(a), b > zero(b))
        new(p, a, b)
    end
end

#@distr_support GeneralizedInverseGaussian 0.0 Inf

#### Parameters

params(d::GeneralizedInverseGaussian) = (d.p, d.a, d.b)

#### Statistics

function mean(d::GeneralizedInverseGaussian)
    (p, a, b) = params(d)
    ω = sqrt(a*b)
    sqrt(b/a) * besselk(p + 1, ω) / besselk(p, ω)
end

function var(d::GeneralizedInverseGaussian)
    (p, a, b) = params(d)
    ω = sqrt(a*b)
    num = besselk(p, ω)*besselk(p+2, ω) - besselk(p+1, ω)*besselk(p+1, ω)
    denom = besselk(p, ω) * besselk(p, ω)
    b/a * num / denom
end

function mode(d::GeneralizedInverseGaussian)
    (p, a, b) = params(d)
    m = (p - 1) + sqrt((p - 1)*(p - 1) + a*b)
    if m < 0
        error("Distribution has no mode!")
    end
    m / a
end

function mgf(d::GeneralizedInverseGaussian, t::Real)
    (p, a, b) = params(d)
    at = a - 2t
    ω = sqrt(a*b)
    ωt = sqrt(at*b)
    (a/at)^(p/2) * besselk(p, ωt)/besselk(p, ω)
end

function cf(d::GeneralizedInverseGaussian, t::Real)
    (p, a, b) = params(d)
    at = a - 2im * t
    ω = sqrt(d.a*d.b)
    ωt = sqrt(at*b)
    (a/at)^(p/2) * besselk(p, ωt)/besselk(p, ω)
end


#### Evaluation & Sampling

# Simple rejection sampling using a gamma distribution as a proposal distribution.
# Probably not the best solution out there. However, for p >= 1, the expected number of
# rejection is lower than sqrt(2).
function rand(d::GeneralizedInverseGaussian)
    (p, a, b) = params(d)
    ω = sqrt(a*b)

    r = p/2 + sqrt(p*p + a*b)/2
    dproposal = Gamma(r, 2/a)
    daccept = Uniform()

    xm = b / ( r - p) / 2

    acceptrate(x::Real) = (x./xm)^(p - r) * exp(-b/2*(1./x - 1./xm))

    while true
        sproposal = rand(dproposal)
        saccept = rand(daccept)

        if saccept <= acceptrate(sproposal)
            return sproposal
        end
    end
end

end # module
