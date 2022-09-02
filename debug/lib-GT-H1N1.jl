# This program was written by Kimihito Ito.
using Distributions
const l = 7 
#----- generation time distribution
# Lessler 2009 (https://www.nejm.org/doi/full/10.1056/nejmoa0906089)
const alpha = 4.5
const theta = 0.60
g1(a) = cdf(Gamma(alpha, theta), a) - cdf(Gamma(alpha, theta), a-1)
function g2(a, c_GT)
    if(a < 1 || a > l )
        return 0
    else
        return (cdf(Gamma(alpha, c_GT * theta), a) -
            cdf(Gamma(alpha, c_GT * theta), a-1))/
            cdf(Gamma(alpha, c_GT * theta),l)
    end
end
function pmf_g(c_GT)
    return map(v -> g2(v,c_GT), 1:l)
end
