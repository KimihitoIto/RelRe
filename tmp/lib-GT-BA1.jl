using Distributions

const l = 20

#----- serial interval distribution
# Hart, 2022 (https://doi.org/10.1016/S1473-3099(22)00001-9)
# Assuming GT_BA1 = 0.60 * GT_Delta

const gt_mean = 4.7
const gt_sd = 3.3

const alpha = 2.03
const theta = 2.32 * 0.60

function g2(a, c_GT)
    if(a < 1 || a > l )
        return 0
    else
        return (cdf(Gamma(alpha, c_GT * theta), a) -
            cdf(Gamma(alpha, c_GT * theta), a-1))/
            cdf(Gamma(alpha, c_GT * theta),l)
    end
end

function g1(a, c_GT)
    if(a == 1)
        return cdf(Gamma(alpha, c_GT * theta), 2)
    elseif(a == l)
        return 1 - cdf(Gamma(alpha, c_GT * theta), a)
    else
        return cdf(Gamma(alpha, c_GT * theta), a + 1) -
            cdf(Gamma(alpha, c_GT * theta), a)
    end
end

function pmf_g(c_SI)
    return map(v -> g2(v,c_SI),1:20)
end

