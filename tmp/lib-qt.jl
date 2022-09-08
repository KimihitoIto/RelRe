using Dates

#----- Relative instantaneous reproduction number
# Ito, 2021 (https://doi.org/10.2807/1560-7917.ES.2021.26.27.2100570)

# l should be defined elsewhere

function model_q(vec_c::Vector{Float64}, vec_k::Vector{Float64},
                 vec_qt::Vector{Float64}, vec_t::Vector{Date}, t_s, t_e)
    num = length(vec_c)
    @assert length(vec_k) == num
    @assert length(vec_qt) == num
    @assert length(vec_t) == num

    duration = (t_e - t_s).value + 1

    g = Matrix{Float64}(undef, l, num + 1)
    for j in 1:num
        g[:,j] = pmf_g(vec_c[j])
    end
    g[:,length(vec_c)+1] = pmf_g(1.0)
    
    q = Matrix{Float64}(undef, duration, num + 1)
    for j in 1:num
        if vec_t[j] <= t_s
            q[1,j] = vec_qt[j]
        else
            q[1,j] = 0.0
        end
    end
    sum_q_subjects =sum(q[1,1:num])
    if(sum_q_subjects > 1)
        map(j -> q[1,j] /= sum_q_subjects,1:num)
        q[1, length(vec_k)+1] = 0
    else
        q[1, length(vec_k)+1] = 1 - sum_q_subjects
    end
    
    vec_sum_nmr=Vector{Float64}(undef, num)
    
    for i in 2:duration
        t = t_s + Day(i-1)
        fill!(vec_sum_nmr, 0.0)
        sum_dnm = 0.0
        for k in 1:l
            t_k = max(1, (t - Day(k) - t_s).value + 1)
            vec_sum_nmr += vec(g[k, 1:num]).* vec_k .* vec(q[t_k, 1:num])
            sum_dnm += g[k, num + 1] * q[t_k, num + 1] +
                sum(vec(g[k, 1:num]).* vec_k .* vec(q[t_k, 1:num]))
        end
        map(j -> q[i,j] = (vec_t[j] == t) ? vec_qt[j] :
            vec_sum_nmr[j] / sum_dnm, 1:num)
        
        sum_q_subjects =sum(q[i,1:num])
        
        if(sum_q_subjects > 1)
            map(j -> q[i,j] /= sum_q_subjects,1:num)
            q[i, num + 1] = 0.0
        else
            q[i, num + 1] = 1 - sum_q_subjects
        end
    end
    return q
end
