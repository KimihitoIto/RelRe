# Copyright 2021 Chayada Piantham, Kimihito Ito

include("lib-GT-BA1.jl")
include("lib-qt.jl")

using CSV, Dates, DataFrames, Distributions, NLopt

println("loading counts")
df_count = DataFrame(CSV.File("Tokyo_BA1_BA2.csv"))

variants = [:Omicron_BA1, :Omicron_BA2]
const baseline= :Omicron_BA1
subjects = [:Omicron_BA2]

dict_index = Dict{Symbol,Int64}()
map(x -> dict_index[subjects[x]] = x, 1:length(subjects))
dict_index[baseline] = length(subjects)+1

const t_start = Date("2022-02-01")
const t_end = maximum(df_count.date)+Day(6)

dict_first = Dict{Symbol,Date}()
map(v -> dict_first[v]=minimum(filter(v => n -> n>0, df_count).date), variants)
dict_first[:Omicron_BA2]=Date("2022-02-01")

dates = filter(row->row.date >= t_start && row.date+Day(6) <= t_end, df_count).date
mat_obs = Matrix(filter(x->x.date in dates, df_count)[:,vcat(subjects,baseline)])

function negLogL(par::Vector, grad::Vector)
    # println(par)
    vec_c = par[1:length(subjects)]
    vec_k = par[length(subjects)+1:length(subjects)*2]
    vec_qt = par[length(subjects)*2+1:length(subjects)*3]
    vec_t = [dict_first[:Omicron_BA2]]
    try
        q = model_q(vec_c, vec_k, vec_qt, vec_t, t_start, t_end)
        sumll = 0.0
        for i in 1:length(dates)
            j = Dates.value(dates[i]-t_start)+1
            rows = map(v -> j + v, 0:6)
            probs = vec(max.(0, mean(q[rows, 1:length(subjects)+1], dims=1)))
            obs = mat_obs[i,:]
            sumll += logpdf(Multinomial(sum(obs), probs), obs)
        end
        if !isfinite(sumll)
            return floatmax(Float64)
        end
        # println(sumll)
        return -sumll
    catch e
        println(e)
        throw(e)
    end
end

vec_c_start = fill(1.0,length(subjects))
vec_c_lb = fill(1.0,length(subjects))
vec_c_ub = fill(1.0,length(subjects))
vec_k_start = fill(1.0,length(subjects))
vec_k_lb = fill(1.0e-10,length(subjects))
vec_k_ub = fill(10.0,length(subjects))
vec_qt_start = fill(0.001,length(subjects))
vec_qt_lb = fill(1.0e-10,length(subjects))
vec_qt_ub = fill(1.0,length(subjects))

par_start = vcat(vec_c_start, vec_k_start, vec_qt_start)
par_lb = vcat(vec_c_lb, vec_k_lb, vec_qt_lb)
par_ub = vcat(vec_c_ub, vec_k_ub, vec_qt_ub)

opt = Opt(:LN_SBPLX, length(par_start))
opt.min_objective = (par, grad) -> negLogL(par, grad)
opt.lower_bounds = par_lb
opt.upper_bounds = par_ub
opt.ftol_abs = 1e-12
opt.maxeval = 500000

println("Maximizing the likehood function")
nmaxll, par_maxll, err = optimize(opt, par_start)
println("Maximization finished")

num = length(subjects)

df_ml = DataFrame()
map(x -> df_ml[!,"t_" * string(x)] = [dict_first[x]],subjects)
map(x -> df_ml[!,"c_" * string(x)] = [par_maxll[dict_index[x]]], subjects)
map(x -> df_ml[!,"k_" * string(x)] = [par_maxll[dict_index[x]+num]], subjects)
map(x -> df_ml[!,"qt_"* string(x)] = [par_maxll[dict_index[x]+num*2]], subjects)
df_ml[!,"maxll"] = [-nmaxll]
df_ml[!,"convergence"] = [err]
CSV.write("estimates_ml.csv", df_ml)

#confidence intervals

function constr(vec_par, vec_grad)
    -1.92 - nmaxll + negLogL(vec_par, vec_grad);
end

function f1(par, grad,i)
    for x = 1:length(grad)
        grad[x]= (x==i) ? 1.0 : 0.0
    end
    par[i]
end

function f2(par, grad,i)
    for x = 1:length(grad)
        grad[x]= (x==i) ? -1.0 : 0.0
    end
    -par[i]
end

mat_par_95 = Matrix{Float64}(undef, length(subjects)*6, length(subjects)*3)

println("Calculating 95% confidence intervals (CIs)")

Threads.@threads for i in 1:length(subjects)
    println("Thread " * string(Threads.threadid()) * " calculates CIs of " *
        string(subjects[i]))
    
    opt_c = Opt(:AUGLAG, length(par_maxll))
    opt_c.lower_bounds = par_lb
    opt_c.upper_bounds = par_ub
    opt_c.maxeval = 500000
    opt_c.ftol_abs = 1e-12

    inequality_constraint!(opt_c, (par,grad) -> constr(par,grad), 1e-6)
    
    opt_l = NLopt.Opt(:LN_SBPLX, length(par_maxll))
    opt_l.xtol_rel = 1e-6
    opt_c.local_optimizer = opt_l

    opt_c.min_objective = (par, grad) -> f1(par, grad, i)
    lb, par_95, err_95 = optimize(opt_c, par_maxll)
    mat_par_95[i, :] = par_95
    
    opt_c.min_objective = (par, grad) -> f2(par, grad, i)
    ub, par_95, err_95 = optimize(opt_c, par_maxll)
    mat_par_95[length(subjects) + i, :] = par_95

    opt_c.min_objective = (par, grad) -> f1(par, grad, num+i)
    lb, par_95, err_95 = optimize(opt_c, par_maxll)
    mat_par_95[2*length(subjects) + i, :] = par_95

    opt_c.min_objective = (par, grad) -> f2(par, grad, num+i)
    ub, par_95, err_95 = optimize(opt_c, par_maxll)
    mat_par_95[3*length(subjects) + i, :] = par_95

    opt_c.min_objective = (par, grad) -> f1(par, grad, 2*num+i)
    lb, par_95, err_95 = optimize(opt_c, par_maxll)
    mat_par_95[4*length(subjects) + i, :] = par_95

    opt_c.min_objective = (par, grad) -> f2(par, grad, 2*num+i)
    ub, par_95, err_95 = optimize(opt_c, par_maxll)
    mat_par_95[5*length(subjects) + i, :] = par_95
    
    println("Calculation for CIs of " * string(subjects[i]) * " finished")
end

df_95 = DataFrame()
map(x -> df_95[!,"t_" * string(x)] = repeat([dict_first[x]],length(subjects)*6),subjects)
map(x -> df_95[!,"c_" * string(x)] = mat_par_95[:,dict_index[x]], subjects)
map(x -> df_95[!,"k_" * string(x)] = mat_par_95[:,dict_index[x]+num], subjects)
map(x -> df_95[!,"qt_"* string(x)] = mat_par_95[:,dict_index[x]+num*2], subjects)
CSV.write("estimates_95CI.csv", df_95)

