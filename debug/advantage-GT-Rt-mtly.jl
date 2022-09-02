# This program was written by Kimihito Ito and Chayada Piantham.
# See https://doi.org/10.2807/1560-7917.ES.2021.26.27.2100570
# and https://doi.org/10.3934/mbe.2022418
#julia --threads 5 advantage-GT-Rt-mtly.jl -m yearseason_freq.csv -b 6B.1 -s 2020-01-01 -e 2021-01-01

#Includes and packages
include("lib-GT-H1N1.jl")
include("lib-qt.jl")

using CSV, Dates, DataFrames, Distributions, NLopt, Getopt

#Init variables
global ftol_prec = 1e-4
global matrixFile = ""
global baseline = ""
global startDate = ""
global endDate = ""

#Pull in information
println(ARGS)
for (opt, arg) in getopt(ARGS, "m:s:e:p:b:", ["matrix=", "state=", "end=", "precision=", "baseline="])
	#@show (opt, arg)
	if opt == "-m" 
		global matrixFile = arg
	elseif opt == "-b" 
		global baseline = Symbol(arg)
	elseif opt == "-p" 
		global ftol_prec = arg
	elseif opt == "-s" 
		global startDate = arg #Try YYYY-MM-DD
	elseif opt == "-e" 
		global endDate = arg	
	end
end

#TODO: Add abort/error handeling!

#Load in specified matrix file
println("loading counts")
df_count = DataFrame(CSV.File(matrixFile))
@show(df_count)

#List the baseline
println("\nBaseline")
println(baseline)

#Specify variants as all, and remove baseline from subjects
variants = propertynames(df_count)[2:size(df_count,2)]
subjects = variants
deleteat!(subjects, subjects .== baseline)
println("\nSubject clades")
@show(subjects)

dict_index = Dict{Symbol,Int64}()
map(x -> dict_index[subjects[x]] = x, 1:length(subjects))
dict_index[baseline] = length(subjects)+1

#Remove data after end date
deleteat!(df_count, df_count.date .> Date(endDate))

#Adjust date range 
const t_start = Date(startDate)
const t_end = maximum(df_count.date)+Day(Dates.daysinmonth(maximum(df_count.date))-1)
#const t_end = Date(endDate) -1
println("\nTime range of analysis")
println("Start: " * Dates.format(t_start, dateformat"yyyy-mm-dd"))
println("End: " * Dates.format(t_end, dateformat"yyyy-mm-dd"))

dict_first = Dict{Symbol,Date}()
map(v -> dict_first[v]=t_start, variants)

months = filter(row->row.date >= t_start, df_count).date
mat_obs = Matrix(filter(x->x.date in months, df_count)[:,vcat(subjects,baseline)])

function negLogL(par::Vector, grad::Vector)
    # println(par)
    vec_c = par[1:length(subjects)]
    vec_k = par[length(subjects)+1:length(subjects)*2]
    vec_qt = par[length(subjects)*2+1:length(subjects)*3]
    vec_t = map(v -> dict_first[v], subjects)
    try
        q = model_q(vec_c, vec_k, vec_qt, vec_t, t_start, t_end)
        sumll = 0.0
        for i in 1:length(months)
            j = Dates.value(months[i]-t_start)+1
            rows = map(v -> j + v, 0:Dates.daysinmonth(months[i])-1)	#FIX THIS LINE
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
opt.ftol_abs = ftol_prec
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
CSV.write(endDate * "estimates_ml.csv", df_ml) #Allow multiple runs at once

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
    opt_c.ftol_abs = ftol_prec

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
CSV.write(endDate * "estimates_95CI.csv", df_95) #Allow multiple runs at once

