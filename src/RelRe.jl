using CSV, Dates, DataFrames, Distributions, NLopt, ArgParse
#julia --threads 5 RelRe.jl -c yearseason_freq.csv -b 6B.1 -s 2020-01-01 -e 2021-01-01
#TODO: Add abort/error handeling!

#Pull in information from command line arguments
s = ArgParseSettings()
@add_arg_table s begin
    "--count", "-c"
    help = "an option with an argument"
    arg_type = String
    required = true
    "--start", "-s"
    arg_type = String
    default = ""
    "--end", "-e"
    arg_type = String
    default = ""
    "--precision", "-p"
    arg_type = Float64
    default = 1e-4
    "--baseline", "-b"
    arg_type = Symbol
    required = true
    "--len", "-l"
    arg_type = Int64
    default = 16 # use 7 for flu
    "--delay", "-d"
    arg_type = Int64
    default = 0
    "--alpha", "-a"
    arg_type = Float64
    default = 2.03 # use 4.5 for flu
    "--theta", "-t"
    arg_type = Float64
    default = 2.32 # use 0.60 for flu
    "--unit", "-u"
    arg_type = Symbol
    default = :D
    required = true
    "--estimateRGT", "-g"
    action = :store_true
end

parsed_args = parse_args(ARGS, s)
@show parsed_args

#Init variables
const ftol_prec = parsed_args["precision"]
const matrixFile = parsed_args["count"]
const baseline = parsed_args["baseline"]
const startDate = parsed_args["start"]
const endDate = parsed_args["end"]
const l = parsed_args["len"]
const alpha = parsed_args["alpha"]
const theta = parsed_args["theta"]
const unit = parsed_args["unit"]
const estimateRGT = parsed_args["estimateRGT"]
const detection_delay = parsed_args["delay"]

#Generation time distribution
function g2(a, c_GT)
    if(a < 1 || a > l )
        return 0
    else
        return (cdf(Gamma(alpha, c_GT * theta), a) -
            cdf(Gamma(alpha, c_GT * theta), a-1))/
            cdf(Gamma(alpha, c_GT * theta),l)
    end
end

#Generation time distribution (old version)
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

function pmf_g(c_GT)
    return map(v -> g2(v,c_GT), 1:l)
end

#Renewal model of variant requencies
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

#Main

#Load in specified matrix file
println("Loading counts")
df_count = DataFrame(CSV.File(matrixFile))
@show(df_count)

#Check whether the first column is a vector of dates
if(typeof(df_count[:,1])!=Vector{Date})
    error("The first column is not a vector of dates")
end

#Check consistency with unit
if(unit==:M)
    if(length(unique(Dates.day.(df_count.date))) != 1)
        error("Dates should be the same day of the month")
    end
elseif(unit==:W)
    if(length(unique(Dates.dayofweek.(df_count.date))) != 1)
        error("Dates should be the same day of the week")
    end
end

#List the baseline
println("\nBaseline")
println(baseline)

#Specify variants as all, and remove baseline from subjects
variants = propertynames(df_count)[2:size(df_count,2)]
subjects = filter(x -> x!= baseline, variants)
println("\nSubject clades")
@show(subjects)

dict_index = Dict{Symbol,Int64}()
map(x -> dict_index[subjects[x]] = x, 1:length(subjects))
dict_index[baseline] = length(subjects)+1

#Remove data after end date
if(endDate!="")
    deleteat!(df_count, df_count.date .> Date(endDate))
end

#Adjust date range 
if(unit==:M)
    const t_end = maximum(df_count.date)+Day(Dates.daysinmonth(maximum(df_count.date))-1)
elseif(unit==:W)
    const t_end = maximum(df_count.date)+Day(6)
elseif(unit==:D)
    const t_end = maximum(df_count.date)
else
    println("error: the unit is not specified")
    exit()
end

#Record the date of variant's first observation
dict_first = Dict{Symbol,Date}()
map(v -> dict_first[v]=minimum(filter(v => n -> n>0, df_count).date)-Day(detection_delay), variants)

#Remove data before start date
if(startDate!="")
    const t_start = Date(startDate)
    deleteat!(df_count, df_count.date .< Date(startDate))
else
    const t_start = minimum(df_count.date)
end

println("\nTime range of analysis")
println("Start: " * Dates.format(t_start, dateformat"yyyy-mm-dd"))
println("End: " * Dates.format(t_end, dateformat"yyyy-mm-dd"))

dates = df_count.date
mat_obs = Matrix(filter(x->x.date in dates, df_count)[:,vcat(subjects,baseline)])

function negLogL(par::Vector, grad::Vector)
    # println(par)
    vec_c = par[1:length(subjects)]
    vec_k = par[length(subjects)+1:length(subjects)*2]
    vec_qt = par[length(subjects)*2+1:length(subjects)*3]
    vec_t = map(v -> dict_first[v], subjects)
    try
        q = model_q(vec_c, vec_k, vec_qt, vec_t, t_start, t_end)
        sumll = 0.0
        for i in 1:length(dates)
            j = Dates.value(dates[i]-t_start)+1
            if(unit==:M)
                rows = map(v -> j + v, 0:Dates.daysinmonth(dates[i])-1)
            elseif(unit==:W)
                rows = map(v -> j + v, 0:6)
            elseif(unit==:D)
                rows = [j]
            end
            probs = vec(max.(0, mean(q[rows, 1:length(subjects)+1], dims=1)))
            obs = mat_obs[i,:]
            sumll += logpdf(Multinomial(sum(obs), probs), obs)
        end
        if !isfinite(sumll)
            return floatmax(Float64)
        end
        return -sumll
    catch e
        println(e)
        throw(e)
    end
end

vec_c_start = fill(1.0,length(subjects))
if estimateRGT
    if unit == :D
        vec_c_lb = fill(1.0e-10,length(subjects))
        vec_c_ub = fill(10.0,length(subjects))
    else
        println("error: The unit should be a day for the -g option")
        exit(1)
    end
else
    vec_c_lb = fill(1.0,length(subjects))
    vec_c_ub = fill(1.0,length(subjects))
end
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

