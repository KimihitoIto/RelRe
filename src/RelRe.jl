using CSV, Dates, DataFrames, Distributions, NLopt, ArgParse
#julia --threads 5 RelRe.jl -b baseline -s 2020-01-01 -e 2021-01-01 variant_freq.csv 
#TODO: Add abort/error handeling!
#TODO: An option for use of the first date.

#Pull in information from command line arguments
s = ArgParseSettings()
@add_arg_table s begin
    "--count_file", "-c"
    arg_type = String
    required = true
    "--start", "-s"
    arg_type = String
    default = ""
    "--end", "-e"
    arg_type = String
    default = ""
    "--baseline", "-b"
    arg_type = Symbol
    required = true
    "--subjects", "-j"
    arg_type = Symbol
    nargs = '*'
    default = []
    "--breaks", "-w"
    arg_type = Date
    nargs = '*'
    default = []
    "--precision", "-p"
    arg_type = Float64
    default = 1e-4
    "--len", "-l"
    arg_type = Int64
    default = 16 # use 7 for flu
    "--alpha", "-a"      # TODO: estimate from mean and var 
    arg_type = Float64
    default = 2.03 # use 4.5 for flu
    "--theta", "-t"      # TODO: estimate from mean and var
    arg_type = Float64
    default = 2.32 # use 0.60 for flu
    "--unit", "-u"
    arg_type = Symbol
    default = :D
    required = true
    "--estimate_GT", "-g"
    action = :store_true
    "--estimate_CI", "-i" # no need. set always true
    action = :store_true
end

parsed_args = parse_args(ARGS, s)
@show parsed_args

#Init variables
const ftol_prec = parsed_args["precision"]
const matrixFile = parsed_args["count_file"]
const baseline = parsed_args["baseline"]
const start_date = parsed_args["start"]
const end_date = parsed_args["end"]
const len_tr = parsed_args["len"]
const alpha = parsed_args["alpha"]
const theta = parsed_args["theta"]
const unit = parsed_args["unit"]
const breaks = parsed_args["breaks"]
const arg_subjects = parsed_args["subjects"]
const estimate_GT = parsed_args["estimate_GT"]
const estimate_CI = parsed_args["estimate_CI"]

#Generation time distribution
function g2(a, c_GT)
    if(a < 1 || a > len_tr )
        return 0
    else
        return (cdf(Gamma(alpha, c_GT * theta), a) -
            cdf(Gamma(alpha, c_GT * theta), a-1))/
            cdf(Gamma(alpha, c_GT * theta),len_tr)
    end
end

#Generation time distribution (old version)
function g1(a, c_GT)
    if(a == 1)
        return cdf(Gamma(alpha, c_GT * theta), 2)
    elseif(a == len_tr)
        return 1 - cdf(Gamma(alpha, c_GT * theta), a)
    else
        return cdf(Gamma(alpha, c_GT * theta), a + 1) -
            cdf(Gamma(alpha, c_GT * theta), a)
    end
end

function pmf_g(c_GT)
    return map(v -> g2(v,c_GT), 1:len_tr)
    #return map(v -> g1(v,c_GT), 1:len_tr)
end

#Renewal model of variant requencies
function model_q(vec_c::Vector{Float64}, vec_k::Vector{Float64},
                 vec_qt::Vector{Float64}, vec_t::Vector{Date},
                 vec_b::Vector{Date},t_s, t_e, l)
    @assert length(vec_c) == num_subjects
    @assert length(vec_k) == num_subjects * (length(vec_b) + 1)
    @assert length(vec_qt) == num_subjects
    @assert length(vec_t) == num_subjects
    
    duration = (t_e - t_s).value + 1

    g = Matrix{Float64}(undef, l, num_subjects + 1)
    for j in 1:num_subjects
        g[:,j] = pmf_g(vec_c[j])
    end
    g[: , num_subjects + 1] = pmf_g(1.0)
    
    q = Matrix{Float64}(undef, duration, num_subjects + 1)
    for j in 1:num_subjects
        if vec_t[j] <= t_s
            q[1,j] = vec_qt[j]
        else
            q[1,j] = 0.0
        end
    end
    sum_q_subjects =sum(q[1,1:num_subjects])
    if(sum_q_subjects > 1)
        map(j -> q[1,j] /= sum_q_subjects,1:num_subjects)
        q[1, num_subjects + 1] = 0
    else
        q[1, num_subjects + 1] = 1 - sum_q_subjects
    end
    
    vec_sum_nmr=Vector{Float64}(undef, num_subjects)
    
    for i in 2:duration
        t = t_s + Day(i-1)
        window_index = 1
        for b in vec_b
            if t < b
                break
            else
                window_index += 1
            end
        end
        subvec_k = map(j -> vec_k[(window_index - 1) * num_subjects + j],1:num_subjects)
        
        fill!(vec_sum_nmr, 0.0)
        sum_dnm = 0.0
        for k in 1:l
            t_k = max(1, (t - Day(k) - t_s).value + 1)
            vec_sum_nmr += vec(g[k, 1:num_subjects]).* subvec_k .* vec(q[t_k, 1:num_subjects])
            sum_dnm += g[k, num_subjects + 1] * q[t_k, num_subjects + 1] +
                sum(vec(g[k, 1:num_subjects]).* subvec_k .* vec(q[t_k, 1:num_subjects]))
        end
        map(j -> q[i,j] = (vec_t[j] == t) ? vec_qt[j] :
            vec_sum_nmr[j] / sum_dnm, 1:num_subjects)
        
        sum_q_subjects =sum(q[i,1:num_subjects])
        
        if(sum_q_subjects > 1)
            map(j -> q[i,j] /= sum_q_subjects,1:num_subjects)
            q[i, num_subjects + 1] = 0.0
        else
            q[i, num_subjects + 1] = 1 - sum_q_subjects
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
if(arg_subjects == [])
    variants = propertynames(df_count)[2:size(df_count,2)]
    subjects = filter(x -> x!= baseline, variants)
else
    subjects = arg_subjects # elements need to be checked
    variants = vcat(subjects,baseline)
end
println("\nSubject clades")
const num_subjects = length(subjects)
@show(subjects)

dict_index = Dict{Symbol,Int64}()
map(x -> dict_index[subjects[x]] = x, 1:length(subjects))
dict_index[baseline] = length(subjects)+1

#Remove data after end date
if(end_date!="")
    deleteat!(df_count, df_count.date .> Date(end_date))
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

#Remove data before start date
if(start_date!="")
    const t_start = Date(start_date)
    deleteat!(df_count, df_count.date .< Date(start_date))
else
    const t_start = minimum(df_count.date)
end

#Record the date of variant's first observation during the period 
dict_first = Dict{Symbol,Date}()
map(v -> dict_first[v]=minimum(filter(v => n -> n>0, df_count).date), variants)



println("\nTime range of analysis")
println("Start: " * Dates.format(t_start, ISODateFormat))
println("End: " * Dates.format(t_end, ISODateFormat))
println("Breakpoints: [" *
    join(map(d -> Dates.format(d, ISODateFormat),breaks),",") * "]")

#Check the breakpoints
for d in breaks
    if !(d > t_start + Day(1) && d < t_end)
        error("Breakpoints should be > t_start + 1 and < t_end")
    end
end

const num_windows = length(breaks) + 1 
if (num_windows > 1)
    println("The analysis period is divided into "
            * string(num_windows) * " windows.")
end

dates = df_count.date
mat_obs = Matrix(filter(x->x.date in dates, df_count)[:,vcat(subjects,baseline)])

function negLogL(par::Vector, grad::Vector)
    vec_c = par[1:num_subjects]
    vec_k = par[num_subjects+1:num_subjects * (num_windows + 1)]
    vec_qt = par[num_subjects*(num_windows+1)+1:num_subjects*(num_windows+2)]
    vec_t = map(v -> dict_first[v], subjects)
    vec_b = breaks
    try
        q = model_q(vec_c, vec_k, vec_qt, vec_t, vec_b, t_start, t_end, len_tr)
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
            probs = vec(max.(0, mean(q[rows, 1:num_subjects+1], dims=1)))
            obs = mat_obs[i,:]
            sumll += logpdf(Multinomial(sum(obs), probs), obs)
        end
        if !isfinite(sumll)
            return floatmax(Float64)
        end
        return -sumll
    catch e
        println(e)
        eixt(1)
    end
end

vec_c_start = fill(1.0,num_subjects)
if estimate_GT
    if unit == :D || unit == :W
        vec_c_lb = fill(1.0e-10,num_subjects)
        vec_c_ub = fill(10.0,num_subjects)
    else
        println("error: The unit should be a day or week for the -g option")
        exit(1)
    end
else
    vec_c_lb = fill(1.0,num_subjects)
    vec_c_ub = fill(1.0,num_subjects)
end
vec_k_start = fill(1.0,num_subjects * num_windows)
vec_k_lb = fill(1.0e-10,num_subjects * num_windows)
vec_k_ub = fill(10.0,num_subjects * num_windows)
vec_qt_start = fill(0.001,num_subjects)
vec_qt_lb = fill(1.0e-10,num_subjects)
vec_qt_ub = fill(1.0,num_subjects)

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

println(par_maxll)
println(-nmaxll) #todo: output to maxll.csv
println(err)

#Confidence Intervals
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

if estimate_CI
    mat_95CI = Matrix{Float64}(undef,
                                 num_subjects * (2 + 2 * num_windows + 2),
                                 num_subjects * (1 + 1 * num_windows + 1))
    
    println("Calculating 95% confidence intervals (CIs)")
    
    Threads.@threads for i in 1:(num_subjects * (2 + 2 * num_windows + 2))
        println("Thread " * string(Threads.threadid()) * " is working on " *
            string(i) * " th loop of the CI calculation")
        
        opt_c = Opt(:AUGLAG, length(par_maxll))
        opt_c.lower_bounds = par_lb
        opt_c.upper_bounds = par_ub
        opt_c.maxeval = 500000
        opt_c.ftol_abs = ftol_prec
        
        inequality_constraint!(opt_c, (par,grad) -> constr(par,grad), 1e-6)
        
        opt_l = NLopt.Opt(:LN_SBPLX, length(par_maxll))
        #opt_l.xtol_rel = 1e-6 #todo precision        
        opt_l.xtol_rel = 1e-4 #todo precision
        opt_c.local_optimizer = opt_l
        
        if(i % 2 == 1) # Lower bound
            opt_c.min_objective = (par, grad) -> f1(par, grad, i)
        else # Upper bound
            opt_c.min_objective = (par, grad) -> f2(par, grad, i)
        end
        lb, par_95, err_95 = optimize(opt_c, par_maxll)
        mat_95CI[i, :] = par_95
        println("Calculation of " * string(i) * " th loop finished")
    end
end

df_estimates = DataFrame()
df_estimates[!,"variant"] = Vector{Symbol}()
df_estimates[!,"date"] = Vector{Date}()
df_estimates[!,"c"] = Vector{Float64}()
df_estimates[!,"c_lb"] = Vector{Float64}()
df_estimates[!,"c_ub"] = Vector{Float64}()
if num_windows > 1
    for w in 1:num_windows
        df_estimates[!,"k_w"*string(w)] = Vector{Float64}()
        df_estimates[!,"k_w"*string(w)*"_lb"] = Vector{Float64}()
        df_estimates[!,"k_w"*string(w)*"_ub"] = Vector{Float64}()
    end
else
    df_estimates[!,"k"] = Vector{Float64}()
    df_estimates[!,"k_lb"] = Vector{Float64}()
    df_estimates[!,"k_ub"] = Vector{Float64}()
end
df_estimates[!,"qt"] = Vector{Float64}()
df_estimates[!,"qt_lb"] = Vector{Float64}()
df_estimates[!,"qt_ub"] = Vector{Float64}()

for j in 1:num_subjects
    row=vcat(subjects[j],
             dict_first[subjects[j]],
             par_maxll[j], #c
             minimum(mat_95CI[:,j]), #lb
             maximum(mat_95CI[:,j])) #ub
    if num_windows > 1
        for w in 1:num_windows
            row=vcat(row,par_maxll[num_subjects+(w-1)*num_subjects + j],#k
                     minimum(mat_95CI[:,num_subjects+(w-1)*num_subjects+j]),#lb
                     maximum(mat_95CI[:,num_subjects+(w-1)*num_subjects+j]))#ub
        end
        row=vcat(row,par_maxll[(1+num_windows)*num_subjects+j], #qt
                 minimum(mat_95CI[:,(1+num_windows)*num_subjects+j]),#lb
                 maximum(mat_95CI[:,(1+num_windows)*num_subjects+j]))#ub
    else
        row=vcat(row,par_maxll[num_subjects+j], #k
                 minimum(mat_95CI[:,num_subjects + j]), #lb
                 maximum(mat_95CI[:,num_subjects + j])) #ub
        row=vcat(row,par_maxll[2 * num_subjects+j], #qt
                 minimum(mat_95CI[:,2 * num_subjects + j]), #lb
                 maximum(mat_95CI[:,2 * num_subjects + j])) #ub 
    end
    push!(df_estimates,row)
end

CSV.write("estimates.csv", df_estimates)
