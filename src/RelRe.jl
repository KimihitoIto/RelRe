using CSV, Dates, DataFrames, Distributions, NLopt, ArgParse
#julia --threads 5 RelRe.jl -b baseline -s 2020-01-01 -e 2021-01-01 variant_freq.csv 
#TODO: Add abort/error handeling!
#TODO: An option for use of the first date.

#Pull in information from command line arguments
s = ArgParseSettings()
@add_arg_table s begin
    "--in", "-i"
    arg_type = String
    required = true
    help = "input file containing temporal count data of variants"
    "--out", "-o"
    arg_type = String
    default = ""
    help = "prefix of output files"
    "--start", "-s"
    arg_type = String
    default = ""
    help = "start date of the analysis"
    "--end", "-e"
    arg_type = String
    default = ""
    help = "end date of the analysis"
    "--future", "-f"
    arg_type = Int64
    default = 0
    help = "duration in days for predicting variant frequencies"
    "--baseline", "-b"
    arg_type = Symbol
    required = true
    help = "variant used as the baseline of relative reproduction numbers"
    "--subjects", "-j"
    arg_type = Symbol
    nargs = '*'
    default = []
    help = "list of variants to calculate relative reproduction numbers"
    "--precision", "-p"
    arg_type = Float64
    default = 1e-4
    help = "stopping criterion used as ftol_abs in NLopt"
    "--len", "-l"
    arg_type = Int64
    default = 16 # use 7 for flu
    help = "trancation point of gamma distribution for generation time"
    "--alpha", "-a"      # TODO: estimate from mean and var 
    arg_type = Float64
    default = 2.03 # use 4.5 for flu
    help = "shape parameter of gamma distribution for generation time"
    "--theta", "-t"      # TODO: estimate from mean and var
    arg_type = Float64
    default = 2.32 # use 0.60 for flu
    help = "scale parameter of gamma distribution for generation time"
    "--delta", "-d"
    arg_type = Float64
    default = 0.5
    help = "unit time of calculation (in days)"
    "--unit", "-u"
    arg_type = Symbol
    default = :D
    help = "unit time of observations: D (Daily),W (Weekly),or M (Monthly)"
    required = true
    "--Dirichlet", "-D"
    action = :store_true
    help = "use Dirichlet multinomial as the observation model"
    "--frequency", "-q"   
    action = :store_true
    help = "calculate the time course of variant frequencies"
    "--estimate_GT", "-g"
    action = :store_true
    help = "estimate relative generation times of variants"
    "--estimate_CI", "-c"
    action = :store_true
    help = "estimate 95% confidence intervals"
    "--undetected", "-n"
    action = :store_true
    help = "assume all variants exist undetected from the start date"
end

parsed_args = parse_args(ARGS, s)
@show parsed_args

#Init variables
const infile = parsed_args["in"]
const outfile_prefix = parsed_args["out"]
const ftol_prec = parsed_args["precision"]
const baseline = parsed_args["baseline"]
const start_date = parsed_args["start"]
const end_date = parsed_args["end"]
const days_to_predict = parsed_args["future"]
const len_tr = parsed_args["len"]
const alpha = parsed_args["alpha"]
const theta = parsed_args["theta"]
const unit = parsed_args["unit"]
const arg_subjects = parsed_args["subjects"]
const estimate_GT = parsed_args["estimate_GT"]
const estimate_CI = parsed_args["estimate_CI"]
const assume_undetected = parsed_args["undetected"]
const calculate_q = parsed_args["frequency"]
const delta = parsed_args["delta"]
const dirichlet = parsed_args["Dirichlet"]

#Generation time distribution
function gt(a, c_GT)
    if(a < delta || a > len_tr )
        return 0
    else
        return (cdf(Gamma(alpha, c_GT * theta), a) -
            cdf(Gamma(alpha, c_GT * theta), a - delta))/
            cdf(Gamma(alpha, c_GT * theta),len_tr)
    end
end

function pmf_gt(c_GT)
    return map(v -> gt(v,c_GT), delta:delta:len_tr)
end

#Renewal model of variant requencies
function model_q(vec_c::Vector{Float64}, vec_k::Vector{Float64},
                 vec_qt::Vector{Float64}, vec_t::Vector{Date},
                 t_s, t_e, l)
    @assert length(vec_c) == num_subjects
    @assert length(vec_k) == num_subjects
    @assert length(vec_qt) == num_subjects
    @assert length(vec_t) == num_subjects
    
    duration = (t_e - t_s).value + 1
    #todo: calculate l automatically
    
    g = Matrix{Float64}(undef, length(delta:delta:l), num_subjects + 1)
    for j in 1:num_subjects
        g[:,j] = pmf_gt(vec_c[j])
    end
    g[: , num_subjects + 1] = pmf_gt(1.0)
    
    q = Matrix{Float64}(undef, length(delta:delta:duration), num_subjects + 1)
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
    
    for i in 2:length(delta:delta:duration)
        t = t_s + Day(Int64(floor((delta:delta:duration)[i])))
        
        fill!(vec_sum_nmr, 0.0)
        sum_dnm = 0.0
        for k in 1:length(delta:delta:l)#here
            t_k = max(1, i - k)
            vec_sum_nmr += vec(g[k, 1:num_subjects]).* vec_k .* vec(q[t_k, 1:num_subjects])
            sum_dnm += g[k, num_subjects + 1] * q[t_k, num_subjects + 1] +
                sum(vec(g[k, 1:num_subjects]).* vec_k .* vec(q[t_k, 1:num_subjects]))
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
    q_day = Matrix{Float64}(undef, duration, num_subjects + 1)
    for i in 1:duration
        row = minimum(findall(y -> y>=(i-1), 0:delta:duration))
        q_day[i, :] = q[row,:]
    end
    return q_day
end

#Main

#Check for values provided in program options
if(delta > 1.0 || !isinteger(1.0/delta))
    error("The delta should be a number obtained by deviding one by an integer")
end

#Load in specified matrix file
println("Loading counts")
df_count = DataFrame(CSV.File(infile))

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

#Remove unnecessary columns
@show(df_count)

#Record the date of variant's first observation during the period 
dict_first = Dict{Symbol,Date}()
if assume_undetected
    map(v -> dict_first[v]=t_start, variants)
else
    map(v -> dict_first[v]=minimum(filter(v => n -> n>0, df_count).date), variants)
end

println("\nTime range of analysis")
println("Start: " * Dates.format(t_start, ISODateFormat))
println("End: " * Dates.format(t_end, ISODateFormat))

dates = df_count.date
mat_obs = Matrix(filter(x->x.date in dates, df_count)[:,vcat(subjects,baseline)])

function negLogL(par::Vector, grad::Vector)
    if(dirichlet)
        @assert length(par) == 3 * num_subjects + 1
    else
        @assert length(par) == 3 * num_subjects
    end
        
    vec_c = par[1:num_subjects]
    vec_k = par[num_subjects+1:num_subjects * 2]
    vec_qt = par[num_subjects*2+1:num_subjects*3]
    if(dirichlet)
        M = par[num_subjects * 3 + 1]
    end
    
    vec_t = map(v -> dict_first[v], subjects)
    try
        q = model_q(vec_c, vec_k, vec_qt, vec_t, t_start, t_end, len_tr)
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
            
            if(dirichlet)
                alphas = probs * M
                sumll += logpdf(DirichletMultinomial(sum(obs), alphas), obs)
            else
                sumll += logpdf(Multinomial(sum(obs), probs), obs)
            end
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
vec_k_start = fill(1.0,num_subjects)
vec_k_lb = fill(1.0e-10,num_subjects)
vec_k_ub = fill(10.0,num_subjects)
vec_qt_start = fill(0.001,num_subjects)
vec_qt_lb = fill(1.0e-10,num_subjects)
vec_qt_ub = fill(1.0,num_subjects)

if(dirichlet)
    M_start = 10
    M_lb = 1.0e-10
    M_ub = 1.0e+5
    par_start = vcat(vec_c_start, vec_k_start, vec_qt_start, M_start)
    par_lb = vcat(vec_c_lb, vec_k_lb, vec_qt_lb, M_lb)
    par_ub = vcat(vec_c_ub, vec_k_ub, vec_qt_ub, M_ub)
else
    par_start = vcat(vec_c_start, vec_k_start, vec_qt_start)
    par_lb = vcat(vec_c_lb, vec_k_lb, vec_qt_lb)
    par_ub = vcat(vec_c_ub, vec_k_ub, vec_qt_ub)
end

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
    -quantile(Chisq(1),0.95)/2 - nmaxll + negLogL(vec_par, vec_grad);
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

if(dirichlet)
    mat_95CI = Matrix{Float64}(undef,
                               num_subjects * 6 + 2,
                               num_subjects * 3 + 1)
    err_95CI = Vector{Symbol}(undef,num_subjects * 6 + 2)
else
    mat_95CI = Matrix{Float64}(undef,
                               num_subjects * 6,
                               num_subjects * 3)
    err_95CI = Vector{Symbol}(undef,num_subjects * 6)
end

if estimate_CI
    println("\nCalculating 95% confidence intervals (CIs)")
    if(dirichlet)
        num_loop = (num_subjects * 6 + 2)
    else
        num_loop = (num_subjects * 6)
    end
    
    Threads.@threads for i in 1:num_loop
        println("Thread " * string(Threads.threadid()) * " is working on the " *
            string(i) * " th loop of the CI calculation")
        
        opt_c = Opt(:AUGLAG, length(par_maxll))
        opt_c.lower_bounds = par_lb
        opt_c.upper_bounds = par_ub
        opt_c.maxeval = 500000
        opt_c.ftol_abs = ftol_prec
        
        inequality_constraint!(opt_c, (par,grad) -> constr(par,grad), 1e-6)
        
        opt_l = NLopt.Opt(:LN_SBPLX, length(par_maxll))
        opt_l.ftol_abs = ftol_prec
        opt_c.local_optimizer = opt_l
        
        if(i % 2 == 1) # Lower bound
            opt_c.min_objective = (par, grad) -> f1(par, grad, Int64((i+1)/2))
        else # Upper bound
            opt_c.min_objective = (par, grad) -> f2(par, grad, Int64(i/2))
        end
        lb, par_95, err_95 = optimize(opt_c, par_maxll)
        mat_95CI[i, :] = par_95
        err_95CI[i] = err_95
        println("Calculation of the " * string(i) * " th loop finished")
    end
    println(err_95CI)
end

df_estimates = DataFrame()
df_estimates[!,"variant"] = Vector{Symbol}()
df_estimates[!,"date"] = Vector{Date}()
df_estimates[!,"c"] = Vector{Float64}()
if estimate_CI
    df_estimates[!,"c_lb"] = Vector{Float64}()
    df_estimates[!,"c_ub"] = Vector{Float64}()
end

df_estimates[!,"k"] = Vector{Float64}()
if estimate_CI
    df_estimates[!,"k_lb"] = Vector{Float64}()
    df_estimates[!,"k_ub"] = Vector{Float64}()
end

df_estimates[!,"qt"] = Vector{Float64}()
if estimate_CI
    df_estimates[!,"qt_lb"] = Vector{Float64}()
    df_estimates[!,"qt_ub"] = Vector{Float64}()
end

for j in 1:num_subjects
    row = vcat(subjects[j], dict_first[subjects[j]], par_maxll[j]) #c
    if estimate_CI
        row = vcat(row, minimum(mat_95CI[:,j]), maximum(mat_95CI[:,j])) #lb, ub
    end

    row=vcat(row, par_maxll[num_subjects+j])#k
    if estimate_CI
        row=vcat(row, minimum(mat_95CI[:,num_subjects + j]), #lb
                 maximum(mat_95CI[:,num_subjects + j])) #ub
    end
    row=vcat(row, par_maxll[2 * num_subjects+j])#qt
    if estimate_CI
        row=vcat(row, minimum(mat_95CI[:,2 * num_subjects + j]), #lb
                 maximum(mat_95CI[:,2 * num_subjects + j])) #ub
    end
    push!(df_estimates,row)
end

if outfile_prefix==""
    outfile_estimate = "estimates.csv"
else
    outfile_estimate = outfile_prefix * "_estimates.csv"
end

CSV.write(outfile_estimate, df_estimates)

if(dirichlet)
    df_M = DataFrame()
    df_M[!,"M"] = [par_maxll[3 * num_subjects + 1]]
    if estimate_CI
        df_M[!,"M_lb"] = [minimum(mat_95CI[:,3 * num_subjects + 1])]
        df_M[!,"M_ub"] = [maximum(mat_95CI[:,3 * num_subjects + 1])]
    end
    if outfile_prefix==""
        outfile_M = "Dirichlet.csv"
    else
        outfile_M = outfile_prefix * "_Dirichlet.csv"
    end
    CSV.write(outfile_M, df_M)
end

#Calculate Trajectory

if estimate_GT
    freedom = num_subjects * 3
else 
    freedom = num_subjects * 2
end

if dirichlet
    freedom = freedom + 1
end

function constr_trajectory(vec_par, vec_grad)
    -quantile(Chisq(freedom),0.95)/2 - nmaxll + negLogL(vec_par, vec_grad);
end

if(dirichlet)
    mat_95CR = Matrix{Float64}(undef,
                               num_subjects * 6 + 2,
                               num_subjects * 3 + 1)
    err_95CR = Vector{Symbol}(undef,num_subjects * 6 + 2)
else
    mat_95CR = Matrix{Float64}(undef,
                               num_subjects * 6,
                               num_subjects * 3)
    err_95CR = Vector{Symbol}(undef,num_subjects * 6)
end

if calculate_q
    println("\nCalculating trajectory")
    vec_c_ml = par_maxll[1 : num_subjects] #c
    vec_k_ml = par_maxll[num_subjects+1 : num_subjects*2] #k
    vec_qt_ml = par_maxll[num_subjects*2+1 : num_subjects*3] #qt
    vec_t_ml = map(v -> dict_first[v], subjects)

    t_future = t_end + Day(days_to_predict)
    q_ml = model_q(vec_c_ml, vec_k_ml, vec_qt_ml, vec_t_ml,
                   t_start, t_future, len_tr)
    q_lb = q_ml
    q_ub = q_ml
    vec_average_c_ml = q_ml * vcat(vec_c_ml, 1.0)
    vec_average_c_lb = vec_average_c_ml
    vec_average_c_ub = vec_average_c_ml

    vec_average_k_ml = q_ml * vcat(vec_k_ml, 1.0)
    vec_average_k_lb = vec_average_k_ml
    vec_average_k_ub = vec_average_k_ml
    
    if estimate_CI
        println("Calculating 95% confidence region (CR)")
        if(dirichlet)
            num_loop = (num_subjects * 6 + 2)
        else
            num_loop = (num_subjects * 6)
        end
        Threads.@threads for i in 1:num_loop
            println("Thread " * string(Threads.threadid()) *
                " is working on the " *
                string(i) * " th loop of the CR calculation")
            opt_c = Opt(:AUGLAG, length(par_maxll))
            opt_c.lower_bounds = par_lb
            opt_c.upper_bounds = par_ub
            opt_c.maxeval = 500000
            opt_c.ftol_abs = ftol_prec
            
            inequality_constraint!(opt_c, (par,grad) -> constr_trajectory(par,grad),0)
            opt_l = NLopt.Opt(:LN_SBPLX, length(par_maxll))
            opt_l.ftol_abs = ftol_prec
            opt_c.local_optimizer = opt_l
            
            if(i % 2 == 1) # Lower bound
                opt_c.min_objective = (par, grad) -> f1(par, grad, Int64((i+1)/2))
            else # Upper bound
                opt_c.min_objective = (par, grad) -> f2(par, grad, Int64(i/2))
            end
            lb, par_95, err_95 = optimize(opt_c, par_maxll)
            mat_95CR[i, :] = par_95
            err_95CR[i] = err_95
            println("Calculation of the " * string(i) * " th loop finished")
        end
        println(err_95CR)
        for i in 1:num_loop
            par_cr = mat_95CR[i, :]
            vec_c_cr = par_cr[1 : num_subjects] #c
            vec_k_cr = par_cr[num_subjects+1 : num_subjects*2] #k
            vec_qt_cr= par_cr[num_subjects*2+1 : num_subjects*3] #qt
            vec_t_cr = map(v -> dict_first[v], subjects)

            q_cr = model_q(vec_c_cr, vec_k_cr, vec_qt_cr, vec_t_cr, 
                           t_start, t_future, len_tr)
            local vec_average_c_cr = q_cr * vcat(vec_c_cr, 1.0) #local?
            vec_average_k_cr = q_cr * vcat(vec_k_cr, 1.0)
            
            global q_lb = min.(q_lb, q_cr)
            global q_ub = max.(q_ub, q_cr)
            global vec_average_c_lb = min.(vec_average_c_lb, vec_average_c_cr)
            global vec_average_c_ub = max.(vec_average_c_ub, vec_average_c_cr)
            global vec_average_k_lb = min.(vec_average_k_lb, vec_average_k_cr)
            global vec_average_k_ub = max.(vec_average_k_ub, vec_average_k_cr)
        end
    end
    println("\nWriting frequencies")
    df_freq = DataFrame()
    df_freq[!,"date"] = collect(t_start:Day(1):t_future)
    map(x -> df_freq[!,string(x)] = q_ml[:,dict_index[x]], subjects)
    df_freq[!,string(baseline)] = q_ml[:,length(subjects)+1]
    df_freq[!,"average_c"] = vec_average_c_ml
    df_freq[!,"average_k"] = vec_average_k_ml
    if estimate_CI
        map(x -> df_freq[!,string(x) * "_lb"] = q_lb[:,dict_index[x]], subjects)
        map(x -> df_freq[!,string(x) * "_ub"] = q_ub[:,dict_index[x]], subjects)
        df_freq[!,string(baseline) * "_lb"] = q_lb[:,length(subjects)+1]
        df_freq[!,string(baseline) * "_ub"] = q_ub[:,length(subjects)+1]
        df_freq[!,"average_c_lb"] = vec_average_c_lb
        df_freq[!,"average_c_ub"] = vec_average_c_ub
        df_freq[!,"average_k_lb"] = vec_average_k_lb
        df_freq[!,"average_k_ub"] = vec_average_k_ub
    end
    if outfile_prefix==""
        outfile_frequency = "frequencies.csv"
    else
        outfile_frequency = outfile_prefix * "_frequencies.csv"
    end
    CSV.write(outfile_frequency, df_freq)
end

df_loglikelihood = DataFrame()
df_loglikelihood[!,"maxll"] = [-nmaxll]
df_loglikelihood[!,"num_pars"] = [freedom]
df_loglikelihood[!,"AIC"] = [2*nmaxll + 2*freedom]

if outfile_prefix==""
    outfile_loglikelihood = "loglikelihood.csv"
else
    outfile_loglikelihood = outfile_prefix * "_loglikelihood.csv"
end
CSV.write(outfile_loglikelihood, df_loglikelihood)

println("done")
