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
    "--out", "-o"
    arg_type = String
    default = ""
    "--start", "-s"
    arg_type = String
    default = ""
    "--end", "-e"
    arg_type = String
    default = ""
    "--future", "-f"
    arg_type = Int64
    default = 0
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
    "--delta", "-d"
    arg_type = Float64
    default = 0.5
    "--unit", "-u"
    arg_type = Symbol
    default = :D
    required = true
    "--frequency", "-q"
    action = :store_true
    "--estimate_GT", "-g"
    action = :store_true
    "--estimate_CI", "-c"
    action = :store_true
    "--undetected", "-n"
    action = :store_true
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
const breaks = parsed_args["breaks"]
const arg_subjects = parsed_args["subjects"]
const estimate_GT = parsed_args["estimate_GT"]
const estimate_CI = parsed_args["estimate_CI"]
const assume_undetected = parsed_args["undetected"]
const calculate_q = parsed_args["frequency"]
const delta = parsed_args["delta"]

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
                 vec_b::Vector{Date},t_s, t_e, l)
    @assert length(vec_c) == num_subjects
    @assert length(vec_k) == num_subjects * (length(vec_b) + 1)
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
        window_index = 1
        for b in vec_b # the first days of next windows
            if t < b
                break
            else
                window_index += 1
            end
        end
        subvec_k = map(j -> vec_k[(window_index - 1) * num_subjects + j],1:num_subjects)
        
        fill!(vec_sum_nmr, 0.0)
        sum_dnm = 0.0
        for k in 1:length(delta:delta:l)#here
            t_k = max(1, i - k)
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
    println("The analysis period was divided into "
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

mat_95CI = Matrix{Float64}(undef,
                           num_subjects * (2 + 2 * num_windows + 2),
                           num_subjects * (1 + 1 * num_windows + 1))
err_95CI = Vector{Symbol}(undef,num_subjects * (2 + 2 * num_windows + 2))

if estimate_CI
    println("\nCalculating 95% confidence intervals (CIs)")
    Threads.@threads for i in 1:(num_subjects * (2 + 2 * num_windows + 2))
        println("Thread " * string(Threads.threadid()) * " is working on the " *
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
if num_windows > 1
    for w in 1:num_windows
        df_estimates[!,"k_w"*string(w)] = Vector{Float64}()
        if estimate_CI
            df_estimates[!,"k_w"*string(w)*"_lb"] = Vector{Float64}()
            df_estimates[!,"k_w"*string(w)*"_ub"] = Vector{Float64}()
        end
    end
else
    df_estimates[!,"k"] = Vector{Float64}()
    if estimate_CI
        df_estimates[!,"k_lb"] = Vector{Float64}()
        df_estimates[!,"k_ub"] = Vector{Float64}()
    end
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
    if num_windows > 1
        for w in 1:num_windows
            row=vcat(row, par_maxll[num_subjects+(w-1)*num_subjects + j])#k
            if estimate_CI
                row=vcat(row,
                         minimum(mat_95CI[:,num_subjects+(w-1)*num_subjects+j]),#lb
                         maximum(mat_95CI[:,num_subjects+(w-1)*num_subjects+j]))#ub
            end
        end
        row=vcat(row, par_maxll[(1+num_windows)*num_subjects+j]) #qt
        if estimate_CI
            row=vcat(row, minimum(mat_95CI[:,(1+num_windows)*num_subjects+j]),#lb
                     maximum(mat_95CI[:,(1+num_windows)*num_subjects+j]))#ub
        end
    else
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
    end
    push!(df_estimates,row)
end

if outfile_prefix==""
    outfile_estimate = "estimates.csv"
else
    outfile_estimate = outfile_prefix * "_estimates.csv"
end

CSV.write(outfile_estimate, df_estimates)

#Calculate Trajectory

if estimate_GT
    freedom = num_subjects * (1 + num_windows + 1)
else 
    freedom = num_subjects * (num_windows + 1)
end

function constr_trajectory(vec_par, vec_grad)
    -quantile(Chisq(freedom),0.95)/2 - nmaxll + negLogL(vec_par, vec_grad);
end

mat_95CR = Matrix{Float64}(undef,
                           num_subjects * (2 + 2 * num_windows + 2),
                           num_subjects * (1 + 1 * num_windows + 1))
err_95CR = Vector{Symbol}(undef,num_subjects * (2 + 2 * num_windows + 2))

if calculate_q
    println("\nCalculating trajectory")
    vec_c_ml = par_maxll[1:num_subjects] #c
    vec_k_ml = par_maxll[num_subjects+1:num_subjects*(1+num_windows)] #k
    vec_qt_ml = par_maxll[num_subjects*(1+num_windows)+1:num_subjects*(2+num_windows)] #qt
    vec_t_ml = map(v -> dict_first[v], subjects)
    vec_b_ml = breaks

    t_future = t_end + Day(days_to_predict)
    q_ml = model_q(vec_c_ml, vec_k_ml, vec_qt_ml, vec_t_ml, vec_b_ml,
                   t_start, t_future, len_tr)
    q_lb = q_ml
    q_ub = q_ml
    vec_average_c_ml = q_ml * vcat(vec_c_ml, 1.0)
    vec_average_c_lb = vec_average_c_ml
    vec_average_c_ub = vec_average_c_ml

    if num_windows > 1
        vec_average_k_ml = Vector{Float64}(undef,0)
        row_from = 1
        for w in 1:length(breaks)
            local subvec_k = map(j -> vec_k_ml[(w - 1) * num_subjects + j],
                           1:num_subjects)
            global row_to = (breaks[w]-t_start).value
            global vec_average_k_ml =
                vcat(vec_average_k_ml,
                     q_ml[row_from:row_to,:] * vcat(subvec_k, 1.0))
            global row_from = row_to +1
        end
        row_to= (t_future - t_start).value + 1
        local subvec_k = map(j -> vec_k_ml[length(breaks) * num_subjects + j],
                             1:num_subjects)
        vec_average_k_ml = vcat(vec_average_k_ml,
                                q_ml[row_from:row_to,:] * vcat(subvec_k, 1.0))        else
        vec_average_k_ml = q_ml * vcat(vec_k_ml, 1.0)
    end
    vec_average_k_lb = vec_average_k_ml
    vec_average_k_ub = vec_average_k_ml
        
    if estimate_CI
        println("Calculating 95% confidence region (CR)")
        Threads.@threads for i in 1:(num_subjects * (2 + 2 * num_windows + 2))
            println("Thread " * string(Threads.threadid()) *
                " is working on the " *
                string(i) * " th loop of the CR calculation")
            opt_c = Opt(:AUGLAG, length(par_maxll))
            opt_c.lower_bounds = par_lb
            opt_c.upper_bounds = par_ub
            opt_c.maxeval = 500000
            opt_c.ftol_abs = ftol_prec
            
            inequality_constraint!(opt_c, (par,grad) -> constr_trajectory(par,grad),
                                   1e-6)
            
            opt_l = NLopt.Opt(:LN_SBPLX, length(par_maxll))
            #opt_l.xtol_rel = 1e-6 #todo precision 
            opt_l.xtol_rel = 1e-4 #todo precision
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

        for i in 1:(num_subjects * (2 + 2 * num_windows + 2))
            par_cr = mat_95CR[i, :]
            vec_c_cr = par_cr[1:num_subjects] #c
            vec_k_cr = par_cr[num_subjects+1:num_subjects*(1+num_windows)] #k
            vec_qt_cr= par_cr[num_subjects*(1+num_windows)+1:num_subjects*(2+num_windows)] #qt
            vec_t_cr = map(v -> dict_first[v], subjects)
            vec_b_cr = breaks

            q_cr = model_q(vec_c_cr, vec_k_cr, vec_qt_cr, vec_t_cr, vec_b_cr,
                           t_start, t_future, len_tr)
            local vec_average_c_cr = q_cr * vcat(vec_c_cr, 1.0)

            if num_windows > 1
                vec_average_k_cr = Vector{Float64}(undef,0)
                global row_from = 1
                for w in 1:length(breaks)
                    local subvec_k =
                        map(j -> vec_k_cr[(w - 1) * num_subjects + j],
                            1:num_subjects)
                    global row_to = (breaks[w]-t_start).value
                    vec_average_k_cr =
                        vcat(vec_average_k_cr,
                             q_cr[row_from:row_to,:] * vcat(subvec_k, 1.0))
                    global row_from = row_to +1
                end
                global row_to= (t_future - t_start).value + 1
                local subvec_k =
                    map(j -> vec_k_cr[length(breaks) * num_subjects + j],
                        1:num_subjects)
                vec_average_k_cr =
                    vcat(vec_average_k_cr,
                         q_cr[row_from:row_to,:] * vcat(subvec_k, 1.0))
            else
                vec_average_k_cr = q_cr * vcat(vec_k_cr, 1.0)
            end
            
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
df_loglikelihood[!,"AIC"] = [2*freedom + 2*nmaxll]

if outfile_prefix==""
    outfile_loglikelihood = "loglikelihood.csv"
else
    outfile_loglikelihood = outfile_prefix * "_loglikelihood.csv"
end
CSV.write(outfile_loglikelihood, df_loglikelihood)

println("done")
