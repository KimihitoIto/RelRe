using CSV, Dates, DataFrames, Distributions, NLopt, ArgParse
#julia --threads 5 RelRe.jl -b baseline -s 2020-01-01 -e 2021-01-01 variant_freq.csv 

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
    # add --values, -v option for giving initial values of optimisation. Set the default value to "".
    "--values", "-v"
    arg_type = String
    default = ""
    help = "for giving initial values of optimisation"
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
    "--ftol_abs"
    arg_type = Float64
    default = 0.0
    help = "stopping criterion used as ftol_abs in NLopt"
    "--ftol_rel"
    arg_type = Float64
    default = 1e-8
    help = "stopping criterion used as ftol_rel in NLopt"
    "--maxeval"
    arg_type = Int64
    default = 5_000_000
    help = "stopping criterion used as maxeval in NLopt"
    "--len", "-l"
    arg_type = Int64
    default = -1 # use 7 for flu; automatically calculated if not given
    help = "truncation point of gamma distribution for generation time"
    "--alpha", "-a"      # TODO: estimate from mean and var 
    arg_type = Float64
    default = 2.03 # use 4.5 for flu
    help = "shape parameter of gamma distribution for generation time"
    "--theta", "-t"      # TODO: estimate from mean and var
    arg_type = Float64
    default = 2.32 # use 0.60 for flu
    help = "scale parameter of gamma distribution for generation time"
    "--division", "-d"
    arg_type = Int64
    default = 1
    help = "the number of calculation time steps in one day"
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
    "--undetected", "-u"
    action = :store_true
    help = "assume all variants exist undetected from the start date"
    "--blocks", "-r"
    arg_type = Int64 
    nargs = '*'
    default = []
    help = "a vector indicating variants sharing same reproduction numbers"
end
parsed_args = parse_args(ARGS, s)
@show parsed_args
#Init variables
const infile = parsed_args["in"]
const outfile_prefix = parsed_args["out"]
const valuefile = parsed_args["values"]
const ftol_abs = parsed_args["ftol_abs"]
const ftol_rel = parsed_args["ftol_rel"]
const maxeval = parsed_args["maxeval"]
const baseline = parsed_args["baseline"]
const start_date = parsed_args["start"]
const end_date = parsed_args["end"]
const days_to_predict = parsed_args["future"]
const alpha = parsed_args["alpha"]
const theta = parsed_args["theta"]
const arg_subjects = parsed_args["subjects"]
const estimate_GT = parsed_args["estimate_GT"]
const estimate_CI = parsed_args["estimate_CI"]
const assume_undetected = parsed_args["undetected"]
const calculate_q = parsed_args["frequency"]
const division = parsed_args["division"]
const dirichlet = parsed_args["Dirichlet"]
if parsed_args["len"] == -1
    const len_tr = Int64(ceil(quantile(Gamma(alpha, theta),0.99)))
else
    const len_tr = parsed_args["len"]
end
const delta = 1.0/division
const arg_blocks = parsed_args["blocks"]

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

#Renewal model of variant frequencies
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
        q[1, num_subjects + 1] = 0.0
    else
        q[1, num_subjects + 1] = 1.0 - sum_q_subjects
    end
    
    vec_sum_nmr=Vector{Float64}(undef, num_subjects)
    
    for i in 2:length(delta:delta:duration)
        if(delta == 1.0)
            t = t_s + Day(i-1)
        else
            t = t_s + Day(Int64(floor((delta:delta:duration)[i])))
        end
        
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
        
        if(sum_q_subjects > 1.0)
            map(j -> q[i,j] /= sum_q_subjects,1:num_subjects)
            q[i, num_subjects + 1] = 0.0
        else
            q[i, num_subjects + 1] = 1.0 - sum_q_subjects
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

#Check whether the first two columns are vectors of dates
if(typeof(df_count[:,1])!=Vector{Date})
    error("The first column is not a vector of dates")
end
if(typeof(df_count[:,2])!=Vector{Date})
    error("The second column is not a vector of dates")
end

#Check the names of the first two columns
if(propertynames(df_count)[1]!=:date_from)
    error("The name of the first column must be date_from")
end
if(propertynames(df_count)[2]!=:date_till)
    error("The name of the second column must be date_till")
end

#Check whether the input table has counts of the baseline variant
if !(baseline in propertynames(df_count)[3:size(df_count,2)])
    error("Counts of the baseline variant can't be found in the table")
end
println("\nBaseline")
println(baseline) # todo: baseline need to be checked

#Specify variants as all, and remove baseline from subjects
if(arg_subjects == [])
    variants = propertynames(df_count)[3:size(df_count,2)]
    subjects = filter(x -> x!= baseline, variants)
else
    if !issubset(arg_subjects,propertynames(df_count)[3:size(df_count,2)])
        error("Counts of a subject can't be found in the table")
    else
        subjects = arg_subjects
        variants = vcat(subjects,baseline)
    end
end
println("\nSubject clades")
const num_subjects = length(subjects)
@show(subjects)

dict_index = Dict{Symbol,Int64}()
map(x -> dict_index[subjects[x]] = x, 1:length(subjects))
dict_index[baseline] = length(subjects)+1

#Determine the date for t_end
if(end_date!="")
    const t_end = Date(end_date)
    if !(t_end in df_count.date_till)
        println("error: The end date should be a day in the date_till column")
        exit(1)
    end
else
    const t_end = maximum(df_count.date_till)
end

#Determine the date for t_start
if(start_date!="")
    const t_start = Date(start_date)
    if !(t_start in df_count.date_from)
        println("error: The start date should be a day in the date_from column")
        exit(1)
    end
else
    const t_start = minimum(df_count.date_from)
end

#Check whether the start date is earlier than the end date
if(t_start >= t_end)
    println("error: The start date should be earlier than the end date")
    exit(1)
end

#Remove data before start date or after end date
deleteat!(df_count, df_count.date_till .< t_start)
deleteat!(df_count, df_count.date_from .> t_end)

@show(df_count)

#Record the date of variant's first observation during the period 
dict_first = Dict{Symbol,Date}()
if assume_undetected
    map(v -> dict_first[v]=t_start, variants)
else
    map(v -> dict_first[v]=minimum(filter(v => n -> n>0, df_count).date_from), variants)
end

println("\nTime range of analysis")
println("Start: " * Dates.format(t_start, ISODateFormat))
println("End: " * Dates.format(t_end, ISODateFormat))

dates_from = df_count.date_from
dates_till = df_count.date_till
mat_obs = Matrix(df_count[:,vcat(subjects,baseline)])

#read initial values from the file if the file name is given 
if valuefile != ""
    @assert isfile(valuefile) == true 
    println("Loading initial values")
    df_values = DataFrame(CSV.File(valuefile))
end

function negLogL(par::Vector, grad::Vector)
    if(dirichlet)
        @assert length(par) == 2 * num_subjects + num_k + 1
    else
        @assert length(par) == 2 * num_subjects + num_k
    end
    vec_c = par[1:num_subjects]
    vec_k_tmp = par[num_subjects+1:num_subjects+num_k]
    vec_k = map(x -> x==0 ? 1.0 : vec_k_tmp[x], vec_blocks)
    vec_qt = par[num_subjects+num_k+1:num_subjects*2+num_k]
    if(dirichlet)
        D = par[num_subjects * 2 + num_k + 1]
    end
    
    vec_t = map(v -> dict_first[v], subjects)
    try
        q = model_q(vec_c, vec_k, vec_qt, vec_t, t_start, t_end, len_tr)
        sumll = 0.0
        for i in 1:length(dates_from)
            obs = mat_obs[i,:]
            if(sum(obs)==0)
                continue
            end
            j1 = Dates.value(dates_from[i]-t_start)+1
            j2 = Dates.value(dates_till[i]-t_start)+1
            rows = collect(j1:j2)
            probs = vec(max.(0, mean(q[rows, 1:num_subjects+1], dims=1)))
            if(dirichlet)
                indices = 1:(num_subjects+1)
                non_zero_indices = indices[map(x->probs[x]!=0.0,indices)]
                alphas = probs[non_zero_indices] * D
                sumll += logpdf(DirichletMultinomial(sum(obs[non_zero_indices]),
                                                     alphas),
                                obs[non_zero_indices])
            else
                sumll += logpdf(Multinomial(sum(obs), probs), obs)
            end
        end
        if !isfinite(sumll)
            println("Warning: sumll is not finite")
            return floatmax(Float64)
        end
        return -sumll
    catch e
        println(e)
        exit(1)
    end
end

function is_restricted_growth(vec::Vector{Int64})
  if vec[1] != 0 
    return false
  end
  for i in 2:length(vec)
    if (vec[i] > 1+maximum(vec[1:(i-1)]) || vec[i] < 0)
      return false
    end
  end
    return true
end

if length(arg_blocks) == 0
  vec_blocks = collect(1:num_subjects)
else
  if length(arg_blocks) != num_subjects + 1
    println("Error: The length of the blocks vector should be the same as the number of variants")
    exit(1)
  end 
  if maximum(arg_blocks) > num_subjects
    println("Error:The maximum number in the blocks vector should be smaller than or equal to the number of subjects")
    exit(1)
  end
  if minimum(arg_blocks) != 0 
    println("Error:The minimum number in the blocks vector should be 0")
    exit(1)
  end
  if is_restricted_growth(arg_blocks) == false
    println("Error: The blocks vector should be a restricted growth string")
    exit(1)
  end
  vec_blocks = arg_blocks[2:length(arg_blocks)]
end
@show vec_blocks
num_k = maximum(vec_blocks)

vec_c_start = fill(1.0,num_subjects)

if estimate_GT
    if Dates.value(maximum(dates_till-dates_from)) < 7
        vec_c_lb = fill(1.0e-10,num_subjects)
        vec_c_ub = fill(10.0,num_subjects)
    else
        println("error: The bin width should be a day or week for the -g option")
        exit(1)
    end
else
    vec_c_lb = fill(1.0,num_subjects)
    vec_c_ub = fill(1.0,num_subjects)
end
#modify the following code to use values read from the value file
#use mean values for initial values of a group
# note that the length of the vector is num_k and not num_subjects
if valuefile != ""
    k_vals = df_values.k[2:end]
    vec_k_start = fill(1.0, num_k)
    for i in 1:num_k
        indices = findall(x -> x == i, vec_blocks)
        mean_val = mean(k_vals[indices])
        vec_k_start[i] = mean_val
    end
    vec_k_lb = fill(1.0e-10,num_k) # no need to change 
    vec_k_ub = fill(10.0,num_k)# no need to change 
    #modify the following code to use values read from the value file
    vec_qt_start = df_values.qt[2:end]
    vec_qt_lb = fill(1.0e-10,num_subjects)# no need to change 
    vec_qt_ub = fill(1.0,num_subjects)# no need to change 
    #Check if the vec_k_start and vec_qt_start were correct.
    @show vec_k_start vec_qt_start
else
    vec_k_start = fill(1.0,num_k)
    vec_k_lb = fill(1.0e-10,num_k)
    vec_k_ub = fill(10.0,num_k)
    vec_qt_start = fill(0.001,num_subjects)
    vec_qt_lb = fill(1.0e-10,num_subjects)
    vec_qt_ub = fill(1.0,num_subjects)
end

if(dirichlet)
    D_start = 10
    D_lb = 1.0e-10
    D_ub = 1.0e+5
    par_start = vcat(vec_c_start, vec_k_start, vec_qt_start, D_start)
    par_lb = vcat(vec_c_lb, vec_k_lb, vec_qt_lb, D_lb)
    par_ub = vcat(vec_c_ub, vec_k_ub, vec_qt_ub, D_ub)
else
    par_start = vcat(vec_c_start, vec_k_start, vec_qt_start)
    par_lb = vcat(vec_c_lb, vec_k_lb, vec_qt_lb)
    par_ub = vcat(vec_c_ub, vec_k_ub, vec_qt_ub)
end

opt = Opt(:LN_SBPLX, length(par_start))
opt.min_objective = (par, grad) -> negLogL(par, grad)
opt.lower_bounds = par_lb
opt.upper_bounds = par_ub
opt.ftol_abs = ftol_abs
opt.ftol_rel = ftol_rel
opt.maxeval = maxeval

println("Maximizing the likelihood function")
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
                               num_subjects * 4 + num_k * 2 + 2,
                               num_subjects * 2 + num_k + 1)
    err_95CI = Vector{Symbol}(undef,num_subjects * 4 + num_k * 2 + 2)
else
    mat_95CI = Matrix{Float64}(undef,
                               num_subjects * 4 + num_k * 2,
                               num_subjects * 2 + num_k)
    err_95CI = Vector{Symbol}(undef,num_subjects * 4 + num_k * 2)
end

if estimate_CI
    println("\nCalculating 95% confidence intervals (CIs)")
    if(dirichlet)
        num_loop = (num_subjects * 4 + num_k * 2 + 2)
    else
        num_loop = (num_subjects * 4 + num_k * 2)
    end
    
    Threads.@threads for i in 1:num_loop
        println("Thread " * string(Threads.threadid()) * " is working on the " *
            string(i) * " th loop of the CI calculation")
        
        opt_c = Opt(:AUGLAG, length(par_maxll))
        opt_c.lower_bounds = par_lb
        opt_c.upper_bounds = par_ub
        opt_c.maxeval = maxeval
        opt_c.ftol_abs = ftol_abs
        opt_c.ftol_rel = ftol_rel
        
        inequality_constraint!(opt_c, (par,grad) -> constr(par,grad), 1e-6)
        
        opt_l = NLopt.Opt(:LN_SBPLX, length(par_maxll))
        opt_l.ftol_abs = ftol_abs
        opt_l.ftol_rel = ftol_rel
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
println("\nWriting estimates")
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

#baseline
row_bl = vcat(baseline, t_start, 1.0) #c
if estimate_CI
    row_bl = vcat(row_bl, 1.0, 1.0) #lb, ub
end
row_bl=vcat(row_bl, 1.0)#k
if estimate_CI
    row_bl = vcat(row_bl, 1.0, 1.0) #lb, ub
end

subjects_at_start = (1:num_subjects)[map(x->(dict_first[subjects[x]]==t_start), 1:num_subjects)]

row_bl=vcat(row_bl,
            1.0-sum(map(j->par_maxll[num_subjects+num_k+j],subjects_at_start)))#qt
if estimate_CI
    sum_qs = sum(mat_95CI[:,(num_subjects + num_k) .+ subjects_at_start],dims=2)
    lb = 1.0 - maximum(sum_qs)
    ub = 1.0 - minimum(sum_qs)
    row_bl=vcat(row_bl,lb, ub)
end
push!(df_estimates,row_bl)

#subjects
for j in 1:num_subjects
    row = vcat(subjects[j], dict_first[subjects[j]], par_maxll[j]) #c
    if estimate_CI
        row = vcat(row, minimum(mat_95CI[:,j]), maximum(mat_95CI[:,j])) #lb, ub
    end
    row=vcat(row, vec_blocks[j]==0 ? 1.0 : par_maxll[num_subjects+vec_blocks[j]])#k
    if estimate_CI
        row=vcat(row, minimum(mat_95CI[:,num_subjects + vec_blocks[j]]), #lb
                 maximum(mat_95CI[:,num_subjects + vec_blocks[j]])) #ub
    end
    row=vcat(row, par_maxll[num_subjects+num_k+j])#qt
    if estimate_CI
        row=vcat(row, minimum(mat_95CI[:,num_subjects + num_k + j]), #lb
                 maximum(mat_95CI[:,num_subjects + num_k + j])) #ub
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
    df_D = DataFrame()
    df_D[!,"D"] = [par_maxll[2 * num_subjects + num_k + 1]]
    if estimate_CI
        df_D[!,"D_lb"] = [minimum(mat_95CI[:,2 * num_subjects + num_k + 1])]
        df_D[!,"D_ub"] = [maximum(mat_95CI[:,2 * num_subjects + num_k + 1])]
    end
    if outfile_prefix==""
        outfile_D = "Dirichlet.csv"
    else
        outfile_D = outfile_prefix * "_Dirichlet.csv"
    end
    CSV.write(outfile_D, df_D)
end

#Calculate Trajectory

if estimate_GT
    freedom = num_subjects * 2 + num_k
else 
    freedom = num_subjects + num_k
end

if dirichlet
    freedom = freedom + 1
end

function constr_trajectory(vec_par, vec_grad)
    -quantile(Chisq(freedom),0.95)/2 - nmaxll + negLogL(vec_par, vec_grad);
end

if(dirichlet)
    mat_95CR = Matrix{Float64}(undef,
                               num_subjects * 4 + num_k * 2 + 2,
                               num_subjects * 2 + num_k + 1)
    err_95CR = Vector{Symbol}(undef,num_subjects * 4 + num_k * 2 + 2)
else
    mat_95CR = Matrix{Float64}(undef,
                               num_subjects * 4 + num_k * 2,
                               num_subjects * 2 + num_k)
    err_95CR = Vector{Symbol}(undef,num_subjects * 4 + num_k * 2)
end

if calculate_q
    println("\nCalculating trajectory")
    vec_c_ml = par_maxll[1 : num_subjects] #c
    vec_k_ml = par_maxll[num_subjects .+ vec_blocks] #k
    vec_qt_ml = par_maxll[(num_subjects + num_k +1) : (num_subjects*2 + num_k)] #qt
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
            num_loop = (num_subjects * 4 + num_k * 2 + 2)
        else
            num_loop = (num_subjects * 4 + num_k * 2)
        end
        Threads.@threads for i in 1:num_loop
            println("Thread " * string(Threads.threadid()) *
                " is working on the " *
                string(i) * " th loop of the CR calculation")
            opt_c = Opt(:AUGLAG, length(par_maxll))
            opt_c.lower_bounds = par_lb
            opt_c.upper_bounds = par_ub
            opt_c.maxeval = maxeval
            opt_c.ftol_abs = ftol_abs
            opt_c.ftol_rel = ftol_rel
            
            inequality_constraint!(opt_c, (par,grad) -> constr_trajectory(par,grad),0)
            opt_l = NLopt.Opt(:LN_SBPLX, length(par_maxll))
            opt_l.ftol_abs = ftol_abs
            opt_l.ftol_rel = ftol_rel
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
            vec_k_cr = par_cr[num_subjects .+ vec_blocks] #k
            vec_qt_cr= par_cr[num_subjects + num_k + 1 : num_subjects * 2 + num_k] #qt
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

df_blocks = DataFrame()
df_blocks[!,"variant"] = vcat(baseline,subjects)
df_blocks[!, "block"] = vcat(0,vec_blocks)

if outfile_prefix==""
    outfile_blocks = "blocks.csv"
else
    outfile_blocks = outfile_prefix * "_blocks.csv"
end
CSV.write(outfile_blocks, df_blocks)

println("done")
