using CSV, Dates, DataFrames, Distributions, NLopt, ArgParse
#julia --threads 5 RelRe.jl -b baseline -s 2020-01-01 -e 2021-01-01 variant_freq.csv 

#Pull in information from command line arguments
s = ArgParseSettings()
@add_arg_table s begin
    "--in", "-i"
    arg_type = String
    required = true
    help = "input file containing model parameters"
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
    "--num", "-n"
    arg_type = Int64
    default = 100
    help = "total number of variants sampled"
end

parsed_args = parse_args(ARGS, s)
@show parsed_args

#Init variables
const infile = parsed_args["in"]
const outfile_prefix = parsed_args["out"]
const start_date = parsed_args["start"]
const end_date = parsed_args["end"]
const alpha = parsed_args["alpha"]
const theta = parsed_args["theta"]
const division = parsed_args["division"]
const dirichlet = parsed_args["Dirichlet"]
const num = parsed_args["num"]

if parsed_args["len"] == -1
    const len_tr = Int64(ceil(quantile(Gamma(alpha, theta),0.99)))
else
    const len_tr = parsed_args["len"]
end
const delta = 1.0/division

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
println("Loading parameters")
df_params = DataFrame(CSV.File(infile))

#Check whether the second column is a vector of dates
if(typeof(df_params[:,2])!=Vector{Date})
    error("The second column is not a vector of dates")
end

#Check the names of the columns
if(propertynames(df_params)[1]!=:variant)
    error("The name of the first column must be variant")
end
if(propertynames(df_params)[2]!=:date)
    error("The name of the second column must be date")
end
if(propertynames(df_params)[3]!=:c)
    error("The name of the third column must be c")
end
if(propertynames(df_params)[4]!=:k)
    error("The name of the fourth column must be k")
end
if(propertynames(df_params)[5]!=:qt)
    error("The name of the fifth column must be qt")
end

#Check whether the input table has counts of the baseline variant
baseline = Symbol(df_params.variant[1])
println("\nBaseline")
println(baseline) 

#Specify variants as all, and remove baseline from subjects
variants = map(x-> Symbol(x),df_params.variant) 
subjects = filter(x -> x!= baseline, variants)
println("\nSubject clades")
const num_subjects = length(subjects)
@show(subjects)
vec_c = df_params.c[2:(num_subjects + 1)] 
vec_k = df_params.k[2:(num_subjects + 1)] 
vec_qt = df_params.qt[2:(num_subjects + 1)] 

dict_index = Dict{Symbol,Int64}()
map(x -> dict_index[subjects[x]] = x, 1:length(subjects))
dict_index[baseline] = length(subjects)+1

#Determine the date for t_end
if(end_date!="")
    const t_end = Date(end_date)
else
    const t_end = maximum(df_params.date)
end

#Determine the date for t_start
if(start_date!="")
    const t_start = Date(start_date)
    if !(t_start in df_params.date)
        println("error: The start date should be a day in the date column")
        exit(1)
    end
else
    const t_start = minimum(df_params.date)
end

#Check whether the start date is earlier than the end date
if(t_start >= t_end)
    println("error: The start date should be earlier than the end date")
    exit(1)
end

@show(df_params)

#Record the date of variant's first observation during the period 
dict_first = Dict{Symbol,Date}()
map(i -> dict_first[variants[i]]=df_params.date[i], 1:(num_subjects + 1))
vec_t = map(v -> dict_first[v], subjects)

println("\nTime range of analysis")
println("Start: " * Dates.format(t_start, ISODateFormat))
println("End: " * Dates.format(t_end, ISODateFormat))

#Calculate Trajectory
subjects_at_start = (1:num_subjects)[map(x->(dict_first[subjects[x]]==t_start), 1:num_subjects)]
println("\nCalculating trajectory")

q_simul = model_q(vec_c, vec_k, vec_qt, vec_t,
               t_start, t_end, len_tr)

vec_average_k = q_simul * vcat(vec_k, 1.0)
vec_average_c = q_simul * vcat(vec_c, 1.0)
    
println("\nWriting frequencies")
df_freq = DataFrame()
df_freq[!,"date"] = collect(t_start:Day(1):t_end)
map(x -> df_freq[!,string(x)] = q_simul[:,dict_index[x]], subjects)
df_freq[!,string(baseline)] = q_simul[:,length(subjects)+1]
df_freq[!,"average_c"] = vec_average_c
df_freq[!,"average_k"] = vec_average_k

if outfile_prefix==""
    outfile_frequency = "pop_freq.csv"
else
    outfile_frequency = outfile_prefix * "_pop_freq.csv"
end
CSV.write(outfile_frequency, df_freq)

num_rows = size(q_simul, 1) 
num_cols = size(q_simul, 2)

new_matrix = zeros(Int64, num_rows, num_cols)

for i in 1:num_rows
    new_matrix[i,:] = rand(Multinomial(num,q_simul[i,:]))
end

df_sim_counts = DataFrame(new_matrix, :auto)
new_names = names(df_freq)[2:(end - 2)]
rename!(df_sim_counts, Symbol.(names(df_sim_counts)) .=> new_names)
date = df_freq[:, 1]
date = DataFrame(date_from=date,date_till=date)
df_sim_counts = hcat(date, df_sim_counts)
CSV.write("counts.csv", df_sim_counts)
println("done")
