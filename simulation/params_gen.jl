using Distributions,ArgParse ,CSV , DataFrames, Printf
const n = parse(Int64, ARGS[1])
const min_k = 0.9
const max_k = 2
const min_qt = 0.01
const max_qt = 0.2
const qt_baseline = 0.8 
start_date = "2024-01-01"

date = fill(start_date,n)
c = fill(1.0,n)
k = vcat(1.0,rand(Uniform(min_k, max_k), n-1))
qt_1 = rand(Uniform(min_qt,max_qt), n-1)
qt_2 = qt_1/sum(qt_1) * (1.0-qt_baseline)
qt = vcat(qt_baseline,qt_2)

variants = map(x -> @sprintf("v%03i", x), 0:(n-1))

df = DataFrame(variant=variants,date=date,c=c,k=k,qt=qt)
CSV.write("parameters.csv",df)
