# RelRe

## About


## Required packages
The RelRe program is written in the Julia language, and it uses CSV, Dates, DataFrames, Distributions, NLopt, ArgParse packages. Install these packages if your julia do not have them.
```julia
import Pkg; Pkg.add(["CSV", "Dates", "DataFrames", "Distributions", "NLopt", "ArgParse"])
```

## The format of input files
The RelRe program calculates relative instantaneous reproduction numbers of subject variants with respect to a baseline variant. It uses observed counts of these variants in time periods, and these count data are assumed to be given to the program in the CSV format. The below is an example of input files. 

```csv
date_from, date_till, Delta, BA.1, BA.2 
2021-12-01, 2021-12-03, 7399, 106, 0
2021-12-04, 2021-12-06, 10630, 368, 1
2021-12-07, 2021-12-09, 4276, 484, 0
2021-12-10, 2021-12-12, 3330, 618, 4
2021-12-13, 2021-12-15, 2055, 1233, 36
2021-12-16, 2021-12-18, 1145, 1159, 67
...
```
The first line of the CSV file should indicate the column names. The first two columns should be date_from and date_till. Other columns should be names of variants. Each line except the first represent observed counts of variants in the period starting at the date in the date_from column and ending at the date in the date_till column. The duration of each observation period can be set freely. Observation periods may have different durations. The minimum observation period is one day, and in this case, date_from and date_till have the same date. Dates are assumed to be given in the YYYY-MM-DD format. 


## Params

Currently the following parameters are recognized

| parameter             | variable      | description                                                                                  |
|-----------------------|---------------|----------------------------------------------------------------------------------------------|
| `-a`, `--alpha`       | ALPHA         | shape parameter of gamma distribution for generation time (type: Float64, default: 2.03)     |
| `-b`, `--baseline`    | BASELINE      | variant used as the baseline of relative reproduction numbers                                |
| `-c`, `--estimate_CI` |               | estimate 95% confidence intervals                                                            |
| `-d`, `--delta`       | DELTA    	    | unit time of calculation (in days) (type:Float64, default: 0.5)                              |
| `-D`, `--Dirichlet`   |  	            | use Dirichlet multinomial as the observation model                                           |
| `-e`, `--end`         | END           | end date of the analysis (default: "")                                                       |
| `-f`, `--future` 	    | FUTURE        | duration in days for predicting variant frequencies (type: Int64, default: 0)                |
| `-g`, `--estimate_GT` |               | estimate relative generation times of variants                                               |
| `-h`, `--help`        |               | show this help message and exit                                                              |
| `-i`, `--in`          | IN            | input file containing temporal count data of variants                                        |
| `-j`, `--subjects`    | [SUBJECTS...]	| list of variants to calculate relative reproduction numbers (type: Symbol)                   |
| `-l`, `--len`         | LEN         	| trancation point of gamma distribution for generation time (type: Int64, default: automatic) |
| `-u`, `--undetected`  |               | assume all variants exist undetected from the start date                                     |
| `-o`, `--out`         | OUT         	| prefix of output files (default: "")                                                         |
| `--ftol_abs`          | FTOL_ABS      |  stopping criterion used as ftol_abs in NLopt (type: Float64, default: 0.0)                  |
| `--ftol_rel`          | FTOL_REL      | stopping criterion used as ftol_rel in NLopt (type: Float64, default: 1.0e-8)                |
| `-q`, `--frequency`   |               | calculate the time course of variant frequencies                                             |
| `-s`, `--start`       | START         | start date of the analysis (default: "")                                                     |
| `-t`, `--theta`       | THETA         | scale parameter of gamma distribution for generation time (type: Float64, default: 2.32)     |

## Usage example

```sh
julia --threads 10 RelRe.jl -b Omicron_BA1 -a 2.03 -t 1.392 -i Tokyo_BA1_BA2.csv -c -q -d 1.0 -f 30 -D -u
```

> **Note**
> This version of the code is optimized for Sars-Cov-2. Parameters should be adjusted as needed.