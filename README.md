# RelRe

## About


## Required packages
The RelRe program is written in the Julia language, and it uses CSV, Dates, DataFrames, Distributions, NLopt, ArgParse packages. Install these packages using install_packages.jl provided with the RelRe code if your julia do not have them.
```sh
julia install_packages.jl
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
| `-d`, `--division`    | DIVISION 	    | divide a day into the given number of equal periods (type: Int64, default: 1) 
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
| `--maxeval`          | MAXEVAL      | stopping criterion used as maxeval in NLopt (type: Float64, default: 5000000)                |
| `-q`, `--frequency`   |               | calculate the time course of variant frequencies                                             |
| `-s`, `--start`       | START         | start date of the analysis (default: "")                                                     |
| `-t`, `--theta`       | THETA         | scale parameter of gamma distribution for generation time (type: Float64, default: 2.32)     |

## Output files
The RelRe program writes the results into four CSV files, namely, estimates.csv, Dirichlet.csv, loglikelihood.csv, and frequencies.csv. Table below shows the summary of each output file.

|File name        | Description        |
|-----------------|--------------------|
|estimates.csv    |The file gives the maximum likelihood estimates of parameters of c_1,...,c_n,k_1,...,k_n,q_(A_1)(t_(A_1)),...,q_(A_n)(t_(A_n)). Parameters c_1,...,c_n are estimated only when the -g option is given. The 95% confidence intervals of parameters are calculated when the -c option is given. The file also contains the date when each variant was assumed to be introduced into the population.|
|Dirichlet.csv    |The value of D, i.e., the sum of the parameters of the Dirichlet distribution is contained. The file is created only when the -D option is given. The 95% confidence intervals are calculated when the -c option is given.|
|loglikelihood.csv|The file contains the logarithm of the maximum likelihood, the number of free parameters, and the value of AIC.|
|frequencies.csv  |The file contains the maximum likelihood estimates of the relative frequency of each variant for each day in the analysis period. The population average of the relative instantaneous reproduction numbers for each day is also given. The lower and upper bound of the 95% confidence intervals are given when the -c option is given.|

## Sample Datasets

The RelRe package contains sample datasets in the “sample_data” folder. Each sample dataset has a Makefile, and one can run the RelRe program by typing “make” in a terminal after moving into the directories containing the sample dataset.

```sh
cd sample_data/01-SARS-CoV-2-Delta-Japan
make
```

## Usage example

The following is an example command to run the RelRe program. 

```sh
julia --threads 10 RelRe.jl -i Japan_Linenage_Counts.csv -b other -a 3.42 -t 1.36 -c -q -f 90
```
The command runs the RelRe program with 10 threads (--threads) based on the input data (-i) provided in "Japan_Lineage_Counts.csv" using "other" as the baseline variant (-b) assuming generation times follows the gamma distribution with a shape parameter (-a) of 3.42 and a scale parameter (-t) of 1.36. The program outputs estimates and their 95% confidence intervals (-c) of parameters as well as the relative frequencies (-q) of variants for each day up to 90 days (-f) after the final date of the observations.

> **Note**
> The default values for the generation time parameters (alpha, theta) are optimized for SARS-CoV-2. These parameters should be adjusted as needed.
